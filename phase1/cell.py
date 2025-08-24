import random
import math
import statistics
import json
import gzip
import os
import copy
from typing import List, Dict, Any, Tuple, Optional

# Try to import numpy, but make it optional
try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False


def KL_div(P: List[float], Q: List[float]) -> float:
    """Calculate Kullback-Leibler divergence."""
    if len(P) != len(Q):
        raise ValueError("Both distributions must have the same length.") 
    
    n = len(P)
    
    for i in range(n):
        if Q[i] == 0 and P[i] > 0:
            return float('inf')
    
    return sum([P[i] * math.log2(P[i]/Q[i]) for i in range(n) if P[i] > 0])


def JS_div(P: List[float], Q: List[float]) -> float:
    """Calculate Jensen-Shannon divergence."""
    if len(P) != len(Q):
        raise ValueError("Both distributions must have the same length.")

    n = len(P)
    M = [(P[i] + Q[i]) / 2 for i in range(n)]
    return (KL_div(P, M) + KL_div(Q, M)) / 2


# Default parameters
N = 1000
RATE = 0.005
GENE_SIZE = 5
BASELINE_METHYLATION_DISTRIBUTION = [1.0] + [0.0 for _ in range(GENE_SIZE)]
T_MAX = 100
DEFAULT_GROWTH_PHASE = 13  # Default growth phase duration in years


class Cell:
    """
    Represents a cell with methylation sites.
    
    Key methods:
        methylate(): Apply stochastic methylation (formerly age_1_year)
        create_daughter_cell(): Create identical copy (mitosis)
        to_dict(): Convert to dictionary for serialization
    """
    
    def __init__(self, n: int = N, rate: float = None, gene_size: int = GENE_SIZE, 
                 baseline_methylation_distribution: List[float] = None,
                 gene_rate_groups: List[Tuple[int, float]] = None) -> None:
        
        # Validation: Can't specify both rate types
        if rate is not None and gene_rate_groups is not None:
            raise ValueError(
                "Cannot specify both 'rate' and 'gene_rate_groups'. "
                "Use 'rate' for uniform methylation or 'gene_rate_groups' for gene-specific rates."
            )
        
        # Validation: Must specify at least one
        if rate is None and gene_rate_groups is None:
            raise ValueError(
                "Must specify either 'rate' (uniform) or 'gene_rate_groups' (gene-specific)"
            )
        
        self.n = n
        self.gene_size = gene_size
        self.rate = rate  # Keep for backward compatibility
        self.gene_rate_groups = gene_rate_groups  # NEW attribute
        
        # Validate gene_rate_groups if provided
        if gene_rate_groups is not None:
            total_genes = sum(n_genes for n_genes, _ in gene_rate_groups)
            expected_genes = n // gene_size
            if total_genes != expected_genes:
                raise ValueError(
                    f"gene_rate_groups specifies {total_genes} genes, but cell has {expected_genes} genes "
                    f"(n={n}, gene_size={gene_size}). Adjust gene counts to match."
                )
            
            # Validate each group
            for i, (n_genes, rate_val) in enumerate(gene_rate_groups):
                if n_genes <= 0:
                    raise ValueError(f"Group {i}: number of genes must be positive, got {n_genes}")
                if rate_val <= 0:
                    raise ValueError(f"Group {i}: rate must be positive, got {rate_val}")
        
        # Validate n % gene_size
        if self.n % self.gene_size != 0:
            raise ValueError(f"Gene size ({gene_size}) must divide the number of sites ({n}) evenly.")
        
        # Build site_rates array (NEW)
        self._build_site_rates()
        
        # Create baseline distribution matching gene_size if not provided
        if baseline_methylation_distribution is None:
            self.baseline_methylation_distribution = [1.0] + [0.0] * gene_size
        else:
            self.baseline_methylation_distribution = baseline_methylation_distribution
        
        # Initialize methylation pattern
        self.cpg_sites = [0 for _ in range(n)]  
        self.age = 0
        self.methylation_proportion = 0.0
        
        self.methylation_distribution = [0.0 for _ in range(0, self.gene_size + 1)]
        self.methylation_distribution[0] = 1.0  # initially all genes are unmethylated
        self.cell_JSD = JS_div(self.methylation_distribution, self.baseline_methylation_distribution)
    
    def _build_site_rates(self) -> None:
        """Build array of per-site methylation rates based on configuration."""
        if self.rate is not None:
            # Uniform rate (backward compatible)
            if HAS_NUMPY:
                self.site_rates = np.full(self.n, self.rate, dtype=np.float64)
            else:
                self.site_rates = [self.rate] * self.n
        else:
            # Gene-specific rates
            site_rates_list = []
            for n_genes, rate_val in self.gene_rate_groups:
                sites_for_group = n_genes * self.gene_size
                site_rates_list.extend([rate_val] * sites_for_group)
            
            if HAS_NUMPY:
                self.site_rates = np.array(site_rates_list, dtype=np.float64)
            else:
                self.site_rates = site_rates_list
        
    def methylate(self) -> None:
        """
        Apply stochastic methylation to unmethylated sites using per-site rates.
        This is the core aging mechanism (formerly age_1_year).
        """
        self.age += 1
        
        # Early exit if already fully methylated
        if self.methylation_proportion >= 1.0:
            return
        
        methylated_count = 0
        
        if HAS_NUMPY and isinstance(self.site_rates, np.ndarray):
            # Numpy vectorized implementation
            unmethylated_mask = np.array(self.cpg_sites) == 0
            random_vals = np.random.random(self.n)
            new_methylations = unmethylated_mask & (random_vals < self.site_rates)
            
            for i in np.where(new_methylations)[0]:
                self.cpg_sites[i] = 1
            
            methylated_count = sum(self.cpg_sites)
        else:
            # Pure Python implementation
            for i in range(self.n):
                if self.cpg_sites[i] == 0:
                    if random.random() < self.site_rates[i]:  # Use per-site rate
                        self.cpg_sites[i] = 1
                        methylated_count += 1
                else:
                    methylated_count += 1
                
        self.methylation_proportion = methylated_count / self.n
        self.compute_methylation_distribution()
        self.cell_JSD = JS_div(self.methylation_distribution, self.baseline_methylation_distribution)
    
    def create_daughter_cell(self) -> 'Cell':
        """
        Create an identical copy of this cell (mitosis).
        Both daughter cells inherit the parent's methylation pattern.
        """
        return copy.deepcopy(self)
            
    def compute_methylation_distribution(self) -> None:
        """Calculate the distribution of methylation levels across genes."""
        distribution = [0.0 for _ in range(0, self.gene_size + 1)]
        
        # Count methylated sites per gene
        for i in range(0, self.n, self.gene_size):
            methylated_count = 0
            for j in range(i, min(i + self.gene_size, self.n)):
                if self.cpg_sites[j] == 1:
                    methylated_count += 1
            distribution[methylated_count] += 1
            
        total_mass = sum(distribution)
        if total_mass > 0:
            distribution = [x / total_mass for x in distribution]
        
        self.methylation_distribution = distribution
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert cell state to dictionary for serialization."""
        data = {
            'cpg_sites': self.cpg_sites[:],
            'methylation_proportion': self.methylation_proportion,
            'methylation_distribution': self.methylation_distribution[:],
            'cell_jsd': self.cell_JSD,
            'age': self.age,
            'gene_size': self.gene_size
        }
        
        # Save rate configuration
        if self.rate is not None:
            data['rate'] = self.rate
        else:
            data['gene_rate_groups'] = self.gene_rate_groups
        
        # Save site_rates for faster loading
        if HAS_NUMPY and isinstance(self.site_rates, np.ndarray):
            data['site_rates'] = self.site_rates.tolist()
        else:
            data['site_rates'] = self.site_rates
        
        return data


class PetriDish:
    """
    Manages a population of cells with division, methylation, and homeostasis.
    
    Life cycle:
        - Years 0-13: Growth phase (1 → 8192 cells)
        - Years 14+: Steady state (maintain ~8192 cells)
    
    Key methods:
        divide_cells(): Each cell divides into two
        methylate_cells(): Apply methylation to all cells
        random_cull_cells(): Randomly remove ~50% of cells
        simulate_year(): Run one year of simulation
        run_simulation(): Run complete simulation
    """
    
    def __init__(self, rate: float = None, gene_rate_groups: List[Tuple[int, float]] = None,
                 n: int = N, gene_size: int = GENE_SIZE, 
                 seed: int = None, growth_phase: int = DEFAULT_GROWTH_PHASE, 
                 cells: List['Cell'] = None, calculate_jsds: bool = True) -> None:
        """
        Initialize petri dish with a single unmethylated cell or provided cells.
        
        Args:
            rate: Uniform methylation rate per site per year
            gene_rate_groups: Gene-specific rates as [(n_genes, rate), ...]
            n: Number of CpG sites per cell
            gene_size: Number of sites per gene
            seed: Random seed for reproducibility
            growth_phase: Duration of growth phase in years (target = 2^growth_phase cells)
            cells: Optional list of cells to start with (for phase2 compatibility)
            calculate_jsds: Whether to calculate cell and gene JSDs (for performance)
        """
        # Validate rate specification
        if rate is not None and gene_rate_groups is not None:
            raise ValueError(
                "Cannot specify both 'rate' and 'gene_rate_groups'. "
                "Use 'rate' for uniform methylation or 'gene_rate_groups' for gene-specific rates."
            )
        
        if rate is None and gene_rate_groups is None:
            # Use default RATE if neither specified
            rate = RATE
        
        # Validate gene_size divides n evenly
        if n % gene_size != 0:
            raise ValueError(f"n ({n}) must be divisible by gene_size ({gene_size})")
        
        # Validate growth_phase
        if growth_phase < 1:
            raise ValueError(f"growth_phase must be >= 1, got {growth_phase}")
        if growth_phase > 20:
            raise ValueError(f"growth_phase must be <= 20 (max 1M cells), got {growth_phase}")
        
        if seed is not None:
            random.seed(seed)
            
        self.rate = rate
        self.gene_rate_groups = gene_rate_groups
        self.n = n
        self.gene_size = gene_size
        self.n_genes = n // gene_size
        self.seed = seed
        self.growth_phase = growth_phase
        self.target_population = 2 ** growth_phase  # Calculate from growth_phase
        self.calculate_jsds = calculate_jsds
        
        # Initialize cells - either provided or single unmethylated cell
        if cells is not None:
            self.cells = cells
            # Validate all cells have same rate configuration
            self._validate_cell_rate_consistency()
        else:
            # Create initial cell with appropriate rate structure
            if self.rate is not None:
                self.cells = [Cell(n=n, rate=self.rate, gene_size=gene_size)]
            else:
                self.cells = [Cell(n=n, gene_rate_groups=self.gene_rate_groups, gene_size=gene_size)]
        
        self.year = 0
        
        # Cell history tracking (renamed for clarity)
        self.absolute_year = 0  # For proper year labeling in phase2
        self.track_cell_history = True  # Renamed from history_enabled
        
        # Gene JSD tracking
        self.track_gene_jsd = False
        self.gene_jsd_history = {}
        # Baseline must match gene_size + 1 bins
        self.BASELINE_GENE_DISTRIBUTION = [1.0] + [0.0] * self.gene_size
        
        # Store initial state if we have cells (renamed for clarity)
        if self.cells:
            self.cell_history = {'0': [cell.to_dict() for cell in self.cells]}
            # Initialize gene JSD at year 0 if tracking
            if self.track_gene_jsd and self.calculate_jsds:
                self.gene_jsd_history['0'] = [0.0] * self.n_genes  # All zeros initially
        else:
            self.cell_history = {}
        
    def _validate_cell_rate_consistency(self) -> None:
        """Ensure all cells have the same rate configuration."""
        if not self.cells:
            return
        
        first_cell = self.cells[0]
        has_uniform = first_cell.rate is not None
        first_groups = first_cell.gene_rate_groups
        
        for i, cell in enumerate(self.cells[1:], 1):
            cell_has_uniform = cell.rate is not None
            
            if has_uniform != cell_has_uniform:
                raise ValueError(
                    f"Cell rate configuration mismatch: cell 0 has {'uniform' if has_uniform else 'gene-specific'} "
                    f"rates, but cell {i} has {'uniform' if cell_has_uniform else 'gene-specific'} rates. "
                    f"All cells must use the same rate configuration."
                )
            
            if has_uniform:
                if cell.rate != first_cell.rate:
                    raise ValueError(
                        f"Cell rate mismatch: cell 0 has rate={first_cell.rate}, "
                        f"but cell {i} has rate={cell.rate}"
                    )
            else:
                if cell.gene_rate_groups != first_groups:
                    raise ValueError(
                        f"Cell gene_rate_groups mismatch: cell 0 has {first_groups}, "
                        f"but cell {i} has {cell.gene_rate_groups}"
                    )
    
    def divide_cells(self) -> None:
        """
        Cell division: Each cell divides into two identical daughters.
        Population doubles: N → 2N
        """
        new_cells = []
        for cell in self.cells:
            daughter1 = cell.create_daughter_cell()
            daughter2 = cell.create_daughter_cell()
            new_cells.extend([daughter1, daughter2])
        self.cells = new_cells
        print(f"  Division: {len(self.cells)//2} → {len(self.cells)} cells")
        
    def methylate_cells(self) -> None:
        """
        Apply stochastic methylation to all cells.
        Each cell independently accumulates methylation.
        """
        for cell in self.cells:
            cell.methylate()
            # Only calculate JSD if enabled
            if not self.calculate_jsds:
                cell.cell_JSD = 0.0
        print(f"  Methylation applied to {len(self.cells)} cells")
        
    def random_cull_cells(self) -> None:
        """
        Randomly kill cells (each has 50% survival chance).
        Used in steady state to maintain population around target.
        """
        initial_count = len(self.cells)
        survivors = []
        for cell in self.cells:
            if random.random() < 0.5:
                survivors.append(cell)
        self.cells = survivors
        print(f"  Random cull: {initial_count} → {len(self.cells)} cells")
        
    def age_cells_growth_phase(self) -> None:
        """
        Age cells during growth phase.
        Process: divide → methylate
        """
        self.divide_cells()
        self.methylate_cells()
        
    def age_cells_steady_state(self) -> None:
        """
        Age cells during steady state (after growth phase).
        Process: divide → cull → methylate
        """
        self.divide_cells()      # 8192 → 16384
        self.random_cull_cells() # 16384 → ~8192
        self.methylate_cells()   # Apply methylation to survivors
        
    def simulate_year(self) -> None:
        """
        Simulate one year of population dynamics.
        Chooses appropriate aging strategy based on current year.
        """
        self.year += 1
        self.absolute_year += 1  # Keep absolute year in sync
        print(f"\nYear {self.year}:")
        
        if self.year <= self.growth_phase:
            # Growth phase: simple division and methylation
            print(f"  Growth phase (year {self.year} of {self.growth_phase})")
            self.age_cells_growth_phase()
            # Show it's predictable during growth
            expected = 2 ** self.year
            actual = len(self.cells)
            if actual == expected:
                print(f"  Final count: {actual} cells (predictable: 2^{self.year})")
            else:
                print(f"  Final count: {actual} cells (expected: {expected})")
        else:
            # Steady state: division, culling, methylation
            print(f"  Steady state phase")
            self.age_cells_steady_state()
            # Show it's random during steady state
            print(f"  Final count: {len(self.cells)} cells (random ~{self.target_population})")
            
        # Store current state after all operations
        self.cell_history[str(self.year)] = [cell.to_dict() for cell in self.cells]
        
        # Report statistics
        jsd_values = [cell.cell_JSD for cell in self.cells]
        mean_jsd = statistics.mean(jsd_values) if jsd_values else 0.0
        print(f"  Mean JSD: {mean_jsd:.4f}")
        
    def run_simulation(self, t_max: int = T_MAX) -> None:
        """
        Run complete simulation from year 0 to t_max.
        
        Args:
            t_max: Maximum simulation time in years
        """
        print("="*60)
        print("STEP1-PRIME SIMULATION")
        print("="*60)
        print(f"Parameters:")
        if self.rate is not None:
            print(f"  Methylation rate: {self.rate:.3%}")
        else:
            print(f"  Gene-specific rates: {len(self.gene_rate_groups)} groups")
            for i, (n_genes, rate) in enumerate(self.gene_rate_groups):
                print(f"    Group {i+1}: {n_genes} genes at {rate:.3%}")
        print(f"  CpG sites per cell: {self.n}")
        print(f"  Gene size: {self.gene_size}")
        print(f"  Growth phase: {self.growth_phase} years")
        print(f"  Target population: {self.target_population} (2^{self.growth_phase})")
        print(f"  Max years: {t_max}")
        print(f"  Random seed: {self.seed}")
        print("="*60)
        
        print(f"\nYear 0: Starting with 1 unmethylated cell")
        
        while self.year < t_max:
            self.simulate_year()
            
        print("\n" + "="*60)
        print("Simulation complete!")
        print(f"Final population: {len(self.cells)} cells")
        print("="*60)
        
    def save_history(self, filename: str = None, directory: str = "data") -> str:
        """
        Save simulation history to compressed JSON file using hierarchical structure.
        
        Args:
            filename: Output filename (ignored, kept for compatibility)
            directory: Base output directory
            
        Returns:
            Path to saved file
        """
        import hashlib
        
        # Generate hierarchical path
        # Level 1: Rate with 5 decimal places or gene rate groups
        if self.rate is not None:
            level1 = f"rate_{self.rate:.5f}"
            rate_str = f"{self.rate:.5f}"
        else:
            # For gene rate groups, create a descriptive name
            groups_str = "_".join([f"{n}x{rate:.5f}" for n, rate in self.gene_rate_groups])
            level1 = f"gene_rates_{groups_str}"[:50]  # Limit length
            rate_str = groups_str
        
        # Level 2: Parameters with hyphen separators
        seed_str = f"seed{self.seed}" if self.seed is not None else "noseed"
        params_str = f"grow{self.growth_phase}-sites{self.n}-years{self.year}-{seed_str}"
        
        # Add 4-char hash for uniqueness
        hash_input = f"{rate_str}-{self.growth_phase}-{self.n}-{self.year}-{seed_str}"
        hash_str = hashlib.md5(hash_input.encode()).hexdigest()[:4]
        level2 = f"{params_str}-{hash_str}"
        
        # Create full directory path
        dir_path = os.path.join(directory, level1, level2)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)
            print(f"Created directory: {dir_path}")
        
        # Fixed filename
        filepath = os.path.join(dir_path, "simulation.json.gz")
        print(f"\nSaving compressed history to {filepath}")
        
        save_start = statistics.mean([0])  # Just to have time module if needed
        import time
        save_start = time.time()
        
        # Use gzip compression and native JSON serialization
        with gzip.open(filepath, 'wt', encoding='utf-8', compresslevel=1) as f:
            json.dump(self.cell_history, f, separators=(',', ':'))
        
        save_time = time.time() - save_start
        file_size_mb = os.path.getsize(filepath) / (1024 * 1024)
        
        print(f"  Save time: {save_time:.2f} seconds")
        print(f"  File size: {file_size_mb:.2f} MB (compressed)")
        
        return filepath
    
    # ==================== Enhanced History Tracking Methods ====================
    
    def enable_history_tracking(self, start_year: int = 0, clear_history: bool = True, track_gene_jsd: bool = True) -> 'PetriDish':
        """
        Enable history tracking from a specific year.
        
        Args:
            start_year: The absolute year to start tracking from
            clear_history: Whether to clear existing history
            track_gene_jsd: Whether to also track gene JSD history
            
        Returns:
            Self for method chaining
        """
        self.track_cell_history = True
        self.track_gene_jsd = track_gene_jsd
        self.absolute_year = start_year
        if clear_history:
            self.cell_history = {}
            self.gene_jsd_history = {}
        self._record_history()  # Record initial state
        return self
    
    def disable_history_tracking(self) -> 'PetriDish':
        """Disable history tracking to save memory."""
        self.track_cell_history = False
        self.track_gene_jsd = False
        return self
    
    def _record_history(self, year: int = None) -> None:
        """
        Internal method to record current state in history.
        
        Args:
            year: Optional year to record at. If None, uses absolute_year
        """
        if year is None:
            year = self.absolute_year
        
        # Record cell history if enabled
        if self.track_cell_history:
            self.cell_history[str(year)] = [cell.to_dict() for cell in self.cells]
        
        # Record gene JSD if enabled
        if self.track_gene_jsd and self.calculate_jsds:
            self.gene_jsd_history[str(year)] = self.calculate_gene_jsd()
    
    def set_absolute_year(self, year: int) -> 'PetriDish':
        """
        Set the absolute year (for phase2 compatibility).
        
        Args:
            year: The absolute year to set
            
        Returns:
            Self for method chaining
        """
        self.absolute_year = year
        self.year = 0  # Reset relative year
        return self
    
    def increment_year(self, record_history: bool = True) -> 'PetriDish':
        """
        Advance both year counters and optionally record history.
        
        Args:
            record_history: Whether to record the state after incrementing
            
        Returns:
            Self for method chaining
        """
        self.year += 1
        self.absolute_year += 1
        if record_history:
            self._record_history()
        return self
    
    def calculate_gene_jsd(self) -> List[float]:
        """
        Calculate JSD for each gene across all cells.
        
        Returns:
            List of JSD values, one per gene
        """
        if not self.cells:
            return [0.0] * self.n_genes
        
        n_cells = len(self.cells)
        
        if HAS_NUMPY:
            # Numpy implementation (faster)
            methylation_matrix = np.array([cell.cpg_sites for cell in self.cells])
            genes_matrix = methylation_matrix.reshape(n_cells, self.n_genes, self.gene_size)
            methylation_counts = genes_matrix.sum(axis=2)
            
            gene_jsds = []
            for gene_idx in range(self.n_genes):
                gene_counts = methylation_counts[:, gene_idx]
                distribution = np.zeros(self.gene_size + 1)
                for count in range(self.gene_size + 1):
                    distribution[count] = np.sum(gene_counts == count) / n_cells
                jsd = JS_div(distribution.tolist(), self.BASELINE_GENE_DISTRIBUTION)
                gene_jsds.append(jsd)
        else:
            # Pure Python implementation (slower but no dependencies)
            gene_jsds = []
            for gene_idx in range(self.n_genes):
                # Count methylation levels for this gene across all cells
                distribution = [0] * (self.gene_size + 1)
                
                for cell in self.cells:
                    # Get sites for this gene
                    start_idx = gene_idx * self.gene_size
                    end_idx = start_idx + self.gene_size
                    gene_sites = cell.cpg_sites[start_idx:end_idx]
                    
                    # Count methylated sites
                    methylated_count = sum(gene_sites)
                    distribution[methylated_count] += 1
                
                # Normalize to get probabilities
                distribution = [count / n_cells for count in distribution]
                
                # Calculate JSD with baseline
                jsd = JS_div(distribution, self.BASELINE_GENE_DISTRIBUTION)
                gene_jsds.append(jsd)
        
        return gene_jsds
    
    def grow_with_homeostasis(self, years: int, growth_phase: int = None,
                             verbose: bool = False, record_history: bool = True) -> 'PetriDish':
        """
        Phase2-style growth with explicit growth phase control.
        Used by phase2's grow_petri_for_years().
        
        Args:
            years: Number of years to simulate
            growth_phase: Years of exponential growth before homeostasis (None uses self.growth_phase)
            verbose: Whether to print progress
            record_history: Whether to track history during growth
            
        Returns:
            Self for method chaining
        """
        if growth_phase is None:
            growth_phase = self.growth_phase
            
        for year_idx in range(years):
            current_year = year_idx + 1
            initial_count = len(self.cells)
            
            if current_year <= growth_phase:
                # Growth phase: divide and methylate
                self.divide_cells()
                self.methylate_cells()
                phase = "growth"
            else:
                # Homeostasis phase: divide, cull, and methylate
                self.divide_cells()
                self.random_cull_cells()
                self.methylate_cells()
                phase = "homeostasis"
            
            final_count = len(self.cells)
            
            if verbose:
                abs_year = self.absolute_year + 1  # Show what year we're moving to
                print(f"      Year {current_year} ({phase}): {initial_count} → {final_count} cells")
            
            # Increment year and record history
            if record_history:
                self.increment_year(record_history=True)
            else:
                self.year += 1
                self.absolute_year += 1
        
        return self


class PetriDishPlotter:
    """
    Handles all plotting for PetriDish objects.
    Reuses the logic from plot_history.py but in an object-oriented way.
    """
    
    def __init__(self, petri: PetriDish):
        """
        Initialize with a PetriDish that has history.
        
        Args:
            petri: PetriDish object with history to plot
            
        Raises:
            ValueError: If PetriDish has no history
        """
        self.petri = petri
        self.stats = None  # Cache calculated statistics
        
        # Validate that petri has history
        if not hasattr(petri, 'cell_history') or not petri.cell_history:
            raise ValueError("PetriDish has no cell history to plot")
    
    def calculate_statistics(self) -> Dict:
        """Calculate statistics from petri's history."""
        stats = {
            'years': [],
            'population_size': [],
            'cell_jsd': {
                'mean': [], 'median': [], 'p5': [], 'p25': [],
                'p75': [], 'p95': [], 'min': [], 'max': []
            },
            'methylation': {
                'mean': [], 'median': [], 'p5': [], 'p25': [],
                'p75': [], 'p95': []
            }
        }
        
        # Sort years numerically
        sorted_years = sorted([int(year) for year in self.petri.cell_history.keys()])
        
        for year in sorted_years:
            year_data = self.petri.cell_history[str(year)]
            
            # Skip if no cells (edge case)
            if not year_data:
                continue
            
            # Extract values
            jsd_values = [cell['cell_jsd'] for cell in year_data]
            meth_values = [cell['methylation_proportion'] * 100 for cell in year_data]
            
            stats['years'].append(year)
            stats['population_size'].append(len(year_data))
            
            # JSD statistics
            if jsd_values:  # Ensure we have data
                stats['cell_jsd']['mean'].append(statistics.mean(jsd_values))
                stats['cell_jsd']['median'].append(statistics.median(jsd_values))
                
                # For single values, all percentiles are the same
                if len(jsd_values) == 1:
                    val = jsd_values[0]
                    stats['cell_jsd']['p5'].append(val)
                    stats['cell_jsd']['p25'].append(val)
                    stats['cell_jsd']['p75'].append(val)
                    stats['cell_jsd']['p95'].append(val)
                    stats['cell_jsd']['min'].append(val)
                    stats['cell_jsd']['max'].append(val)
                else:
                    # Use statistics.quantiles for percentiles
                    import numpy as np
                    stats['cell_jsd']['p5'].append(np.percentile(jsd_values, 5))
                    stats['cell_jsd']['p25'].append(np.percentile(jsd_values, 25))
                    stats['cell_jsd']['p75'].append(np.percentile(jsd_values, 75))
                    stats['cell_jsd']['p95'].append(np.percentile(jsd_values, 95))
                    stats['cell_jsd']['min'].append(min(jsd_values))
                    stats['cell_jsd']['max'].append(max(jsd_values))
            
            # Methylation statistics
            if meth_values:
                stats['methylation']['mean'].append(statistics.mean(meth_values))
                stats['methylation']['median'].append(statistics.median(meth_values))
                
                if len(meth_values) == 1:
                    val = meth_values[0]
                    stats['methylation']['p5'].append(val)
                    stats['methylation']['p25'].append(val)
                    stats['methylation']['p75'].append(val)
                    stats['methylation']['p95'].append(val)
                else:
                    import numpy as np
                    stats['methylation']['p5'].append(np.percentile(meth_values, 5))
                    stats['methylation']['p25'].append(np.percentile(meth_values, 25))
                    stats['methylation']['p75'].append(np.percentile(meth_values, 75))
                    stats['methylation']['p95'].append(np.percentile(meth_values, 95))
        
        self.stats = stats
        return stats
    
    def detect_growth_phase(self) -> int:
        """Auto-detect growth phase from population dynamics."""
        if not self.stats:
            self.calculate_statistics()
        
        pop_sizes = self.stats['population_size']
        
        for i in range(1, len(pop_sizes)):
            if i > 1:
                expected = 2 ** i
                actual = pop_sizes[i]
                if actual != expected:
                    return i - 1
        
        # Fallback to petri's growth_phase or default
        return getattr(self.petri, 'growth_phase', 13)
    
    def plot_jsd(self, title: str = None, output_path: str = None,
                 width: int = 1200, height: int = 500):
        """
        Create JSD plot with secondary y-axis for cell count.
        
        Args:
            title: Optional title for the plot
            output_path: Optional path to save the plot
            width: Plot width in pixels
            height: Plot height in pixels
            
        Returns:
            Plotly figure object
        """
        # Import plotly here to avoid dependency issues
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not self.stats:
            self.calculate_statistics()
        
        # Create figure with secondary y-axis
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        years = self.stats['years']
        growth_phase = self.detect_growth_phase()
        
        # Add vertical line for growth phase
        if growth_phase > 0 and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text=f"End of growth phase",
                annotation_position="top"
            )
        
        # Add 5-95 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p95'] + self.stats['cell_jsd']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add 25-75 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p75'] + self.stats['cell_jsd']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add mean line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['cell_jsd']['mean'],
                mode='lines',
                name='Mean JSD',
                line=dict(color='rgb(99, 110, 250)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean JSD: %{y:.4f}<br>Population: %{customdata}<extra></extra>',
                customdata=self.stats['population_size']
            ),
            secondary_y=False
        )
        
        # Add cell count on secondary y-axis
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>'
            ),
            secondary_y=True
        )
        
        # Update layout
        if title is None:
            title = f"JSD Score vs Time"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            height=height,
            margin=dict(t=120),
            hovermode='x unified',
            showlegend=True,
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            ),
            template='plotly_white'
        )
        
        # Set axis titles
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="JSD Score", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
        
        # Add annotation
        final_idx = -1
        final_year = years[final_idx]
        final_pop = self.stats['population_size'][final_idx]
        
        annotation_text = (
            f"<b>Final Statistics (Year {final_year}):</b><br>"
            f"Population: {final_pop} cells<br>"
            f"JSD Mean: {self.stats['cell_jsd']['mean'][final_idx]:.4f}<br>"
            f"JSD 25-75%: [{self.stats['cell_jsd']['p25'][final_idx]:.4f}, {self.stats['cell_jsd']['p75'][final_idx]:.4f}]<br>"
            f"JSD 5-95%: [{self.stats['cell_jsd']['p5'][final_idx]:.4f}, {self.stats['cell_jsd']['p95'][final_idx]:.4f}]"
        )
        
        if growth_phase > 0:
            annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
        
        fig.add_annotation(
            text=annotation_text,
            xref="paper", yref="paper",
            x=0.5, y=1.15,
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1,
            font=dict(size=10),
            align="center",
            xanchor="center",
            yanchor="top"
        )
        
        if output_path:
            fig.write_image(output_path, width=width, height=height, scale=2)
            
        return fig
    
    def plot_methylation(self, title: str = None, output_path: str = None,
                        width: int = 1200, height: int = 500):
        """Create methylation plot with secondary y-axis for cell count."""
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not self.stats:
            self.calculate_statistics()
        
        # Create figure with secondary y-axis
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        years = self.stats['years']
        growth_phase = self.detect_growth_phase()
        
        # Add vertical line for growth phase
        if growth_phase > 0 and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text=f"End of growth phase",
                annotation_position="top"
            )
        
        # Add 5-95 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p95'] + self.stats['methylation']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add 25-75 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p75'] + self.stats['methylation']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add mean line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['methylation']['mean'],
                mode='lines',
                name='Mean Methylation',
                line=dict(color='rgb(239, 85, 59)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Methylation: %{y:.2f}%<br>Population: %{customdata}<extra></extra>',
                customdata=self.stats['population_size']
            ),
            secondary_y=False
        )
        
        # Add cell count on secondary y-axis
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>'
            ),
            secondary_y=True
        )
        
        # Update layout
        if title is None:
            title = f"Methylation Proportion vs Time"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            height=height,
            margin=dict(t=120),
            hovermode='x unified',
            showlegend=True,
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            ),
            template='plotly_white'
        )
        
        # Set axis titles
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Methylation (%)", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
        
        # Add annotation
        final_idx = -1
        final_year = years[final_idx]
        final_pop = self.stats['population_size'][final_idx]
        
        annotation_text = (
            f"<b>Final Statistics (Year {final_year}):</b><br>"
            f"Population: {final_pop} cells<br>"
            f"Methylation Mean: {self.stats['methylation']['mean'][final_idx]:.2f}%<br>"
            f"Methylation 25-75%: [{self.stats['methylation']['p25'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p75'][final_idx]:.2f}%]<br>"
            f"Methylation 5-95%: [{self.stats['methylation']['p5'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p95'][final_idx]:.2f}%]"
        )
        
        if growth_phase > 0:
            annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
        
        fig.add_annotation(
            text=annotation_text,
            xref="paper", yref="paper",
            x=0.5, y=1.15,
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1,
            font=dict(size=10),
            align="center",
            xanchor="center",
            yanchor="top"
        )
        
        if output_path:
            fig.write_image(output_path, width=width, height=height, scale=2)
            
        return fig
    
    def plot_combined(self, title: str = None, output_path: str = None,
                     width: int = 1200, height: int = 800):
        """Create combined JSD and methylation plot."""
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not self.stats:
            self.calculate_statistics()
        
        # Create subplots with secondary y-axes
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('JSD Score vs Time', 'Methylation Proportion vs Time'),
            vertical_spacing=0.12,
            row_heights=[0.5, 0.5],
            specs=[[{"secondary_y": True}], [{"secondary_y": True}]]
        )
        
        years = self.stats['years']
        growth_phase = self.detect_growth_phase()
        
        # Add vertical lines
        if growth_phase > 0 and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text=f"End of growth phase",
                annotation_position="top",
                row=1, col=1
            )
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                row=2, col=1
            )
        
        # ===== JSD Plot (Row 1) =====
        # Add percentile bands
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p95'] + self.stats['cell_jsd']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            row=1, col=1, secondary_y=False
        )
        
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p75'] + self.stats['cell_jsd']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            row=1, col=1, secondary_y=False
        )
        
        # Add mean JSD line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['cell_jsd']['mean'],
                mode='lines',
                name='Mean JSD',
                line=dict(color='rgb(99, 110, 250)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean JSD: %{y:.4f}<extra></extra>'
            ),
            row=1, col=1, secondary_y=False
        )
        
        # Cell count for JSD plot
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
                showlegend=True
            ),
            row=1, col=1, secondary_y=True
        )
        
        # ===== Methylation Plot (Row 2) =====
        # Add percentile bands
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p95'] + self.stats['methylation']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=False,
                hoverinfo='skip'
            ),
            row=2, col=1, secondary_y=False
        )
        
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p75'] + self.stats['methylation']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=False,
                hoverinfo='skip'
            ),
            row=2, col=1, secondary_y=False
        )
        
        # Add mean methylation line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['methylation']['mean'],
                mode='lines',
                name='Mean Methylation',
                line=dict(color='rgb(239, 85, 59)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Methylation: %{y:.2f}%<extra></extra>'
            ),
            row=2, col=1, secondary_y=False
        )
        
        # Cell count for methylation plot
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
                showlegend=False
            ),
            row=2, col=1, secondary_y=True
        )
        
        # Update layout
        if title is None:
            title = f"Simulation Results"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            height=height,
            margin=dict(t=120),
            hovermode='x unified',
            showlegend=True,
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            ),
            template='plotly_white'
        )
        
        # Update axes
        fig.update_xaxes(title_text="", row=1, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_xaxes(title_text="Age (years)", row=2, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="JSD Score", row=1, col=1, secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", row=1, col=1, secondary_y=True, showgrid=False)
        fig.update_yaxes(title_text="Methylation (%)", row=2, col=1, secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", row=2, col=1, secondary_y=True, showgrid=False)
        
        # Add annotation
        final_idx = -1
        final_year = years[final_idx]
        final_pop = self.stats['population_size'][final_idx]
        
        annotation_text = (
            f"<b>Final Statistics (Year {final_year}):</b><br>"
            f"Population: {final_pop} cells<br>"
            f"JSD Mean: {self.stats['cell_jsd']['mean'][final_idx]:.4f}<br>"
            f"JSD 25-75%: [{self.stats['cell_jsd']['p25'][final_idx]:.4f}, {self.stats['cell_jsd']['p75'][final_idx]:.4f}]<br>"
            f"JSD 5-95%: [{self.stats['cell_jsd']['p5'][final_idx]:.4f}, {self.stats['cell_jsd']['p95'][final_idx]:.4f}]<br>"
            f"Methylation Mean: {self.stats['methylation']['mean'][final_idx]:.2f}%<br>"
            f"Methylation 25-75%: [{self.stats['methylation']['p25'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p75'][final_idx]:.2f}%]<br>"
            f"Methylation 5-95%: [{self.stats['methylation']['p5'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p95'][final_idx]:.2f}%]"
        )
        
        if growth_phase > 0:
            annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
        
        fig.add_annotation(
            text=annotation_text,
            xref="paper", yref="paper",
            x=0.5, y=1.15,
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1,
            font=dict(size=10),
            align="center",
            xanchor="center",
            yanchor="top"
        )
        
        if output_path:
            fig.write_image(output_path, width=width, height=height, scale=2)
            
        return fig
    
    def plot_all(self, base_path: str, title_prefix: str = None):
        """
        Generate all three plot types.
        
        Args:
            base_path: Base path for output files (without extension)
            title_prefix: Optional prefix for plot titles
            
        Returns:
            Self for method chaining
        """
        title = title_prefix or "PetriDish Simulation"
        
        self.plot_jsd(title, f"{base_path}_jsd.png")
        self.plot_methylation(title, f"{base_path}_methylation.png")
        self.plot_combined(title, f"{base_path}_combined.png")
        
        return self