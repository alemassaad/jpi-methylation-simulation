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


def rate_to_gene_rate_groups(rate: float, n: int, gene_size: int) -> List[Tuple[int, float]]:
    """
    Convert uniform methylation rate to gene_rate_groups format.
    
    Args:
        rate: Uniform methylation rate
        n: Number of CpG sites
        gene_size: Sites per gene
    
    Returns:
        List with single tuple: [(n_genes, rate)]
    """
    n_genes = n // gene_size
    return [(n_genes, rate)]


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
        
        # Convert rate to gene_rate_groups if provided
        if rate is not None:
            if gene_rate_groups is not None:
                raise ValueError(
                    "Cannot specify both 'rate' and 'gene_rate_groups'. "
                    "Use 'rate' for uniform methylation or 'gene_rate_groups' for gene-specific rates."
                )
            gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
        
        # Must have gene_rate_groups by now
        if gene_rate_groups is None:
            raise ValueError(
                "Must specify either 'rate' (uniform) or 'gene_rate_groups' (gene-specific)"
            )
        
        self.n = n
        self.gene_size = gene_size
        self.gene_rate_groups = gene_rate_groups  # Only store this
        
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
        """Build array of per-site methylation rates from gene_rate_groups."""
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
    
    @property
    def rate(self) -> Optional[float]:
        """
        Get uniform rate if applicable (backward compatibility).
        Returns None if gene-specific rates are used.
        """
        # Check if all rates are the same
        rates = set(rate for _, rate in self.gene_rate_groups)
        if len(rates) == 1:
            return rates.pop()
        return None
    
    def to_dict(self) -> Dict[str, Any]:
        """
        Convert cell state to dictionary for serialization.
        
        New lean format: Only stores essential cell-specific data.
        Shared parameters (rate, gene_size) are stored once at the PetriDish level.
        """
        return {
            'methylated': self.cpg_sites[:],  # The methylation state (0s and 1s)
            'cell_JSD': self.cell_JSD,        # The JSD value
            'age': self.age                    # Cell age in years
            # Note: rate, gene_size are stored in parameters
            # methylation_proportion and distribution can be calculated from methylated array
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any], rate: float = None, 
                  gene_rate_groups: List[Tuple[int, float]] = None,
                  gene_size: int = GENE_SIZE) -> 'Cell':
        """
        Create a Cell from dictionary data (new lean format only).
        
        Args:
            data: Dictionary with 'methylated' and 'cell_JSD' keys
            rate: Uniform methylation rate (will be converted to gene_rate_groups)
            gene_rate_groups: Gene-specific rates
            gene_size: Sites per gene
            
        Returns:
            Reconstructed Cell object
        """
        n = len(data['methylated'])
        
        # Let __init__ handle the conversion
        cell = cls(n=n, rate=rate, gene_rate_groups=gene_rate_groups, gene_size=gene_size)
        cell.cpg_sites = data['methylated'][:]
        cell.cell_JSD = data.get('cell_JSD', 0.0)
        cell.age = data.get('age', 0)  # Restore age (default to 0 for old data)
        
        return cell


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
        validate_cell_consistency(): Ensure all cells have same rate configuration
        validate_cells_compatible(): Static method to pre-check cell compatibility
    
    Important:
        All cells in a PetriDish must have identical gene_rate_groups configuration.
        This is validated automatically when creating a PetriDish with existing cells.
    """
    
    def __init__(self, rate: float = None, gene_rate_groups: List[Tuple[int, float]] = None,
                 n: int = N, gene_size: int = GENE_SIZE, 
                 seed: int = None, growth_phase: int = DEFAULT_GROWTH_PHASE, 
                 cells: List['Cell'] = None, calculate_cell_jsds: bool = True) -> None:
        """
        Initialize petri dish with a single unmethylated cell or provided cells.
        
        Args:
            rate: Uniform methylation rate per site per year (will be converted to gene_rate_groups)
            gene_rate_groups: Gene-specific rates as [(n_genes, rate), ...]
            n: Number of CpG sites per cell
            gene_size: Number of sites per gene
            seed: Random seed for reproducibility
            growth_phase: Duration of growth phase in years (target = 2^growth_phase cells)
            cells: Optional list of cells to start with (for phase2 compatibility).
                   All cells must have identical gene_rate_groups configuration.
            calculate_cell_jsds: Whether to calculate cell JSDs (for performance)
            
        Raises:
            ValueError: If provided cells have different gene_rate_groups configurations
        """
        # Convert rate to gene_rate_groups if provided
        if rate is not None:
            if gene_rate_groups is not None:
                raise ValueError(
                    "Cannot specify both 'rate' and 'gene_rate_groups'. "
                    "Use 'rate' for uniform methylation or 'gene_rate_groups' for gene-specific rates."
                )
            gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
        
        # Use default if neither specified
        if gene_rate_groups is None:
            gene_rate_groups = rate_to_gene_rate_groups(RATE, n, gene_size)
        
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
        
        self.gene_rate_groups = gene_rate_groups  # Only store this
        self.n = n
        self.gene_size = gene_size
        self.n_genes = n // gene_size
        self.seed = seed
        self.growth_phase = growth_phase
        self.target_population = 2 ** growth_phase  # Calculate from growth_phase
        self.calculate_cell_jsds = calculate_cell_jsds
        
        # Initialize cells - either provided or single unmethylated cell
        if cells is not None:
            self.cells = cells
            # Validate all cells have same rate configuration
            self.validate_cell_consistency()
        else:
            # Create initial cell - single path now
            self.cells = [Cell(n=n, gene_rate_groups=gene_rate_groups, gene_size=gene_size)]
        
        self.year = 0
        
        # Cell history tracking (renamed for clarity)
        self.track_cell_history = True  # Renamed from history_enabled
        
        # Gene JSD tracking
        self.track_gene_jsd = False
        self.gene_jsd_history = {}
        self.mean_gene_jsd_history = {}    # Track mean gene JSD over time
        self.median_gene_jsd_history = {}  # Track median gene JSD over time
        # Baseline must match gene_size + 1 bins
        self.BASELINE_GENE_DISTRIBUTION = [1.0] + [0.0] * self.gene_size
        
        # Store initial state if we have cells (renamed for clarity)
        if self.cells:
            self.cell_history = {'0': [cell.to_dict() for cell in self.cells]}
            # Initialize gene JSD at year 0 if tracking
            if self.track_gene_jsd and self.calculate_cell_jsds:
                self.gene_jsd_history['0'] = [0.0] * self.n_genes  # All zeros initially
                self.mean_gene_jsd_history['0'] = 0.0
                self.median_gene_jsd_history['0'] = 0.0
        else:
            self.cell_history = {}
        
    @property
    def rate(self) -> Optional[float]:
        """
        Get uniform rate if applicable (backward compatibility).
        Returns None if gene-specific rates are used.
        """
        # Check if all rates are the same
        rates = set(rate for _, rate in self.gene_rate_groups)
        if len(rates) == 1:
            return rates.pop()
        return None
    
    @staticmethod
    def validate_cells_compatible(cells: List['Cell']) -> None:
        """
        Validate that a list of cells have identical gene_rate_groups.
        
        This static method can be called before creating a PetriDish to ensure
        cell compatibility. Useful for pre-validation in phase2 or when combining
        cells from different sources.
        
        Args:
            cells: List of Cell objects to validate
            
        Raises:
            ValueError: If cells have different gene_rate_groups configurations
            
        Example:
            >>> cells = [cell1, cell2, cell3]
            >>> PetriDish.validate_cells_compatible(cells)  # Raises if incompatible
            >>> petri = PetriDish(cells=cells, ...)  # Safe to create
        """
        if not cells:
            return
        
        first_groups = cells[0].gene_rate_groups
        
        for i, cell in enumerate(cells[1:], 1):
            if cell.gene_rate_groups != first_groups:
                raise ValueError(
                    f"Rate configuration mismatch: All cells must have "
                    f"identical gene_rate_groups.\n"
                    f"  Cell 0: gene_rate_groups={first_groups}\n"
                    f"  Cell {i}: gene_rate_groups={cell.gene_rate_groups}\n"
                    f"Ensure all cells are created with the same rate configuration."
                )
    
    def validate_cell_consistency(self) -> None:
        """
        Validate that all cells have identical gene_rate_groups configuration.
        
        This ensures rate consistency across the population, which is required
        for accurate simulation and analysis. Called automatically during
        PetriDish initialization when cells are provided.
        
        Raises:
            ValueError: If cells have different gene_rate_groups configurations
        """
        if not self.cells:
            return
        
        first_groups = self.cells[0].gene_rate_groups
        
        for i, cell in enumerate(self.cells[1:], 1):
            if cell.gene_rate_groups != first_groups:
                raise ValueError(
                    f"Rate configuration mismatch: All cells in a PetriDish must have "
                    f"identical gene_rate_groups.\n"
                    f"  Cell 0: gene_rate_groups={first_groups}\n"
                    f"  Cell {i}: gene_rate_groups={cell.gene_rate_groups}\n"
                    f"Ensure all cells are created with the same rate configuration."
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
        
        # Update metadata if it exists
        if hasattr(self, 'metadata'):
            self.metadata['num_cells'] = len(self.cells)
        
    def methylate_cells(self) -> None:
        """
        Apply stochastic methylation to all cells.
        Each cell independently accumulates methylation.
        """
        for cell in self.cells:
            cell.methylate()
            # Only calculate JSD if enabled
            if not self.calculate_cell_jsds:
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
        
        # Update metadata if it exists
        if hasattr(self, 'metadata'):
            self.metadata['num_cells'] = len(self.cells)
        
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
        if self.track_cell_history:
            self._record_history(self.year)
        else:
            # Always store in cell_history for backward compatibility
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
        print("PHASE 1 SIMULATION")
        print("="*60)
        print(f"Parameters:")
        # Always show as gene groups (even if uniform)
        print(f"  Gene-specific rates: {len(self.gene_rate_groups)} group(s)")
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
        
    def save_history(self, filename: str = None, directory: str = "data", compress: bool = True) -> str:
        """
        Save simulation history to JSON file (compressed or uncompressed) using hierarchical structure.
        
        Args:
            filename: Output filename (ignored, kept for compatibility)
            directory: Base output directory
            compress: If True, save as .json.gz; if False, save as .json
            
        Returns:
            Path to saved file
        """
        from datetime import datetime
        
        # Generate hierarchical path
        # Always use gene_rates format for consistency
        groups_str = "_".join([f"{n}x{rate:.5f}" for n, rate in self.gene_rate_groups])
        level1 = f"gene_rates_{groups_str}"[:50]  # Limit length
        
        # Level 2: Parameters with hyphen separators
        seed_str = f"seed{self.seed}" if self.seed is not None else "noseed"
        # Use size instead of grow for clarity (size = 2^growth_phase)
        population_size = 2 ** self.growth_phase
        params_str = f"size{population_size}-sites{self.n}-genesize{self.gene_size}-years{self.year}-{seed_str}"
        
        # Add timestamp for uniqueness (YYYYMMDD-HHMMSS format)
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        level2 = f"{params_str}-{timestamp}"
        
        # Create full directory path
        dir_path = os.path.join(directory, level1, level2)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path, exist_ok=True)
            print(f"Created directory: {dir_path}")
        
        # Fixed filename with appropriate extension
        extension = ".json.gz" if compress else ".json"
        filepath = os.path.join(dir_path, f"simulation{extension}")
        format_type = "compressed" if compress else "uncompressed"
        print(f"\nSaving {format_type} history to {filepath}")
        
        save_start = statistics.mean([0])  # Just to have time module if needed
        import time
        save_start = time.time()
        
        # Prepare data to save with new format
        save_data = {
            'parameters': {
                'gene_rate_groups': self.gene_rate_groups,  # Only store this
                'n': self.n,
                'gene_size': self.gene_size,
                'growth_phase': self.growth_phase,
                'years': self.year,
                'seed': self.seed,
                'track_cell_history': self.track_cell_history,
                'track_gene_jsd': self.track_gene_jsd
            },
            'history': {}
        }
        
        # Add cell history
        if hasattr(self, 'cell_history') and self.cell_history:
            for year_str, cells in self.cell_history.items():
                save_data['history'][year_str] = {
                    'cells': cells  # Already in dict format from to_dict()
                }
        
        # Add gene_jsd_history if it exists
        if hasattr(self, 'gene_jsd_history') and self.gene_jsd_history:
            for year_str, gene_jsds in self.gene_jsd_history.items():
                if year_str not in save_data['history']:
                    save_data['history'][year_str] = {}
                save_data['history'][year_str]['gene_jsd'] = gene_jsds
        
        # Add mean and median gene JSD histories
        if hasattr(self, 'mean_gene_jsd_history') and self.mean_gene_jsd_history:
            for year_str, mean_jsd in self.mean_gene_jsd_history.items():
                if year_str not in save_data['history']:
                    save_data['history'][year_str] = {}
                save_data['history'][year_str]['mean_gene_jsd'] = mean_jsd
        
        if hasattr(self, 'median_gene_jsd_history') and self.median_gene_jsd_history:
            for year_str, median_jsd in self.median_gene_jsd_history.items():
                if year_str not in save_data['history']:
                    save_data['history'][year_str] = {}
                save_data['history'][year_str]['median_gene_jsd'] = median_jsd
        
        # Save with or without compression
        if compress:
            # Use gzip compression and compact JSON
            with gzip.open(filepath, 'wt', encoding='utf-8', compresslevel=1) as f:
                json.dump(save_data, f, separators=(',', ':'))
        else:
            # Save uncompressed with indentation for readability
            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(save_data, f, indent=2)
        
        save_time = time.time() - save_start
        file_size_mb = os.path.getsize(filepath) / (1024 * 1024)
        
        print(f"  Save time: {save_time:.2f} seconds")
        format_note = " (compressed)" if compress else " (uncompressed)"
        print(f"  File size: {file_size_mb:.2f} MB{format_note}")
        
        return filepath
    
    # ==================== Enhanced History Tracking Methods ====================
    
    def enable_history_tracking(self, clear_history: bool = True, track_gene_jsd: bool = True) -> 'PetriDish':
        """
        Enable history tracking.
        
        Args:
            clear_history: Whether to clear existing history
            track_gene_jsd: Whether to also track gene JSD history
            
        Returns:
            Self for method chaining
        """
        self.track_cell_history = True
        self.track_gene_jsd = track_gene_jsd
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
            year: Optional year to record at. If None, uses self.year
        """
        if year is None:
            year = self.year
        
        # Record cell history if enabled
        if self.track_cell_history:
            self.cell_history[str(year)] = [cell.to_dict() for cell in self.cells]
        
        # Record gene JSD if enabled
        if self.track_gene_jsd and self.calculate_cell_jsds:
            gene_jsds = self.calculate_gene_jsd()
            self.gene_jsd_history[str(year)] = gene_jsds
            
            # Calculate and store summary statistics
            if gene_jsds:
                self.mean_gene_jsd_history[str(year)] = sum(gene_jsds) / len(gene_jsds)
                self.median_gene_jsd_history[str(year)] = statistics.median(gene_jsds)
            else:
                self.mean_gene_jsd_history[str(year)] = 0.0
                self.median_gene_jsd_history[str(year)] = 0.0
    
    
    def increment_year(self, record_history: bool = True) -> 'PetriDish':
        """
        Advance year counter and optionally record history.
        
        Args:
            record_history: Whether to record the state after incrementing
            
        Returns:
            Self for method chaining
        """
        self.year += 1
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
    
    def calculate_gene_jsds(self) -> List[float]:
        """
        Calculate Jensen-Shannon Divergence for each gene across all cells.
        This is identical to calculate_gene_jsd but with a clearer name.
        
        Returns:
            List of JSD values, one per gene
        """
        return self.calculate_gene_jsd()
    
    def calculate_gene_mean_methylation(self) -> List[float]:
        """
        Calculate mean methylation proportion for each gene across all cells.
        
        Returns:
            List of mean methylation proportions (0.0 to 1.0), one per gene.
            0.0 = completely unmethylated, 1.0 = fully methylated
        """
        if not self.cells:
            return []
        
        gene_means = []
        
        for gene_idx in range(self.n_genes):
            start = gene_idx * self.gene_size
            end = start + self.gene_size
            
            # Calculate mean methylation for this gene across all cells
            total_methylation = 0
            for cell in self.cells:
                # Sum methylation for this gene in this cell
                gene_methylation = sum(cell.cpg_sites[start:end])
                total_methylation += gene_methylation
            
            # Calculate mean count across all cells
            mean_count = total_methylation / len(self.cells)
            # Convert to proportion (0.0 to 1.0)
            mean_proportion = mean_count / self.gene_size
            gene_means.append(mean_proportion)
        
        return gene_means
    
    @property
    def mean_gene_jsd(self) -> float:
        """Current mean JSD across all genes."""
        if not self.cells:
            return 0.0
        gene_jsds = self.calculate_gene_jsd()
        return sum(gene_jsds) / len(gene_jsds) if gene_jsds else 0.0
    
    @property
    def median_gene_jsd(self) -> float:
        """Current median JSD across all genes."""
        if not self.cells:
            return 0.0
        gene_jsds = self.calculate_gene_jsd()
        return statistics.median(gene_jsds) if gene_jsds else 0.0
    
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
                print(f"      Year {current_year} ({phase}): {initial_count} → {final_count} cells")
            
            # Increment year and record history
            if record_history:
                self.increment_year(record_history=True)
            else:
                self.year += 1
        
        return self
    
    # ==================== Professional Pipeline Methods ====================
    
    @classmethod
    def from_snapshot_cell(cls, cell: 'Cell', growth_phase: int = 7, 
                          calculate_cell_jsds: bool = True,
                          metadata: Dict = None) -> 'PetriDish':
        """
        Create a PetriDish from a single snapshot cell.
        This is the professional way to create individuals in phase2 pipeline.
        
        Args:
            cell: The cell to use as the founding population
            growth_phase: The growth phase for this PetriDish
            calculate_cell_jsds: Whether to calculate cell JSDs
            metadata: Optional metadata to attach to the PetriDish
            
        Returns:
            A properly initialized PetriDish with the given cell
        """
        # Extract rate configuration from the cell
        if cell.rate is not None:
            # Uniform rate
            petri = cls(
                n=cell.n,
                rate=cell.rate,
                growth_phase=growth_phase,
                seed=None,  # Don't set seed for individual dishes
                calculate_cell_jsds=calculate_cell_jsds
            )
        else:
            # Gene-specific rates
            petri = cls(
                n=cell.n,
                gene_rate_groups=cell.gene_rate_groups,
                growth_phase=growth_phase,
                seed=None,
                calculate_cell_jsds=calculate_cell_jsds
            )
        
        # Replace the initial unmethylated cell with our snapshot cell
        petri.cells = [cell]
        
        # Set metadata if provided
        if metadata:
            petri.metadata = metadata.copy()
        else:
            petri.metadata = {}
        
        # Ensure critical metadata is set
        petri.update_metadata({
            'num_cells': 1,
            'creation_method': 'from_snapshot_cell'
        })
        
        return petri
    
    @classmethod
    def from_cells(cls, cells, growth_phase: int = 7, 
                   calculate_cell_jsds: bool = True,
                   metadata: Dict = None) -> 'PetriDish':
        """
        Create a PetriDish from one or more cells without intermediate steps.
        This is cleaner than from_snapshot_cell as it doesn't create/replace cells.
        
        Args:
            cells: Single cell or list of cells to use as founding population
            growth_phase: The growth phase for this PetriDish
            calculate_cell_jsds: Whether to calculate cell JSDs
            metadata: Optional metadata to attach to the PetriDish
            
        Returns:
            A properly initialized PetriDish with the given cells
        """
        # Handle single cell or list
        cells_list = cells if isinstance(cells, list) else [cells]
        first_cell = cells_list[0]
        
        # Extract rate configuration from first cell
        if first_cell.rate is not None:
            # Uniform rate
            petri = cls(
                cells=cells_list,  # Pass cells directly - no replacement needed!
                n=first_cell.n,
                rate=first_cell.rate,
                gene_size=first_cell.gene_size,
                growth_phase=growth_phase,
                seed=None,
                calculate_cell_jsds=calculate_cell_jsds
            )
        else:
            # Gene-specific rates
            petri = cls(
                cells=cells_list,  # Pass cells directly - no replacement needed!
                n=first_cell.n,
                gene_rate_groups=first_cell.gene_rate_groups,
                gene_size=first_cell.gene_size,
                growth_phase=growth_phase,
                seed=None,
                calculate_cell_jsds=calculate_cell_jsds
            )
        
        # Set metadata if provided
        if metadata:
            petri.metadata = metadata.copy()
        else:
            petri.metadata = {}
        
        # Ensure critical metadata is set
        petri.update_metadata({
            'num_cells': len(cells_list),
            'creation_method': 'from_cells'
        })
        
        return petri
    
    def update_metadata(self, updates: Dict) -> 'PetriDish':
        """
        Update the PetriDish metadata with new values.
        
        Args:
            updates: Dictionary of metadata updates
            
        Returns:
            Self for method chaining
        """
        if not hasattr(self, 'metadata'):
            self.metadata = {}
        self.metadata.update(updates)
        
        # Always keep num_cells in sync
        self.metadata['num_cells'] = len(self.cells)
        
        return self
    
    def _calculate_mean_methylation_proportion(self) -> float:
        """
        Calculate the mean methylation proportion across all cells.
        
        Returns:
            Mean proportion of methylated sites (0.0 to 1.0)
        """
        if not self.cells:
            return 0.0
        
        total_methylation = sum(sum(cell.cpg_sites) for cell in self.cells)
        total_sites = len(self.cells) * self.n
        return total_methylation / total_sites if total_sites > 0 else 0.0
    
    def get_current_stats(self) -> Dict:
        """
        Get current statistics about the PetriDish state.
        Always returns fresh, accurate values.
        
        Returns:
            Dictionary with current statistics
        """
        stats = {
            'num_cells': len(self.cells),
            'year': self.year,
            'mean_methylation': self._calculate_mean_methylation_proportion(),
            'mean_cell_jsd': statistics.mean([c.cell_JSD for c in self.cells]) if self.cells else 0.0
        }
        
        # Add gene metrics if we can calculate them
        if self.calculate_cell_jsds and self.cells:
            stats['gene_jsds'] = self.calculate_gene_jsds()
            stats['gene_mean_methylation'] = self.calculate_gene_mean_methylation()
            stats['n_genes'] = self.n_genes
        
        return stats
    
    def merge_metadata(self, other_metadata: Dict, preserve_keys: List[str] = None) -> 'PetriDish':
        """
        Merge metadata from another source while preserving critical values.
        
        Args:
            other_metadata: Metadata to merge in
            preserve_keys: Keys that should not be overwritten
            
        Returns:
            Self for method chaining
        """
        if not hasattr(self, 'metadata'):
            self.metadata = {}
        
        if preserve_keys is None:
            # By default, preserve critical calculated values
            preserve_keys = ['num_cells', 'gene_jsds', 'gene_mean_methylation']
        
        # Save values to preserve
        preserved = {}
        for key in preserve_keys:
            if key in self.metadata:
                preserved[key] = self.metadata[key]
        
        # Merge other metadata
        if other_metadata:
            self.metadata.update(other_metadata)
        
        # Restore preserved values
        self.metadata.update(preserved)
        
        # Always ensure num_cells is accurate
        self.metadata['num_cells'] = len(self.cells)
        
        return self
    
    def prepare_for_save(self, include_gene_metrics: bool = True) -> Dict:
        """
        Prepare the PetriDish data for saving, ensuring all metadata is current.
        
        Args:
            include_gene_metrics: Whether to include gene-level metrics
            
        Returns:
            Dictionary ready for JSON serialization
        """
        # Get fresh stats
        current_stats = self.get_current_stats()
        
        # Build the save data
        save_data = {
            'cells': [cell.to_dict() for cell in self.cells],
            'metadata': {}
        }
        
        # Start with existing metadata
        if hasattr(self, 'metadata'):
            save_data['metadata'] = self.metadata.copy()
        
        # Update with fresh values
        save_data['metadata']['num_cells'] = current_stats['num_cells']
        save_data['metadata']['year'] = self.year
        
        # Add gene metrics if requested
        if include_gene_metrics and self.calculate_cell_jsds:
            save_data['metadata']['gene_jsds'] = current_stats.get('gene_jsds', [])
            save_data['metadata']['gene_mean_methylation'] = current_stats.get('gene_mean_methylation', [])
            save_data['metadata']['n_genes'] = self.n_genes
        
        # Add rate configuration and parameters
        save_data['metadata']['gene_rate_groups'] = self.gene_rate_groups
        save_data['metadata']['gene_size'] = self.gene_size
        save_data['metadata']['n'] = self.n
        save_data['metadata']['growth_phase'] = self.growth_phase
        
        return save_data


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
            
            # Extract values (handle both naming conventions)
            jsd_values = []
            meth_values = []
            for cell in year_data:
                # Handle both 'cell_jsd' and 'cell_JSD' for compatibility
                if isinstance(cell, dict):
                    jsd = cell.get('cell_jsd', cell.get('cell_JSD', 0.0))
                    jsd_values.append(jsd)
                    # Handle missing methylation_proportion
                    if 'methylation_proportion' in cell:
                        meth_values.append(cell['methylation_proportion'] * 100)
                    else:
                        # Calculate from methylated array if available
                        if 'methylated' in cell:
                            meth_prop = sum(cell['methylated']) / len(cell['methylated'])
                            meth_values.append(meth_prop * 100)
                        else:
                            meth_values.append(0.0)
                else:
                    # Handle Cell objects
                    jsd_values.append(getattr(cell, 'cell_JSD', 0.0))
                    meth_values.append(getattr(cell, 'methylation_proportion', 0.0) * 100)
            
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
        Create Cell JSD plot with secondary y-axis for cell count.
        
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
                name='Mean Cell JSD',
                line=dict(color='rgb(99, 110, 250)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Cell JSD: %{y:.4f}<br>Population: %{customdata}<extra></extra>',
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
            title = f"Cell JSD Score vs Time"
        
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
        fig.update_yaxes(title_text="Cell JSD Score", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
        
        # Add annotation
        final_idx = -1
        final_year = years[final_idx]
        final_pop = self.stats['population_size'][final_idx]
        
        annotation_text = (
            f"<b>Final Statistics (Year {final_year}):</b><br>"
            f"Population: {final_pop} cells<br>"
            f"Cell JSD Mean: {self.stats['cell_jsd']['mean'][final_idx]:.4f}<br>"
            f"Cell JSD 25-75%: [{self.stats['cell_jsd']['p25'][final_idx]:.4f}, {self.stats['cell_jsd']['p75'][final_idx]:.4f}]<br>"
            f"Cell JSD 5-95%: [{self.stats['cell_jsd']['p5'][final_idx]:.4f}, {self.stats['cell_jsd']['p95'][final_idx]:.4f}]"
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
            subplot_titles=('Cell JSD Score vs Time', 'Methylation Proportion vs Time'),
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
                name='Mean Cell JSD',
                line=dict(color='rgb(99, 110, 250)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Cell JSD: %{y:.4f}<extra></extra>'
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
        fig.update_yaxes(title_text="Cell JSD Score", row=1, col=1, secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
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
            f"Cell JSD Mean: {self.stats['cell_jsd']['mean'][final_idx]:.4f}<br>"
            f"Cell JSD 25-75%: [{self.stats['cell_jsd']['p25'][final_idx]:.4f}, {self.stats['cell_jsd']['p75'][final_idx]:.4f}]<br>"
            f"Cell JSD 5-95%: [{self.stats['cell_jsd']['p5'][final_idx]:.4f}, {self.stats['cell_jsd']['p95'][final_idx]:.4f}]<br>"
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
        
        # If gene JSD history exists, plot basic gene JSD (keep the existing simple plot)
        if hasattr(self.petri, 'gene_jsd_history') and self.petri.gene_jsd_history:
            self.plot_gene_jsd(title, f"{base_path}_gene_jsd.png")
            # Note: Advanced gene JSD plots (heatmap, rate comparison) are now generated in phase2
        
        return self
    
    def plot_gene_jsd(self, title: str = None, output_path: str = None,
                      width: int = 1200, height: int = 600):
        """
        Create Gene JSD plot showing per-gene JSD evolution over time.
        
        Args:
            title: Optional title for the plot
            output_path: Optional path to save the plot
            width: Plot width in pixels
            height: Plot height in pixels
        """
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            print("Warning: Plotly not installed. Skipping gene JSD plot.")
            return
        
        # Check if gene JSD history exists
        if not hasattr(self.petri, 'gene_jsd_history') or not self.petri.gene_jsd_history:
            print("No gene JSD history available. Enable with track_gene_jsd=True")
            return
        
        # Extract data
        years = sorted([int(y) for y in self.petri.gene_jsd_history.keys()])
        n_genes = len(next(iter(self.petri.gene_jsd_history.values())))
        
        # Create figure with secondary y-axis for cell count
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        # Calculate statistics for each year
        mean_gene_jsds = []
        median_gene_jsds = []
        p25_gene_jsds = []
        p75_gene_jsds = []
        p5_gene_jsds = []
        p95_gene_jsds = []
        
        for year in years:
            gene_jsds = self.petri.gene_jsd_history[str(year)]
            mean_gene_jsds.append(np.mean(gene_jsds))
            median_gene_jsds.append(np.median(gene_jsds))
            p25_gene_jsds.append(np.percentile(gene_jsds, 25))
            p75_gene_jsds.append(np.percentile(gene_jsds, 75))
            p5_gene_jsds.append(np.percentile(gene_jsds, 5))
            p95_gene_jsds.append(np.percentile(gene_jsds, 95))
        
        # Add percentile bands
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=p95_gene_jsds + p5_gene_jsds[::-1],
                fill='toself',
                fillcolor='rgba(31, 119, 180, 0.1)',
                line=dict(width=0),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=p75_gene_jsds + p25_gene_jsds[::-1],
                fill='toself',
                fillcolor='rgba(31, 119, 180, 0.25)',
                line=dict(width=0),
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
                y=mean_gene_jsds,
                mode='lines',
                name='Mean',
                line=dict(color='#1f77b4', width=2),
                hovertemplate='Year: %{x}<br>Mean Gene JSD: %{y:.4f}<extra></extra>'
            ),
            secondary_y=False
        )
        
        # Add median line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=median_gene_jsds,
                mode='lines',
                name='Median',
                line=dict(color='#ff7f0e', width=2, dash='dash'),
                hovertemplate='Year: %{x}<br>Median Gene JSD: %{y:.4f}<extra></extra>'
            ),
            secondary_y=False
        )
        
        # Add cell count on secondary axis
        if self.stats and 'population_size' in self.stats:
            fig.add_trace(
                go.Scatter(
                    x=self.stats['years'],
                    y=self.stats['population_size'],
                    mode='lines',
                    name='Cell Count',
                    line=dict(color='gray', width=1, dash='dot'),
                    hovertemplate='Year: %{x}<br>Cells: %{y}<extra></extra>'
                ),
                secondary_y=True
            )
        
        # Add growth phase indicator
        growth_phase = self.detect_growth_phase()
        if growth_phase > 0 and years and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text="End of growth phase",
                annotation_position="top"
            )
        
        # Update layout
        if title is None:
            title = f"Gene JSD Score vs Time (Population-level, {n_genes} genes)"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            template='plotly_white',
            width=width,
            height=height,
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            )
        )
        
        # Set axis titles
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Gene JSD Score", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
        
        # Add final statistics annotation
        if years:
            final_year = years[-1]
            final_mean = mean_gene_jsds[-1]
            final_median = median_gene_jsds[-1]
            
            annotation_text = (
                f"<b>Final Gene JSD Statistics (Year {final_year}):</b><br>"
                f"Mean: {final_mean:.4f}<br>"
                f"Median: {final_median:.4f}<br>"
                f"25-75%: [{p25_gene_jsds[-1]:.4f}, {p75_gene_jsds[-1]:.4f}]<br>"
                f"5-95%: [{p5_gene_jsds[-1]:.4f}, {p95_gene_jsds[-1]:.4f}]<br>"
                f"Genes tracked: {n_genes}"
            )
            
            fig.add_annotation(
                text=annotation_text,
                xref="paper", yref="paper",
                x=0.98, y=0.02,
                showarrow=False,
                font=dict(size=10, family="monospace"),
                align="right",
                xanchor="right",
                yanchor="bottom",
                bgcolor="rgba(255, 255, 255, 0.9)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            )
        
        # Save if path provided
        if output_path:
            os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
            fig.write_image(output_path, width=width, height=height, scale=2)
            print(f"Gene JSD plot saved to {output_path}")
        else:
            fig.show()
        
        return fig

    def plot_gene_jsd_heatmap(self, title: str = None, output_path: str = None,
                              width: int = 1200, height: int = 800):
        """
        Create heatmap showing each gene's JSD evolution over time.
        
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
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not hasattr(self.petri, 'gene_jsd_history') or not self.petri.gene_jsd_history:
            raise ValueError("No gene JSD history found. Enable gene JSD tracking with --track-gene-jsd flag.")
        
        # Extract data from gene_jsd_history
        years = sorted(self.petri.gene_jsd_history.keys())
        n_genes = len(self.petri.gene_jsd_history[years[0]])
        
        # Create 2D array: rows = genes, columns = years
        heatmap_data = []
        for gene_idx in range(n_genes):
            gene_row = []
            for year in years:
                gene_jsd_value = self.petri.gene_jsd_history[year][gene_idx]
                gene_row.append(gene_jsd_value)
            heatmap_data.append(gene_row)
        
        # Create heatmap
        fig = go.Figure(data=go.Heatmap(
            z=heatmap_data,
            x=years,
            y=list(range(n_genes)),
            colorscale='Blues',
            colorbar=dict(title="Gene JSD Score"),
            hovertemplate='Year: %{x}<br>Gene: %{y}<br>JSD: %{z:.4f}<extra></extra>'
        ))
        
        # Add horizontal lines every 50 genes if using gene-rate-groups
        if hasattr(self.petri.cells[0], 'gene_rate_groups') and self.petri.cells[0].gene_rate_groups:
            gene_group_boundaries = []
            cumulative = 0
            for n_genes_in_group, _ in self.petri.cells[0].gene_rate_groups:
                cumulative += n_genes_in_group
                if cumulative < n_genes:  # Don't add line after last group
                    gene_group_boundaries.append(cumulative - 0.5)
            
            # Add horizontal lines
            for boundary in gene_group_boundaries:
                fig.add_hline(
                    y=boundary,
                    line_dash="solid",
                    line_color="white",
                    line_width=2,
                    opacity=0.8
                )
        
        # Calculate max JSD for annotation
        max_jsd = max(max(row) for row in heatmap_data)
        
        # Add annotation box
        stats_text = (f"<b>Heatmap Info</b><br>"
                     f"Genes: {n_genes}<br>"
                     f"Years: {len(years)}<br>"
                     f"Max JSD: {max_jsd:.3f}")
        
        fig.add_annotation(
            text=stats_text,
            xref="paper", yref="paper",
            x=0.98,
            y=0.97,
            showarrow=False,
            font=dict(size=11, family="Arial"),
            align="right",
            xanchor="right",
            yanchor="top",
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="#333333",
            borderwidth=1
        )
        
        # Update layout
        if title is None:
            title = "Gene JSD Evolution Heatmap"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            yaxis_title="Gene Index",
            width=width,
            height=height,
            template='plotly_white',
            showlegend=False
        )
        
        # Update axes with grid
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        
        # Save if path provided
        if output_path:
            os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
            fig.write_image(output_path, width=width, height=height, scale=2)
            print(f"Gene JSD heatmap saved to {output_path}")
        else:
            fig.show()
        
        return fig

    def plot_gene_jsd_by_rate_group(self, title: str = None, output_path: str = None,
                                    width: int = 1200, height: int = 600):
        """
        Line plot showing mean JSD for each gene rate group over time.
        
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
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not hasattr(self.petri, 'gene_jsd_history') or not self.petri.gene_jsd_history:
            raise ValueError("No gene JSD history found. Enable gene JSD tracking with --track-gene-jsd flag.")
        
        # Check if we have gene rate groups
        if not hasattr(self.petri.cells[0], 'gene_rate_groups') or not self.petri.cells[0].gene_rate_groups:
            raise ValueError("No gene rate groups found. This plot requires --gene-rate-groups parameter.")
        
        # Extract gene rate groups information
        gene_rate_groups = self.petri.cells[0].gene_rate_groups
        
        # Extract data from gene_jsd_history
        years = sorted(self.petri.gene_jsd_history.keys())
        
        # Calculate mean JSD for each rate group at each time point
        fig = go.Figure()
        
        # Define colors for each group (blue to red gradient)
        n_groups = len(gene_rate_groups)
        if n_groups == 4:
            colors = ['#0066cc', '#3399ff', '#ff6600', '#cc0000']
        else:
            # Generate color gradient
            import colorsys
            colors = []
            for i in range(n_groups):
                hue = 240 - (240 * i / max(1, n_groups - 1)) / 360  # Blue to red
                rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
                colors.append(f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})')
        
        # Calculate statistics for each group over time
        group_stats = []
        gene_start_idx = 0
        
        for group_idx, (n_genes, rate) in enumerate(gene_rate_groups):
            gene_end_idx = gene_start_idx + n_genes
            group_means = []
            group_p25 = []
            group_p75 = []
            
            for year in years:
                gene_jsds_for_group = self.petri.gene_jsd_history[year][gene_start_idx:gene_end_idx]
                group_means.append(np.mean(gene_jsds_for_group))
                group_p25.append(np.percentile(gene_jsds_for_group, 25))
                group_p75.append(np.percentile(gene_jsds_for_group, 75))
            
            group_stats.append({
                'mean': group_means,
                'p25': group_p25,
                'p75': group_p75,
                'rate': rate,
                'n_genes': n_genes
            })
            gene_start_idx = gene_end_idx
        
        # Add confidence bands (25-75 percentile) for each group
        for group_idx, stats in enumerate(group_stats):
            # Add confidence band
            fig.add_trace(
                go.Scatter(
                    x=years + years[::-1],
                    y=stats['p75'] + stats['p25'][::-1],
                    fill='toself',
                    fillcolor=colors[group_idx].replace('rgb', 'rgba').replace(')', ', 0.2)'),
                    line=dict(color='rgba(255,255,255,0)'),
                    name=f'Group {group_idx+1} (25-75%)',
                    showlegend=False,
                    hoverinfo='skip'
                )
            )
        
        # Add mean lines for each group
        for group_idx, stats in enumerate(group_stats):
            rate_pct = stats['rate'] * 100
            fig.add_trace(
                go.Scatter(
                    x=years,
                    y=stats['mean'],
                    mode='lines',
                    name=f'Group {group_idx+1}: {rate_pct:.1f}% rate',
                    line=dict(color=colors[group_idx], width=2),
                    hovertemplate=f'Year: %{{x}}<br>Mean Gene JSD: %{{y:.4f}}<br>Group {group_idx+1} ({stats["n_genes"]} genes @ {rate_pct:.1f}%)<extra></extra>'
                )
            )
        
        # Add growth phase indicator if available
        growth_phase = self.detect_growth_phase()
        if growth_phase > 0 and years and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text="End of growth phase",
                annotation_position="top"
            )
        
        # Add annotation box with group information
        group_info_lines = ["<b>Gene Groups:</b>"]
        for i, (n_genes, rate) in enumerate(gene_rate_groups):
            rate_pct = rate * 100
            group_info_lines.append(f"Group {i+1}: {n_genes} genes @ {rate_pct:.1f}%")
        
        group_info_text = "<br>".join(group_info_lines)
        
        fig.add_annotation(
            text=group_info_text,
            xref="paper", yref="paper",
            x=0.98,
            y=0.97,
            showarrow=False,
            font=dict(size=11, family="Arial"),
            align="right",
            xanchor="right",
            yanchor="top",
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="#333333",
            borderwidth=1
        )
        
        # Update layout
        if title is None:
            title = "Gene JSD by Methylation Rate Group"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            yaxis_title="Mean Gene JSD Score",
            width=width,
            height=height,
            template='plotly_white',
            showlegend=True,
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            ),
            hovermode='x unified'
        )
        
        # Update axes with grid
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        
        # Save if path provided
        if output_path:
            os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
            fig.write_image(output_path, width=width, height=height, scale=2)
            print(f"Gene rate group comparison saved to {output_path}")
        else:
            fig.show()
        
        return fig