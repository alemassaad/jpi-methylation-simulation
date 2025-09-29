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
        self.cell_methylation_proportion = 0.0
        
        self.methylation_distribution = [0.0 for _ in range(0, self.gene_size + 1)]
        self.methylation_distribution[0] = 1.0  # initially all genes are unmethylated
        self.cell_jsd = JS_div(self.methylation_distribution, self.baseline_methylation_distribution)
    
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
        if self.cell_methylation_proportion >= 1.0:
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
                
        self.cell_methylation_proportion = methylated_count / self.n
        self.compute_methylation_distribution()
        self.cell_jsd = JS_div(self.methylation_distribution, self.baseline_methylation_distribution)
    
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
        
        Stores essential cell data including rate configuration for verification.
        While gene_rate_groups is redundant (stored in parameters), including it
        makes each cell self-contained and enables verification of consistency.
        """
        # Calculate methylation stats
        n_methylated = sum(self.cpg_sites)
        
        return {
            'cpg_sites': self.cpg_sites[:],            # The methylation state (0s and 1s)
            'cell_jsd': self.cell_jsd,                 # The JSD value
            'age': self.age,                            # Cell age in years
            'gene_rate_groups': self.gene_rate_groups, # Rate configuration (for verification)
            'cell_methylation_proportion': n_methylated / self.n,  # Proportion methylated
            'n_methylated': n_methylated               # Count of methylated sites
        }
    
    @classmethod
    def from_dict(cls, data: Dict[str, Any], rate: float = None, 
                  gene_rate_groups: List[Tuple[int, float]] = None,
                  gene_size: int = GENE_SIZE) -> 'Cell':
        """
        Create a Cell from dictionary data.
        
        Args:
            data: Dictionary with cell data
            rate: Uniform methylation rate (will be converted to gene_rate_groups)
            gene_rate_groups: Gene-specific rates (can be overridden from data)
            gene_size: Sites per gene
            
        Returns:
            Reconstructed Cell object
        """
        # Get n from cpg_sites
        n = len(data['cpg_sites'])
        
        # Use gene_rate_groups from data if available (new format)
        # Otherwise fall back to provided parameters (backward compatibility)
        if 'gene_rate_groups' in data:
            actual_groups = data['gene_rate_groups']
            # Convert list of lists to list of tuples if needed
            if actual_groups and isinstance(actual_groups[0], list):
                actual_groups = [tuple(group) for group in actual_groups]
            gene_rate_groups = actual_groups
        
        # Let __init__ handle the conversion
        cell = cls(n=n, rate=rate, gene_rate_groups=gene_rate_groups, gene_size=gene_size)
        # Set cpg_sites from data
        cell.cpg_sites = data['cpg_sites'][:]
        # Handle both old and new naming for backward compatibility
        cell.cell_jsd = data.get('cell_jsd', data.get('cell_JSD', 0.0))
        cell.age = data.get('age', 0)  # Restore age (default to 0 for old data)
        
        # Verify stored stats if available (for validation)
        # Handle both old and new keys for backward compatibility during transition
        if 'cell_methylation_proportion' in data:
            actual_prop = sum(cell.cpg_sites) / n
            stored_prop = data['cell_methylation_proportion']
            if abs(actual_prop - stored_prop) > 0.0001:
                print(f"Warning: Cell methylation proportion mismatch: stored={stored_prop:.4f}, actual={actual_prop:.4f}")
        elif 'methylation_proportion' in data:
            # Legacy support - can remove after transition
            actual_prop = sum(cell.cpg_sites) / n
            stored_prop = data['methylation_proportion']
            if abs(actual_prop - stored_prop) > 0.0001:
                print(f"Warning: Cell methylation proportion mismatch: stored={stored_prop:.4f}, actual={actual_prop:.4f}")
        
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
                 seed: int = None, growth_phase: Optional[int] = DEFAULT_GROWTH_PHASE, 
                 cells: List['Cell'] = None) -> None:
        """
        Initialize petri dish with a single unmethylated cell or provided cells.
        
        Args:
            rate: Uniform methylation rate per site per year (will be converted to gene_rate_groups)
            gene_rate_groups: Gene-specific rates as [(n_genes, rate), ...]
            n: Number of CpG sites per cell
            gene_size: Number of sites per gene
            seed: Random seed for reproducibility
            growth_phase: Duration of growth phase in years (target = 2^growth_phase cells).
                         Use None for static populations (e.g., snapshots) that won't be aged.
            cells: Optional list of cells to start with (for phase2 compatibility).
                   All cells must have identical gene_rate_groups configuration and ages.
            
        Raises:
            ValueError: If provided cells have different gene_rate_groups configurations or ages
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
        
        # Validate growth_phase if provided
        if growth_phase is not None:
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
        self.target_population = 2 ** growth_phase if growth_phase is not None else None
        
        # Initialize cells - either provided or single unmethylated cell
        if cells is not None:
            self.cells = cells
            # Validate all cells have same rate configuration
            self.validate_cell_consistency()
            
            # NEW: Validate and set year from cell ages
            if self.cells:
                ages = [cell.age for cell in self.cells]
                unique_ages = set(ages)
                if len(unique_ages) > 1:
                    raise ValueError(
                        f"All cells must have the same age. Found ages: {sorted(unique_ages)}"
                    )
                # Set PetriDish year to match cell age
                self.year = ages[0] if ages else 0
                self.initial_year = self.year  # Track when this PetriDish started
            else:
                self.year = 0
                self.initial_year = 0  # Starting from year 0
        else:
            # Create initial cell - single path now
            self.cells = [Cell(n=n, gene_rate_groups=gene_rate_groups, gene_size=gene_size)]
            self.year = 0  # Fresh cell starts at age 0
            self.initial_year = 0  # Starting from year 0
        
        # Cell history tracking (always enabled)
        # Gene JSD tracking (always enabled)
        self.gene_jsd_history = {}
        # Baseline must match gene_size + 1 bins
        self.BASELINE_GENE_DISTRIBUTION = [1.0] + [0.0] * self.gene_size
        
        # Store initial state if we have cells (renamed for clarity)
        if self.cells:
            # Use actual year as key, not always '0'
            year_key = str(self.year)
            self.cell_history = {year_key: [cell.to_dict() for cell in self.cells]}
            
            # Calculate gene JSD for all years (including year 0)
            gene_jsds = self.calculate_gene_jsd()
            self.gene_jsd_history[year_key] = gene_jsds
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
        
    def methylate_cells(self) -> None:
        """
        Apply stochastic methylation to all cells.
        Each cell independently accumulates methylation.
        """
        for cell in self.cells:
            cell.methylate()
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
        Chooses appropriate aging strategy based on years since initial_year.
        
        Raises:
            ValueError: If growth_phase is None (static population)
        """
        if self.growth_phase is None:
            raise ValueError(
                "Cannot age a static PetriDish (created with growth_phase=None). "
                "Static PetriDishes are snapshots and cannot undergo growth or homeostasis."
            )
        
        # Determine if we should grow or maintain based on years since initial_year
        if not hasattr(self, 'initial_year'):
            self.initial_year = 0  # Backward compatibility
        
        years_since_start = self.year - self.initial_year
        
        if years_since_start < self.growth_phase:
            # Use exponential growth for one year
            self.grow_exponentially(1, verbose=True)
        else:
            # Use homeostasis for one year
            self.maintain_homeostasis(1, target_population=self.target_population, verbose=True)
        
    def grow_exponentially(self, n_years: int, verbose: bool = True) -> None:
        """
        Grow the population exponentially for n_years.
        Cells divide and methylate each year without culling.
        
        Args:
            n_years: Number of years to grow exponentially
            verbose: Whether to print progress
        """
        for i in range(n_years):
            self.year += 1
            if verbose:
                print(f"\nYear {self.year}:")
                print(f"  Exponential growth (year {i+1} of {n_years})")
            
            # Growth phase operations
            initial_count = len(self.cells)
            self.divide_cells()
            self.methylate_cells()
            final_count = len(self.cells)
            
            if verbose:
                print(f"  Final count: {final_count} cells (2^{i+1} = {2**(i+1)} expected)")
            
            # Record history (always enabled now)
            self._record_history(self.year)
            
            # Report statistics
            if verbose:
                jsd_values = [cell.cell_jsd for cell in self.cells]
                mean_jsd = statistics.mean(jsd_values) if jsd_values else 0.0
                print(f"  Mean cell JSD: {mean_jsd:.4f}")
    
    def maintain_homeostasis(self, n_years: int, target_population: int = None, verbose: bool = True) -> None:
        """
        Maintain population in homeostasis for n_years.
        Cells divide, then ~50% are culled, then methylate.
        
        Args:
            n_years: Number of years to maintain homeostasis
            target_population: Expected population size (for display only)
            verbose: Whether to print progress
        """
        if target_population is None:
            target_population = len(self.cells)
        
        for i in range(n_years):
            self.year += 1
            if verbose:
                print(f"\nYear {self.year}:")
                print(f"  Homeostasis (year {i+1} of {n_years})")
            
            # Homeostasis operations
            initial_count = len(self.cells)
            self.divide_cells()
            intermediate_count = len(self.cells)
            self.random_cull_cells()
            self.methylate_cells()
            final_count = len(self.cells)
            
            if verbose:
                print(f"  Final count: {final_count} cells (target ~{target_population})")
            
            # Check for extinction
            if final_count == 0:
                if verbose:
                    print(f"  WARNING: Population extinct at year {self.year}")
                # Still record the history even if extinct
            
            # Record history
            self._record_history(self.year)
            
            # Report statistics
            if verbose and self.cells:
                jsd_values = [cell.cell_jsd for cell in self.cells]
                mean_jsd = statistics.mean(jsd_values) if jsd_values else 0.0
                print(f"  Mean cell JSD: {mean_jsd:.4f}")
    
    def run_simulation(self, t_max: int = T_MAX) -> None:
        """
        Run complete simulation from year 0 to t_max.
        
        Args:
            t_max: Maximum simulation time in years
            
        Raises:
            ValueError: If growth_phase is None (static population)
        """
        if self.growth_phase is None:
            raise ValueError(
                "Cannot run simulation on a static PetriDish (created with growth_phase=None). "
                "Static PetriDishes are snapshots and cannot be simulated."
            )
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
        
        print(f"\nYear {self.year}: Starting with {len(self.cells)} cell(s)")
        
        # Use the new methods based on growth_phase
        if not hasattr(self, 'initial_year'):
            self.initial_year = 0  # Backward compatibility
            
        # Calculate remaining growth years
        years_grown = self.year - self.initial_year
        growth_years_remaining = max(0, self.growth_phase - years_grown)
        total_years_to_simulate = t_max - self.year
        homeostasis_years = max(0, total_years_to_simulate - growth_years_remaining)
        
        if growth_years_remaining > 0:
            years_to_grow = min(growth_years_remaining, total_years_to_simulate)
            print(f"\nGrowth phase: {years_to_grow} years")
            self.grow_exponentially(years_to_grow)
        
        if homeostasis_years > 0 and self.year < t_max:
            print(f"\nHomeostasis phase: {homeostasis_years} years")
            self.maintain_homeostasis(homeostasis_years, target_population=self.target_population)
        
        print("\n" + "="*60)
        print("Simulation complete!")
        print(f"Final population: {len(self.cells)} cells")
        print("="*60)
        
    def save_history(self, filename: str = None, directory: str = "data", compress: bool = True) -> str:
        """
        Save simulation history to JSON file (compressed or uncompressed) using hierarchical structure.
        
        Args:
            filename: If provided as absolute path, use it directly.
                     Otherwise, auto-generate filename.
            directory: Base output directory (ignored if filename is absolute)
            compress: If True, save as .json.gz; if False, save as .json
            
        Returns:
            Path to saved file
        """
        from datetime import datetime
        
        # Check if filename is an absolute path
        if filename and os.path.isabs(filename):
            # Use the provided absolute path directly
            filepath = filename
            dir_path = os.path.dirname(filepath)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path, exist_ok=True)
        else:
            # Original hierarchical path generation
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
            'config': {
                'gene_rate_groups': self.gene_rate_groups,  # Only store this
                'n': self.n,
                'gene_size': self.gene_size,
                'growth_phase': self.growth_phase,
                'years': self.year,
                'seed': self.seed
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
    
    
    def _record_history(self, year: int = None) -> None:
        """
        Internal method to record current state in history.
        
        Args:
            year: Optional year to record at. If None, uses self.year
        """
        if year is None:
            year = self.year
        
        # Record cell history (always enabled)
        self.cell_history[str(year)] = [cell.to_dict() for cell in self.cells]
        
        # Record gene JSD (always enabled)
        gene_jsds = self.calculate_gene_jsd()
        self.gene_jsd_history[str(year)] = gene_jsds
    
    
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
        Calculate Jensen-Shannon Divergence for each gene across all cells.
        
        Returns:
            List of JSD values, one per gene (0.0 to 1.0)
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
    
    @property
    def gene_methylation_proportions(self) -> List[float]:
        """
        Get current gene methylation proportions for the population.
        Calculates on-demand, providing clean OOP access.
        
        Returns:
            List of methylation proportions (0.0 to 1.0), one per gene
        """
        return self.calculate_gene_mean_methylation()
    
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
    def from_cells(cls, cells, growth_phase: int = 7,
                   metadata: Dict = None) -> 'PetriDish':
        """
        Create a PetriDish from cell(s).
        
        Args:
            cells: Single Cell object or list of Cell objects to use as founding population
            growth_phase: The growth phase for this PetriDish
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
                seed=None
            )
        else:
            # Gene-specific rates
            petri = cls(
                cells=cells_list,  # Pass cells directly - no replacement needed!
                n=first_cell.n,
                gene_rate_groups=first_cell.gene_rate_groups,
                gene_size=first_cell.gene_size,
                growth_phase=growth_phase,
                seed=None
            )
        
        # Set metadata if provided
        if metadata:
            petri.metadata = metadata.copy()
        else:
            petri.metadata = {}
        
        # Don't auto-add metadata - let phase2 control what's stored
        
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
    
    def _calculate_mean_cell_methylation_proportion(self) -> float:
        """
        Calculate the mean cell methylation proportion across all cells.
        
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
            'mean_cell_methylation_proportion': self._calculate_mean_cell_methylation_proportion(),
            'mean_cell_jsd': statistics.mean([c.cell_jsd for c in self.cells]) if self.cells else 0.0
        }
        
        # Add gene metrics
        if self.cells:
            stats['gene_jsds'] = self.calculate_gene_jsd()
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
            # By default, preserve nothing - let caller decide
            preserve_keys = []
        
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
        
        return self
    
    def prepare_for_save(self, include_gene_metrics: bool = False) -> Dict:
        """
        Prepare the PetriDish data for saving with minimal metadata.
        
        Args:
            include_gene_metrics: Whether to include gene-level metrics (deprecated, ignored)
            
        Returns:
            Dictionary ready for JSON serialization
        """
        # Build the save data
        save_data = {
            'cells': [cell.to_dict() for cell in self.cells],
            'metadata': {}
        }
        
        # Only include explicitly set metadata
        if hasattr(self, 'metadata'):
            save_data['metadata'] = self.metadata.copy()
        
        return save_data


# Plotting functionality is available in phase3
# To visualize simulation data:
#   cd phase3
#   python run_pipeline.py --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/
