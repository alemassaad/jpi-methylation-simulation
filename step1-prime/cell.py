import random
import math
import statistics
import json
import gzip
import os
import copy
from typing import List, Dict, Any


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
    
    def __init__(self, n: int = N, rate: float = RATE, gene_size: int = GENE_SIZE, 
                 baseline_methylation_distribution: List[float] = BASELINE_METHYLATION_DISTRIBUTION) -> None:
        self.n = n
        self.rate = rate
        self.gene_size = gene_size
        self.baseline_methylation_distribution = baseline_methylation_distribution
        
        if self.n % self.gene_size != 0:
            raise ValueError("Gene size must divide the number of cells evenly.")
    
        self.cpg_sites = [0 for _ in range(n)]  
        self.age = 0
        self.methylation_proportion = 0.0
        
        self.methylation_distribution = [0.0 for _ in range(0, self.gene_size + 1)]
        self.methylation_distribution[0] = 1.0  # initially all genes are unmethylated
        self.JSD = JS_div(self.methylation_distribution, self.baseline_methylation_distribution)
        
    def methylate(self) -> None:
        """
        Apply stochastic methylation to unmethylated sites.
        This is the core aging mechanism (formerly age_1_year).
        """
        self.age += 1
        
        # Early exit if already fully methylated
        if self.methylation_proportion >= 1.0:
            return
        
        # Optimized: count methylated sites while updating
        methylated_count = 0
        for i in range(self.n):
            if self.cpg_sites[i] == 0:
                if random.random() < self.rate:
                    self.cpg_sites[i] = 1
                    methylated_count += 1
            else:
                methylated_count += 1
                
        self.methylation_proportion = methylated_count / self.n
        self.compute_methylation_distribution()
        self.JSD = JS_div(self.methylation_distribution, self.baseline_methylation_distribution)
    
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
        return {
            'cpg_sites': self.cpg_sites[:],
            'methylation_proportion': self.methylation_proportion,
            'methylation_distribution': self.methylation_distribution[:],
            'jsd': self.JSD,
            'age': self.age,
            'gene_size': self.gene_size,
            'rate': self.rate
        }


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
    
    def __init__(self, rate: float = RATE, n: int = N, gene_size: int = GENE_SIZE, 
                 seed: int = None, growth_phase: int = DEFAULT_GROWTH_PHASE) -> None:
        """
        Initialize petri dish with a single unmethylated cell.
        
        Args:
            rate: Methylation rate per site per year
            n: Number of CpG sites per cell
            gene_size: Number of sites per gene
            seed: Random seed for reproducibility
            growth_phase: Duration of growth phase in years (target = 2^growth_phase cells)
        """
        # Validate growth_phase
        if growth_phase < 1:
            raise ValueError(f"growth_phase must be >= 1, got {growth_phase}")
        if growth_phase > 20:
            raise ValueError(f"growth_phase must be <= 20 (max 1M cells), got {growth_phase}")
        
        if seed is not None:
            random.seed(seed)
            
        self.rate = rate
        self.n = n
        self.gene_size = gene_size
        self.seed = seed
        self.growth_phase = growth_phase
        self.target_population = 2 ** growth_phase  # Calculate from growth_phase
        
        # Start with single unmethylated cell
        self.cells = [Cell(n=n, rate=rate, gene_size=gene_size)]
        self.year = 0
        
        # Store initial state
        self.history = {'0': [self.cells[0].to_dict()]}
        
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
        self.history[str(self.year)] = [cell.to_dict() for cell in self.cells]
        
        # Report statistics
        jsd_values = [cell.JSD for cell in self.cells]
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
        print(f"  Methylation rate: {self.rate:.3%}")
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
        
    def save_history(self, filename: str = None, directory: str = "data") -> None:
        """
        Save simulation history to compressed JSON file.
        
        Args:
            filename: Output filename (auto-generated if None)
            directory: Output directory
        """
        # Create directory if it doesn't exist
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")
        
        # Generate filename if not provided
        if filename is None:
            # Always include seed info in filename
            seed_str = f"_seed{self.seed}" if self.seed is not None else "_noseed"
            # Use actual final population count, not target
            final_pop = len(self.cells)
            # Include growth_phase in filename
            filename = f"simulation_rate_{self.rate:.6f}_g{self.growth_phase}_m{final_pop}_n{self.n}_t{self.year}{seed_str}"
        
        filepath = os.path.join(directory, filename + ".json.gz")
        print(f"\nSaving compressed history to {filepath}")
        
        save_start = statistics.mean([0])  # Just to have time module if needed
        import time
        save_start = time.time()
        
        # Use gzip compression and native JSON serialization
        with gzip.open(filepath, 'wt', encoding='utf-8', compresslevel=1) as f:
            json.dump(self.history, f, separators=(',', ':'))
        
        save_time = time.time() - save_start
        file_size_mb = os.path.getsize(filepath) / (1024 * 1024)
        
        print(f"  Save time: {save_time:.2f} seconds")
        print(f"  File size: {file_size_mb:.2f} MB (compressed)")
        
        return filepath