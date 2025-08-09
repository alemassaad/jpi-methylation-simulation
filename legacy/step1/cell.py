import random
import math
import statistics
import json
import os
from typing import List, Dict, Any


def KL_div(P: List[float], Q: List[float]) -> float:
    if len(P) != len(Q):
        raise ValueError("Both distributions must have the same length.") 
    
    n = len(P)
    
    for i in range(n):
        if Q[i] == 0 and P[i] > 0:
                return float('inf')
    
    return sum([P[i] * math.log2(P[i]/Q[i]) for i in range(n) if P[i] > 0])


def JS_div(P: List[float], Q: List[float]) -> float:
    if len(P) != len(Q):
        raise ValueError("Both distributions must have the same length.")

    n = len(P)

    M = [(P[i] + Q[i]) / 2 for i in range(n)]

    return (KL_div(P, M) + KL_div(Q, M)) / 2


N = 1000
RATE = 0.01
GENE_SIZE = 5
BASELINE_METHYLATION_DISTRIBUTION = [1.0] + [0.0 for _ in range(GENE_SIZE)]

M = 10_000 
T_MAX = 100


class Cell:
    """Represents a cell with methylation sites that age over time.
    
    Attributes:
        n: Number of CpG sites
        rate: Methylation rate per year
        gene_size: Number of CpG sites per gene
        cpg_sites: List of methylation states (0 or 1)
        age: Current age in years
        methylation_proportion: Fraction of methylated sites
        methylation_distribution: Distribution of methylation across genes
        JSD: Jensen-Shannon divergence from baseline
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
        self.methylation_proportion = 0
        
        self.methylation_distribution = [0.0 for _ in range(0, self.gene_size + 1)]
        self.methylation_distribution[0] = 1.0  # initially all cells are unmethylated
        self.JSD = JS_div(self.methylation_distribution, self.baseline_methylation_distribution)
        
    def age_1_year(self) -> None:
        """Age the cell by one year, randomly methylating unmethylated sites.""" 
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
        
        # print(f"Cell age: {self.age}, Percentage of Methylation: {self.methylation_proportion:.2%}, Methylation Distribution: {self.methylation_distribution},  JSD: {self.JSD:.4f}")

    def age_multiple_years(self, years: int) -> None:
        """Age the cell for multiple years."""
        for _ in range(years):
            self.age_1_year()
            
    def compute_methylation_distribution(self) -> None:
        """Calculate the distribution of methylation levels across genes."""
        distribution = [0.0 for _ in range(0, self.gene_size + 1)]
        
        # Optimized: avoid creating sublists, count methylated sites directly
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
            'cpg_sites': self.cpg_sites[:],  # Faster slicing instead of copy()
            'methylation_proportion': self.methylation_proportion,
            'methylation_distribution': self.methylation_distribution[:],  # Faster slicing instead of copy()
            'jsd': self.JSD,
            'age': self.age,
            'gene_size': self.gene_size,
            'rate': self.rate
        }


class History:
    """Manages a population of cells and tracks their state over time.
    
    Attributes:
        m: Number of cells in the population
        n: Number of CpG sites per cell
        rate: Methylation rate
        gene_size: Sites per gene
        t_max: Maximum simulation time
        cells: List of Cell objects
        history: Dictionary mapping years to cell states
    """
    def __init__(self, m: int = M, n: int = N, rate: float = RATE, 
                 gene_size: int = GENE_SIZE, t_max: int = T_MAX) -> None:
        self.m = m
        self.n = n
        self.rate = rate
        self.gene_size = gene_size
        self.t_max = t_max
        
        # Create cells directly
        if self.n % self.gene_size != 0:
            raise ValueError("Gene size must divide the number of cells evenly.")
        
        self.cells = [Cell(self.n, self.rate, self.gene_size) for _ in range(self.m)]
        
        # Store initial state
        self.history = {0: [cell.to_dict() for cell in self.cells]}
        
    def age_multiple_years(self, years: int) -> None:
        """Age all cells in the population for the specified number of years."""
        for year in range(years):
            print(f"Aging year {year + 1}...")
            
            # Age all cells
            for cell in self.cells:
                cell.age_1_year()
            
            # Store current state
            current_age = self.cells[0].age  # All cells have same age
            self.history[current_age] = [cell.to_dict() for cell in self.cells]

    def save_history(self, filename: str = "simulation_history", directory: str = "history") -> None:
        """Save the simulation history to a JSON file in the specified directory."""
        # Create directory if it doesn't exist
        if not os.path.exists(directory):
            os.makedirs(directory)
            print(f"Created directory: {directory}")
        
        # Construct full path
        filepath = os.path.join(directory, filename + ".json")
        print(f"Saving history to {filepath}")
        
        # Custom formatting to keep cpg_sites compact
        output = "{\n"
        years = sorted(self.history.keys())
        for i, year in enumerate(years):
            output += f'  "{year}": [\n'
            cells = self.history[year]
            for j, cell in enumerate(cells):
                output += "    {\n"
                output += f'      "cpg_sites": {cell["cpg_sites"]},\n'
                output += f'      "methylation_proportion": {cell["methylation_proportion"]},\n'
                output += f'      "methylation_distribution": {cell["methylation_distribution"]},\n'
                output += f'      "jsd": {cell["jsd"]},\n'
                output += f'      "age": {cell["age"]},\n'
                output += f'      "gene_size": {cell["gene_size"]},\n'
                output += f'      "rate": {cell["rate"]}\n'
                output += "    }"
                if j < len(cells) - 1:
                    output += ","
                output += "\n"
            output += "  ]"
            if i < len(years) - 1:
                output += ","
            output += "\n"
        output += "}"
        
        with open(filepath, 'w') as f:
            f.write(output)
        print("Save complete.")