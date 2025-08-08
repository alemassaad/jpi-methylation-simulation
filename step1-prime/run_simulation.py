#!/usr/bin/env python3
"""
Step1-Prime: Biologically realistic methylation simulation with growth and homeostasis.

This simulation starts with a single unmethylated cell and models:
1. Growth phase (years 0-13): Exponential growth to 8192 cells
2. Steady state (years 14+): Maintain population through division and culling

Usage:
    python run_simulation.py --rate 0.005 --years 100 --seed 42
"""

import argparse
import time
from cell import PetriDish, RATE, N, GENE_SIZE, DEFAULT_GROWTH_YEARS, T_MAX


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run Step1-Prime methylation simulation with growth and homeostasis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Simulation parameters
    parser.add_argument('-r', '--rate', type=float, default=RATE,
                        help='Methylation rate per site per year (0.01 = 1%%, 0.005 = 0.5%%)')
    parser.add_argument('-n', '--sites', type=int, default=N,
                        help='Number of CpG sites per cell')
    parser.add_argument('-t', '--years', type=int, default=T_MAX,
                        help='Number of years to simulate')
    parser.add_argument('-g', '--gene-size', type=int, default=GENE_SIZE,
                        help='Size of genes (sites per gene)')
    parser.add_argument('--growth-years', type=int, default=DEFAULT_GROWTH_YEARS,
                        help='Years of exponential growth phase (target = 2^growth_years cells). '
                             'Default: 13 (8192 cells). Range: 1-20')
    
    # Optional parameters
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output filename (auto-generated if not specified)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility')
    parser.add_argument('--no-save', action='store_true',
                        help='Run simulation without saving results')
    
    return parser.parse_args()


def main():
    """Main entry point for the simulation."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Extract parameters
    rate = args.rate
    n = args.sites
    t_max = args.years
    gene_size = args.gene_size
    growth_years = args.growth_years
    seed = args.seed
    
    # Validate parameters
    if n % gene_size != 0:
        print(f"Error: Number of sites ({n}) must be divisible by gene size ({gene_size})")
        return 1
    
    # Validate growth_years
    if growth_years < 1 or growth_years > 20:
        print(f"Error: growth-years must be between 1 and 20, got {growth_years}")
        return 1
    
    # Start timing
    start_time = time.time()
    
    # Create and run simulation
    print("\nInitializing PetriDish simulation...")
    petri_dish = PetriDish(
        rate=rate,
        n=n,
        gene_size=gene_size,
        seed=seed,
        growth_years=growth_years
    )
    
    # Run simulation
    petri_dish.run_simulation(t_max=t_max)
    
    # Calculate runtime
    runtime = time.time() - start_time
    print(f"\nSimulation runtime: {runtime:.2f} seconds")
    
    # Save results
    if not args.no_save:
        filepath = petri_dish.save_history(filename=args.output)
        print(f"\nResults saved to: {filepath}")
    else:
        print("\nResults not saved (--no-save flag)")
    
    # Print summary statistics
    print("\nFinal Statistics:")
    print(f"  Total cells: {len(petri_dish.cells)}")
    
    # Calculate JSD statistics
    jsd_values = [cell.JSD for cell in petri_dish.cells]
    if jsd_values:
        import statistics
        print(f"  Mean JSD: {statistics.mean(jsd_values):.4f}")
        print(f"  Median JSD: {statistics.median(jsd_values):.4f}")
        print(f"  Std Dev JSD: {statistics.stdev(jsd_values):.4f}")
        print(f"  Min JSD: {min(jsd_values):.4f}")
        print(f"  Max JSD: {max(jsd_values):.4f}")
    
    # Calculate methylation statistics
    meth_props = [cell.methylation_proportion for cell in petri_dish.cells]
    if meth_props:
        print(f"  Mean methylation: {statistics.mean(meth_props):.2%}")
        print(f"  Median methylation: {statistics.median(meth_props):.2%}")
    
    print("\nSimulation complete!")
    return 0


if __name__ == "__main__":
    exit(main())