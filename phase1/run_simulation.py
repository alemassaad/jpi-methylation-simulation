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
from cell import PetriDish, RATE, N, GENE_SIZE, DEFAULT_GROWTH_PHASE, T_MAX


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run Step1-Prime methylation simulation with growth and homeostasis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Simulation parameters
    parser.add_argument('-r', '--rate', type=float, default=None, nargs='?', const=RATE,
                        help='Uniform methylation rate per site per year (default: 0.005 if flag used without value)')
    parser.add_argument('-n', '--sites', type=int, default=N,
                        help='Number of CpG sites per cell')
    parser.add_argument('-t', '--years', type=int, default=T_MAX,
                        help='Number of years to simulate')
    parser.add_argument('-g', '--gene-size', type=int, default=GENE_SIZE,
                        help='Size of genes (sites per gene)')
    parser.add_argument('--growth-phase', type=int, default=DEFAULT_GROWTH_PHASE,
                        help='Duration of growth phase in years (cells double each year). '
                             'Final population = 2^growth_phase. Default: 13 (8192 cells). Range: 1-20')
    parser.add_argument('--gene-rate-groups', type=str, default=None,
                        help='Gene-specific rates as "n1:rate1,n2:rate2,..." '
                             'Example: "50:0.004,50:0.0045,50:0.005,50:0.0055" '
                             'Total genes must equal n_sites/gene_size')
    
    # Optional parameters
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output filename (auto-generated if not specified)')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility (default: 42, use -1 for no seed)')
    parser.add_argument('--no-save', action='store_true',
                        help='Run simulation without saving results')
    
    # Performance and tracking options
    parser.add_argument('--no-gene-jsd', action='store_true',
                        help='Disable gene JSD tracking (enabled by default)')
    parser.add_argument('--no-jsds', action='store_true',
                        help='Disable all JSD calculations for maximum performance')
    
    return parser.parse_args()


def parse_gene_rate_groups(groups_str: str, n_sites: int, gene_size: int):
    """
    Parse gene rate groups string into list of tuples.
    
    Args:
        groups_str: String like "50:0.004,50:0.0045,50:0.005,50:0.0055"
        n_sites: Total number of CpG sites
        gene_size: Sites per gene
        
    Returns:
        List of (n_genes, rate) tuples
    """
    groups = []
    
    try:
        for i, group in enumerate(groups_str.split(',')):
            if ':' not in group:
                raise ValueError(f"Group {i} missing ':' separator: '{group}'")
            
            parts = group.split(':')
            if len(parts) != 2:
                raise ValueError(f"Group {i} should have format 'n:rate', got: '{group}'")
            
            n_genes = int(parts[0])
            rate = float(parts[1])
            
            if n_genes <= 0:
                raise ValueError(f"Group {i}: number of genes must be positive, got {n_genes}")
            if rate <= 0:
                raise ValueError(f"Group {i}: rate must be positive, got {rate}")
            
            groups.append((n_genes, rate))
            
    except ValueError as e:
        raise ValueError(f"Invalid gene-rate-groups format: {e}")
    
    # Validate total
    total_genes = sum(n for n, _ in groups)
    expected_genes = n_sites // gene_size
    if total_genes != expected_genes:
        raise ValueError(
            f"gene-rate-groups specifies {total_genes} genes, but simulation has {expected_genes} genes "
            f"(n_sites={n_sites}, gene_size={gene_size})"
        )
    
    return groups


def main():
    """Main entry point for the simulation."""
    # Parse command-line arguments
    args = parse_arguments()
    
    # Extract parameters
    n = args.sites
    t_max = args.years
    gene_size = args.gene_size
    growth_phase = args.growth_phase
    seed = args.seed if args.seed != -1 else None  # -1 means no seed
    
    # Validate rate specification
    if args.rate is not None and args.gene_rate_groups is not None:
        print("Error: Cannot specify both --rate and --gene-rate-groups")
        return 1
    
    if args.rate is None and args.gene_rate_groups is None:
        print("Error: Must specify either --rate (uniform) or --gene-rate-groups (gene-specific)")
        return 1
    
    # Parse gene rate groups if specified
    gene_rate_groups = None
    if args.gene_rate_groups:
        try:
            gene_rate_groups = parse_gene_rate_groups(
                args.gene_rate_groups, n, gene_size
            )
        except ValueError as e:
            print(f"Error: {e}")
            return 1
    
    # Validate parameters
    if n % gene_size != 0:
        print(f"Error: Number of sites ({n}) must be divisible by gene size ({gene_size})")
        return 1
    
    # Validate growth_phase
    if growth_phase < 1 or growth_phase > 20:
        print(f"Error: growth-phase must be between 1 and 20, got {growth_phase}")
        return 1
    
    # Start timing
    start_time = time.time()
    
    # Create and run simulation with appropriate rate configuration
    print("\nInitializing PetriDish simulation...")
    if args.rate is not None:
        petri_dish = PetriDish(
            rate=args.rate,
            n=n,
            gene_size=gene_size,
            seed=seed,
            growth_phase=growth_phase,
            calculate_jsds=not args.no_jsds
        )
    else:
        petri_dish = PetriDish(
            gene_rate_groups=gene_rate_groups,
            n=n,
            gene_size=gene_size,
            seed=seed,
            growth_phase=growth_phase,
            calculate_jsds=not args.no_jsds
        )
    
    # Enable history tracking with gene JSD by default (unless disabled)
    if not hasattr(args, 'no_gene_jsd') or not args.no_gene_jsd:
        petri_dish.enable_history_tracking(start_year=0, track_gene_jsd=True)
        print("Gene JSD tracking enabled (default behavior)")
    else:
        print("Gene JSD tracking disabled by --no-gene-jsd flag")
    
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
    jsd_values = [cell.cell_JSD for cell in petri_dish.cells]
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