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
import os
import sys
import copy
from typing import Dict, Any, Optional
from cell import PetriDish, RATE, N, GENE_SIZE, DEFAULT_GROWTH_PHASE, T_MAX, rate_to_gene_rate_groups

# Try to import yaml, provide helpful error if not available
try:
    import yaml
except ImportError:
    yaml = None
    print("Warning: PyYAML not installed. Config file support disabled.")
    print("Install with: pip install pyyaml")


def deep_merge(base: Dict, override: Dict) -> Dict:
    """Deep merge two dictionaries, with override taking precedence."""
    result = copy.deepcopy(base)
    for key, value in override.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """
    Load configuration from YAML file(s).
    
    Args:
        config_path: Path to user config file (optional)
        
    Returns:
        Merged configuration dictionary
    """
    config = {}
    
    if yaml is None:
        return config
    
    # Load default config if it exists
    default_path = os.path.join(os.path.dirname(__file__), "config_default.yaml")
    if os.path.exists(default_path):
        try:
            with open(default_path, 'r') as f:
                config = yaml.safe_load(f) or {}
        except Exception as e:
            print(f"Warning: Could not load default config: {e}")
    
    # Override with user config if provided
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f) or {}
            config = deep_merge(config, user_config)
        except Exception as e:
            print(f"Error loading config file {config_path}: {e}")
            sys.exit(1)
    elif config_path:
        print(f"Error: Config file not found: {config_path}")
        sys.exit(1)
    
    return config


def merge_config_and_args(config: Dict[str, Any], args: argparse.Namespace) -> Dict[str, Any]:
    """
    Merge configuration with command-line arguments.
    CLI arguments take precedence over config file.
    """
    result = copy.deepcopy(config)
    
    # Ensure nested dictionaries exist
    if 'simulation' not in result:
        result['simulation'] = {}
    if 'output' not in result:
        result['output'] = {}
    if 'performance' not in result:
        result['performance'] = {}
    
    # Override with CLI arguments if provided
    if args.rate is not None:
        result['simulation']['rate'] = args.rate
        # Clear gene_rate_groups if rate is specified
        if 'gene_rate_groups' in result['simulation']:
            result['simulation']['gene_rate_groups'] = None
    
    if args.gene_rate_groups is not None:
        result['simulation']['gene_rate_groups'] = args.gene_rate_groups
        # Clear rate if gene_rate_groups is specified
        if 'rate' in result['simulation']:
            result['simulation']['rate'] = None
    
    if args.sites is not None:
        result['simulation']['sites'] = args.sites
    
    if args.years is not None:
        result['simulation']['years'] = args.years
    
    if args.gene_size is not None:
        result['simulation']['gene_size'] = args.gene_size
    
    if args.growth_phase is not None:
        result['simulation']['growth_phase'] = args.growth_phase
    
    if args.output is not None:
        result['output']['filename'] = args.output
    
    if args.seed is not None:
        result['seed'] = args.seed
    
    # Handle boolean flags (only override if explicitly set)
    if args.no_save:
        result['output']['save'] = False
    
    if args.no_compress:
        result['output']['compress'] = False
    
    if args.no_gene_jsd:
        result['performance']['track_gene_jsd'] = False
    
    if args.no_cell_jsds:
        result['performance']['calculate_cell_jsds'] = False
    
    if args.no_jsds:
        # Disable both cell and gene JSDs
        result['performance']['calculate_cell_jsds'] = False
        result['performance']['track_gene_jsd'] = False
    
    return result


def validate_config(config: Dict[str, Any]) -> None:
    """
    Validate configuration for consistency and required fields.
    
    Raises:
        ValueError: If configuration is invalid
    """
    sim = config.get('simulation', {})
    
    # Check for rate configuration
    has_rate = sim.get('rate') is not None
    has_gene_rates = sim.get('gene_rate_groups') is not None and sim.get('gene_rate_groups')
    
    if has_rate and has_gene_rates:
        raise ValueError("Cannot specify both 'rate' and 'gene_rate_groups'")
    
    # Validate growth phase
    growth_phase = sim.get('growth_phase', DEFAULT_GROWTH_PHASE)
    if growth_phase < 1 or growth_phase > 20:
        raise ValueError(f"growth_phase must be between 1 and 20, got {growth_phase}")
    
    # Validate sites and gene_size compatibility
    sites = sim.get('sites', N)
    gene_size = sim.get('gene_size', GENE_SIZE)
    if sites % gene_size != 0:
        raise ValueError(f"sites ({sites}) must be divisible by gene_size ({gene_size})")
    
    # Validate gene_rate_groups if specified
    if has_gene_rates:
        groups_str = sim['gene_rate_groups']
        try:
            groups = parse_gene_rate_groups(groups_str, sites, gene_size)
        except ValueError as e:
            raise ValueError(f"Invalid gene_rate_groups: {e}")


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run Step1-Prime methylation simulation with growth and homeostasis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Config file argument (highest priority)
    parser.add_argument('-c', '--config', type=str, default=None,
                        help='Path to configuration YAML file')
    
    # Simulation parameters - all now optional for config override
    parser.add_argument('-r', '--rate', type=float, default=None,
                        help='Uniform methylation rate per site per year')
    parser.add_argument('-n', '--sites', type=int, default=None,
                        help='Number of CpG sites per cell')
    parser.add_argument('-t', '--years', type=int, default=None,
                        help='Number of years to simulate')
    parser.add_argument('-g', '--gene-size', type=int, default=None,
                        help='Size of genes (sites per gene)')
    parser.add_argument('--growth-phase', type=int, default=None,
                        help='Duration of growth phase in years (cells double each year). '
                             'Final population = 2^growth_phase. Range: 1-20')
    parser.add_argument('--gene-rate-groups', type=str, default=None,
                        help='Gene-specific rates as "n1:rate1,n2:rate2,..." '
                             'Example: "50:0.004,50:0.0045,50:0.005,50:0.0055" '
                             'Total genes must equal n_sites/gene_size')
    
    # Optional parameters
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output filename (auto-generated if not specified)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility (use -1 for no seed)')
    parser.add_argument('--no-save', action='store_true',
                        help='Run simulation without saving results')
    
    # Performance and tracking options
    parser.add_argument('--no-gene-jsd', action='store_true',
                        help='Disable gene JSD tracking (population-level heterogeneity per gene)')
    parser.add_argument('--no-cell-jsds', action='store_true',
                        help='Disable cell JSD calculations (individual cell divergence from baseline)')
    parser.add_argument('--no-jsds', action='store_true',
                        help='Disable ALL JSD calculations (both cell and gene) for maximum performance')
    parser.add_argument('--no-compress', action='store_true',
                        help='Save output as uncompressed JSON instead of .json.gz (larger files but easier to inspect)')
    
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
    
    # Load configuration (defaults + user config + CLI overrides)
    config = load_config(args.config)
    config = merge_config_and_args(config, args)
    
    # Validate configuration
    try:
        validate_config(config)
    except ValueError as e:
        print(f"Error: {e}")
        return 1
    
    # Extract parameters from config
    sim_config = config.get('simulation', {})
    out_config = config.get('output', {})
    perf_config = config.get('performance', {})
    
    # Get simulation parameters with defaults
    n = sim_config.get('sites', N)
    t_max = sim_config.get('years', T_MAX)
    gene_size = sim_config.get('gene_size', GENE_SIZE)
    growth_phase = sim_config.get('growth_phase', DEFAULT_GROWTH_PHASE)
    
    # Get seed (handle -1 as None)
    seed = config.get('seed', 42)
    if seed == -1 or seed is None:
        seed = None
    
    # Get rate configuration
    rate = sim_config.get('rate')
    gene_rate_groups_str = sim_config.get('gene_rate_groups')
    
    # Parse gene rate groups if specified as string
    if gene_rate_groups_str:
        try:
            gene_rate_groups = parse_gene_rate_groups(
                gene_rate_groups_str, n, gene_size
            )
        except ValueError as e:
            print(f"Error: {e}")
            return 1
    # Convert rate to gene_rate_groups if specified
    elif rate is not None:
        gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
    else:
        gene_rate_groups = None  # Will use default in PetriDish
    
    # Start timing
    start_time = time.time()
    
    # Create and run simulation with appropriate rate configuration
    print("\nInitializing PetriDish simulation...")
    
    # Get performance flags
    calculate_cell_jsds = perf_config.get('calculate_cell_jsds', True)
    track_gene_jsd = perf_config.get('track_gene_jsd', True)
    
    # Create PetriDish with unified interface
    petri_dish = PetriDish(
        gene_rate_groups=gene_rate_groups,  # Always use this
        n=n,
        gene_size=gene_size,
        seed=seed,
        growth_phase=growth_phase,
        calculate_cell_jsds=calculate_cell_jsds
    )
    
    # Enable history tracking with gene JSD based on config
    if track_gene_jsd:
        petri_dish.enable_history_tracking(track_gene_jsd=True)
        print("Gene JSD tracking enabled")
    else:
        petri_dish.enable_history_tracking(track_gene_jsd=False)
        print("Gene JSD tracking disabled")
    
    # Run simulation
    petri_dish.run_simulation(t_max=t_max)
    
    # Calculate runtime
    runtime = time.time() - start_time
    print(f"\nSimulation runtime: {runtime:.2f} seconds")
    
    # Save results based on config
    save_results = out_config.get('save', True)
    compress = out_config.get('compress', True)
    output_filename = out_config.get('filename')
    
    if save_results:
        filepath = petri_dish.save_history(filename=output_filename, compress=compress)
        print(f"\nResults saved to: {filepath}")
    else:
        print("\nResults not saved (save disabled in config)")
    
    # Print summary statistics
    print("\nFinal Statistics:")
    print(f"  Total cells: {len(petri_dish.cells)}")
    
    # Calculate JSD statistics
    jsd_values = [cell.cell_jsd for cell in petri_dish.cells]
    if jsd_values:
        import statistics
        print(f"  Mean JSD: {statistics.mean(jsd_values):.4f}")
        print(f"  Median JSD: {statistics.median(jsd_values):.4f}")
        print(f"  Std Dev JSD: {statistics.stdev(jsd_values):.4f}")
        print(f"  Min JSD: {min(jsd_values):.4f}")
        print(f"  Max JSD: {max(jsd_values):.4f}")
    
    # Calculate methylation statistics
    meth_props = [cell.cell_methylation_proportion for cell in petri_dish.cells]
    if meth_props:
        print(f"  Mean cell methylation proportion: {statistics.mean(meth_props):.2%}")
        print(f"  Median cell methylation proportion: {statistics.median(meth_props):.2%}")
    
    print("\nSimulation complete!")
    return 0


if __name__ == "__main__":
    exit(main())