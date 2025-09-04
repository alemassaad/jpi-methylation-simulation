#!/usr/bin/env python3
"""
Phase 2 Pipeline: Refactored to use PetriDish and Cell classes directly.

This pipeline:
1. Extracts year 50 and 60 snapshots as Cell objects
2. Creates individuals using PetriDish objects
3. Grows individuals using PetriDish division/methylation
4. Mixes populations
5. Analyzes results

Key improvement: No dictionary conversions during processing!
"""

import argparse
import os
import sys
import time
import json
import gzip
import numpy as np
import math
from typing import List, Dict, Optional, Tuple, Any

# Try to import YAML support
try:
    import yaml
except ImportError:
    yaml = None

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish, Cell, PetriDishPlotter, rate_to_gene_rate_groups
from core.individual_helpers import (
    create_individual, process_batch_growth, process_batch_mixing
)
from core.pipeline_utils import (
    load_snapshot_as_cells, save_snapshot_cells, load_snapshot_cells,
    save_petri_dish, load_petri_dish, load_all_petri_dishes,
    sample_by_quantiles, sample_uniform,
    create_pure_snapshot_petri,
    create_control2_with_uniform_base,
    check_petri_files_state,
    calculate_population_statistics, print_mixing_statistics,
    create_uniform_mixing_pool, mix_petri_uniform, normalize_populations,
    normalize_individuals_for_uniform_mixing
)
from core.pipeline_analysis import (
    plot_cell_jsd_distribution,
    analyze_populations_from_dishes,
    plot_gene_jsd_distribution_comparison,
    plot_top_variable_genes,
    plot_gene_jsd_distributions,
    plot_gene_jsd_individual_comparison
)
from core.path_utils import parse_step1_simulation_path, generate_step23_output_dir
# Note: validate_timeline_parameters, check_timeline_warnings, print_timeline_breakdown are defined locally below


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
    default_path = os.path.join(os.path.dirname(__file__), "configs", "config_default.yaml")
    if os.path.exists(default_path):
        try:
            with open(default_path, 'r') as f:
                default_config = yaml.safe_load(f) or {}
                config.update(default_config)
        except Exception as e:
            print(f"Warning: Could not load default config: {e}")
    
    # Load user config if provided
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f) or {}
                # Deep merge user config into default config
                deep_merge(config, user_config)
        except Exception as e:
            print(f"Error loading config file {config_path}: {e}")
            sys.exit(1)
    
    return config


def deep_merge(base: dict, update: dict) -> None:
    """
    Deep merge update dictionary into base dictionary.
    
    Args:
        base: Base dictionary to update
        update: Dictionary with updates
    """
    for key, value in update.items():
        if key in base and isinstance(base[key], dict) and isinstance(value, dict):
            deep_merge(base[key], value)
        else:
            base[key] = value


def merge_config_and_args(config: Dict[str, Any], args: argparse.Namespace) -> argparse.Namespace:
    """
    Merge config file values with command-line arguments.
    CLI arguments take precedence over config file values.
    
    Args:
        config: Configuration dictionary from YAML
        args: Command-line arguments
    
    Returns:
        Updated args namespace
    """
    # Map config structure to args
    # Priority: CLI args > config file > defaults
    
    # Input configuration
    if 'input' in config:
        if not args.simulation and config['input'].get('simulation'):
            args.simulation = config['input']['simulation']
        if not args.rate and not args.gene_rate_groups:
            if config['input'].get('rate'):
                args.rate = config['input']['rate']
            elif config['input'].get('gene_rate_groups'):
                args.gene_rate_groups = config['input']['gene_rate_groups']
        if hasattr(args, 'gene_size') and args.gene_size == 5:  # default value
            args.gene_size = config['input'].get('gene_size', 5)
    
    # Snapshots configuration
    if 'snapshots' in config:
        if args.first_snapshot == 50:  # default value
            args.first_snapshot = config['snapshots'].get('first', 50)
        if args.second_snapshot == 60:  # default value
            args.second_snapshot = config['snapshots'].get('second', 60)
        if hasattr(args, 'force_reload') and not args.force_reload:
            args.force_reload = config['snapshots'].get('force_reload', False)
    
    # Individuals configuration
    if 'individuals' in config:
        if args.individual_growth_phase == 7:  # default value
            args.individual_growth_phase = config['individuals'].get('growth_phase', 7)
        if args.n_quantiles == 10:  # default value
            args.n_quantiles = config['individuals'].get('n_quantiles', 10)
        if args.cells_per_quantile == 3:  # default value
            args.cells_per_quantile = config['individuals'].get('cells_per_quantile', 3)
    
    # Mixing configuration
    if 'mixing' in config:
        if args.mix_ratio == 80:  # default value
            args.mix_ratio = config['mixing'].get('ratio', 80)
        if not args.uniform_mixing:
            args.uniform_mixing = config['mixing'].get('uniform', False)
        if not args.normalize_size:
            args.normalize_size = config['mixing'].get('normalize_size', False)
    
    # Visualization configuration
    if 'visualization' in config:
        if args.bins == 200:  # default value
            args.bins = config['visualization'].get('bins', 200)
        if not args.plot_individuals:
            args.plot_individuals = config['visualization'].get('plot_individuals', False)
    
    # Output configuration
    if 'output' in config:
        if hasattr(args, 'output_dir'):
            args.output_dir = config['output'].get('directory', args.output_dir)
        if hasattr(args, 'no_compress') and not args.no_compress:
            # Config uses 'compress: true/false', CLI uses '--no-compress'
            compress = config['output'].get('compress', True)
            args.no_compress = not compress
    
    # Other settings
    if args.seed == 42:  # default value
        args.seed = config.get('seed', 42)
    
    if 'verbose' in config and not hasattr(args, 'verbose'):
        args.verbose = config.get('verbose', False)
    
    return args


def validate_pipeline_config(args: argparse.Namespace) -> None:
    """
    Validate the final merged configuration.
    
    Args:
        args: Merged configuration
    
    Raises:
        ValueError: If configuration is invalid
    """
    # Required fields
    if not args.simulation:
        raise ValueError("--simulation is required (or set in config file)")
    
    # Rate specification is now optional - can be inferred from simulation
    if args.rate and args.gene_rate_groups:
        raise ValueError("Cannot specify both --rate and --gene-rate-groups")
    
    # Value ranges
    if args.mix_ratio < 0 or args.mix_ratio > 100:
        raise ValueError(f"--mix-ratio must be between 0 and 100, got {args.mix_ratio}")
    
    if args.n_quantiles < 1:
        raise ValueError(f"--n-quantiles must be >= 1, got {args.n_quantiles}")
    
    if args.cells_per_quantile < 1:
        raise ValueError(f"--cells-per-quantile must be >= 1, got {args.cells_per_quantile}")
    
    if args.bins < 10:
        raise ValueError(f"--bins must be >= 10, got {args.bins}")
    
    # Timeline validation is done elsewhere
    if args.second_snapshot <= args.first_snapshot:
        raise ValueError(f"--second-snapshot ({args.second_snapshot}) must be > --first-snapshot ({args.first_snapshot})")


def parse_gene_rate_groups(gene_rate_str: str) -> List[Tuple[int, float]]:
    """
    Parse gene rate groups string into list of tuples.
    
    Args:
        gene_rate_str: String like "50:0.004,50:0.0045,50:0.005,50:0.0055"
    
    Returns:
        List of (n_genes, rate) tuples
    """
    groups = []
    for group in gene_rate_str.split(','):
        n_genes, rate = group.strip().split(':')
        groups.append((int(n_genes), float(rate)))
    return groups


def validate_gene_rate_groups(groups: List[Tuple[int, float]], n_sites: int, gene_size: int) -> None:
    """
    Validate gene rate groups configuration.
    
    Args:
        groups: List of (n_genes, rate) tuples
        n_sites: Total number of CpG sites
        gene_size: Sites per gene
    
    Raises:
        ValueError: If configuration is invalid
    """
    total_genes = sum(n for n, _ in groups)
    expected_genes = n_sites // gene_size
    
    if total_genes != expected_genes:
        raise ValueError(
            f"Gene rate groups total {total_genes} genes, but need {expected_genes} "
            f"({n_sites} sites / {gene_size} sites per gene)"
        )
    
    for n_genes, rate in groups:
        if n_genes <= 0:
            raise ValueError(f"Number of genes must be positive, got {n_genes}")
        if rate < 0 or rate > 1:
            raise ValueError(f"Rate must be between 0 and 1, got {rate}")


def validate_all_gene_rate_groups(
    expected: List[Tuple[int, float]], 
    petri_dishes: List[PetriDish] = None,
    cells: List[Cell] = None,
    label: str = ""
) -> None:
    """Validate all objects have expected gene_rate_groups.
    
    Args:
        expected: Expected gene rate groups
        petri_dishes: Optional list of PetriDish objects to validate
        cells: Optional list of Cell objects to validate
        label: Description for error messages
    """
    if petri_dishes:
        for i, petri in enumerate(petri_dishes):
            # Check PetriDish itself
            if petri.gene_rate_groups != expected:
                raise ValueError(
                    f"{label} PetriDish {i} has wrong gene_rate_groups!\n"
                    f"  Expected: {expected}\n"
                    f"  Found: {petri.gene_rate_groups}"
                )
            
            # Check all cells in the dish
            for j, cell in enumerate(petri.cells):
                if cell.gene_rate_groups != expected:
                    raise ValueError(
                        f"{label} PetriDish {i}, Cell {j} has wrong gene_rate_groups!\n"
                        f"  Expected: {expected}\n"
                        f"  Found: {cell.gene_rate_groups}"
                    )
    
    if cells:
        for i, cell in enumerate(cells):
            if cell.gene_rate_groups != expected:
                raise ValueError(
                    f"{label} Cell {i} has wrong gene_rate_groups!\n"
                    f"  Expected: {expected}\n"
                    f"  Found: {cell.gene_rate_groups}"
                )
    
    print(f"✓ {label}: All gene_rate_groups verified as {expected}")


def calculate_timeline_metrics(first_snapshot: int, second_snapshot: int, individual_growth_phase: int) -> dict:
    """Calculate all timeline-related metrics."""
    timeline_duration = second_snapshot - first_snapshot
    exponential_years = individual_growth_phase
    homeostasis_years = timeline_duration - exponential_years
    target_cells = 2 ** exponential_years
    
    return {
        'timeline_duration': timeline_duration,
        'exponential_years': exponential_years,
        'homeostasis_years': homeostasis_years,
        'target_cells': target_cells
    }


def create_timeline_error_message(first_snapshot: int, second_snapshot: int, individual_growth_phase: int) -> str:
    """Create detailed error message for timeline validation failures."""
    timeline_duration = second_snapshot - first_snapshot
    target_cells = 2 ** individual_growth_phase
    max_allowed_growth = timeline_duration
    min_second_snapshot = first_snapshot + individual_growth_phase + 1
    
    return f"""Timeline validation failed:
  Individual growth phase: {individual_growth_phase} years (exponential growth to 2^{individual_growth_phase} = {target_cells:,} cells)
  Available time: {timeline_duration} years (year {first_snapshot} → year {second_snapshot})
  Problem: Need time for both exponential growth AND homeostasis phases
  
  Fix option 1: Reduce --individual-growth-phase to {max_allowed_growth} or less
  Fix option 2: Increase --second-snapshot to {min_second_snapshot} or higher
  
  Example working commands:
    # Option 1: Smaller populations
    --individual-growth-phase {max_allowed_growth} --first-snapshot {first_snapshot} --second-snapshot {second_snapshot}
    
    # Option 2: Longer timeline  
    --individual-growth-phase {individual_growth_phase} --first-snapshot {first_snapshot} --second-snapshot {min_second_snapshot}"""


def check_timeline_warnings(metrics: dict, individual_growth_phase: int) -> List[str]:
    """Check for edge cases and return warning messages."""
    warnings = []
    timeline_duration = metrics['timeline_duration']
    homeostasis_years = metrics['homeostasis_years']
    target_cells = metrics['target_cells']
    
    if homeostasis_years == 0:
        warnings.append(f"⚠️  No homeostasis phase: All {timeline_duration} years used for exponential growth")
        warnings.append("   Population may be unstable. Consider adding 1-2 more years to second-snapshot.")
    
    elif homeostasis_years == 1:
        warnings.append(f"⚠️  Only 1 year for homeostasis phase")
        warnings.append("   Consider more time for population stability.")
    
    if individual_growth_phase > 10:
        warnings.append(f"⚠️  Very large target population: {target_cells:,} cells")
        warnings.append("   This may be slow and memory-intensive.")
    
    if homeostasis_years > timeline_duration * 0.8:
        warnings.append(f"⚠️  Very long homeostasis phase: {homeostasis_years}/{timeline_duration} years")
        warnings.append("   Most time spent in steady state, minimal exponential growth.")
    
    return warnings


def validate_timeline_parameters(first_snapshot: int, second_snapshot: int, individual_growth_phase: int) -> dict:
    """Validate timeline parameters and raise detailed errors if invalid."""
    # Basic validation
    if second_snapshot <= first_snapshot:
        raise ValueError(f"second-snapshot ({second_snapshot}) must be > first-snapshot ({first_snapshot})")
    
    if individual_growth_phase < 1 or individual_growth_phase > 15:
        raise ValueError(f"individual-growth-phase must be between 1 and 15, got {individual_growth_phase}")
    
    # Calculate metrics
    metrics = calculate_timeline_metrics(first_snapshot, second_snapshot, individual_growth_phase)
    
    # Timeline validation
    if individual_growth_phase > metrics['timeline_duration']:
        error_message = create_timeline_error_message(first_snapshot, second_snapshot, individual_growth_phase)
        raise ValueError(error_message)
    
    return metrics


def print_timeline_breakdown(metrics: dict, first_snapshot: int, second_snapshot: int, individual_growth_phase: int):
    """Display visual timeline breakdown."""
    timeline_duration = metrics['timeline_duration']
    exponential_years = metrics['exponential_years']
    homeostasis_years = metrics['homeostasis_years']
    target_cells = metrics['target_cells']
    
    print(f"\n{'='*60}")
    print("TIMELINE BREAKDOWN")
    print(f"{'='*60}")
    print(f"📊 Individual Simulation Timeline:")
    print(f"   Year {first_snapshot}: Sample individual cells")
    print(f"   Years {first_snapshot}→{first_snapshot + exponential_years}: Exponential growth (1 → {target_cells:,} cells)")
    if homeostasis_years > 0:
        print(f"   Years {first_snapshot + exponential_years}→{second_snapshot}: Homeostasis (~{target_cells:,} cells)")
    print(f"   Year {second_snapshot}: Extract for mixing")
    print(f"")
    print(f"📈 Growth Summary:")
    print(f"   Total aging time: {timeline_duration} years")
    print(f"   Exponential phase: {exponential_years} years ({exponential_years/timeline_duration*100:.1f}%)")
    print(f"   Homeostasis phase: {homeostasis_years} years ({homeostasis_years/timeline_duration*100:.1f}%)")
    print(f"   Target population: {target_cells:,} cells (2^{exponential_years})")


def run_pipeline(args):
    """
    Main pipeline orchestrator using PetriDish objects.
    Now uses gene_rate_groups as the single internal representation.
    
    Args:
        args: Command-line arguments
    """
    # Calculate derived values and validate using new helper functions
    metrics = validate_timeline_parameters(args.first_snapshot, args.second_snapshot, args.individual_growth_phase)
    
    # Extract metrics for backward compatibility
    timeline_duration = metrics['timeline_duration']
    homeostasis_years = metrics['homeostasis_years']
    expected_population = metrics['target_cells']
    
    # Check for warnings and display if any
    warnings = check_timeline_warnings(metrics, args.individual_growth_phase)
    if warnings:
        print(f"\n⚠️  TIMELINE WARNINGS:")
        for warning in warnings:
            print(warning)
    
    # Display timeline breakdown
    print_timeline_breakdown(metrics, args.first_snapshot, args.second_snapshot, args.individual_growth_phase)
    
    start_time = time.time()
    
    # Set global random seed for reproducibility
    import random
    import numpy as np
    random.seed(args.seed)
    np.random.seed(args.seed)
    print(f"Random seed set to {args.seed} for reproducibility")
    
    # ========================================================================
    # SETUP
    # ========================================================================
    
    # Determine compression strategy once, based on priority:
    # 1. Explicit CLI flag (--no-compress)
    # 2. Input file format (.json vs .json.gz)
    # 3. Config file setting (already applied)
    # 4. Default (compress)
    
    if hasattr(args, 'no_compress') and args.no_compress:
        # CLI flag explicitly says no compression
        args.use_compression = False
        print(f"Compression strategy: uncompressed (.json) [from --no-compress flag]")
    elif args.simulation.endswith('.gz'):
        # Match input format - compressed
        args.use_compression = True
        print(f"Compression strategy: compressed (.json.gz) [matching input format]")
    elif args.simulation.endswith('.json'):
        # Match input format - uncompressed
        args.use_compression = False
        print(f"Compression strategy: uncompressed (.json) [matching input format]")
    else:
        # Default to compressed
        args.use_compression = True
        print(f"Compression strategy: compressed (.json.gz) [default]")
    
    # Load simulation and get its parameters
    print("Loading simulation parameters...")
    with gzip.open(args.simulation, 'rt') if args.simulation.endswith('.gz') else open(args.simulation, 'r') as f:
        sim_data = json.load(f)
    
    sim_params_raw = sim_data['parameters']
    n_sites = sim_params_raw.get('n', 1000)
    gene_size = sim_params_raw.get('gene_size', 5)
    
    # Get gene_rate_groups from simulation
    sim_gene_rate_groups = sim_params_raw.get('gene_rate_groups')
    if not sim_gene_rate_groups:
        # Convert rate to gene_rate_groups for old simulations
        if 'rate' in sim_params_raw:
            sim_gene_rate_groups = rate_to_gene_rate_groups(
                sim_params_raw['rate'], n_sites, gene_size
            )
        else:
            raise ValueError("Simulation has no rate configuration!")
    
    # Convert to tuples
    sim_gene_rate_groups = [tuple(group) for group in sim_gene_rate_groups]
    
    # Get gene_rate_groups - either from arguments or infer from simulation
    if args.gene_rate_groups:
        args_gene_rate_groups = parse_gene_rate_groups(args.gene_rate_groups)
        # Validate they're valid for the n_sites
        validate_gene_rate_groups(args_gene_rate_groups, n_sites, gene_size)
        
        # Verify simulation matches arguments
        if sim_gene_rate_groups != args_gene_rate_groups:
            raise ValueError(
                f"FATAL: Rate configuration mismatch!\n"
                f"  Simulation has: gene_rate_groups={sim_gene_rate_groups}\n"
                f"  You specified: gene_rate_groups={args_gene_rate_groups}\n"
                f"  These must match exactly.\n"
                f"  Please use the same rate configuration as the simulation."
            )
    elif args.rate:
        args_gene_rate_groups = rate_to_gene_rate_groups(args.rate, n_sites, gene_size)
        
        # Verify simulation matches arguments
        if sim_gene_rate_groups != args_gene_rate_groups:
            raise ValueError(
                f"FATAL: Rate configuration mismatch!\n"
                f"  Simulation has: gene_rate_groups={sim_gene_rate_groups}\n"
                f"  You specified: rate={args.rate} (converted to {args_gene_rate_groups})\n"
                f"  These must match exactly.\n"
                f"  Please use the same rate configuration as the simulation."
            )
    else:
        # No rate specified - infer from simulation
        print(f"  ℹ No rate specified, inferring from simulation...")
        args_gene_rate_groups = sim_gene_rate_groups
    
    # This is THE gene_rate_groups for the entire pipeline
    gene_rate_groups = sim_gene_rate_groups
    print(f"✓ Pipeline gene_rate_groups: {gene_rate_groups}")
    
    # Parse simulation parameters from input file (for path generation)
    sim_params = parse_step1_simulation_path(args.simulation)
    if not sim_params:
        print(f"Error: Could not parse simulation parameters from {args.simulation}")
        sys.exit(1)
    
    # Calculate total expected individuals
    expected_individuals = args.n_quantiles * args.cells_per_quantile
    
    # Generate hierarchical output directory
    base_dir = generate_step23_output_dir(args, sim_params)
    
    snapshots_dir = os.path.join(base_dir, "snapshots")
    individuals_dir = os.path.join(base_dir, "individuals")
    mutant_dir = os.path.join(individuals_dir, "mutant")
    control1_dir = os.path.join(individuals_dir, "control1")
    control2_dir = os.path.join(individuals_dir, "control2")
    results_dir = os.path.join(base_dir, "results")
    
    # Create all directories (removed plots_dir)
    for dir_path in [snapshots_dir, mutant_dir, control1_dir, control2_dir, results_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    print("=" * 80)
    print("PHASE 2 PIPELINE")
    print("=" * 80)
    print(f"Gene rate groups: {gene_rate_groups}")
    print(f"Simulation: {args.simulation}")
    print(f"Output: {base_dir}")
    print(f"Quantiles: {args.n_quantiles} ({'quartiles' if args.n_quantiles == 4 else 'deciles' if args.n_quantiles == 10 else f'{args.n_quantiles}-tiles'})")
    print(f"Cells per quantile: {args.cells_per_quantile}")
    print(f"Total individuals: {expected_individuals} per group")
    print(f"Snapshots: year {args.first_snapshot} → year {args.second_snapshot}")
    print(f"Growth: {timeline_duration} years total ({args.individual_growth_phase} growth + {homeostasis_years} homeostasis)")
    print(f"Target population: {expected_population} cells (2^{args.individual_growth_phase})")
    print(f"Mix ratio: {args.mix_ratio}% year {args.second_snapshot}, {100-args.mix_ratio}% grown")
    print(f"Seed: {args.seed}")
    print("=" * 80)
    
    # ========================================================================
    # STAGE 1: Extract First Snapshot
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 1: Extract Year {args.first_snapshot} Snapshot")
    print(f"{'='*60}")
    
    # File extensions based on compression setting
    snapshot_ext = ".json.gz" if args.use_compression else ".json"
    individual_ext = ".json.gz" if args.use_compression else ".json"
    first_snapshot_path = os.path.join(snapshots_dir, f"year{args.first_snapshot}_snapshot{snapshot_ext}")
    
    if os.path.exists(first_snapshot_path) and not args.force_reload:
        print(f"  ⏭ Loading existing snapshot from {first_snapshot_path}...")
        first_snapshot_cells = load_snapshot_cells(first_snapshot_path)
        print(f"  Loaded {len(first_snapshot_cells)} Cell objects")
    else:
        print(f"  ✓ Extracting year {args.first_snapshot} from simulation...")
        first_snapshot_cells = load_snapshot_as_cells(
            args.simulation, args.first_snapshot, 
            expected_gene_rate_groups=gene_rate_groups
        )
        save_snapshot_cells(first_snapshot_cells, first_snapshot_path, compress=args.use_compression)
        print(f"  Cached {len(first_snapshot_cells)} cells for future use")
    
    # ========================================================================
    # VALIDATION: Snapshot cells
    # ========================================================================
    from core.validation import PipelineValidator, ValidationError
    validator = PipelineValidator(verbose=not args.quiet if hasattr(args, 'quiet') else True)
    
    try:
        min_cells_needed = args.n_quantiles * args.cells_per_quantile * 2  # Need enough for sampling
        validator.validate_snapshot(
            cells=first_snapshot_cells,
            expected_year=args.first_snapshot,
            expected_gene_rate_groups=gene_rate_groups,
            min_cells=min_cells_needed
        )
    except ValidationError as e:
        print(f"\n❌ Snapshot validation failed: {e}")
        sys.exit(1)
    
    # ========================================================================
    # STAGE 2: Plot JSD Distribution
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 2: Plot JSD Distribution")
    print(f"{'='*60}")
    
    plot_path = os.path.join(results_dir, f"year{args.first_snapshot}_jsd_distribution_{args.bins}bins.png")
    # Pass gene_rate_groups for plot labeling
    plot_cell_jsd_distribution(first_snapshot_cells, args.bins, plot_path, 
                                    gene_rate_groups=gene_rate_groups, year=args.first_snapshot)
    
    # ========================================================================
    # STAGE 3: Create Initial Individuals (as PetriDish objects)
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 3: Create Initial Individuals")
    print(f"{'='*60}")
    
    # Create mutant individuals in memory
    mutant_dishes = []
    print(f"\n  ✓ Creating {expected_individuals} mutant individuals...")
    print(f"  Sampling by quantiles ({args.n_quantiles} quantiles, {args.cells_per_quantile} each)...")
    
    sampled = sample_by_quantiles(first_snapshot_cells, 
                                 n_quantiles=args.n_quantiles,
                                 cells_per_quantile=args.cells_per_quantile,
                                 seed=args.seed)
    
    for i, (cell, quantile) in enumerate(sampled):
        file_index = i + 1  # Start IDs at 1 instead of 0
        
        # Use helper function for cleaner code
        additional_metadata = {
            'source_quantile': quantile,
            'n_quantiles': args.n_quantiles
        }
        
        petri = create_individual(
            cell=cell,
            individual_type='mutant',
            individual_id=file_index,
            growth_phase=args.individual_growth_phase,
            additional_metadata=additional_metadata
        )
        
        mutant_dishes.append(petri)
    
    print(f"  Created {len(mutant_dishes)} mutant individuals (in memory)")
    
    # Create control1 individuals in memory
    control1_dishes = []
    print(f"\n  ✓ Creating {expected_individuals} control1 individuals...")
    print(f"  Sampling uniformly...")
    
    sampled_cells = sample_uniform(first_snapshot_cells, n_samples=expected_individuals, 
                                  seed=args.seed + 1000)
    
    for i, cell in enumerate(sampled_cells):
        file_index = i + 1  # Start IDs at 1 instead of 0
        
        # Use helper function for cleaner code
        additional_metadata = {
            'source': 'uniform'
        }
        
        petri = create_individual(
            cell=cell,
            individual_type='control1',
            individual_id=file_index,
            growth_phase=args.individual_growth_phase,
            additional_metadata=additional_metadata
        )
        
        control1_dishes.append(petri)
    
    print(f"  Created {len(control1_dishes)} control1 individuals (in memory)")
    
    # Validate all initial individuals have correct gene_rate_groups
    validate_all_gene_rate_groups(
        expected=gene_rate_groups,
        petri_dishes=mutant_dishes,
        label="Mutant individuals"
    )
    validate_all_gene_rate_groups(
        expected=gene_rate_groups,
        petri_dishes=control1_dishes,
        label="Control1 individuals"
    )
    
    # ========================================================================
    # VALIDATION: Initial individuals
    # ========================================================================
    try:
        validator.validate_initial_individuals(
            mutant_dishes=mutant_dishes,
            control1_dishes=control1_dishes,
            expected_count=expected_individuals,
            expected_gene_rate_groups=gene_rate_groups,
            snapshot_year=args.first_snapshot
        )
    except ValidationError as e:
        print(f"\n❌ Initial individuals validation failed: {e}")
        sys.exit(1)
    
    # ========================================================================
    # STAGE 4: Grow Individuals (using PetriDish methods with homeostasis)
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 4: Grow Individuals ({timeline_duration} years)")
    print(f"{'='*60}")
    
    # Grow mutant individuals in memory
    process_batch_growth(
        dishes=mutant_dishes,
        batch_name='mutant',
        years=timeline_duration,
        growth_phase=args.individual_growth_phase,
        expected_population=expected_population,
        start_year=args.first_snapshot,
        verbose=True
    )
    
    # Grow control1 individuals in memory
    process_batch_growth(
        dishes=control1_dishes,
        batch_name='control1',
        years=timeline_duration,
        growth_phase=args.individual_growth_phase,
        expected_population=expected_population,
        start_year=args.first_snapshot,
        verbose=True
    )
    
    # ========================================================================
    # VALIDATION: Grown individuals
    # ========================================================================
    extinction_report = validator.validate_grown_individuals(
        mutant_dishes=mutant_dishes,
        control1_dishes=control1_dishes,
        expected_population=expected_population,
        growth_years=timeline_duration,
        allow_extinction=True  # Allow but track extinction
    )
    
    if extinction_report['total_extinct'] > 0:
        print(f"\n  ⚠️  Warning: {extinction_report['total_extinct']} individuals went extinct during growth")
        if extinction_report['complete_batch_extinction']:
            print("\n❌ Fatal: Entire batch extinct. Consider different parameters.")
            sys.exit(1)
    
    # ========================================================================
    # STAGE 5: Extract Second Snapshot
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 5: Extract Year {args.second_snapshot} Snapshot")
    print(f"{'='*60}")
    
    second_snapshot_path = os.path.join(snapshots_dir, f"year{args.second_snapshot}_snapshot{snapshot_ext}")
    
    if os.path.exists(second_snapshot_path) and not args.force_reload:
        print(f"  ⏭ Loading existing snapshot from {second_snapshot_path}...")
        second_snapshot_cells = load_snapshot_cells(second_snapshot_path)
        print(f"  Loaded {len(second_snapshot_cells)} Cell objects")
    else:
        print(f"  ✓ Extracting year {args.second_snapshot} from simulation...")
        second_snapshot_cells = load_snapshot_as_cells(
            args.simulation, args.second_snapshot,
            expected_gene_rate_groups=gene_rate_groups
        )
        save_snapshot_cells(second_snapshot_cells, second_snapshot_path, compress=args.use_compression)
        print(f"  Cached {len(second_snapshot_cells)} cells for future use")
    
    # ========================================================================
    # STAGE 5.5: Plot Second Snapshot JSD Distribution
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 5.5: Plot Year {args.second_snapshot} JSD Distribution")
    print(f"{'='*60}")
    
    second_plot_path = os.path.join(results_dir, f"year{args.second_snapshot}_jsd_distribution_{args.bins}bins.png")
    # Pass gene_rate_groups for plot labeling
    plot_cell_jsd_distribution(second_snapshot_cells, args.bins, second_plot_path, 
                                    gene_rate_groups=gene_rate_groups, year=args.second_snapshot)
    
    # ========================================================================
    # STAGE 6: Mix Populations
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 6: Mix Populations")
    print(f"{'='*60}")
    
    # Use in-memory dishes from previous stages (already have history)
    
    # Validate compatibility before mixing
    print("  Validating cell compatibility for mixing...")
    if mutant_dishes and second_snapshot_cells:
        sample_grown = mutant_dishes[0].cells[:5]
        sample_snapshot = second_snapshot_cells[:5]
        PetriDish.validate_cells_compatible(sample_grown + sample_snapshot)
        print("  ✓ Grown and snapshot cells are compatible")
    
    # Apply normalization if requested
    normalization_threshold = None
    if args.normalize_size:
        print("\n  === APPLYING SIZE NORMALIZATION ===")
        print("  Using median - 0.5σ threshold")
        
        # Normalize populations to same size (in memory)
        mutant_dishes, control1_dishes, normalization_threshold = normalize_populations(
            mutant_dishes, control1_dishes, 
            seed=args.seed + 5000
        )
        
        # Update metadata but don't save yet
        for dish in mutant_dishes:
            if not hasattr(dish, 'metadata'):
                dish.metadata = {}
            dish.metadata['normalized'] = True
            dish.metadata['normalization_threshold'] = normalization_threshold
        
        for dish in control1_dishes:
            if not hasattr(dish, 'metadata'):
                dish.metadata = {}
            dish.metadata['normalized'] = True
            dish.metadata['normalization_threshold'] = normalization_threshold
        
        print(f"  Normalized all individuals to {normalization_threshold} cells (in memory)")
        
        # Validate normalization
        try:
            validator.validate_normalized_populations(
                mutant_after=mutant_dishes,
                control1_after=control1_dishes,
                threshold=normalization_threshold
            )
        except ValidationError as e:
            print(f"\n❌ Normalization validation failed: {e}")
            sys.exit(1)
    
    if args.uniform_mixing:
        print("\n  === UNIFORM MIXING MODE (Fixed Normalization) ===")
        
        # Step 1: Normalize all individuals to same size first
        print(f"  Step 1: Normalizing individuals for uniform mixing...")
        normalized_mutant, normalized_control1, normalized_size = normalize_individuals_for_uniform_mixing(
            mutant_dishes, control1_dishes, seed=args.seed + 999
        )
        
        # Update the dish lists to use normalized versions
        mutant_dishes = normalized_mutant
        control1_dishes = normalized_control1
        
        # Step 2: Create uniform pool based on normalized size
        print(f"\n  Step 2: Creating uniform mixing pool from year {args.second_snapshot}...")
        uniform_pool, uniform_indices = create_uniform_mixing_pool(
            second_snapshot_cells,
            normalized_size,  # Now using the exact normalized size
            args.mix_ratio / 100,
            seed=args.seed + 1000
        )
        
        # Calculate statistics after normalization
        mutant_stats = calculate_population_statistics(mutant_dishes, "Mutant (normalized)")
        control1_stats = calculate_population_statistics(control1_dishes, "Control1 (normalized)")
        print_mixing_statistics(mutant_stats, control1_stats, normalized_size)
        
        # Save statistics to file
        stats_path = os.path.join(results_dir, "mixing_statistics.json")
        mixing_stats = {
            'mutant': mutant_stats,
            'control1': control1_stats,
            'combined_median': normalized_size,
            'uniform_mixing': True,
            'normalization_applied': True
        }
        with open(stats_path, 'w') as f:
            json.dump(mixing_stats, f, indent=2, default=str)
        print(f"\n  Saved mixing statistics to {stats_path}")
        
        # Step 3: Mix all individuals (simplified since they're all the same size)
        print(f"\n  Step 3: Mixing normalized individuals with uniform pool...")
        
        # Mix mutant individuals
        print(f"  Processing {len(mutant_dishes)} mutant individuals...")
        for i, petri in enumerate(mutant_dishes):
            # Ensure metadata exists before accessing it
            if not hasattr(petri, 'metadata'):
                petri.metadata = {}
            
            # Use individual_id from metadata for consistent file naming
            individual_id = petri.metadata.get('individual_id', i + 1)
            initial_size = len(petri.cells)
            final_size = mix_petri_uniform(petri, uniform_pool, args.mix_ratio / 100)
            print(f"    Mutant {individual_id:02d}: {initial_size} → {final_size} cells")
            
            # Update metadata
            petri.metadata['mixed'] = True
            petri.metadata['mix_mode'] = 'uniform'
            petri.metadata['mix_ratio'] = args.mix_ratio
            petri.metadata['normalized'] = True
            petri.metadata['final_cells'] = final_size
            # Don't save yet - will save all at the end
        
        # Mix control1 individuals (using SAME pool)
        print(f"  Processing {len(control1_dishes)} control1 individuals...")
        for i, petri in enumerate(control1_dishes):
            # Ensure metadata exists before accessing it
            if not hasattr(petri, 'metadata'):
                petri.metadata = {}
            
            # Use individual_id from metadata for consistent file naming
            individual_id = petri.metadata.get('individual_id', i + 1)
            initial_size = len(petri.cells)
            final_size = mix_petri_uniform(petri, uniform_pool, args.mix_ratio / 100)
            print(f"    Control1 {individual_id:02d}: {initial_size} → {final_size} cells")
            
            # Update metadata
            petri.metadata['mixed'] = True
            petri.metadata['mix_mode'] = 'uniform'
            petri.metadata['mix_ratio'] = args.mix_ratio
            petri.metadata['normalized'] = True
            petri.metadata['final_cells'] = final_size
            # Don't save yet - will save all at the end
                
    else:
        print("\n  === INDEPENDENT MIXING MODE (default) ===")
        
        # Calculate expected final size after mixing
        # If grown cells are 20% of final, what's 100%?
        grown_cells = expected_population  # e.g., 128 cells
        grown_fraction = (100 - args.mix_ratio) / 100  # e.g., 0.2 for 20% grown
        expected_final_cells = int(grown_cells / grown_fraction)
        
        print(f"  Target size after mixing: {expected_final_cells} cells")
        print(f"  Mix ratio: {args.mix_ratio}% from year {args.second_snapshot}")
        
        # Mix mutant individuals in memory
        process_batch_mixing(
            dishes=mutant_dishes,
            batch_name='mutant',
            snapshot_cells=second_snapshot_cells,
            mix_ratio=args.mix_ratio,
            expected_population=expected_population,
            expected_final_cells=expected_final_cells,
            base_seed=args.seed + 100
        )
    
        # Mix control1 individuals in memory
        process_batch_mixing(
            dishes=control1_dishes,
            batch_name='control1',
            snapshot_cells=second_snapshot_cells,
            mix_ratio=args.mix_ratio,
            expected_population=expected_population,
            expected_final_cells=expected_final_cells,
            base_seed=args.seed + 200
        )
    
    # Store uniform pool data for Control2 creation (if using uniform mixing)
    uniform_pool_for_control2 = None
    uniform_indices_for_control2 = None
    if args.uniform_mixing:
        uniform_pool_for_control2 = uniform_pool
        uniform_indices_for_control2 = uniform_indices
    
    # ========================================================================
    # VALIDATION: Mixed populations
    # ========================================================================
    try:
        validator.validate_mixed_populations(
            mutant_dishes=mutant_dishes,
            control1_dishes=control1_dishes,
            mix_ratio=args.mix_ratio,
            uniform_mixing=args.uniform_mixing,
            second_snapshot_year=args.second_snapshot
        )
    except ValidationError as e:
        print(f"\n❌ Mixing validation failed: {e}")
        sys.exit(1)
    
    # ========================================================================
    # VALIDATION: Pre-save check
    # ========================================================================
    try:
        validator.validate_before_save(
            mutant_dishes=mutant_dishes,
            control1_dishes=control1_dishes
        )
    except ValidationError as e:
        print(f"\n❌ Pre-save validation failed: {e}")
        sys.exit(1)
    
    # ========================================================================
    # SINGLE SAVE AT END OF STAGE 6 - Save all mutant and control1 individuals
    # ========================================================================
    print("\n  === SAVING ALL INDIVIDUALS TO DISK ===")
    
    # Import the new save function
    from core.pipeline_utils import save_all_individuals
    
    # Save mutant individuals
    save_all_individuals(
        dishes=mutant_dishes,
        output_dir=mutant_dir,
        batch_name='mutant',
        compress=args.use_compression,
        include_cell_history=True,
        include_gene_metrics=True
    )
    
    # Save control1 individuals
    save_all_individuals(
        dishes=control1_dishes,
        output_dir=control1_dir,
        batch_name='control1',
        compress=args.use_compression,
        include_cell_history=True,
        include_gene_metrics=True
    )
    
    # ========================================================================
    # STAGE 7: Create Control2 Individuals (Pure Second Snapshot)
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 7: Create Control2 Individuals")
    print(f"{'='*60}")
    
    # Use in-memory dish counts (already have current sizes after mixing)
    # Determine how many control2 individuals to create
    # Should match the average of mutant and control1 after normalization
    if args.normalize_size:
        # After normalization, we may have fewer individuals
        num_control2 = (len(mutant_dishes) + len(control1_dishes)) // 2
        print(f"  Adjusted control2 count after normalization: {num_control2}")
        print(f"    (Based on {len(mutant_dishes)} mutant + {len(control1_dishes)} control1)")
    else:
        # Use original expected count
        num_control2 = expected_individuals
    
    # Calculate expected final cells (needed for both mixing modes)
    if args.uniform_mixing:
        # With new normalization: all individuals have same final size
        # Calculate based on normalized size + uniform pool
        if 'normalized_size' in locals():
            # Use our normalized size + pool calculation
            expected_final_cells = int(normalized_size / (1 - args.mix_ratio / 100))
        else:
            # Fallback: Use median of actual mixed sizes
            all_sizes = [len(dish.cells) for dish in mutant_dishes] + [len(dish.cells) for dish in control1_dishes]
            expected_final_cells = int(np.median(all_sizes))
        print(f"  Using uniform mixing final size: {expected_final_cells} cells")
    else:
        # Use the standard calculation
        grown_cells = expected_population
        grown_fraction = (100 - args.mix_ratio) / 100
        expected_final_cells = int(grown_cells / grown_fraction)
    
    control2_state = check_petri_files_state(control2_dir, expected_final_cells)
    
    control2_dishes = []
    if control2_state['total_files'] >= num_control2 and not args.force_recreate:
        print(f"  ⏭ Control2 individuals exist ({control2_state['total_files']} files)")
        control2_dishes = load_all_petri_dishes(control2_dir)
    else:
        print(f"  ✓ Creating {num_control2} control2 individuals (pure year {args.second_snapshot})...")
        
        # Clean up any existing control2 files first
        import glob
        for old_file in glob.glob(os.path.join(control2_dir, f"individual_*{individual_ext}")):
            os.remove(old_file)
        
        # Adjust target size if snapshot has fewer cells
        actual_control2_size = min(expected_final_cells, len(second_snapshot_cells))
        if actual_control2_size < expected_final_cells:
            print(f"    Warning: Snapshot has only {len(second_snapshot_cells)} cells, adjusting control2 size")
        print(f"    Each with {actual_control2_size} pure year {args.second_snapshot} cells")
        
        for i in range(num_control2):
            print(f"    Creating individual {i+1}/{num_control2}")
            
            # IMPORTANT: individual_id must match the filename index to ensure consistency
            # across multiple runs (even if files are partially created/deleted)
            file_index = i + 1  # Start IDs at 1 instead of 0
            
            # Create PetriDish based on mixing mode
            if args.uniform_mixing and uniform_pool_for_control2 is not None:
                # Use uniform base + additional sampling
                print(f"      Using uniform base + additional cells")
                petri = create_control2_with_uniform_base(
                    second_snapshot_cells,
                    uniform_pool_for_control2,
                    uniform_indices_for_control2,
                    actual_control2_size,
                    seed=args.seed + 300 + i
                )
                
                # Add metadata indicating uniform base was used
                if not hasattr(petri, 'metadata'):
                    petri.metadata = {}
                petri.metadata.update({
                    'individual_id': file_index,  # Must match filename for consistency
                    'individual_type': 'control2',
                    'source': f'uniform_base_plus_snapshot_year{args.second_snapshot}',
                    'uniform_base': True
                    # Note: year is now correctly set from cell ages (50)
                })
            else:
                # Original random sampling
                petri = create_pure_snapshot_petri(second_snapshot_cells, n_cells=actual_control2_size,
                                                  seed=args.seed + 300 + i)
                
                # Add metadata
                if not hasattr(petri, 'metadata'):
                    petri.metadata = {}
                petri.metadata.update({
                    'individual_id': file_index,  # Must match filename for consistency
                    'individual_type': 'control2',
                    'source': f'pure_year{args.second_snapshot}'
                    # Note: year is now correctly set from cell ages (50)
                })
            
            control2_dishes.append(petri)
            
            # Save without history - control2 doesn't need it (pure snapshot, no growth)
            filepath = os.path.join(control2_dir, f"individual_{file_index:02d}{individual_ext}")
            save_petri_dish(petri, filepath, include_cell_history=False,  # No history needed!
                          include_gene_metrics=True, compress=args.use_compression)
        
        print(f"  Created and saved {len(control2_dishes)} control2 individuals")
    
    # Validate control2 has correct gene_rate_groups
    validate_all_gene_rate_groups(
        expected=gene_rate_groups,
        petri_dishes=control2_dishes,
        label="Control2 individuals"
    )
    
    # ========================================================================
    # VALIDATION: Control2 individuals
    # ========================================================================
    try:
        validator.validate_control2(
            control2_dishes=control2_dishes,
            second_snapshot_year=args.second_snapshot,
            expected_count=num_control2
        )
    except ValidationError as e:
        print(f"\n❌ Control2 validation failed: {e}")
        sys.exit(1)
    
    # ========================================================================
    # Generate Individual Growth Trajectories
    # ========================================================================
    print(f"\n{'='*60}")
    print("Generating Individual Growth Trajectories")
    print(f"{'='*60}")
    
    # Create output directory for individual trajectories
    individual_trajectories_dir = os.path.join(results_dir, 'individual_trajectories')
    os.makedirs(individual_trajectories_dir, exist_ok=True)
    
    # Plot mutant individuals
    print("\n  Creating plots for mutant individuals...")
    mutant_plot_count = 0
    for i, filename in enumerate(sorted(os.listdir(mutant_dir)), 1):
        if filename.endswith(('.json', '.json.gz')):
            filepath = os.path.join(mutant_dir, filename)
            try:
                # Load with history
                petri = load_petri_dish(filepath, include_cell_history=True)
            
                # Check if history exists
                if hasattr(petri, 'cell_history') and petri.cell_history:
                    from cell import PetriDishPlotter
                    plotter = PetriDishPlotter(petri)
                
                    # Generate plots
                    jsd_path = os.path.join(individual_trajectories_dir, f"mutant_{i:02d}_jsd.png")
                    meth_path = os.path.join(individual_trajectories_dir, f"mutant_{i:02d}_methylation.png")
                
                    plotter.plot_jsd(f"Mutant Individual {i:02d}", jsd_path)
                    plotter.plot_methylation(f"Mutant Individual {i:02d}", meth_path)
                
                    mutant_plot_count += 1
                    print(f"    ✓ Plotted mutant_{i:02d} trajectories")
            except Exception as e:
                print(f"    ⚠️  Could not plot mutant_{i:02d}: {str(e)}")
    
    # Plot control1 individuals
    print("\n  Creating plots for control1 individuals...")
    control1_plot_count = 0
    for i, filename in enumerate(sorted(os.listdir(control1_dir)), 1):
        if filename.endswith(('.json', '.json.gz')):
            filepath = os.path.join(control1_dir, filename)
            try:
                # Load with history
                petri = load_petri_dish(filepath, include_cell_history=True)
            
                # Check if history exists
                if hasattr(petri, 'cell_history') and petri.cell_history:
                    from cell import PetriDishPlotter
                    plotter = PetriDishPlotter(petri)
                
                    # Generate plots
                    jsd_path = os.path.join(individual_trajectories_dir, f"control1_{i:02d}_jsd.png")
                    meth_path = os.path.join(individual_trajectories_dir, f"control1_{i:02d}_methylation.png")
                
                    plotter.plot_jsd(f"Control1 Individual {i:02d}", jsd_path)
                    plotter.plot_methylation(f"Control1 Individual {i:02d}", meth_path)
                
                    control1_plot_count += 1
                    print(f"    ✓ Plotted control1_{i:02d} trajectories")
            except Exception as e:
                print(f"    ⚠️  Could not plot control1_{i:02d}: {str(e)}")
    
    print(f"\n  Skipping control2 (pure snapshots, no growth trajectories)")
    print(f"  ✓ Generated {(mutant_plot_count + control1_plot_count) * 2} trajectory plots in {individual_trajectories_dir}")
    
    # ========================================================================
    # STAGE 8: Analysis
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 8: Analysis and Visualization")
    print(f"{'='*60}")
    
    # Use in-memory dishes for mutant and control1, load control2 from disk
    print("\n  Using final populations for analysis...")
    # mutant_dishes and control1_dishes are already in memory from Stage 6
    # Only need to load control2 from disk
    control2_dishes = load_all_petri_dishes(control2_dir)
    
    print(f"    Loaded {len(mutant_dishes)} mutant individuals")
    print(f"    Loaded {len(control1_dishes)} control1 individuals")
    print(f"    Loaded {len(control2_dishes)} control2 individuals")
    
    # Use the enhanced analysis function for cell-level JSD
    analysis_results = analyze_populations_from_dishes(
        mutant_dishes, control1_dishes, control2_dishes, results_dir
    )
    
    # Generate gene-level JSD analysis
    from core.pipeline_analysis import generate_gene_jsd_analysis
    gene_analysis_results = generate_gene_jsd_analysis(
        mutant_dishes, control1_dishes, control2_dishes, results_dir
    )
    
    # Generate gene-level JSD distribution plots
    if 'gene_jsd_analysis' in gene_analysis_results:
        plot_gene_jsd_distributions(gene_analysis_results['gene_jsd_analysis'], results_dir, 
                                   max_genes=args.max_gene_plots)
    
        # Generate individual-averaged gene JSD comparison plot
        gene_jsd_path = os.path.join(results_dir, 'gene_jsd_analysis.json')
        if os.path.exists(gene_jsd_path):
            plot_gene_jsd_individual_comparison(
                gene_jsd_path=gene_jsd_path,
                output_dir=results_dir,
                verbose=True  # Always verbose for now
            )
    
    # Generate additional gene JSD plots if data is available
    print(f"\nGenerating additional gene JSD plots...")
    
    # Check if we have snapshots for gene JSD distribution comparison
    try:
        # Load snapshots for gene JSD comparison
        snapshot1_path = os.path.join(snapshots_dir, f"year{args.first_snapshot}_snapshot.json.gz")
        snapshot2_path = os.path.join(snapshots_dir, f"year{args.second_snapshot}_snapshot.json.gz")
        snapshot1_cells = load_snapshot_cells(snapshot1_path)
        snapshot2_cells = load_snapshot_cells(snapshot2_path)
    
        # Gene JSD distribution comparison (year vs year)
        gene_dist_path = os.path.join(results_dir, "gene_jsd_distribution.png")
        plot_gene_jsd_distribution_comparison(
            snapshot1_cells, snapshot2_cells, 
            args.first_snapshot, args.second_snapshot, gene_dist_path
        )
    
    except Exception as e:
        print(f"  Skipping gene JSD distribution comparison: {e}")
    
    # Gene vs Cell JSD scatter plot - removed per user request
    # (This plot was not providing useful insights)
    
    # Top variable genes plot
    try:
        # Combine all dishes for gene variability analysis
        all_dishes = mutant_dishes + control1_dishes + control2_dishes
        top_genes_path = os.path.join(results_dir, "top_variable_genes.png")
        plot_top_variable_genes(all_dishes, n_top=20, output_path=top_genes_path)
    except Exception as e:
        print(f"  Skipping top variable genes plot: {e}")
    
    # Generate phase1-style gene JSD plots using original simulation data
    try:
        print(f"\nGenerating phase1-style gene JSD plots from simulation data...")
    
        # Load the original simulation data
        print(f"  Loading original simulation: {args.simulation}")
        with gzip.open(args.simulation, 'rt') as f:
            sim_data = json.load(f)
    
        # Check if gene_jsd_history exists in simulation
        if 'gene_jsd_history' in sim_data:
            # Create a temporary PetriDish with gene_jsd_history for plotting
            temp_petri = PetriDish()
            temp_petri.gene_jsd_history = {int(year): values for year, values in sim_data['gene_jsd_history'].items()}
            
            # Extract mean and median gene JSD histories from the history structure if available
            temp_petri.mean_gene_jsd_history = {}
            temp_petri.median_gene_jsd_history = {}
            if 'history' in sim_data:
                for year_str, year_data in sim_data['history'].items():
                    if 'mean_gene_jsd' in year_data:
                        temp_petri.mean_gene_jsd_history[int(year_str)] = year_data['mean_gene_jsd']
                    if 'median_gene_jsd' in year_data:
                        temp_petri.median_gene_jsd_history[int(year_str)] = year_data['median_gene_jsd']
        
            # Add some cells for gene rate group detection (if available)
            if sim_data.get('history') and len(sim_data['history']) > 0:
                # Get cells from the last year for gene rate group info
                last_year_data = list(sim_data['history'].values())[-1]
                if 'cells' in last_year_data and len(last_year_data['cells']) > 0:
                    # Convert first cell to get gene rate group info
                    from core.pipeline_utils import dict_to_cell
                    sample_cell = dict_to_cell(last_year_data['cells'][0])
                    temp_petri.cells = [sample_cell]  # Just need one for gene rate group detection
        
            # Create plotter and generate phase1-style plots
            from cell import PetriDishPlotter
            plotter = PetriDishPlotter(temp_petri)
        
            # Gene JSD Heatmap
            heatmap_path = os.path.join(results_dir, "simulation_gene_jsd_heatmap.png")
            plotter.plot_gene_jsd_heatmap(
                title="Gene JSD Evolution Heatmap (from Phase1 simulation)",
                output_path=heatmap_path
            )
        
            # Gene Rate Group Comparison (only if gene rate groups exist)
            if (temp_petri.cells and hasattr(temp_petri.cells[0], 'gene_rate_groups') 
                and temp_petri.cells[0].gene_rate_groups):
                rate_comparison_path = os.path.join(results_dir, "simulation_gene_rate_comparison.png")
                plotter.plot_gene_jsd_by_rate_group(
                    title="Gene JSD by Rate Group (from Phase1 simulation)",
                    output_path=rate_comparison_path
                )
            else:
                print(f"    No gene rate groups found - skipping rate comparison plot")
            
        else:
            print(f"    No gene_jsd_history found in simulation")
            print(f"    This simulation was likely run with --no-gene-jsd or --no-jsds flags")
            print(f"    Gene JSD tracking is now enabled by default in new simulations")
        
    except Exception as e:
        print(f"  Skipping phase1-style gene JSD plots: {e}")
    
    # ========================================================================
    # Generate Original Simulation Timeline Plots
    # ========================================================================
    print(f"\n{'='*60}")
    print("Generating Original Simulation Timeline Plots")
    print(f"{'='*60}")
    
    try:
        print(f"  Loading simulation with full history...")
    
        # Determine if simulation file is compressed
        if args.simulation.endswith('.gz'):
            with gzip.open(args.simulation, 'rt') as f:
                sim_data = json.load(f)
        else:
            with open(args.simulation, 'r') as f:
                sim_data = json.load(f)
    
        # Convert history to PetriDish format
        if 'history' in sim_data and sim_data['history']:
            # Create a PetriDish with history
            from core.pipeline_utils import dict_to_cell
        
            # Get parameters from simulation
            params = sim_data.get('parameters', {})
            rate = params.get('rate')
            gene_rate_groups = params.get('gene_rate_groups')
            n_sites = params.get('n', 1000)
            gene_size = params.get('gene_size', 5)
        
            # Create PetriDish with appropriate configuration
            if gene_rate_groups:
                # Convert gene_rate_groups format
                rate_groups = [(g[0], g[1]) for g in gene_rate_groups]
                original_petri = PetriDish(gene_rate_groups=rate_groups, n=n_sites, gene_size=gene_size)
            else:
                original_petri = PetriDish(rate=rate, n=n_sites, gene_size=gene_size)
        
            # Build cell history from simulation data
            original_petri.cell_history = {}
            original_petri.years_simulated = len(sim_data['history']) - 1
        
            for year_str, year_data in sim_data['history'].items():
                if 'cells' in year_data:
                    # Keep year as string for PetriDishPlotter compatibility
                    # Store cells as dictionaries for plotter compatibility
                    cells = []
                    for cell_dict in year_data['cells']:
                        # Handle new lean format - convert to dict format expected by plotter
                        if 'methylated' in cell_dict:
                            cell_data = {
                                'methylated': cell_dict['methylated'],
                                'cell_JSD': cell_dict.get('cell_JSD', 0.0),
                                'cell_jsd': cell_dict.get('cell_JSD', 0.0),  # Include both naming conventions
                                'methylation_proportion': sum(cell_dict['methylated']) / len(cell_dict['methylated']),
                                'age': cell_dict.get('age', int(year_str))
                            }
                            cells.append(cell_data)
                        else:
                            # Old format - pass through
                            cells.append(cell_dict)
                    original_petri.cell_history[year_str] = cells
        
            # Also copy gene_jsd_history if available
            if 'gene_jsd_history' in sim_data:
                original_petri.gene_jsd_history = {
                    int(year): values for year, values in sim_data['gene_jsd_history'].items()
                }
        
            # Create plotter and generate timeline plots
            plotter = PetriDishPlotter(original_petri)
        
            # JSD timeline
            jsd_timeline_path = os.path.join(results_dir, "original_simulation_jsd_timeline.png")
            plotter.plot_jsd("Original Simulation JSD Timeline", jsd_timeline_path)
            print(f"  ✓ Generated original_simulation_jsd_timeline.png")
        
            # Methylation timeline
            meth_timeline_path = os.path.join(results_dir, "original_simulation_methylation_timeline.png")
            plotter.plot_methylation("Original Simulation Methylation Timeline", meth_timeline_path)
            print(f"  ✓ Generated original_simulation_methylation_timeline.png")
        
            # Combined plot if desired
            combined_timeline_path = os.path.join(results_dir, "original_simulation_combined_timeline.png")
            plotter.plot_combined("Original Simulation Combined Timeline", combined_timeline_path)
            print(f"  ✓ Generated original_simulation_combined_timeline.png")
        
        else:
            print(f"  ⚠️  No history found in simulation file")
            print(f"      The simulation may have been run without --track-cell-history")
        
    except Exception as e:
        import traceback
        print(f"  ⚠️  Could not generate original simulation timeline plots: {str(e)}")
        if hasattr(args, 'verbose') and args.verbose:
            print(f"      Full error: {traceback.format_exc()}")
    
    # Get statistics for summary from new consolidated format
    cell_jsd_data = analysis_results['cell_jsd_analysis']
    summary_stats = cell_jsd_data['summary_statistics']
    statistical_tests = cell_jsd_data['statistical_tests']
    
    # ========================================================================
    # Summary
    # ========================================================================
    elapsed_time = time.time() - start_time
    
    # Print validation summary
    validation_summary = validator.get_summary()
    if validation_summary['warnings_count'] > 0:
        print(f"\n{'='*80}")
        print("VALIDATION SUMMARY")
        print(f"{'='*80}")
        print(f"Total warnings: {validation_summary['warnings_count']}")
        if validation_summary['warnings_count'] <= 10:
            for warning in validation_summary['warnings']:
                print(f"  ⚠️  {warning}")
        else:
            print(f"  (Showing first 10 of {validation_summary['warnings_count']} warnings)")
            for warning in validation_summary['warnings'][:10]:
                print(f"  ⚠️  {warning}")
    
    print(f"\n{'='*80}")
    print("PIPELINE COMPLETE")
    print(f"{'='*80}")
    print(f"Total time: {elapsed_time/60:.1f} minutes")
    print(f"Output directory: {base_dir}")
    
    print("\nKey results:")
    print(f"  Mutant mean JSD: {summary_stats['mutant']['mean']:.6f} ± {summary_stats['mutant']['std']:.6f}")
    print(f"  Control1 mean JSD: {summary_stats['control1']['mean']:.6f} ± {summary_stats['control1']['std']:.6f}")
    print(f"  Control2 mean JSD: {summary_stats['control2']['mean']:.6f} ± {summary_stats['control2']['std']:.6f}")
    
    print("\nStatistical tests:")
    for comparison, values in statistical_tests.items():
        print(f"  {comparison}: p={values['p_value']:.6f}")
    
    # Save pipeline metadata (parameters only, no results)
    metadata = {
        "pipeline_version": "phase2",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_time_minutes": elapsed_time / 60,
        "parameters": {
            "rate": args.rate,
                "first_snapshot": args.first_snapshot,
                "second_snapshot": args.second_snapshot,
                "individual_growth_phase": args.individual_growth_phase,
                "timeline_duration": timeline_duration,
                "homeostasis_years": homeostasis_years,
                "expected_population": expected_population,
                "n_quantiles": args.n_quantiles,
                "cells_per_quantile": args.cells_per_quantile,
                "total_individuals": expected_individuals,
                "mix_ratio": args.mix_ratio,
                "uniform_mixing": args.uniform_mixing,
                "normalize_size": args.normalize_size,
                "normalization_threshold": normalization_threshold,
                "seed": args.seed,
                "bins": args.bins,
                "plot_individuals": args.plot_individuals
            }
        }
    
    metadata_path = os.path.join(results_dir, "pipeline_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\nPipeline metadata saved to: {metadata_path}")
    
    # Generate individual plots if requested
    if args.plot_individuals:
        print(f"\n{'='*60}")
        print("Generating Individual Growth Trajectory Plots")
        print(f"{'='*60}")
        
        from visualization.plot_individuals import plot_all_individuals
        plot_all_individuals(base_dir, plot_combined=True)
    
    return summary_stats


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 Pipeline using PetriDish and Cell classes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Config file argument
    parser.add_argument("--config", type=str,
                       help="Path to YAML configuration file")
    
    # Required arguments (can be set via config)
    # Rate configuration - mutually exclusive group
    rate_group = parser.add_mutually_exclusive_group(required=False)
    rate_group.add_argument("--rate", type=float,
                           help="Uniform methylation rate (must match simulation)")
    rate_group.add_argument("--gene-rate-groups", type=str,
                           help="Gene-specific methylation rates (format: 'n1:rate1,n2:rate2,...')")
    parser.add_argument("--gene-size", type=int, default=5,
                       help="Sites per gene (only used with --gene-rate-groups, default: 5)")
    
    parser.add_argument("--simulation", type=str,
                       help="Path to phase1 simulation file")
    
    # Quantile sampling parameters
    parser.add_argument("--n-quantiles", type=int, default=10,
                       help="Number of quantiles for sampling (e.g., 4 for quartiles, 10 for deciles)")
    parser.add_argument("--cells-per-quantile", type=int, default=3,
                       help="Number of cells to sample per quantile")
    
    # Snapshot and growth parameters
    parser.add_argument("--first-snapshot", type=int, default=50,
                       help="Year to extract initial cells from simulation")
    parser.add_argument("--second-snapshot", type=int, default=60,
                       help="Year to extract mixing cells from simulation")
    parser.add_argument("--individual-growth-phase", type=int, default=7,
                       help="Years of exponential growth before homeostasis (7=128 cells, 8=256 cells)")
    parser.add_argument("--mix-ratio", type=int, default=80,
                       help="Percentage of second snapshot cells in mix (0-100)")
    
    # Visualization parameters
    parser.add_argument("--bins", type=int, default=200,
                       help="Number of bins for JSD histograms")
    parser.add_argument("--plot-individuals", action='store_true',
                       help="Generate growth trajectory plots for each individual")
    parser.add_argument("--max-gene-plots", type=int, default=None,
                       help="Maximum number of gene JSD plots to generate (default: all)")
    parser.add_argument("--no-compress", action='store_true',
                       help="Save output as uncompressed JSON instead of .json.gz (larger files but easier to inspect)")
    
    # Mixing parameters
    parser.add_argument("--uniform-mixing", action='store_true',
                       help="Use same snapshot cells for all individuals (default: independent random sampling)")
    parser.add_argument("--normalize-size", action='store_true',
                       help="Normalize all individuals to same size before mixing (uses median - 0.5σ threshold)")
    
    # Other parameters
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducibility")
    parser.add_argument("--output-dir", type=str, default="data",
                       help="Output directory")
    
    # Force options
    parser.add_argument("--force-reload", action='store_true',
                       help="Force reload of snapshots even if cached")
    parser.add_argument("--force-recreate", action='store_true',
                       help="Force recreation of individuals even if they exist")
    
    args = parser.parse_args()
    
    # Load config file if provided
    config = load_config(args.config)
    
    # Merge config with CLI arguments
    args = merge_config_and_args(config, args)
    
    # Validate final configuration
    validate_pipeline_config(args)
    
    # Validate rate configuration (rate or gene_rate_groups must be specified)
    # Actual parsing and validation happens in run_pipeline()
    
    # Handle glob patterns in simulation file path
    import glob
    if '*' in args.simulation:
        matching_files = glob.glob(args.simulation)
        if not matching_files:
            print(f"Error: No files matching pattern: {args.simulation}")
            sys.exit(1)
        # Sort to prefer .json.gz over .json if both exist
        matching_files.sort(key=lambda x: (0 if x.endswith('.json.gz') else 1))
        args.simulation = matching_files[0]
        if len(matching_files) > 1:
            print(f"Multiple files found, using: {args.simulation}")
    
    # Check simulation file exists and has valid extension
    if not os.path.exists(args.simulation):
        print(f"Error: Simulation file not found: {args.simulation}")
        sys.exit(1)
    
    if not (args.simulation.endswith('.json') or args.simulation.endswith('.json.gz')):
        print(f"Error: Simulation file must be .json or .json.gz, got: {args.simulation}")
        sys.exit(1)
    
    # Run pipeline (no rate_config needed anymore)
    run_pipeline(args)


if __name__ == "__main__":
    main()