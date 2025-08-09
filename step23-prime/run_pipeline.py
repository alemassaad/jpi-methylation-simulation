#!/usr/bin/env python3
"""
Step23-Prime Pipeline: Refactored to use PetriDish and Cell classes directly.

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
from typing import List, Dict, Optional

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'step1-prime'))
from cell import PetriDish, Cell
from pipeline_utils import (
    load_snapshot_as_cells, save_snapshot_cells, load_snapshot_cells,
    save_petri_dish, load_petri_dish, load_all_petri_dishes,
    sample_by_quantiles, sample_uniform,
    mix_petri_with_snapshot, create_pure_snapshot_petri,
    grow_petri_for_years, get_petri_statistics, check_petri_files_state,
    get_jsd_array
)
from pipeline_analysis import (
    plot_jsd_distribution_from_cells,
    analyze_populations_from_dishes,
    plot_cell_level_distributions
)
from path_utils import parse_step1_simulation_path, generate_step23_output_dir


def run_pipeline(args):
    """
    Main pipeline orchestrator using PetriDish objects.
    """
    # Calculate snapshot years (always dynamic now)
    first_snapshot_year = args.snapshot_year
    second_snapshot_year = args.snapshot_year + args.growth_years
    
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
    
    # Parse simulation parameters from input file
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
    plots_dir = os.path.join(base_dir, "plots")
    results_dir = os.path.join(base_dir, "results")
    
    # Create all directories
    for dir_path in [snapshots_dir, mutant_dir, control1_dir, control2_dir, plots_dir, results_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    print("=" * 80)
    print("STEP23-PRIME PIPELINE")
    print("=" * 80)
    print(f"Rate: {args.rate}")
    print(f"Simulation: {args.simulation}")
    print(f"Output: {base_dir}")
    print(f"Quantiles: {args.n_quantiles} ({'quartiles' if args.n_quantiles == 4 else 'deciles' if args.n_quantiles == 10 else f'{args.n_quantiles}-tiles'})")
    print(f"Cells per quantile: {args.cells_per_quantile}")
    print(f"Total individuals: {expected_individuals} per group")
    print(f"Growth years: {args.growth_years} (target: {2**args.growth_years} cells)")
    print(f"Snapshot years: {first_snapshot_year} → {second_snapshot_year} (growth: {args.growth_years} years)")
    print(f"Mix ratio: {args.mix_ratio}% year {second_snapshot_year}, {100-args.mix_ratio}% grown")
    print(f"Seed: {args.seed}")
    print("=" * 80)
    
    # ========================================================================
    # STAGE 1: Extract First Snapshot
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 1: Extract Year {first_snapshot_year} Snapshot")
    print(f"{'='*60}")
    
    first_snapshot_path = os.path.join(snapshots_dir, f"year{first_snapshot_year}_snapshot.json.gz")
    
    if os.path.exists(first_snapshot_path) and not args.force_reload:
        print(f"  ⏭ Loading existing snapshot from {first_snapshot_path}...")
        first_snapshot_cells = load_snapshot_cells(first_snapshot_path)
        print(f"  Loaded {len(first_snapshot_cells)} Cell objects")
    else:
        print(f"  ✓ Extracting year {first_snapshot_year} from simulation...")
        first_snapshot_cells = load_snapshot_as_cells(args.simulation, first_snapshot_year)
        save_snapshot_cells(first_snapshot_cells, first_snapshot_path)
        print(f"  Cached {len(first_snapshot_cells)} cells for future use")
    
    # ========================================================================
    # STAGE 2: Plot JSD Distribution
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 2: Plot JSD Distribution")
    print(f"{'='*60}")
    
    plot_path = os.path.join(plots_dir, f"year{first_snapshot_year}_jsd_distribution_{args.bins}bins.png")
    plot_jsd_distribution_from_cells(first_snapshot_cells, args.bins, plot_path, rate=args.rate)
    
    # ========================================================================
    # STAGE 3: Create Initial Individuals (as PetriDish objects)
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 3: Create Initial Individuals")
    print(f"{'='*60}")
    
    # Check existing mutant state
    mutant_state = check_petri_files_state(mutant_dir, expected_cells=1)
    
    mutant_dishes = []
    if mutant_state['total_files'] >= expected_individuals and not args.force_recreate:
        print(f"\n  ⏭ Mutant individuals exist ({mutant_state['total_files']} files), loading...")
        mutant_dishes = load_all_petri_dishes(mutant_dir)
        print(f"  Loaded {len(mutant_dishes)} mutant PetriDish objects")
    else:
        print(f"\n  ✓ Creating {expected_individuals} mutant individuals...")
        print(f"  Sampling by quantiles ({args.n_quantiles} quantiles, {args.cells_per_quantile} each)...")
        
        sampled = sample_by_quantiles(first_snapshot_cells, 
                                     n_quantiles=args.n_quantiles,
                                     cells_per_quantile=args.cells_per_quantile,
                                     seed=args.seed)
        
        for i, (cell, quantile) in enumerate(sampled):
            # Create PetriDish with single cell
            petri = PetriDish(rate=args.rate, n=len(cell.cpg_sites), 
                            gene_size=cell.gene_size, seed=None)
            petri.cells = [cell]  # Start with 1 cell
            petri.year = 50
            
            # Add metadata
            if not hasattr(petri, 'metadata'):
                petri.metadata = {}
            petri.metadata.update({
                'individual_id': i,
                'individual_type': 'mutant',
                'source_quantile': quantile,
                'initial_year': 50,
                'n_quantiles': args.n_quantiles
            })
            
            mutant_dishes.append(petri)
            
            # Save immediately
            filepath = os.path.join(mutant_dir, f"individual_{i:02d}.json.gz")
            save_petri_dish(petri, filepath)
        
        print(f"  Created and saved {len(mutant_dishes)} mutant individuals")
    
    # Check existing control1 state
    control1_state = check_petri_files_state(control1_dir, expected_cells=1)
    
    control1_dishes = []
    if control1_state['total_files'] >= expected_individuals and not args.force_recreate:
        print(f"\n  ⏭ Control1 individuals exist ({control1_state['total_files']} files), loading...")
        control1_dishes = load_all_petri_dishes(control1_dir)
        print(f"  Loaded {len(control1_dishes)} control1 PetriDish objects")
    else:
        print(f"\n  ✓ Creating {expected_individuals} control1 individuals...")
        print(f"  Sampling uniformly...")
        
        sampled_cells = sample_uniform(first_snapshot_cells, n_samples=expected_individuals, 
                                      seed=args.seed + 1000)
        
        for i, cell in enumerate(sampled_cells):
            # Create PetriDish with single cell
            petri = PetriDish(rate=args.rate, n=len(cell.cpg_sites),
                            gene_size=cell.gene_size, seed=None)
            petri.cells = [cell]
            petri.year = 50
            
            # Add metadata
            if not hasattr(petri, 'metadata'):
                petri.metadata = {}
            petri.metadata.update({
                'individual_id': i,
                'individual_type': 'control1',
                'source': 'uniform',
                'initial_year': 50
            })
            
            control1_dishes.append(petri)
            
            # Save immediately
            filepath = os.path.join(control1_dir, f"individual_{i:02d}.json.gz")
            save_petri_dish(petri, filepath)
        
        print(f"  Created and saved {len(control1_dishes)} control1 individuals")
    
    # ========================================================================
    # STAGE 4: Grow Individuals (10 years using PetriDish methods)
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 4: Grow Individuals ({args.growth_years} years)")
    print(f"{'='*60}")
    
    expected_cells_after_growth = 2 ** args.growth_years
    
    # Check mutant growth state
    mutant_state = check_petri_files_state(mutant_dir, expected_cells_after_growth)
    
    if mutant_state['all_expected'] or mutant_state['all_above']:
        print(f"\n  ⏭ Mutant individuals already grown/mixed:")
        print(f"    At target ({expected_cells_after_growth} cells): {mutant_state['expected_count']}")
        print(f"    Above target (mixed): {mutant_state['above_count']}")
    else:
        print(f"\n  ✓ Growing mutant individuals to {expected_cells_after_growth} cells...")
        
        # Reload to get current state
        mutant_dishes = load_all_petri_dishes(mutant_dir)
        
        for i, petri in enumerate(mutant_dishes):
            current_cells = len(petri.cells)
            
            if current_cells < expected_cells_after_growth:
                # Calculate remaining growth needed
                current_doublings = math.log2(current_cells) if current_cells > 0 else 0
                target_doublings = args.growth_years
                remaining_years = int(target_doublings - current_doublings)
                
                if remaining_years > 0:
                    print(f"    Individual {i:02d}: {current_cells} → {expected_cells_after_growth} cells")
                    grow_petri_for_years(petri, remaining_years, verbose=True)
                    
                    # Save updated state
                    filepath = os.path.join(mutant_dir, f"individual_{i:02d}.json.gz")
                    save_petri_dish(petri, filepath)
            elif current_cells == expected_cells_after_growth:
                print(f"    Individual {i:02d}: Already at {current_cells} cells")
            else:
                print(f"    Individual {i:02d}: Already mixed ({current_cells} cells)")
    
    # Check control1 growth state
    control1_state = check_petri_files_state(control1_dir, expected_cells_after_growth)
    
    if control1_state['all_expected'] or control1_state['all_above']:
        print(f"\n  ⏭ Control1 individuals already grown/mixed:")
        print(f"    At target ({expected_cells_after_growth} cells): {control1_state['expected_count']}")
        print(f"    Above target (mixed): {control1_state['above_count']}")
    else:
        print(f"\n  ✓ Growing control1 individuals to {expected_cells_after_growth} cells...")
        
        # Reload to get current state
        control1_dishes = load_all_petri_dishes(control1_dir)
        
        for i, petri in enumerate(control1_dishes):
            current_cells = len(petri.cells)
            
            if current_cells < expected_cells_after_growth:
                # Calculate remaining growth needed
                current_doublings = math.log2(current_cells) if current_cells > 0 else 0
                target_doublings = args.growth_years
                remaining_years = int(target_doublings - current_doublings)
                
                if remaining_years > 0:
                    print(f"    Individual {i:02d}: {current_cells} → {expected_cells_after_growth} cells")
                    grow_petri_for_years(petri, remaining_years, verbose=True)
                    
                    # Save updated state
                    filepath = os.path.join(control1_dir, f"individual_{i:02d}.json.gz")
                    save_petri_dish(petri, filepath)
    
    # ========================================================================
    # STAGE 5: Extract Second Snapshot (Year 60 or dynamic)
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 5: Extract Year {second_snapshot_year} Snapshot")
    print(f"{'='*60}")
    
    second_snapshot_path = os.path.join(snapshots_dir, f"year{second_snapshot_year}_snapshot.json.gz")
    
    if os.path.exists(second_snapshot_path) and not args.force_reload:
        print(f"  ⏭ Loading existing snapshot from {second_snapshot_path}...")
        second_snapshot_cells = load_snapshot_cells(second_snapshot_path)
        print(f"  Loaded {len(second_snapshot_cells)} Cell objects")
    else:
        print(f"  ✓ Extracting year {second_snapshot_year} from simulation...")
        second_snapshot_cells = load_snapshot_as_cells(args.simulation, second_snapshot_year)
        save_snapshot_cells(second_snapshot_cells, second_snapshot_path)
        print(f"  Cached {len(second_snapshot_cells)} cells for future use")
    
    # ========================================================================
    # STAGE 6: Mix Populations
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 6: Mix Populations")
    print(f"{'='*60}")
    
    # Calculate expected final size after mixing
    grown_cells = expected_cells_after_growth
    expected_final_cells = int(grown_cells / ((100 - args.mix_ratio) / 100))
    
    print(f"  Target size after mixing: {expected_final_cells} cells")
    print(f"  Mix ratio: {args.mix_ratio}% from year {second_snapshot_year}")
    
    # Check mutant mix state
    mutant_state = check_petri_files_state(mutant_dir, expected_cells_after_growth)
    
    if mutant_state['all_above']:
        print(f"\n  ⏭ Mutant individuals already mixed")
    else:
        print(f"\n  ✓ Mixing mutant individuals with year {second_snapshot_year} cells...")
        
        # Reload current state
        mutant_dishes = load_all_petri_dishes(mutant_dir)
        
        for i, petri in enumerate(mutant_dishes):
            if len(petri.cells) == expected_cells_after_growth:
                print(f"    Individual {i:02d}: Mixing {len(petri.cells)} → {expected_final_cells} cells")
                
                total_cells = mix_petri_with_snapshot(petri, second_snapshot_cells,
                                                     mix_ratio=args.mix_ratio / 100,
                                                     seed=args.seed + 100 + i)
                
                # Update metadata
                if not hasattr(petri, 'metadata'):
                    petri.metadata = {}
                petri.metadata['mixed'] = True
                petri.metadata['mix_ratio'] = args.mix_ratio
                petri.metadata['final_cells'] = total_cells
                
                # Save updated state
                filepath = os.path.join(mutant_dir, f"individual_{i:02d}.json.gz")
                save_petri_dish(petri, filepath)
            elif len(petri.cells) > expected_cells_after_growth:
                print(f"    Individual {i:02d}: Already mixed ({len(petri.cells)} cells)")
    
    # Check control1 mix state
    control1_state = check_petri_files_state(control1_dir, expected_cells_after_growth)
    
    if control1_state['all_above']:
        print(f"\n  ⏭ Control1 individuals already mixed")
    else:
        print(f"\n  ✓ Mixing control1 individuals with year {second_snapshot_year} cells...")
        
        # Reload current state
        control1_dishes = load_all_petri_dishes(control1_dir)
        
        for i, petri in enumerate(control1_dishes):
            if len(petri.cells) == expected_cells_after_growth:
                print(f"    Individual {i:02d}: Mixing {len(petri.cells)} → {expected_final_cells} cells")
                
                total_cells = mix_petri_with_snapshot(petri, second_snapshot_cells,
                                                     mix_ratio=args.mix_ratio / 100,
                                                     seed=args.seed + 200 + i)
                
                # Update metadata
                if not hasattr(petri, 'metadata'):
                    petri.metadata = {}
                petri.metadata['mixed'] = True
                petri.metadata['mix_ratio'] = args.mix_ratio
                petri.metadata['final_cells'] = total_cells
                
                # Save updated state
                filepath = os.path.join(control1_dir, f"individual_{i:02d}.json.gz")
                save_petri_dish(petri, filepath)
    
    # ========================================================================
    # STAGE 7: Create Control2 Individuals (Pure Year 60)
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 7: Create Control2 Individuals")
    print(f"{'='*60}")
    
    control2_state = check_petri_files_state(control2_dir, expected_final_cells)
    
    control2_dishes = []
    if control2_state['total_files'] >= expected_individuals and not args.force_recreate:
        print(f"  ⏭ Control2 individuals exist ({control2_state['total_files']} files)")
        control2_dishes = load_all_petri_dishes(control2_dir)
    else:
        print(f"  ✓ Creating {expected_individuals} control2 individuals (pure year {second_snapshot_year})...")
        print(f"    Each with {expected_final_cells} pure year {second_snapshot_year} cells")
        
        for i in range(expected_individuals):
            print(f"    Creating individual {i+1}/{expected_individuals}")
            
            # Create PetriDish with pure year 60 cells
            petri = create_pure_snapshot_petri(second_snapshot_cells, n_cells=expected_final_cells,
                                              rate=args.rate, seed=args.seed + 300 + i)
            
            # Add metadata
            if not hasattr(petri, 'metadata'):
                petri.metadata = {}
            petri.metadata.update({
                'individual_id': i,
                'individual_type': 'control2',
                'source': 'pure_year60',
                'year': 60
            })
            
            control2_dishes.append(petri)
            
            # Save
            filepath = os.path.join(control2_dir, f"individual_{i:02d}.json.gz")
            save_petri_dish(petri, filepath)
        
        print(f"  Created and saved {len(control2_dishes)} control2 individuals")
    
    # ========================================================================
    # STAGE 8: Analysis
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 8: Analysis and Visualization")
    print(f"{'='*60}")
    
    # Reload all final states
    print("\n  Loading final populations for analysis...")
    mutant_dishes = load_all_petri_dishes(mutant_dir)
    control1_dishes = load_all_petri_dishes(control1_dir)
    control2_dishes = load_all_petri_dishes(control2_dir)
    
    print(f"    Loaded {len(mutant_dishes)} mutant individuals")
    print(f"    Loaded {len(control1_dishes)} control1 individuals")
    print(f"    Loaded {len(control2_dishes)} control2 individuals")
    
    # Use the enhanced analysis function
    analysis_results = analyze_populations_from_dishes(
        mutant_dishes, control1_dishes, control2_dishes, results_dir
    )
    
    # Create additional cell-level distribution plot
    cell_plot_path = os.path.join(plots_dir, "cell_level_jsd_distributions.png")
    plot_cell_level_distributions(mutant_dishes, control1_dishes, control2_dishes, cell_plot_path)
    
    # Get statistics for summary
    stats = analysis_results['statistics']
    
    # ========================================================================
    # Summary
    # ========================================================================
    elapsed_time = time.time() - start_time
    
    print(f"\n{'='*80}")
    print("PIPELINE COMPLETE")
    print(f"{'='*80}")
    print(f"Total time: {elapsed_time/60:.1f} minutes")
    print(f"Output directory: {base_dir}")
    
    print("\nKey results:")
    print(f"  Mutant mean JSD: {stats['mutant']['mean']:.6f} ± {stats['mutant']['std']:.6f}")
    print(f"  Control1 mean JSD: {stats['control1']['mean']:.6f} ± {stats['control1']['std']:.6f}")
    print(f"  Control2 mean JSD: {stats['control2']['mean']:.6f} ± {stats['control2']['std']:.6f}")
    
    print("\nStatistical tests:")
    for comparison, values in stats['comparisons'].items():
        print(f"  {comparison}: p={values['p_value']:.6f}")
    
    # Save pipeline metadata
    metadata = {
        "pipeline_version": "step23-prime",
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "elapsed_time_minutes": elapsed_time / 60,
        "parameters": {
            "rate": args.rate,
            "n_quantiles": args.n_quantiles,
            "cells_per_quantile": args.cells_per_quantile,
            "total_individuals": expected_individuals,
            "growth_years": args.growth_years,
            "mix_ratio": args.mix_ratio,
            "seed": args.seed,
            "bins": args.bins
        },
        "results_summary": stats
    }
    
    metadata_path = os.path.join(results_dir, "pipeline_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\nPipeline metadata saved to: {metadata_path}")
    
    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Step23-Prime Pipeline using PetriDish and Cell classes",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--rate", type=float, required=True,
                       help="Methylation rate (must match simulation)")
    parser.add_argument("--simulation", type=str, required=True,
                       help="Path to step1/step1-prime simulation file")
    
    # Quantile sampling parameters
    parser.add_argument("--n-quantiles", type=int, default=10,
                       help="Number of quantiles for sampling (e.g., 4 for quartiles, 10 for deciles)")
    parser.add_argument("--cells-per-quantile", type=int, default=3,
                       help="Number of cells to sample per quantile")
    
    # Snapshot and growth parameters
    parser.add_argument("--snapshot-year", type=int, default=50,
                       help="Year for first snapshot (default: 50)")
    parser.add_argument("--growth-years", type=int, default=10,
                       help="Number of years to grow individuals")
    parser.add_argument("--mix-ratio", type=int, default=80,
                       help="Percentage of second snapshot cells in mix (0-100)")
    
    # Visualization parameters
    parser.add_argument("--bins", type=int, default=200,
                       help="Number of bins for JSD histograms")
    
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
    
    # Validate arguments
    if not os.path.exists(args.simulation):
        print(f"Error: Simulation file not found: {args.simulation}")
        sys.exit(1)
    
    if not (0 <= args.mix_ratio <= 100):
        print(f"Error: mix-ratio must be between 0 and 100, got {args.mix_ratio}")
        sys.exit(1)
    
    if args.n_quantiles < 1:
        print(f"Error: n-quantiles must be at least 1, got {args.n_quantiles}")
        sys.exit(1)
    
    if args.cells_per_quantile < 1:
        print(f"Error: cells-per-quantile must be at least 1, got {args.cells_per_quantile}")
        sys.exit(1)
    
    # Validate snapshot years
    if args.snapshot_year < 0:
        print(f"Error: snapshot-year must be non-negative, got {args.snapshot_year}")
        sys.exit(1)
    
    if args.growth_years < 0:
        print(f"Error: growth-years must be non-negative, got {args.growth_years}")
        sys.exit(1)
    
    # Check if simulation is long enough (basic check - more detailed check in run_pipeline)
    second_snapshot = args.snapshot_year + args.growth_years
    # We can't check the actual simulation length here without loading it,
    # but we can warn about obviously problematic values
    if second_snapshot > 200:
        print(f"Warning: Second snapshot year {second_snapshot} seems very high. Make sure your simulation runs that long.")
    
    # Run pipeline
    run_pipeline(args)


if __name__ == "__main__":
    main()