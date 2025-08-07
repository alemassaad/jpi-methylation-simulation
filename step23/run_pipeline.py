#!/usr/bin/env python3
"""
Step23 Unified Pipeline: Combines steps 2 and 3 into a single streamlined process.

This pipeline:
1. Extracts year 50 snapshot
2. Plots JSD distribution
3. Creates initial individuals (single cells)
4. Grows individuals in-place for 10 years
5. Extracts year 60 snapshot
6. Mixes populations
7. Creates control2 group
8. Analyzes and compares all three groups
"""

import argparse
import os
import sys
import time
import json
import gzip
import glob
import random
import numpy as np

from pipeline_utils import (
    extract_snapshot, save_snapshot, load_individual, save_individual,
    sample_by_deciles, sample_uniform, grow_individual_inplace,
    mix_with_year60_inplace, create_control2_individual
)
from pipeline_analysis import (
    plot_jsd_distribution, analyze_populations
)
from pipeline_checkpoint import PipelineCheckpoint


def run_pipeline(args):
    """Main pipeline orchestrator."""
    
    start_time = time.time()
    
    # Setup output directories
    rate_str = f"rate_{args.rate:.6f}"
    base_dir = os.path.join(args.output_dir, rate_str)
    
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
    print("STEP23 UNIFIED PIPELINE")
    print("=" * 80)
    print(f"Rate: {args.rate}")
    print(f"Simulation: {args.simulation}")
    print(f"Output: {base_dir}")
    print(f"Mix ratio: {args.mix_ratio}% year 60, {100-args.mix_ratio}% grown")
    print(f"Seed: {args.seed}")
    print("=" * 80)
    
    # Set random seed for reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    # Initialize checkpoint system
    checkpoint_file = os.path.join(base_dir, "pipeline_checkpoint.json")
    checkpoint = PipelineCheckpoint(checkpoint_file)
    
    # Store pipeline parameters
    checkpoint.set_parameter("rate", args.rate)
    checkpoint.set_parameter("n_individuals", args.n_individuals)
    checkpoint.set_parameter("growth_years", args.growth_years)
    checkpoint.set_parameter("mix_ratio", args.mix_ratio)
    checkpoint.set_parameter("seed", args.seed)
    
    # ========================================================================
    # STAGE 1: Extract year 50 snapshot
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 1: Extract Year 50 Snapshot")
    print(f"{'='*60}")
    
    year50_path = os.path.join(snapshots_dir, "year50_snapshot.json.gz")
    if os.path.exists(year50_path) and checkpoint.is_stage_complete("extract_year50"):
        print(f"  ⏭ Snapshot already exists at {year50_path}, loading...")
        with gzip.open(year50_path, 'rt') as f:
            snapshot_data = json.load(f)
        year50_cells = snapshot_data['cells']
        print(f"  Loaded {len(year50_cells)} cells from existing snapshot")
    else:
        print(f"  ✓ Extracting year 50 snapshot...")
        year50_cells = extract_snapshot(args.simulation, 50)
        save_snapshot(year50_cells, year50_path)
        checkpoint.mark_stage_complete("extract_year50", {"cells": len(year50_cells)})
    
    # ========================================================================
    # STAGE 2: Plot JSD distribution
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 2: Plot JSD Distribution")
    print(f"{'='*60}")
    
    plot_path = os.path.join(plots_dir, f"year50_jsd_distribution_{args.bins}bins.png")
    plot_jsd_distribution(year50_cells, args.bins, plot_path, rate=args.rate)
    
    # ========================================================================
    # STAGE 3: Create initial individuals
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 3: Create Initial Individuals")
    print(f"{'='*60}")
    
    # Check if mutant individuals already exist
    mutant_files = sorted(glob.glob(os.path.join(mutant_dir, "individual_*.json.gz")))
    skip_mutant = len(mutant_files) == args.n_individuals
    
    if skip_mutant:
        print(f"\nMutant individuals already exist ({len(mutant_files)} files), skipping creation...")
    else:
        print("\nSampling mutant cells (decile-based)...")
        mutant_cells = sample_by_deciles(year50_cells, cells_per_decile=3, seed=args.seed)
        
        # Save as single-cell individuals
        print("Saving mutant individuals...")
        for i, cell in enumerate(mutant_cells):
            filepath = os.path.join(mutant_dir, f"individual_{i:02d}.json.gz")
            metadata = {
                "individual_id": i,
                "individual_type": "mutant",
                "source_decile": cell.get('source_decile', 0),
                "initial_year": 50
            }
            save_individual([cell], filepath, metadata)
        print(f"  Created {len(mutant_cells)} mutant individuals")
    
    # Check if control1 individuals already exist
    control1_files = sorted(glob.glob(os.path.join(control1_dir, "individual_*.json.gz")))
    skip_control1 = len(control1_files) == args.n_individuals
    
    if skip_control1:
        print(f"\nControl1 individuals already exist ({len(control1_files)} files), skipping creation...")
    else:
        print("\nSampling control1 cells (uniform)...")
        control1_cells = sample_uniform(year50_cells, n_cells=args.n_individuals, seed=args.seed + 1)
        
        # Save as single-cell individuals
        print("Saving control1 individuals...")
        for i, cell in enumerate(control1_cells):
            filepath = os.path.join(control1_dir, f"individual_{i:02d}.json.gz")
            metadata = {
                "individual_id": i,
                "individual_type": "control1",
                "source_decile": 0,  # Uniform sampling
                "initial_year": 50
            }
            save_individual([cell], filepath, metadata)
        print(f"  Created {len(control1_cells)} control1 individuals")
    
    # ========================================================================
    # STAGE 4: Grow individuals
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 4: Grow Individuals (10 years)")
    print(f"{'='*60}")
    
    # Always get current file lists (they may exist from previous runs)
    mutant_files = sorted(glob.glob(os.path.join(mutant_dir, "*.json.gz")))
    control1_files = sorted(glob.glob(os.path.join(control1_dir, "*.json.gz")))
    
    # Check if mutant individuals need growing
    skip_mutant_growth = False
    if mutant_files:
        with gzip.open(mutant_files[0], 'rt') as f:
            sample_data = json.load(f)
            current_cells = len(sample_data['cells'])
            expected_cells = 2**args.growth_years
            skip_mutant_growth = (current_cells == expected_cells)
    
    if skip_mutant_growth:
        print(f"\nMutant individuals already grown ({current_cells} cells each), skipping growth...")
    else:
        print("\nGrowing mutant individuals...")
        for i, filepath in enumerate(mutant_files):
            print(f"  Individual {i+1}/{len(mutant_files)}")
            final_count = grow_individual_inplace(filepath, args.growth_years, rate=args.rate)
            if final_count != 2**args.growth_years:
                print(f"    WARNING: Expected {2**args.growth_years} cells, got {final_count}")
    
    # Check if control1 individuals need growing
    skip_control1_growth = False
    if control1_files:
        with gzip.open(control1_files[0], 'rt') as f:
            sample_data = json.load(f)
            current_cells = len(sample_data['cells'])
            expected_cells = 2**args.growth_years
            skip_control1_growth = (current_cells == expected_cells)
    
    if skip_control1_growth:
        print(f"\nControl1 individuals already grown ({current_cells} cells each), skipping growth...")
    else:
        print("\nGrowing control1 individuals...")
        for i, filepath in enumerate(control1_files):
            print(f"  Individual {i+1}/{len(control1_files)}")
            final_count = grow_individual_inplace(filepath, args.growth_years, rate=args.rate)
            if final_count != 2**args.growth_years:
                print(f"    WARNING: Expected {2**args.growth_years} cells, got {final_count}")
    
    # ========================================================================
    # STAGE 5: Extract year 60 snapshot
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 5: Extract Year 60 Snapshot")
    print(f"{'='*60}")
    
    year60_path = os.path.join(snapshots_dir, "year60_snapshot.json.gz")
    if os.path.exists(year60_path):
        print(f"  Snapshot already exists at {year60_path}, loading...")
        with gzip.open(year60_path, 'rt') as f:
            snapshot_data = json.load(f)
        year60_cells = snapshot_data['cells']
        print(f"  Loaded {len(year60_cells)} cells from existing snapshot")
    else:
        year60_cells = extract_snapshot(args.simulation, 60)
        save_snapshot(year60_cells, year60_path)
    
    # Also plot year 60 distribution
    plot_path = os.path.join(plots_dir, f"year60_jsd_distribution_{args.bins}bins.png")
    plot_jsd_distribution(year60_cells, args.bins, plot_path, rate=args.rate)
    
    # ========================================================================
    # STAGE 6: Mix populations
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 6: Mix Populations")
    print(f"{'='*60}")
    
    # Calculate expected final size after mixing
    grown_cells = 2**args.growth_years
    expected_final_cells = int(grown_cells / ((100 - args.mix_ratio) / 100))
    
    # Check if mutant individuals need mixing
    skip_mutant_mixing = False
    if mutant_files:
        with gzip.open(mutant_files[0], 'rt') as f:
            sample_data = json.load(f)
            current_cells = len(sample_data['cells'])
            # Skip if already mixed (has more cells than just grown)
            skip_mutant_mixing = (current_cells > grown_cells)
    
    if skip_mutant_mixing:
        print(f"\nMutant individuals already mixed ({current_cells} cells each), skipping mixing...")
    else:
        print("\nMixing mutant individuals with year 60 cells...")
        for i, filepath in enumerate(mutant_files):
            print(f"  Individual {i+1}/{len(mutant_files)}")
            total_cells = mix_with_year60_inplace(filepath, year60_cells, 
                                                 mix_ratio=args.mix_ratio, 
                                                 seed=args.seed + 100 + i)
    
    # Check if control1 individuals need mixing
    skip_control1_mixing = False
    if control1_files:
        with gzip.open(control1_files[0], 'rt') as f:
            sample_data = json.load(f)
            current_cells = len(sample_data['cells'])
            # Skip if already mixed (has more cells than just grown)
            skip_control1_mixing = (current_cells > grown_cells)
    
    if skip_control1_mixing:
        print(f"\nControl1 individuals already mixed ({current_cells} cells each), skipping mixing...")
    else:
        print("\nMixing control1 individuals with year 60 cells...")
        for i, filepath in enumerate(control1_files):
            print(f"  Individual {i+1}/{len(control1_files)}")
            total_cells = mix_with_year60_inplace(filepath, year60_cells, 
                                                 mix_ratio=args.mix_ratio,
                                                 seed=args.seed + 200 + i)
    
    # ========================================================================
    # STAGE 7: Create control2 individuals
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 7: Create Control2 Individuals")
    print(f"{'='*60}")
    
    # Check if control2 individuals already exist
    control2_files = sorted(glob.glob(os.path.join(control2_dir, "individual_*.json.gz")))
    skip_control2 = len(control2_files) == args.n_individuals
    
    if skip_control2:
        print(f"\nControl2 individuals already exist ({len(control2_files)} files), skipping creation...")
    else:
        # Calculate target cell count (same as mixed individuals)
        with gzip.open(mutant_files[0], 'rt') as gf:
            sample_data = json.load(gf)
            target_cell_count = len(sample_data['cells'])
        
        print(f"\nCreating {args.n_individuals} control2 individuals ({target_cell_count} cells each)...")
        for i in range(args.n_individuals):
            print(f"  Individual {i+1}/{args.n_individuals}")
            cells = create_control2_individual(year60_cells, target_cell_count, 
                                              seed=args.seed + 300 + i)
            filepath = os.path.join(control2_dir, f"individual_{i:02d}.json.gz")
            metadata = {
                "individual_id": i,
                "individual_type": "control2",
                "source": "pure_year60",
                "year": 60
            }
            save_individual(cells, filepath, metadata)
    
    # ========================================================================
    # STAGE 8: Analysis
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 8: Analysis and Visualization")
    print(f"{'='*60}")
    
    results = analyze_populations(mutant_dir, control1_dir, control2_dir, results_dir)
    
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
    stats = results['statistics']
    print(f"  Mutant mean JSD: {stats['mutant']['mean']:.6f} ± {stats['mutant']['std']:.6f}")
    print(f"  Control1 mean JSD: {stats['control1']['mean']:.6f} ± {stats['control1']['std']:.6f}")
    print(f"  Control2 mean JSD: {stats['control2']['mean']:.6f} ± {stats['control2']['std']:.6f}")
    
    if 'comparisons' in stats:
        print("\nStatistical tests:")
        for comparison, values in stats['comparisons'].items():
            print(f"  {comparison}: p={values['p_value']:.6f}")
    
    # Save pipeline metadata
    metadata = {
        "rate": args.rate,
        "simulation_file": args.simulation,
        "mix_ratio": args.mix_ratio,
        "n_individuals": args.n_individuals,
        "growth_years": args.growth_years,
        "bins": args.bins,
        "seed": args.seed,
        "elapsed_time_seconds": elapsed_time,
        "output_directory": base_dir
    }
    metadata_path = os.path.join(results_dir, "pipeline_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\nPipeline metadata saved to {metadata_path}")
    print("=" * 80)


def main():
    parser = argparse.ArgumentParser(
        description="Step23 Unified Pipeline: Streamlined cell sampling, growth, and analysis",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument('--rate', type=float, required=True,
                       help='Methylation rate (must match simulation file)')
    parser.add_argument('--simulation', type=str, required=True,
                       help='Path to simulation file from step1')
    
    # Optional arguments
    parser.add_argument('--bins', type=int, default=200,
                       help='Number of bins for JSD distribution plots')
    parser.add_argument('--mix-ratio', type=int, default=80,
                       help='Percentage of year 60 cells in final mixture')
    parser.add_argument('--n-individuals', type=int, default=30,
                       help='Number of individuals per group')
    parser.add_argument('--growth-years', type=int, default=10,
                       help='Years to grow individuals (50->60)')
    parser.add_argument('--seed', type=int, default=42,
                       help='Random seed for reproducibility')
    parser.add_argument('--output-dir', type=str, default='data',
                       help='Base output directory')
    
    args = parser.parse_args()
    
    # Validate inputs
    if not os.path.exists(args.simulation):
        print(f"Error: Simulation file not found: {args.simulation}")
        sys.exit(1)
    
    if args.mix_ratio < 0 or args.mix_ratio > 100:
        print(f"Error: mix-ratio must be between 0 and 100")
        sys.exit(1)
    
    # Run pipeline
    try:
        run_pipeline(args)
    except Exception as e:
        print(f"\n{'='*80}")
        print(f"ERROR: Pipeline failed")
        print(f"{'='*80}")
        print(f"{str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()