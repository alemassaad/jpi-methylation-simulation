#!/usr/bin/env python3
"""
Step23 Unified Pipeline V2: Improved version with flexible quantiles and better skip logic.

Key improvements:
1. Flexible quantile system (--n-quantiles, --cells-per-quantile)
2. Better skip logic that checks ALL files, not just the first
3. Handles corrupted files gracefully
4. Supports quick testing with small datasets
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
from typing import List, Dict, Tuple

from pipeline_utils import (
    extract_snapshot, save_snapshot, load_individual, save_individual,
    sample_uniform, grow_individual_inplace,
    mix_with_year60_inplace, create_control2_individual
)
from pipeline_analysis import (
    plot_jsd_distribution, analyze_populations
)
from pipeline_checkpoint import PipelineCheckpoint


def sample_by_quantiles(cells: List[Dict], n_quantiles: int = 10, 
                        cells_per_quantile: int = 3, seed: int = 42) -> List[Dict]:
    """
    Sample cells from each quantile based on JSD values.
    
    Args:
        cells: List of cell dictionaries
        n_quantiles: Number of quantiles (e.g., 4 for quartiles, 10 for deciles)
        cells_per_quantile: Number of cells to sample from each quantile
        seed: Random seed for reproducibility
    
    Returns:
        List of sampled cells (n_quantiles * cells_per_quantile total)
    """
    random.seed(seed)
    
    # Sort cells by JSD
    sorted_cells = sorted(cells, key=lambda c: c['jsd'])
    
    # Calculate quantile boundaries
    n_cells = len(sorted_cells)
    quantile_size = n_cells // n_quantiles
    
    sampled_cells = []
    
    for q in range(n_quantiles):
        # Define quantile boundaries
        start_idx = q * quantile_size
        if q == n_quantiles - 1:
            # Last quantile includes any remainder
            end_idx = n_cells
        else:
            end_idx = (q + 1) * quantile_size
        
        # Get cells in this quantile
        quantile_cells = sorted_cells[start_idx:end_idx]
        
        # Sample from this quantile
        if len(quantile_cells) >= cells_per_quantile:
            sampled = random.sample(quantile_cells, cells_per_quantile)
        else:
            # If not enough cells, take all
            sampled = quantile_cells
        
        # Add quantile information
        for cell in sampled:
            cell = cell.copy()
            cell['source_quantile'] = q
            sampled_cells.append(cell)
    
    print(f"  Sampled {len(sampled_cells)} cells from {n_quantiles} quantiles")
    return sampled_cells


def check_individual_states(directory: str, expected_grown: int = 1024) -> Dict:
    """
    Check the state of all individuals in a directory.
    
    Returns dict with:
        - total_files: number of files
        - cell_counts: dict of filename -> cell count
        - corrupted: list of corrupted files
        - all_grown: bool, True if all valid files have expected_grown cells
        - all_mixed: bool, True if all valid files have > expected_grown cells
        - mixed_count: number of files that are mixed
        - grown_count: number of files that are grown
    """
    files = sorted(glob.glob(os.path.join(directory, "*.json.gz")))
    
    result = {
        'total_files': len(files),
        'cell_counts': {},
        'corrupted': [],
        'all_grown': False,
        'all_mixed': False,
        'mixed_count': 0,
        'grown_count': 0,
        'valid_files': []
    }
    
    if not files:
        return result
    
    for f in files:
        fname = os.path.basename(f)
        try:
            with gzip.open(f, 'rt') as file:
                data = json.load(file)
                cell_count = len(data.get('cells', []))
                result['cell_counts'][fname] = cell_count
                result['valid_files'].append(f)
                
                if cell_count == expected_grown:
                    result['grown_count'] += 1
                elif cell_count > expected_grown:
                    result['mixed_count'] += 1
                    
        except Exception as e:
            result['corrupted'].append((fname, str(e)))
            result['cell_counts'][fname] = None
    
    # Check if all valid files are in the same state
    valid_counts = [c for c in result['cell_counts'].values() if c is not None]
    if valid_counts:
        result['all_grown'] = all(c == expected_grown for c in valid_counts)
        result['all_mixed'] = all(c > expected_grown for c in valid_counts)
    
    return result


def run_pipeline(args):
    """Main pipeline orchestrator with improved skip logic."""
    
    start_time = time.time()
    
    # Calculate total expected individuals
    expected_individuals = args.n_quantiles * args.cells_per_quantile
    
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
    print("STEP23 UNIFIED PIPELINE V2")
    print("=" * 80)
    print(f"Rate: {args.rate}")
    print(f"Simulation: {args.simulation}")
    print(f"Output: {base_dir}")
    print(f"Quantiles: {args.n_quantiles} (e.g., {'quartiles' if args.n_quantiles == 4 else 'deciles' if args.n_quantiles == 10 else f'{args.n_quantiles}-tiles'})")
    print(f"Cells per quantile: {args.cells_per_quantile}")
    print(f"Total individuals: {expected_individuals} per group")
    print(f"Growth years: {args.growth_years} (target: {2**args.growth_years} cells)")
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
    checkpoint.set_parameter("n_quantiles", args.n_quantiles)
    checkpoint.set_parameter("cells_per_quantile", args.cells_per_quantile)
    checkpoint.set_parameter("expected_individuals", expected_individuals)
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
    if os.path.exists(year50_path):
        print(f"  ⏭ Snapshot already exists at {year50_path}, loading...")
        with gzip.open(year50_path, 'rt') as f:
            snapshot_data = json.load(f)
        year50_cells = snapshot_data['cells']
        print(f"  Loaded {len(year50_cells)} cells from existing snapshot")
        # Mark as complete if not already done
        if not checkpoint.is_stage_complete("extract_year50"):
            checkpoint.mark_stage_complete("extract_year50", {"cells": len(year50_cells)})
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
    
    # Check state of mutant individuals
    mutant_state = check_individual_states(mutant_dir, 2**args.growth_years)
    
    if mutant_state['total_files'] >= expected_individuals and mutant_state['corrupted']:
        print(f"\n⚠️  Warning: {len(mutant_state['corrupted'])} corrupted mutant files detected")
        for fname, error in mutant_state['corrupted'][:3]:
            print(f"    {fname}: {error[:50]}...")
    
    skip_mutant_creation = (mutant_state['total_files'] >= expected_individuals and 
                           len(mutant_state['valid_files']) == expected_individuals)
    
    if skip_mutant_creation:
        print(f"\nMutant individuals already exist ({mutant_state['total_files']} files), skipping creation...")
    else:
        print(f"\nSampling mutant cells (quantile-based: {args.n_quantiles} quantiles, {args.cells_per_quantile} cells each)...")
        mutant_cells = sample_by_quantiles(year50_cells, 
                                          n_quantiles=args.n_quantiles,
                                          cells_per_quantile=args.cells_per_quantile,
                                          seed=args.seed)
        
        # Save as single-cell individuals
        print("Saving mutant individuals...")
        for i, cell in enumerate(mutant_cells):
            filepath = os.path.join(mutant_dir, f"individual_{i:02d}.json.gz")
            metadata = {
                "individual_id": i,
                "individual_type": "mutant",
                "source_quantile": cell.get('source_quantile', 0),
                "n_quantiles": args.n_quantiles,
                "initial_year": 50
            }
            save_individual([cell], filepath, metadata)
        print(f"  Created {len(mutant_cells)} mutant individuals")
    
    # Check state of control1 individuals
    control1_state = check_individual_states(control1_dir, 2**args.growth_years)
    
    skip_control1_creation = (control1_state['total_files'] >= expected_individuals and
                             len(control1_state['valid_files']) == expected_individuals)
    
    if skip_control1_creation:
        print(f"\nControl1 individuals already exist ({control1_state['total_files']} files), skipping creation...")
    else:
        print(f"\nSampling control1 cells (uniform: {expected_individuals} cells)...")
        control1_cells = sample_uniform(year50_cells, n_cells=expected_individuals, seed=args.seed + 1)
        
        # Save as single-cell individuals
        print("Saving control1 individuals...")
        for i, cell in enumerate(control1_cells):
            filepath = os.path.join(control1_dir, f"individual_{i:02d}.json.gz")
            metadata = {
                "individual_id": i,
                "individual_type": "control1",
                "source_quantile": 0,  # Uniform sampling
                "initial_year": 50
            }
            save_individual([cell], filepath, metadata)
        print(f"  Created {len(control1_cells)} control1 individuals")
    
    # ========================================================================
    # STAGE 4: Grow individuals (with improved skip logic)
    # ========================================================================
    print(f"\n{'='*60}")
    print(f"STAGE 4: Grow Individuals ({args.growth_years} years)")
    print(f"{'='*60}")
    
    expected_cells_after_growth = 2**args.growth_years
    
    # Re-check mutant state after potential creation
    mutant_state = check_individual_states(mutant_dir, expected_cells_after_growth)
    
    if mutant_state['all_grown'] or mutant_state['all_mixed']:
        print(f"\nMutant individuals already processed:")
        print(f"  Grown (={expected_cells_after_growth} cells): {mutant_state['grown_count']}")
        print(f"  Mixed (>{expected_cells_after_growth} cells): {mutant_state['mixed_count']}")
        print(f"  Skipping growth stage for mutants...")
    else:
        print(f"\nGrowing mutant individuals:")
        print(f"  Current state: {mutant_state['grown_count']} grown, {mutant_state['mixed_count']} mixed")
        
        for filepath in mutant_state['valid_files']:
            fname = os.path.basename(filepath)
            current_cells = mutant_state['cell_counts'][fname]
            
            if current_cells and current_cells < expected_cells_after_growth:
                print(f"  Growing {fname}: {current_cells} → {expected_cells_after_growth} cells")
                final_count = grow_individual_inplace(filepath, args.growth_years, rate=args.rate)
                if final_count != expected_cells_after_growth:
                    print(f"    WARNING: Expected {expected_cells_after_growth} cells, got {final_count}")
            elif current_cells == expected_cells_after_growth:
                print(f"  Skipping {fname}: already at {current_cells} cells")
            elif current_cells and current_cells > expected_cells_after_growth:
                print(f"  Skipping {fname}: already mixed ({current_cells} cells)")
    
    # Re-check control1 state
    control1_state = check_individual_states(control1_dir, expected_cells_after_growth)
    
    if control1_state['all_grown'] or control1_state['all_mixed']:
        print(f"\nControl1 individuals already processed:")
        print(f"  Grown (={expected_cells_after_growth} cells): {control1_state['grown_count']}")
        print(f"  Mixed (>{expected_cells_after_growth} cells): {control1_state['mixed_count']}")
        print(f"  Skipping growth stage for control1...")
    else:
        print(f"\nGrowing control1 individuals:")
        print(f"  Current state: {control1_state['grown_count']} grown, {control1_state['mixed_count']} mixed")
        
        for filepath in control1_state['valid_files']:
            fname = os.path.basename(filepath)
            current_cells = control1_state['cell_counts'][fname]
            
            if current_cells and current_cells < expected_cells_after_growth:
                print(f"  Growing {fname}: {current_cells} → {expected_cells_after_growth} cells")
                final_count = grow_individual_inplace(filepath, args.growth_years, rate=args.rate)
                if final_count != expected_cells_after_growth:
                    print(f"    WARNING: Expected {expected_cells_after_growth} cells, got {final_count}")
    
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
    # STAGE 6: Mix populations (with improved skip logic)
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 6: Mix Populations")
    print(f"{'='*60}")
    
    # Calculate expected final size after mixing
    grown_cells = expected_cells_after_growth
    expected_final_cells = int(grown_cells / ((100 - args.mix_ratio) / 100))
    
    # Re-check states after growth
    mutant_state = check_individual_states(mutant_dir, grown_cells)
    
    if mutant_state['all_mixed']:
        print(f"\nMutant individuals already mixed ({mutant_state['mixed_count']} files), skipping mixing...")
    else:
        print(f"\nMixing mutant individuals with year 60 cells:")
        print(f"  Target: {expected_final_cells} cells per individual")
        
        for i, filepath in enumerate(mutant_state['valid_files']):
            fname = os.path.basename(filepath)
            current_cells = mutant_state['cell_counts'][fname]
            
            if current_cells == grown_cells:
                print(f"  Mixing {fname}: {current_cells} → {expected_final_cells} cells")
                total_cells = mix_with_year60_inplace(filepath, year60_cells,
                                                     mix_ratio=args.mix_ratio,
                                                     seed=args.seed + 100 + i)
            elif current_cells and current_cells > grown_cells:
                print(f"  Skipping {fname}: already mixed ({current_cells} cells)")
    
    control1_state = check_individual_states(control1_dir, grown_cells)
    
    if control1_state['all_mixed']:
        print(f"\nControl1 individuals already mixed ({control1_state['mixed_count']} files), skipping mixing...")
    else:
        print(f"\nMixing control1 individuals with year 60 cells:")
        
        for i, filepath in enumerate(control1_state['valid_files']):
            fname = os.path.basename(filepath)
            current_cells = control1_state['cell_counts'][fname]
            
            if current_cells == grown_cells:
                print(f"  Mixing {fname}: {current_cells} → {expected_final_cells} cells")
                total_cells = mix_with_year60_inplace(filepath, year60_cells,
                                                     mix_ratio=args.mix_ratio,
                                                     seed=args.seed + 200 + i)
            elif current_cells and current_cells > grown_cells:
                print(f"  Skipping {fname}: already mixed ({current_cells} cells)")
    
    # ========================================================================
    # STAGE 7: Create control2 individuals
    # ========================================================================
    print(f"\n{'='*60}")
    print("STAGE 7: Create Control2 Individuals")
    print(f"{'='*60}")
    
    control2_state = check_individual_states(control2_dir, grown_cells)
    skip_control2 = control2_state['total_files'] >= expected_individuals
    
    if skip_control2:
        print(f"\nControl2 individuals already exist ({control2_state['total_files']} files), skipping creation...")
    else:
        # Calculate target cell count (same as mixed individuals)
        print(f"\nCreating {expected_individuals} control2 individuals ({expected_final_cells} cells each)...")
        for i in range(expected_individuals):
            print(f"  Individual {i+1}/{expected_individuals}")
            cells = create_control2_individual(year60_cells, expected_final_cells,
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
        "pipeline_version": "v2",
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
        "results_summary": {
            "mutant_mean_jsd": stats['mutant']['mean'],
            "control1_mean_jsd": stats['control1']['mean'],
            "control2_mean_jsd": stats['control2']['mean']
        }
    }
    
    metadata_path = os.path.join(results_dir, "pipeline_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    print(f"\nPipeline metadata saved to: {metadata_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Step23 Unified Pipeline V2 with flexible quantile sampling",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--rate", type=float, required=True,
                       help="Methylation rate (must match simulation)")
    parser.add_argument("--simulation", type=str, required=True,
                       help="Path to step1 simulation file")
    
    # Quantile sampling parameters (NEW)
    parser.add_argument("--n-quantiles", type=int, default=10,
                       help="Number of quantiles for sampling (e.g., 4 for quartiles, 10 for deciles)")
    parser.add_argument("--cells-per-quantile", type=int, default=3,
                       help="Number of cells to sample per quantile")
    
    # Growth and mixing parameters
    parser.add_argument("--growth-years", type=int, default=10,
                       help="Number of years to grow individuals")
    parser.add_argument("--mix-ratio", type=int, default=80,
                       help="Percentage of year 60 cells in mix (0-100)")
    
    # Visualization parameters
    parser.add_argument("--bins", type=int, default=200,
                       help="Number of bins for JSD histograms")
    
    # Other parameters
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducibility")
    parser.add_argument("--output-dir", type=str, default="data",
                       help="Output directory")
    
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
    
    # Run pipeline
    run_pipeline(args)


if __name__ == "__main__":
    main()