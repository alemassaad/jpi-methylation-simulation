#!/usr/bin/env python3
"""
Create cell lineages with division from a snapshot, supporting both mutant and control types.

Mutant lineages: Sample 3 cells from each JSD decile (30 total), age with division
Control lineages: Sample 30 cells uniformly from population, age with division

Each year: cells divide (creating identical copies), then both daughters methylate.
After 10 years, each lineage contains 1,024 cells (2^10).
"""

import argparse
import json
import gzip
import numpy as np
import random
import copy
import time
import os
import sys
sys.path.append('../..')  # To import from parent directory

from step1.cell import Cell, GENE_SIZE, BASELINE_METHYLATION_DISTRIBUTION


def load_snapshot(filepath):
    """Load snapshot data from json.gz file."""
    print(f"Loading snapshot from {filepath}...")
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    print(f"  Loaded {data['metadata']['num_cells']} cells from year {data['metadata']['extracted_year']}")
    return data


def sample_cells_by_deciles(snapshot_data, cells_per_decile=3):
    """
    Sample cells uniformly from each JSD decile.
    
    Args:
        snapshot_data: The snapshot data containing all cells
        cells_per_decile: Number of cells to sample from each decile (default: 3)
    
    Returns:
        List of tuples: (cell_data, decile_number, within_decile_index)
    """
    # Extract all cells with their indices
    all_cells = snapshot_data['cells']
    cells_with_jsd = [(i, cell, cell['jsd']) for i, cell in enumerate(all_cells)]
    
    # Sort by JSD
    cells_with_jsd.sort(key=lambda x: x[2])
    
    # Calculate decile boundaries
    n_cells = len(cells_with_jsd)
    decile_size = n_cells // 10
    
    sampled_cells_with_info = []
    
    for decile in range(10):
        # Define the range for this decile
        start_idx = decile * decile_size
        if decile == 9:  # Last decile includes any remainder
            end_idx = n_cells
        else:
            end_idx = (decile + 1) * decile_size
        
        # Get cells in this decile
        decile_cells = cells_with_jsd[start_idx:end_idx]
        
        # Sample uniformly from this decile
        if len(decile_cells) >= cells_per_decile:
            sampled_indices = random.sample(range(len(decile_cells)), cells_per_decile)
            for idx, sample_idx in enumerate(sampled_indices):
                original_idx, cell_data, jsd = decile_cells[sample_idx]
                sampled_cells_with_info.append((cell_data, decile + 1, idx))  # decile is 1-indexed
        else:
            print(f"Warning: Decile {decile+1} has only {len(decile_cells)} cells, sampling all")
            for idx, (original_idx, cell_data, jsd) in enumerate(decile_cells):
                sampled_cells_with_info.append((cell_data, decile + 1, idx))
    
    print(f"Sampled {len(sampled_cells_with_info)} cells total across {10} deciles")
    return sampled_cells_with_info


def sample_cells_uniformly(snapshot_data, n_cells=30):
    """
    Sample cells uniformly from the entire population.
    
    Args:
        snapshot_data: The snapshot data containing all cells
        n_cells: Number of cells to sample uniformly (default: 30)
    
    Returns:
        List of tuples: (cell_data, 0, index) - decile is 0 for uniform sampling
    """
    all_cells = snapshot_data['cells']
    
    if n_cells > len(all_cells):
        print(f"Warning: Requested {n_cells} cells but only {len(all_cells)} available, sampling all")
        n_cells = len(all_cells)
    
    # Sample uniformly from all cells
    sampled_indices = random.sample(range(len(all_cells)), n_cells)
    sampled_cells_with_info = []
    
    for idx, cell_idx in enumerate(sampled_indices):
        cell_data = all_cells[cell_idx]
        sampled_cells_with_info.append((cell_data, 0, idx))  # decile=0 for uniform sampling
    
    print(f"Sampled {len(sampled_cells_with_info)} cells uniformly from entire population")
    return sampled_cells_with_info


def age_lineage_with_division(initial_cell_dict, years=10, rate=None):
    """
    Age a single cell lineage with division.
    Each year: divide (double the cells), then age all cells.
    
    Args:
        initial_cell_dict: Dictionary representation of the initial cell
        years: Number of years to age
        rate: Methylation rate (if None, uses cell's existing rate)
    
    Returns:
        List of cell dictionaries after aging
    """
    # Convert dict to Cell object
    initial_cell = Cell(
        n=len(initial_cell_dict['cpg_sites']),
        rate=initial_cell_dict['rate'] if rate is None else rate,
        gene_size=initial_cell_dict.get('gene_size', GENE_SIZE)
    )
    initial_cell.cpg_sites = initial_cell_dict['cpg_sites'].copy()
    initial_cell.age = initial_cell_dict['age']
    # Set the existing computed values
    initial_cell.methylation_proportion = initial_cell_dict['methylation_proportion']
    initial_cell.methylation_distribution = initial_cell_dict['methylation_distribution'].copy()
    initial_cell.JSD = initial_cell_dict['jsd']
    
    # Start with one cell
    cells = [initial_cell]
    
    # Age with division
    for year in range(years):
        # Division: each cell creates an identical copy
        new_cells = []
        for cell in cells:
            # Create daughter cells (identical copies)
            daughter1 = copy.deepcopy(cell)
            daughter2 = copy.deepcopy(cell)
            new_cells.extend([daughter1, daughter2])
        
        cells = new_cells
        
        # Age all cells
        for cell in cells:
            cell.age_1_year()
    
    # Convert back to dictionaries
    return [cell.to_dict() for cell in cells]


def save_lineage(lineage_data, output_dir="lineages", lineage_type="mutant"):
    """Save a single lineage to a compressed JSON file."""
    lineage_id = lineage_data['metadata']['lineage_id']
    
    # Create filename based on type
    if lineage_type == "control":
        filename = f"control_{lineage_id:02d}.json.gz"
    else:  # mutant
        decile = lineage_data['metadata']['source_decile']
        filename = f"lineage_{lineage_id:02d}_decile_{decile:02d}.json.gz"
    
    filepath = os.path.join(output_dir, filename)
    
    # Save compressed
    with gzip.open(filepath, 'wt') as f:
        json.dump(lineage_data, f, indent=2)
    
    return filepath


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Create cell lineages with division from a snapshot",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('snapshot_file', 
                        help='Input snapshot file (e.g., data/snapshots/year50_snapshot.json.gz)')
    parser.add_argument('-t', '--type', choices=['mutant', 'control', 'both'], default='both',
                        help='Type of lineages to create')
    parser.add_argument('-c', '--cells-per-decile', type=int, default=3,
                        help='Number of cells to sample from each JSD decile (for mutant)')
    parser.add_argument('-n', '--n-control', type=int, default=30,
                        help='Number of control cells to sample uniformly')
    parser.add_argument('-y', '--years', type=int, default=10,
                        help='Number of years to age cells with division')
    parser.add_argument('-o', '--output-dir', default='../data/lineages',
                        help='Output directory for lineage files')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')
    
    return parser.parse_args()


def process_lineages(sampled_cells, snapshot_data, args, lineage_type="mutant"):
    """Process a set of lineages (either mutant or control)."""
    output_dir = os.path.join(args.output_dir, lineage_type)
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\nProcessing {len(sampled_cells)} {lineage_type} lineages...")
    print(f"Each lineage will have {2**args.years} cells after {args.years} years of division")
    print(f"Output directory: {output_dir}/")
    
    start_time = time.time()
    
    for lineage_id, (cell_data, decile, within_idx) in enumerate(sampled_cells):
        # Progress update
        if lineage_type == "mutant":
            print(f"\nProcessing {lineage_type} lineage {lineage_id + 1}/{len(sampled_cells)} "
                  f"(from decile {decile}, cell {within_idx + 1})")
        else:
            print(f"\nProcessing {lineage_type} lineage {lineage_id + 1}/{len(sampled_cells)}")
        
        lineage_start = time.time()
        
        # Age this lineage with division
        aged_cells = age_lineage_with_division(cell_data, years=args.years)
        
        # Create lineage data structure
        lineage_data = {
            "metadata": {
                "lineage_id": lineage_id,
                "lineage_type": lineage_type,
                "source_decile": decile,
                "within_decile_index": within_idx,
                "source_year": snapshot_data['metadata']['extracted_year'],
                "final_year": snapshot_data['metadata']['extracted_year'] + args.years,
                "generations": args.years,
                "final_cell_count": len(aged_cells),
                "expected_cell_count": 2**args.years,
                "source_cell": copy.deepcopy(cell_data)  # Store original cell
            },
            "cells": aged_cells
        }
        
        # Save lineage
        filepath = save_lineage(lineage_data, output_dir, lineage_type)
        
        lineage_time = time.time() - lineage_start
        print(f"  Saved {len(aged_cells)} cells to {filepath}")
        print(f"  Time: {lineage_time:.1f}s")
    
    # Summary for this type
    total_time = time.time() - start_time
    total_cells = len(sampled_cells) * (2**args.years)
    
    print(f"\n{lineage_type.upper()} SUMMARY:")
    print(f"  Lineages processed: {len(sampled_cells)}")
    print(f"  Total cells created: {total_cells:,}")
    print(f"  Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"  Average time per lineage: {total_time/len(sampled_cells):.1f}s")
    
    return total_time, total_cells


def main():
    """Main function to process all lineages."""
    # Parse arguments
    args = parse_arguments()
    
    # Set random seed for reproducibility
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    # Load snapshot
    snapshot_data = load_snapshot(args.snapshot_file)
    
    print("=" * 60)
    print("CELL LINEAGE CREATION WITH DIVISION")
    print("=" * 60)
    
    grand_start_time = time.time()
    total_stats = {"cells": 0, "time": 0, "lineages": 0}
    
    # Process mutant lineages if requested
    if args.type in ['mutant', 'both']:
        print("\n### CREATING MUTANT LINEAGES ###")
        sampled_cells = sample_cells_by_deciles(snapshot_data, args.cells_per_decile)
        time_taken, cells_created = process_lineages(sampled_cells, snapshot_data, args, "mutant")
        total_stats["cells"] += cells_created
        total_stats["time"] += time_taken
        total_stats["lineages"] += len(sampled_cells)
    
    # Process control lineages if requested
    if args.type in ['control', 'both']:
        print("\n### CREATING CONTROL LINEAGES ###")
        sampled_cells = sample_cells_uniformly(snapshot_data, args.n_control)
        time_taken, cells_created = process_lineages(sampled_cells, snapshot_data, args, "control")
        total_stats["cells"] += cells_created
        total_stats["time"] += time_taken
        total_stats["lineages"] += len(sampled_cells)
    
    # Grand summary
    grand_total_time = time.time() - grand_start_time
    print("\n" + "=" * 60)
    print("GRAND SUMMARY")
    print("=" * 60)
    print(f"Total lineages created: {total_stats['lineages']}")
    print(f"Total cells created: {total_stats['cells']:,}")
    print(f"Total processing time: {grand_total_time:.1f}s ({grand_total_time/60:.1f} minutes)")
    print(f"Output directory: {args.output_dir}/")
    
    if args.type == 'both':
        print(f"  - Mutant lineages: {args.output_dir}/mutant/")
        print(f"  - Control lineages: {args.output_dir}/control/")


if __name__ == "__main__":
    main()