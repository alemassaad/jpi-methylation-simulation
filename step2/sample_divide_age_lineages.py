#!/usr/bin/env python3
"""
Sample cells from JSD deciles, then age each cell separately with division to track lineages.
Process: Sample 3 cells per decile → Age each separately for 10 years → Save 30 lineage files.
Each year: cells divide (creating identical copies), then both daughters methylate.
"""

import json
import gzip
import numpy as np
import random
import copy
import time
import os
import sys
sys.path.append('..')  # To import from parent directory

from cell import Cell, GENE_SIZE, BASELINE_METHYLATION_DISTRIBUTION


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


def save_lineage(lineage_data, output_dir="lineages"):
    """Save a single lineage to a compressed JSON file."""
    lineage_id = lineage_data['metadata']['lineage_id']
    decile = lineage_data['metadata']['source_decile']
    
    # Create filename
    filename = f"lineage_{lineage_id:02d}_decile_{decile:02d}.json.gz"
    filepath = os.path.join(output_dir, filename)
    
    # Save compressed
    with gzip.open(filepath, 'wt') as f:
        json.dump(lineage_data, f, indent=2)
    
    return filepath


def main():
    """Main function to process all lineages."""
    # Parameters
    snapshot_file = "year50_snapshot.json.gz"
    cells_per_decile = 3
    years_to_age = 10
    output_dir = "lineages"
    
    # Set random seed for reproducibility
    random.seed(42)
    np.random.seed(42)
    
    # Load snapshot
    snapshot_data = load_snapshot(snapshot_file)
    
    # Sample cells by deciles
    sampled_cells = sample_cells_by_deciles(snapshot_data, cells_per_decile)
    
    print(f"\nAging {len(sampled_cells)} lineages separately...")
    print(f"Each lineage will have {2**years_to_age} cells after {years_to_age} years of division")
    
    # Process each lineage
    start_time = time.time()
    
    for lineage_id, (cell_data, decile, within_decile_idx) in enumerate(sampled_cells):
        # Progress update
        print(f"\nProcessing lineage {lineage_id + 1}/{len(sampled_cells)} "
              f"(from decile {decile}, cell {within_decile_idx + 1}/3)")
        
        lineage_start = time.time()
        
        # Age this lineage with division
        aged_cells = age_lineage_with_division(cell_data, years=years_to_age)
        
        # Create lineage data structure
        lineage_data = {
            "metadata": {
                "lineage_id": lineage_id,
                "source_decile": decile,
                "within_decile_index": within_decile_idx,
                "source_year": snapshot_data['metadata']['extracted_year'],
                "final_year": snapshot_data['metadata']['extracted_year'] + years_to_age,
                "generations": years_to_age,
                "final_cell_count": len(aged_cells),
                "expected_cell_count": 2**years_to_age,
                "source_cell": copy.deepcopy(cell_data)  # Store original cell
            },
            "cells": aged_cells
        }
        
        # Save lineage
        filepath = save_lineage(lineage_data, output_dir)
        
        lineage_time = time.time() - lineage_start
        print(f"  Saved {len(aged_cells)} cells to {filepath}")
        print(f"  Time: {lineage_time:.1f}s")
    
    # Summary
    total_time = time.time() - start_time
    total_cells = len(sampled_cells) * (2**years_to_age)
    
    print(f"\n{'='*60}")
    print(f"COMPLETED!")
    print(f"{'='*60}")
    print(f"Total lineages processed: {len(sampled_cells)}")
    print(f"Total cells created: {total_cells:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"Average time per lineage: {total_time/len(sampled_cells):.1f}s")
    print(f"Output directory: {output_dir}/")


if __name__ == "__main__":
    main()