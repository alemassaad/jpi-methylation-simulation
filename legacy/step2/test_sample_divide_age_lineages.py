#!/usr/bin/env python3
"""
Test version of sample_divide_age_lineages.py with smaller parameters.
- Only 1 cell per decile (10 total instead of 30)
- Only 3 years of aging (8 cells per lineage instead of 1024)
- Faster execution for testing
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


def sample_cells_by_deciles(snapshot_data, cells_per_decile=1):
    """
    Sample cells uniformly from each JSD decile.
    TEST VERSION: Only 1 cell per decile
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
    
    print(f"TEST MODE: Sampled {len(sampled_cells_with_info)} cells total (1 per decile)")
    return sampled_cells_with_info


def age_lineage_with_division(initial_cell_dict, years=3, rate=None):
    """
    Age a single cell lineage with division.
    TEST VERSION: Only 3 years
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
        
        print(f"    Year {year + 1}: {len(cells)} cells")
    
    # Convert back to dictionaries
    return [cell.to_dict() for cell in cells]


def save_lineage(lineage_data, output_dir="test_lineages"):
    """Save a single lineage to a compressed JSON file."""
    lineage_id = lineage_data['metadata']['lineage_id']
    decile = lineage_data['metadata']['source_decile']
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Create filename
    filename = f"lineage_{lineage_id:02d}_decile_{decile:02d}.json.gz"
    filepath = os.path.join(output_dir, filename)
    
    # Save compressed
    with gzip.open(filepath, 'wt') as f:
        json.dump(lineage_data, f, indent=2)
    
    return filepath


def verify_lineage_file(filepath):
    """Verify a lineage file was created correctly."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    metadata = data['metadata']
    cells = data['cells']
    
    print(f"\n  Verification for {os.path.basename(filepath)}:")
    print(f"    Lineage ID: {metadata['lineage_id']}")
    print(f"    Source decile: {metadata['source_decile']}")
    print(f"    Cell count: {len(cells)} (expected: {metadata['expected_cell_count']})")
    print(f"    Source cell JSD: {metadata['source_cell']['jsd']:.4f}")
    print(f"    Final cells JSD range: {min(c['jsd'] for c in cells):.4f} - {max(c['jsd'] for c in cells):.4f}")
    
    return len(cells) == metadata['expected_cell_count']


def main():
    """Main function to process all lineages in TEST MODE."""
    # TEST Parameters
    snapshot_file = "year50_snapshot.json.gz"
    cells_per_decile = 1  # TEST: Only 1 cell per decile
    years_to_age = 3      # TEST: Only 3 years (2^3 = 8 cells)
    output_dir = "test_lineages"
    
    print("="*60)
    print("RUNNING IN TEST MODE")
    print("- 1 cell per decile (10 total)")
    print("- 3 years of aging (8 cells per lineage)")
    print("="*60)
    
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
    created_files = []
    
    for lineage_id, (cell_data, decile, within_decile_idx) in enumerate(sampled_cells):
        # Progress update
        print(f"\nProcessing lineage {lineage_id + 1}/{len(sampled_cells)} "
              f"(from decile {decile})")
        
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
        created_files.append(filepath)
        
        lineage_time = time.time() - lineage_start
        print(f"  Saved {len(aged_cells)} cells to {filepath}")
        print(f"  Time: {lineage_time:.1f}s")
    
    # Verify all files
    print(f"\n{'='*60}")
    print("VERIFYING ALL FILES")
    print(f"{'='*60}")
    
    all_valid = True
    for filepath in created_files:
        valid = verify_lineage_file(filepath)
        if not valid:
            all_valid = False
    
    # Summary
    total_time = time.time() - start_time
    total_cells = len(sampled_cells) * (2**years_to_age)
    
    print(f"\n{'='*60}")
    print(f"TEST COMPLETED!")
    print(f"{'='*60}")
    print(f"Total lineages processed: {len(sampled_cells)}")
    print(f"Total cells created: {total_cells}")
    print(f"Total time: {total_time:.1f}s")
    print(f"Average time per lineage: {total_time/len(sampled_cells):.2f}s")
    print(f"Output directory: {output_dir}/")
    print(f"All files valid: {all_valid}")
    
    # List all created files
    print(f"\nCreated files:")
    for i, filepath in enumerate(created_files):
        print(f"  {i+1}. {os.path.basename(filepath)}")


if __name__ == "__main__":
    main()