#!/usr/bin/env python3
"""
Sample cells from JSD deciles, then age them with cell division.
Process: Sample 3 cells per decile → Age 10 years with division → Save results.
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
        List of sampled cells (30 cells total if cells_per_decile=3)
    """
    # Extract all cells with their indices
    all_cells = snapshot_data['cells']
    cells_with_jsd = [(i, cell, cell['jsd']) for i, cell in enumerate(all_cells)]
    
    # Sort by JSD
    cells_with_jsd.sort(key=lambda x: x[2])
    
    # Calculate decile boundaries
    n_cells = len(cells_with_jsd)
    decile_size = n_cells // 10
    
    sampled_cells = []
    
    for decile in range(10):
        # Define decile boundaries
        start_idx = decile * decile_size
        if decile == 9:  # Last decile includes any remainder
            end_idx = n_cells
        else:
            end_idx = (decile + 1) * decile_size
        
        # Get cells in this decile
        decile_cells = cells_with_jsd[start_idx:end_idx]
        
        # Sample cells_per_decile cells randomly from this decile
        if len(decile_cells) >= cells_per_decile:
            sampled_indices = random.sample(range(len(decile_cells)), cells_per_decile)
            for idx in sampled_indices:
                original_idx, cell, jsd = decile_cells[idx]
                sampled_cells.append((original_idx, cell, jsd, decile + 1))
        
        print(f"  Decile {decile + 1}: JSD range [{decile_cells[0][2]:.4f}, {decile_cells[-1][2]:.4f}], "
              f"sampled {len(sampled_indices)} cells")
    
    return sampled_cells


def create_cell_from_dict(cell_dict):
    """Create a Cell object from dictionary representation."""
    # Create new cell with same parameters
    cell = Cell(
        n=len(cell_dict['cpg_sites']),
        rate=cell_dict['rate'],
        gene_size=cell_dict['gene_size'],
        baseline_methylation_distribution=BASELINE_METHYLATION_DISTRIBUTION
    )
    
    # Set the state from the dictionary
    cell.cpg_sites = copy.deepcopy(cell_dict['cpg_sites'])
    cell.age = cell_dict['age']
    
    # Recompute derived properties
    cell.compute_methylation_distribution()
    
    return cell


def age_cells_with_division(initial_cells, years_to_age):
    """
    Age cells for specified years with division.
    Each year: divide (copy), then both daughters methylate.
    
    Args:
        initial_cells: List of (original_idx, cell_dict, jsd, decile) tuples
        years_to_age: Number of years to simulate
    
    Returns:
        List of all final cells
    """
    # Convert initial cells to Cell objects
    current_population = []
    for orig_idx, cell_dict, jsd, decile in initial_cells:
        cell_obj = create_cell_from_dict(cell_dict)
        # Store metadata for tracking
        cell_obj.original_idx = orig_idx
        cell_obj.original_decile = decile
        cell_obj.generation = 0  # Original cells are generation 0
        current_population.append(cell_obj)
    
    print(f"\nStarting with {len(current_population)} cells")
    
    # Age for specified years
    for year in range(1, years_to_age + 1):
        new_population = []
        
        # Each cell divides then both daughters age
        for parent_cell in current_population:
            # Create two daughter cells (exact copies)
            daughter1 = copy.deepcopy(parent_cell)
            daughter2 = copy.deepcopy(parent_cell)
            
            # Update generation
            daughter1.generation = parent_cell.generation + 1
            daughter2.generation = parent_cell.generation + 1
            
            # Both daughters age independently
            daughter1.age_1_year()
            daughter2.age_1_year()
            
            # Add to new population
            new_population.extend([daughter1, daughter2])
        
        current_population = new_population
        print(f"  Year {year}: {len(current_population)} cells")
    
    return current_population


def save_aged_cells(cells, original_metadata, years_aged, output_filename):
    """Save aged cells to json.gz file."""
    # Convert cells back to dictionary format
    cells_data = []
    for i, cell in enumerate(cells):
        cell_dict = cell.to_dict()
        # Add tracking metadata
        cell_dict['original_idx'] = cell.original_idx
        cell_dict['original_decile'] = cell.original_decile
        cell_dict['generation'] = cell.generation
        cells_data.append(cell_dict)
    
    # Create output data structure
    output_data = {
        "metadata": {
            "source_file": original_metadata.get('source_file', 'unknown'),
            "original_year": original_metadata.get('extracted_year', 50),
            "final_year": original_metadata.get('extracted_year', 50) + years_aged,
            "years_aged_with_division": years_aged,
            "initial_cells": 30,
            "final_cells": len(cells_data),
            "cells_per_decile": 3,
            "division_model": "divide_then_methylate"
        },
        "cells": cells_data
    }
    
    # Save to file
    print(f"\nSaving {len(cells_data)} cells to {output_filename}...")
    with gzip.open(output_filename, 'wt') as f:
        json.dump(output_data, f, indent=2)
    
    # Report file size
    file_size_mb = os.path.getsize(output_filename) / (1024 * 1024)
    print(f"Saved successfully! File size: {file_size_mb:.2f} MB")


def main():
    # Set random seed for reproducibility
    random.seed(42)
    np.random.seed(42)
    
    # Parameters
    years_to_age = 10
    cells_per_decile = 3
    
    print("=" * 60)
    print("CELL AGING WITH DIVISION SIMULATION")
    print("=" * 60)
    print(f"Parameters:")
    print(f"  Starting year: 50")
    print(f"  Years to age: {years_to_age}")
    print(f"  Final year: 60")
    print(f"  Cells per decile: {cells_per_decile}")
    print(f"  Total starting cells: {10 * cells_per_decile}")
    print(f"  Expected final cells: {30 * (2 ** years_to_age)} ({30 * (2 ** years_to_age):,})")
    print("=" * 60)
    
    # Load snapshot
    snapshot_data = load_snapshot('year50_snapshot.json.gz')
    
    # Sample cells by deciles
    print("\nSampling cells by JSD deciles...")
    sampled_cells = sample_cells_by_deciles(snapshot_data, cells_per_decile)
    print(f"\nSampled {len(sampled_cells)} cells total")
    
    # Age cells with division
    print("\nAging cells with division...")
    start_time = time.time()
    final_cells = age_cells_with_division(sampled_cells, years_to_age)
    aging_time = time.time() - start_time
    
    print(f"\nAging complete in {aging_time:.1f} seconds")
    print(f"Final population: {len(final_cells)} cells")
    
    # Calculate final statistics
    final_methylation = [cell.methylation_proportion for cell in final_cells]
    final_jsd = [cell.JSD for cell in final_cells]
    
    print(f"\nFinal statistics:")
    print(f"  Average methylation: {np.mean(final_methylation):.2%}")
    print(f"  Average JSD: {np.mean(final_jsd):.4f}")
    
    # Save results
    output_filename = 'year60_from_decile_sampling_30cells.json.gz'
    save_aged_cells(final_cells, snapshot_data['metadata'], years_to_age, output_filename)
    
    print("\nSimulation complete!")


if __name__ == "__main__":
    main()