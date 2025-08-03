#!/usr/bin/env python3
"""
Create 30 control individuals using only original year 60 cells.
Each control individual has 5,120 cells (same size as mixed individuals).
Sampling is done without replacement within each individual.
"""

import json
import gzip
import numpy as np
import random
import os
import time


def load_snapshot(filepath):
    """Load snapshot data from json.gz file."""
    print(f"Loading {filepath}...")
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    return data


def create_control_individual(original_cells, control_id, n_cells=5120, seed=None):
    """
    Create one control individual by sampling from original cells.
    
    Args:
        original_cells: List of all original year 60 cells
        control_id: ID for this control individual
        n_cells: Number of cells to sample (default 5120 to match mixed individuals)
        seed: Random seed for reproducibility
    
    Returns:
        Dictionary with control individual data
    """
    if seed is not None:
        random.seed(seed + control_id)  # Different seed per individual
        np.random.seed(seed + control_id)
    
    # Sample cells WITHOUT replacement (uniformly from all 10,000 cells)
    if n_cells > len(original_cells):
        raise ValueError(f"Cannot sample {n_cells} cells from {len(original_cells)} available")
    
    sampled_cells = random.sample(original_cells, k=n_cells)
    
    # Create copies and tag as control
    control_cells = []
    for cell in sampled_cells:
        cell_copy = cell.copy()
        cell_copy['origin'] = 'control'
        control_cells.append(cell_copy)
    
    # Calculate summary statistics
    jsds = [cell['jsd'] for cell in control_cells]
    meths = [cell['methylation_proportion'] for cell in control_cells]
    
    mean_jsd = np.mean(jsds)
    mean_meth = np.mean(meths)
    
    # Create control individual data
    control_data = {
        "metadata": {
            "control_id": control_id,
            "type": "control",
            "n_cells": len(control_cells),
            "mean_jsd": mean_jsd,
            "std_jsd": np.std(jsds),
            "mean_methylation": mean_meth,
            "std_methylation": np.std(meths),
            "random_seed": seed + control_id if seed else None,
            "sampling": "without_replacement"
        },
        "cells": control_cells
    }
    
    return control_data


def save_control_individual(control_data, output_dir="control_individuals"):
    """Save control individual data to compressed JSON file."""
    os.makedirs(output_dir, exist_ok=True)
    
    control_id = control_data['metadata']['control_id']
    filename = f"control_individual_{control_id:02d}.json.gz"
    filepath = os.path.join(output_dir, filename)
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(control_data, f, indent=2)
    
    return filepath


def main():
    """Create 30 control individuals."""
    # Paths
    original_snapshot = "year60_original_snapshot.json.gz"
    output_dir = "control_individuals"
    
    # Parameters
    n_control_individuals = 30
    n_cells_per_individual = 5120  # Same as mixed individuals
    seed = 123  # Different from mixed individuals seed
    
    print("="*60)
    print("Creating 30 control individuals")
    print(f"Cells per individual: {n_cells_per_individual}")
    print(f"Sampling: WITHOUT replacement within each individual")
    print("="*60)
    
    # Load original cells
    print(f"\nLoading original year 60 cells...")
    original_data = load_snapshot(original_snapshot)
    original_cells = original_data['cells']
    print(f"Loaded {len(original_cells)} original cells")
    
    # Verify we have enough cells
    if n_cells_per_individual > len(original_cells):
        print(f"ERROR: Cannot sample {n_cells_per_individual} cells from {len(original_cells)}")
        return
    
    print(f"Will sample {n_cells_per_individual} cells ({n_cells_per_individual/len(original_cells)*100:.1f}%) for each individual")
    
    # Create control individuals
    start_time = time.time()
    created_files = []
    all_mean_jsds = []
    
    for i in range(n_control_individuals):
        print(f"\nCreating control individual {i+1}/{n_control_individuals}")
        
        # Create control individual
        control_data = create_control_individual(
            original_cells,
            control_id=i,
            n_cells=n_cells_per_individual,
            seed=seed
        )
        
        # Save control individual
        filepath = save_control_individual(control_data, output_dir)
        created_files.append(filepath)
        all_mean_jsds.append(control_data['metadata']['mean_jsd'])
        
        print(f"  Saved to: {filepath}")
        print(f"  Mean JSD: {control_data['metadata']['mean_jsd']:.4f}")
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"COMPLETED!")
    print(f"{'='*60}")
    print(f"Created {len(created_files)} control individuals")
    print(f"Total time: {total_time:.1f}s")
    print(f"Output directory: {output_dir}/")
    
    # Display summary statistics
    print(f"\nControl distribution summary ({len(all_mean_jsds)} individuals):")
    print(f"  Mean of means: {np.mean(all_mean_jsds):.4f}")
    print(f"  Std of means:  {np.std(all_mean_jsds):.4f}")
    print(f"  Min mean JSD:  {np.min(all_mean_jsds):.4f}")
    print(f"  Max mean JSD:  {np.max(all_mean_jsds):.4f}")
    print(f"  Range:         {np.max(all_mean_jsds) - np.min(all_mean_jsds):.4f}")
    
    # Show all values
    print(f"\nAll control mean JSD values:")
    for i, jsd in enumerate(all_mean_jsds):
        print(f"  Control {i:02d}: {jsd:.4f}")


if __name__ == "__main__":
    main()