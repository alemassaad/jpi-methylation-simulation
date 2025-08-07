#!/usr/bin/env python3
"""
Test version of create_control_individuals.py with smaller parameters.
Creates 5 control individuals with 1000 cells each.
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


def create_control_individual(original_cells, control_id, n_cells=1000, seed=None):
    """
    Create one control individual by sampling from original cells.
    TEST VERSION: Smaller sample size
    """
    if seed is not None:
        random.seed(seed + control_id)
        np.random.seed(seed + control_id)
    
    # Sample cells WITHOUT replacement
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
    
    # Verify no duplicates (should be true with random.sample)
    cell_ids = [id(cell) for cell in sampled_cells]
    n_unique = len(set(cell_ids))
    if n_unique != n_cells:
        print(f"  WARNING: Found duplicates! {n_cells} cells but only {n_unique} unique")
    
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


def save_control_individual(control_data, output_dir="test_control_individuals"):
    """Save control individual data to compressed JSON file."""
    os.makedirs(output_dir, exist_ok=True)
    
    control_id = control_data['metadata']['control_id']
    filename = f"control_individual_{control_id:02d}.json.gz"
    filepath = os.path.join(output_dir, filename)
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(control_data, f, indent=2)
    
    return filepath


def main():
    """TEST: Create 5 control individuals with 1000 cells each."""
    # Paths
    original_snapshot = "year60_original_snapshot.json.gz"
    output_dir = "test_control_individuals"
    
    # TEST Parameters
    n_control_individuals = 5    # Fewer individuals
    n_cells_per_individual = 1000  # Fewer cells
    seed = 123
    
    print("="*60)
    print("TEST MODE: Creating control individuals")
    print(f"Individuals: {n_control_individuals}")
    print(f"Cells per individual: {n_cells_per_individual}")
    print("="*60)
    
    # Load original cells
    print(f"\nLoading original year 60 cells...")
    original_data = load_snapshot(original_snapshot)
    original_cells = original_data['cells']
    print(f"Loaded {len(original_cells)} original cells")
    
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
        print(f"  Cells sampled: {control_data['metadata']['n_cells']}")
    
    # Verify files
    print(f"\nVerifying saved files...")
    for filepath in created_files:
        with gzip.open(filepath, 'rt') as f:
            data = json.load(f)
        n_cells = len(data['cells'])
        mean_jsd = data['metadata']['mean_jsd']
        print(f"  {os.path.basename(filepath)}: {n_cells} cells, mean JSD = {mean_jsd:.4f}")
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"TEST COMPLETED!")
    print(f"{'='*60}")
    print(f"Created {len(created_files)} control individuals")
    print(f"Total time: {total_time:.1f}s")
    print(f"\nControl mean JSD values:")
    for i, jsd in enumerate(all_mean_jsds):
        print(f"  Control {i}: {jsd:.4f}")
    print(f"\nMean: {np.mean(all_mean_jsds):.4f} ± {np.std(all_mean_jsds):.4f}")


if __name__ == "__main__":
    main()