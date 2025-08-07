#!/usr/bin/env python3
"""
Create 30 individuals by mixing each lineage with original year 60 cells.
Each individual represents a mixed population with 80% original cells and 20% lineage cells.
"""

import json
import gzip
import numpy as np
import random
import os
import glob
import time


def load_snapshot(filepath):
    """Load snapshot data from json.gz file."""
    print(f"Loading {filepath}...")
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    return data


def create_individual(lineage_data, original_cells, individual_id, ratio=0.8, seed=None):
    """
    Create one individual by mixing lineage cells with original cells.
    
    Args:
        lineage_data: Data from one lineage file
        original_cells: List of all original year 60 cells
        individual_id: ID for this individual
        ratio: Fraction of original cells (default 0.8 = 80%)
        seed: Random seed for reproducibility
    
    Returns:
        Dictionary with individual data
    """
    if seed is not None:
        random.seed(seed + individual_id)  # Different seed per individual
        np.random.seed(seed + individual_id)
    
    # Get lineage cells
    lineage_cells = lineage_data['cells']
    n_lineage = len(lineage_cells)
    
    # Calculate how many original cells we need
    # If lineage is 20%, then n_lineage = 0.2 * total
    # So total = n_lineage / (1 - ratio)
    total_cells = int(n_lineage / (1 - ratio))
    n_original = total_cells - n_lineage
    
    print(f"  Lineage cells: {n_lineage}")
    print(f"  Original cells to sample: {n_original}")
    print(f"  Total cells: {total_cells}")
    
    # Sample original cells (with replacement since we need more than available)
    sampled_original = random.choices(original_cells, k=n_original)
    
    # Create copies and tag cells with their origin
    tagged_lineage = []
    for cell in lineage_cells:
        cell_copy = cell.copy()
        cell_copy['origin'] = 'lineage'
        tagged_lineage.append(cell_copy)
    
    tagged_original = []
    for cell in sampled_original:
        cell_copy = cell.copy()
        cell_copy['origin'] = 'original'
        tagged_original.append(cell_copy)
    
    # Combine all cells
    all_cells = tagged_lineage + tagged_original
    
    # Calculate summary statistics
    jsds = [cell['jsd'] for cell in all_cells]
    meths = [cell['methylation_proportion'] for cell in all_cells]
    
    mean_jsd = np.mean(jsds)
    mean_meth = np.mean(meths)
    
    # Create individual data
    individual_data = {
        "metadata": {
            "individual_id": individual_id,
            "source_lineage_id": lineage_data['metadata']['lineage_id'],
            "source_decile": lineage_data['metadata']['source_decile'],
            "source_lineage_file": lineage_data['metadata'].get('source_file', 'unknown'),
            "mixing_ratio": ratio,
            "n_original": n_original,
            "n_lineage": n_lineage,
            "total_cells": len(all_cells),
            "mean_jsd": mean_jsd,
            "std_jsd": np.std(jsds),
            "mean_methylation": mean_meth,
            "std_methylation": np.std(meths),
            "random_seed": seed + individual_id if seed else None
        },
        "cells": all_cells
    }
    
    return individual_data


def save_individual(individual_data, output_dir="individuals"):
    """Save individual data to compressed JSON file."""
    os.makedirs(output_dir, exist_ok=True)
    
    ind_id = individual_data['metadata']['individual_id']
    lineage_id = individual_data['metadata']['source_lineage_id']
    decile = individual_data['metadata']['source_decile']
    
    filename = f"individual_{ind_id:02d}_from_lineage_{lineage_id:02d}_decile_{decile:02d}.json.gz"
    filepath = os.path.join(output_dir, filename)
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(individual_data, f, indent=2)
    
    return filepath


def main():
    """Create 30 individuals from lineage files."""
    # Paths
    original_snapshot = "year60_original_snapshot.json.gz"
    lineages_dir = "../step2/lineages"
    output_dir = "individuals"
    
    # Parameters
    ratio = 0.8  # 80% original, 20% lineage
    seed = 42    # For reproducibility
    
    print("="*60)
    print("Creating 30 individuals from lineages")
    print(f"Mixing ratio: {ratio*100:.0f}% original, {(1-ratio)*100:.0f}% lineage")
    print("="*60)
    
    # Load original cells
    print(f"\nLoading original year 60 cells...")
    original_data = load_snapshot(original_snapshot)
    original_cells = original_data['cells']
    print(f"Loaded {len(original_cells)} original cells")
    
    # Get all lineage files
    lineage_files = sorted(glob.glob(os.path.join(lineages_dir, "lineage_*.json.gz")))
    print(f"\nFound {len(lineage_files)} lineage files")
    
    if len(lineage_files) != 30:
        print(f"WARNING: Expected 30 lineage files, found {len(lineage_files)}")
    
    # Process each lineage
    start_time = time.time()
    created_files = []
    all_mean_jsds = []
    
    for i, lineage_file in enumerate(lineage_files):
        print(f"\n{'='*40}")
        print(f"Processing lineage {i+1}/{len(lineage_files)}: {os.path.basename(lineage_file)}")
        
        # Load lineage
        lineage_data = load_snapshot(lineage_file)
        
        # Store source file info
        lineage_data['metadata']['source_file'] = os.path.basename(lineage_file)
        
        # Create individual
        individual_data = create_individual(
            lineage_data, 
            original_cells, 
            individual_id=i,
            ratio=ratio,
            seed=seed
        )
        
        # Save individual
        filepath = save_individual(individual_data, output_dir)
        created_files.append(filepath)
        all_mean_jsds.append(individual_data['metadata']['mean_jsd'])
        
        print(f"  Saved to: {filepath}")
        print(f"  Mean JSD: {individual_data['metadata']['mean_jsd']:.4f}")
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n{'='*60}")
    print(f"COMPLETED!")
    print(f"{'='*60}")
    print(f"Created {len(created_files)} individuals")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"Output directory: {output_dir}/")
    
    # Display summary statistics
    print(f"\nSummary of mean JSD values across {len(all_mean_jsds)} individuals:")
    print(f"  Overall mean: {np.mean(all_mean_jsds):.4f}")
    print(f"  Overall std:  {np.std(all_mean_jsds):.4f}")
    print(f"  Min JSD:      {np.min(all_mean_jsds):.4f}")
    print(f"  Max JSD:      {np.max(all_mean_jsds):.4f}")
    
    # Group by decile
    print(f"\nMean JSD by source decile:")
    decile_jsds = {}
    for filepath in created_files:
        with gzip.open(filepath, 'rt') as f:
            data = json.load(f)
        meta = data['metadata']
        decile = meta['source_decile']
        if decile not in decile_jsds:
            decile_jsds[decile] = []
        decile_jsds[decile].append(meta['mean_jsd'])
    
    for decile in sorted(decile_jsds.keys()):
        jsds = decile_jsds[decile]
        print(f"  Decile {decile:02d}: {np.mean(jsds):.4f} ± {np.std(jsds):.4f} (n={len(jsds)})")


if __name__ == "__main__":
    main()