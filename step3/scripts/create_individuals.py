#!/usr/bin/env python3
"""
Create individuals by mixing lineage cells with original year 60 cells.
Works for both mutant lineages (decile-based) and control lineages (uniform).
Each individual: 80% original cells + 20% lineage cells.
"""

import argparse
import json
import gzip
import numpy as np
import random
import os
import glob
import time
import sys
sys.path.append('../..')  # To import from parent directory

from step1.cell import Cell


def load_snapshot(filepath):
    """Load snapshot data from json.gz file."""
    print(f"Loading {filepath}...")
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    return data


def load_lineage(filepath):
    """Load a single lineage file."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    return data


def create_individual_from_lineage(lineage_data, original_cells, individual_id, 
                                 individual_type="mutant", ratio=0.8, seed=None):
    """
    Create one individual by mixing lineage cells with original cells.
    
    Args:
        lineage_data: Data from one lineage file
        original_cells: List of all original year 60 cells
        individual_id: ID for this individual
        individual_type: "mutant" or "control"
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
    
    # Sample original cells WITHOUT replacement
    if n_original > len(original_cells):
        raise ValueError(f"Need {n_original} original cells but only {len(original_cells)} available")
    
    sampled_original_indices = random.sample(range(len(original_cells)), n_original)
    sampled_original_cells = [original_cells[i] for i in sampled_original_indices]
    
    # Combine cells
    all_cells = []
    
    # Add lineage cells with source tracking
    for cell in lineage_cells:
        cell_copy = cell.copy()
        cell_copy['source'] = 'lineage'
        cell_copy['lineage_id'] = lineage_data['metadata']['lineage_id']
        all_cells.append(cell_copy)
    
    # Add original cells with source tracking
    for cell in sampled_original_cells:
        cell_copy = cell.copy()
        cell_copy['source'] = 'original'
        all_cells.append(cell_copy)
    
    # Shuffle to mix them
    random.shuffle(all_cells)
    
    # Calculate statistics
    jsd_values = [cell['jsd'] for cell in all_cells]
    mean_jsd = np.mean(jsd_values)
    std_jsd = np.std(jsd_values)
    
    # Get lineage metadata
    lineage_metadata = lineage_data['metadata']
    
    # Create individual data structure
    individual_data = {
        "metadata": {
            "individual_id": individual_id,
            "individual_type": individual_type,
            "lineage_id": lineage_metadata['lineage_id'],
            "source_decile": lineage_metadata.get('source_decile', 0),
            "n_lineage_cells": n_lineage,
            "n_original_cells": n_original,
            "total_cells": len(all_cells),
            "ratio_original": ratio,
            "mean_jsd": mean_jsd,
            "std_jsd": std_jsd,
            "creation_time": time.time()
        },
        "cells": all_cells
    }
    
    return individual_data


def process_lineages(lineage_dir, original_cells, output_dir, individual_type="mutant", 
                    ratio=0.8, seed=42):
    """Process all lineages in a directory to create individuals."""
    
    # Find all lineage files
    pattern = os.path.join(lineage_dir, "*.json.gz")
    lineage_files = sorted(glob.glob(pattern))
    
    print(f"\nFound {len(lineage_files)} {individual_type} lineage files")
    
    if not lineage_files:
        print(f"Warning: No lineage files found in {lineage_dir}")
        return
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Process each lineage
    start_time = time.time()
    
    for idx, lineage_file in enumerate(lineage_files):
        print(f"\nProcessing {individual_type} individual {idx + 1}/{len(lineage_files)}")
        print(f"  From: {os.path.basename(lineage_file)}")
        
        # Load lineage
        lineage_data = load_lineage(lineage_file)
        
        # Create individual
        individual_data = create_individual_from_lineage(
            lineage_data, original_cells, idx, individual_type, ratio, seed
        )
        
        # Save individual
        filename = f"individual_{idx:02d}.json.gz"
        filepath = os.path.join(output_dir, filename)
        
        with gzip.open(filepath, 'wt') as f:
            json.dump(individual_data, f, indent=2)
        
        print(f"  Saved to: {filepath}")
    
    # Summary
    total_time = time.time() - start_time
    print(f"\n{individual_type.upper()} SUMMARY:")
    print(f"  Individuals created: {len(lineage_files)}")
    print(f"  Time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"  Output directory: {output_dir}/")


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Create individuals by mixing lineage cells with original cells",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-t', '--type', choices=['mutant', 'control', 'both'], 
                        default='both', help='Type of individuals to create')
    parser.add_argument('-r', '--ratio', type=float, default=0.8,
                        help='Fraction of original cells (0.8 = 80%% original, 20%% lineage)')
    parser.add_argument('--year60-snapshot', 
                        default='../data/snapshots/year60_original_snapshot.json.gz',
                        help='Path to year 60 snapshot file')
    parser.add_argument('--lineage-base-dir', 
                        default='../../step2/data/lineages',
                        help='Base directory containing mutant/ and control/ subdirectories')
    parser.add_argument('--test-mode', action='store_true',
                        help='Use test lineage directories (test_lineages and test_lineages_control)')
    parser.add_argument('--output-base-dir',
                        default='../data/individuals',
                        help='Base output directory')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed for reproducibility')
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()
    
    print("=" * 60)
    print("INDIVIDUAL CREATION FROM LINEAGES")
    print("=" * 60)
    print(f"Parameters:")
    print(f"  Type: {args.type}")
    print(f"  Ratio: {args.ratio:.1%} original, {1-args.ratio:.1%} lineage")
    print(f"  Seed: {args.seed}")
    print("=" * 60)
    
    # Load year 60 snapshot
    snapshot_data = load_snapshot(args.year60_snapshot)
    original_cells = snapshot_data['cells']
    print(f"\nLoaded {len(original_cells)} original year 60 cells")
    
    # Process mutant individuals if requested
    if args.type in ['mutant', 'both']:
        print("\n### CREATING MUTANT INDIVIDUALS ###")
        if args.test_mode:
            mutant_lineage_dir = os.path.join(args.lineage_base_dir, 'test_lineages')
        else:
            mutant_lineage_dir = os.path.join(args.lineage_base_dir, 'mutant')
        mutant_output_dir = os.path.join(args.output_base_dir, 'mutant')
        
        process_lineages(mutant_lineage_dir, original_cells, mutant_output_dir,
                        'mutant', args.ratio, args.seed)
    
    # Process control individuals if requested
    if args.type in ['control', 'both']:
        print("\n### CREATING CONTROL INDIVIDUALS ###")
        if args.test_mode:
            control_lineage_dir = os.path.join(args.lineage_base_dir, 'test_lineages_control')
        else:
            control_lineage_dir = os.path.join(args.lineage_base_dir, 'control')
        control_output_dir = os.path.join(args.output_base_dir, 'control')
        
        process_lineages(control_lineage_dir, original_cells, control_output_dir,
                        'control', args.ratio, args.seed)
    
    print("\n" + "=" * 60)
    print("COMPLETED!")
    print("=" * 60)


if __name__ == "__main__":
    main()