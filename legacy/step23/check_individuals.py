#!/usr/bin/env python3
"""Check the state of individual files."""

import json
import gzip
import glob
import os
import sys

def check_individuals(base_dir):
    """Check all individual files in a directory."""
    files = sorted(glob.glob(os.path.join(base_dir, '*.json.gz')))
    
    if not files:
        print(f"No files found in {base_dir}")
        return
    
    cell_counts = {}
    corrupted = []
    
    for f in files:
        fname = os.path.basename(f)
        try:
            with gzip.open(f, 'rt') as file:
                data = json.load(file)
                cells = len(data.get('cells', []))
                cell_counts[fname] = cells
        except Exception as e:
            corrupted.append((fname, str(e)))
            cell_counts[fname] = 'CORRUPTED'
    
    # Summary
    valid_counts = [c for c in cell_counts.values() if c != 'CORRUPTED']
    if valid_counts:
        print(f"Files found: {len(files)}")
        print(f"Unique cell counts: {sorted(set(valid_counts))}")
        print(f"Files with 1024 cells: {valid_counts.count(1024)}")
        print(f"Files with 5120 cells: {valid_counts.count(5120)}")
    
    if corrupted:
        print(f"\nCorrupted files: {len(corrupted)}")
        for fname, error in corrupted[:3]:
            print(f"  {fname}: {error}")
    
    # Show sample
    print("\nSample files:")
    for i, fname in enumerate(sorted(cell_counts.keys())):
        if i < 3 or i >= len(cell_counts) - 3:
            print(f"  {fname}: {cell_counts[fname]} cells")
        elif i == 3:
            print("  ...")

if __name__ == "__main__":
    base = '/Users/alessandromassaad/jpi-methylation-simulation/step23/data/rate_0.005000/individuals'
    
    for group in ['mutant', 'control1', 'control2']:
        group_dir = os.path.join(base, group)
        if os.path.exists(group_dir):
            print(f"\n{'='*60}")
            print(f"Checking {group.upper()} individuals")
            print(f"{'='*60}")
            check_individuals(group_dir)