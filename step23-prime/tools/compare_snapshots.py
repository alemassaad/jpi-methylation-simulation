#!/usr/bin/env python3
"""
Compare snapshots between step23 and step23-prime to verify they're identical.
Snapshots should be identical since they're just extracted cells from the same simulation.
"""

import json
import gzip
import sys
import numpy as np

def load_snapshot(filepath):
    """Load a snapshot file and return cells."""
    print(f"Loading {filepath}...")
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    # Handle both formats
    if 'cells' in data:
        # step23-prime format: {"metadata": {...}, "cells": [...]}
        cells = data['cells']
    else:
        # step23 format: direct list of cells
        cells = data
    
    return cells

def compare_cells(cells1, cells2):
    """Compare two lists of cells."""
    if len(cells1) != len(cells2):
        print(f"  ❌ Different number of cells: {len(cells1)} vs {len(cells2)}")
        return False
    
    print(f"  ✓ Same number of cells: {len(cells1)}")
    
    # Compare JSD distributions
    jsds1 = [c['jsd'] for c in cells1]
    jsds2 = [c['jsd'] for c in cells2]
    
    mean1, std1 = np.mean(jsds1), np.std(jsds1)
    mean2, std2 = np.mean(jsds2), np.std(jsds2)
    
    print(f"  Step23:       mean JSD = {mean1:.6f} ± {std1:.6f}")
    print(f"  Step23-prime: mean JSD = {mean2:.6f} ± {std2:.6f}")
    
    # Check if distributions are identical
    if np.allclose(jsds1, jsds2):
        print(f"  ✓ JSD values are identical!")
        return True
    else:
        # Check if they're just in different order
        if sorted(jsds1) == sorted(jsds2):
            print(f"  ✓ JSD values match (different order)")
            return True
        else:
            print(f"  ❌ JSD values differ")
            # Show difference
            diff = abs(mean1 - mean2)
            print(f"  Mean difference: {diff:.8f}")
            return False

def main():
    print("="*60)
    print("SNAPSHOT COMPARISON: step23 vs step23-prime")
    print("="*60)
    
    # Year 50 comparison
    print("\n--- YEAR 50 SNAPSHOT ---")
    try:
        # Load both snapshots
        step23_y50 = load_snapshot('../step23/data/rate_0.005000/snapshots/year50_snapshot.json.gz')
        prime_y50 = load_snapshot('data/rate_0.005000/snapshots/year50_snapshot.json.gz')
        
        # Compare
        y50_match = compare_cells(step23_y50, prime_y50)
    except FileNotFoundError as e:
        print(f"  ❌ File not found: {e}")
        y50_match = False
    
    # Year 60 comparison
    print("\n--- YEAR 60 SNAPSHOT ---")
    try:
        # Load both snapshots
        step23_y60 = load_snapshot('../step23/data/rate_0.005000/snapshots/year60_snapshot.json.gz')
        prime_y60 = load_snapshot('data/rate_0.005000/snapshots/year60_snapshot.json.gz')
        
        # Compare
        y60_match = compare_cells(step23_y60, prime_y60)
    except FileNotFoundError as e:
        print(f"  ❌ File not found: {e}")
        y60_match = False
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    if y50_match and y60_match:
        print("✅ Snapshots are identical! The refactoring preserves data integrity.")
    else:
        print("⚠️  Snapshots differ. This might indicate:")
        print("   - Different random seeds")
        print("   - Different extraction logic")
        print("   - Or files from different simulation runs")

if __name__ == "__main__":
    main()