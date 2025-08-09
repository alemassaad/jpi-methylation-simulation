#!/usr/bin/env python3
"""
Detailed cell-by-cell comparison of snapshots between step23 and step23-prime.
This checks if the exact same cells were extracted in the exact same order.
"""

import json
import gzip
import sys
import numpy as np
from typing import List, Dict, Any

def load_snapshot(filepath: str) -> List[Dict[str, Any]]:
    """Load a snapshot file and return cells."""
    print(f"Loading {filepath}...")
    try:
        with gzip.open(filepath, 'rt') as f:
            data = json.load(f)
        
        # Handle both formats
        if isinstance(data, dict) and 'cells' in data:
            # step23-prime format: {"metadata": {...}, "cells": [...]}
            cells = data['cells']
            print(f"  Format: step23-prime (dict with 'cells' key)")
        elif isinstance(data, list):
            # step23 format: direct list of cells
            cells = data
            print(f"  Format: step23 (direct list)")
        else:
            # Might be a dict without 'cells' key
            cells = data
            print(f"  Format: unknown (assuming dict)")
        
        print(f"  Loaded {len(cells)} cells")
        return cells
    except Exception as e:
        print(f"  ERROR loading file: {e}")
        return []

def compare_single_cell(cell1: Dict, cell2: Dict, index: int) -> bool:
    """Compare two cells and report differences."""
    differences = []
    
    # Check all fields
    fields_to_check = ['cpg_sites', 'age', 'rate', 'gene_size', 'jsd', 
                       'methylation_proportion', 'methylation_distribution']
    
    for field in fields_to_check:
        if field not in cell1 and field not in cell2:
            continue  # Field missing in both, skip
        
        if field not in cell1:
            differences.append(f"{field} missing in step23")
            continue
            
        if field not in cell2:
            differences.append(f"{field} missing in step23-prime")
            continue
        
        val1 = cell1[field]
        val2 = cell2[field]
        
        if field == 'cpg_sites':
            # For CpG sites, check if arrays are identical
            if val1 != val2:
                # Count differences
                n_diff = sum(1 for a, b in zip(val1, val2) if a != b)
                differences.append(f"cpg_sites: {n_diff}/{len(val1)} sites differ")
        
        elif field in ['jsd', 'methylation_proportion', 'rate']:
            # For floats, check if close enough
            if not np.isclose(val1, val2, rtol=1e-9):
                differences.append(f"{field}: {val1:.9f} vs {val2:.9f}")
        
        elif field == 'methylation_distribution':
            # Check if distributions match
            if not all(np.isclose(a, b, rtol=1e-9) for a, b in zip(val1, val2)):
                differences.append(f"methylation_distribution differs")
        
        else:
            # Direct comparison for other fields
            if val1 != val2:
                differences.append(f"{field}: {val1} vs {val2}")
    
    if differences:
        print(f"\n  Cell {index:4d} DIFFERS:")
        for diff in differences:
            print(f"    - {diff}")
        return False
    
    return True

def compare_snapshots_detailed(cells1: List[Dict], cells2: List[Dict], 
                              max_cells_to_show: int = 10) -> Dict:
    """Perform detailed cell-by-cell comparison."""
    
    results = {
        'n_cells_match': len(cells1) == len(cells2),
        'n_cells1': len(cells1),
        'n_cells2': len(cells2),
        'cells_identical': False,
        'n_different_cells': 0,
        'different_indices': []
    }
    
    if not results['n_cells_match']:
        print(f"\n❌ Different number of cells: {len(cells1)} vs {len(cells2)}")
        return results
    
    print(f"\n✓ Same number of cells: {len(cells1)}")
    print(f"\nComparing {len(cells1)} cells one by one...")
    
    # Compare each cell
    identical_count = 0
    different_count = 0
    
    for i in range(len(cells1)):
        if compare_single_cell(cells1[i], cells2[i], i):
            identical_count += 1
        else:
            different_count += 1
            results['different_indices'].append(i)
            
            # Stop showing individual differences after max_cells_to_show
            if different_count >= max_cells_to_show:
                print(f"\n  ... stopping detailed output after {max_cells_to_show} differences")
                # Continue counting but don't print
                for j in range(i + 1, len(cells1)):
                    if not compare_single_cell(cells1[j], cells2[j], j):
                        different_count += 1
                        results['different_indices'].append(j)
                break
    
    results['n_different_cells'] = different_count
    results['cells_identical'] = (different_count == 0)
    
    # Summary statistics
    print(f"\n{'='*40}")
    print(f"CELL-BY-CELL COMPARISON RESULTS:")
    print(f"  Identical cells: {identical_count}/{len(cells1)}")
    print(f"  Different cells: {different_count}/{len(cells1)}")
    
    if different_count == 0:
        print(f"\n✅ ALL CELLS ARE IDENTICAL!")
    else:
        print(f"\n⚠️  {different_count} cells differ")
        
        # Analyze JSD distribution even if cells differ
        jsds1 = [c['jsd'] for c in cells1]
        jsds2 = [c['jsd'] for c in cells2]
        
        print(f"\nJSD Distribution Comparison:")
        print(f"  Step23:       mean={np.mean(jsds1):.6f}, std={np.std(jsds1):.6f}")
        print(f"  Step23-prime: mean={np.mean(jsds2):.6f}, std={np.std(jsds2):.6f}")
        print(f"  Difference:   {abs(np.mean(jsds1) - np.mean(jsds2)):.9f}")
        
        # Check if sorted JSDs match
        if sorted(jsds1) == sorted(jsds2):
            print(f"\n  Note: JSD values match when sorted (just different order)")
        
    return results

def main():
    print("="*60)
    print("DETAILED SNAPSHOT COMPARISON: step23 vs step23-prime")
    print("="*60)
    print("\nThis script does CELL-BY-CELL comparison to check if")
    print("the exact same cells were extracted in the exact same order.\n")
    
    # File paths
    step23_y50 = '../step23/data/rate_0.005000/snapshots/year50_snapshot.json.gz'
    prime_y50 = 'data/rate_0.005000/snapshots/year50_snapshot.json.gz'
    step23_y60 = '../step23/data/rate_0.005000/snapshots/year60_snapshot.json.gz'
    prime_y60 = 'data/rate_0.005000/snapshots/year60_snapshot.json.gz'
    
    # Year 50 comparison
    print("\n" + "="*60)
    print("YEAR 50 SNAPSHOT")
    print("="*60)
    
    cells1_y50 = load_snapshot(step23_y50)
    cells2_y50 = load_snapshot(prime_y50)
    
    if cells1_y50 and cells2_y50:
        results_y50 = compare_snapshots_detailed(cells1_y50, cells2_y50)
    else:
        print("❌ Could not load year 50 snapshots")
        results_y50 = None
    
    # Year 60 comparison
    print("\n" + "="*60)
    print("YEAR 60 SNAPSHOT")
    print("="*60)
    
    cells1_y60 = load_snapshot(step23_y60)
    cells2_y60 = load_snapshot(prime_y60)
    
    if cells1_y60 and cells2_y60:
        results_y60 = compare_snapshots_detailed(cells1_y60, cells2_y60)
    else:
        print("❌ Could not load year 60 snapshots")
        results_y60 = None
    
    # Final summary
    print("\n" + "="*60)
    print("FINAL SUMMARY")
    print("="*60)
    
    if results_y50 and results_y60:
        if results_y50['cells_identical'] and results_y60['cells_identical']:
            print("\n✅ PERFECT MATCH! Both snapshots are identical cell-by-cell!")
            print("   This means step23-prime extracts exactly the same cells")
            print("   in exactly the same order as the original step23.")
        elif results_y50['n_cells_match'] and results_y60['n_cells_match']:
            print("\n⚠️  Snapshots have same number of cells but some cells differ.")
            print("   Year 50: {} different cells".format(results_y50['n_different_cells']))
            print("   Year 60: {} different cells".format(results_y60['n_different_cells']))
            print("\n   Possible reasons:")
            print("   - Different extraction logic")
            print("   - Different file format handling")
            print("   - Floating point precision differences")
        else:
            print("\n❌ Snapshots have different numbers of cells!")
            print("   This indicates a fundamental difference in extraction logic.")
    else:
        print("\n❌ Could not complete comparison due to missing files.")
        print("   Make sure both pipelines have been run with rate 0.005")

if __name__ == "__main__":
    main()