#!/usr/bin/env python3
"""
Compare two runs of step23-prime to check if they're identical.
"""

import sys
import os
import json
import gzip
import hashlib
import glob
import numpy as np

def get_file_hash(filepath):
    """Get MD5 hash of a file."""
    if not os.path.exists(filepath):
        return None
    with open(filepath, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()[:8]

def compare_directories(dir1, dir2):
    """Compare two pipeline output directories."""
    
    print(f"Comparing:")
    print(f"  Run 1: {dir1}")
    print(f"  Run 2: {dir2}")
    print("="*50)
    
    results = {
        'snapshots_identical': True,
        'individuals_identical': True,
        'stats_identical': True
    }
    
    # 1. Compare snapshots (should be IDENTICAL)
    print("\n1. SNAPSHOTS:")
    for year in [50, 60]:
        file1 = f"{dir1}/snapshots/year{year}_snapshot.json.gz"
        file2 = f"{dir2}/snapshots/year{year}_snapshot.json.gz"
        
        hash1 = get_file_hash(file1)
        hash2 = get_file_hash(file2)
        
        if hash1 and hash2:
            if hash1 == hash2:
                print(f"  Year {year}: ✅ Identical (hash: {hash1})")
            else:
                print(f"  Year {year}: ❌ Different! ({hash1} vs {hash2})")
                results['snapshots_identical'] = False
    
    # 2. Compare individuals
    print("\n2. INDIVIDUALS:")
    for group in ['mutant', 'control1', 'control2']:
        print(f"\n  {group.upper()}:")
        
        # Get file lists
        files1 = sorted(glob.glob(f"{dir1}/individuals/{group}/individual_*.json.gz"))
        files2 = sorted(glob.glob(f"{dir2}/individuals/{group}/individual_*.json.gz"))
        
        if len(files1) != len(files2):
            print(f"    ❌ Different file counts: {len(files1)} vs {len(files2)}")
            results['individuals_identical'] = False
            continue
        
        # Sample first few individuals for detailed comparison
        identical_count = 0
        different_count = 0
        
        for i in range(min(3, len(files1))):  # Check first 3
            with gzip.open(files1[i], 'rt') as f:
                data1 = json.load(f)
            with gzip.open(files2[i], 'rt') as f:
                data2 = json.load(f)
            
            cells1 = data1['cells'] if 'cells' in data1 else data1
            cells2 = data2['cells'] if 'cells' in data2 else data2
            
            # Compare JSD distributions
            jsds1 = sorted([c['jsd'] for c in cells1])
            jsds2 = sorted([c['jsd'] for c in cells2])
            
            if jsds1 == jsds2:
                identical_count += 1
            else:
                different_count += 1
                mean_diff = abs(np.mean(jsds1) - np.mean(jsds2))
                print(f"    Individual {i:02d}: Mean diff = {mean_diff:.9f}")
        
        if identical_count == min(3, len(files1)):
            print(f"    ✅ First {identical_count} individuals identical")
        else:
            print(f"    ⚠️  {different_count} individuals differ")
            results['individuals_identical'] = False
    
    # 3. Compare statistics
    print("\n3. STATISTICS:")
    stats1_file = f"{dir1}/results/statistics.json"
    stats2_file = f"{dir2}/results/statistics.json"
    
    if os.path.exists(stats1_file) and os.path.exists(stats2_file):
        with open(stats1_file) as f:
            stats1 = json.load(f)
        with open(stats2_file) as f:
            stats2 = json.load(f)
        
        for group in ['mutant', 'control1', 'control2']:
            if group in stats1 and group in stats2:
                mean1 = stats1[group]['mean']
                mean2 = stats2[group]['mean']
                diff = abs(mean1 - mean2)
                
                if diff < 1e-9:
                    print(f"  {group}: ✅ Identical (mean = {mean1:.6f})")
                else:
                    print(f"  {group}: ⚠️  Difference = {diff:.9f}")
                    if diff > 0.001:
                        results['stats_identical'] = False
    
    # Final verdict
    print("\n" + "="*50)
    print("VERDICT:")
    
    if results['snapshots_identical'] and results['individuals_identical']:
        print("✅ RUNS ARE IDENTICAL - Perfect reproducibility!")
    elif results['snapshots_identical'] and not results['individuals_identical']:
        print("⚠️  Snapshots match but individuals differ")
        print("    This suggests randomness in growth/mixing")
    else:
        print("❌ RUNS DIFFER - Reproducibility issue detected")
    
    return results

def main():
    if len(sys.argv) != 3:
        print("Usage: python compare_two_runs.py <dir1> <dir2>")
        sys.exit(1)
    
    dir1 = sys.argv[1]
    dir2 = sys.argv[2]
    
    if not os.path.exists(dir1):
        print(f"Error: {dir1} not found")
        sys.exit(1)
    if not os.path.exists(dir2):
        print(f"Error: {dir2} not found")
        sys.exit(1)
    
    results = compare_directories(dir1, dir2)
    
    # Exit code
    if results['snapshots_identical'] and results['individuals_identical']:
        sys.exit(0)  # Success
    else:
        sys.exit(1)  # Differences found

if __name__ == "__main__":
    main()