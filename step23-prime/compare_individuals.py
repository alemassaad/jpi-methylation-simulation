#!/usr/bin/env python3
"""
Compare individuals between step23 and step23-prime.
These might differ slightly due to random sampling, but statistics should be similar.
"""

import json
import gzip
import sys
import os
import numpy as np
import glob

def load_individual(filepath):
    """Load an individual file."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    # Handle both formats
    if 'cells' in data:
        cells = data['cells']
        metadata = data.get('metadata', {})
    else:
        # Older format might just be a list
        cells = data if isinstance(data, list) else data.get('cells', [])
        metadata = {}
    
    return cells, metadata

def compare_individual_files(step23_file, prime_file):
    """Compare two individual files."""
    cells1, meta1 = load_individual(step23_file)
    cells2, meta2 = load_individual(prime_file)
    
    # Basic checks
    n_cells1 = len(cells1)
    n_cells2 = len(cells2)
    
    # Get JSD statistics
    jsds1 = [c['jsd'] for c in cells1]
    jsds2 = [c['jsd'] for c in cells2]
    
    mean1 = np.mean(jsds1)
    mean2 = np.mean(jsds2)
    
    return {
        'n_cells': (n_cells1, n_cells2),
        'mean_jsd': (mean1, mean2),
        'cells_match': n_cells1 == n_cells2,
        'jsd_diff': abs(mean1 - mean2)
    }

def compare_group(group_name):
    """Compare all individuals in a group."""
    print(f"\n--- {group_name.upper()} INDIVIDUALS ---")
    
    step23_dir = f'../step23/data/rate_0.005000/individuals/{group_name}'
    prime_dir = f'data/rate_0.005000/individuals/{group_name}'
    
    # Get file lists
    step23_files = sorted(glob.glob(f'{step23_dir}/individual_*.json.gz'))
    prime_files = sorted(glob.glob(f'{prime_dir}/individual_*.json.gz'))
    
    print(f"  Step23 files: {len(step23_files)}")
    print(f"  Step23-prime files: {len(prime_files)}")
    
    if len(step23_files) != len(prime_files):
        print(f"  ❌ Different number of individuals!")
        return None
    
    # Compare each individual
    all_results = []
    for i, (s23_file, prime_file) in enumerate(zip(step23_files, prime_files)):
        result = compare_individual_files(s23_file, prime_file)
        all_results.append(result)
        
        # Report major differences
        if not result['cells_match']:
            fname = os.path.basename(s23_file)
            print(f"  ⚠️  {fname}: cell count mismatch {result['n_cells'][0]} vs {result['n_cells'][1]}")
    
    # Calculate group statistics
    all_means_s23 = [r['mean_jsd'][0] for r in all_results]
    all_means_prime = [r['mean_jsd'][1] for r in all_results]
    
    group_mean_s23 = np.mean(all_means_s23)
    group_mean_prime = np.mean(all_means_prime)
    group_std_s23 = np.std(all_means_s23)
    group_std_prime = np.std(all_means_prime)
    
    print(f"\n  Group Statistics:")
    print(f"  Step23:       {group_mean_s23:.6f} ± {group_std_s23:.6f}")
    print(f"  Step23-prime: {group_mean_prime:.6f} ± {group_std_prime:.6f}")
    print(f"  Difference:   {abs(group_mean_s23 - group_mean_prime):.6f}")
    
    # Check if all have correct cell count
    expected_cells = 5120  # After mixing
    cells_correct = all(r['n_cells'][0] == expected_cells and r['n_cells'][1] == expected_cells 
                        for r in all_results)
    
    if cells_correct:
        print(f"  ✓ All individuals have {expected_cells} cells")
    else:
        print(f"  ⚠️  Some individuals have wrong cell count (expected {expected_cells})")
    
    return {
        'group': group_name,
        'n_individuals': len(all_results),
        'mean_jsd': (group_mean_s23, group_mean_prime),
        'std_jsd': (group_std_s23, group_std_prime),
        'all_cells_correct': cells_correct
    }

def main():
    print("="*60)
    print("INDIVIDUAL COMPARISON: step23 vs step23-prime")
    print("="*60)
    print("\nNote: Individuals may differ due to random sampling,")
    print("but group statistics should be very similar.")
    
    results = {}
    
    # Compare each group
    for group in ['mutant', 'control1', 'control2']:
        result = compare_group(group)
        if result:
            results[group] = result
    
    # Final summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    
    if len(results) == 3:
        # Check if means are close
        all_close = True
        for group, data in results.items():
            diff = abs(data['mean_jsd'][0] - data['mean_jsd'][1])
            print(f"{group.capitalize():10s}: Δ = {diff:.6f}", end="")
            if diff < 0.01:  # Within 0.01 JSD
                print(" ✓")
            else:
                print(" ⚠️  (>0.01 difference)")
                all_close = False
        
        if all_close:
            print("\n✅ All group means are within 0.01 JSD - statistically equivalent!")
        else:
            print("\n⚠️  Some groups show larger differences.")
            print("This could be due to:")
            print("  - Different random seeds in sampling")
            print("  - Different random number generator state")
            print("  - Actual implementation differences")
    else:
        print("❌ Could not compare all groups")

if __name__ == "__main__":
    main()