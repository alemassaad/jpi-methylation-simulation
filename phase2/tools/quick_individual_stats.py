#!/usr/bin/env python3
"""
Quick statistics for all individuals - faster version that just checks basics.
"""

import json
import gzip
import numpy as np
import glob
import os

def quick_check_group(group_name, step23_dir, prime_dir):
    """Quick check of a group without loading all cells."""
    
    print(f"\n{group_name.upper()}:")
    
    # Count files
    step23_files = glob.glob(f'{step23_dir}/individuals/{group_name}/individual_*.json.gz')
    prime_files = glob.glob(f'{prime_dir}/individuals/{group_name}/individual_*.json.gz')
    
    print(f"  Files: step23={len(step23_files)}, phase2={len(prime_files)}")
    
    if len(step23_files) > 0 and len(prime_files) > 0:
        # Sample first and last files for quick check
        for idx in [0, len(step23_files)-1]:
            if idx < len(step23_files) and idx < len(prime_files):
                s23_file = sorted(step23_files)[idx]
                pr_file = sorted(prime_files)[idx]
                
                # Load just metadata if available, otherwise count cells
                with gzip.open(s23_file, 'rt') as f:
                    s23_data = json.load(f)
                with gzip.open(pr_file, 'rt') as f:
                    pr_data = json.load(f)
                
                s23_cells = s23_data['cells'] if 'cells' in s23_data else s23_data
                pr_cells = pr_data['cells'] if 'cells' in pr_data else pr_data
                
                print(f"  Individual {idx:02d}: step23={len(s23_cells)} cells, prime={len(pr_cells)} cells")

def main():
    step23_dir = "../step23/data/rate_0.005000"
    prime_dir = "data/rate_0.005000"
    
    print("QUICK INDIVIDUAL CHECK")
    print("=" * 40)
    
    for group in ['mutant', 'control1', 'control2']:
        quick_check_group(group, step23_dir, prime_dir)
    
    print("\n" + "=" * 40)
    print("Quick check complete!")

if __name__ == "__main__":
    main()