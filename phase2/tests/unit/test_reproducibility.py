#!/usr/bin/env python3
"""
Test if phase2 produces reproducible results when run twice with the same seed.
"""

import json
import gzip
import numpy as np
import hashlib
import os

def get_individual_signature(filepath):
    """Get a signature from an individual file."""
    try:
        with gzip.open(filepath, 'rt') as f:
            data = json.load(f)
        
        cells = data['cells'] if 'cells' in data else data
        
        # Get sorted JSD values as signature
        jsds = sorted([c['cell_jsd'] for c in cells])
        
        # Create a hash of the JSD values
        jsd_str = ','.join([f"{j:.9f}" for j in jsds[:100]])  # First 100 for speed
        jsd_hash = hashlib.md5(jsd_str.encode()).hexdigest()[:8]
        
        return {
            'n_cells': len(cells),
            'mean_cell_jsd': np.mean(jsds),
            'std_cell_jsd': np.std(jsds),
            'hash': jsd_hash
        }
    except:
        return None

def test_reproducibility():
    """
    Compare two runs of phase2 to see if they're reproducible.
    Assumes you've run the pipeline twice with same parameters.
    """
    
    print("="*60)
    print("PHASE 2 REPRODUCIBILITY TEST")
    print("="*60)
    print("\nThis tests if running phase2 twice with the same seed")
    print("produces identical results.\n")
    
    # Test with first mutant individual
    run1_dir = "data/rate_0.005000"
    run2_dir = "data_run2/rate_0.005000"  # Assume second run saved here
    
    if not os.path.exists(run2_dir):
        print("❌ Second run not found!")
        print("\nTo test reproducibility:")
        print("1. Run pipeline once:")
        print("   python run_pipeline.py --rate 0.005 --simulation ... --seed 42")
        print("2. Move results:")
        print("   mv data/rate_0.005000 data_run2/rate_0.005000")
        print("3. Run pipeline again with SAME parameters:")
        print("   python run_pipeline.py --rate 0.005 --simulation ... --seed 42")
        print("4. Run this test again")
        return
    
    print("Comparing individuals from two runs...\n")
    
    # Compare first few individuals from each group
    all_match = True
    
    for group in ['mutant', 'control1', 'control2']:
        print(f"--- {group.upper()} ---")
        
        for i in range(3):  # Check first 3 individuals
            file1 = f"{run1_dir}/individuals/{group}/individual_{i:02d}.json.gz"
            file2 = f"{run2_dir}/individuals/{group}/individual_{i:02d}.json.gz"
            
            sig1 = get_individual_signature(file1)
            sig2 = get_individual_signature(file2)
            
            if sig1 and sig2:
                match = sig1['hash'] == sig2['hash']
                symbol = "✅" if match else "❌"
                
                print(f"  Individual {i:02d}: {symbol}")
                print(f"    Run 1: {sig1['n_cells']} cells, mean={sig1['mean_cell_jsd']:.6f}, hash={sig1['hash']}")
                print(f"    Run 2: {sig2['n_cells']} cells, mean={sig2['mean_cell_jsd']:.6f}, hash={sig2['hash']}")
                
                if not match:
                    all_match = False
                    diff = abs(sig1['mean_cell_jsd'] - sig2['mean_cell_jsd'])
                    print(f"    Difference: {diff:.9f}")
            else:
                print(f"  Individual {i:02d}: ❌ File not found")
                all_match = False
        print()
    
    # Summary
    print("="*60)
    print("REPRODUCIBILITY RESULT")
    print("="*60)
    
    if all_match:
        print("✅ REPRODUCIBLE! Same seed produces identical results.")
    else:
        print("❌ NOT REPRODUCIBLE! Same seed produces different results.")
        print("\nThis means phase2 has a reproducibility issue.")
        print("The problem is likely:")
        print("1. No global random seed set at start")
        print("2. Random state not controlled between function calls")
        print("3. Some functions not seeding properly")

if __name__ == "__main__":
    test_reproducibility()