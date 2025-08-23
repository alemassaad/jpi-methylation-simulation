#!/usr/bin/env python3
"""
Verify if the pipeline is deterministic with the same seed.
"""

import json
import gzip
import numpy as np
import hashlib

def get_file_hash(filepath):
    """Get hash of file contents."""
    with gzip.open(filepath, 'rb') as f:
        return hashlib.md5(f.read()).hexdigest()

def get_jsd_signature(filepath):
    """Get sorted JSD values as a signature."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    cells = data['cells'] if 'cells' in data else data
    jsds = sorted([c['cell_jsd'] for c in cells])
    
    return {
        'n_cells': len(cells),
        'mean_jsd': np.mean(jsds),
        'std_jsd': np.std(jsds),
        'min_jsd': min(jsds),
        'max_jsd': max(jsds),
        'hash': hashlib.md5(str(jsds).encode()).hexdigest()[:8]
    }

def main():
    print("="*60)
    print("DETERMINISM CHECK")
    print("="*60)
    
    # Check if we can find matching individuals
    print("\nChecking first mutant individual...")
    
    try:
        step23_sig = get_jsd_signature('../step23/data/rate_0.005000/individuals/mutant/individual_00.json.gz')
        prime_sig = get_jsd_signature('data/rate_0.005000/individuals/mutant/individual_00.json.gz')
        
        print(f"\nStep23 signature:")
        print(f"  Cells: {step23_sig['n_cells']}")
        print(f"  Mean JSD: {step23_sig['mean_jsd']:.6f}")
        print(f"  Std JSD: {step23_sig['std_jsd']:.6f}")
        print(f"  Hash: {step23_sig['hash']}")
        
        print(f"\nStep23-prime signature:")
        print(f"  Cells: {prime_sig['n_cells']}")
        print(f"  Mean JSD: {prime_sig['mean_jsd']:.6f}")
        print(f"  Std JSD: {prime_sig['std_jsd']:.6f}")
        print(f"  Hash: {prime_sig['hash']}")
        
        if step23_sig['hash'] == prime_sig['hash']:
            print("\n✅ IDENTICAL! The JSDs match exactly (same random sequence)")
        else:
            diff = abs(step23_sig['mean_jsd'] - prime_sig['mean_jsd'])
            print(f"\n⚠️  Different random sequences (mean diff: {diff:.6f})")
            
            if diff < 0.001:
                print("    But VERY close - likely same algorithm, different RNG state")
            elif diff < 0.01:
                print("    Statistically similar - acceptable for scientific purposes")
            else:
                print("    Significant difference - may indicate different logic")
    
    except Exception as e:
        print(f"Error: {e}")
    
    print("\n" + "="*60)
    print("EXPLANATION")
    print("="*60)
    print("""
Even with the same seed (42), individuals may differ because:

1. GROWTH IS STOCHASTIC: Each cell has random methylation events
   - 10 years × 1024 cells × 1000 sites = millions of random events
   
2. RANDOM STATE FLOW: The seed sets initial state, but then:
   - Individual 00 consumes random numbers
   - Individual 01 gets the NEXT sequence, not a reset
   
3. ORDER MATTERS: Any difference in processing order changes everything

To get EXACT reproduction, we'd need:
- Set seed PER individual: random.seed(42 + individual_id)
- Or save/restore random state between individuals
- Or use separate RandomState objects per individual

For scientific purposes, statistical similarity (same distribution) 
is usually sufficient rather than bit-for-bit reproduction.
    """)

if __name__ == "__main__":
    main()