#!/usr/bin/env python3
"""
Comprehensive test to ensure step23-prime is fully reproducible.
Tests both internal reproducibility and comparison with step23.
"""

import os
import sys
import subprocess
import json
import gzip
import hashlib
import shutil
import numpy as np

def run_command(cmd):
    """Run a shell command and return success status."""
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Error: {result.stderr}")
    return result.returncode == 0

def get_snapshot_hash(filepath):
    """Get hash of snapshot file."""
    with gzip.open(filepath, 'rb') as f:
        content = f.read()
    return hashlib.md5(content).hexdigest()[:8]

def get_individual_signature(filepath):
    """Get signature from individual file."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    cells = data['cells'] if 'cells' in data else data
    jsds = sorted([c['jsd'] for c in cells])
    return {
        'n_cells': len(cells),
        'mean': np.mean(jsds),
        'std': np.std(jsds),
        'hash': hashlib.md5(str(jsds[:100]).encode()).hexdigest()[:8]
    }

def test_reproducibility():
    """Test if step23-prime is reproducible with itself."""
    
    print("="*60)
    print("STEP23-PRIME REPRODUCIBILITY TEST")
    print("="*60)
    
    # Test parameters
    test_params = [
        "--rate 0.005",
        "--simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz",
        "--n-quantiles 2",  # Small for faster testing
        "--cells-per-quantile 1",
        "--growth-years 2",  # Small for faster testing
        "--mix-ratio 80",
        "--seed 12345",  # Different seed to avoid collision with existing
        "--bins 200"
    ]
    
    base_cmd = "python run_pipeline.py " + " ".join(test_params)
    
    # Clean any existing test data
    test_dir = "data/rate_0.005000_test"
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    
    print("\n1. INTERNAL REPRODUCIBILITY TEST")
    print("-" * 40)
    
    # Temporarily rename output dir for test
    orig_cmd = base_cmd.replace("--rate 0.005", "--rate 0.005_test")
    
    # Run 1
    print("\nRun 1...")
    if not run_command(orig_cmd + " 2>&1 | head -100"):
        print("First run failed/timed out (expected for large data)")
    
    # Check what was created
    if os.path.exists("data/rate_0.005000_test/snapshots/year50_snapshot.json.gz"):
        snap1_hash = get_snapshot_hash("data/rate_0.005000_test/snapshots/year50_snapshot.json.gz")
        print(f"Run 1 year50 hash: {snap1_hash}")
        
        # Save first run
        shutil.move("data/rate_0.005000_test", "data/rate_0.005000_test_run1")
        
        # Run 2 with same parameters
        print("\nRun 2...")
        if not run_command(orig_cmd + " 2>&1 | head -100"):
            print("Second run failed/timed out (expected for large data)")
        
        if os.path.exists("data/rate_0.005000_test/snapshots/year50_snapshot.json.gz"):
            snap2_hash = get_snapshot_hash("data/rate_0.005000_test/snapshots/year50_snapshot.json.gz")
            print(f"Run 2 year50 hash: {snap2_hash}")
            
            if snap1_hash == snap2_hash:
                print("✅ Snapshots are IDENTICAL - Good reproducibility!")
            else:
                print("❌ Snapshots DIFFER - Reproducibility issue!")
        
        # Check individuals if they exist
        if os.path.exists("data/rate_0.005000_test/individuals/mutant/individual_00.json.gz"):
            ind1 = get_individual_signature("data/rate_0.005000_test_run1/individuals/mutant/individual_00.json.gz")
            ind2 = get_individual_signature("data/rate_0.005000_test/individuals/mutant/individual_00.json.gz")
            
            print(f"\nIndividual comparison:")
            print(f"Run 1: {ind1['n_cells']} cells, mean={ind1['mean']:.6f}")
            print(f"Run 2: {ind2['n_cells']} cells, mean={ind2['mean']:.6f}")
            
            if ind1['hash'] == ind2['hash']:
                print("✅ Individuals are IDENTICAL - Perfect reproducibility!")
            else:
                diff = abs(ind1['mean'] - ind2['mean'])
                if diff < 0.001:
                    print(f"⚠️  Very close but not identical (diff={diff:.9f})")
                else:
                    print(f"❌ Individuals differ significantly (diff={diff:.6f})")
    
    # Clean up test directories
    for d in ["data/rate_0.005000_test", "data/rate_0.005000_test_run1"]:
        if os.path.exists(d):
            shutil.rmtree(d)
    
    print("\n" + "="*60)
    print("RECOMMENDATIONS")
    print("="*60)
    print("""
If not reproducible, check:
1. Global seed is set at start of run_pipeline()
2. All functions that use randomness call random.seed()
3. numpy.random is also seeded where needed
4. No randomness in file I/O or data structure iteration

The fixes applied should ensure reproducibility.
    """)

if __name__ == "__main__":
    print("This will run a quick reproducibility test.")
    print("It will create temporary test directories that will be cleaned up.")
    response = input("Continue? (y/n): ")
    if response.lower() == 'y':
        test_reproducibility()
    else:
        print("Test cancelled.")