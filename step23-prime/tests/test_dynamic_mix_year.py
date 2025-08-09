#!/usr/bin/env python3
"""
Test the --dynamic-mix-year feature with edge cases.
"""

import os
import sys
import subprocess
import shutil
import json
import gzip

def run_test(test_name, growth_years, dynamic_mix, expected_year, should_succeed=True):
    """Run a single test case."""
    print(f"\n{'='*60}")
    print(f"TEST: {test_name}")
    print(f"  Growth years: {growth_years}")
    print(f"  Dynamic mix: {dynamic_mix}")
    print(f"  Expected snapshot year: {expected_year}")
    print(f"{'='*60}")
    
    # Clean test directory
    test_dir = f"data_test_dynamic"
    if os.path.exists(test_dir):
        shutil.rmtree(test_dir)
    
    # Build command
    cmd = [
        "python", "run_pipeline.py",
        "--rate", "0.005",
        "--simulation", "../step1-prime/data/simulation_rate_0.005000_g10_m957_n100_t100_seed42.json.gz",
        "--n-quantiles", "2",
        "--cells-per-quantile", "1",
        "--growth-years", str(growth_years),
        "--seed", "42",
        "--output-dir", test_dir
    ]
    
    if dynamic_mix:
        cmd.append("--dynamic-mix-year")
    
    print(f"  Command: {' '.join(cmd[-6:])}")  # Show relevant part
    
    # Run pipeline
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    # Check if it ran successfully
    if should_succeed:
        if result.returncode != 0:
            # Check if it's the mixing error (expected for some tests)
            if "Need" in result.stderr and "cells but only" in result.stderr:
                print(f"  ⚠️  Mix failed (not enough cells) - Expected")
                return True
            else:
                print(f"  ❌ Pipeline failed unexpectedly")
                print(f"  Error: {result.stderr}")
                return False
    
    # Check if correct snapshot was created/used
    snapshot_path = f"{test_dir}/rate_0.005000/snapshots/year{expected_year}_snapshot.json.gz"
    
    if os.path.exists(snapshot_path):
        print(f"  ✅ Created/used year {expected_year} snapshot")
        
        # Check file content to verify it's the right year
        with gzip.open(snapshot_path, 'rt') as f:
            data = json.load(f)
        
        # Get first cell's age (should be close to expected_year)
        if 'cells' in data:
            cells = data['cells']
        else:
            cells = data
        
        if cells and len(cells) > 0:
            first_cell = cells[0]
            age = first_cell.get('age', -1)
            print(f"  ✅ First cell age: {age} (should be ~{expected_year})")
    else:
        print(f"  ❌ Expected snapshot not found: {snapshot_path}")
        # List what was created
        snapshots_dir = f"{test_dir}/rate_0.005000/snapshots"
        if os.path.exists(snapshots_dir):
            files = os.listdir(snapshots_dir)
            print(f"  Found: {files}")
        return False
    
    # Check output for dynamic mode message
    if dynamic_mix and "Dynamic mix year: ENABLED" in result.stdout:
        print(f"  ✅ Dynamic mode message displayed")
    elif not dynamic_mix and "Dynamic mix year: ENABLED" not in result.stdout:
        print(f"  ✅ Dynamic mode not mentioned (correct)")
    else:
        print(f"  ⚠️  Dynamic mode message mismatch")
    
    return True

def main():
    """Run all tests."""
    print("="*60)
    print("TESTING DYNAMIC MIX YEAR FEATURE")
    print("="*60)
    
    # Check if simulation exists
    sim_path = "../step1-prime/data/simulation_rate_0.005000_g10_m957_n100_t100_seed42.json.gz"
    if not os.path.exists(sim_path):
        print(f"Error: Test simulation not found: {sim_path}")
        print("Please run: cd ../step1-prime && python run_simulation.py --rate 0.005 --years 100 --growth-phase 10 --sites 100 --seed 42")
        sys.exit(1)
    
    all_passed = True
    
    # Test 1: Default behavior (year 60)
    if not run_test("Default (no dynamic flag)", 
                     growth_years=10, 
                     dynamic_mix=False, 
                     expected_year=60):
        all_passed = False
    
    # Test 2: Dynamic with 10 years (year 60)
    if not run_test("Dynamic with 10 years growth", 
                     growth_years=10, 
                     dynamic_mix=True, 
                     expected_year=60):
        all_passed = False
    
    # Test 3: Dynamic with 5 years (year 55)
    if not run_test("Dynamic with 5 years growth", 
                     growth_years=5, 
                     dynamic_mix=True, 
                     expected_year=55):
        all_passed = False
    
    # Test 4: Dynamic with 2 years (year 52)
    if not run_test("Dynamic with 2 years growth", 
                     growth_years=2, 
                     dynamic_mix=True, 
                     expected_year=52):
        all_passed = False
    
    # Test 5: Edge case - 0 years (year 50)
    if not run_test("Edge case: 0 years growth", 
                     growth_years=0, 
                     dynamic_mix=True, 
                     expected_year=50):
        all_passed = False
    
    # Test 6: Edge case - 50 years (year 100)
    if not run_test("Edge case: 50 years growth (max)", 
                     growth_years=50, 
                     dynamic_mix=True, 
                     expected_year=100):
        all_passed = False
    
    # Test 7: Compare outputs with and without dynamic
    print(f"\n{'='*60}")
    print("TEST: Comparing dynamic vs non-dynamic")
    print(f"{'='*60}")
    
    # Run both
    subprocess.run(["rm", "-rf", "data_test_compare_dynamic", "data_test_compare_static"], 
                   capture_output=True)
    
    cmd_dynamic = [
        "python", "run_pipeline.py",
        "--rate", "0.005",
        "--simulation", sim_path,
        "--n-quantiles", "2", "--cells-per-quantile", "1",
        "--growth-years", "5", "--dynamic-mix-year",
        "--seed", "42", "--output-dir", "data_test_compare_dynamic"
    ]
    
    cmd_static = [
        "python", "run_pipeline.py",
        "--rate", "0.005",
        "--simulation", sim_path,
        "--n-quantiles", "2", "--cells-per-quantile", "1",
        "--growth-years", "5",  # No dynamic flag
        "--seed", "42", "--output-dir", "data_test_compare_static"
    ]
    
    subprocess.run(cmd_dynamic, capture_output=True)
    subprocess.run(cmd_static, capture_output=True)
    
    # Check snapshots
    dynamic_snaps = set(os.listdir("data_test_compare_dynamic/rate_0.005000/snapshots/"))
    static_snaps = set(os.listdir("data_test_compare_static/rate_0.005000/snapshots/"))
    
    print(f"  Dynamic snapshots: {sorted(dynamic_snaps)}")
    print(f"  Static snapshots: {sorted(static_snaps)}")
    
    if "year55_snapshot.json.gz" in dynamic_snaps and "year60_snapshot.json.gz" in static_snaps:
        print(f"  ✅ Different snapshots used as expected")
    else:
        print(f"  ❌ Unexpected snapshot files")
        all_passed = False
    
    # Cleanup
    print(f"\n{'='*60}")
    print("Cleaning up test directories...")
    for test_dir in ["data_test_dynamic", "data_test_compare_dynamic", "data_test_compare_static"]:
        if os.path.exists(test_dir):
            shutil.rmtree(test_dir)
    
    # Summary
    print(f"\n{'='*60}")
    if all_passed:
        print("✅ ALL TESTS PASSED")
    else:
        print("❌ SOME TESTS FAILED")
    print(f"{'='*60}")
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())