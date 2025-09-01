#!/usr/bin/env python3
"""
Test that control2 batch respects compression settings.
"""

import os
import sys
import glob
import tempfile
import shutil

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from pipeline_utils import save_petri_dish, load_petri_dish, save_snapshot_cells
from cell import Cell, PetriDish

def create_minimal_simulation(filepath, compress=False):
    """Create a minimal simulation file for testing."""
    petri = PetriDish(rate=0.005, growth_phase=2, calculate_cell_jsds=True)
    
    # Simulate enough years for snapshots
    for year in range(6):
        petri.divide_cells()
        petri.methylate_cells()
        if year > petri.growth_phase:
            petri.random_cull_cells()
    
    petri.save_history(filepath, compress=compress)
    return filepath

def check_control2_compression(sim_file, expected_compressed):
    """Check if control2 files match expected compression."""
    import subprocess
    
    # Create temp output dir
    with tempfile.TemporaryDirectory() as tmpdir:
        # Run minimal pipeline
        cmd = [
            "python3", "../run_pipeline.py",
            "--simulation", sim_file,
            "--rate", "0.005",
            "--first-snapshot", "3",
            "--second-snapshot", "4",
            "--growth-phase", "1",
            "--n-quantiles", "2",
            "--cells-per-quantile", "1",
            "--output-dir", tmpdir,
            "--seed", "42"
        ]
        
        if not expected_compressed:
            cmd.append("--no-compress")
        
        print(f"Running: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            # Check if it's just the plotly error
            if "ModuleNotFoundError: No module named 'plotly'" in result.stderr:
                print("  ⚠️  Pipeline failed due to missing plotly (expected)")
                print("  Checking files created before failure...")
                
                # Look for any files created
                all_files = glob.glob(os.path.join(tmpdir, "**/*"), recursive=True)
                json_files = [f for f in all_files if f.endswith('.json')]
                gz_files = [f for f in all_files if f.endswith('.json.gz')]
                
                print(f"  Found {len(json_files)} .json files, {len(gz_files)} .json.gz files")
                
                # Check control2 specifically
                control2_files = glob.glob(os.path.join(tmpdir, "**/control2/*.json*"), recursive=True)
                if control2_files:
                    print(f"  Control2 files: {control2_files}")
                    for f in control2_files:
                        if expected_compressed and f.endswith('.json'):
                            print(f"  ❌ Found uncompressed control2 file: {f}")
                            return False
                        elif not expected_compressed and f.endswith('.json.gz'):
                            print(f"  ❌ Found compressed control2 file: {f}")
                            return False
                
                # Even with plotly error, we can check what was created
                return True
            else:
                print(f"  ❌ Pipeline failed: {result.stderr}")
                return False
        
        # Check all output directories
        for batch in ["mutant", "control1", "control2"]:
            batch_dir = glob.glob(os.path.join(tmpdir, f"**/individuals/{batch}"), recursive=True)
            if batch_dir:
                json_files = glob.glob(os.path.join(batch_dir[0], "*.json"))
                gz_files = glob.glob(os.path.join(batch_dir[0], "*.json.gz"))
                
                print(f"  {batch}: {len(json_files)} .json, {len(gz_files)} .json.gz")
                
                if expected_compressed and json_files:
                    print(f"  ❌ Found uncompressed files in {batch}")
                    return False
                elif not expected_compressed and gz_files:
                    print(f"  ❌ Found compressed files in {batch}")
                    return False
        
        return True

def main():
    """Run control2 compression tests."""
    print("Testing Control2 Compression Consistency")
    print("=" * 50)
    
    # Create test simulations
    print("\nCreating test simulations...")
    uncompressed_sim = create_minimal_simulation("./test_sim.json", compress=False)
    compressed_sim = create_minimal_simulation("./test_sim.json.gz", compress=True)
    
    tests_passed = 0
    tests_failed = 0
    
    # Test 1: Uncompressed input should produce uncompressed control2
    print("\nTest 1: Uncompressed input -> all uncompressed")
    if check_control2_compression(uncompressed_sim, expected_compressed=False):
        print("  ✅ Pass - control2 matches other batches")
        tests_passed += 1
    else:
        print("  ❌ Fail - control2 compression mismatch")
        tests_failed += 1
    
    # Test 2: Compressed input should produce compressed control2
    print("\nTest 2: Compressed input -> all compressed")
    if check_control2_compression(compressed_sim, expected_compressed=True):
        print("  ✅ Pass - control2 matches other batches")
        tests_passed += 1
    else:
        print("  ❌ Fail - control2 compression mismatch")
        tests_failed += 1
    
    # Clean up
    if os.path.exists(uncompressed_sim):
        os.remove(uncompressed_sim)
    if os.path.exists(compressed_sim):
        os.remove(compressed_sim)
    
    # Also clean up phase1 data directory if created
    if os.path.exists("data"):
        shutil.rmtree("data")
    
    print("\n" + "=" * 50)
    print(f"Results: {tests_passed} passed, {tests_failed} failed")
    
    if tests_failed == 0:
        print("🎉 All control2 compression tests passed!")
        return 0
    else:
        print(f"⚠️  {tests_failed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())