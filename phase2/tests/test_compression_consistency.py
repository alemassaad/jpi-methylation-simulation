#!/usr/bin/env python3
"""
Test compression consistency across all pipeline outputs.
"""

import os
import sys
import json
import gzip
import glob
import tempfile
import shutil
import subprocess

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish

def create_test_simulation(filepath, compress=False):
    """Create a minimal test simulation file."""
    # Create a small PetriDish with history
    petri = PetriDish(rate=0.005, growth_phase=2, calculate_cell_jsds=True)
    
    # Simulate a few years
    for year in range(5):
        petri.divide_cells()
        petri.methylate_cells()
        petri.random_cull_cells()
    
    # Save history
    petri.save_history(filepath, compress=compress)
    print(f"Created test simulation: {filepath}")

def check_file_formats(directory, expected_extension):
    """Check that all files in directory tree have expected extension."""
    issues = []
    
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.startswith('.'):  # Skip hidden files
                continue
            if not file.endswith(expected_extension):
                if file.endswith('.json') or file.endswith('.json.gz'):
                    issues.append(os.path.join(root, file))
    
    return issues

def run_pipeline_test(sim_file, test_name, extra_args=""):
    """Run the pipeline and check output consistency."""
    print(f"\n{'='*60}")
    print(f"Test: {test_name}")
    print(f"Input: {sim_file}")
    print(f"{'='*60}")
    
    # Create temp output dir
    output_dir = f"test_compression_{test_name}"
    
    # Clean up if exists
    if os.path.exists(output_dir):
        shutil.rmtree(output_dir)
    
    # Run pipeline
    cmd = f"python3 ../run_pipeline.py --simulation {sim_file} --rate 0.005 " \
          f"--first-snapshot 2 --second-snapshot 3 " \
          f"--growth-phase 1 --n-quantiles 2 --cells-per-quantile 2 " \
          f"--output-dir {output_dir} --seed 42 {extra_args}"
    
    print(f"Running: {cmd}")
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.returncode != 0:
        print(f"❌ Pipeline failed!")
        print(f"Error: {result.stderr}")
        return False
    
    # Determine expected extension
    if "--no-compress" in extra_args:
        expected_ext = ".json"
        unexpected_ext = ".json.gz"
    elif sim_file.endswith('.gz'):
        expected_ext = ".json.gz"
        unexpected_ext = ".json"
    else:
        expected_ext = ".json"
        unexpected_ext = ".json.gz"
    
    print(f"\nExpected extension: {expected_ext}")
    
    # Find the actual output directory (with timestamp)
    pattern = os.path.join(output_dir, "rate_*/snap*")
    output_dirs = glob.glob(pattern)
    
    if not output_dirs:
        print(f"❌ No output directory found matching {pattern}")
        return False
    
    actual_output = output_dirs[0]
    print(f"Output directory: {actual_output}")
    
    # Check all files have consistent format
    wrong_format = check_file_formats(actual_output, expected_ext)
    
    if wrong_format:
        print(f"\n❌ Found files with wrong format:")
        for file in wrong_format:
            print(f"  - {file}")
        return False
    
    # Count files of each type
    all_json = glob.glob(os.path.join(actual_output, "**/*.json"), recursive=True)
    all_gz = glob.glob(os.path.join(actual_output, "**/*.json.gz"), recursive=True)
    
    print(f"\nFile counts:")
    print(f"  .json files: {len(all_json)}")
    print(f"  .json.gz files: {len(all_gz)}")
    
    # Check specific directories
    for subdir in ["snapshots", "individuals/mutant", "individuals/control1", "individuals/control2"]:
        dir_path = os.path.join(actual_output, subdir)
        if os.path.exists(dir_path):
            json_files = glob.glob(os.path.join(dir_path, "*.json"))
            gz_files = glob.glob(os.path.join(dir_path, "*.json.gz"))
            print(f"\n  {subdir}:")
            print(f"    .json: {len(json_files)}")
            print(f"    .json.gz: {len(gz_files)}")
            
            # Verify consistency
            if expected_ext == ".json" and gz_files:
                print(f"    ❌ Unexpected .json.gz files in {subdir}")
                return False
            elif expected_ext == ".json.gz" and json_files:
                print(f"    ❌ Unexpected .json files in {subdir}")
                return False
    
    # Clean up
    shutil.rmtree(output_dir)
    
    print(f"\n✅ Test '{test_name}' passed - all files use {expected_ext}")
    return True

def main():
    """Run all compression consistency tests."""
    print("Testing Compression Consistency")
    print("================================")
    
    # Create test simulations
    uncompressed_sim = "test_sim_uncompressed.json"
    compressed_sim = "test_sim_compressed.json.gz"
    
    create_test_simulation(uncompressed_sim, compress=False)
    create_test_simulation(compressed_sim, compress=True)
    
    tests_passed = 0
    tests_failed = 0
    
    # Test 1: Uncompressed input -> uncompressed output
    if run_pipeline_test(uncompressed_sim, "uncompressed_input"):
        tests_passed += 1
    else:
        tests_failed += 1
    
    # Test 2: Compressed input -> compressed output
    if run_pipeline_test(compressed_sim, "compressed_input"):
        tests_passed += 1
    else:
        tests_failed += 1
    
    # Test 3: Compressed input with --no-compress -> uncompressed output
    if run_pipeline_test(compressed_sim, "compressed_with_flag", "--no-compress"):
        tests_passed += 1
    else:
        tests_failed += 1
    
    # Test 4: Uncompressed input with --no-compress -> uncompressed output
    if run_pipeline_test(uncompressed_sim, "uncompressed_with_flag", "--no-compress"):
        tests_passed += 1
    else:
        tests_failed += 1
    
    # Clean up test files
    if os.path.exists(uncompressed_sim):
        os.remove(uncompressed_sim)
    if os.path.exists(compressed_sim):
        os.remove(compressed_sim)
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Test Summary")
    print(f"{'='*60}")
    print(f"✅ Passed: {tests_passed}")
    print(f"❌ Failed: {tests_failed}")
    
    if tests_failed == 0:
        print("\n🎉 All compression consistency tests passed!")
        return 0
    else:
        print(f"\n⚠️  {tests_failed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())