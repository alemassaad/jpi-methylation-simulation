#!/usr/bin/env python3
"""
Test script for --no-compress functionality.
Tests both phase1 and phase2 with compressed and uncompressed outputs.
"""

import os
import sys
import subprocess
import json
import gzip
import glob
import tempfile
import shutil

def run_command(cmd, cwd=None):
    """Run a command and return stdout, stderr, and return code."""
    print(f"\n  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return result.stdout, result.stderr, result.returncode

def check_file_format(filepath, should_be_compressed):
    """Check if a file is in the expected format (compressed or uncompressed)."""
    if not os.path.exists(filepath):
        return False, f"File does not exist: {filepath}"
    
    # Try to open as compressed
    try:
        with gzip.open(filepath, 'rt') as f:
            data = json.load(f)
        is_compressed = True
    except (gzip.BadGzipFile, UnicodeDecodeError):
        # Not compressed, try regular JSON
        try:
            with open(filepath, 'r') as f:
                data = json.load(f)
            is_compressed = False
        except json.JSONDecodeError:
            return False, f"File is not valid JSON: {filepath}"
    except json.JSONDecodeError:
        return False, f"File is not valid JSON (compressed): {filepath}"
    
    if is_compressed != should_be_compressed:
        return False, f"File compression mismatch: {filepath} (compressed={is_compressed}, expected={should_be_compressed})"
    
    return True, f"File format correct: {filepath} (compressed={is_compressed})"

def test_phase1_compressed():
    """Test phase1 with default compressed output."""
    print("\n" + "="*60)
    print("TEST 1: Phase 1 - Compressed Output (default)")
    print("="*60)
    
    # Run simulation with default compression
    cmd = [
        "python3", "run_simulation.py",
        "--rate", "0.005",
        "--years", "10",
        "--growth-phase", "3",
        "--seed", "123"
    ]
    
    stdout, stderr, returncode = run_command(cmd, cwd="phase1")
    if returncode != 0:
        print(f"  ❌ Command failed with return code {returncode}")
        print(f"  Error: {stderr}")
        return False
    
    # Find the output file
    pattern = "phase1/data/rate_*/grow3-sites1000-years10-seed123-*/simulation.json.gz"
    files = glob.glob(pattern)
    if not files:
        print(f"  ❌ No output file found matching pattern: {pattern}")
        return False
    
    filepath = files[0]
    success, msg = check_file_format(filepath, should_be_compressed=True)
    print(f"  {msg}")
    
    # Check file size
    size_mb = os.path.getsize(filepath) / (1024 * 1024)
    print(f"  File size: {size_mb:.3f} MB (compressed)")
    
    # Clean up
    shutil.rmtree(os.path.dirname(filepath))
    
    return success

def test_phase1_uncompressed():
    """Test phase1 with --no-compress flag."""
    print("\n" + "="*60)
    print("TEST 2: Phase 1 - Uncompressed Output (--no-compress)")
    print("="*60)
    
    # Run simulation with --no-compress
    cmd = [
        "python3", "run_simulation.py",
        "--rate", "0.005",
        "--years", "10",
        "--growth-phase", "3",
        "--seed", "124",
        "--no-compress"
    ]
    
    stdout, stderr, returncode = run_command(cmd, cwd="phase1")
    if returncode != 0:
        print(f"  ❌ Command failed with return code {returncode}")
        print(f"  Error: {stderr}")
        return False
    
    # Find the output file
    pattern = "phase1/data/rate_*/grow3-sites1000-years10-seed124-*/simulation.json"
    files = glob.glob(pattern)
    if not files:
        print(f"  ❌ No output file found matching pattern: {pattern}")
        return False
    
    filepath = files[0]
    success, msg = check_file_format(filepath, should_be_compressed=False)
    print(f"  {msg}")
    
    # Check file size
    size_mb = os.path.getsize(filepath) / (1024 * 1024)
    print(f"  File size: {size_mb:.3f} MB (uncompressed)")
    
    # Verify we can read the JSON directly
    try:
        with open(filepath, 'r') as f:
            data = json.load(f)
        print(f"  ✓ Successfully loaded uncompressed JSON directly")
        print(f"  ✓ Contains {len(data)} years of data")
    except Exception as e:
        print(f"  ❌ Failed to read uncompressed JSON: {e}")
        success = False
    
    # Clean up
    shutil.rmtree(os.path.dirname(filepath))
    
    return success

def test_phase2_compressed():
    """Test phase2 with compressed output."""
    print("\n" + "="*60)
    print("TEST 3: Phase 2 - Compressed Output (default)")
    print("="*60)
    
    # First create a phase1 simulation
    print("\n  Creating phase1 simulation for testing...")
    cmd = [
        "python3", "run_simulation.py",
        "--rate", "0.005",
        "--years", "20",
        "--growth-phase", "4",
        "--seed", "125"
    ]
    stdout, stderr, returncode = run_command(cmd, cwd="phase1")
    if returncode != 0:
        print(f"  ❌ Failed to create phase1 simulation")
        return False
    
    # Find the simulation file
    sim_pattern = "phase1/data/rate_*/grow4-sites1000-years20-seed125-*/simulation.json.gz"
    sim_files = glob.glob(sim_pattern)
    if not sim_files:
        print(f"  ❌ No simulation file found")
        return False
    sim_file = sim_files[0]
    
    # Run phase2 pipeline with default compression
    cmd = [
        "python3", "run_pipeline.py",
        "--rate", "0.005",
        "--simulation", f"../{sim_file}",
        "--first-snapshot", "10",
        "--second-snapshot", "15",
        "--individual-growth-phase", "2",
        "--n-quantiles", "2",
        "--cells-per-quantile", "1",
        "--seed", "126"
    ]
    
    stdout, stderr, returncode = run_command(cmd, cwd="phase2")
    if returncode != 0:
        print(f"  ❌ Pipeline failed with return code {returncode}")
        print(f"  Error: {stderr}")
        # Clean up
        shutil.rmtree(os.path.dirname(sim_file))
        return False
    
    # Check output files
    success = True
    
    # Check snapshots
    snapshot_pattern = "phase2/data/*/snap10to15-*/snapshots/*.json.gz"
    snapshot_files = glob.glob(snapshot_pattern)
    print(f"\n  Found {len(snapshot_files)} snapshot files")
    for sf in snapshot_files:
        file_success, msg = check_file_format(sf, should_be_compressed=True)
        print(f"    {msg}")
        success = success and file_success
    
    # Check individuals
    individual_pattern = "phase2/data/*/snap10to15-*/individuals/*/*.json.gz"
    individual_files = glob.glob(individual_pattern)
    print(f"\n  Found {len(individual_files)} individual files")
    if individual_files:
        # Just check first few
        for sf in individual_files[:3]:
            file_success, msg = check_file_format(sf, should_be_compressed=True)
            print(f"    {msg}")
            success = success and file_success
    
    # Clean up
    shutil.rmtree(os.path.dirname(sim_file))
    output_dirs = glob.glob("phase2/data/*/snap10to15-*")
    for d in output_dirs:
        shutil.rmtree(d)
    
    return success

def test_phase2_uncompressed():
    """Test phase2 with --no-compress flag."""
    print("\n" + "="*60)
    print("TEST 4: Phase 2 - Uncompressed Output (--no-compress)")
    print("="*60)
    
    # First create a phase1 simulation
    print("\n  Creating phase1 simulation for testing...")
    cmd = [
        "python3", "run_simulation.py",
        "--rate", "0.005",
        "--years", "20",
        "--growth-phase", "4",
        "--seed", "127",
        "--no-compress"  # Use uncompressed for phase1 too
    ]
    stdout, stderr, returncode = run_command(cmd, cwd="phase1")
    if returncode != 0:
        print(f"  ❌ Failed to create phase1 simulation")
        return False
    
    # Find the simulation file
    sim_pattern = "phase1/data/rate_*/grow4-sites1000-years20-seed127-*/simulation.json"
    sim_files = glob.glob(sim_pattern)
    if not sim_files:
        print(f"  ❌ No simulation file found")
        return False
    sim_file = sim_files[0]
    
    # Run phase2 pipeline with --no-compress
    cmd = [
        "python3", "run_pipeline.py",
        "--rate", "0.005",
        "--simulation", f"../{sim_file}",
        "--first-snapshot", "10",
        "--second-snapshot", "15",
        "--individual-growth-phase", "2",
        "--n-quantiles", "2",
        "--cells-per-quantile", "1",
        "--seed", "128",
        "--no-compress"
    ]
    
    stdout, stderr, returncode = run_command(cmd, cwd="phase2")
    if returncode != 0:
        print(f"  ❌ Pipeline failed with return code {returncode}")
        print(f"  Error: {stderr}")
        # Clean up
        shutil.rmtree(os.path.dirname(sim_file))
        return False
    
    # Check output files
    success = True
    
    # Check snapshots
    snapshot_pattern = "phase2/data/*/snap10to15-*/snapshots/*.json"
    snapshot_files = glob.glob(snapshot_pattern)
    print(f"\n  Found {len(snapshot_files)} snapshot files")
    for sf in snapshot_files:
        file_success, msg = check_file_format(sf, should_be_compressed=False)
        print(f"    {msg}")
        success = success and file_success
        
        # Verify we can read directly
        try:
            with open(sf, 'r') as f:
                data = json.load(f)
            print(f"      ✓ Can read {os.path.basename(sf)} directly as JSON")
        except Exception as e:
            print(f"      ❌ Cannot read {os.path.basename(sf)}: {e}")
            success = False
    
    # Check individuals
    individual_pattern = "phase2/data/*/snap10to15-*/individuals/*/*.json"
    individual_files = glob.glob(individual_pattern)
    print(f"\n  Found {len(individual_files)} individual files")
    if individual_files:
        # Just check first few
        for sf in individual_files[:3]:
            file_success, msg = check_file_format(sf, should_be_compressed=False)
            print(f"    {msg}")
            success = success and file_success
    
    # Clean up
    shutil.rmtree(os.path.dirname(sim_file))
    output_dirs = glob.glob("phase2/data/*/snap10to15-*")
    for d in output_dirs:
        shutil.rmtree(d)
    
    return success

def test_mixed_mode():
    """Test phase2 reading compressed phase1 data but saving uncompressed."""
    print("\n" + "="*60)
    print("TEST 5: Mixed Mode - Compressed input, Uncompressed output")
    print("="*60)
    
    # Create compressed phase1 simulation
    print("\n  Creating compressed phase1 simulation...")
    cmd = [
        "python3", "run_simulation.py",
        "--rate", "0.005",
        "--years", "20",
        "--growth-phase", "4",
        "--seed", "129"
        # No --no-compress, so it will be compressed
    ]
    stdout, stderr, returncode = run_command(cmd, cwd="phase1")
    if returncode != 0:
        print(f"  ❌ Failed to create phase1 simulation")
        return False
    
    # Find the simulation file
    sim_pattern = "phase1/data/rate_*/grow4-sites1000-years20-seed129-*/simulation.json.gz"
    sim_files = glob.glob(sim_pattern)
    if not sim_files:
        print(f"  ❌ No simulation file found")
        return False
    sim_file = sim_files[0]
    print(f"  ✓ Created compressed simulation: {sim_file}")
    
    # Run phase2 with uncompressed output
    cmd = [
        "python3", "run_pipeline.py",
        "--rate", "0.005",
        "--simulation", f"../{sim_file}",
        "--first-snapshot", "10",
        "--second-snapshot", "15",
        "--individual-growth-phase", "2",
        "--n-quantiles", "2",
        "--cells-per-quantile", "1",
        "--seed", "130",
        "--no-compress"  # Save uncompressed
    ]
    
    stdout, stderr, returncode = run_command(cmd, cwd="phase2")
    if returncode != 0:
        print(f"  ❌ Pipeline failed with return code {returncode}")
        print(f"  Error: {stderr}")
        # Clean up
        shutil.rmtree(os.path.dirname(sim_file))
        return False
    
    print(f"  ✓ Pipeline completed successfully")
    print(f"  ✓ Read compressed input, saved uncompressed output")
    
    # Verify output is uncompressed
    snapshot_pattern = "phase2/data/*/snap10to15-*/snapshots/*.json"
    snapshot_files = glob.glob(snapshot_pattern)
    success = len(snapshot_files) > 0
    
    if snapshot_files:
        file_success, msg = check_file_format(snapshot_files[0], should_be_compressed=False)
        print(f"  {msg}")
        success = success and file_success
    
    # Clean up
    shutil.rmtree(os.path.dirname(sim_file))
    output_dirs = glob.glob("phase2/data/*/snap10to15-*")
    for d in output_dirs:
        shutil.rmtree(d)
    
    return success

def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("TESTING --no-compress FUNCTIONALITY")
    print("="*60)
    
    tests = [
        ("Phase 1 Compressed", test_phase1_compressed),
        ("Phase 1 Uncompressed", test_phase1_uncompressed),
        ("Phase 2 Compressed", test_phase2_compressed),
        ("Phase 2 Uncompressed", test_phase2_uncompressed),
        ("Mixed Mode", test_mixed_mode)
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"\n  ❌ Test failed with exception: {e}")
            results.append((name, False))
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    all_passed = True
    for name, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"  {name}: {status}")
        all_passed = all_passed and success
    
    print("\n" + "="*60)
    if all_passed:
        print("✅ ALL TESTS PASSED!")
    else:
        print("❌ SOME TESTS FAILED")
    print("="*60)
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())