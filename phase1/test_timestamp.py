#!/usr/bin/env python3
"""
Test that phase1 uses timestamp instead of hash in directory names.
"""

import os
import sys
import re
import time
import tempfile
from datetime import datetime

from cell import PetriDish


def test_timestamp_in_save_history():
    """Test that save_history uses timestamp in directory name."""
    print("\nTest: Phase1 Timestamp in save_history")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create a simple PetriDish and run simulation
        petri = PetriDish(rate=0.005, growth_phase=2, seed=42, n=100)
        petri.run_simulation(t_max=5)  # Run simulation for a few years
        
        # Save history
        before_time = datetime.now().strftime("%Y%m%d%H%M%S")
        filepath = petri.save_history(directory=tmpdir, compress=False)
        after_time = datetime.now().strftime("%Y%m%d%H%M%S")
        
        print(f"Saved to: {filepath}")
        
        # Extract timestamp from path
        # Pattern: ...seed42-YYYYMMDDHHMMSS/simulation.json
        pattern = r"seed\d+-(\d{14})"
        match = re.search(pattern, filepath)
        
        assert match, f"Path doesn't contain timestamp: {filepath}"
        timestamp = match.group(1)
        
        print(f"  Extracted timestamp: {timestamp}")
        print(f"  Before time: {before_time}")
        print(f"  After time: {after_time}")
        
        # Verify timestamp format (14 digits: YYYYMMDDHHMMSS)
        assert len(timestamp) == 14, f"Timestamp wrong length: {timestamp}"
        assert timestamp.isdigit(), f"Timestamp not all digits: {timestamp}"
        
        # Verify timestamp is in valid range
        assert before_time <= timestamp <= after_time, f"Timestamp out of range"
        
        # Parse timestamp to verify it's valid
        try:
            dt = datetime.strptime(timestamp, "%Y%m%d%H%M%S")
            print(f"  Parsed as: {dt.strftime('%Y-%m-%d %H:%M:%S')}")
        except ValueError as e:
            assert False, f"Invalid timestamp format: {e}"
        
        print("  ✓ Timestamp format correct (YYYYMMDDHHMMSS)")
        print("  ✓ Timestamp is valid datetime")
        print("  ✓ File saved successfully")
        
        # Verify file exists
        assert os.path.exists(filepath), f"File not created: {filepath}"
        
    return True


def test_uniqueness_phase1():
    """Test that successive saves generate unique paths."""
    print("\nTest: Phase1 Path Uniqueness")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        paths = []
        
        for i in range(3):
            petri = PetriDish(rate=0.005, growth_phase=2, seed=42, n=100)
            petri.run_simulation(t_max=5)
            filepath = petri.save_history(directory=tmpdir, compress=False)
            paths.append(filepath)
            print(f"  Path {i+1}: {os.path.dirname(filepath).split('/')[-1]}")
            time.sleep(1.1)  # Sleep to ensure different timestamp
        
        # Extract directory names
        dirs = [os.path.dirname(p) for p in paths]
        
        # Check all directories are unique
        assert len(dirs) == len(set(dirs)), "Directories are not unique!"
        print("  ✓ All paths are unique")
        
    return True


def test_gene_rates_timestamp():
    """Test timestamp with gene-specific rates."""
    print("\nTest: Gene Rates with Timestamp")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create PetriDish with gene-specific rates
        # With n=100 and gene_size=5, we have 20 genes total
        gene_rate_groups = [(10, 0.004), (10, 0.006)]
        petri = PetriDish(
            gene_rate_groups=gene_rate_groups,
            growth_phase=2,
            seed=42,
            n=100
        )
        petri.run_simulation(t_max=5)
        
        filepath = petri.save_history(directory=tmpdir, compress=False)
        print(f"Saved to: {filepath}")
        
        # Check path contains gene_rates
        assert "gene_rates_" in filepath, "Path should contain gene_rates"
        
        # Extract timestamp
        pattern = r"seed\d+-(\d{14})"
        match = re.search(pattern, filepath)
        assert match, f"Path doesn't contain timestamp: {filepath}"
        
        timestamp = match.group(1)
        print(f"  Extracted timestamp: {timestamp}")
        
        # Verify timestamp format
        assert len(timestamp) == 14, f"Timestamp wrong length: {timestamp}"
        print("  ✓ Gene rates path uses timestamp")
        
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing Phase1 Timestamp-based Path Generation")
    print("=" * 60)
    
    tests = [
        test_timestamp_in_save_history,
        test_uniqueness_phase1,
        test_gene_rates_timestamp
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"\n  ❌ Test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Test Summary: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("\n🎉 All Phase1 timestamp tests passed!")
        print("   ✅ save_history uses YYYYMMDDHHMMSS format")
        print("   ✅ Paths are unique")
        print("   ✅ Works with gene-specific rates")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())