#!/usr/bin/env python3
"""
Test that timestamp-based path generation works correctly.
"""

import os
import sys
import re
import time
from datetime import datetime

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from path_utils import generate_step23_output_dir


def test_timestamp_format():
    """Test that paths use timestamp instead of hash."""
    print("\nTest: Timestamp Path Generation")
    print("=" * 50)
    
    # Create mock args object
    class Args:
        rate = 0.005
        gene_rate_groups = None
        first_snapshot = 50
        second_snapshot = 60
        individual_growth_phase = 7
        n_quantiles = 10
        cells_per_quantile = 3
        mix_ratio = 80
        uniform_mixing = True
        normalize_size = False
        seed = 42
        output_dir = "data"
    
    # Mock sim_params
    sim_params = {
        'growth_phase': 13,
        'n_sites': 1000,
        'sim_years': 100
    }
    
    # Generate path
    before_time = datetime.now().strftime("%Y%m%d%H%M%S")
    output_path = generate_step23_output_dir(Args(), sim_params)
    time.sleep(0.01)  # Small delay to ensure different timestamp if called again
    after_time = datetime.now().strftime("%Y%m%d%H%M%S")
    
    print(f"Generated path: {output_path}")
    
    # Extract timestamp from path
    # Pattern: ...seed42-YYYYMMDDHHMMSS
    pattern = r"seed\d+-(\d{14})$"
    match = re.search(pattern, output_path)
    
    assert match, f"Path doesn't contain timestamp: {output_path}"
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
    print("  ✓ Timestamp in expected range")
    
    return True


def test_uniqueness():
    """Test that successive calls generate unique paths."""
    print("\nTest: Path Uniqueness")
    print("=" * 50)
    
    class Args:
        rate = 0.005
        gene_rate_groups = None
        first_snapshot = 50
        second_snapshot = 60
        individual_growth_phase = 7
        n_quantiles = 10
        cells_per_quantile = 3
        mix_ratio = 80
        uniform_mixing = False
        normalize_size = False
        seed = 42
        output_dir = "data"
    
    sim_params = {
        'growth_phase': 13,
        'n_sites': 1000,
        'sim_years': 100
    }
    
    # Generate multiple paths
    paths = []
    for i in range(3):
        path = generate_step23_output_dir(Args(), sim_params)
        paths.append(path)
        print(f"  Path {i+1}: {path}")
        time.sleep(1.1)  # Sleep > 1 second to ensure different timestamp
    
    # Check all paths are unique
    assert len(paths) == len(set(paths)), "Paths are not unique!"
    print("  ✓ All paths are unique")
    
    return True


def test_sorting():
    """Test that timestamp-based paths sort chronologically."""
    print("\nTest: Chronological Sorting")
    print("=" * 50)
    
    # Create sample paths with timestamps
    base = "data/rate_0.00500-grow13-sites1000-years100/snap50to60-growth7-quant10x3-mix80-seed42-"
    
    paths = [
        base + "20240115143022",  # 2024-01-15 14:30:22
        base + "20240115092515",  # 2024-01-15 09:25:15
        base + "20240116083045",  # 2024-01-16 08:30:45
        base + "20240114235959",  # 2024-01-14 23:59:59
    ]
    
    sorted_paths = sorted(paths)
    
    print("Original order:")
    for p in paths:
        print(f"  {p.split('-')[-1]}")
    
    print("\nSorted order:")
    for p in sorted_paths:
        print(f"  {p.split('-')[-1]}")
    
    # Verify chronological order
    expected = [
        base + "20240114235959",
        base + "20240115092515",
        base + "20240115143022",
        base + "20240116083045",
    ]
    
    assert sorted_paths == expected, "Paths don't sort chronologically"
    print("\n  ✓ Paths sort chronologically")
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing Timestamp-based Path Generation")
    print("=" * 60)
    
    tests = [
        test_timestamp_format,
        test_uniqueness,
        test_sorting
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
        print("\n🎉 All timestamp tests passed!")
        print("   ✅ YYYYMMDDHHMMSS format working")
        print("   ✅ Paths are unique")
        print("   ✅ Chronological sorting works")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())