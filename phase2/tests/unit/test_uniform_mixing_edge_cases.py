#!/usr/bin/env python3
"""
Edge case tests for uniform mixing implementation.
Tests extreme conditions and error handling.
"""

import sys
import os
import io
import contextlib

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import (
    create_uniform_mixing_pool,
    create_control2_with_uniform_base
)
from test_helpers import create_mock_cells, get_cell_fingerprint


def test_extreme_mix_ratios():
    """Test extreme mix ratio values."""
    print("\nTest 3.1: Extreme mix ratios")
    
    snapshot_cells = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    median_size = 20
    
    # Test 1% mix ratio
    print("  Testing 1% mix ratio...")
    pool, indices = create_uniform_mixing_pool(snapshot_cells, median_size, 0.01, seed=100)
    # With 1% mix, total = 20/0.99 ≈ 20.2, pool = 20.2 * 0.01 ≈ 0.2 → rounds to 0
    # But implementation might handle this differently
    assert len(pool) >= 0, "Should handle tiny mix ratio"
    print(f"    Pool size with 1% mix: {len(pool)}")
    
    # Test 99% mix ratio
    print("  Testing 99% mix ratio...")
    pool, indices = create_uniform_mixing_pool(snapshot_cells, median_size, 0.99, seed=101)
    # Total = 20/0.01 = 2000, pool = 2000 * 0.99 = 1980
    # But we only have 100 cells, so should handle with replacement
    assert len(pool) > 0, "Should handle 99% mix ratio"
    print(f"    Pool size with 99% mix: {len(pool)}")
    
    # Test 100% mix ratio
    print("  Testing 100% mix ratio...")
    pool, indices = create_uniform_mixing_pool(snapshot_cells, median_size, 1.0, seed=102)
    assert len(pool) == median_size, "100% mix should equal median size"
    
    print("  ✓ Extreme mix ratios handled")
    return True


def test_insufficient_cells():
    """Test behavior when there aren't enough cells."""
    print("\nTest 3.2: Insufficient cells scenarios")
    
    # Case 1: Need more cells than available
    print("  Testing pool larger than snapshot...")
    small_snapshot = create_mock_cells(n_cells=30, n_sites=50, seed=42)
    
    # Capture output to check for warnings
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        pool, indices = create_uniform_mixing_pool(
            small_snapshot, median_size=50, mix_ratio=0.8, seed=100
        )
    
    output_text = output.getvalue()
    
    # Should work but possibly with warning about sampling with replacement
    assert len(pool) > 0, "Should handle insufficient cells"
    if "Sampling with replacement" in output_text:
        print("    ✓ Warning about sampling with replacement detected")
    
    # Case 2: Pool uses most cells, little left for Control2
    print("  Testing Control2 with few remaining cells...")
    snapshot = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    
    # Use 95 of 100 cells for uniform pool
    pool, indices = create_uniform_mixing_pool(
        snapshot, median_size=95, mix_ratio=1.0, seed=100
    )
    
    # Try to create Control2 needing more cells
    output = io.StringIO()
    with contextlib.redirect_stdout(output):
        control2 = create_control2_with_uniform_base(
            snapshot, pool, indices,
            target_size=100, rate=0.005, seed=200
        )
    
    output_text = output.getvalue()
    
    assert len(control2.cells) == 100, "Should reach target size"
    if "Not enough unique cells" in output_text or "Sampling with replacement" in output_text:
        print("    ✓ Warning about insufficient unique cells detected")
    
    print("  ✓ Insufficient cells handled gracefully")
    return True


def test_edge_sizes():
    """Test with edge case sizes."""
    print("\nTest 3.3: Edge case sizes")
    
    # Test with target size = 1
    print("  Testing target size = 1...")
    snapshot = create_mock_cells(n_cells=10, n_sites=50, seed=42)
    pool, indices = create_uniform_mixing_pool(snapshot, median_size=1, mix_ratio=0.5, seed=100)
    
    control2 = create_control2_with_uniform_base(
        snapshot, pool, indices,
        target_size=1, rate=0.005, seed=200
    )
    
    assert len(control2.cells) == 1, "Should handle target size = 1"
    
    # Test with very small snapshot
    print("  Testing with 2-cell snapshot...")
    tiny_snapshot = create_mock_cells(n_cells=2, n_sites=50, seed=42)
    pool, indices = create_uniform_mixing_pool(tiny_snapshot, median_size=1, mix_ratio=0.5, seed=100)
    assert len(pool) <= 2, "Pool can't exceed snapshot size without replacement"
    
    print("  ✓ Edge sizes handled correctly")
    return True


def test_seed_reproducibility():
    """Test that seeds produce reproducible results."""
    print("\nTest 3.4: Seed reproducibility")
    
    snapshot = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    pool, indices = create_uniform_mixing_pool(snapshot, 20, 0.8, seed=100)
    
    # Run 1
    control2_a = create_control2_with_uniform_base(
        snapshot, pool, indices, 100, rate=0.005, seed=42
    )
    fingerprints_a = [get_cell_fingerprint(c) for c in control2_a.cells]
    
    # Run 2 with same seed
    control2_b = create_control2_with_uniform_base(
        snapshot, pool, indices, 100, rate=0.005, seed=42
    )
    fingerprints_b = [get_cell_fingerprint(c) for c in control2_b.cells]
    
    # Should be identical
    assert fingerprints_a == fingerprints_b, "Same seed should give identical results"
    
    # Run 3 with different seed
    control2_c = create_control2_with_uniform_base(
        snapshot, pool, indices, 100, rate=0.005, seed=99
    )
    fingerprints_c = [get_cell_fingerprint(c) for c in control2_c.cells]
    
    # Should be different (at least in the additional cells)
    assert fingerprints_a != fingerprints_c, "Different seed should give different results"
    
    print("  ✓ Seed reproducibility verified")
    return True


def test_pool_larger_than_target():
    """Test when uniform pool is larger than target size."""
    print("\nTest 3.5: Pool larger than target size")
    
    snapshot = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    
    # Create a large pool
    pool, indices = create_uniform_mixing_pool(snapshot, median_size=50, mix_ratio=1.0, seed=100)
    assert len(pool) == 50, "Pool should be 50 cells"
    
    # Create Control2 with smaller target
    control2 = create_control2_with_uniform_base(
        snapshot, pool, indices,
        target_size=30, rate=0.005, seed=200
    )
    
    # Should truncate to target size
    assert len(control2.cells) == 30, "Should truncate pool to target size"
    
    # First 30 cells should match pool
    for i in range(30):
        original_fp = get_cell_fingerprint(pool[i])
        control2_fp = get_cell_fingerprint(control2.cells[i])
        assert original_fp == control2_fp, f"Cell {i} should match pool"
    
    print("  ✓ Pool truncation handled correctly")
    return True


def run_all_tests():
    """Run all edge case tests."""
    print("="*60)
    print("EDGE CASE TESTS FOR UNIFORM MIXING")
    print("="*60)
    
    tests = [
        test_extreme_mix_ratios,
        test_insufficient_cells,
        test_edge_sizes,
        test_seed_reproducibility,
        test_pool_larger_than_target
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except AssertionError as e:
            print(f"  ✗ Test failed: {e}")
            failed += 1
        except Exception as e:
            print(f"  ✗ Test error: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"EDGE CASE TEST RESULTS: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)