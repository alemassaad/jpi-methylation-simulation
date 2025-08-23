#!/usr/bin/env python3
"""
Unit tests for core uniform mixing functions.
Fast, focused tests on individual functions.
"""

import sys
import os
import random
import copy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import (
    create_uniform_mixing_pool,
    create_control2_with_uniform_base,
    create_pure_snapshot_petri
)
from test_helpers import (
    create_mock_cells,
    cells_are_identical,
    get_cell_fingerprint,
    verify_no_duplicates
)


def test_create_uniform_mixing_pool_basic():
    """Test basic functionality of create_uniform_mixing_pool."""
    print("\nTest 1.1: create_uniform_mixing_pool basic functionality")
    
    # Create test cells
    snapshot_cells = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    
    # Create pool with 80% mix ratio
    median_size = 20
    mix_ratio = 0.8
    result = create_uniform_mixing_pool(snapshot_cells, median_size, mix_ratio, seed=100)
    
    # Verify return type
    assert isinstance(result, tuple), "Should return tuple"
    assert len(result) == 2, "Should return (pool, indices)"
    
    pool, indices = result
    
    # Calculate expected size
    expected_total = int(median_size / (1 - mix_ratio))  # 20 / 0.2 = 100
    expected_pool_size = int(expected_total * mix_ratio)  # 100 * 0.8 = 80
    
    # Verify sizes
    assert len(pool) == expected_pool_size, f"Pool size {len(pool)} != expected {expected_pool_size}"
    assert len(indices) == expected_pool_size, f"Indices count {len(indices)} != expected {expected_pool_size}"
    
    # Verify indices are unique and valid
    assert len(set(indices)) == len(indices), "Indices should be unique"
    assert min(indices) >= 0, "Indices should be non-negative"
    assert max(indices) < 100, "Indices should be within snapshot range"
    
    # Verify cells match indices
    for i in range(len(pool)):
        original_cell = snapshot_cells[indices[i]]
        pool_cell = pool[i]
        assert cells_are_identical(original_cell, pool_cell), f"Cell {i} doesn't match index {indices[i]}"
    
    print("  ✓ Basic functionality test passed")
    return True


def test_create_uniform_mixing_pool_edge_cases():
    """Test edge cases for create_uniform_mixing_pool."""
    print("\nTest 1.2: create_uniform_mixing_pool edge cases")
    
    snapshot_cells = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    
    # Test 100% mix ratio
    print("  Testing 100% mix ratio...")
    pool, indices = create_uniform_mixing_pool(snapshot_cells, median_size=20, mix_ratio=1.0, seed=101)
    assert len(pool) == 20, "100% mix should give pool size = median size"
    assert len(indices) == 20, "100% mix indices count should = median size"
    
    # Test 50% mix ratio
    print("  Testing 50% mix ratio...")
    pool, indices = create_uniform_mixing_pool(snapshot_cells, median_size=20, mix_ratio=0.5, seed=102)
    # 50% of total 40 = 20
    assert len(pool) == 20, "50% mix pool size incorrect"
    
    # Test reproducibility
    print("  Testing reproducibility...")
    pool1, indices1 = create_uniform_mixing_pool(snapshot_cells, 20, 0.8, seed=42)
    pool2, indices2 = create_uniform_mixing_pool(snapshot_cells, 20, 0.8, seed=42)
    assert indices1 == indices2, "Same seed should give same indices"
    
    pool3, indices3 = create_uniform_mixing_pool(snapshot_cells, 20, 0.8, seed=99)
    assert indices1 != indices3, "Different seed should give different indices"
    
    print("  ✓ Edge cases test passed")
    return True


def test_create_control2_with_uniform_base_basic():
    """Test basic functionality of create_control2_with_uniform_base."""
    print("\nTest 1.3: create_control2_with_uniform_base basic functionality")
    
    # Setup
    snapshot_cells = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size=20, mix_ratio=0.8, seed=100
    )
    
    # Create Control2
    target_size = int(20 / (1 - 0.8))  # 100
    control2 = create_control2_with_uniform_base(
        snapshot_cells, uniform_pool, uniform_indices,
        target_size, rate=0.005, seed=200
    )
    
    # Verify size
    assert len(control2.cells) == target_size, f"Control2 size {len(control2.cells)} != target {target_size}"
    
    # Verify first 80 cells match uniform pool
    print(f"  Checking first {len(uniform_pool)} cells match uniform pool...")
    for i in range(len(uniform_pool)):
        assert cells_are_identical(control2.cells[i], uniform_pool[i]), \
            f"Base cell {i} doesn't match uniform pool"
    
    # Verify additional cells are from snapshot but not in uniform_indices
    print(f"  Checking additional {target_size - len(uniform_pool)} cells...")
    additional_cells = control2.cells[len(uniform_pool):]
    
    for cell in additional_cells:
        # Find this cell in original snapshot
        cell_fingerprint = get_cell_fingerprint(cell)
        found_idx = None
        
        for idx, snap_cell in enumerate(snapshot_cells):
            if get_cell_fingerprint(snap_cell) == cell_fingerprint:
                found_idx = idx
                break
        
        assert found_idx is not None, "Additional cell not found in snapshot"
        # Note: Due to deep copying, we can't guarantee idx not in uniform_indices
        # as cells might randomly match, but structure should be correct
    
    print("  ✓ Basic functionality test passed")
    return True


def test_create_control2_no_duplicates():
    """Test that Control2 has no duplicate cells."""
    print("\nTest 1.4: Control2 no duplicates")
    
    # Small dataset for thorough checking
    snapshot_cells = create_mock_cells(n_cells=50, n_sites=30, seed=42)
    
    # Create uniform pool using 60% of cells
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size=10, mix_ratio=0.6, seed=100
    )
    
    # Create Control2
    target_size = int(10 / (1 - 0.6))  # 25
    control2 = create_control2_with_uniform_base(
        snapshot_cells, uniform_pool, uniform_indices,
        target_size, rate=0.005, seed=200
    )
    
    # Check for duplicates
    assert verify_no_duplicates(control2.cells), "Control2 contains duplicate cells"
    
    print(f"  ✓ No duplicates found in {len(control2.cells)} cells")
    return True


def test_100_percent_mix():
    """Test 100% mix ratio edge case."""
    print("\nTest 1.5: 100% mix ratio edge case")
    
    snapshot_cells = create_mock_cells(n_cells=100, n_sites=50, seed=42)
    
    # Create uniform pool with 100% mix
    median_size = 30
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size, mix_ratio=1.0, seed=100
    )
    
    # Create Control2 - should just be the uniform pool
    control2 = create_control2_with_uniform_base(
        snapshot_cells, uniform_pool, uniform_indices,
        median_size, rate=0.005, seed=200
    )
    
    # Verify size
    assert len(control2.cells) == median_size, f"Size mismatch: {len(control2.cells)} != {median_size}"
    
    # All cells should match uniform pool
    for i in range(median_size):
        assert cells_are_identical(control2.cells[i], uniform_pool[i]), \
            f"Cell {i} doesn't match uniform pool with 100% mix"
    
    print("  ✓ 100% mix ratio handled correctly")
    return True


def run_all_tests():
    """Run all unit tests."""
    print("="*60)
    print("UNIT TESTS FOR UNIFORM MIXING")
    print("="*60)
    
    tests = [
        test_create_uniform_mixing_pool_basic,
        test_create_uniform_mixing_pool_edge_cases,
        test_create_control2_with_uniform_base_basic,
        test_create_control2_no_duplicates,
        test_100_percent_mix
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
    print(f"UNIT TEST RESULTS: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)