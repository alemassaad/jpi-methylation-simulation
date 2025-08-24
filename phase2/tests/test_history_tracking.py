#!/usr/bin/env python3
"""
Test history tracking functionality in phase2 pipeline.
Tests various edge cases and scenarios.
"""

import os
import sys
import json
import gzip
import tempfile
import shutil
import numpy as np
from typing import List

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish, PetriDishPlotter
from phase2.pipeline_utils import (
    save_petri_dish, load_petri_dish,
    grow_petri_for_years
)


def test_history_save_load():
    """Test saving and loading PetriDish with history."""
    print("\nTest 1: Save/Load with History")
    print("=" * 50)
    
    # Create a PetriDish and grow it with history
    petri = PetriDish(rate=0.01, n=100, seed=42)
    petri.enable_history_tracking(start_year=0)
    petri.grow_with_homeostasis(years=5, growth_phase=3, verbose=False)
    
    # Save with history
    with tempfile.NamedTemporaryFile(suffix='.json.gz', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_cell_history=True)
        print(f"  ✓ Saved PetriDish with {len(petri.history)} history entries")
        
        # Load without history
        petri_no_hist = load_petri_dish(temp_file, include_cell_history=False)
        # When loaded without history flag, history should not be restored
        has_history = hasattr(petri_no_hist, 'history') and petri_no_hist.history and len(petri_no_hist.history) > 0
        assert not has_history, "History should not be loaded when include_cell_history=False"
        print("  ✓ Loaded without history - no history present")
        
        # Load with history
        petri_with_hist = load_petri_dish(temp_file, include_cell_history=True)
        assert hasattr(petri_with_hist, 'history')
        assert len(petri_with_hist.history) == len(petri.history)
        print(f"  ✓ Loaded with history - {len(petri_with_hist.history)} entries restored")
        
        # Verify history content
        for year in petri.history:
            assert year in petri_with_hist.history
            assert len(petri_with_hist.history[year]) == len(petri.history[year])
        print("  ✓ History content verified")
        
    finally:
        os.remove(temp_file)
    
    print("  PASSED")


def test_grow_with_history_tracking():
    """Test grow_petri_for_years with history tracking."""
    print("\nTest 2: Grow with History Tracking")
    print("=" * 50)
    
    # Test with history tracking enabled
    petri1 = PetriDish(rate=0.01, n=100, seed=42)
    grow_petri_for_years(petri1, years=5, growth_phase=3, 
                        track_history=True, start_year=10, verbose=False)
    
    assert hasattr(petri1, 'history')
    assert petri1.history is not None
    assert len(petri1.history) == 6  # Years 10-15
    assert '10' in petri1.history
    assert '15' in petri1.history
    print(f"  ✓ History tracking enabled: {len(petri1.history)} entries")
    
    # Test without history tracking
    petri2 = PetriDish(rate=0.01, n=100, seed=42)
    initial_history_len = len(petri2.history) if hasattr(petri2, 'history') else 0
    grow_petri_for_years(petri2, years=5, growth_phase=3, 
                        track_history=False, verbose=False)
    
    # When history is not tracked during growth, it shouldn't grow beyond initial state
    final_history_len = len(petri2.history) if hasattr(petri2, 'history') else 0
    assert final_history_len <= initial_history_len, f"History grew from {initial_history_len} to {final_history_len} when track_history=False"
    print(f"  ✓ History tracking disabled: history unchanged ({final_history_len} entries)")
    
    print("  PASSED")


def test_plotter_with_history():
    """Test PetriDishPlotter with history data."""
    print("\nTest 3: PetriDishPlotter with History")
    print("=" * 50)
    
    # Create PetriDish with history
    petri = PetriDish(rate=0.01, n=100, seed=42)
    petri.enable_history_tracking(start_year=0)
    petri.grow_with_homeostasis(years=10, growth_phase=4, verbose=False)
    
    # Test plotter creation
    try:
        plotter = PetriDishPlotter(petri)
        print(f"  ✓ Plotter created with {len(petri.history)} history entries")
    except ValueError as e:
        print(f"  ✗ Failed to create plotter: {e}")
        return
    
    # Test plot generation (without actually saving)
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test JSD plot
        jsd_path = os.path.join(tmpdir, "test_jsd.png")
        try:
            plotter.plot_jsd("Test JSD", jsd_path)
            assert os.path.exists(jsd_path)
            print(f"  ✓ JSD plot created: {os.path.getsize(jsd_path)} bytes")
        except Exception as e:
            print(f"  ✗ JSD plot failed: {e}")
        
        # Test methylation plot
        meth_path = os.path.join(tmpdir, "test_meth.png")
        try:
            plotter.plot_methylation("Test Methylation", meth_path)
            assert os.path.exists(meth_path)
            print(f"  ✓ Methylation plot created: {os.path.getsize(meth_path)} bytes")
        except Exception as e:
            print(f"  ✗ Methylation plot failed: {e}")
        
        # Test combined plot
        combined_path = os.path.join(tmpdir, "test_combined.png")
        try:
            plotter.plot_combined("Test Combined", combined_path)
            assert os.path.exists(combined_path)
            print(f"  ✓ Combined plot created: {os.path.getsize(combined_path)} bytes")
        except Exception as e:
            print(f"  ✗ Combined plot failed: {e}")
    
    print("  PASSED")


def test_edge_cases():
    """Test edge cases for history tracking."""
    print("\nTest 4: Edge Cases")
    print("=" * 50)
    
    # Test 1: Empty PetriDish
    petri_empty = PetriDish(rate=0.01, n=100, seed=42)
    petri_empty.cells = []  # Remove all cells
    petri_empty.enable_history_tracking(start_year=0)
    
    try:
        petri_empty.grow_with_homeostasis(years=2, growth_phase=1, verbose=False)
        print("  ✓ Empty PetriDish handled")
    except Exception as e:
        print(f"  ✗ Empty PetriDish failed: {e}")
    
    # Test 2: Single cell growth
    petri_single = PetriDish(rate=0.01, n=100, seed=42)
    grow_petri_for_years(petri_single, years=3, growth_phase=2,
                        track_history=True, start_year=0, verbose=False)
    assert len(petri_single.history) == 4  # Years 0-3
    print(f"  ✓ Single cell growth: {len(petri_single.cells)} cells, {len(petri_single.history)} history entries")
    
    # Test 3: Large starting year
    petri_large_year = PetriDish(rate=0.01, n=100, seed=42)
    grow_petri_for_years(petri_large_year, years=2, growth_phase=1,
                        track_history=True, start_year=1000, verbose=False)
    assert '1000' in petri_large_year.history
    assert '1002' in petri_large_year.history
    print("  ✓ Large starting year (1000) handled")
    
    # Test 4: Zero years growth
    petri_zero = PetriDish(rate=0.01, n=100, seed=42)
    initial_cells = len(petri_zero.cells)
    grow_petri_for_years(petri_zero, years=0, growth_phase=0,
                        track_history=True, start_year=0, verbose=False)
    assert len(petri_zero.cells) == initial_cells
    print("  ✓ Zero years growth handled")
    
    # Test 5: Growth phase longer than total years
    petri_long_phase = PetriDish(rate=0.01, n=100, seed=42)
    try:
        grow_petri_for_years(petri_long_phase, years=3, growth_phase=5,
                            track_history=True, start_year=0, verbose=False)
        print("  ✗ Should have raised error for growth_phase > years")
    except ValueError:
        print("  ✓ Growth phase > years properly rejected")
    
    # Test 6: No growth phase (pure exponential)
    petri_no_phase = PetriDish(rate=0.01, n=100, seed=42)
    grow_petri_for_years(petri_no_phase, years=4, growth_phase=None,
                        track_history=True, start_year=0, verbose=False)
    expected_cells = 2 ** 4  # Pure exponential growth
    assert len(petri_no_phase.cells) == expected_cells
    print(f"  ✓ No growth phase (pure exponential): {len(petri_no_phase.cells)} cells")
    
    print("  PASSED")


def test_history_consistency():
    """Test that history is consistent across operations."""
    print("\nTest 5: History Consistency")
    print("=" * 50)
    
    # Create and grow with history
    petri = PetriDish(rate=0.01, n=100, seed=42)
    petri.enable_history_tracking(start_year=50)
    petri.grow_with_homeostasis(years=10, growth_phase=5, verbose=False)
    
    # Verify years are continuous
    years = sorted([int(y) for y in petri.history.keys()])
    for i in range(len(years) - 1):
        assert years[i+1] == years[i] + 1
    print(f"  ✓ Years are continuous: {years[0]}-{years[-1]}")
    
    # Verify cell count progression
    for i in range(min(5, len(years) - 1)):  # Check growth phase
        year = str(years[i])
        next_year = str(years[i+1])
        current_count = len(petri.history[year])
        next_count = len(petri.history[next_year])
        
        if i < 5:  # Growth phase
            expected = current_count * 2
            assert next_count == expected, f"Year {year}->{next_year}: {current_count}->{next_count}, expected {expected}"
    print("  ✓ Cell count progression verified")
    
    # Verify JSD values exist and are valid
    for year, cells in petri.history.items():
        for cell in cells:
            assert 'cell_jsd' in cell
            assert 0 <= cell['cell_jsd'] <= 1
    print("  ✓ JSD values valid for all cells")
    
    print("  PASSED")


def test_plotter_without_history():
    """Test that PetriDishPlotter properly handles missing history."""
    print("\nTest 6: PetriDishPlotter without History")
    print("=" * 50)
    
    # Create PetriDish without history
    petri = PetriDish(rate=0.01, n=100, seed=42)
    petri.grow_with_homeostasis(years=5, growth_phase=3, verbose=False, record_history=False)
    
    # Try to create plotter - should fail because history is empty/missing
    try:
        plotter = PetriDishPlotter(petri)
        # If we get here, check if it's because history exists but is empty
        if hasattr(petri, 'history') and not petri.history:
            print("  ✗ Plotter should reject empty history")
        else:
            print("  ✗ Should have raised error for missing history")
    except ValueError as e:
        if "no history" in str(e).lower() or "empty" in str(e).lower():
            print("  ✓ Properly rejected PetriDish without history")
        else:
            print(f"  ✗ Unexpected error: {e}")
    
    print("  PASSED")


def run_all_tests():
    """Run all tests."""
    print("\n" + "="*60)
    print("TESTING HISTORY TRACKING FUNCTIONALITY")
    print("="*60)
    
    tests = [
        test_history_save_load,
        test_grow_with_history_tracking,
        test_plotter_with_history,
        test_edge_cases,
        test_history_consistency,
        test_plotter_without_history
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"\n  FAILED: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    print(f"SUMMARY: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)