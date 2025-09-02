#!/usr/bin/env python3
"""
Test script for mean and median gene JSD tracking.
Tests the new dish-level gene JSD summary statistics.
"""

import sys
import os
import json
import gzip
import tempfile
import statistics

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cell import Cell, PetriDish

def test_basic_tracking():
    """Test that mean and median gene JSD are tracked correctly."""
    print("Testing basic mean/median gene JSD tracking...")
    
    # Create a small petri dish
    petri = PetriDish(rate=0.01, n=20, gene_size=5, seed=42, growth_phase=2)
    
    # Enable tracking
    petri.enable_history_tracking(track_gene_jsd=True)
    
    # Verify initial values
    assert hasattr(petri, 'mean_gene_jsd_history'), "mean_gene_jsd_history not initialized"
    assert hasattr(petri, 'median_gene_jsd_history'), "median_gene_jsd_history not initialized"
    assert '0' in petri.mean_gene_jsd_history, "Year 0 not in mean history"
    assert '0' in petri.median_gene_jsd_history, "Year 0 not in median history"
    assert petri.mean_gene_jsd_history['0'] == 0.0, "Initial mean should be 0"
    assert petri.median_gene_jsd_history['0'] == 0.0, "Initial median should be 0"
    
    print("  ✓ Initial histories created correctly")
    
    # Simulate a few years
    for year in range(1, 5):
        petri.simulate_year()
    
    # Check that histories are populated
    assert len(petri.mean_gene_jsd_history) == 5, f"Should have 5 years, got {len(petri.mean_gene_jsd_history)}"
    assert len(petri.median_gene_jsd_history) == 5, f"Should have 5 years, got {len(petri.median_gene_jsd_history)}"
    
    # Verify values increase over time (should generally increase as methylation accumulates)
    mean_values = [petri.mean_gene_jsd_history[str(y)] for y in range(5)]
    median_values = [petri.median_gene_jsd_history[str(y)] for y in range(5)]
    
    print(f"  Mean gene JSD over time: {[f'{v:.4f}' for v in mean_values]}")
    print(f"  Median gene JSD over time: {[f'{v:.4f}' for v in median_values]}")
    
    # Check that values are non-negative and bounded
    for year in range(5):
        mean_val = petri.mean_gene_jsd_history[str(year)]
        median_val = petri.median_gene_jsd_history[str(year)]
        assert 0 <= mean_val <= 1, f"Mean JSD out of bounds at year {year}: {mean_val}"
        assert 0 <= median_val <= 1, f"Median JSD out of bounds at year {year}: {median_val}"
    
    print("  ✓ Values tracked correctly over time")
    
    # Test properties
    current_mean = petri.mean_gene_jsd
    current_median = petri.median_gene_jsd
    
    # These should match the last recorded values
    assert abs(current_mean - mean_values[-1]) < 0.0001, "Property doesn't match history"
    assert abs(current_median - median_values[-1]) < 0.0001, "Property doesn't match history"
    
    print("  ✓ Properties work correctly")
    
    return True

def test_save_and_load():
    """Test that mean/median histories are saved and loaded correctly."""
    print("\nTesting save and load of mean/median gene JSD...")
    
    # Create and simulate
    petri = PetriDish(rate=0.01, n=20, gene_size=5, seed=123, growth_phase=2)
    petri.enable_history_tracking(track_gene_jsd=True)
    
    for year in range(1, 4):
        petri.simulate_year()
    
    # Save original values
    original_mean_history = dict(petri.mean_gene_jsd_history)
    original_median_history = dict(petri.median_gene_jsd_history)
    
    # Save to file
    with tempfile.TemporaryDirectory() as temp_dir:
        saved_path = petri.save_history(directory=temp_dir, compress=True)
        
        # Load the file
        with gzip.open(saved_path, 'rt') as f:
            loaded_data = json.load(f)
    
    # Check that mean/median are in the history
    for year in range(4):
        year_str = str(year)
        assert year_str in loaded_data['history'], f"Year {year} not in history"
        year_data = loaded_data['history'][year_str]
        
        if year_str in original_mean_history:
            assert 'mean_gene_jsd' in year_data, f"mean_gene_jsd missing for year {year}"
            assert abs(year_data['mean_gene_jsd'] - original_mean_history[year_str]) < 0.0001, \
                f"Mean mismatch at year {year}"
        
        if year_str in original_median_history:
            assert 'median_gene_jsd' in year_data, f"median_gene_jsd missing for year {year}"
            assert abs(year_data['median_gene_jsd'] - original_median_history[year_str]) < 0.0001, \
                f"Median mismatch at year {year}"
    
    print("  ✓ Mean/median histories saved correctly")
    
    return True

def test_consistency_with_gene_jsds():
    """Test that mean/median are consistent with actual gene JSDs."""
    print("\nTesting consistency with gene JSDs...")
    
    petri = PetriDish(rate=0.01, n=20, gene_size=5, seed=999, growth_phase=2)
    petri.enable_history_tracking(track_gene_jsd=True)
    
    # Simulate and check consistency at each year
    for year in range(1, 4):
        petri.simulate_year()
        
        # Get the stored values
        stored_mean = petri.mean_gene_jsd_history[str(year)]
        stored_median = petri.median_gene_jsd_history[str(year)]
        
        # Calculate from gene JSDs
        gene_jsds = petri.gene_jsd_history[str(year)]
        calculated_mean = sum(gene_jsds) / len(gene_jsds) if gene_jsds else 0.0
        calculated_median = statistics.median(gene_jsds) if gene_jsds else 0.0
        
        # Compare
        assert abs(stored_mean - calculated_mean) < 0.0001, \
            f"Mean mismatch at year {year}: stored={stored_mean}, calculated={calculated_mean}"
        assert abs(stored_median - calculated_median) < 0.0001, \
            f"Median mismatch at year {year}: stored={stored_median}, calculated={calculated_median}"
    
    print("  ✓ Mean/median values consistent with gene JSDs")
    
    return True

def test_disabled_tracking():
    """Test that tracking can be disabled."""
    print("\nTesting disabled tracking...")
    
    # Create without gene JSD tracking
    petri = PetriDish(rate=0.01, n=20, gene_size=5, seed=42, growth_phase=2, 
                      calculate_cell_jsds=False)
    
    # Simulate
    for year in range(1, 3):
        petri.simulate_year()
    
    # Histories should be empty
    assert len(petri.mean_gene_jsd_history) == 0, "Should have no mean history when disabled"
    assert len(petri.median_gene_jsd_history) == 0, "Should have no median history when disabled"
    
    print("  ✓ Tracking correctly disabled when calculate_cell_jsds=False")
    
    return True

def main():
    """Run all tests."""
    print("="*60)
    print("Testing Mean/Median Gene JSD Implementation")
    print("="*60)
    
    all_passed = True
    
    try:
        test_basic_tracking()
    except AssertionError as e:
        print(f"  ✗ Basic tracking test failed: {e}")
        all_passed = False
    
    try:
        test_save_and_load()
    except AssertionError as e:
        print(f"  ✗ Save/load test failed: {e}")
        all_passed = False
    
    try:
        test_consistency_with_gene_jsds()
    except AssertionError as e:
        print(f"  ✗ Consistency test failed: {e}")
        all_passed = False
    
    try:
        test_disabled_tracking()
    except AssertionError as e:
        print(f"  ✗ Disabled tracking test failed: {e}")
        all_passed = False
    
    print("\n" + "="*60)
    if all_passed:
        print("All tests PASSED! ✓")
        print("="*60)
        return 0
    else:
        print("Some tests FAILED! ✗")
        print("="*60)
        return 1

if __name__ == "__main__":
    sys.exit(main())