#!/usr/bin/env python3
"""
Test the simplified year tracking system.
Verifies that absolute_year is gone and year tracking is simple and clear.
"""

import os
import sys
import json
import tempfile

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, grow_petri_for_years


def test_no_absolute_year():
    """Test that absolute_year is completely removed."""
    print("\nTest 1: No absolute_year")
    print("=" * 50)
    
    # Create PetriDish
    petri = PetriDish(n=100, rate=0.005, growth_phase=5)
    
    # Verify no absolute_year
    assert not hasattr(petri, 'absolute_year'), "absolute_year should not exist!"
    
    # Verify year starts at 0
    assert petri.year == 0, "Year should start at 0"
    
    # Simulate a year
    petri.simulate_year()
    assert petri.year == 1, "Year should be 1 after simulation"
    assert not hasattr(petri, 'absolute_year'), "absolute_year should still not exist!"
    
    print("  ✅ absolute_year is completely removed")
    print(f"  ✅ Simple year tracking: year = {petri.year}")
    
    return True


def test_year_tracking_during_growth():
    """Test that year increments properly during growth."""
    print("\nTest 2: Year Tracking During Growth")
    print("=" * 50)
    
    # Create individual from snapshot
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=7,
        metadata={
            'individual_id': 0,
            'individual_type': 'test',
            'initial_year': 50  # Snapshot from year 50
        }
    )
    
    # Verify starts at 0
    assert petri.year == 0, "PetriDish age should start at 0"
    assert petri.metadata['initial_year'] == 50, "Should track snapshot year in metadata"
    
    print(f"  Created from year {petri.metadata['initial_year']} snapshot")
    print(f"  PetriDish age: {petri.year}")
    
    # Grow for 20 years
    grow_petri_for_years(petri, years=20, growth_phase=7, verbose=False)
    
    # Year should be 20 (age of PetriDish)
    assert petri.year == 20, f"After 20 years growth, year should be 20, got {petri.year}"
    
    print(f"  After 20 years growth: year = {petri.year}")
    print(f"  Simulation year = initial_year + year = {petri.metadata['initial_year']} + {petri.year} = {petri.metadata['initial_year'] + petri.year}")
    
    print("  ✅ Year tracking is simple and clear")
    
    return True


def test_metadata_clarity():
    """Test that metadata is clear about year tracking."""
    print("\nTest 3: Metadata Clarity")
    print("=" * 50)
    
    # Create and grow a PetriDish
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=7,
        metadata={
            'individual_id': 0,
            'individual_type': 'mutant',
            'initial_year': 50,
            'source_quantile': 0.9
        }
    )
    
    # Grow for 20 years
    for _ in range(20):
        petri.divide_cells()
        petri.methylate_cells()
        petri.year += 1  # Manual increment for this test
    
    # Update metadata with growth info
    petri.update_metadata({
        'growth_complete': True,
        'growth_years': 20
    })
    
    # Save and check
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_gene_metrics=True, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        metadata = data['metadata']
        
        print("  Saved metadata:")
        print(f"    year: {metadata['year']} (PetriDish age)")
        print(f"    initial_year: {metadata.get('initial_year', 'N/A')} (snapshot origin)")
        print(f"    growth_years: {metadata.get('growth_years', 'N/A')} (growth duration)")
        
        # Should NOT have absolute_year
        assert 'absolute_year' not in metadata, "Should not have absolute_year in metadata!"
        
        # Should have clear year tracking
        assert metadata['year'] == 20, "Year should be 20"
        assert metadata['initial_year'] == 50, "Initial year should be preserved"
        assert metadata['growth_years'] == 20, "Growth years should be recorded"
        
        # Can calculate simulation year if needed
        sim_year = metadata['initial_year'] + metadata['year']
        print(f"    Calculated simulation year: {sim_year}")
        
        print("\n  ✅ Metadata is clear and simple")
        print("  ✅ No confusing absolute_year")
        print("  ✅ Can derive any needed info from simple fields")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def test_increment_year():
    """Test that increment_year works without absolute_year."""
    print("\nTest 4: increment_year() Method")
    print("=" * 50)
    
    petri = PetriDish(n=100, rate=0.005, growth_phase=5)
    
    # Test increment_year
    assert petri.year == 0, "Should start at 0"
    
    petri.increment_year()
    assert petri.year == 1, "Should be 1 after increment"
    
    petri.increment_year()
    assert petri.year == 2, "Should be 2 after second increment"
    
    # No absolute_year to worry about
    assert not hasattr(petri, 'absolute_year'), "No absolute_year to track"
    
    print(f"  Year increments simply: 0 → 1 → 2")
    print("  ✅ increment_year() is simple without absolute_year")
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Simplified Year Tracking Tests")
    print("=" * 60)
    
    tests = [
        test_no_absolute_year,
        test_year_tracking_during_growth,
        test_metadata_clarity,
        test_increment_year
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
        print("\n🎉 Year Tracking Simplified Successfully!")
        print("   ✅ No more absolute_year")
        print("   ✅ Simple year counter (PetriDish age)")
        print("   ✅ Clear metadata (initial_year for origin)")
        print("   ✅ Easy to understand and maintain")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())