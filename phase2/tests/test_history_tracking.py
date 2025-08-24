#!/usr/bin/env python3
"""
Test that history tracking works correctly in the pipeline.
"""

import os
import sys
import tempfile
import gzip
import json

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, load_petri_dish, grow_petri_for_years

def test_history_tracking():
    """Test that history is tracked and saved correctly."""
    print("Testing history tracking...")
    
    # Create a small PetriDish
    petri = PetriDish(rate=0.005, n=100, growth_phase=2)
    
    # Grow with history tracking
    print("  Growing PetriDish with history tracking...")
    grow_petri_for_years(
        petri, 
        years=3,
        growth_phase=2,
        verbose=False,
        track_history=True,
        start_year=0
    )
    
    # Check that history was created
    assert hasattr(petri, 'cell_history'), "PetriDish should have cell_history attribute"
    assert petri.cell_history is not None, "cell_history should not be None"
    assert len(petri.cell_history) > 0, "cell_history should have entries"
    print(f"  ✓ History created with {len(petri.cell_history)} time points")
    
    # Save with history
    with tempfile.NamedTemporaryFile(suffix='.json.gz', delete=False) as tmp:
        tmp_path = tmp.name
    
    print("  Saving PetriDish with history...")
    save_petri_dish(petri, tmp_path, include_cell_history=True)
    
    # Load back with history
    print("  Loading PetriDish with history...")
    loaded_petri = load_petri_dish(tmp_path, include_cell_history=True)
    
    # Verify history was loaded
    assert hasattr(loaded_petri, 'cell_history'), "Loaded PetriDish should have cell_history"
    assert loaded_petri.cell_history is not None, "Loaded cell_history should not be None"
    assert len(loaded_petri.cell_history) == len(petri.cell_history), "History should be preserved"
    print(f"  ✓ History loaded with {len(loaded_petri.cell_history)} time points")
    
    # Clean up
    os.remove(tmp_path)
    
    return True


def test_history_plotting():
    """Test that plotting can work with history."""
    print("\nTesting history plotting compatibility...")
    
    # Create a PetriDish with history
    petri = PetriDish(rate=0.005, n=100, growth_phase=2)
    
    # Enable history and grow
    petri.enable_history_tracking(start_year=0)
    petri.grow_with_homeostasis(years=3, growth_phase=2, verbose=False, record_history=True)
    
    # Check that PetriDishPlotter can work with it
    from cell import PetriDishPlotter
    
    try:
        plotter = PetriDishPlotter(petri)
        print("  ✓ PetriDishPlotter initialized successfully")
        
        # Check that plotter can access history
        assert hasattr(plotter.petri, 'cell_history'), "Plotter should see cell_history"
        print(f"  ✓ Plotter can access {len(plotter.petri.cell_history)} history points")
        
        return True
    except ValueError as e:
        print(f"  ✗ Failed to create plotter: {e}")
        return False


def test_gene_rate_history():
    """Test history tracking with gene-specific rates."""
    print("\nTesting history with gene-specific rates...")
    
    # Create PetriDish with gene rates
    gene_groups = [(10, 0.004), (10, 0.006)]
    petri = PetriDish(
        gene_rate_groups=gene_groups,
        n=100,
        gene_size=5,
        growth_phase=2
    )
    
    # Grow with history
    grow_petri_for_years(
        petri,
        years=2,
        growth_phase=1,
        verbose=False,
        track_history=True,
        start_year=0
    )
    
    # Verify history exists
    assert hasattr(petri, 'cell_history'), "Gene-rate PetriDish should have history"
    assert len(petri.cell_history) > 0, "Should have history entries"
    print(f"  ✓ Gene-rate history created with {len(petri.cell_history)} points")
    
    # Save and load
    with tempfile.NamedTemporaryFile(suffix='.json.gz', delete=False) as tmp:
        tmp_path = tmp.name
    
    save_petri_dish(petri, tmp_path, include_cell_history=True)
    loaded = load_petri_dish(tmp_path, include_cell_history=True)
    
    assert hasattr(loaded, 'cell_history'), "Loaded gene-rate PetriDish should have history"
    print(f"  ✓ Gene-rate history preserved after save/load")
    
    # Clean up
    os.remove(tmp_path)
    
    return True


def run_all_tests():
    """Run all history tracking tests."""
    print("="*60)
    print("Testing History Tracking Fix")
    print("="*60)
    
    tests = [
        test_history_tracking,
        test_history_plotting,
        test_gene_rate_history
    ]
    
    failed = 0
    for test_func in tests:
        try:
            if not test_func():
                failed += 1
        except Exception as e:
            print(f"\n✗ {test_func.__name__} FAILED:")
            print(f"  {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    if failed == 0:
        print("✅ All history tracking tests passed!")
        print("\nHistory tracking is now working correctly:")
        print("  • History is created during growth when track_history=True")
        print("  • History is saved as 'cell_history' attribute")
        print("  • History is loaded correctly from files")
        print("  • PetriDishPlotter can access the history")
        print("  • Works with both uniform and gene-specific rates")
    else:
        print(f"❌ {failed} test(s) failed")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(run_all_tests())