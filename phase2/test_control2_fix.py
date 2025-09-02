#!/usr/bin/env python3
"""
Quick test to verify control2 creation doesn't have bogus history.
"""

import sys
import os
import json
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from core.pipeline_utils import create_pure_snapshot_petri, create_control2_with_uniform_base
import copy

def test_pure_snapshot():
    """Test pure snapshot creation (no uniform mixing)."""
    print("Testing pure snapshot creation...")
    
    # Create some test cells (age 50)
    test_cells = []
    for i in range(10):
        cell = Cell(n=20, rate=0.005, gene_size=5)
        cell.age = 50  # Simulate year 50 cells
        test_cells.append(cell)
    
    # Create control2 with pure snapshot
    petri = create_pure_snapshot_petri(test_cells, n_cells=5, rate=0.005, seed=42)
    
    # Check results
    print(f"  Number of cells: {len(petri.cells)}")
    print(f"  Cell ages: {set([c.age for c in petri.cells])}")
    print(f"  Has cell_history? {hasattr(petri, 'cell_history') and bool(petri.cell_history)}")
    if hasattr(petri, 'cell_history') and petri.cell_history:
        print(f"  History years: {list(petri.cell_history.keys())}")
        if '0' in petri.cell_history:
            print(f"  Year 0 cells: {len(petri.cell_history['0'])}")
            if petri.cell_history['0']:
                print(f"  Year 0 cell age: {petri.cell_history['0'][0].get('age', 'N/A')}")
    
    assert len(petri.cells) == 5, "Should have 5 cells"
    assert all(c.age == 50 for c in petri.cells), "All cells should be age 50"
    
    # Check for bogus history
    if hasattr(petri, 'cell_history') and '0' in petri.cell_history:
        if petri.cell_history['0']:
            first_cell_age = petri.cell_history['0'][0].get('age', -1)
            assert first_cell_age != 0, f"FAILED: Found bogus age-0 cell in history!"
    
    print("  ✓ Pure snapshot test passed!\n")
    return True

def test_uniform_base():
    """Test uniform base creation."""
    print("Testing uniform base creation...")
    
    # Create test cells
    test_cells = []
    for i in range(20):
        cell = Cell(n=20, rate=0.005, gene_size=5)
        cell.age = 50
        test_cells.append(cell)
    
    # Create uniform pool (first 8 cells)
    uniform_pool = [copy.deepcopy(test_cells[i]) for i in range(8)]
    uniform_indices = list(range(8))
    
    # Create control2 with uniform base
    petri = create_control2_with_uniform_base(
        test_cells, uniform_pool, uniform_indices,
        target_size=10, rate=0.005, seed=42
    )
    
    # Check results
    print(f"  Number of cells: {len(petri.cells)}")
    print(f"  Cell ages: {set([c.age for c in petri.cells])}")
    print(f"  Has cell_history? {hasattr(petri, 'cell_history') and bool(petri.cell_history)}")
    if hasattr(petri, 'cell_history') and petri.cell_history:
        print(f"  History years: {list(petri.cell_history.keys())}")
        if '0' in petri.cell_history:
            print(f"  Year 0 cells: {len(petri.cell_history['0'])}")
            if petri.cell_history['0']:
                print(f"  Year 0 cell age: {petri.cell_history['0'][0].get('age', 'N/A')}")
    
    assert len(petri.cells) == 10, "Should have 10 cells"
    assert all(c.age == 50 for c in petri.cells), "All cells should be age 50"
    
    # Check for bogus history
    if hasattr(petri, 'cell_history') and '0' in petri.cell_history:
        if petri.cell_history['0']:
            first_cell_age = petri.cell_history['0'][0].get('age', -1)
            assert first_cell_age != 0, f"FAILED: Found bogus age-0 cell in history!"
    
    print("  ✓ Uniform base test passed!\n")
    return True

def main():
    print("="*60)
    print("Testing Control2 Creation Fix")
    print("="*60)
    print()
    
    try:
        test_pure_snapshot()
        test_uniform_base()
        
        print("="*60)
        print("All tests PASSED! ✓")
        print("Control2 creation no longer creates bogus age-0 cells!")
        print("="*60)
        return 0
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())