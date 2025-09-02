#!/usr/bin/env python3
"""
Test that control2 is saved without cell_history.
"""

import sys
import os
import json
import tempfile
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from core.pipeline_utils import create_pure_snapshot_petri, save_petri_dish

def test_control2_no_history():
    """Test that control2 is saved without history."""
    print("Testing control2 saved without history...")
    
    # Create test cells
    test_cells = []
    for i in range(10):
        cell = Cell(n=20, rate=0.005, gene_size=5)
        cell.age = 50
        test_cells.append(cell)
    
    # Create control2
    petri = create_pure_snapshot_petri(test_cells, n_cells=5, rate=0.005, seed=42)
    
    # Save WITHOUT history (as control2 should be)
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        save_petri_dish(petri, f.name, include_cell_history=False, compress=False)
        
        # Read back and check
        with open(f.name) as rf:
            data = json.load(rf)
        
        print(f"  Top-level keys: {list(data.keys())}")
        print(f"  Has cell_history? {'cell_history' in data}")
        print(f"  Number of cells: {len(data['cells'])}")
        print(f"  Metadata present? {'metadata' in data}")
        if 'metadata' in data:
            print(f"  Metadata keys: {list(data['metadata'].keys())[:5]}...")
        
        assert 'cell_history' not in data, "cell_history should NOT be saved!"
        assert 'cells' in data, "cells should be saved"
        assert 'metadata' in data, "metadata should be saved"
        assert len(data['cells']) == 5, "Should have 5 cells"
        
        os.unlink(f.name)
    
    print("  ✓ Test passed! No cell_history in saved file\n")
    return True

def test_control1_has_history():
    """Test that control1/mutant still have history (for comparison)."""
    print("Testing that grown individuals still have history...")
    
    # This is just to confirm the difference - control1/mutant SHOULD have history
    # because they have growth trajectories we want to plot
    
    # Create a PetriDish with one cell (like mutant/control1 start)
    cell = Cell(n=20, rate=0.005, gene_size=5)
    cell.age = 30
    petri = PetriDish(cells=[cell], rate=0.005, n=20, gene_size=5)
    
    # Save WITH history (as mutant/control1 should be)
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        save_petri_dish(petri, f.name, include_cell_history=True, compress=False)
        
        with open(f.name) as rf:
            data = json.load(rf)
        
        print(f"  Has cell_history? {'cell_history' in data}")
        assert 'cell_history' in data, "Grown individuals should have history"
        
        os.unlink(f.name)
    
    print("  ✓ Grown individuals still have history\n")
    return True

def main():
    print("="*60)
    print("Testing Control2 No History")
    print("="*60)
    print()
    
    try:
        test_control2_no_history()
        test_control1_has_history()
        
        print("="*60)
        print("All tests PASSED! ✓")
        print("Control2: No history (pure snapshot)")
        print("Control1/Mutant: Has history (for trajectory plots)")
        print("="*60)
        return 0
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())