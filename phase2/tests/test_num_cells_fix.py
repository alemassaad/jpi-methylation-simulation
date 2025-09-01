#!/usr/bin/env python3
"""
Test that num_cells is correctly saved even when petri.metadata contains old values.
"""

import os
import sys
import json
import tempfile

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish

def test_num_cells_preservation():
    """Test that num_cells reflects actual cell count, not metadata."""
    print("Testing num_cells preservation...")
    
    # Create PetriDish with multiple cells
    petri = PetriDish(rate=0.005, growth_phase=3)
    
    # Grow to get multiple cells
    for _ in range(4):
        petri.divide_cells()
        petri.methylate_cells()
    
    actual_cells = len(petri.cells)
    print(f"  Actual cells after growth: {actual_cells}")
    
    # Simulate the bug: add metadata with wrong num_cells
    petri.metadata = {
        'individual_id': 0,
        'individual_type': 'test',
        'num_cells': 1,  # Wrong value!
        'source': 'test'
    }
    
    # Save with the fix
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_gene_metrics=True, compress=False)
        
        # Load and check
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        saved_num_cells = data['metadata']['num_cells']
        actual_cells_in_file = len(data['cells'])
        
        print(f"  Saved num_cells in metadata: {saved_num_cells}")
        print(f"  Actual cells in array: {actual_cells_in_file}")
        
        # Verify they match
        assert saved_num_cells == actual_cells, f"Metadata num_cells ({saved_num_cells}) doesn't match actual ({actual_cells})"
        assert saved_num_cells == actual_cells_in_file, f"Metadata ({saved_num_cells}) doesn't match array ({actual_cells_in_file})"
        
        # Verify other metadata was preserved
        assert data['metadata']['individual_id'] == 0, "individual_id not preserved"
        assert data['metadata']['individual_type'] == 'test', "individual_type not preserved"
        
        print("  ✅ num_cells correctly reflects actual cell count!")
        print("  ✅ Other metadata preserved correctly!")
        
        return True
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)

def main():
    """Run the test."""
    print("="*60)
    print("Testing num_cells Bug Fix")
    print("="*60)
    
    try:
        if test_num_cells_preservation():
            print("\n🎉 Test passed! The fix works correctly.")
            return 0
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())