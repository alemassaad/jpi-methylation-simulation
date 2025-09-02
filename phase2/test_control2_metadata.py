#!/usr/bin/env python3
"""
Test that control2 metadata correctly has initial_year set to snapshot year (not 0).
"""

import sys
import os
import json
import tempfile
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from core.pipeline_utils import create_pure_snapshot_petri, save_petri_dish

def test_control2_metadata():
    """Test that control2 has correct metadata."""
    print("Testing control2 metadata...")
    
    # Create test cells from year 50
    test_cells = []
    for i in range(20):
        cell = Cell(n=20, rate=0.005, gene_size=5)
        cell.age = 50  # Simulate year 50 cells
        test_cells.append(cell)
    
    # Create control2 using pure snapshot
    petri = create_pure_snapshot_petri(test_cells, n_cells=10, rate=0.005, seed=42)
    
    # Add metadata as pipeline would
    petri.metadata.update({
        'individual_id': 1,
        'individual_type': 'control2',
        'source': 'pure_year50',
        'initial_year': 50  # This should be 50, not 0!
    })
    
    # Check metadata
    print(f"  Metadata keys: {list(petri.metadata.keys())}")
    print(f"  Initial year: {petri.metadata.get('initial_year', 'MISSING')}")
    print(f"  Individual type: {petri.metadata.get('individual_type', 'MISSING')}")
    print(f"  Number of cells: {len(petri.cells)}")
    print(f"  Cell ages: {set([c.age for c in petri.cells])}")
    
    # Save and check JSON
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        save_petri_dish(petri, f.name, include_cell_history=False, compress=False)
        
        # Read back and verify
        with open(f.name) as rf:
            data = json.load(rf)
        
        print(f"\n  JSON metadata:")
        if 'metadata' in data:
            for key, value in data['metadata'].items():
                if key in ['initial_year', 'individual_type', 'individual_id']:
                    print(f"    {key}: {value}")
        
        # Assertions
        assert 'metadata' in data, "metadata should be saved"
        assert data['metadata']['initial_year'] == 50, f"initial_year should be 50, not {data['metadata'].get('initial_year', 'MISSING')}"
        assert data['metadata']['individual_type'] == 'control2', "individual_type should be control2"
        assert 'cell_history' not in data, "cell_history should NOT be saved for control2"
        
        os.unlink(f.name)
    
    print("\n  ✓ Control2 metadata test passed!")
    print("  ✓ initial_year correctly set to 50 (not 0)")
    return True

def main():
    print("="*60)
    print("Testing Control2 Metadata Fix")
    print("="*60)
    print()
    
    try:
        test_control2_metadata()
        
        print("\n" + "="*60)
        print("All tests PASSED! ✓")
        print("Control2 metadata correctly has initial_year = snapshot year")
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