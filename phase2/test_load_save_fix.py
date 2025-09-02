#!/usr/bin/env python3
"""
Test that the load_petri_dish fix correctly preserves gene rate groups.
"""

import sys
import os
import json
import tempfile
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from core.pipeline_utils import save_petri_dish, load_petri_dish

def test_save_load_cycle():
    """Test that gene rate groups survive a save/load cycle."""
    print("Testing save/load cycle with gene rate groups...")
    
    # Create PetriDish with gene rate groups
    gene_rate_groups = [(5, 0.004), (5, 0.005), (5, 0.006), (5, 0.007)]
    cell = Cell(n=100, gene_rate_groups=gene_rate_groups, gene_size=5)
    
    original_petri = PetriDish.from_cells(cell, growth_phase=7)
    original_petri.metadata = {
        'individual_id': 1,
        'individual_type': 'test'
    }
    
    print(f"  Original PetriDish.rate: {original_petri.rate}")
    print(f"  Original PetriDish.gene_rate_groups: {original_petri.gene_rate_groups}")
    
    # Save to temp file
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        save_petri_dish(original_petri, f.name, compress=False)
        temp_path = f.name
    
    # Check what was actually saved
    with open(temp_path, 'r') as f:
        saved_data = json.load(f)
    
    print(f"\n  Saved metadata['rate']: {saved_data['metadata'].get('rate')}")
    print(f"  Saved metadata['gene_rate_groups']: {saved_data['metadata'].get('gene_rate_groups')}")
    
    # Load back
    loaded_petri = load_petri_dish(temp_path)
    
    print(f"\n  Loaded PetriDish.rate: {loaded_petri.rate}")
    print(f"  Loaded PetriDish.gene_rate_groups: {loaded_petri.gene_rate_groups}")
    
    # Clean up
    os.unlink(temp_path)
    
    # Verify
    assert loaded_petri.rate is None, f"Loaded rate should be None, got {loaded_petri.rate}"
    assert loaded_petri.gene_rate_groups == gene_rate_groups, f"Gene rate groups not preserved! Got {loaded_petri.gene_rate_groups}"
    
    print("\n  ✓ Gene rate groups correctly preserved through save/load cycle!")
    return True

def test_with_growth():
    """Test that gene rate groups survive growth and re-saving."""
    print("\nTesting with growth cycle...")
    
    # Create PetriDish with gene rate groups
    gene_rate_groups = [(5, 0.004), (5, 0.005), (5, 0.006), (5, 0.007)]
    cell = Cell(n=100, gene_rate_groups=gene_rate_groups, gene_size=5)
    
    petri = PetriDish.from_cells(cell, growth_phase=2)
    
    # Simulate growth
    petri.divide_cells()
    petri.methylate_cells()
    
    print(f"  After growth - cells: {len(petri.cells)}")
    print(f"  PetriDish.rate: {petri.rate}")
    print(f"  PetriDish.gene_rate_groups: {petri.gene_rate_groups}")
    
    # Save
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        save_petri_dish(petri, f.name, compress=False)
        temp_path = f.name
    
    # Load back
    loaded = load_petri_dish(temp_path)
    
    print(f"\n  After load - cells: {len(loaded.cells)}")
    print(f"  Loaded.rate: {loaded.rate}")
    print(f"  Loaded.gene_rate_groups: {loaded.gene_rate_groups}")
    
    os.unlink(temp_path)
    
    assert loaded.gene_rate_groups == gene_rate_groups, "Gene rate groups lost after growth!"
    print("\n  ✓ Gene rate groups preserved after growth and reload!")
    return True

def main():
    print("="*60)
    print("Testing Load/Save Fix for Gene Rate Groups")
    print("="*60)
    print()
    
    try:
        test_save_load_cycle()
        test_with_growth()
        
        print("\n" + "="*60)
        print("All tests PASSED! ✓")
        print("Fix successfully preserves gene rate groups")
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