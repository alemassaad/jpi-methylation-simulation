#!/usr/bin/env python3
"""
Test that verifies the num_cells bug is completely fixed with the professional implementation.
This simulates the exact scenario that exposed the bug originally.
"""

import os
import sys
import json
import tempfile

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, load_petri_dish


def test_num_cells_through_pipeline():
    """Test that num_cells stays accurate through the entire pipeline."""
    print("\nTest: num_cells Through Complete Pipeline")
    print("=" * 50)
    
    # Stage 1: Create individual from snapshot (simulating Stage 3)
    print("\n  Stage 1: Create individual from snapshot cell")
    
    # Create a snapshot cell with some methylation
    snapshot_cell = Cell(n=100, rate=0.005)
    for _ in range(20):
        snapshot_cell.methylate()
    
    # Create individual using professional method
    individual = PetriDish.from_snapshot_cell(
        cell=snapshot_cell,
        growth_phase=7,
        metadata={
            'individual_id': 0,
            'individual_type': 'control1',
            'source': 'uniform',
            'initial_year': 50
        }
    )
    
    assert individual.metadata['num_cells'] == 1, "Should start with 1 cell"
    print(f"    Created individual with {individual.metadata['num_cells']} cell")
    
    # Save and verify
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(individual, temp_file, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        assert data['metadata']['num_cells'] == 1, "Saved file should show 1 cell"
        assert len(data['cells']) == 1, "Should have 1 cell in array"
        print("    ✅ Saved correctly with num_cells = 1")
        
        # Stage 2: Growth (simulating Stage 4)
        print("\n  Stage 2: Grow individual")
        
        # Reload to simulate pipeline
        individual = load_petri_dish(temp_file)
        
        # Grow for 7 divisions
        for i in range(7):
            individual.divide_cells()
            individual.methylate_cells()
        
        expected_cells = 2 ** 7  # 128 cells
        assert len(individual.cells) == expected_cells, f"Should have {expected_cells} cells"
        
        # Update metadata and save
        individual.update_metadata({'growth_complete': True})
        assert individual.metadata['num_cells'] == expected_cells, f"Metadata should show {expected_cells}"
        
        save_petri_dish(individual, temp_file, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        assert data['metadata']['num_cells'] == expected_cells, f"Saved file should show {expected_cells} cells"
        assert len(data['cells']) == expected_cells, f"Should have {expected_cells} cells in array"
        print(f"    ✅ Grown and saved correctly with num_cells = {expected_cells}")
        
        # Stage 3: Mixing (simulating Stage 6)
        print("\n  Stage 3: Mix with additional cells")
        
        # Reload to simulate pipeline
        individual = load_petri_dish(temp_file)
        
        # Simulate mixing by adding cells
        mix_ratio = 0.8  # 80% from snapshot
        target_size = int(len(individual.cells) / (1 - mix_ratio))
        cells_to_add = target_size - len(individual.cells)
        
        # Add cells from "second snapshot"
        for _ in range(cells_to_add):
            mix_cell = Cell(n=100, rate=0.005)
            for _ in range(30):  # More methylation for year 60
                mix_cell.methylate()
            individual.cells.append(mix_cell)
        
        final_cells = len(individual.cells)
        
        # Update metadata
        individual.update_metadata({
            'mixed': True,
            'mix_ratio': 80,
            'final_cells': final_cells
        })
        
        assert individual.metadata['num_cells'] == final_cells, f"Metadata should show {final_cells}"
        
        # Save with gene metrics
        save_petri_dish(individual, temp_file, include_gene_metrics=True, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        assert data['metadata']['num_cells'] == final_cells, f"Saved file should show {final_cells} cells"
        assert len(data['cells']) == final_cells, f"Should have {final_cells} cells in array"
        assert 'gene_jsds' in data['metadata'], "Should have gene metrics"
        print(f"    ✅ Mixed and saved correctly with num_cells = {final_cells}")
        
        # Final verification
        print("\n  Final Verification:")
        print(f"    Metadata num_cells: {data['metadata']['num_cells']}")
        print(f"    Actual cell count: {len(data['cells'])}")
        print(f"    Individual type: {data['metadata']['individual_type']}")
        print(f"    Growth complete: {data['metadata'].get('growth_complete', False)}")
        print(f"    Mixed: {data['metadata'].get('mixed', False)}")
        
        # Verify the original bug is fixed
        assert data['metadata']['num_cells'] == len(data['cells']), \
            "CRITICAL: num_cells doesn't match actual cell count!"
        
        print("\n  ✅ num_cells bug is FIXED - metadata always matches actual count!")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def test_metadata_preservation_with_wrong_values():
    """Test that wrong num_cells values get corrected automatically."""
    print("\nTest: Auto-correction of Wrong num_cells")
    print("=" * 50)
    
    # Create PetriDish with multiple cells
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(cell, growth_phase=5)
    
    # Grow to get multiple cells
    for _ in range(5):
        petri.divide_cells()
    
    actual_cells = len(petri.cells)  # Should be 32
    print(f"  Actual cells after growth: {actual_cells}")
    
    # Intentionally set wrong num_cells (simulating the bug)
    petri.metadata['num_cells'] = 1  # Wrong!
    petri.metadata['some_other_field'] = 'preserved'
    
    print(f"  Intentionally set wrong num_cells: {petri.metadata['num_cells']}")
    
    # Any update should fix it
    petri.update_metadata({'test': 'value'})
    
    assert petri.metadata['num_cells'] == actual_cells, "Should auto-correct"
    assert petri.metadata['some_other_field'] == 'preserved', "Other fields preserved"
    print(f"  After update_metadata(): num_cells = {petri.metadata['num_cells']}")
    
    # Test with save
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        # Set wrong value again
        petri.metadata['num_cells'] = 999
        
        # Save should fix it
        save_petri_dish(petri, temp_file, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        assert data['metadata']['num_cells'] == actual_cells, "Save should correct num_cells"
        print(f"  After save: num_cells = {data['metadata']['num_cells']}")
        
        print("\n  ✅ Auto-correction works - wrong values get fixed!")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing num_cells Bug Fix")
    print("=" * 60)
    
    tests = [
        test_num_cells_through_pipeline,
        test_metadata_preservation_with_wrong_values
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
        print("🎉 The num_cells bug is completely FIXED!")
        print("   Metadata always reflects the actual cell count.")
        return 0
    else:
        print(f"⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())