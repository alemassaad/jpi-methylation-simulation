#!/usr/bin/env python3
"""
Test cell age tracking system.
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


def test_cell_age_basics():
    """Test basic cell age functionality."""
    print("\nTest 1: Basic Cell Age")
    print("=" * 50)
    
    # Create new cell
    cell = Cell(n=100, rate=0.005)
    assert cell.age == 0, "New cell should have age 0"
    print(f"  ✅ New cell age: {cell.age}")
    
    # Methylate increments age
    cell.methylate()
    assert cell.age == 1, "Age should be 1 after methylation"
    
    cell.methylate()
    assert cell.age == 2, "Age should be 2 after second methylation"
    print(f"  ✅ Age after 2 methylations: {cell.age}")
    
    # Daughter cell inherits age
    daughter = cell.create_daughter_cell()
    assert daughter.age == 2, "Daughter should inherit parent's age"
    print(f"  ✅ Daughter inherits age: {daughter.age}")
    
    return True


def test_age_serialization():
    """Test that age is saved and loaded correctly."""
    print("\nTest 2: Age Serialization")
    print("=" * 50)
    
    # Create cell with some age
    cell = Cell(n=100, rate=0.005)
    for _ in range(50):
        cell.methylate()
    
    assert cell.age == 50, "Cell should be 50 years old"
    print(f"  Cell age before save: {cell.age}")
    
    # Test to_dict
    cell_dict = cell.to_dict()
    assert 'age' in cell_dict, "to_dict should include age"
    assert cell_dict['age'] == 50, "Saved age should be 50"
    print(f"  ✅ to_dict includes age: {cell_dict['age']}")
    
    # Test from_dict
    restored = Cell.from_dict(cell_dict, rate=0.005)
    assert restored.age == 50, "Restored cell should have same age"
    print(f"  ✅ from_dict restores age: {restored.age}")
    
    return True


def test_petri_dish_with_aged_cells():
    """Test that PetriDish preserves cell ages."""
    print("\nTest 3: PetriDish with Aged Cells")
    print("=" * 50)
    
    # Create aged cell (simulating year 50 snapshot)
    snapshot_cell = Cell(n=100, rate=0.005)
    for _ in range(50):
        snapshot_cell.methylate()
    
    print(f"  Snapshot cell age: {snapshot_cell.age}")
    
    # Create PetriDish from snapshot
    petri = PetriDish.from_snapshot_cell(
        cell=snapshot_cell,
        growth_phase=7,
        metadata={'initial_year': 50}
    )
    
    # Check initial state
    assert petri.year == 0, "PetriDish year should be 0"
    assert petri.cells[0].age == 50, "Cell should maintain age 50"
    print(f"  ✅ PetriDish year: {petri.year}, Cell age: {petri.cells[0].age}")
    
    # Grow for 20 years
    for _ in range(20):
        for cell in petri.cells:
            cell.methylate()
        if petri.year < 7:
            petri.divide_cells()
        petri.year += 1
    
    # Check ages after growth
    print(f"  After 20 years growth:")
    print(f"    PetriDish year: {petri.year}")
    print(f"    Cell ages: {petri.cells[0].age} (should be 70)")
    
    assert petri.year == 20, "PetriDish should be 20 years old"
    assert all(cell.age == 70 for cell in petri.cells), "All cells should be 70 years old"
    
    print("  ✅ Cell ages tracked correctly through growth")
    
    return True


def test_save_load_with_ages():
    """Test that ages are preserved through save/load cycle."""
    print("\nTest 4: Save/Load with Ages")
    print("=" * 50)
    
    # Create PetriDish with aged cells
    cell = Cell(n=100, rate=0.005)
    for _ in range(30):
        cell.methylate()
    
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=5,
        metadata={'initial_year': 30, 'test': 'age_tracking'}
    )
    
    # Grow a bit (5 years)
    for year in range(5):
        petri.divide_cells()
        petri.methylate_cells()  # This increments cell age
        petri.year += 1
    
    # Save
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, compress=False)
        
        # Check saved data
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        # Verify cells have age field
        assert all('age' in cell for cell in data['cells']), "All cells should have age"
        ages = [cell['age'] for cell in data['cells']]
        assert all(age == 35 for age in ages), f"All cells should be 35 years old, got {set(ages)}"
        
        print(f"  ✅ Saved {len(data['cells'])} cells with age {ages[0]}")
        
        # Load and verify
        loaded_petri = load_petri_dish(temp_file)
        
        assert all(cell.age == 35 for cell in loaded_petri.cells), "Loaded cells should maintain age"
        print(f"  ✅ Loaded cells maintain age: {loaded_petri.cells[0].age}")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def test_mixed_age_populations():
    """Test mixing cells of different ages."""
    print("\nTest 5: Mixed Age Populations")
    print("=" * 50)
    
    # Create PetriDish with young cells
    young_cell = Cell(n=100, rate=0.005)
    for _ in range(20):
        young_cell.methylate()
    
    petri = PetriDish.from_snapshot_cell(young_cell, growth_phase=5)
    print(f"  Initial cell age: {petri.cells[0].age}")
    
    # Add older cells (simulating mixing)
    old_cells = []
    for _ in range(3):
        old_cell = Cell(n=100, rate=0.005)
        for _ in range(60):
            old_cell.methylate()
        old_cells.append(old_cell)
    
    petri.cells.extend(old_cells)
    petri.update_metadata({'mixed': True})
    
    # Check age distribution
    ages = [cell.age for cell in petri.cells]
    young_count = sum(1 for age in ages if age == 20)
    old_count = sum(1 for age in ages if age == 60)
    
    print(f"  After mixing:")
    print(f"    Young cells (age 20): {young_count}")
    print(f"    Old cells (age 60): {old_count}")
    print(f"    All ages: {sorted(set(ages))}")
    
    assert 20 in ages and 60 in ages, "Should have mixed ages"
    print("  ✅ Mixed population maintains different cell ages")
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Cell Age Tracking Tests")
    print("=" * 60)
    
    tests = [
        test_cell_age_basics,
        test_age_serialization,
        test_petri_dish_with_aged_cells,
        test_save_load_with_ages,
        test_mixed_age_populations
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
        print("\n🎉 Cell Age Tracking Complete!")
        print("   ✅ Cells track their biological age")
        print("   ✅ Age increments with methylation")
        print("   ✅ Daughters inherit parent age")
        print("   ✅ Age preserved through save/load")
        print("   ✅ Mixed populations maintain individual ages")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())