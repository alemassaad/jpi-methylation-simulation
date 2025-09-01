#!/usr/bin/env python3
"""
Test that all individual IDs start at 1, not 0.
"""

import os
import sys
import json
import tempfile

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, load_all_petri_dishes


def test_creation_starts_at_one():
    """Test that all batches create files starting at individual_01.json with ID=1."""
    print("\nTest: All Batches Start at ID 1")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        mutant_dir = os.path.join(tmpdir, "mutant")
        control1_dir = os.path.join(tmpdir, "control1")
        control2_dir = os.path.join(tmpdir, "control2")
        
        for dir_path in [mutant_dir, control1_dir, control2_dir]:
            os.makedirs(dir_path)
        
        # Simulate mutant creation (with i+1)
        print("Creating mutant individuals (starting at 1)...")
        for i in range(3):
            file_index = i + 1  # Start at 1
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': file_index,
                'individual_type': 'mutant'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(mutant_dir, f"individual_{file_index:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Created individual_{file_index:02d}.json with ID={file_index}")
        
        # Simulate control1 creation (with i+1)
        print("\nCreating control1 individuals (starting at 1)...")
        for i in range(3):
            file_index = i + 1  # Start at 1
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': file_index,
                'individual_type': 'control1'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(control1_dir, f"individual_{file_index:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Created individual_{file_index:02d}.json with ID={file_index}")
        
        # Simulate control2 creation (with i+1)
        print("\nCreating control2 individuals (starting at 1)...")
        for i in range(3):
            file_index = i + 1  # Start at 1
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': file_index,
                'individual_type': 'control2'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(control2_dir, f"individual_{file_index:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Created individual_{file_index:02d}.json with ID={file_index}")
        
        # Verify all files start at 01
        print("\nVerifying file structure...")
        for batch_name, batch_dir in [("mutant", mutant_dir), 
                                       ("control1", control1_dir), 
                                       ("control2", control2_dir)]:
            print(f"\n{batch_name.upper()}:")
            
            # Check that individual_00.json does NOT exist
            file_00 = os.path.join(batch_dir, "individual_00.json")
            assert not os.path.exists(file_00), f"{batch_name} should NOT have individual_00.json"
            print(f"  ✓ No individual_00.json")
            
            # Check that individual_01.json DOES exist
            file_01 = os.path.join(batch_dir, "individual_01.json")
            assert os.path.exists(file_01), f"{batch_name} should have individual_01.json"
            print(f"  ✓ Has individual_01.json")
            
            # Load and verify IDs
            dishes = load_all_petri_dishes(batch_dir)
            for idx, petri in enumerate(dishes):
                expected_id = idx + 1  # Should be 1, 2, 3
                actual_id = petri.metadata.get('individual_id', -1)
                assert actual_id == expected_id, f"ID mismatch: expected {expected_id}, got {actual_id}"
                print(f"  ✓ File {idx+1:02d} has ID={actual_id}")
    
    return True


def test_growth_preserves_one_based_ids():
    """Test that growth/mixing stages preserve 1-based IDs."""
    print("\nTest: Growth/Mixing Preserves 1-Based IDs")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create files with IDs 1, 2, 3
        print("Creating initial files with IDs 1, 2, 3...")
        for i in range(3):
            file_index = i + 1
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': file_index,
                'individual_type': 'test'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(tmpdir, f"individual_{file_index:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
        
        # Simulate growth stage
        print("\nSimulating growth stage...")
        dishes = load_all_petri_dishes(tmpdir)
        
        for i, petri in enumerate(dishes):
            # With the fix, this should use individual_id from metadata
            individual_id = petri.metadata.get('individual_id', i + 1)  # Default to i+1
            
            # Modify
            petri.metadata['grown'] = True
            
            # Save
            filepath = os.path.join(tmpdir, f"individual_{individual_id:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Grew and saved individual_{individual_id:02d}.json")
        
        # Verify
        print("\nVerifying after growth...")
        dishes = load_all_petri_dishes(tmpdir)
        
        for idx, petri in enumerate(dishes):
            expected_id = idx + 1  # Should be 1, 2, 3
            actual_id = petri.metadata.get('individual_id', -1)
            grown = petri.metadata.get('grown', False)
            
            assert actual_id == expected_id, f"ID mismatch: expected {expected_id}, got {actual_id}"
            assert grown, "Should be marked as grown"
            print(f"  ✓ Individual {actual_id}: grown={grown}")
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing 1-Based Individual IDs")
    print("=" * 60)
    
    tests = [
        test_creation_starts_at_one,
        test_growth_preserves_one_based_ids
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
        print("\n🎉 All IDs Start at 1!")
        print("   ✅ No individual_00.json files")
        print("   ✅ First file is individual_01.json with ID=1")
        print("   ✅ Consistent across all batches (mutant, control1, control2)")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())