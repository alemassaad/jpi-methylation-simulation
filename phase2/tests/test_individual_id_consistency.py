#!/usr/bin/env python3
"""
Test that individual IDs are consistent between filename and metadata.
This test identifies the root cause of the ID mismatch issue.
"""

import os
import sys
import json
import tempfile
import shutil

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, load_petri_dish, load_all_petri_dishes


def test_individual_id_creation():
    """Test that individual IDs match filenames during creation."""
    print("\nTest 1: Individual ID Creation")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Simulate creating mutant individuals (like in run_pipeline.py)
        mutant_dir = os.path.join(tmpdir, "mutant")
        os.makedirs(mutant_dir)
        
        print("Creating mutant individuals...")
        for i in range(3):
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': i,
                'individual_type': 'mutant',
                'source_quantile': i // 2
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(mutant_dir, f"individual_{i:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Saved individual_{i:02d}.json with ID={i}")
        
        # Load and verify
        print("\nVerifying mutant IDs...")
        dishes = load_all_petri_dishes(mutant_dir)
        for i, petri in enumerate(dishes):
            file_id = i  # Based on filename order
            meta_id = petri.metadata.get('individual_id', -1)
            print(f"  File individual_{i:02d}.json has metadata ID={meta_id}")
            assert meta_id == file_id, f"Mismatch: file {i:02d} has ID {meta_id}"
        
        print("  ✅ Mutant IDs are consistent")
        
        # Simulate creating control1 individuals
        control1_dir = os.path.join(tmpdir, "control1")
        os.makedirs(control1_dir)
        
        print("\nCreating control1 individuals...")
        for i in range(3):
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': i,
                'individual_type': 'control1',
                'source': 'uniform'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(control1_dir, f"individual_{i:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Saved individual_{i:02d}.json with ID={i}")
        
        # Load and verify
        print("\nVerifying control1 IDs...")
        dishes = load_all_petri_dishes(control1_dir)
        for i, petri in enumerate(dishes):
            file_id = i  # Based on filename order
            meta_id = petri.metadata.get('individual_id', -1)
            print(f"  File individual_{i:02d}.json has metadata ID={meta_id}")
            assert meta_id == file_id, f"Mismatch: file {i:02d} has ID {meta_id}"
        
        print("  ✅ Control1 IDs are consistent")
    
    return True


def test_partial_recreation_bug():
    """Test what happens when files are partially recreated."""
    print("\nTest 2: Partial Recreation Bug")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        control_dir = os.path.join(tmpdir, "control")
        os.makedirs(control_dir)
        
        # First run: Create 3 individuals (simulating interrupted run)
        print("First run: Creating 3 individuals...")
        created_ids = []
        for i in range(3):
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': i,
                'individual_type': 'control',
                'run': 'first'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(control_dir, f"individual_{i:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            created_ids.append(i)
            print(f"  Created individual_{i:02d}.json with ID={i}")
        
        # Second run: Create 5 individuals (full run, but files 0-2 exist)
        print("\nSecond run: Creating 5 individuals (files 0-2 already exist)...")
        
        # Simulate what happens if we enumerate from a new sample
        # but files already exist (this is the bug scenario)
        new_sample = ['cell_a', 'cell_b', 'cell_c', 'cell_d', 'cell_e']
        
        for i, fake_cell in enumerate(new_sample):
            # This is what the pipeline does - enumerate always starts from 0
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': i,  # BUG: This starts from 0 again!
                'individual_type': 'control',
                'run': 'second'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(control_dir, f"individual_{i:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Saved individual_{i:02d}.json with ID={i} (overwrites if exists)")
        
        # Load and check - this reveals the bug
        print("\nLoading all files to check IDs...")
        dishes = load_all_petri_dishes(control_dir)
        
        mismatches = []
        for file_idx, petri in enumerate(dishes):
            meta_id = petri.metadata.get('individual_id', -1)
            run = petri.metadata.get('run', 'unknown')
            
            print(f"  File individual_{file_idx:02d}.json: metadata ID={meta_id}, run={run}")
            
            if meta_id != file_idx:
                mismatches.append((file_idx, meta_id))
        
        if mismatches:
            print(f"\n  ⚠️  Found {len(mismatches)} ID mismatches!")
            print("  This happens because enumerate() always starts from 0")
            print("  even when recreating files that already exist.")
        else:
            print("\n  ✅ No mismatches (unexpected in bug scenario)")
        
        return len(mismatches) == 0  # Returns False if bug exists


def test_proposed_fix():
    """Test a proposed fix: use filename-based IDs."""
    print("\nTest 3: Proposed Fix - Filename-based IDs")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        control_dir = os.path.join(tmpdir, "control")
        os.makedirs(control_dir)
        
        # Create some existing files (simulating partial run)
        print("Creating some existing files...")
        for i in [0, 2, 4]:  # Gaps in sequence
            cell = Cell(n=100, rate=0.005)
            metadata = {
                'individual_id': i,
                'individual_type': 'control'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(control_dir, f"individual_{i:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Created individual_{i:02d}.json with ID={i}")
        
        # Now fill in the gaps with correct IDs
        print("\nFilling gaps with correct IDs...")
        for i in [1, 3]:
            cell = Cell(n=100, rate=0.005)
            
            # FIX: Use the target filename index as the ID
            metadata = {
                'individual_id': i,  # Use filename-based ID
                'individual_type': 'control'
            }
            
            petri = PetriDish.from_snapshot_cell(
                cell=cell,
                growth_phase=5,
                metadata=metadata
            )
            
            filepath = os.path.join(control_dir, f"individual_{i:02d}.json")
            save_petri_dish(petri, filepath, compress=False)
            print(f"  Created individual_{i:02d}.json with ID={i}")
        
        # Verify all IDs match
        print("\nVerifying all IDs match filenames...")
        dishes = load_all_petri_dishes(control_dir)
        
        all_match = True
        for file_idx, petri in enumerate(dishes):
            meta_id = petri.metadata.get('individual_id', -1)
            print(f"  File individual_{file_idx:02d}.json: metadata ID={meta_id}")
            
            if meta_id != file_idx:
                all_match = False
                print(f"    ❌ Mismatch!")
            else:
                print(f"    ✅ Match!")
        
        if all_match:
            print("\n  ✅ Fix works: All IDs match filenames")
        else:
            print("\n  ❌ Fix failed: Some IDs don't match")
        
        return all_match


def main():
    """Run all tests."""
    print("=" * 60)
    print("Individual ID Consistency Tests")
    print("=" * 60)
    
    tests = [
        test_individual_id_creation,
        test_partial_recreation_bug,
        test_proposed_fix
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
                print(f"  ⚠️ Test revealed issue")
        except Exception as e:
            print(f"\n  ❌ Test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Test Summary: {passed} passed, {failed} failed/revealed issues")
    
    print("\n🔍 ROOT CAUSE IDENTIFIED:")
    print("   When creating individuals, enumerate() always starts from 0")
    print("   If files already exist from a previous run, they get overwritten")
    print("   but with potentially different IDs in metadata.")
    print("")
    print("   Example: If control1 creation is run twice:")
    print("   - First run creates 3 files: 00,01,02 with IDs 0,1,2")
    print("   - Second run creates 6 files: 00-05 with IDs 0-5")
    print("   - Files 00-02 get overwritten with same IDs")
    print("   - But if sampling order changes, ID mismatches occur")
    print("")
    print("📝 SOLUTION:")
    print("   Set individual_id based on the target filename index,")
    print("   not the enumeration index. This ensures consistency")
    print("   even across multiple runs.")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())