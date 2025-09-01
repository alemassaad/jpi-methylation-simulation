#!/usr/bin/env python3
"""
Test automatic file format detection for both compressed and uncompressed files.
"""

import os
import sys
import tempfile
import json
import gzip

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

from cell import Cell, PetriDish
from pipeline_utils import smart_open, save_petri_dish, load_petri_dish, save_snapshot_cells, load_snapshot_cells


def test_smart_open():
    """Test smart_open utility function."""
    print("\nTest 1: Smart Open Utility")
    print("=" * 50)
    
    # Test with .json file
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        json_file = f.name
    
    try:
        # Write with smart_open
        data = {"test": "uncompressed"}
        with smart_open(json_file, 'w') as f:
            json.dump(data, f)
        
        # Read with smart_open
        with smart_open(json_file, 'r') as f:
            loaded = json.load(f)
        
        assert loaded == data, "Uncompressed JSON failed"
        print("  ✓ Uncompressed .json files work")
        
    finally:
        os.remove(json_file)
    
    # Test with .json.gz file
    with tempfile.NamedTemporaryFile(suffix='.json.gz', delete=False) as f:
        gz_file = f.name
    
    try:
        # Write with smart_open
        data = {"test": "compressed"}
        with smart_open(gz_file, 'w') as f:
            json.dump(data, f)
        
        # Read with smart_open
        with smart_open(gz_file, 'r') as f:
            loaded = json.load(f)
        
        assert loaded == data, "Compressed JSON failed"
        print("  ✓ Compressed .json.gz files work")
        
    finally:
        os.remove(gz_file)
    
    # Test error on wrong extension
    try:
        with smart_open("test.txt", 'r') as f:
            pass
        print("  ✗ Should have raised error for .txt file")
    except ValueError as e:
        print(f"  ✓ Correctly rejected .txt file: {e}")
    
    print("  PASSED")


def test_petri_dish_save_load():
    """Test PetriDish saving and loading with both formats."""
    print("\nTest 2: PetriDish Save/Load")
    print("=" * 50)
    
    # Create a PetriDish
    petri = PetriDish(rate=0.01, n=100, seed=42)
    petri.grow_with_homeostasis(years=2, growth_phase=1, verbose=False)
    original_cells = len(petri.cells)
    
    # Test compressed save/load
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        gz_file = f.name + '.gz'
    
    try:
        # Save compressed
        save_petri_dish(petri, gz_file[:-3], compress=True)  # Pass without .gz
        assert os.path.exists(gz_file), "Compressed file not created"
        print(f"  ✓ Saved compressed: {os.path.basename(gz_file)}")
        
        # Load compressed
        loaded = load_petri_dish(gz_file)
        assert len(loaded.cells) == original_cells, "Cell count mismatch"
        print(f"  ✓ Loaded compressed: {len(loaded.cells)} cells")
        
    finally:
        if os.path.exists(gz_file):
            os.remove(gz_file)
    
    # Test uncompressed save/load
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        json_file = f.name
    
    try:
        # Save uncompressed
        save_petri_dish(petri, json_file, compress=False)
        assert os.path.exists(json_file), "Uncompressed file not created"
        assert not json_file.endswith('.gz'), "Should not have .gz extension"
        print(f"  ✓ Saved uncompressed: {os.path.basename(json_file)}")
        
        # Load uncompressed
        loaded = load_petri_dish(json_file)
        assert len(loaded.cells) == original_cells, "Cell count mismatch"
        print(f"  ✓ Loaded uncompressed: {len(loaded.cells)} cells")
        
    finally:
        if os.path.exists(json_file):
            os.remove(json_file)
    
    print("  PASSED")


def test_snapshot_save_load():
    """Test snapshot saving and loading with both formats."""
    print("\nTest 3: Snapshot Save/Load")
    print("=" * 50)
    
    # Create some cells
    cells = []
    for i in range(10):
        cell = Cell(n=100, rate=0.01)
        cell.methylate()
        cells.append(cell)
    
    # Test compressed
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        gz_file = f.name + '.gz'
    
    try:
        # Save compressed
        save_snapshot_cells(cells, gz_file[:-3], compress=True)
        assert os.path.exists(gz_file), "Compressed snapshot not created"
        print(f"  ✓ Saved compressed snapshot: {os.path.basename(gz_file)}")
        
        # Load compressed
        loaded = load_snapshot_cells(gz_file)
        assert len(loaded) == len(cells), "Cell count mismatch"
        print(f"  ✓ Loaded compressed snapshot: {len(loaded)} cells")
        
    finally:
        if os.path.exists(gz_file):
            os.remove(gz_file)
    
    # Test uncompressed
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        json_file = f.name
    
    try:
        # Save uncompressed
        save_snapshot_cells(cells, json_file, compress=False)
        assert os.path.exists(json_file), "Uncompressed snapshot not created"
        print(f"  ✓ Saved uncompressed snapshot: {os.path.basename(json_file)}")
        
        # Load uncompressed
        loaded = load_snapshot_cells(json_file)
        assert len(loaded) == len(cells), "Cell count mismatch"
        print(f"  ✓ Loaded uncompressed snapshot: {len(loaded)} cells")
        
    finally:
        if os.path.exists(json_file):
            os.remove(json_file)
    
    print("  PASSED")


def test_mixed_formats():
    """Test that we can load either format transparently."""
    print("\nTest 4: Mixed Format Loading")
    print("=" * 50)
    
    # Create test data
    petri = PetriDish(rate=0.01, n=100, seed=42)
    petri.grow_with_homeostasis(years=2, growth_phase=1, verbose=False)
    
    # Save in both formats
    with tempfile.TemporaryDirectory() as tmpdir:
        json_path = os.path.join(tmpdir, "test.json")
        gz_path = os.path.join(tmpdir, "test.json.gz")
        
        # Save both
        save_petri_dish(petri, json_path, compress=False)
        save_petri_dish(petri, json_path, compress=True)  # This creates .json.gz
        
        assert os.path.exists(json_path), "JSON file not created"
        assert os.path.exists(gz_path), "GZ file not created"
        
        # Load both and verify they're identical
        loaded_json = load_petri_dish(json_path)
        loaded_gz = load_petri_dish(gz_path)
        
        assert len(loaded_json.cells) == len(loaded_gz.cells), "Different cell counts"
        assert len(loaded_json.cells) == len(petri.cells), "Cell count changed"
        
        # Check file sizes
        json_size = os.path.getsize(json_path)
        gz_size = os.path.getsize(gz_path)
        compression_ratio = gz_size / json_size
        
        print(f"  ✓ JSON size: {json_size:,} bytes")
        print(f"  ✓ GZ size: {gz_size:,} bytes")
        print(f"  ✓ Compression ratio: {compression_ratio:.1%}")
        assert gz_size < json_size, "Compressed should be smaller"
    
    print("  PASSED")


def test_phase1_simulation_format():
    """Test loading Phase 1 simulation files in both formats."""
    print("\nTest 5: Phase 1 Simulation Format")
    print("=" * 50)
    
    # Create a minimal simulation-like structure
    sim_data = {
        "parameters": {
            "rate": 0.005,
            "n": 100,
            "gene_size": 5,
            "years": 10,
            "seed": 42
        },
        "history": {
            "5": {
                "cells": [
                    {"methylated": [0] * 100, "cell_JSD": 0.01},
                    {"methylated": [1, 0] * 50, "cell_JSD": 0.02}
                ]
            }
        }
    }
    
    # Test compressed format
    with tempfile.NamedTemporaryFile(suffix='.json.gz', delete=False) as f:
        gz_file = f.name
    
    try:
        # Save compressed
        with smart_open(gz_file, 'w') as f:
            json.dump(sim_data, f)
        
        # Load and verify
        from pipeline_utils import load_snapshot_as_cells
        cells = load_snapshot_as_cells(gz_file, year=5)
        assert len(cells) == 2, f"Expected 2 cells, got {len(cells)}"
        print(f"  ✓ Loaded compressed simulation: {len(cells)} cells")
        
    finally:
        os.remove(gz_file)
    
    # Test uncompressed format
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        json_file = f.name
    
    try:
        # Save uncompressed
        with smart_open(json_file, 'w') as f:
            json.dump(sim_data, f)
        
        # Load and verify
        cells = load_snapshot_as_cells(json_file, year=5)
        assert len(cells) == 2, f"Expected 2 cells, got {len(cells)}"
        print(f"  ✓ Loaded uncompressed simulation: {len(cells)} cells")
        
    finally:
        os.remove(json_file)
    
    print("  PASSED")


def run_all_tests():
    """Run all file format detection tests."""
    print("\n" + "="*60)
    print("TESTING FILE FORMAT DETECTION")
    print("="*60)
    
    tests = [
        test_smart_open,
        test_petri_dish_save_load,
        test_snapshot_save_load,
        test_mixed_formats,
        test_phase1_simulation_format
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"\n  FAILED: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    print(f"SUMMARY: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)