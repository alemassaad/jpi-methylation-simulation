#!/usr/bin/env python3
"""
Minimal test to verify compression consistency without plotly dependency.
"""

import os
import sys
import glob
import json
import gzip

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from pipeline_utils import save_petri_dish, load_petri_dish, save_snapshot_cells, load_snapshot_cells
from cell import Cell, PetriDish

def test_save_functions():
    """Test that save functions respect compression parameter."""
    print("Testing save functions with compression parameter")
    print("=" * 50)
    
    # Create a test PetriDish
    petri = PetriDish(rate=0.005, growth_phase=2)
    petri.divide_cells()
    petri.methylate_cells()
    
    # Test 1: Save with compression=False
    print("\nTest 1: save_petri_dish with compress=False")
    save_petri_dish(petri, "./test_uncompressed.json", compress=False)
    assert os.path.exists("test_uncompressed.json"), "File should exist"
    assert not os.path.exists("test_uncompressed.json.gz"), "Compressed file should not exist"
    
    # Verify it's actually uncompressed
    with open("test_uncompressed.json", 'r') as f:
        data = json.load(f)
        assert 'cells' in data, "Should be valid JSON"
    print("  ✅ Created uncompressed .json file")
    
    # Test 2: Save with compression=True
    print("\nTest 2: save_petri_dish with compress=True")
    save_petri_dish(petri, "./test_compressed.json.gz", compress=True)
    assert os.path.exists("test_compressed.json.gz"), "Compressed file should exist"
    
    # Verify it's actually compressed
    with gzip.open("test_compressed.json.gz", 'rt') as f:
        data = json.load(f)
        assert 'cells' in data, "Should be valid JSON"
    print("  ✅ Created compressed .json.gz file")
    
    # Test 3: Load both formats
    print("\nTest 3: Loading both formats")
    petri1 = load_petri_dish("test_uncompressed.json")
    petri2 = load_petri_dish("test_compressed.json.gz")
    assert len(petri1.cells) == len(petri2.cells), "Should load same data"
    print("  ✅ Both formats load correctly")
    
    # Test 4: Save snapshot cells
    print("\nTest 4: save_snapshot_cells with compression")
    cells = [Cell(n=100, rate=0.005) for _ in range(2)]
    
    save_snapshot_cells(cells, "./snapshot_uncompressed.json", compress=False)
    assert os.path.exists("snapshot_uncompressed.json")
    print("  ✅ Created uncompressed snapshot")
    
    save_snapshot_cells(cells, "./snapshot_compressed.json.gz", compress=True)
    assert os.path.exists("snapshot_compressed.json.gz")
    print("  ✅ Created compressed snapshot")
    
    # Clean up
    for f in ["test_uncompressed.json", "test_compressed.json.gz", 
              "snapshot_uncompressed.json", "snapshot_compressed.json.gz"]:
        if os.path.exists(f):
            os.remove(f)
    
    print("\n" + "=" * 50)
    print("All save function tests passed! 🎉")

def check_directory_compression(directory, expected_compressed):
    """Check if all JSON files in directory match expected compression."""
    json_files = glob.glob(os.path.join(directory, "**/*.json"), recursive=True)
    gz_files = glob.glob(os.path.join(directory, "**/*.json.gz"), recursive=True)
    
    if expected_compressed:
        if json_files:
            print(f"  ❌ Found uncompressed files when expecting compressed:")
            for f in json_files[:3]:  # Show first 3
                print(f"     - {f}")
            return False
        print(f"  ✅ All {len(gz_files)} files are compressed (.json.gz)")
    else:
        if gz_files:
            print(f"  ❌ Found compressed files when expecting uncompressed:")
            for f in gz_files[:3]:  # Show first 3
                print(f"     - {f}")
            return False
        print(f"  ✅ All {len(json_files)} files are uncompressed (.json)")
    
    return True

if __name__ == "__main__":
    test_save_functions()