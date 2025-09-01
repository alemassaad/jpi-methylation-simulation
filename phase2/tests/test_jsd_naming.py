#!/usr/bin/env python3
"""
Test that the new JSD naming convention is working correctly.
"""

import os
import sys
import json
import tempfile

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, load_petri_dish, get_petri_statistics, get_cell_jsd_array

def test_metadata_keys():
    """Test that save_petri_dish creates correct metadata keys."""
    print("Testing metadata key names...")
    
    # Create a test PetriDish
    petri = PetriDish(rate=0.005, growth_phase=2)
    petri.divide_cells()
    petri.methylate_cells()
    
    # Save to temp file
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        # Save without compression
        save_petri_dish(petri, temp_file, compress=False)
        
        # Load the JSON directly to check keys
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        # Check metadata keys
        metadata = data['metadata']
        expected_keys = ['mean_cell_jsd', 'std_cell_jsd', 'min_cell_jsd', 
                        'max_cell_jsd', 'median_cell_jsd']
        old_keys = ['mean_jsd', 'std_jsd', 'min_jsd', 'max_jsd', 'median_jsd']
        
        print("\n  Checking for new cell JSD keys:")
        for key in expected_keys:
            if key in metadata:
                print(f"    ✓ Found {key}: {metadata[key]:.4f}")
            else:
                print(f"    ✗ Missing {key}")
                return False
        
        print("\n  Checking old keys are gone:")
        for key in old_keys:
            if key in metadata:
                print(f"    ✗ Found old key {key} (should be removed)")
                return False
            else:
                print(f"    ✓ Old key {key} not present")
        
        print("\n  ✅ All metadata keys correct!")
        return True
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)

def test_function_names():
    """Test that renamed functions work correctly."""
    print("\nTesting function names...")
    
    # Create test PetriDish
    petri = PetriDish(rate=0.005, growth_phase=2)
    petri.divide_cells()
    petri.methylate_cells()
    
    # Test get_cell_jsd_array
    try:
        jsd_array = get_cell_jsd_array(petri)
        print(f"  ✓ get_cell_jsd_array() works: {len(jsd_array)} values")
    except AttributeError as e:
        if "get_jsd_array" in str(e):
            print(f"  ✗ Old function name still being used: {e}")
            return False
        raise
    
    # Test get_petri_statistics
    stats = get_petri_statistics(petri)
    
    # Check that stats has new keys
    expected_keys = ['mean_cell_jsd', 'std_cell_jsd', 'min_cell_jsd', 
                    'max_cell_jsd', 'median_cell_jsd']
    
    print("\n  Checking get_petri_statistics() keys:")
    for key in expected_keys:
        if key in stats:
            print(f"    ✓ Found {key}: {stats[key]:.4f}")
        else:
            print(f"    ✗ Missing {key}")
            return False
    
    print("\n  ✅ All function names correct!")
    return True

def test_pipeline_analysis_import():
    """Test that pipeline_analysis functions are renamed."""
    print("\nTesting pipeline_analysis imports...")
    
    try:
        from pipeline_analysis import plot_cell_jsd_distribution, get_mean_cell_jsds_from_petri_dishes
        print("  ✓ plot_cell_jsd_distribution imported")
        print("  ✓ get_mean_cell_jsds_from_petri_dishes imported")
        
        # Try to import old names (should fail)
        try:
            from pipeline_analysis import plot_jsd_distribution_from_cells
            print("  ✗ Old function plot_jsd_distribution_from_cells still exists")
            return False
        except ImportError:
            print("  ✓ Old function plot_jsd_distribution_from_cells removed")
        
        try:
            from pipeline_analysis import get_mean_jsds_from_petri_dishes
            print("  ✗ Old function get_mean_jsds_from_petri_dishes still exists")
            return False
        except ImportError:
            print("  ✓ Old function get_mean_jsds_from_petri_dishes removed")
        
        print("\n  ✅ Pipeline analysis functions renamed correctly!")
        return True
        
    except ImportError as e:
        print(f"  ✗ Import error: {e}")
        return False

def main():
    """Run all naming convention tests."""
    print("="*60)
    print("Testing JSD Naming Convention Changes")
    print("="*60)
    
    tests_passed = 0
    tests_failed = 0
    
    # Test 1: Metadata keys
    if test_metadata_keys():
        tests_passed += 1
    else:
        tests_failed += 1
    
    # Test 2: Function names
    if test_function_names():
        tests_passed += 1
    else:
        tests_failed += 1
    
    # Test 3: Pipeline analysis imports
    if test_pipeline_analysis_import():
        tests_passed += 1
    else:
        tests_failed += 1
    
    # Summary
    print("\n" + "="*60)
    print(f"Test Summary: {tests_passed} passed, {tests_failed} failed")
    
    if tests_failed == 0:
        print("🎉 All JSD naming convention tests passed!")
        return 0
    else:
        print(f"⚠️  {tests_failed} test(s) failed")
        return 1

if __name__ == "__main__":
    sys.exit(main())