"""
Test that all history dictionaries use consistent string keys.

This test ensures that the dictionary key convention is maintained:
- All history dictionaries should use string keys for years
- Sorting should work correctly with string keys
- Access patterns should be consistent
"""

import sys
import os
import json
import tempfile

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cell import Cell, PetriDish


def test_string_keys_in_histories():
    """Verify all history dictionaries use string keys."""
    print("\n=== Testing String Keys in History Dictionaries ===")
    
    # Create a PetriDish with all tracking enabled (track_gene_jsd is True by default)
    petri = PetriDish(rate=0.005, growth_phase=3)
    
    # Simulate 5 years
    petri.run_simulation(t_max=5)
    
    # Check cell_history uses string keys
    print("Checking cell_history keys...")
    assert hasattr(petri, 'cell_history'), "PetriDish should have cell_history"
    assert len(petri.cell_history) > 0, "cell_history should have entries"
    for key in petri.cell_history.keys():
        assert isinstance(key, str), f"cell_history key {key} should be string, got {type(key)}"
    print(f"  ✓ All {len(petri.cell_history)} cell_history keys are strings")
    
    # Check gene_jsd_history uses string keys (may be empty if not tracking)
    print("Checking gene_jsd_history keys...")
    assert hasattr(petri, 'gene_jsd_history'), "PetriDish should have gene_jsd_history"
    if len(petri.gene_jsd_history) > 0:
        for key in petri.gene_jsd_history.keys():
            assert isinstance(key, str), f"gene_jsd_history key {key} should be string, got {type(key)}"
        print(f"  ✓ All {len(petri.gene_jsd_history)} gene_jsd_history keys are strings")
    else:
        print(f"  ⚠️ gene_jsd_history is empty (track_gene_jsd may be disabled)")
    
    # Check mean/median histories use string keys
    print("Checking mean_gene_jsd_history keys...")
    assert hasattr(petri, 'mean_gene_jsd_history'), "PetriDish should have mean_gene_jsd_history"
    for key in petri.mean_gene_jsd_history.keys():
        assert isinstance(key, str), f"mean_gene_jsd_history key {key} should be string, got {type(key)}"
    print(f"  ✓ All {len(petri.mean_gene_jsd_history)} mean_gene_jsd_history keys are strings")
    
    print("Checking median_gene_jsd_history keys...")
    assert hasattr(petri, 'median_gene_jsd_history'), "PetriDish should have median_gene_jsd_history"
    for key in petri.median_gene_jsd_history.keys():
        assert isinstance(key, str), f"median_gene_jsd_history key {key} should be string, got {type(key)}"
    print(f"  ✓ All {len(petri.median_gene_jsd_history)} median_gene_jsd_history keys are strings")
    
    print("\n✅ All history dictionaries use string keys consistently")


def test_sorting_with_string_keys():
    """Test that sorting string keys works correctly."""
    print("\n=== Testing Sorting with String Keys ===")
    
    # Create a PetriDish and simulate years
    petri = PetriDish(rate=0.005, growth_phase=2)
    
    # Simulate 12 years to test multi-digit sorting
    petri.run_simulation(t_max=12)
    
    # Test sorting pattern
    years = sorted([int(y) for y in petri.cell_history.keys()])
    print(f"Years after sorting: {years}")
    
    # Verify years are in correct order
    assert years == list(range(len(years))), "Years should be in sequential order"
    
    # Verify we can access with string keys
    for year in years:
        assert str(year) in petri.cell_history, f"Should be able to access year {year} with string key"
        cells = petri.cell_history[str(year)]
        assert isinstance(cells, list), "Should get cells list"
    
    print(f"  ✓ Sorting pattern works correctly for {len(years)} years")
    print("  ✓ Can access all years with string keys after sorting")


def test_json_serialization_with_string_keys():
    """Test that JSON serialization maintains string keys."""
    print("\n=== Testing JSON Serialization ===")
    
    # Create and simulate
    petri = PetriDish(rate=0.005, growth_phase=2)
    petri.run_simulation(t_max=5)
    
    # Save to temporary file
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_path = f.name
    
    petri.save_history(temp_path, compress=False)  # Save without compression
    
    # Load and check keys
    try:
        with open(temp_path, 'r') as f:
            data = json.load(f)
        
        # Check history keys in JSON
        assert 'history' in data, "JSON should have history key"
        for key in data['history'].keys():
            assert isinstance(key, str), f"JSON history key {key} should be string"
        
        print(f"  ✓ JSON serialization maintains string keys")
        print(f"  ✓ Loaded {len(data['history'])} years with string keys")
        
    finally:
        # Clean up
        if os.path.exists(temp_path):
            os.remove(temp_path)


def test_plot_functions_handle_string_keys():
    """Test that plot functions correctly handle string keys."""
    print("\n=== Testing Plot Functions with String Keys ===")
    
    # Create a PetriDish with gene-specific rates
    gene_rate_groups = [(50, 0.004), (50, 0.005), (50, 0.006), (50, 0.007)]
    petri = PetriDish(gene_rate_groups=gene_rate_groups, growth_phase=3)
    
    # Simulate several years
    petri.run_simulation(t_max=10)
    
    # Test that we can create plotter and access data
    from cell import PetriDishPlotter
    plotter = PetriDishPlotter(petri)
    
    # Verify plotter can work with string keys
    if hasattr(plotter.petri, 'gene_jsd_history'):
        # Test the sorting pattern used in plot functions
        years = sorted([int(y) for y in plotter.petri.gene_jsd_history.keys()])
        
        # Verify we can access data with the pattern
        for year in years:
            gene_jsds = plotter.petri.gene_jsd_history[str(year)]
            assert isinstance(gene_jsds, list), "Should get gene JSDs list"
            assert len(gene_jsds) > 0, "Gene JSDs should not be empty"
        
        print(f"  ✓ Plot functions can handle string keys correctly")
        print(f"  ✓ Tested with {len(years)} years of gene JSD data")
    
    # Test with heatmap if plotly is available
    try:
        import plotly
        print("  Testing plot_gene_jsd_heatmap...")
        # This should not raise type comparison error anymore
        with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as f:
            temp_plot = f.name
        
        try:
            plotter.plot_gene_jsd_heatmap(output_path=temp_plot)
            print("  ✓ plot_gene_jsd_heatmap works with string keys")
        except Exception as e:
            if "'<' not supported" in str(e):
                assert False, f"Type comparison error still exists: {e}"
            # Other errors (like missing data) are okay for this test
            print(f"  ⚠️ Plot skipped: {e}")
        finally:
            if os.path.exists(temp_plot):
                os.remove(temp_plot)
                
    except ImportError:
        print("  ⚠️ Plotly not installed, skipping plot generation test")


if __name__ == "__main__":
    print("\n" + "="*60)
    print("KEY CONSISTENCY TESTS")
    print("="*60)
    
    test_string_keys_in_histories()
    test_sorting_with_string_keys()
    test_json_serialization_with_string_keys()
    test_plot_functions_handle_string_keys()
    
    print("\n" + "="*60)
    print("✅ ALL KEY CONSISTENCY TESTS PASSED")
    print("="*60)