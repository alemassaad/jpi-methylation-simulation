#!/usr/bin/env python3
"""
Simple test to verify imports and basic functionality after cleanup.
"""

import os
import sys

# Add paths
phase2_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
phase1_dir = os.path.join(os.path.dirname(phase2_dir), 'phase1')
sys.path.insert(0, phase2_dir)
sys.path.append(phase1_dir)

def test_imports():
    """Test that all necessary imports still work."""
    print("Testing imports after cleanup...")
    try:
        # Core imports from phase1
        from cell import PetriDish, Cell, PetriDishPlotter, rate_to_gene_rate_groups
        
        # Helper functions
        from core.individual_helpers import (
            create_individual, process_batch_growth, process_batch_mixing
        )
        
        # Pipeline utilities - with removed unused imports
        from core.pipeline_utils import (
            load_snapshot_as_cells, save_snapshot_cells, load_snapshot_cells,
            save_petri_dish, load_petri_dish, load_all_petri_dishes,
            sample_by_quantiles, sample_uniform,
            create_pure_snapshot_petri,
            create_control2_with_uniform_base,
            check_petri_files_state,
            calculate_population_statistics, print_mixing_statistics,
            create_uniform_mixing_pool, mix_petri_uniform, normalize_populations,
            normalize_individuals_for_uniform_mixing
        )
        
        # Analysis functions
        from core.pipeline_analysis import (
            plot_cell_jsd_distribution,
            analyze_populations_from_dishes,
            plot_gene_jsd_distributions,
            plot_gene_jsd_individual_comparison
        )
        
        # Path utilities
        from core.path_utils import parse_step1_simulation_path, generate_step23_output_dir
        
        # Validation
        from core.validation import PipelineValidator, ValidationError
        
        print("  ✓ All imports successful")
        return True
    except ImportError as e:
        print(f"  ✗ Import failed: {e}")
        return False

def test_helper_functions_still_work():
    """Test that the helper functions that use the 'removed' functions still work."""
    print("\nTesting that individual_helpers still work...")
    try:
        from cell import Cell, PetriDish
        from core.individual_helpers import create_individual, process_batch_growth
        
        # Create a test cell
        cell = Cell(n=20, rate=0.005)
        cell.age = 10
        
        # Create an individual (uses internal functions)
        petri = create_individual(
            cell=cell,
            individual_type='test',
            individual_id=1,
            growth_phase=2,
            additional_metadata={'test': True}
        )
        
        if petri and len(petri.cells) == 1:
            print("  ✓ create_individual works")
        else:
            print("  ✗ create_individual failed")
            return False
        
        # Test growth (uses grow_petri_for_years internally)
        dishes = [petri]
        process_batch_growth(
            dishes=dishes,
            batch_name='test',
            years=2,
            growth_phase=2,
            expected_population=4,
            start_year=10,
            verbose=False
        )
        
        if len(dishes[0].cells) > 1:
            print("  ✓ process_batch_growth works (uses grow_petri_for_years internally)")
        else:
            print("  ✗ process_batch_growth failed")
            return False
        
        return True
        
    except Exception as e:
        print(f"  ✗ Test failed: {e}")
        return False

def main():
    """Run all tests."""
    print("="*60)
    print("Testing Pipeline After Cleanup")
    print("="*60)
    
    all_passed = True
    
    # Test imports
    if not test_imports():
        all_passed = False
    
    # Test that helper functions still work
    if not test_helper_functions_still_work():
        all_passed = False
    
    print("\n" + "="*60)
    if all_passed:
        print("✅ All tests passed! Cleanup was successful.")
        print("\nRemoved files:")
        print("  - diagnose_control1.py")
        print("  - test_control2_fix.py")
        print("  - test_control2_metadata.py")
        print("  - test_control2_no_history.py")
        print("  - test_gene_rate_preservation.py")
        print("  - test_load_save_fix.py")
        print("\nRemoved unused imports from run_pipeline.py:")
        print("  - get_cell_jsd_array")
        print("  - get_petri_statistics")
        print("  - grow_petri_for_years")
        print("  - mix_petri_with_snapshot")
        print("\nNote: These functions are still defined in pipeline_utils.py")
        print("because they're used internally by other functions.")
        return 0
    else:
        print("❌ Some tests failed. Check cleanup changes.")
        return 1

if __name__ == "__main__":
    exit(main())