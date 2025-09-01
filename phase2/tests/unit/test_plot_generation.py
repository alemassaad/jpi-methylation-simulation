#!/usr/bin/env python3
"""
Test that plot_individuals.py can now generate plots with the fix.
"""

import os
import sys
import tempfile
import shutil

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, grow_petri_for_years
from plot_individuals import plot_individual

def test_plot_generation():
    """Test that plots can be generated with the fixed attribute name."""
    print("Testing plot generation with fixed history attribute...")
    
    # Create temporary directory structure
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create directories
        individuals_dir = os.path.join(tmpdir, "individuals")
        mutant_dir = os.path.join(individuals_dir, "mutant")
        plots_dir = os.path.join(tmpdir, "individual_plots")
        os.makedirs(mutant_dir)
        os.makedirs(plots_dir)
        
        # Create a PetriDish with history
        print("  Creating PetriDish with history...")
        petri = PetriDish(rate=0.005, n=100, growth_phase=2)
        
        # Grow with history tracking
        grow_petri_for_years(
            petri,
            years=5,
            growth_phase=2,
            verbose=False,
            track_history=True,
            start_year=0
        )
        
        # Verify history exists
        assert hasattr(petri, 'cell_history'), "PetriDish should have cell_history"
        assert len(petri.cell_history) > 0, "History should have entries"
        print(f"  ✓ Created PetriDish with {len(petri.cell_history)} history points")
        
        # Save the individual
        filepath = os.path.join(mutant_dir, "individual_00.json.gz")
        save_petri_dish(petri, filepath, include_cell_history=True)
        print(f"  ✓ Saved individual with history")
        
        # Try to generate plots
        print("  Attempting to generate plots...")
        success = plot_individual(filepath, plots_dir, "mutant", 0)
        
        if success:
            print("  ✓ Plots generated successfully!")
            
            # Check that plot files were created
            expected_files = [
                "mutant_00_jsd.png",
                "mutant_00_methylation.png",
                "mutant_00_combined.png"
            ]
            
            for filename in expected_files:
                plot_path = os.path.join(plots_dir, filename)
                if os.path.exists(plot_path):
                    print(f"    ✓ Created: {filename}")
                else:
                    print(f"    ✗ Missing: {filename}")
                    success = False
        else:
            print("  ✗ Failed to generate plots")
        
        return success


def run_test():
    """Run the plot generation test."""
    print("="*60)
    print("Testing Plot Generation with History Fix")
    print("="*60)
    
    try:
        if test_plot_generation():
            print("\n" + "="*60)
            print("✅ Plot generation test passed!")
            print("\nThe history attribute fix is working:")
            print("  • plot_individuals.py now correctly looks for 'cell_history'")
            print("  • Plots can be generated from individuals with history")
            print("  • All three plot types are created (JSD, methylation, combined)")
            return 0
        else:
            print("\n" + "="*60)
            print("❌ Plot generation test failed")
            return 1
    except Exception as e:
        print(f"\n✗ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return 1


if __name__ == "__main__":
    sys.exit(run_test())