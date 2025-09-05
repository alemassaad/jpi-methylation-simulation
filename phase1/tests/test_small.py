#!/usr/bin/env python3
"""
Quick test of Phase 1 simulation with small parameters.
Tests the growth phase and steady state transitions.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cell import PetriDish


def test_small_simulation():
    """Run a small test simulation to verify functionality."""
    print("Running small test simulation...")
    print("-" * 40)
    
    # Small parameters for quick testing
    petri_dish = PetriDish(
        rate=0.01,      # Higher rate for visible changes
        n=100,          # Fewer CpG sites
        gene_size=5,    # Standard gene size
        seed=42,        # Fixed seed for reproducibility
        growth_phase=4  # Small target (2^4 = 16 cells)
    )
    
    # Run for 10 years to test both phases
    petri_dish.run_simulation(t_max=10)
    
    # Check growth phase worked correctly
    print("\n" + "="*40)
    print("TEST RESULTS:")
    print("="*40)
    
    # Year-by-year population check
    expected_pops = {
        0: 1,    # Initial
        1: 2,    # 2^1
        2: 4,    # 2^2
        3: 8,    # 2^3
        4: 16,   # 2^4 - reaches target
        # After year 4, should maintain ~16 cells
    }
    
    print("\nPopulation verification:")
    for year in range(11):
        if str(year) in petri_dish.cell_history:
            actual = len(petri_dish.cell_history[str(year)])
            if year in expected_pops:
                expected = expected_pops[year]
                status = "✓" if actual == expected else "✗"
                print(f"  Year {year}: {actual} cells (expected: {expected}) {status}")
            else:
                # Steady state - should be around 16
                in_range = 8 <= actual <= 24  # Allow variation
                status = "✓" if in_range else "✗"
                print(f"  Year {year}: {actual} cells (steady state ~16) {status}")
    
    # Test methylation is happening
    print("\nMethylation verification:")
    year_10_cells = petri_dish.cell_history.get('10', [])
    if year_10_cells:
        # In lean format, check if any sites are methylated
        methylated_cells = sum(1 for cell in year_10_cells if any(cell['methylated']))
        print(f"  Cells with methylation at year 10: {methylated_cells}/{len(year_10_cells)}")
        
        # Check JSD is increasing
        mean_jsd_year_0 = petri_dish.cell_history['0'][0]['cell_jsd']
        mean_jsd_year_10 = sum(c['cell_jsd'] for c in year_10_cells) / len(year_10_cells)
        print(f"  Mean JSD year 0: {mean_jsd_year_0:.4f}")
        print(f"  Mean JSD year 10: {mean_jsd_year_10:.4f}")
        
        if mean_jsd_year_10 > mean_jsd_year_0:
            print("  ✓ JSD is increasing over time")
        else:
            print("  ✗ JSD should increase over time")
    
    print("\n" + "="*40)
    print("Test complete!")
    return petri_dish


if __name__ == "__main__":
    petri_dish = test_small_simulation()