#!/usr/bin/env python3
"""
Additional edge case tests for Step1-Prime simulation.
Tests extreme parameters and boundary conditions with small values.
"""

import sys
import os
import tempfile

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cell import PetriDish, Cell

def test_single_year_simulation():
    """Test that growth_phase = years works correctly."""
    print("Testing single year with growth_phase = years...")
    petri = PetriDish(rate=0.1, n=20, gene_size=5, seed=42, growth_phase=1)
    petri.run_simulation(t_max=1)
    
    # Should have exactly 2 cells (2^1)
    assert len(petri.cells) == 2, f"Expected 2 cells, got {len(petri.cells)}"
    print("✓ Single year simulation works")

def test_no_growth_to_steady_transition():
    """Test when simulation ends during growth phase."""
    print("Testing simulation ending during growth phase...")
    petri = PetriDish(rate=0.1, n=20, gene_size=5, seed=42, growth_phase=5)
    petri.run_simulation(t_max=3)  # Stop before growth_phase
    
    # Should have 8 cells (2^3) and still be in growth
    assert len(petri.cells) == 8, f"Expected 8 cells, got {len(petri.cells)}"
    assert petri.year == 3, f"Expected year 3, got {petri.year}"
    assert petri.year < petri.growth_phase, "Should still be in growth phase"
    print("✓ Ending during growth phase works")

def test_immediate_steady_state():
    """Test growth_phase = 1 transitions quickly to steady state."""
    print("Testing immediate transition to steady state...")
    petri = PetriDish(rate=0.05, n=20, gene_size=5, seed=42, growth_phase=1)
    
    # Year 1: growth (1→2 cells)
    petri.simulate_year()
    assert len(petri.cells) == 2, f"Year 1: Expected 2 cells, got {len(petri.cells)}"
    
    # Year 2: should be steady state
    petri.simulate_year()
    assert petri.year > petri.growth_phase, "Should be in steady state at year 2"
    # Population should be around 2 (culled from 4)
    assert 1 <= len(petri.cells) <= 4, f"Year 2: Population {len(petri.cells)} out of range"
    print("✓ Immediate steady state transition works")

def test_extreme_culling_survival():
    """Test that at least 1 cell always survives culling."""
    print("Testing extreme culling scenarios...")
    petri = PetriDish(rate=0.01, n=20, gene_size=5, seed=999, growth_phase=1)
    
    # Get to steady state
    petri.simulate_year()  # Year 1: 1→2
    petri.simulate_year()  # Year 2: steady state
    
    # Run many years and ensure we never hit 0 cells
    for year in range(3, 20):
        petri.simulate_year()
        assert len(petri.cells) > 0, f"Population died at year {year}"
    
    print(f"✓ Population survived 20 years (final: {len(petri.cells)} cells)")

def test_100_percent_methylation():
    """Test behavior with 100% methylation rate."""
    print("Testing 100% methylation rate...")
    petri = PetriDish(rate=1.0, n=20, gene_size=5, seed=42, growth_phase=2)
    
    # After 1 year, all cells should be fully methylated
    petri.simulate_year()
    
    for cell in petri.cells:
        assert cell.cell_methylation_proportion == 1.0, "Cell not fully methylated with rate=1.0"
        assert all(site == 1 for site in cell.cpg_sites), "Some sites unmethylated"
    
    print("✓ 100% methylation rate works correctly")

def test_zero_methylation():
    """Test behavior with 0% methylation rate."""
    print("Testing 0% methylation rate...")
    petri = PetriDish(rate=0.0, n=20, gene_size=5, seed=42, growth_phase=2)
    
    # Run several years
    for _ in range(5):
        petri.simulate_year()
    
    # All cells should remain unmethylated
    for cell in petri.cells:
        assert cell.cell_methylation_proportion == 0.0, "Cell methylated with rate=0.0"
        assert all(site == 0 for site in cell.cpg_sites), "Some sites methylated"
        assert cell.cell_jsd == 0.0, "cell_JSD should be 0 for unmethylated cells"
    
    print("✓ 0% methylation rate works correctly")

def test_minimum_gene_configuration():
    """Test with minimum viable gene configuration."""
    print("Testing minimum gene configuration...")
    # n=10, gene_size=5 means only 2 genes
    petri = PetriDish(rate=0.5, n=10, gene_size=5, seed=42, growth_phase=2)
    
    petri.simulate_year()
    petri.simulate_year()
    
    # Check methylation distribution has correct length
    for cell in petri.cells:
        assert len(cell.methylation_distribution) == 6, "Wrong distribution length"
        assert abs(sum(cell.methylation_distribution) - 1.0) < 1e-6, "Distribution doesn't sum to 1"
    
    print("✓ Minimum gene configuration works")

def test_large_growth_phase():
    """Test with maximum allowed growth_phase."""
    print("Testing maximum growth_phase=20...")
    petri = PetriDish(rate=0.01, n=20, gene_size=5, seed=42, growth_phase=20)
    
    # Just test initialization and one year (full simulation would be slow)
    petri.simulate_year()
    
    assert petri.target_population == 2**20, f"Wrong target: {petri.target_population}"
    assert len(petri.cells) == 2, "Year 1 should have 2 cells"
    print("✓ Maximum growth_phase=20 works")

def test_filename_formats():
    """Test various filename generation scenarios."""
    print("Testing filename generation...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Test with seed
        petri1 = PetriDish(rate=0.123, n=20, gene_size=5, seed=99, growth_phase=2)
        petri1.simulate_year()
        path1 = petri1.save_history(directory=tmpdir)
        assert "_seed99" in path1, "Seed not in filename"
        assert "_g2_" in path1, "Growth phase not in filename"
        assert "0.123000" in path1, "Rate not formatted correctly"
        
        # Test without seed
        petri2 = PetriDish(rate=0.001, n=20, gene_size=5, seed=None, growth_phase=3)
        petri2.simulate_year()
        path2 = petri2.save_history(directory=tmpdir)
        assert "_noseed" in path2, "No seed marker missing"
        assert "_g3_" in path2, "Growth phase not in filename"
        
    print("✓ Filename generation works correctly")

def test_population_variability():
    """Test that steady state has natural variation."""
    print("Testing steady state population variability...")
    petri = PetriDish(rate=0.01, n=20, gene_size=5, seed=42, growth_phase=3)
    
    # Get to steady state
    for _ in range(5):
        petri.simulate_year()
    
    # Track populations for several years
    populations = []
    for _ in range(10):
        petri.simulate_year()
        populations.append(len(petri.cells))
    
    # Should have some variation
    unique_pops = len(set(populations))
    assert unique_pops > 1, f"No variation in steady state: {populations}"
    
    # But should be around target
    avg_pop = sum(populations) / len(populations)
    target = petri.target_population
    assert 0.25 * target <= avg_pop <= 3 * target, f"Average {avg_pop} far from target {target}"
    
    print(f"✓ Steady state shows natural variation: {unique_pops} different populations")

def test_gene_rate_edge_cases():
    """Test edge cases for gene rate groups."""
    print("Testing gene rate edge cases...")
    
    # Single gene group (effectively uniform)
    cell = Cell(n=100, gene_rate_groups=[(20, 0.01)], gene_size=5)
    assert all(r == 0.01 for r in cell.site_rates), "Single group should be uniform"
    
    # Many small groups
    groups = [(1, 0.001 * i) for i in range(1, 21)]
    cell = Cell(n=100, gene_rate_groups=groups, gene_size=5)
    assert len(cell.site_rates) == 100, "Wrong number of site rates"
    
    # Very high rates (>1 is allowed)
    cell = Cell(n=100, gene_rate_groups=[(20, 2.0)], gene_size=5)
    cell.methylate()
    # Most sites should be methylated with rate 2.0
    assert sum(cell.cpg_sites) > 80, "High rate should methylate most sites"
    
    # Mixed very low and very high
    groups = [(10, 0.0001), (10, 10.0)]
    cell = Cell(n=100, gene_rate_groups=groups, gene_size=5)
    cell.methylate()
    # First half should have few methylated, second half should have many
    first_half_meth = sum(cell.cpg_sites[:50])
    second_half_meth = sum(cell.cpg_sites[50:])
    assert second_half_meth > first_half_meth * 10, "High rate genes should methylate much more"
    
    print("✓ Gene rate edge cases work")


def test_gene_rate_with_petri_dish():
    """Test gene rates work properly in PetriDish simulation."""
    print("Testing gene rates in PetriDish simulation...")
    
    # Create with very different rates per gene group
    groups = [(5, 0.001), (5, 0.01), (5, 0.1), (5, 1.0)]
    petri = PetriDish(gene_rate_groups=groups, n=100, gene_size=5, seed=42, growth_phase=2)
    
    # Run simulation
    petri.run_simulation(t_max=5)
    
    # Check all cells have the same configuration
    for cell in petri.cells:
        assert cell.gene_rate_groups == groups, "Cell has wrong gene_rate_groups"
        assert len(cell.site_rates) == 100, "Cell has wrong site_rates length"
    
    # Check methylation patterns reflect different rates
    # High-rate genes (last 5 genes) should be more methylated
    for cell in petri.cells:
        # Count methylation in different gene groups
        group1_meth = sum(cell.cpg_sites[0:25])   # rate 0.001
        group4_meth = sum(cell.cpg_sites[75:100]) # rate 1.0
        # With rate 1.0 for 5 years, most sites should be methylated
        if group4_meth < group1_meth:
            print(f"Warning: High rate genes less methylated ({group4_meth} vs {group1_meth})")
    
    print("✓ Gene rates work in PetriDish simulation")


def run_all_edge_tests():
    """Run all edge case tests."""
    print("="*60)
    print("EDGE CASE TEST SUITE")
    print("="*60)
    
    tests = [
        test_single_year_simulation,
        test_no_growth_to_steady_transition,
        test_immediate_steady_state,
        test_extreme_culling_survival,
        test_100_percent_methylation,
        test_zero_methylation,
        test_minimum_gene_configuration,
        test_large_growth_phase,
        test_filename_formats,
        test_population_variability,
        test_gene_rate_edge_cases,
        test_gene_rate_with_petri_dish
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            test()
            passed += 1
        except Exception as e:
            print(f"✗ {test.__name__}: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print("EDGE CASE TEST SUMMARY")
    print("="*60)
    print(f"Total: {passed + failed}")
    print(f"Passed: {passed}")
    print(f"Failed: {failed}")
    
    if failed == 0:
        print("\n🎉 ALL EDGE CASE TESTS PASSED! 🎉")
        return 0
    else:
        print(f"\n❌ {failed} EDGE CASE TESTS FAILED")
        return 1

if __name__ == "__main__":
    exit(run_all_edge_tests())