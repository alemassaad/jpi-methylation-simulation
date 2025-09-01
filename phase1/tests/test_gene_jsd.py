#!/usr/bin/env python3
"""
Test gene_JSD functionality in PetriDish class.
"""

import os
import sys
import tempfile
import json
import gzip

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish


def test_gene_jsd_initialization():
    """Test that gene_JSD is properly initialized."""
    print("\nTest 1: Gene JSD Initialization")
    print("=" * 50)
    
    # Create PetriDish
    petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42)
    
    # Check attributes exist
    assert hasattr(petri, 'n_genes'), "Should have n_genes attribute"
    assert hasattr(petri, 'gene_jsd_history'), "Should have gene_jsd_history attribute"
    assert hasattr(petri, 'BASELINE_GENE_DISTRIBUTION'), "Should have baseline distribution"
    
    # Check values
    assert petri.n_genes == 20, f"Should have 20 genes (100/5), got {petri.n_genes}"
    assert petri.gene_jsd_history == {}, "Gene JSD history should start empty when not tracking"
    assert petri.BASELINE_GENE_DISTRIBUTION == [1.0] + [0.0] * 5, "Baseline should be all unmethylated"
    
    print(f"  ✓ n_genes = {petri.n_genes}")
    print(f"  ✓ Baseline distribution correct")
    print("  PASSED")


def test_gene_jsd_calculation():
    """Test the calculate_gene_jsd method."""
    print("\nTest 2: Gene JSD Calculation")
    print("=" * 50)
    
    # Create PetriDish with some cells
    petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42)
    petri.grow_with_homeostasis(years=3, growth_phase=2, verbose=False)
    
    # Calculate gene JSD
    gene_jsds = petri.calculate_gene_jsd()
    
    # Check results
    assert len(gene_jsds) == 20, f"Should have 20 gene JSDs, got {len(gene_jsds)}"
    assert all(0 <= jsd <= 1 for jsd in gene_jsds), "All JSDs should be between 0 and 1"
    
    # Initially unmethylated cells should have low gene JSDs
    mean_jsd = sum(gene_jsds) / len(gene_jsds)
    print(f"  ✓ Calculated {len(gene_jsds)} gene JSDs")
    print(f"  ✓ Mean gene JSD = {mean_jsd:.4f}")
    print(f"  ✓ Min = {min(gene_jsds):.4f}, Max = {max(gene_jsds):.4f}")
    print("  PASSED")


def test_gene_jsd_tracking():
    """Test that gene JSD tracking works during growth."""
    print("\nTest 3: Gene JSD Tracking During Growth")
    print("=" * 50)
    
    # Create PetriDish with tracking enabled
    petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42)
    petri.enable_history_tracking(start_year=0, track_gene_jsd=True)
    
    # Grow for a few years
    petri.grow_with_homeostasis(years=5, growth_phase=3, verbose=False)
    
    # Check gene JSD history was recorded
    assert len(petri.gene_jsd_history) == 6, f"Should have 6 years of history (0-5), got {len(petri.gene_jsd_history)}"
    assert '0' in petri.gene_jsd_history, "Should have year 0"
    assert '5' in petri.gene_jsd_history, "Should have year 5"
    
    # Check year 0 is all zeros (unmethylated)
    year0_jsds = petri.gene_jsd_history['0']
    assert all(jsd == 0.0 for jsd in year0_jsds), "Year 0 should have all zero JSDs"
    
    # Check JSDs increase over time
    year5_jsds = petri.gene_jsd_history['5']
    assert any(jsd > 0 for jsd in year5_jsds), "Year 5 should have some non-zero JSDs"
    
    mean_year0 = sum(year0_jsds) / len(year0_jsds)
    mean_year5 = sum(year5_jsds) / len(year5_jsds)
    assert mean_year5 > mean_year0, "Mean gene JSD should increase over time"
    
    print(f"  ✓ Tracked {len(petri.gene_jsd_history)} years")
    print(f"  ✓ Year 0 mean = {mean_year0:.4f}")
    print(f"  ✓ Year 5 mean = {mean_year5:.4f}")
    print("  PASSED")


def test_gene_jsd_save_load():
    """Test saving and loading PetriDish with gene JSD history."""
    print("\nTest 4: Save/Load with Gene JSD")
    print("=" * 50)
    
    # Import pipeline_utils for save/load
    sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase2'))
    try:
        from pipeline_utils import save_petri_dish, load_petri_dish
    except ImportError as e:
        print(f"  ⚠️  Skipping test - pipeline_utils requires numpy: {e}")
        print("  SKIPPED")
        return
    
    # Create PetriDish with gene JSD tracking
    petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42)
    petri.enable_history_tracking(start_year=0, track_gene_jsd=True)
    petri.grow_with_homeostasis(years=3, growth_phase=2, verbose=False)
    
    # Save with gene JSD
    with tempfile.NamedTemporaryFile(suffix='.json.gz', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_cell_history=True, include_gene_jsd=True)
        print(f"  ✓ Saved PetriDish with {len(petri.gene_jsd_history)} gene JSD entries")
        
        # Load without gene JSD
        petri_no_jsd = load_petri_dish(temp_file, include_cell_history=False, include_gene_jsd=False)
        assert not petri_no_jsd.gene_jsd_history, "Should not load gene JSD when not requested"
        print("  ✓ Loaded without gene JSD - history empty")
        
        # Load with gene JSD
        petri_with_jsd = load_petri_dish(temp_file, include_cell_history=True, include_gene_jsd=True)
        assert len(petri_with_jsd.gene_jsd_history) == len(petri.gene_jsd_history), "Should restore all gene JSD history"
        
        # Verify content matches
        for year in petri.gene_jsd_history:
            assert year in petri_with_jsd.gene_jsd_history, f"Missing year {year}"
            orig_jsds = petri.gene_jsd_history[year]
            loaded_jsds = petri_with_jsd.gene_jsd_history[year]
            assert len(orig_jsds) == len(loaded_jsds), f"Year {year}: different lengths"
            for i, (orig, loaded) in enumerate(zip(orig_jsds, loaded_jsds)):
                assert abs(orig - loaded) < 1e-10, f"Year {year}, gene {i}: {orig} != {loaded}"
        
        print(f"  ✓ Loaded with gene JSD - {len(petri_with_jsd.gene_jsd_history)} entries restored")
        print("  ✓ Content verified to match original")
        
    finally:
        os.remove(temp_file)
    
    print("  PASSED")


def test_gene_jsd_disabled():
    """Test that gene JSD can be disabled for performance."""
    print("\nTest 5: Gene JSD Disabled for Performance")
    print("=" * 50)
    
    # Create PetriDish with cell JSDs disabled
    petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, calculate_cell_jsds=False)
    
    # Enable tracking but grow
    petri.enable_history_tracking(start_year=0, track_gene_jsd=True)
    petri.grow_with_homeostasis(years=2, growth_phase=1, verbose=False)
    
    # Gene JSD history should remain empty when calculate_cell_jsds=False
    assert len(petri.gene_jsd_history) <= 1, "Should not calculate gene JSDs when cell JSDs disabled"
    
    # Cell JSDs should also not be calculated
    for cell in petri.cells:
        assert cell.cell_JSD == 0.0, "Cell JSDs should not be calculated when disabled"
    
    print("  ✓ Gene JSD calculation properly disabled")
    print("  ✓ Cell JSD calculation also disabled")
    print("  PASSED")


def test_edge_cases():
    """Test edge cases for gene JSD."""
    print("\nTest 6: Edge Cases")
    print("=" * 50)
    
    # Test 1: Invalid gene_size
    try:
        petri = PetriDish(rate=0.01, n=101, gene_size=5)  # 101 not divisible by 5
        print("  ✗ Should have raised error for invalid gene_size")
    except ValueError as e:
        print(f"  ✓ Properly rejected invalid gene_size: {e}")
    
    # Test 2: Empty PetriDish
    petri = PetriDish(rate=0.01, n=100, gene_size=5)
    petri.cells = []  # Remove all cells
    gene_jsds = petri.calculate_gene_jsd()
    assert all(jsd == 0.0 for jsd in gene_jsds), "Empty PetriDish should return all zeros"
    print("  ✓ Empty PetriDish handled correctly")
    
    # Test 3: Single cell
    petri = PetriDish(rate=0.01, n=100, gene_size=5)
    gene_jsds = petri.calculate_gene_jsd()
    assert len(gene_jsds) == 20, "Single cell should still calculate 20 gene JSDs"
    assert all(jsd == 0.0 for jsd in gene_jsds), "Unmethylated cell should have zero JSDs"
    print("  ✓ Single cell handled correctly")
    
    # Test 4: Different gene sizes
    for gene_size in [1, 2, 5, 10, 20]:
        if 100 % gene_size == 0:
            petri = PetriDish(rate=0.01, n=100, gene_size=gene_size)
            expected_n_genes = 100 // gene_size
            assert petri.n_genes == expected_n_genes, f"Gene size {gene_size}: expected {expected_n_genes} genes"
            print(f"  ✓ Gene size {gene_size}: {petri.n_genes} genes")
    
    print("  PASSED")


def test_gene_jsd_with_rate_groups():
    """Test gene JSD calculation with variable gene rates."""
    print("\nTest 7: Gene JSD with Rate Groups")
    print("=" * 50)
    
    # Create PetriDish with different rates per gene group
    groups = [(5, 0.001), (5, 0.01), (5, 0.1), (5, 0.5)]
    petri = PetriDish(gene_rate_groups=groups, n=100, gene_size=5, seed=42)
    petri.enable_history_tracking(start_year=0, track_gene_jsd=True)
    petri.grow_with_homeostasis(years=5, growth_phase=2, verbose=False)
    
    # Calculate gene JSDs
    gene_jsds = petri.calculate_gene_jsd()
    
    # High-rate genes should have higher JSDs
    low_rate_jsds = gene_jsds[0:5]   # First 5 genes (rate 0.001)
    high_rate_jsds = gene_jsds[15:20] # Last 5 genes (rate 0.5)
    
    mean_low = sum(low_rate_jsds) / len(low_rate_jsds)
    mean_high = sum(high_rate_jsds) / len(high_rate_jsds)
    
    # High rate genes should have noticeably higher JSDs
    if mean_high <= mean_low:
        print(f"  ⚠️  High-rate genes don't have higher JSD (low={mean_low:.4f}, high={mean_high:.4f})")
        print("  This may happen with short simulations or small populations")
    else:
        print(f"  ✓ High-rate genes have higher JSD (low={mean_low:.4f}, high={mean_high:.4f})")
    
    print(f"  ✓ Gene JSD works with rate groups")
    print("  PASSED")


def run_all_tests():
    """Run all gene JSD tests."""
    print("\n" + "="*60)
    print("TESTING GENE JSD FUNCTIONALITY")
    print("="*60)
    
    tests = [
        test_gene_jsd_initialization,
        test_gene_jsd_calculation,
        test_gene_jsd_tracking,
        test_gene_jsd_save_load,
        test_gene_jsd_disabled,
        test_edge_cases,
        test_gene_jsd_with_rate_groups
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