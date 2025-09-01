#!/usr/bin/env python3
"""
Test that gene_mean_methylation returns proportions (0.0 to 1.0) not counts.
"""

import os
import sys
import json
import tempfile

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish


def test_proportions_not_counts():
    """Test that gene_mean_methylation returns proportions."""
    print("\nTest 1: Proportions Not Counts")
    print("=" * 50)
    
    # Create PetriDish with known methylation
    petri = PetriDish(n=20, rate=0.005, growth_phase=2)
    petri.n_genes = 4
    petri.gene_size = 5
    
    # Create cells with specific methylation patterns
    # Gene 0: 0/5, 2/5, 5/5 methylated
    cell1 = Cell(n=20, rate=0.005)
    cell2 = Cell(n=20, rate=0.005)
    cell3 = Cell(n=20, rate=0.005)
    
    # Set methylation patterns
    cell1.cpg_sites = [0,0,0,0,0] + [0]*15  # Gene 0: 0/5 = 0.0
    cell2.cpg_sites = [1,1,0,0,0] + [0]*15  # Gene 0: 2/5 = 0.4
    cell3.cpg_sites = [1,1,1,1,1] + [0]*15  # Gene 0: 5/5 = 1.0
    
    petri.cells = [cell1, cell2, cell3]
    
    # Calculate gene means
    gene_means = petri.calculate_gene_mean_methylation()
    
    # Check Gene 0
    # Mean count: (0 + 2 + 5) / 3 = 2.33
    # Mean proportion: 2.33 / 5 = 0.467
    expected_proportion = (0 + 2 + 5) / 3 / 5
    
    print(f"  Gene 0 methylation patterns: 0/5, 2/5, 5/5")
    print(f"  Mean count: {(0 + 2 + 5) / 3:.2f} sites")
    print(f"  Mean proportion: {gene_means[0]:.3f}")
    print(f"  Expected: {expected_proportion:.3f}")
    
    assert abs(gene_means[0] - expected_proportion) < 0.001, \
        f"Gene 0 should be {expected_proportion:.3f}, got {gene_means[0]:.3f}"
    
    # All values should be between 0 and 1
    for i, mean in enumerate(gene_means):
        assert 0 <= mean <= 1.0, f"Gene {i} proportion {mean} out of range [0, 1]"
    
    print("  ✅ Returns proportions (0.0 to 1.0)")
    return True


def test_edge_cases():
    """Test edge cases: fully methylated and unmethylated."""
    print("\nTest 2: Edge Cases")
    print("=" * 50)
    
    petri = PetriDish(n=20, rate=0.005, growth_phase=2)
    
    # Test 1: Completely unmethylated
    cell = Cell(n=20, rate=0.005)
    cell.cpg_sites = [0] * 20  # All unmethylated
    petri.cells = [cell]
    
    gene_means = petri.calculate_gene_mean_methylation()
    
    assert all(mean == 0.0 for mean in gene_means), "Unmethylated genes should be 0.0"
    print("  ✅ Unmethylated genes = 0.0")
    
    # Test 2: Fully methylated
    cell.cpg_sites = [1] * 20  # All methylated
    petri.cells = [cell]
    
    gene_means = petri.calculate_gene_mean_methylation()
    
    assert all(mean == 1.0 for mean in gene_means), "Fully methylated genes should be 1.0"
    print("  ✅ Fully methylated genes = 1.0")
    
    # Test 3: Half methylated (uniform)
    cell.cpg_sites = [1,1,0,0,0] * 4  # 2/5 = 0.4 for each gene
    petri.cells = [cell]
    
    gene_means = petri.calculate_gene_mean_methylation()
    
    assert all(abs(mean - 0.4) < 0.001 for mean in gene_means), \
        f"Half methylated genes should be 0.4, got {gene_means}"
    print("  ✅ Partially methylated genes = 0.4 (2/5)")
    
    return True


def test_saved_as_proportions():
    """Test that saved JSON contains proportions."""
    print("\nTest 3: Saved as Proportions")
    print("=" * 50)
    
    # Create PetriDish with some methylation
    petri = PetriDish(n=100, rate=0.005, growth_phase=3)
    
    # Add some methylation
    for _ in range(5):
        petri.methylate_cells()
        if petri.year < 3:
            petri.divide_cells()
    
    # Save with gene metrics
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_gene_metrics=True, compress=False)
        
        # Load and check
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        gene_means = data['metadata']['gene_mean_methylation']
        
        # All values should be proportions (0 to 1)
        assert all(0 <= mean <= 1.0 for mean in gene_means), \
            f"Some values out of proportion range: {[m for m in gene_means if m > 1.0]}"
        
        # Print some examples
        print(f"  Sample gene proportions:")
        for i in range(min(5, len(gene_means))):
            percentage = gene_means[i] * 100
            print(f"    Gene {i}: {gene_means[i]:.3f} ({percentage:.1f}%)")
        
        print(f"  ✅ All {len(gene_means)} genes saved as proportions")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def test_interpretation():
    """Test that proportions are easy to interpret."""
    print("\nTest 4: Interpretation")
    print("=" * 50)
    
    petri = PetriDish(n=20, rate=0.005, growth_phase=2)
    
    # Create cells with interpretable patterns
    cells = []
    
    # 10 cells with Gene 0 at 20% methylated (1/5)
    for _ in range(10):
        cell = Cell(n=20, rate=0.005)
        cell.cpg_sites = [1,0,0,0,0] + [0]*15
        cells.append(cell)
    
    # 10 cells with Gene 0 at 80% methylated (4/5)
    for _ in range(10):
        cell = Cell(n=20, rate=0.005)
        cell.cpg_sites = [1,1,1,1,0] + [0]*15
        cells.append(cell)
    
    petri.cells = cells
    gene_means = petri.calculate_gene_mean_methylation()
    
    # Gene 0 should average to 50% (0.5)
    # (10 cells * 0.2 + 10 cells * 0.8) / 20 = 0.5
    expected = 0.5
    
    print(f"  Population: 10 cells at 20% + 10 cells at 80%")
    print(f"  Gene 0 mean: {gene_means[0]:.1%}")  # Format as percentage
    print(f"  Expected: {expected:.1%}")
    
    assert abs(gene_means[0] - expected) < 0.001, \
        f"Should be {expected}, got {gene_means[0]}"
    
    print("  ✅ Proportions are intuitive (0.5 = 50%)")
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Gene Mean Methylation Proportion Tests")
    print("=" * 60)
    
    tests = [
        test_proportions_not_counts,
        test_edge_cases,
        test_saved_as_proportions,
        test_interpretation
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"\n  ❌ Test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Test Summary: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("\n🎉 Gene Methylation Proportions Working!")
        print("   ✅ Returns proportions (0.0 to 1.0)")
        print("   ✅ Not counts (0 to 5)")
        print("   ✅ Easy to interpret as percentages")
        print("   ✅ Consistent with other proportion metrics")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())