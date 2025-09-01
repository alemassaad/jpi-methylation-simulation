#!/usr/bin/env python3
"""
Comprehensive test suite for gene-level JSD and mean methylation metrics.
Tests calculation, storage, and consistency of gene metrics.
"""

import os
import sys
import json
import tempfile
import numpy as np

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, load_petri_dish


def test_gene_metrics_calculation():
    """Test that gene metrics are calculated correctly."""
    print("\nTest 1: Gene Metrics Calculation")
    print("=" * 50)
    
    # Create a PetriDish with known methylation patterns
    petri = PetriDish(rate=0.005, growth_phase=2, calculate_cell_jsds=True)
    
    # Manually create cells with specific methylation patterns
    # 20 sites total (4 genes × 5 sites/gene)
    petri.n = 20
    petri.n_genes = 4
    petri.gene_size = 5
    
    # Create 3 cells with different methylation patterns
    cell1 = Cell(n=20, rate=0.005)
    cell2 = Cell(n=20, rate=0.005)
    cell3 = Cell(n=20, rate=0.005)
    
    # Gene 0 (sites 0-4): Variable methylation
    cell1.cpg_sites = [0, 0, 0, 0, 0] + [0]*15  # 0/5 methylated
    cell2.cpg_sites = [1, 1, 0, 0, 0] + [0]*15  # 2/5 methylated
    cell3.cpg_sites = [1, 1, 1, 1, 1] + [0]*15  # 5/5 methylated
    
    petri.cells = [cell1, cell2, cell3]
    
    # Calculate gene metrics
    gene_jsds = petri.calculate_gene_jsds()
    gene_means = petri.calculate_gene_mean_methylation()
    
    print(f"  Gene JSDs: {[f'{jsd:.4f}' for jsd in gene_jsds]}")
    print(f"  Gene means: {[f'{mean:.2f}' for mean in gene_means]}")
    
    # Check Gene 0 metrics
    # Mean methylation: (0 + 2 + 5) / 3 = 2.33
    expected_mean_gene0 = 7/3
    assert abs(gene_means[0] - expected_mean_gene0) < 0.01, f"Gene 0 mean incorrect: {gene_means[0]} vs {expected_mean_gene0}"
    
    # Gene 0 should have high JSD (heterogeneous)
    assert gene_jsds[0] > 0.1, f"Gene 0 JSD should be high (heterogeneous): {gene_jsds[0]}"
    
    # Genes 1-3 should have JSD = 0 (all unmethylated)
    for i in range(1, 4):
        assert gene_jsds[i] < 0.001, f"Gene {i} JSD should be ~0 (homogeneous): {gene_jsds[i]}"
        assert gene_means[i] < 0.001, f"Gene {i} mean should be 0: {gene_means[i]}"
    
    print("  ✅ Gene metrics calculated correctly")
    return True


def test_gene_metrics_with_uniform_pattern():
    """Test gene metrics with uniform methylation across genes."""
    print("\nTest 2: Uniform Methylation Pattern")
    print("=" * 50)
    
    # Create PetriDish with n=20 (4 genes)
    petri = PetriDish(n=20, rate=0.005, growth_phase=2)
    petri.cells = []  # Clear initial cell
    
    # Create cells with uniform 40% methylation
    for _ in range(10):
        cell = Cell(n=20, rate=0.005)
        # Each gene has 2/5 sites methylated
        cell.cpg_sites = [1, 1, 0, 0, 0] * 4
        petri.cells.append(cell)
    
    gene_jsds = petri.calculate_gene_jsds()
    gene_means = petri.calculate_gene_mean_methylation()
    
    print(f"  Gene JSDs: {[f'{jsd:.4f}' for jsd in gene_jsds]}")
    print(f"  Gene means: {[f'{mean:.2f}' for mean in gene_means]}")
    
    # All genes should have same JSD (all cells identical)
    for i in range(4):
        assert abs(gene_means[i] - 2.0) < 0.001, f"Gene {i} mean should be 2.0"
        # JSD should be > 0 because different from baseline (all unmethylated)
        assert gene_jsds[i] > 0.1, f"Gene {i} JSD should be > 0"
    
    # All genes should have similar JSD values
    jsd_std = np.std(gene_jsds)
    assert jsd_std < 0.001, f"JSDs should be uniform: std={jsd_std}"
    
    print("  ✅ Uniform pattern metrics correct")
    return True


def test_gene_metrics_in_saved_file():
    """Test that gene metrics are saved correctly in JSON."""
    print("\nTest 3: Gene Metrics in Saved Files")
    print("=" * 50)
    
    # Create PetriDish with some cells
    petri = PetriDish(rate=0.005, growth_phase=3)
    for _ in range(3):
        petri.divide_cells()
        petri.methylate_cells()
    
    # Save with gene metrics
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_gene_metrics=True, compress=False)
        
        # Load JSON directly
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        metadata = data['metadata']
        
        # Check that gene metrics are present
        assert 'gene_jsds' in metadata, "gene_jsds missing from metadata"
        assert 'gene_mean_methylation' in metadata, "gene_mean_methylation missing"
        assert 'n_genes' in metadata, "n_genes missing from metadata"
        
        # Check dimensions
        n_genes = metadata['n_genes']
        assert len(metadata['gene_jsds']) == n_genes, f"Wrong number of gene JSDs"
        assert len(metadata['gene_mean_methylation']) == n_genes, f"Wrong number of gene means"
        
        # Check value ranges
        for jsd in metadata['gene_jsds']:
            assert 0 <= jsd <= 1, f"JSD out of range: {jsd}"
        
        for mean in metadata['gene_mean_methylation']:
            assert 0 <= mean <= 5, f"Mean methylation out of range: {mean}"
        
        print(f"  ✅ Saved {n_genes} genes with metrics")
        print(f"     Mean gene JSD: {np.mean(metadata['gene_jsds']):.4f}")
        print(f"     Mean gene methylation: {np.mean(metadata['gene_mean_methylation']):.4f}")
        
        return True
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)


def test_gene_metrics_consistency():
    """Test that gene metrics are consistent across saves/loads."""
    print("\nTest 4: Gene Metrics Consistency")
    print("=" * 50)
    
    # Create PetriDish with known seed
    np.random.seed(42)
    petri = PetriDish(rate=0.005, growth_phase=4)
    
    # Grow to get some methylation
    for _ in range(5):
        petri.divide_cells()
        petri.methylate_cells()
        if petri.year > petri.growth_phase:
            petri.random_cull_cells()
    
    # Calculate metrics directly
    direct_jsds = petri.calculate_gene_jsds()
    direct_means = petri.calculate_gene_mean_methylation()
    
    # Save and reload
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_gene_metrics=True, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        saved_jsds = data['metadata']['gene_jsds']
        saved_means = data['metadata']['gene_mean_methylation']
        
        # Compare
        for i in range(len(direct_jsds)):
            assert abs(direct_jsds[i] - saved_jsds[i]) < 0.0001, \
                f"Gene {i} JSD mismatch: {direct_jsds[i]} vs {saved_jsds[i]}"
            assert abs(direct_means[i] - saved_means[i]) < 0.0001, \
                f"Gene {i} mean mismatch: {direct_means[i]} vs {saved_means[i]}"
        
        print("  ✅ Metrics consistent between calculation and save")
        return True
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)


def test_gene_metrics_with_gene_rates():
    """Test gene metrics with gene-specific methylation rates."""
    print("\nTest 5: Gene-Specific Methylation Rates")
    print("=" * 50)
    
    # Create PetriDish with different rates per gene
    # Need to specify groups that sum to actual gene count
    # For n=20, gene_size=5 -> 4 genes total
    gene_rate_groups = [(1, 0.01), (1, 0.001), (1, 0.005), (1, 0.005)]
    
    petri = PetriDish(n=20, gene_rate_groups=gene_rate_groups, growth_phase=3)
    
    # Simulate to accumulate methylation
    for _ in range(10):
        for cell in petri.cells:
            cell.methylate()
        if petri.year < petri.growth_phase:
            petri.divide_cells()
    
    gene_jsds = petri.calculate_gene_jsds()
    gene_means = petri.calculate_gene_mean_methylation()
    
    print(f"  Gene rates: {[g[1] for g in gene_rate_groups]}")
    print(f"  Gene JSDs: {[f'{jsd:.4f}' for jsd in gene_jsds]}")
    print(f"  Gene means: {[f'{mean:.3f}' for mean in gene_means]}")
    
    # Gene 0 (high rate) should have higher mean methylation than Gene 1 (low rate)
    assert gene_means[0] > gene_means[1], \
        f"High-rate gene should have more methylation: {gene_means[0]} vs {gene_means[1]}"
    
    print("  ✅ Gene-specific rates affect metrics correctly")
    return True


def test_edge_cases():
    """Test edge cases: empty dish, single cell, etc."""
    print("\nTest 6: Edge Cases")
    print("=" * 50)
    
    # Empty PetriDish
    empty_petri = PetriDish(rate=0.005, growth_phase=1)
    empty_petri.cells = []
    
    empty_jsds = empty_petri.calculate_gene_jsds()
    empty_means = empty_petri.calculate_gene_mean_methylation()
    
    # Empty dish returns zeros for each gene, not empty list
    assert len(empty_jsds) == empty_petri.n_genes, "Empty dish should return JSDs for all genes"
    assert all(jsd == 0.0 for jsd in empty_jsds), "Empty dish JSDs should all be 0"
    assert empty_means == [], "Empty dish should return empty gene means"
    print("  ✅ Empty dish handled correctly")
    
    # Single cell
    single_petri = PetriDish(n=100, rate=0.005, growth_phase=1)
    # Already has one cell with correct size
    
    single_jsds = single_petri.calculate_gene_jsds()
    single_means = single_petri.calculate_gene_mean_methylation()
    
    assert len(single_jsds) == 20, "Single cell should still return 20 gene JSDs"
    assert len(single_means) == 20, "Single cell should still return 20 gene means"
    
    # Single cell means no heterogeneity, but JSD compares to baseline
    print(f"  Single cell JSDs: min={min(single_jsds):.4f}, max={max(single_jsds):.4f}")
    print("  ✅ Single cell handled correctly")
    
    return True


def test_gene_metrics_without_flag():
    """Test that gene metrics are NOT saved when flag is False."""
    print("\nTest 7: Gene Metrics Flag Control")
    print("=" * 50)
    
    petri = PetriDish(rate=0.005, growth_phase=2)
    petri.divide_cells()
    petri.methylate_cells()
    
    # Save WITHOUT gene metrics
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, include_gene_metrics=False, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        metadata = data['metadata']
        
        # Check that gene metrics are NOT present
        assert 'gene_jsds' not in metadata, "gene_jsds should not be saved when flag is False"
        assert 'gene_mean_methylation' not in metadata, "gene_mean_methylation should not be saved"
        
        print("  ✅ Gene metrics correctly excluded when flag is False")
        return True
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)


def main():
    """Run all gene metrics tests."""
    print("="*60)
    print("Testing Gene-Level Metrics Implementation")
    print("="*60)
    
    tests = [
        test_gene_metrics_calculation,
        test_gene_metrics_with_uniform_pattern,
        test_gene_metrics_in_saved_file,
        test_gene_metrics_consistency,
        test_gene_metrics_with_gene_rates,
        test_edge_cases,
        test_gene_metrics_without_flag
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"  ❌ Test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    print(f"Test Summary: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All gene metrics tests passed!")
        return 0
    else:
        print(f"⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())