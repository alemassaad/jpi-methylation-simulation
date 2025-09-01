#!/usr/bin/env python3
"""
Test gene JSD analysis generation.
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
from pipeline_analysis import generate_gene_jsd_analysis


def test_gene_jsd_analysis():
    """Test that gene_jsd_analysis.json is created with correct structure."""
    print("\nTest: Gene JSD Analysis Generation")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create mock PetriDish objects with known gene structure
        print("Creating mock PetriDish objects...")
        
        n_sites = 100
        gene_size = 5
        n_genes = n_sites // gene_size  # 20 genes
        
        # Create mutant dishes
        mutant_dishes = []
        for i in range(3):
            petri = PetriDish(rate=0.005, n=n_sites, gene_size=gene_size)
            # Add some cells
            for j in range(10):
                cell = Cell(n=n_sites, rate=0.005, gene_size=gene_size)
                # Apply some methylation for testing
                cell.methylate()  # Apply random methylation
                petri.cells.append(cell)
            mutant_dishes.append(petri)
        
        # Create control1 dishes
        control1_dishes = []
        for i in range(3):
            petri = PetriDish(rate=0.005, n=n_sites, gene_size=gene_size)
            for j in range(10):
                cell = Cell(n=n_sites, rate=0.005, gene_size=gene_size)
                # Apply methylation
                cell.methylate()
                petri.cells.append(cell)
            control1_dishes.append(petri)
        
        # Create control2 dishes
        control2_dishes = []
        for i in range(3):
            petri = PetriDish(rate=0.005, n=n_sites, gene_size=gene_size)
            for j in range(10):
                cell = Cell(n=n_sites, rate=0.005, gene_size=gene_size)
                # Apply methylation
                cell.methylate()
                petri.cells.append(cell)
            control2_dishes.append(petri)
        
        # Run analysis
        print("Running gene JSD analysis...")
        results_dir = os.path.join(tmpdir, "results")
        results = generate_gene_jsd_analysis(
            mutant_dishes, control1_dishes, control2_dishes, results_dir
        )
        
        # Check that file exists
        gene_jsd_path = os.path.join(results_dir, "gene_jsd_analysis.json")
        assert os.path.exists(gene_jsd_path), "gene_jsd_analysis.json not created"
        print("  ✓ gene_jsd_analysis.json created")
        
        # Load and verify structure
        with open(gene_jsd_path, 'r') as f:
            gene_jsd_data = json.load(f)
        
        # Check top-level keys
        assert "gene_jsd_distributions" in gene_jsd_data, "Missing gene_jsd_distributions"
        assert "gene_metadata" in gene_jsd_data, "Missing gene_metadata"
        print("  ✓ Top-level keys present")
        
        # Check gene_jsd_distributions structure
        distributions = gene_jsd_data["gene_jsd_distributions"]
        assert len(distributions) == n_genes, f"Expected {n_genes} genes, got {len(distributions)}"
        
        for gene_idx in range(n_genes):
            gene_key = f"gene_{gene_idx}"
            assert gene_key in distributions, f"Missing {gene_key}"
            
            gene_data = distributions[gene_key]
            assert "mutant" in gene_data, f"Missing mutant data for {gene_key}"
            assert "control1" in gene_data, f"Missing control1 data for {gene_key}"
            assert "control2" in gene_data, f"Missing control2 data for {gene_key}"
            
            # Check that each batch has correct number of values
            assert len(gene_data["mutant"]) == 3, f"Wrong number of mutant values for {gene_key}"
            assert len(gene_data["control1"]) == 3, f"Wrong number of control1 values for {gene_key}"
            assert len(gene_data["control2"]) == 3, f"Wrong number of control2 values for {gene_key}"
            
            # Check that values are floats and in valid range
            for batch in ["mutant", "control1", "control2"]:
                for val in gene_data[batch]:
                    assert isinstance(val, float), f"Non-float value in {batch} for {gene_key}"
                    assert 0 <= val <= 1, f"Invalid JSD value {val} in {batch} for {gene_key}"
        
        print(f"  ✓ All {n_genes} genes have correct structure")
        
        # Check metadata
        metadata = gene_jsd_data["gene_metadata"]
        assert metadata["n_genes"] == n_genes, f"Wrong n_genes in metadata"
        assert metadata["n_individuals_per_batch"]["mutant"] == 3
        assert metadata["n_individuals_per_batch"]["control1"] == 3
        assert metadata["n_individuals_per_batch"]["control2"] == 3
        print("  ✓ Metadata correct")
        
        print("\n✅ Gene JSD analysis test passed!")
        return True


def test_with_gene_rates():
    """Test gene JSD analysis with gene-specific rates."""
    print("\nTest: Gene JSD Analysis with Gene-Specific Rates")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create dishes with gene-specific rates
        gene_rate_groups = [(10, 0.004), (10, 0.006)]  # 20 genes total
        
        mutant_dishes = []
        for i in range(2):
            petri = PetriDish(gene_rate_groups=gene_rate_groups, n=100, gene_size=5)
            for j in range(5):
                cell = Cell(n=100, gene_rate_groups=gene_rate_groups, gene_size=5)
                petri.cells.append(cell)
            mutant_dishes.append(petri)
        
        control1_dishes = mutant_dishes  # Reuse for simplicity
        control2_dishes = mutant_dishes
        
        # Run analysis
        results_dir = os.path.join(tmpdir, "results")
        results = generate_gene_jsd_analysis(
            mutant_dishes, control1_dishes, control2_dishes, results_dir
        )
        
        # Load and check for gene_rate_groups in metadata
        with open(os.path.join(results_dir, "gene_jsd_analysis.json"), 'r') as f:
            gene_jsd_data = json.load(f)
        
        metadata = gene_jsd_data["gene_metadata"]
        assert "gene_rate_groups" in metadata, "Missing gene_rate_groups in metadata"
        
        rate_groups = metadata["gene_rate_groups"]
        assert len(rate_groups) == 2, "Wrong number of rate groups"
        assert rate_groups[0]["start"] == 0
        assert rate_groups[0]["end"] == 9
        assert rate_groups[0]["rate"] == 0.004
        assert rate_groups[1]["start"] == 10
        assert rate_groups[1]["end"] == 19
        assert rate_groups[1]["rate"] == 0.006
        
        print("  ✓ Gene rate groups correctly recorded in metadata")
        print("\n✅ Gene rate test passed!")
        return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing Gene JSD Analysis Generation")
    print("=" * 60)
    
    tests = [
        test_gene_jsd_analysis,
        test_with_gene_rates
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
        print("\n🎉 All gene JSD analysis tests passed!")
        print("   ✅ gene_jsd_analysis.json created successfully")
        print("   ✅ Structure matches specification")
        print("   ✅ Gene-specific rates handled correctly")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())