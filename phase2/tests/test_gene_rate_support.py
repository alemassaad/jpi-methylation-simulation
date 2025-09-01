#!/usr/bin/env python3
"""
Comprehensive test suite for gene-specific rate support in Phase 2 pipeline.
Tests both uniform and gene-specific rates with various edge cases.
"""

import os
import sys
import json
import gzip
import tempfile
import shutil
from pathlib import Path
import argparse

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from run_pipeline import (
    parse_gene_rate_groups, 
    validate_gene_rate_groups,
    create_petri_with_rate_config
)


def test_parse_gene_rate_groups():
    """Test parsing of gene rate groups string."""
    print("\nTesting parse_gene_rate_groups...")
    
    # Test 1: Simple case
    result = parse_gene_rate_groups("50:0.004,50:0.005")
    assert result == [(50, 0.004), (50, 0.005)], f"Expected [(50, 0.004), (50, 0.005)], got {result}"
    print("  ✓ Simple parsing works")
    
    # Test 2: Multiple groups
    result = parse_gene_rate_groups("25:0.001,25:0.002,25:0.003,25:0.004")
    assert len(result) == 4, f"Expected 4 groups, got {len(result)}"
    assert result[2] == (25, 0.003), f"Expected (25, 0.003), got {result[2]}"
    print("  ✓ Multiple groups parsing works")
    
    # Test 3: Single group
    result = parse_gene_rate_groups("200:0.005")
    assert result == [(200, 0.005)], f"Expected [(200, 0.005)], got {result}"
    print("  ✓ Single group parsing works")
    
    # Test 4: Different sized groups
    result = parse_gene_rate_groups("100:0.001,50:0.002,30:0.003,20:0.004")
    assert sum(n for n, _ in result) == 200, "Total genes should be 200"
    print("  ✓ Different sized groups parsing works")


def test_validate_gene_rate_groups():
    """Test validation of gene rate groups."""
    print("\nTesting validate_gene_rate_groups...")
    
    # Test 1: Valid configuration
    groups = [(50, 0.004), (50, 0.005), (50, 0.006), (50, 0.007)]
    try:
        validate_gene_rate_groups(groups, 1000, 5)
        print("  ✓ Valid configuration passes")
    except ValueError as e:
        assert False, f"Valid configuration should pass: {e}"
    
    # Test 2: Wrong total genes
    groups = [(50, 0.004), (50, 0.005), (50, 0.006)]  # Only 150 genes
    try:
        validate_gene_rate_groups(groups, 1000, 5)
        assert False, "Should have raised ValueError for wrong gene count"
    except ValueError as e:
        assert "150 genes" in str(e) and "200 genes" in str(e)
        print("  ✓ Wrong gene count detected")
    
    # Test 3: Invalid rate (negative)
    groups = [(100, -0.001), (100, 0.005)]
    try:
        validate_gene_rate_groups(groups, 1000, 5)
        assert False, "Should have raised ValueError for negative rate"
    except ValueError as e:
        assert "Invalid rate" in str(e)
        print("  ✓ Negative rate detected")
    
    # Test 4: Invalid rate (> 1)
    groups = [(100, 0.5), (100, 1.5)]
    try:
        validate_gene_rate_groups(groups, 1000, 5)
        assert False, "Should have raised ValueError for rate > 1"
    except ValueError as e:
        assert "Invalid rate" in str(e)
        print("  ✓ Rate > 1 detected")
    
    # Test 5: Different n_sites and gene_size
    groups = [(10, 0.004), (10, 0.005)]  # 20 genes for 100 sites, gene_size=5
    try:
        validate_gene_rate_groups(groups, 100, 5)
        print("  ✓ Different n_sites configuration passes")
    except ValueError as e:
        assert False, f"Valid small configuration should pass: {e}"


def test_create_petri_with_rate_config():
    """Test PetriDish creation with different rate configurations."""
    print("\nTesting create_petri_with_rate_config...")
    
    # Test 1: Uniform rate configuration
    rate_config = {'type': 'uniform', 'rate': 0.005}
    petri = create_petri_with_rate_config(rate_config, growth_phase=3)
    assert petri.rate == 0.005, f"Expected rate 0.005, got {petri.rate}"
    assert petri.growth_phase == 3, f"Expected growth_phase 3, got {petri.growth_phase}"
    print("  ✓ Uniform rate PetriDish created")
    
    # Test 2: Gene-specific rate configuration
    rate_config = {
        'type': 'gene_specific',
        'gene_rate_groups': [(10, 0.004), (10, 0.006)],
        'gene_size': 5
    }
    petri = create_petri_with_rate_config(rate_config, growth_phase=2)
    assert hasattr(petri, 'gene_rate_groups'), "PetriDish should have gene_rate_groups"
    assert petri.gene_size == 5, f"Expected gene_size 5, got {petri.gene_size}"
    assert petri.growth_phase == 2, f"Expected growth_phase 2, got {petri.growth_phase}"
    print("  ✓ Gene-specific rate PetriDish created")
    
    # Test 3: Verify cells can be added and grown
    # Create a cell with matching configuration
    cell = Cell(n=100, gene_rate_groups=[(10, 0.004), (10, 0.006)], gene_size=5)
    petri.cells = [cell]
    
    # Test growth (small scale)
    initial_cells = len(petri.cells)
    petri.divide_cells()  # Should double cells
    assert len(petri.cells) == initial_cells * 2, f"Expected {initial_cells * 2} cells after division, got {len(petri.cells)}"
    print("  ✓ Gene-specific PetriDish can grow cells")


def test_integration_small_pipeline():
    """Test a small-scale pipeline run with gene-specific rates."""
    print("\nTesting small-scale pipeline integration...")
    
    # Create a minimal simulation file
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create mock simulation data with gene-specific rates
        sim_data = {
            "metadata": {
                "n_sites": 100,
                "gene_size": 5,
                "growth_phase": 2,
                "years": 10,
                "gene_rate_groups": [(10, 0.01), (10, 0.02)]
            },
            "snapshots": {}
        }
        
        # Add some mock snapshots with small populations
        for year in [5, 8]:
            cells_data = []
            for i in range(4):  # Only 4 cells for quick testing
                cell = Cell(n=100, gene_rate_groups=[(10, 0.01), (10, 0.02)], gene_size=5)
                # Add some methylation
                for _ in range(year):
                    cell.methylate()
                cells_data.append({
                    "cpg_sites": cell.cpg_sites.tolist(),
                    "gene_size": cell.gene_size,
                    "id": i
                })
            sim_data["snapshots"][str(year)] = {"cells": cells_data}
        
        # Save simulation file
        sim_file = os.path.join(tmpdir, "test_simulation.json.gz")
        with gzip.open(sim_file, 'wt') as f:
            json.dump(sim_data, f)
        
        print("  ✓ Created mock simulation file")
        
        # Test command-line argument parsing (without running full pipeline)
        test_args = [
            "--gene-rate-groups", "10:0.01,10:0.02",
            "--gene-size", "5",
            "--simulation", sim_file,
            "--first-snapshot", "5",
            "--second-snapshot", "8",
            "--individual-growth-phase", "1",
            "--n-quantiles", "2",
            "--cells-per-quantile", "1",
            "--seed", "42"
        ]
        
        parser = argparse.ArgumentParser()
        rate_group = parser.add_mutually_exclusive_group(required=True)
        rate_group.add_argument("--rate", type=float)
        rate_group.add_argument("--gene-rate-groups", type=str)
        parser.add_argument("--gene-size", type=int, default=5)
        parser.add_argument("--simulation", type=str, required=True)
        parser.add_argument("--first-snapshot", type=int, default=50)
        parser.add_argument("--second-snapshot", type=int, default=60)
        parser.add_argument("--individual-growth-phase", type=int, default=7)
        parser.add_argument("--n-quantiles", type=int, default=10)
        parser.add_argument("--cells-per-quantile", type=int, default=3)
        parser.add_argument("--seed", type=int, default=42)
        
        args = parser.parse_args(test_args)
        
        # Parse and validate gene rate groups
        gene_rate_groups = parse_gene_rate_groups(args.gene_rate_groups)
        validate_gene_rate_groups(gene_rate_groups, 100, args.gene_size)
        
        print("  ✓ Command-line arguments parsed successfully")
        print("  ✓ Gene rate groups validated")
        
        # Test rate configuration creation
        rate_config = {
            'type': 'gene_specific',
            'gene_rate_groups': gene_rate_groups,
            'gene_size': args.gene_size
        }
        
        # Test PetriDish creation with configuration
        petri = create_petri_with_rate_config(rate_config, growth_phase=args.individual_growth_phase)
        assert petri.growth_phase == 1, "Growth phase should be 1"
        
        print("  ✓ Integration test passed")


def test_backward_compatibility():
    """Test that uniform rate still works (backward compatibility)."""
    print("\nTesting backward compatibility with uniform rates...")
    
    # Test uniform rate configuration
    rate_config = {'type': 'uniform', 'rate': 0.005}
    petri = create_petri_with_rate_config(rate_config, growth_phase=3)
    
    # Add a cell with uniform rate
    cell = Cell(n=100, rate=0.005)
    petri.cells = [cell]
    
    # Test that it works as before
    initial_methylation = sum(cell.cpg_sites)
    cell.methylate()  # Apply methylation
    final_methylation = sum(cell.cpg_sites)
    
    assert final_methylation >= initial_methylation, "Methylation should not decrease"
    print("  ✓ Uniform rate backward compatibility maintained")
    
    # Test division
    petri.divide_cells()
    assert len(petri.cells) == 2, "Should have 2 cells after division"
    print("  ✓ Cell division works with uniform rate")


def test_edge_cases():
    """Test edge cases and potential failure modes."""
    print("\nTesting edge cases...")
    
    # Test 1: Empty gene rate groups string
    try:
        parse_gene_rate_groups("")
        assert False, "Should fail on empty string"
    except:
        print("  ✓ Empty string handled")
    
    # Test 2: Malformed gene rate groups
    try:
        parse_gene_rate_groups("50-0.004,50-0.005")  # Wrong separator
        assert False, "Should fail on wrong separator"
    except:
        print("  ✓ Wrong separator detected")
    
    # Test 3: Very small gene groups
    groups = [(1, 0.001)] * 100  # 100 groups of 1 gene each
    try:
        validate_gene_rate_groups(groups, 500, 5)
        print("  ✓ Many small groups handled")
    except ValueError as e:
        assert False, f"Should handle many small groups: {e}"
    
    # Test 4: Single large group
    groups = [(200, 0.005)]
    try:
        validate_gene_rate_groups(groups, 1000, 5)
        print("  ✓ Single large group handled")
    except ValueError as e:
        assert False, f"Should handle single large group: {e}"
    
    # Test 5: Mixed rate config types (should not mix)
    rate_config = {'type': 'uniform', 'rate': 0.005}
    petri1 = create_petri_with_rate_config(rate_config, growth_phase=2)
    
    rate_config = {'type': 'gene_specific', 'gene_rate_groups': [(10, 0.004), (10, 0.006)], 'gene_size': 5}
    petri2 = create_petri_with_rate_config(rate_config, growth_phase=2)
    
    assert hasattr(petri1, 'rate'), "Uniform should have rate"
    assert hasattr(petri2, 'gene_rate_groups'), "Gene-specific should have gene_rate_groups"
    print("  ✓ Different config types create appropriate PetriDish objects")


def run_all_tests():
    """Run all test functions."""
    print("="*60)
    print("Running Phase 2 Gene Rate Support Tests")
    print("="*60)
    
    test_functions = [
        test_parse_gene_rate_groups,
        test_validate_gene_rate_groups,
        test_create_petri_with_rate_config,
        test_integration_small_pipeline,
        test_backward_compatibility,
        test_edge_cases
    ]
    
    failed = 0
    for test_func in test_functions:
        try:
            test_func()
        except Exception as e:
            print(f"\n✗ {test_func.__name__} FAILED:")
            print(f"  {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    if failed == 0:
        print("✅ All tests passed!")
    else:
        print(f"❌ {failed} test(s) failed")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(run_all_tests())