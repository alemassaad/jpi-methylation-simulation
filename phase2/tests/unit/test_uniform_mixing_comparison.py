#!/usr/bin/env python3
"""
Comparison tests for uniform mixing vs regular Control2.
Verifies statistical properties of the implementation.
"""

import sys
import os
import random

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import (
    create_uniform_mixing_pool,
    create_control2_with_uniform_base,
    create_pure_snapshot_petri,
    mix_petri_uniform
)
from test_helpers import (
    create_mock_cells,
    create_mock_petri_dish,
    get_cell_fingerprint,
    cells_are_identical
)


def test_uniform_vs_regular_statistics():
    """Compare statistical properties of uniform vs regular Control2."""
    print("\nTest 4.1: Uniform vs regular Control2 statistics")
    
    # Create snapshot
    snapshot_cells = create_mock_cells(n_cells=500, n_sites=50, seed=42)
    
    # Parameters
    median_size = 20
    mix_ratio = 0.8
    target_size = int(median_size / (1 - mix_ratio))  # 100
    base_size = int(target_size * mix_ratio)  # 80
    n_individuals = 10
    
    # Create uniform pool
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size, mix_ratio, seed=100
    )
    
    # Create uniform-based Control2 individuals
    print(f"  Creating {n_individuals} uniform-based Control2 individuals...")
    uniform_individuals = []
    for i in range(n_individuals):
        control2 = create_control2_with_uniform_base(
            snapshot_cells, uniform_pool, uniform_indices,
            target_size, rate=0.005, seed=200 + i
        )
        uniform_individuals.append(control2)
    
    # Create regular Control2 individuals
    print(f"  Creating {n_individuals} regular Control2 individuals...")
    regular_individuals = []
    for i in range(n_individuals):
        control2 = create_pure_snapshot_petri(
            snapshot_cells, n_cells=target_size,
            rate=0.005, seed=300 + i
        )
        regular_individuals.append(control2)
    
    # Analyze uniform-based individuals
    print("  Analyzing uniform-based variance...")
    uniform_base_variance = []
    for pos in range(base_size):
        # Get fingerprints at this position across all individuals
        fingerprints_at_pos = []
        for ind in uniform_individuals:
            fingerprints_at_pos.append(get_cell_fingerprint(ind.cells[pos]))
        
        # Count unique fingerprints
        unique_count = len(set(fingerprints_at_pos))
        uniform_base_variance.append(unique_count)
    
    # All positions should have exactly 1 unique fingerprint (no variance)
    assert all(v == 1 for v in uniform_base_variance), \
        f"Uniform base should have no variance, got {uniform_base_variance[:5]}..."
    
    print(f"    ✓ All {base_size} base positions have identical cells across individuals")
    
    # Analyze regular individuals
    print("  Analyzing regular variance...")
    regular_variance = []
    for pos in range(base_size):
        fingerprints_at_pos = []
        for ind in regular_individuals:
            fingerprints_at_pos.append(get_cell_fingerprint(ind.cells[pos]))
        
        unique_count = len(set(fingerprints_at_pos))
        regular_variance.append(unique_count)
    
    # Most positions should have variance
    positions_with_variance = sum(v > 1 for v in regular_variance)
    variance_ratio = positions_with_variance / base_size
    
    print(f"    Regular: {positions_with_variance}/{base_size} positions vary ({variance_ratio:.1%})")
    assert variance_ratio > 0.8, "Regular Control2 should have variance in most positions"
    
    print("  ✓ Statistical properties verified")
    return True


def test_three_group_alignment():
    """Test that all three groups share the same base with uniform mixing."""
    print("\nTest 4.2: Three-group alignment with uniform mixing")
    
    # Create snapshot
    snapshot_cells = create_mock_cells(n_cells=200, n_sites=50, seed=42)
    
    # Parameters
    median_size = 20
    mix_ratio = 0.8
    target_size = int(median_size / (1 - mix_ratio))  # 100
    base_size = int(target_size * mix_ratio)  # 80
    
    # Create uniform pool
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size, mix_ratio, seed=100
    )
    
    # Simulate mutant (with grown cells)
    print("  Creating mutant individual...")
    mutant = PetriDish(rate=0.005, n=50, gene_size=5, seed=None)
    # Start with uniform pool
    mutant.cells = [c for c in uniform_pool]
    # Add "grown" cells (just different mock cells for testing)
    grown_cells = create_mock_cells(n_cells=20, n_sites=50, seed=500)
    mutant.cells.extend(grown_cells)
    
    # Simulate control1 (with different grown cells)
    print("  Creating control1 individual...")
    control1 = PetriDish(rate=0.005, n=50, gene_size=5, seed=None)
    control1.cells = [c for c in uniform_pool]
    grown_cells2 = create_mock_cells(n_cells=20, n_sites=50, seed=600)
    control1.cells.extend(grown_cells2)
    
    # Create control2 with uniform base
    print("  Creating control2 individual...")
    control2 = create_control2_with_uniform_base(
        snapshot_cells, uniform_pool, uniform_indices,
        target_size, rate=0.005, seed=700
    )
    
    # Verify base alignment
    print(f"  Checking base alignment (first {base_size} cells)...")
    mismatches = 0
    for pos in range(base_size):
        mut_fp = get_cell_fingerprint(mutant.cells[pos])
        con1_fp = get_cell_fingerprint(control1.cells[pos])
        con2_fp = get_cell_fingerprint(control2.cells[pos])
        
        if not (mut_fp == con1_fp == con2_fp):
            mismatches += 1
    
    assert mismatches == 0, f"Found {mismatches} mismatches in base cells"
    print(f"    ✓ All {base_size} base cells identical across groups")
    
    # Verify additional portions differ
    print(f"  Checking additional cells (last {target_size - base_size})...")
    
    # Mutant vs Control1 grown cells should differ
    mut_grown = mutant.cells[base_size:]
    con1_grown = control1.cells[base_size:]
    
    identical_count = sum(
        cells_are_identical(mut_grown[i], con1_grown[i])
        for i in range(len(mut_grown))
    )
    
    # Should be mostly different (allowing for random matches)
    assert identical_count < len(mut_grown) * 0.2, \
        "Mutant and Control1 grown cells should differ"
    
    print(f"    ✓ Mutant and Control1 have different grown cells ({identical_count}/{len(mut_grown)} match)")
    
    # Control2 additional cells should also differ from mutant/control1
    con2_additional = control2.cells[base_size:]
    
    # Check against mutant grown
    identical_to_mut = sum(
        cells_are_identical(con2_additional[i], mut_grown[i])
        for i in range(len(con2_additional))
    )
    
    assert identical_to_mut < len(con2_additional) * 0.2, \
        "Control2 additional cells should differ from mutant grown"
    
    print(f"    ✓ Control2 additional cells differ from mutant ({identical_to_mut}/{len(con2_additional)} match)")
    
    print("  ✓ Three-group alignment verified")
    return True


def test_mixing_consistency():
    """Test that uniform mixing produces consistent results."""
    print("\nTest 4.3: Mixing consistency")
    
    # Create snapshot
    snapshot_cells = create_mock_cells(n_cells=200, n_sites=50, seed=42)
    
    # Create uniform pool
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size=20, mix_ratio=0.8, seed=100
    )
    
    # Create multiple Control2 individuals
    individuals = []
    for i in range(5):
        control2 = create_control2_with_uniform_base(
            snapshot_cells, uniform_pool, uniform_indices,
            target_size=100, rate=0.005, seed=200 + i
        )
        individuals.append(control2)
    
    # Verify all share the same base (first 80 cells)
    print("  Verifying shared base across 5 individuals...")
    base_size = len(uniform_pool)
    
    for i in range(1, len(individuals)):
        for pos in range(base_size):
            assert cells_are_identical(
                individuals[0].cells[pos],
                individuals[i].cells[pos]
            ), f"Individual {i} differs at base position {pos}"
    
    print(f"    ✓ All individuals share identical {base_size}-cell base")
    
    # Verify additional cells differ
    print("  Verifying unique additional cells...")
    for i in range(len(individuals)):
        for j in range(i + 1, len(individuals)):
            # Compare additional cells
            add_i = individuals[i].cells[base_size:]
            add_j = individuals[j].cells[base_size:]
            
            # Count matches
            matches = sum(
                cells_are_identical(add_i[k], add_j[k])
                for k in range(len(add_i))
            )
            
            # Should be mostly different
            match_ratio = matches / len(add_i)
            assert match_ratio < 0.3, \
                f"Individuals {i} and {j} have too similar additional cells ({match_ratio:.1%})"
    
    print("    ✓ Each individual has unique additional cells")
    
    print("  ✓ Mixing consistency verified")
    return True


def run_all_tests():
    """Run all comparison tests."""
    print("="*60)
    print("COMPARISON TESTS FOR UNIFORM MIXING")
    print("="*60)
    
    tests = [
        test_uniform_vs_regular_statistics,
        test_three_group_alignment,
        test_mixing_consistency
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except AssertionError as e:
            print(f"  ✗ Test failed: {e}")
            failed += 1
        except Exception as e:
            print(f"  ✗ Test error: {e}")
            failed += 1
    
    print("\n" + "="*60)
    print(f"COMPARISON TEST RESULTS: {passed} passed, {failed} failed")
    print("="*60)
    
    return failed == 0


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)