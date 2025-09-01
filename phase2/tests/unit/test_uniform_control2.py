#!/usr/bin/env python3
"""
Test that Control2 shares the same base with mutant/control1 when using --uniform-mixing.
"""

import os
import sys
import json
import gzip
import random
import numpy as np
from typing import List, Dict, Tuple

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import (
    create_uniform_mixing_pool, 
    create_control2_with_uniform_base,
    create_pure_snapshot_petri,
    mix_petri_uniform
)


def create_test_cells(n_cells: int, n_sites: int = 100, seed: int = 42) -> List[Cell]:
    """Create test cells with varying methylation patterns."""
    random.seed(seed)
    np.random.seed(seed)
    
    cells = []
    for i in range(n_cells):
        cell = Cell(n=n_sites, rate=0.005, gene_size=5, seed=seed + i)
        # Randomly methylate some sites
        for j in range(n_sites):
            if random.random() < 0.1 * (i / n_cells):  # Gradient of methylation
                cell.cpg_sites[j] = 1
                cell.methylated[j] = True
        cell.age = 60  # Simulate year 60 snapshot
        cells.append(cell)
    
    return cells


def test_uniform_base_sharing():
    """Verify all groups share the same 80% base cells when using uniform mixing."""
    print("Testing uniform base sharing...")
    
    # Create test snapshot
    snapshot_cells = create_test_cells(1000, seed=42)
    
    # Create uniform pool (80% of median size)
    median_size = 100
    mix_ratio = 0.8
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size, mix_ratio, seed=100
    )
    
    # Verify pool size
    expected_pool_size = int(median_size / (1 - mix_ratio) * mix_ratio)
    assert len(uniform_pool) == expected_pool_size, f"Pool size mismatch: {len(uniform_pool)} != {expected_pool_size}"
    assert len(uniform_indices) == expected_pool_size, f"Indices count mismatch"
    
    # Create Control2 individuals with uniform base
    control2_individuals = []
    for i in range(3):
        target_size = int(median_size / (1 - mix_ratio))
        petri = create_control2_with_uniform_base(
            snapshot_cells, uniform_pool, uniform_indices,
            target_size, rate=0.005, seed=200 + i
        )
        control2_individuals.append(petri)
    
    # Verify all Control2 individuals have the same base
    print(f"  Created {len(control2_individuals)} Control2 individuals")
    
    # Check that first 80% of cells are identical across individuals
    base_size = len(uniform_pool)
    for i in range(1, len(control2_individuals)):
        for j in range(base_size):
            cell1 = control2_individuals[0].cells[j]
            cell2 = control2_individuals[i].cells[j]
            
            # Compare methylation patterns
            assert np.array_equal(cell1.methylated, cell2.methylated), \
                f"Base cells differ between individual 0 and {i} at position {j}"
    
    print("  ✓ All Control2 individuals share the same base cells")
    
    # Verify the additional 20% differs between individuals
    additional_start = base_size
    for i in range(len(control2_individuals)):
        for j in range(i + 1, len(control2_individuals)):
            # Check that at least some additional cells differ
            all_identical = True
            for k in range(additional_start, len(control2_individuals[i].cells)):
                cell1 = control2_individuals[i].cells[k]
                cell2 = control2_individuals[j].cells[k]
                if not np.array_equal(cell1.methylated, cell2.methylated):
                    all_identical = False
                    break
            
            assert not all_identical, f"Additional cells identical between individuals {i} and {j}"
    
    print("  ✓ Additional cells differ between individuals")
    

def test_no_duplicate_cells():
    """Ensure no cell appears twice in Control2."""
    print("\nTesting for duplicate cells in Control2...")
    
    # Create test snapshot
    snapshot_cells = create_test_cells(500, seed=42)
    
    # Create uniform pool
    median_size = 50
    mix_ratio = 0.8
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size, mix_ratio, seed=100
    )
    
    # Create Control2 individual
    target_size = int(median_size / (1 - mix_ratio))
    petri = create_control2_with_uniform_base(
        snapshot_cells, uniform_pool, uniform_indices,
        target_size, rate=0.005, seed=200
    )
    
    # Check for duplicates by comparing methylation patterns
    patterns_seen = set()
    for i, cell in enumerate(petri.cells):
        pattern = tuple(cell.methylated)
        if pattern in patterns_seen:
            # Double-check it's really a duplicate (not just same pattern by chance)
            for j in range(i):
                if tuple(petri.cells[j].methylated) == pattern:
                    if petri.cells[j].age == cell.age:  # Same age too
                        assert False, f"Duplicate cell found at positions {j} and {i}"
        patterns_seen.add(pattern)
    
    print(f"  ✓ No duplicate cells found in Control2 ({len(petri.cells)} cells checked)")


def test_edge_case_100_percent_mix():
    """Test edge case where mix ratio is 100%."""
    print("\nTesting 100% mix ratio edge case...")
    
    # Create test snapshot
    snapshot_cells = create_test_cells(500, seed=42)
    
    # Create uniform pool with 100% mix ratio
    median_size = 50
    mix_ratio = 1.0
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size, mix_ratio, seed=100
    )
    
    # With 100% mix, pool should equal median size
    assert len(uniform_pool) == median_size, f"Pool size should equal median with 100% mix"
    
    # Create Control2 individual
    petri = create_control2_with_uniform_base(
        snapshot_cells, uniform_pool, uniform_indices,
        median_size, rate=0.005, seed=200
    )
    
    # Should have exactly median_size cells, all from uniform pool
    assert len(petri.cells) == median_size, f"Control2 size mismatch with 100% mix"
    
    # All cells should match uniform pool
    for i in range(median_size):
        assert np.array_equal(petri.cells[i].methylated, uniform_pool[i].methylated), \
            f"Cell {i} doesn't match uniform pool with 100% mix"
    
    print(f"  ✓ 100% mix ratio handled correctly")


def test_comparison_with_regular_control2():
    """Compare uniform-base Control2 with regular Control2."""
    print("\nComparing uniform-base Control2 with regular Control2...")
    
    # Create test snapshot
    snapshot_cells = create_test_cells(1000, seed=42)
    
    # Create uniform pool
    median_size = 100
    mix_ratio = 0.8
    uniform_pool, uniform_indices = create_uniform_mixing_pool(
        snapshot_cells, median_size, mix_ratio, seed=100
    )
    
    target_size = int(median_size / (1 - mix_ratio))
    
    # Create uniform-base Control2
    uniform_control2 = create_control2_with_uniform_base(
        snapshot_cells, uniform_pool, uniform_indices,
        target_size, rate=0.005, seed=200
    )
    
    # Create regular Control2
    regular_control2 = create_pure_snapshot_petri(
        snapshot_cells, n_cells=target_size,
        rate=0.005, seed=200
    )
    
    # Both should have same number of cells
    assert len(uniform_control2.cells) == len(regular_control2.cells), \
        "Control2 variants have different sizes"
    
    # But they should have different cells (except by chance)
    n_identical = 0
    for i in range(len(uniform_control2.cells)):
        if np.array_equal(uniform_control2.cells[i].methylated, 
                         regular_control2.cells[i].methylated):
            n_identical += 1
    
    # Some overlap is possible by chance, but not all
    assert n_identical < len(uniform_control2.cells), \
        "Uniform and regular Control2 are unexpectedly identical"
    
    print(f"  ✓ Uniform-base and regular Control2 differ as expected")
    print(f"    {n_identical}/{target_size} cells identical by chance")


def main():
    """Run all tests."""
    print("="*60)
    print("Testing Uniform Control2 Implementation")
    print("="*60)
    
    test_uniform_base_sharing()
    test_no_duplicate_cells()
    test_edge_case_100_percent_mix()
    test_comparison_with_regular_control2()
    
    print("\n" + "="*60)
    print("All tests passed! ✅")
    print("="*60)


if __name__ == "__main__":
    main()