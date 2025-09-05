#!/usr/bin/env python3
"""
Test rate consistency validation in Phase 2 pipeline.

Tests ensure that:
1. Simulation gene_rate_groups match arguments
2. All cells have consistent gene_rate_groups
3. Mixing validates compatibility
4. Error messages are clear
"""

import sys
import os
import json
import gzip
import tempfile
import shutil

# Add paths for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish, rate_to_gene_rate_groups
from core.pipeline_utils import (
    dict_to_cell, load_snapshot_as_cells, save_snapshot_cells,
    load_snapshot_cells, create_pure_snapshot_petri
)


def create_mock_simulation(gene_rate_groups, n=100, gene_size=5, years=10):
    """Create a mock simulation file for testing."""
    
    # Create cells for each year
    history = {}
    for year in range(years + 1):
        cells = []
        for i in range(10):  # 10 cells per year
            cell = Cell(n=n, gene_rate_groups=gene_rate_groups, gene_size=gene_size)
            cell.age = year
            # Add some methylation
            for j in range(year):
                if j < n:
                    cell.cpg_sites[j] = 1
            # Set JSD directly for test
            cell.cell_jsd = year * 0.01
            cells.append(cell.to_dict())
        
        history[str(year)] = {'cells': cells}
    
    # Create simulation data
    sim_data = {
        'parameters': {
            'gene_rate_groups': gene_rate_groups,
            'n': n,
            'gene_size': gene_size,
            'years': years,
            'growth_phase': 3,
            'seed': 42
        },
        'history': history
    }
    
    return sim_data


def test_dict_to_cell_requires_gene_rate_groups():
    """Test that dict_to_cell requires gene_rate_groups in cell data."""
    print("\n" + "="*60)
    print("TEST: dict_to_cell requires gene_rate_groups")
    print("="*60)
    
    # Cell dict without gene_rate_groups
    cell_dict = {
        'methylated': [0, 1, 0, 1, 0],
        'age': 10,
        'cell_jsd': 0.123
    }
    
    try:
        cell = dict_to_cell(cell_dict)
        print("✗ Should have failed for missing gene_rate_groups")
        return False
    except ValueError as e:
        if "missing gene_rate_groups" in str(e):
            print("✓ Correctly rejected cell without gene_rate_groups")
            print(f"  Error: {str(e)[:100]}...")
            return True
        else:
            print(f"✗ Wrong error: {e}")
            return False


def test_snapshot_loading_validates_groups():
    """Test that loading snapshots validates gene_rate_groups."""
    print("\n" + "="*60)
    print("TEST: Snapshot loading validates gene_rate_groups")
    print("="*60)
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create simulation with specific gene_rate_groups
        sim_groups = [(10, 0.004), (10, 0.006)]
        sim_data = create_mock_simulation(sim_groups, n=100, gene_size=5)
        
        # Save simulation
        sim_path = os.path.join(tmpdir, "simulation.json.gz")
        with gzip.open(sim_path, 'wt') as f:
            json.dump(sim_data, f)
        
        # Test 1: Loading with correct expected groups should work
        try:
            cells = load_snapshot_as_cells(sim_path, 5, expected_gene_rate_groups=sim_groups)
            print(f"✓ Loaded snapshot with matching gene_rate_groups: {sim_groups}")
        except Exception as e:
            print(f"✗ Failed to load with matching groups: {e}")
            return False
        
        # Test 2: Loading with wrong expected groups should fail
        wrong_groups = [(20, 0.005)]
        try:
            cells = load_snapshot_as_cells(sim_path, 5, expected_gene_rate_groups=wrong_groups)
            print("✗ Should have failed for mismatched gene_rate_groups")
            return False
        except ValueError as e:
            if "mismatch" in str(e).lower():
                print("✓ Correctly rejected mismatched gene_rate_groups")
                return True
            else:
                print(f"✗ Wrong error: {e}")
                return False


def test_mixing_validation():
    """Test that mixing cells validates compatibility."""
    print("\n" + "="*60)
    print("TEST: Mixing validation")
    print("="*60)
    
    # Create cells with same gene_rate_groups
    groups1 = [(10, 0.005)]
    cells1 = [Cell(n=50, gene_rate_groups=groups1, gene_size=5) for _ in range(5)]
    
    # Create cells with different gene_rate_groups
    groups2 = [(10, 0.006)]
    cells2 = [Cell(n=50, gene_rate_groups=groups2, gene_size=5) for _ in range(5)]
    
    # Test 1: Mixing compatible cells should work
    try:
        PetriDish.validate_cells_compatible(cells1 + cells1)
        print("✓ Compatible cells pass validation")
    except Exception as e:
        print(f"✗ Compatible cells failed: {e}")
        return False
    
    # Test 2: Mixing incompatible cells should fail
    try:
        PetriDish.validate_cells_compatible(cells1 + cells2)
        print("✗ Incompatible cells should have failed validation")
        return False
    except ValueError as e:
        if "mismatch" in str(e).lower():
            print("✓ Incompatible cells correctly rejected")
            print(f"  Error shows clear mismatch: {groups1} vs {groups2}")
            return True
        else:
            print(f"✗ Wrong error: {e}")
            return False


def test_control2_creation():
    """Test that Control2 creation maintains gene_rate_groups."""
    print("\n" + "="*60)
    print("TEST: Control2 creation maintains gene_rate_groups")
    print("="*60)
    
    # Create snapshot cells with specific gene_rate_groups
    groups = [(20, 0.004)]
    snapshot_cells = [Cell(n=100, gene_rate_groups=groups, gene_size=5) for _ in range(100)]
    
    # Set age for all cells
    for i, cell in enumerate(snapshot_cells):
        cell.age = 50
        # Add some methylation
        for j in range(25):
            cell.cpg_sites[j] = 1
        # Set JSD directly for test
        cell.cell_jsd = 0.25
    
    # Create Control2 PetriDish
    try:
        control2 = create_pure_snapshot_petri(snapshot_cells, n_cells=50, seed=42)
        
        # Verify all cells have correct gene_rate_groups
        for cell in control2.cells:
            if cell.gene_rate_groups != groups:
                print(f"✗ Cell has wrong gene_rate_groups: {cell.gene_rate_groups}")
                return False
        
        # Verify PetriDish has correct gene_rate_groups
        if control2.gene_rate_groups != groups:
            print(f"✗ PetriDish has wrong gene_rate_groups: {control2.gene_rate_groups}")
            return False
        
        print(f"✓ Control2 maintains gene_rate_groups: {groups}")
        print(f"  Created with {len(control2.cells)} cells")
        return True
        
    except Exception as e:
        print(f"✗ Control2 creation failed: {e}")
        return False


def test_rate_to_gene_rate_groups_conversion():
    """Test conversion from uniform rate to gene_rate_groups."""
    print("\n" + "="*60)
    print("TEST: Rate to gene_rate_groups conversion")
    print("="*60)
    
    # Test conversion
    rate = 0.005
    n = 100
    gene_size = 5
    
    groups = rate_to_gene_rate_groups(rate, n, gene_size)
    
    # Verify format
    expected_groups = [(20, 0.005)]  # 100 sites / 5 sites per gene = 20 genes
    
    if groups == expected_groups:
        print(f"✓ Correctly converted rate {rate} to groups {groups}")
        print(f"  {n} sites / {gene_size} sites per gene = {n//gene_size} genes")
        return True
    else:
        print(f"✗ Wrong conversion: expected {expected_groups}, got {groups}")
        return False


def test_error_messages():
    """Test that error messages are clear and actionable."""
    print("\n" + "="*60)
    print("TEST: Error message quality")
    print("="*60)
    
    # Create cells with different gene_rate_groups
    cell1 = Cell(n=50, gene_rate_groups=[(10, 0.004)], gene_size=5)
    cell2 = Cell(n=50, gene_rate_groups=[(10, 0.006)], gene_size=5)
    
    try:
        PetriDish.validate_cells_compatible([cell1, cell2])
    except ValueError as e:
        error_msg = str(e)
        print("Error message received:")
        print(error_msg)
        print()
        
        # Check error message quality
        checks = [
            ("mismatch" in error_msg.lower(), "Contains 'mismatch'"),
            ("Cell 0:" in error_msg, "Identifies first cell"),
            ("Cell 1:" in error_msg, "Identifies second cell"),
            ("gene_rate_groups=" in error_msg, "Shows parameter name"),
            ("[(10, 0.004)]" in error_msg, "Shows first cell's groups"),
            ("[(10, 0.006)]" in error_msg, "Shows second cell's groups")
        ]
        
        for check, description in checks:
            if check:
                print(f"✓ {description}")
            else:
                print(f"✗ {description}")
        
        return all(check for check, _ in checks)
    
    print("✗ Should have raised an error")
    return False


def main():
    """Run all tests."""
    print("="*60)
    print("PHASE 2 RATE CONSISTENCY VALIDATION TESTS")
    print("="*60)
    
    tests = [
        ("dict_to_cell requires gene_rate_groups", test_dict_to_cell_requires_gene_rate_groups),
        ("Snapshot loading validates groups", test_snapshot_loading_validates_groups),
        ("Mixing validation", test_mixing_validation),
        ("Control2 creation", test_control2_creation),
        ("Rate to gene_rate_groups conversion", test_rate_to_gene_rate_groups_conversion),
        ("Error messages", test_error_messages),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            passed = test_func()
            results.append((test_name, passed))
        except Exception as e:
            print(f"\n✗ Test '{test_name}' crashed: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for test_name, passed in results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{status}: {test_name}")
    
    passed_count = sum(1 for _, passed in results if passed)
    total_count = len(results)
    
    print(f"\nTotal: {passed_count}/{total_count} tests passed")
    
    if passed_count == total_count:
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n❌ {total_count - passed_count} test(s) failed")
        return 1


if __name__ == "__main__":
    exit(main())