#!/usr/bin/env python3
"""
Test rate consistency validation in PetriDish.

Tests ensure that:
1. All cells in a PetriDish have identical gene_rate_groups
2. Validation catches mismatches
3. Cell division maintains consistency
4. Error messages are helpful
"""

import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from cell import Cell, PetriDish


def test_identical_uniform_rates():
    """Test cells with same uniform rate are compatible."""
    print("\n" + "="*60)
    print("TEST: Identical uniform rates")
    print("="*60)
    
    # Create cells with same uniform rate
    cell1 = Cell(n=100, rate=0.005, gene_size=5)
    cell2 = Cell(n=100, rate=0.005, gene_size=5)
    cell3 = Cell(n=100, rate=0.005, gene_size=5)
    
    # Should work with static method
    try:
        PetriDish.validate_cells_compatible([cell1, cell2, cell3])
        print("✓ Static validation passed for identical uniform rates")
    except ValueError as e:
        print(f"✗ Static validation failed: {e}")
        return False
    
    # Should work when creating PetriDish
    try:
        petri = PetriDish(cells=[cell1, cell2, cell3], n=100, gene_size=5)
        print("✓ PetriDish creation passed with identical uniform rates")
        print(f"  All cells have gene_rate_groups: {cell1.gene_rate_groups}")
    except ValueError as e:
        print(f"✗ PetriDish creation failed: {e}")
        return False
    
    return True


def test_different_uniform_rates():
    """Test cells with different uniform rates are incompatible."""
    print("\n" + "="*60)
    print("TEST: Different uniform rates")
    print("="*60)
    
    # Create cells with different uniform rates
    cell1 = Cell(n=100, rate=0.004, gene_size=5)
    cell2 = Cell(n=100, rate=0.005, gene_size=5)
    cell3 = Cell(n=100, rate=0.006, gene_size=5)
    
    print(f"Cell 1 gene_rate_groups: {cell1.gene_rate_groups}")
    print(f"Cell 2 gene_rate_groups: {cell2.gene_rate_groups}")
    print(f"Cell 3 gene_rate_groups: {cell3.gene_rate_groups}")
    
    # Should fail with static method
    try:
        PetriDish.validate_cells_compatible([cell1, cell2, cell3])
        print("✗ Static validation should have failed for different rates")
        return False
    except ValueError as e:
        print(f"✓ Static validation correctly failed: {str(e)[:100]}...")
    
    # Should fail when creating PetriDish
    try:
        petri = PetriDish(cells=[cell1, cell2], n=100, gene_size=5)
        print("✗ PetriDish creation should have failed for different rates")
        return False
    except ValueError as e:
        print(f"✓ PetriDish creation correctly failed")
    
    return True


def test_identical_gene_groups():
    """Test cells with identical gene-specific rates are compatible."""
    print("\n" + "="*60)
    print("TEST: Identical gene-specific rates")
    print("="*60)
    
    # Create cells with identical gene-specific rates
    groups = [(10, 0.004), (10, 0.006)]
    cell1 = Cell(n=100, gene_rate_groups=groups, gene_size=5)
    cell2 = Cell(n=100, gene_rate_groups=groups, gene_size=5)
    cell3 = Cell(n=100, gene_rate_groups=groups, gene_size=5)
    
    print(f"All cells have gene_rate_groups: {groups}")
    
    # Should work with static method
    try:
        PetriDish.validate_cells_compatible([cell1, cell2, cell3])
        print("✓ Static validation passed for identical gene-specific rates")
    except ValueError as e:
        print(f"✗ Static validation failed: {e}")
        return False
    
    # Should work when creating PetriDish
    try:
        petri = PetriDish(cells=[cell1, cell2, cell3], n=100, gene_size=5)
        print("✓ PetriDish creation passed with identical gene-specific rates")
    except ValueError as e:
        print(f"✗ PetriDish creation failed: {e}")
        return False
    
    return True


def test_different_gene_groups():
    """Test cells with different gene-specific rates are incompatible."""
    print("\n" + "="*60)
    print("TEST: Different gene-specific rates")
    print("="*60)
    
    # Create cells with different gene-specific rates
    groups1 = [(10, 0.004), (10, 0.006)]
    groups2 = [(10, 0.003), (10, 0.007)]
    groups3 = [(5, 0.004), (5, 0.005), (10, 0.006)]  # Different structure
    
    cell1 = Cell(n=100, gene_rate_groups=groups1, gene_size=5)
    cell2 = Cell(n=100, gene_rate_groups=groups2, gene_size=5)
    
    print(f"Cell 1 gene_rate_groups: {groups1}")
    print(f"Cell 2 gene_rate_groups: {groups2}")
    
    # Should fail with static method
    try:
        PetriDish.validate_cells_compatible([cell1, cell2])
        print("✗ Static validation should have failed for different gene groups")
        return False
    except ValueError as e:
        print(f"✓ Static validation correctly failed")
    
    # Test mixing uniform with gene-specific
    cell_uniform = Cell(n=100, rate=0.005, gene_size=5)
    cell_specific = Cell(n=100, gene_rate_groups=groups1, gene_size=5)
    
    print(f"\nMixing uniform rate with gene-specific:")
    print(f"Uniform cell gene_rate_groups: {cell_uniform.gene_rate_groups}")
    print(f"Specific cell gene_rate_groups: {cell_specific.gene_rate_groups}")
    
    try:
        PetriDish.validate_cells_compatible([cell_uniform, cell_specific])
        print("✗ Should have failed mixing uniform with gene-specific")
        return False
    except ValueError as e:
        print(f"✓ Correctly failed mixing rate types")
    
    return True


def test_division_maintains_consistency():
    """Test that cell division preserves rate configuration."""
    print("\n" + "="*60)
    print("TEST: Division maintains consistency")
    print("="*60)
    
    # Start with PetriDish with gene-specific rates
    groups = [(5, 0.003), (10, 0.005), (5, 0.007)]
    petri = PetriDish(gene_rate_groups=groups, n=100, gene_size=5, seed=42)
    
    print(f"Initial cell gene_rate_groups: {petri.cells[0].gene_rate_groups}")
    
    # Divide cells multiple times
    for generation in range(3):
        petri.divide_cells()
        print(f"Generation {generation+1}: {len(petri.cells)} cells")
    
    # Check all cells still have same gene_rate_groups
    try:
        petri.validate_cell_consistency()
        print("✓ All cells maintain same gene_rate_groups after division")
        
        # Verify they actually have the right groups
        for i, cell in enumerate(petri.cells[:3]):  # Check first few
            if cell.gene_rate_groups != groups:
                print(f"✗ Cell {i} has wrong groups: {cell.gene_rate_groups}")
                return False
        print(f"✓ All cells have correct gene_rate_groups: {groups}")
        
    except ValueError as e:
        print(f"✗ Validation failed after division: {e}")
        return False
    
    # Test after culling
    petri.random_cull_cells()
    print(f"After culling: {len(petri.cells)} cells")
    
    try:
        petri.validate_cell_consistency()
        print("✓ Consistency maintained after culling")
    except ValueError as e:
        print(f"✗ Validation failed after culling: {e}")
        return False
    
    return True


def test_edge_cases():
    """Test edge cases like empty lists, single cells, None values."""
    print("\n" + "="*60)
    print("TEST: Edge cases")
    print("="*60)
    
    # Empty list
    try:
        PetriDish.validate_cells_compatible([])
        print("✓ Empty cell list handled correctly")
    except Exception as e:
        print(f"✗ Failed on empty list: {e}")
        return False
    
    # Single cell
    cell = Cell(n=100, rate=0.005, gene_size=5)
    try:
        PetriDish.validate_cells_compatible([cell])
        petri = PetriDish(cells=[cell], n=100, gene_size=5)
        print("✓ Single cell handled correctly")
    except Exception as e:
        print(f"✗ Failed on single cell: {e}")
        return False
    
    # Large number of cells
    cells = [Cell(n=100, rate=0.005, gene_size=5) for _ in range(100)]
    try:
        PetriDish.validate_cells_compatible(cells)
        print("✓ Large cell list (100 cells) validated quickly")
    except Exception as e:
        print(f"✗ Failed on large list: {e}")
        return False
    
    # Different n values (should work - only checking gene_rate_groups)
    cell1 = Cell(n=100, rate=0.005, gene_size=5)
    cell2 = Cell(n=200, rate=0.005, gene_size=5)  # Different n but will have different gene_rate_groups!
    
    print(f"\nCells with different n values:")
    print(f"Cell 1 (n=100): gene_rate_groups={cell1.gene_rate_groups}")
    print(f"Cell 2 (n=200): gene_rate_groups={cell2.gene_rate_groups}")
    
    # These will actually be incompatible because they have different numbers of genes
    try:
        PetriDish.validate_cells_compatible([cell1, cell2])
        print("✗ Should have failed for different n values (different gene counts)")
        return False
    except ValueError:
        print("✓ Correctly failed for cells with different n (leads to different gene_rate_groups)")
    
    return True


def test_error_messages():
    """Test that error messages are informative."""
    print("\n" + "="*60)
    print("TEST: Error message quality")
    print("="*60)
    
    cell1 = Cell(n=100, rate=0.004, gene_size=5)
    cell2 = Cell(n=100, rate=0.006, gene_size=5)
    
    try:
        PetriDish.validate_cells_compatible([cell1, cell2])
    except ValueError as e:
        error_msg = str(e)
        print("Error message received:")
        print(error_msg)
        print()
        
        # Check error message quality
        checks = [
            ("Rate configuration mismatch" in error_msg, "Contains clear error description"),
            ("Cell 0:" in error_msg, "Identifies first cell"),
            ("Cell 1:" in error_msg, "Identifies mismatched cell"),
            ("gene_rate_groups=" in error_msg, "Shows actual values"),
            ("[(20, 0.004)]" in error_msg, "Shows first cell's groups"),
            ("[(20, 0.006)]" in error_msg, "Shows second cell's groups"),
            ("same rate configuration" in error_msg, "Provides guidance")
        ]
        
        for check, description in checks:
            if check:
                print(f"✓ {description}")
            else:
                print(f"✗ {description}")
                
        return all(check for check, _ in checks)
    
    print("✗ Should have raised an error")
    return False


def test_static_method():
    """Test that static method works without PetriDish instance."""
    print("\n" + "="*60)
    print("TEST: Static method functionality")
    print("="*60)
    
    # Create cells
    cells = [Cell(n=100, rate=0.005, gene_size=5) for _ in range(3)]
    
    # Should be callable without instance
    try:
        # Call directly on class
        PetriDish.validate_cells_compatible(cells)
        print("✓ Static method callable on class")
    except Exception as e:
        print(f"✗ Failed to call static method on class: {e}")
        return False
    
    # Should also work on instance (though not typical usage)
    petri = PetriDish(rate=0.005, n=100, gene_size=5)
    try:
        petri.validate_cells_compatible(cells)
        print("✓ Static method callable on instance")
    except Exception as e:
        print(f"✗ Failed to call static method on instance: {e}")
        return False
    
    return True


def main():
    """Run all tests."""
    print("="*60)
    print("RATE CONSISTENCY VALIDATION TESTS")
    print("="*60)
    
    tests = [
        ("Identical uniform rates", test_identical_uniform_rates),
        ("Different uniform rates", test_different_uniform_rates),
        ("Identical gene groups", test_identical_gene_groups),
        ("Different gene groups", test_different_gene_groups),
        ("Division maintains consistency", test_division_maintains_consistency),
        ("Edge cases", test_edge_cases),
        ("Error messages", test_error_messages),
        ("Static method", test_static_method),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            passed = test_func()
            results.append((test_name, passed))
        except Exception as e:
            print(f"\n✗ Test '{test_name}' crashed: {e}")
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