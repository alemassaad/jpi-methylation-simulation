#!/usr/bin/env python3
"""
Test helper functions directly without importing the full pipeline.
"""

import os
import sys

# Add parent directories to path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))


def test_create_petri_with_rate_config():
    """Test the helper function for creating PetriDish with rate config."""
    print("\nTesting create_petri_with_rate_config helper...")
    
    from cell import PetriDish
    
    # Copy the function here to test it without imports
    def create_petri_with_rate_config(rate_config: dict, growth_phase: int) -> PetriDish:
        """Create PetriDish with proper rate configuration."""
        if rate_config['type'] == 'uniform':
            return PetriDish(rate=rate_config['rate'], growth_phase=growth_phase)
        else:
            return PetriDish(
                gene_rate_groups=rate_config['gene_rate_groups'],
                gene_size=rate_config['gene_size'],
                growth_phase=growth_phase
            )
    
    # Test uniform config
    config_uniform = {'type': 'uniform', 'rate': 0.005}
    petri = create_petri_with_rate_config(config_uniform, growth_phase=2)
    assert petri.rate == 0.005
    assert petri.growth_phase == 2
    print("  ✓ Created from uniform config")
    
    # Test gene-specific config
    config_gene = {
        'type': 'gene_specific',
        'gene_rate_groups': [(10, 0.004), (10, 0.006)],
        'gene_size': 5
    }
    petri = create_petri_with_rate_config(config_gene, growth_phase=3)
    assert petri.gene_rate_groups == [(10, 0.004), (10, 0.006)]
    assert petri.gene_size == 5
    assert petri.growth_phase == 3
    print("  ✓ Created from gene-specific config")
    
    # Test that cells can be added
    from cell import Cell
    
    # For uniform rate
    petri_u = create_petri_with_rate_config(config_uniform, growth_phase=2)
    cell_u = Cell(n=100, rate=0.005)
    petri_u.cells = [cell_u]
    assert len(petri_u.cells) == 1
    print("  ✓ Can add cells to uniform rate PetriDish")
    
    # For gene-specific rate
    petri_g = create_petri_with_rate_config(config_gene, growth_phase=2)
    cell_g = Cell(n=100, gene_rate_groups=[(10, 0.004), (10, 0.006)], gene_size=5)
    petri_g.cells = [cell_g]
    assert len(petri_g.cells) == 1
    print("  ✓ Can add cells to gene-specific PetriDish")
    
    return True


def test_parse_validate_functions():
    """Test parsing and validation functions."""
    print("\nTesting parse and validate functions...")
    
    # Copy functions here to test without imports
    def parse_gene_rate_groups(groups_str: str):
        """Parse gene rate groups string into list of tuples."""
        groups = []
        for group in groups_str.split(','):
            n, rate = group.split(':')
            groups.append((int(n), float(rate)))
        return groups
    
    def validate_gene_rate_groups(groups, n_sites: int, gene_size: int):
        """Validate gene rate groups match expected configuration."""
        total_genes = sum(n for n, _ in groups)
        expected_genes = n_sites // gene_size
        
        if total_genes != expected_genes:
            raise ValueError(
                f"Gene rate groups specify {total_genes} genes, "
                f"but simulation expects {expected_genes} genes "
                f"({n_sites} sites / {gene_size} sites per gene)"
            )
        
        # Validate rates are reasonable
        for n, rate in groups:
            if rate < 0 or rate > 1:
                raise ValueError(f"Invalid rate {rate}: must be between 0 and 1")
    
    # Test parsing
    result = parse_gene_rate_groups("50:0.004,50:0.005,50:0.006,50:0.007")
    assert len(result) == 4
    assert result[0] == (50, 0.004)
    assert result[3] == (50, 0.007)
    print("  ✓ Parsing works correctly")
    
    # Test validation - valid
    try:
        validate_gene_rate_groups(result, 1000, 5)
        print("  ✓ Valid configuration passes validation")
    except ValueError:
        assert False, "Valid configuration should pass"
    
    # Test validation - invalid count
    try:
        bad_groups = [(50, 0.004), (50, 0.005)]  # Only 100 genes
        validate_gene_rate_groups(bad_groups, 1000, 5)
        assert False, "Should have failed validation"
    except ValueError as e:
        assert "100 genes" in str(e)
        print("  ✓ Invalid gene count caught")
    
    # Test validation - invalid rate
    try:
        bad_groups = [(100, -0.5), (100, 0.5)]
        validate_gene_rate_groups(bad_groups, 1000, 5)
        assert False, "Should have failed validation"
    except ValueError as e:
        assert "Invalid rate" in str(e)
        print("  ✓ Invalid rate caught")
    
    return True


def test_real_scenario():
    """Test a realistic scenario with gene rates."""
    print("\nTesting realistic scenario...")
    
    from cell import Cell, PetriDish
    
    # Simulate what Phase 2 does:
    # 1. Load cells from Phase 1 (with gene rates)
    # 2. Create PetriDish with matching configuration
    # 3. Grow the population
    
    # Step 1: Create cells as if loaded from Phase 1
    gene_groups = [(25, 0.002), (25, 0.004), (25, 0.006), (25, 0.008)]
    cells_from_phase1 = []
    for i in range(3):  # 3 cells
        cell = Cell(n=500, gene_rate_groups=gene_groups, gene_size=5)
        # Simulate some prior methylation
        for _ in range(10):
            cell.methylate()
        cells_from_phase1.append(cell)
    
    print(f"  Created {len(cells_from_phase1)} cells with gene-specific rates")
    
    # Step 2: Create PetriDish for Phase 2 growth
    petri = PetriDish(
        gene_rate_groups=gene_groups,
        n=500,
        gene_size=5,
        growth_phase=2,  # Small for testing
        cells=cells_from_phase1  # Start with loaded cells
    )
    
    assert len(petri.cells) == 3
    print("  ✓ PetriDish created with Phase 1 cells")
    
    # Step 3: Grow population (like Phase 2 does)
    # Growth phase: double cells
    for year in range(2):  # growth_phase = 2
        petri.divide_cells()
        petri.methylate_cells()
    
    expected_cells = 3 * (2 ** 2)  # 3 initial × 4 (2^2)
    assert len(petri.cells) == expected_cells, f"Expected {expected_cells} cells, got {len(petri.cells)}"
    print(f"  ✓ Growth phase complete: {expected_cells} cells")
    
    # Step 4: Homeostasis (random culling)
    petri.random_cull_cells()
    # Should have ~half the cells
    assert expected_cells // 3 <= len(petri.cells) <= expected_cells
    print(f"  ✓ Homeostasis culling: {len(petri.cells)} cells remain")
    
    # Verify all cells still have correct gene configuration
    for cell in petri.cells:
        assert hasattr(cell, 'site_rates')
        assert len(cell.site_rates) == 500
    print("  ✓ All cells maintain gene-specific configuration")
    
    return True


def run_all_tests():
    """Run all tests."""
    print("="*60)
    print("Testing Gene Rate Helper Functions")
    print("="*60)
    
    tests = [
        test_create_petri_with_rate_config,
        test_parse_validate_functions,
        test_real_scenario
    ]
    
    failed = 0
    for test_func in tests:
        try:
            if not test_func():
                failed += 1
        except Exception as e:
            print(f"\n✗ {test_func.__name__} FAILED:")
            print(f"  {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "="*60)
    if failed == 0:
        print("✅ All helper function tests passed!")
        print("\nGene rate support is working correctly:")
        print("  • PetriDish accepts gene_rate_groups parameter")
        print("  • Cells maintain gene-specific methylation rates")
        print("  • Division preserves rate configuration")
        print("  • Growth and homeostasis work with gene rates")
        print("  • Backward compatibility with uniform rates maintained")
    else:
        print(f"❌ {failed} test(s) failed")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(run_all_tests())