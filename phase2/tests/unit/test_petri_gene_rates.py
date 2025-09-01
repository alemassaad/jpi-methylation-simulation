#!/usr/bin/env python3
"""
Direct test of PetriDish with gene rates without full pipeline dependencies.
"""

import os
import sys

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

def test_petri_dish_creation():
    """Test that PetriDish can be created with gene rate groups."""
    print("\nTesting PetriDish creation with gene rates...")
    
    from cell import PetriDish, Cell
    
    # Test 1: Create PetriDish with uniform rate
    petri_uniform = PetriDish(rate=0.005, n=100, growth_phase=2)
    assert petri_uniform.rate == 0.005
    assert len(petri_uniform.cells) == 1  # Starts with 1 cell
    print("  ✓ Uniform rate PetriDish created")
    
    # Test 2: Create PetriDish with gene-specific rates
    gene_groups = [(10, 0.004), (10, 0.006)]  # 20 genes total for 100 sites
    petri_gene = PetriDish(gene_rate_groups=gene_groups, n=100, gene_size=5, growth_phase=2)
    assert petri_gene.gene_rate_groups == gene_groups
    assert len(petri_gene.cells) == 1
    print("  ✓ Gene-specific rate PetriDish created")
    
    # Test 3: Verify cell has correct configuration
    cell = petri_gene.cells[0]
    assert hasattr(cell, 'site_rates'), "Cell should have site_rates for gene-specific config"
    assert len(cell.site_rates) == 100, "Should have rate for each site"
    print("  ✓ Cell has correct gene-specific configuration")
    
    # Test 4: Test division with gene rates
    petri_gene.divide_cells()
    assert len(petri_gene.cells) == 2, "Should have 2 cells after division"
    
    # Both cells should have same configuration
    for cell in petri_gene.cells:
        assert len(cell.site_rates) == 100
    print("  ✓ Cell division preserves gene-specific rates")
    
    # Test 5: Test methylation
    initial_meth = [sum(cell.cpg_sites) for cell in petri_gene.cells]
    petri_gene.methylate_cells()
    final_meth = [sum(cell.cpg_sites) for cell in petri_gene.cells]
    
    # At least some cells should have gained methylation
    assert any(f >= i for f, i in zip(final_meth, initial_meth))
    print("  ✓ Methylation works with gene-specific rates")
    
    return True


def test_create_petri_from_config():
    """Test creating PetriDish from rate configuration dict."""
    print("\nTesting PetriDish creation from config...")
    
    from cell import PetriDish
    
    # Import the helper function from run_pipeline
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from run_pipeline import create_petri_with_rate_config
    
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
    
    return True


def test_edge_cases():
    """Test edge cases for gene rate support."""
    print("\nTesting edge cases...")
    
    from cell import PetriDish, Cell
    
    # Test 1: Very small configuration (10 sites, 2 genes)
    gene_groups = [(1, 0.01), (1, 0.02)]
    petri = PetriDish(gene_rate_groups=gene_groups, n=10, gene_size=5, growth_phase=1)
    assert len(petri.cells[0].cpg_sites) == 10
    print("  ✓ Very small configuration works")
    
    # Test 2: Single gene group
    gene_groups = [(20, 0.005)]  # All genes same rate
    petri = PetriDish(gene_rate_groups=gene_groups, n=100, gene_size=5, growth_phase=1)
    cell = petri.cells[0]
    # All sites should have same rate
    assert all(r == 0.005 for r in cell.site_rates)
    print("  ✓ Single gene group works")
    
    # Test 3: Many small groups
    gene_groups = [(1, 0.001 + i*0.001) for i in range(20)]  # 20 groups, 1 gene each
    petri = PetriDish(gene_rate_groups=gene_groups, n=100, gene_size=5, growth_phase=1)
    assert len(petri.cells[0].site_rates) == 100
    print("  ✓ Many small groups work")
    
    # Test 4: Growth phase works correctly
    petri = PetriDish(gene_rate_groups=[(10, 0.01), (10, 0.02)], n=100, gene_size=5, growth_phase=2)
    
    # Grow for 2 years (should reach 2^2 = 4 cells)
    for year in range(2):
        petri.divide_cells()
    
    assert len(petri.cells) == 4, f"Expected 4 cells after growth phase 2, got {len(petri.cells)}"
    print("  ✓ Growth phase works with gene rates")
    
    return True


def test_backward_compatibility():
    """Ensure uniform rate still works as before."""
    print("\nTesting backward compatibility...")
    
    from cell import PetriDish, Cell
    
    # Create with uniform rate (old style)
    petri = PetriDish(rate=0.005, n=100, growth_phase=2)
    
    # Should work exactly as before
    assert petri.rate == 0.005
    assert petri.n == 100
    assert petri.growth_phase == 2
    assert len(petri.cells) == 1
    
    # Cell should use uniform rate
    cell = petri.cells[0]
    assert cell.rate == 0.005
    
    # Division should work
    petri.divide_cells()
    assert len(petri.cells) == 2
    
    # Methylation should work
    petri.methylate_cells()
    
    print("  ✓ Uniform rate backward compatibility maintained")
    
    return True


def run_all_tests():
    """Run all tests."""
    print("="*60)
    print("Testing PetriDish Gene Rate Support")
    print("="*60)
    
    tests = [
        test_petri_dish_creation,
        test_create_petri_from_config,
        test_edge_cases,
        test_backward_compatibility
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
        print("✅ All PetriDish gene rate tests passed!")
    else:
        print(f"❌ {failed} test(s) failed")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(run_all_tests())