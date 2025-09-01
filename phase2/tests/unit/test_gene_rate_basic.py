#!/usr/bin/env python3
"""
Basic tests for gene rate support without numpy dependency.
"""

import os
import sys

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

def test_parse_gene_rate_groups():
    """Test parsing of gene rate groups string."""
    print("\nTesting parse_gene_rate_groups...")
    
    # Import here to test the function exists
    from run_pipeline import parse_gene_rate_groups
    
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
    
    return True


def test_validate_gene_rate_groups():
    """Test validation of gene rate groups."""
    print("\nTesting validate_gene_rate_groups...")
    
    from run_pipeline import validate_gene_rate_groups
    
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
    
    return True


def test_rate_config_creation():
    """Test rate configuration dictionary creation."""
    print("\nTesting rate configuration creation...")
    
    # Test uniform rate config
    rate_config = {
        'type': 'uniform',
        'rate': 0.005
    }
    assert rate_config['type'] == 'uniform'
    assert rate_config['rate'] == 0.005
    print("  ✓ Uniform rate config created")
    
    # Test gene-specific rate config
    rate_config = {
        'type': 'gene_specific',
        'gene_rate_groups': [(50, 0.004), (50, 0.005), (50, 0.006), (50, 0.007)],
        'gene_size': 5
    }
    assert rate_config['type'] == 'gene_specific'
    assert len(rate_config['gene_rate_groups']) == 4
    assert rate_config['gene_size'] == 5
    print("  ✓ Gene-specific rate config created")
    
    return True


def test_path_generation():
    """Test that path generation handles gene rate groups."""
    print("\nTesting path generation with gene rates...")
    
    # Create a mock args object
    class MockArgs:
        def __init__(self):
            self.gene_rate_groups = "50:0.004,50:0.005,50:0.006,50:0.007"
            self.rate = None
            self.first_snapshot = 50
            self.second_snapshot = 60
            self.individual_growth_phase = 7
            self.n_quantiles = 10
            self.cells_per_quantile = 3
            self.mix_ratio = 80
            self.seed = 42
            self.output_dir = "data"
            self.uniform_mixing = False
            self.normalize_size = False
    
    args = MockArgs()
    sim_params = {
        'growth_phase': 13,
        'n_sites': 1000,
        'sim_years': 100
    }
    
    # Import and test path generation
    from path_utils import generate_step23_output_dir
    
    output_dir = generate_step23_output_dir(args, sim_params)
    
    # Check that gene_rates is in the path
    assert "gene_rates" in output_dir, f"Expected 'gene_rates' in path: {output_dir}"
    assert "50x0.00400" in output_dir, f"Expected '50x0.00400' in path: {output_dir}"
    print(f"  ✓ Gene rates path generated: {output_dir}")
    
    # Test with uniform rate
    args.gene_rate_groups = None
    args.rate = 0.005
    output_dir = generate_step23_output_dir(args, sim_params)
    assert "rate_0.00500" in output_dir, f"Expected 'rate_0.00500' in path: {output_dir}"
    print(f"  ✓ Uniform rate path generated: {output_dir}")
    
    return True


def run_all_tests():
    """Run all basic tests."""
    print("="*60)
    print("Running Basic Gene Rate Support Tests")
    print("="*60)
    
    tests = [
        test_parse_gene_rate_groups,
        test_validate_gene_rate_groups,
        test_rate_config_creation,
        test_path_generation
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
        print("✅ All basic tests passed!")
        print("\nNote: Full integration tests require numpy/scipy installation")
    else:
        print(f"❌ {failed} test(s) failed")
        return 1
    
    return 0


if __name__ == "__main__":
    sys.exit(run_all_tests())