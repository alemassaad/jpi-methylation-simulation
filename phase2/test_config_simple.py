#!/usr/bin/env python3
"""
Simple test for Phase 2 config loading and merging functionality.
Tests the config functions directly without running the full pipeline.
"""

import os
import sys
import tempfile
import yaml
import argparse

# Add path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

# Import the config functions directly
from run_pipeline import load_config, merge_config_and_args, validate_pipeline_config, deep_merge


def test_load_default_config():
    """Test loading default config."""
    print("\n1. Testing load default config...")
    
    config = load_config()
    
    if 'snapshots' in config and 'individuals' in config:
        print(f"✓ Default config loaded successfully")
        print(f"  Snapshots: {config['snapshots']}")
        print(f"  Individuals: {config['individuals']}")
        return True
    else:
        print("✗ Default config missing expected sections")
        return False


def test_load_custom_config():
    """Test loading custom config file."""
    print("\n2. Testing load custom config...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        config_path = os.path.join(tmpdir, 'test.yaml')
        
        test_config = {
            'snapshots': {
                'first': 25,
                'second': 35
            },
            'seed': 999,
            'mixing': {
                'ratio': 70,
                'uniform': True
            }
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(test_config, f)
        
        config = load_config(config_path)
        
        if config.get('snapshots', {}).get('first') == 25:
            print("✓ Custom config loaded successfully")
            print(f"  Loaded snapshots: {config.get('snapshots')}")
            return True
        else:
            print("✗ Custom config values not found")
            return False


def test_deep_merge():
    """Test deep merge functionality."""
    print("\n3. Testing deep merge...")
    
    base = {
        'level1': {
            'level2': {
                'a': 1,
                'b': 2
            },
            'other': 'value'
        },
        'top': 'level'
    }
    
    update = {
        'level1': {
            'level2': {
                'b': 3,
                'c': 4
            }
        },
        'new': 'field'
    }
    
    deep_merge(base, update)
    
    checks = [
        base['level1']['level2']['a'] == 1,  # Original preserved
        base['level1']['level2']['b'] == 3,  # Updated
        base['level1']['level2']['c'] == 4,  # New added
        base['level1']['other'] == 'value',  # Other preserved
        base['top'] == 'level',              # Top preserved
        base['new'] == 'field'                # New top-level added
    ]
    
    if all(checks):
        print("✓ Deep merge works correctly")
        return True
    else:
        print("✗ Deep merge failed")
        print(f"  Result: {base}")
        return False


def test_merge_config_and_args():
    """Test merging config with command-line arguments."""
    print("\n4. Testing merge config and args...")
    
    # Create mock args
    args = argparse.Namespace(
        rate=None,
        gene_rate_groups=None,
        gene_size=5,
        simulation=None,
        first_snapshot=50,  # Default
        second_snapshot=60,  # Default
        individual_growth_phase=7,  # Default
        n_quantiles=10,  # Default
        cells_per_quantile=3,  # Default
        mix_ratio=80,  # Default
        bins=200,  # Default
        uniform_mixing=False,
        normalize_size=False,
        plot_individuals=False,
        no_compress=False,
        seed=42,  # Default
        force_reload=False
    )
    
    config = {
        'input': {
            'rate': 0.006,
            'simulation': 'test.json.gz'
        },
        'snapshots': {
            'first': 30,
            'second': 40
        },
        'individuals': {
            'growth_phase': 5,
            'n_quantiles': 5
        },
        'seed': 123
    }
    
    # Merge config into args
    merged = merge_config_and_args(config, args)
    
    # Check that config values were applied to defaults
    checks = [
        merged.rate == 0.006,  # From config
        merged.simulation == 'test.json.gz',  # From config
        merged.first_snapshot == 30,  # From config (was default)
        merged.second_snapshot == 40,  # From config (was default)
        merged.individual_growth_phase == 5,  # From config
        merged.n_quantiles == 5,  # From config
        merged.seed == 123,  # From config
        merged.cells_per_quantile == 3  # Default preserved
    ]
    
    if all(checks):
        print("✓ Config merged with args correctly")
        return True
    else:
        print("✗ Merge failed")
        print(f"  Rate: {merged.rate} (expected 0.006)")
        print(f"  First snapshot: {merged.first_snapshot} (expected 30)")
        print(f"  Seed: {merged.seed} (expected 123)")
        return False


def test_cli_override():
    """Test that CLI args override config values."""
    print("\n5. Testing CLI override...")
    
    # Create args with non-default values
    args = argparse.Namespace(
        rate=0.008,  # CLI override
        gene_rate_groups=None,
        gene_size=5,
        simulation='cli.json.gz',  # CLI override
        first_snapshot=45,  # CLI override
        second_snapshot=60,  # Default
        individual_growth_phase=7,  # Default
        n_quantiles=10,  # Default
        cells_per_quantile=3,  # Default
        mix_ratio=80,  # Default
        bins=200,  # Default
        uniform_mixing=True,  # CLI override
        normalize_size=False,
        plot_individuals=False,
        no_compress=False,
        seed=500,  # CLI override
        force_reload=False
    )
    
    config = {
        'input': {
            'rate': 0.006,
            'simulation': 'config.json.gz'
        },
        'snapshots': {
            'first': 30,
            'second': 40
        },
        'seed': 123,
        'mixing': {
            'uniform': False
        }
    }
    
    # Merge - CLI should win
    merged = merge_config_and_args(config, args)
    
    checks = [
        merged.rate == 0.008,  # CLI wins
        merged.simulation == 'cli.json.gz',  # CLI wins
        merged.first_snapshot == 45,  # CLI wins
        merged.second_snapshot == 40,  # Config (CLI was default)
        merged.seed == 500,  # CLI wins
        merged.uniform_mixing == True  # CLI wins
    ]
    
    if all(checks):
        print("✓ CLI arguments override config correctly")
        return True
    else:
        print("✗ CLI override failed")
        print(f"  Rate: {merged.rate} (expected 0.008 from CLI)")
        print(f"  Simulation: {merged.simulation} (expected cli.json.gz)")
        return False


def test_validation():
    """Test config validation."""
    print("\n6. Testing validation...")
    
    # Test valid config
    valid_args = argparse.Namespace(
        rate=0.005,
        gene_rate_groups=None,
        simulation='test.json.gz',
        first_snapshot=50,
        second_snapshot=60,
        mix_ratio=80,
        n_quantiles=10,
        cells_per_quantile=3,
        bins=200
    )
    
    try:
        validate_pipeline_config(valid_args)
        print("✓ Valid config passes validation")
        valid_passed = True
    except Exception as e:
        print(f"✗ Valid config failed: {e}")
        valid_passed = False
    
    # Test invalid config - missing rate
    invalid_args1 = argparse.Namespace(
        rate=None,
        gene_rate_groups=None,
        simulation='test.json.gz',
        first_snapshot=50,
        second_snapshot=60,
        mix_ratio=80,
        n_quantiles=10,
        cells_per_quantile=3,
        bins=200
    )
    
    try:
        validate_pipeline_config(invalid_args1)
        print("✗ Should have caught missing rate")
        invalid1_passed = False
    except ValueError as e:
        if "rate" in str(e).lower():
            print("✓ Caught missing rate")
            invalid1_passed = True
        else:
            print(f"✗ Wrong error for missing rate: {e}")
            invalid1_passed = False
    
    # Test invalid config - bad mix ratio
    invalid_args2 = argparse.Namespace(
        rate=0.005,
        gene_rate_groups=None,
        simulation='test.json.gz',
        first_snapshot=50,
        second_snapshot=60,
        mix_ratio=150,  # Invalid
        n_quantiles=10,
        cells_per_quantile=3,
        bins=200
    )
    
    try:
        validate_pipeline_config(invalid_args2)
        print("✗ Should have caught invalid mix ratio")
        invalid2_passed = False
    except ValueError as e:
        if "mix" in str(e).lower():
            print("✓ Caught invalid mix ratio")
            invalid2_passed = True
        else:
            print(f"✗ Wrong error for mix ratio: {e}")
            invalid2_passed = False
    
    return valid_passed and invalid1_passed and invalid2_passed


def test_gene_rate_groups():
    """Test gene rate group configuration."""
    print("\n7. Testing gene rate groups...")
    
    args = argparse.Namespace(
        rate=None,
        gene_rate_groups=None,
        gene_size=5,
        simulation=None,
        first_snapshot=50,
        second_snapshot=60,
        individual_growth_phase=7,
        n_quantiles=10,
        cells_per_quantile=3,
        mix_ratio=80,
        bins=200,
        uniform_mixing=False,
        normalize_size=False,
        plot_individuals=False,
        no_compress=False,
        seed=42,
        force_reload=False
    )
    
    config = {
        'input': {
            'gene_rate_groups': '50:0.004,50:0.005,50:0.006,50:0.007',
            'gene_size': 5,
            'simulation': 'gene_test.json.gz'
        }
    }
    
    merged = merge_config_and_args(config, args)
    
    if merged.gene_rate_groups == '50:0.004,50:0.005,50:0.006,50:0.007':
        print("✓ Gene rate groups loaded from config")
        return True
    else:
        print("✗ Gene rate groups not loaded")
        return False


def test_compression_flag():
    """Test compression flag handling."""
    print("\n8. Testing compression flag...")
    
    args = argparse.Namespace(
        rate=0.005,
        gene_rate_groups=None,
        gene_size=5,
        simulation='test.json.gz',
        first_snapshot=50,
        second_snapshot=60,
        individual_growth_phase=7,
        n_quantiles=10,
        cells_per_quantile=3,
        mix_ratio=80,
        bins=200,
        uniform_mixing=False,
        normalize_size=False,
        plot_individuals=False,
        no_compress=False,  # Default: compress
        seed=42,
        force_reload=False
    )
    
    # Config says don't compress
    config = {
        'output': {
            'compress': False
        }
    }
    
    merged = merge_config_and_args(config, args)
    
    # Config compress=False should set no_compress=True
    if merged.no_compress == True:
        print("✓ Compression flag converted correctly")
        return True
    else:
        print("✗ Compression flag not converted")
        print(f"  no_compress: {merged.no_compress} (expected True)")
        return False


def main():
    """Run all tests."""
    print("=" * 60)
    print("PHASE 2 CONFIG FUNCTIONALITY TEST")
    print("=" * 60)
    
    # Check YAML is available
    try:
        import yaml
        print("✓ YAML module available")
    except ImportError:
        print("✗ YAML module not available")
        print("  Install with: pip install pyyaml")
        return 1
    
    tests = [
        test_load_default_config,
        test_load_custom_config,
        test_deep_merge,
        test_merge_config_and_args,
        test_cli_override,
        test_validation,
        test_gene_rate_groups,
        test_compression_flag
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"✗ Test crashed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)
    
    # Also test that example configs are valid YAML
    print("\nChecking example config files...")
    example_configs = [
        'configs/quick_test.yaml',
        'configs/debug.yaml',
        'configs/uniform_mixing.yaml',
        'configs/full_analysis.yaml'
    ]
    
    for config_file in example_configs:
        if os.path.exists(config_file):
            try:
                with open(config_file, 'r') as f:
                    yaml.safe_load(f)
                print(f"✓ {config_file} is valid YAML")
            except Exception as e:
                print(f"✗ {config_file} has YAML error: {e}")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())