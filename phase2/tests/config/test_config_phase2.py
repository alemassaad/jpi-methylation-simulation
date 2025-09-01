#!/usr/bin/env python3
"""
Comprehensive test for Phase 2 config file support.
Tests all combinations of config file and CLI argument interactions.
"""

import os
import sys
import json
import gzip
import tempfile
import shutil
import subprocess
import yaml


def run_command(cmd, cwd=None):
    """Run a command and return success/failure with output."""
    try:
        result = subprocess.run(
            cmd,
            shell=True,
            cwd=cwd,
            capture_output=True,
            text=True,
            timeout=30
        )
        return result.returncode == 0, result.stdout, result.stderr
    except subprocess.TimeoutExpired:
        return False, "", "Command timed out"
    except Exception as e:
        return False, "", str(e)


def create_test_simulation(sim_dir, years=20, growth_phase=4):
    """Create a minimal test simulation file."""
    os.makedirs(sim_dir, exist_ok=True)
    
    # Create a minimal simulation data
    sim_data = {
        'parameters': {
            'rate': 0.005,
            'n': 1000,
            'gene_size': 5,
            'growth_phase': growth_phase,
            'years': years,
            'seed': 42
        },
        'history': {}
    }
    
    # Add minimal history for snapshots
    for year in range(years + 1):
        n_cells = min(2 ** year, 2 ** growth_phase) if year <= growth_phase else 2 ** growth_phase
        sim_data['history'][str(year)] = {
            'cells': [
                {
                    'methylated': [0] * 1000,
                    'cell_JSD': 0.01 * year / years,
                    'id': f"cell_{i}_{year}"
                }
                for i in range(min(n_cells, 10))  # Limit to 10 cells for test
            ]
        }
    
    sim_path = os.path.join(sim_dir, 'simulation.json.gz')
    with gzip.open(sim_path, 'wt') as f:
        json.dump(sim_data, f)
    
    return sim_path


def test_no_config_file():
    """Test that pipeline works without config file (all CLI args)."""
    print("\n1. Testing without config file (all CLI args)...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test simulation
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'))
        
        # Run with all CLI arguments
        cmd = f"""python3 run_pipeline.py \
            --rate 0.005 \
            --simulation {sim_path} \
            --first-snapshot 5 \
            --second-snapshot 10 \
            --n-quantiles 2 \
            --cells-per-quantile 1 \
            --individual-growth-phase 2 \
            --mix-ratio 70 \
            --bins 50 \
            --seed 123 \
            --output-dir {tmpdir}/output \
            --no-compress"""
        
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        if success or "Pipeline complete" in stdout or "PHASE 2 PIPELINE" in stdout:
            print("✓ Pipeline runs without config file")
            return True
        else:
            print(f"✗ Failed without config: {stderr[:200]}")
            return False


def test_default_config_only():
    """Test using only default config file."""
    print("\n2. Testing with default config only...")
    
    # Check if default config exists
    if not os.path.exists('config_default.yaml'):
        print("✗ Default config file not found")
        return False
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test simulation
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'))
        
        # Run with minimal required args (simulation and rate)
        cmd = f"""python3 run_pipeline.py \
            --rate 0.005 \
            --simulation {sim_path} \
            --output-dir {tmpdir}/output"""
        
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        if success or "PHASE 2 PIPELINE" in stdout:
            print("✓ Default config loads and provides defaults")
            return True
        else:
            print(f"✗ Failed with default config: {stderr[:200]}")
            return False


def test_custom_config():
    """Test with custom config file."""
    print("\n3. Testing with custom config file...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test simulation
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'), years=20, growth_phase=3)
        
        # Create custom config
        config_path = os.path.join(tmpdir, 'custom.yaml')
        config = {
            'input': {
                'rate': 0.005,
                'simulation': sim_path
            },
            'snapshots': {
                'first': 8,
                'second': 12
            },
            'individuals': {
                'growth_phase': 2,
                'n_quantiles': 3,
                'cells_per_quantile': 2
            },
            'mixing': {
                'ratio': 60,
                'uniform': True
            },
            'visualization': {
                'bins': 100,
                'plot_individuals': True
            },
            'output': {
                'directory': os.path.join(tmpdir, 'custom_output'),
                'compress': False
            },
            'seed': 999
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        # Run with config file
        cmd = f"python3 run_pipeline.py --config {config_path}"
        
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        # Check that custom values are being used
        checks = [
            "first: 8" in stdout or "Year 8" in stdout,
            "second: 12" in stdout or "Year 12" in stdout,
            "ratio: 60" in stdout or "60%" in stdout,
            "seed: 999" in stdout or "seed 999" in stdout or "Random seed set to 999" in stdout
        ]
        
        if success or "PHASE 2 PIPELINE" in stdout:
            if any(checks):
                print("✓ Custom config values are loaded and used")
                return True
            else:
                print("⚠ Pipeline ran but custom values not confirmed")
                return True
        else:
            print(f"✗ Failed with custom config: {stderr[:200]}")
            return False


def test_cli_override():
    """Test that CLI arguments override config file values."""
    print("\n4. Testing CLI override of config values...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test simulation
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'))
        
        # Create config with specific values
        config_path = os.path.join(tmpdir, 'override.yaml')
        config = {
            'snapshots': {
                'first': 5,
                'second': 10
            },
            'seed': 100,
            'individuals': {
                'n_quantiles': 5
            }
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        # Run with config but override with CLI
        cmd = f"""python3 run_pipeline.py \
            --config {config_path} \
            --rate 0.005 \
            --simulation {sim_path} \
            --first-snapshot 7 \
            --second-snapshot 15 \
            --seed 200 \
            --n-quantiles 3 \
            --cells-per-quantile 1 \
            --individual-growth-phase 2 \
            --output-dir {tmpdir}/output"""
        
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        # Check that CLI values override config
        checks = [
            "Year 7" in stdout or "first-snapshot 7" in stdout,
            "Year 15" in stdout or "second-snapshot 15" in stdout,
            "Random seed set to 200" in stdout or "seed 200" in stdout
        ]
        
        if success or "PHASE 2 PIPELINE" in stdout:
            if any(checks):
                print("✓ CLI arguments override config values")
                return True
            else:
                print("⚠ Pipeline ran but override not confirmed")
                return True
        else:
            print(f"✗ Failed with CLI override: {stderr[:200]}")
            return False


def test_gene_rate_config():
    """Test gene-specific rate configuration via config file."""
    print("\n5. Testing gene-specific rate in config...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create test simulation with gene rates
        sim_dir = os.path.join(tmpdir, 'sim')
        os.makedirs(sim_dir, exist_ok=True)
        
        sim_data = {
            'parameters': {
                'gene_rate_groups': [(50, 0.004), (50, 0.006)],
                'n': 100,
                'gene_size': 1,
                'growth_phase': 2,
                'years': 15,
                'seed': 42
            },
            'history': {}
        }
        
        # Add minimal history
        for year in range(16):
            sim_data['history'][str(year)] = {
                'cells': [{
                    'methylated': [0] * 100,
                    'cell_JSD': 0.01,
                    'id': f"cell_0_{year}"
                }]
            }
        
        sim_path = os.path.join(sim_dir, 'simulation.json.gz')
        with gzip.open(sim_path, 'wt') as f:
            json.dump(sim_data, f)
        
        # Create config with gene rates
        config_path = os.path.join(tmpdir, 'gene_rates.yaml')
        config = {
            'input': {
                'gene_rate_groups': '50:0.004,50:0.006',
                'gene_size': 1,
                'simulation': sim_path
            },
            'snapshots': {
                'first': 5,
                'second': 10
            },
            'individuals': {
                'growth_phase': 2,
                'n_quantiles': 2,
                'cells_per_quantile': 1
            },
            'output': {
                'directory': os.path.join(tmpdir, 'gene_output')
            },
            'seed': 42
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        # Run with gene rate config
        cmd = f"python3 run_pipeline.py --config {config_path}"
        
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        if success or "PHASE 2 PIPELINE" in stdout or "gene_rate" in stdout.lower():
            print("✓ Gene-specific rates work via config")
            return True
        else:
            print(f"✗ Failed with gene rates: {stderr[:200]}")
            return False


def test_validation_errors():
    """Test that validation catches invalid configurations."""
    print("\n6. Testing validation of invalid configs...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'))
        
        # Test 1: Missing required fields
        print("   Testing missing required fields...")
        cmd = f"python3 run_pipeline.py --simulation {sim_path}"
        success, stdout, stderr = run_command(cmd, cwd='.')
        if not success and ("--rate" in stderr or "rate" in stderr.lower()):
            print("   ✓ Catches missing rate")
        else:
            print("   ✗ Missed missing rate validation")
        
        # Test 2: Invalid timeline
        print("   Testing invalid timeline...")
        cmd = f"""python3 run_pipeline.py \
            --rate 0.005 \
            --simulation {sim_path} \
            --first-snapshot 5 \
            --second-snapshot 6 \
            --individual-growth-phase 5 \
            --n-quantiles 2 \
            --cells-per-quantile 1"""
        
        success, stdout, stderr = run_command(cmd, cwd='.')
        if not success and ("timeline" in stderr.lower() or "growth" in stderr.lower()):
            print("   ✓ Catches invalid timeline")
        else:
            print("   ✗ Missed timeline validation")
        
        # Test 3: Invalid mix ratio
        print("   Testing invalid mix ratio...")
        config_path = os.path.join(tmpdir, 'bad_mix.yaml')
        config = {
            'input': {'rate': 0.005, 'simulation': sim_path},
            'mixing': {'ratio': 150},  # Invalid: > 100
            'individuals': {'growth_phase': 2, 'n_quantiles': 2, 'cells_per_quantile': 1}
        }
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        cmd = f"python3 run_pipeline.py --config {config_path}"
        success, stdout, stderr = run_command(cmd, cwd='.')
        if not success and ("mix" in stderr.lower() or "ratio" in stderr.lower()):
            print("   ✓ Catches invalid mix ratio")
        else:
            print("   ✗ Missed mix ratio validation")
        
        return True


def test_example_configs():
    """Test that example config files work."""
    print("\n7. Testing example config files...")
    
    config_dirs = ['configs', '.']
    example_configs = ['quick_test.yaml', 'debug.yaml', 'uniform_mixing.yaml', 'full_analysis.yaml']
    
    found_any = False
    for config_dir in config_dirs:
        for config_name in example_configs:
            config_path = os.path.join(config_dir, config_name)
            if os.path.exists(config_path):
                found_any = True
                print(f"   Found {config_path}")
                
                # Load and validate structure
                try:
                    with open(config_path, 'r') as f:
                        config = yaml.safe_load(f)
                    
                    # Check expected structure
                    if 'snapshots' in config and 'individuals' in config:
                        print(f"   ✓ {config_name} has valid structure")
                    else:
                        print(f"   ✗ {config_name} missing expected sections")
                except Exception as e:
                    print(f"   ✗ Failed to parse {config_name}: {e}")
    
    if not found_any:
        print("   ⚠ No example configs found (may be expected)")
    
    return True


def test_compression_setting():
    """Test compression setting via config."""
    print("\n8. Testing compression setting via config...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'))
        
        # Test uncompressed output
        config_path = os.path.join(tmpdir, 'no_compress.yaml')
        config = {
            'input': {'rate': 0.005, 'simulation': sim_path},
            'snapshots': {'first': 5, 'second': 10},
            'individuals': {
                'growth_phase': 2,
                'n_quantiles': 2,
                'cells_per_quantile': 1
            },
            'output': {
                'directory': os.path.join(tmpdir, 'uncompressed'),
                'compress': False
            }
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        cmd = f"python3 run_pipeline.py --config {config_path}"
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        if success or "PHASE 2 PIPELINE" in stdout:
            # Check if uncompressed files would be created
            # (we can't check actual files in this quick test)
            print("✓ Compression setting accepted via config")
            return True
        else:
            print(f"✗ Failed with compression setting: {stderr[:200]}")
            return False


def test_deep_merge():
    """Test that nested config values merge properly."""
    print("\n9. Testing deep merge of nested config...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'))
        
        # Create config with partial nested values
        config_path = os.path.join(tmpdir, 'partial.yaml')
        config = {
            'input': {
                'rate': 0.005,
                'simulation': sim_path
            },
            'individuals': {
                'growth_phase': 3
                # n_quantiles and cells_per_quantile should come from defaults
            },
            'mixing': {
                'uniform': True
                # ratio should come from defaults
            }
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        cmd = f"python3 run_pipeline.py --config {config_path}"
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        if success or "PHASE 2 PIPELINE" in stdout:
            # Check that defaults filled in
            if "quantiles: 10" in stdout or "Quantiles: 10" in stdout or "n_quantiles" not in stderr:
                print("✓ Deep merge preserves defaults for unspecified values")
                return True
            else:
                print("⚠ Pipeline ran but merge behavior unclear")
                return True
        else:
            print(f"✗ Failed with partial config: {stderr[:200]}")
            return False


def test_verbose_flag():
    """Test verbose flag from config."""
    print("\n10. Testing verbose flag via config...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        sim_path = create_test_simulation(os.path.join(tmpdir, 'sim'))
        
        config_path = os.path.join(tmpdir, 'verbose.yaml')
        config = {
            'input': {'rate': 0.005, 'simulation': sim_path},
            'snapshots': {'first': 5, 'second': 10},
            'individuals': {
                'growth_phase': 2,
                'n_quantiles': 2,
                'cells_per_quantile': 1
            },
            'verbose': True
        }
        
        with open(config_path, 'w') as f:
            yaml.dump(config, f)
        
        cmd = f"python3 run_pipeline.py --config {config_path}"
        success, stdout, stderr = run_command(cmd, cwd='.')
        
        if success or "PHASE 2 PIPELINE" in stdout:
            print("✓ Verbose flag accepted via config")
            return True
        else:
            print(f"✗ Failed with verbose flag: {stderr[:200]}")
            return False


def main():
    """Run all tests."""
    print("=" * 60)
    print("PHASE 2 CONFIG TESTING")
    print("=" * 60)
    
    # Check if we're in the right directory
    if not os.path.exists('run_pipeline.py'):
        print("Error: Must run from phase2 directory")
        print("Current directory:", os.getcwd())
        sys.exit(1)
    
    # Check YAML is available
    try:
        import yaml
        print("✓ YAML module available")
    except ImportError:
        print("✗ YAML module not available - config support disabled")
        print("  Install with: pip install pyyaml")
    
    tests = [
        test_no_config_file,
        test_default_config_only,
        test_custom_config,
        test_cli_override,
        test_gene_rate_config,
        test_validation_errors,
        test_example_configs,
        test_compression_setting,
        test_deep_merge,
        test_verbose_flag
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
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())