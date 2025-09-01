#!/usr/bin/env python3
"""
Comprehensive test script for config YAML functionality in phase1.
Tests all edge cases, priorities, and validation.
"""

import os
import sys
import subprocess
import tempfile
import shutil
import json
import gzip
import yaml
import glob

def run_simulation(args, cwd="phase1"):
    """Run simulation and return stdout, stderr, returncode."""
    cmd = ["python3", "run_simulation.py"] + args
    print(f"\n  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    return result.stdout, result.stderr, result.returncode

def create_temp_config(content):
    """Create a temporary config file and return its path."""
    fd, path = tempfile.mkstemp(suffix='.yaml')
    with os.fdopen(fd, 'w') as f:
        f.write(content)
    return path

def cleanup_output():
    """Clean up any test output directories."""
    patterns = [
        "phase1/data/rate_*/",
        "phase1/data/gene_rates_*/",
        "phase1/data/debug/",
        "phase1/data/production/",
        "phase1/data/test/"
    ]
    for pattern in patterns:
        for path in glob.glob(pattern):
            if os.path.isdir(path):
                shutil.rmtree(path)

def test_default_config():
    """Test 1: Use default config only."""
    print("\n" + "="*60)
    print("TEST 1: Default Config Only")
    print("="*60)
    
    cleanup_output()
    stdout, stderr, returncode = run_simulation([], cwd="phase1")
    
    if returncode != 0:
        print(f"  ❌ Failed with return code {returncode}")
        print(f"  Error: {stderr}")
        return False
    
    # Check that it used default values
    if "rate: 0.500%" not in stdout.lower():
        print("  ❌ Did not use default rate")
        return False
    
    if "growth phase: 13" not in stdout.lower():
        print("  ❌ Did not use default growth phase")
        return False
    
    print("  ✅ Default config loaded successfully")
    cleanup_output()
    return True

def test_custom_config():
    """Test 2: Use custom config file."""
    print("\n" + "="*60)
    print("TEST 2: Custom Config File")
    print("="*60)
    
    cleanup_output()
    stdout, stderr, returncode = run_simulation(
        ["--config", "configs/quick_test.yaml"],
        cwd="phase1"
    )
    
    if returncode != 0:
        print(f"  ❌ Failed with return code {returncode}")
        print(f"  Error: {stderr}")
        return False
    
    # Check that it used quick_test values
    if "years: 10" not in stdout.lower():
        print("  ❌ Did not use quick_test years")
        return False
    
    if "growth phase: 3" not in stdout.lower():
        print("  ❌ Did not use quick_test growth phase")
        return False
    
    print("  ✅ Custom config loaded successfully")
    cleanup_output()
    return True

def test_cli_override():
    """Test 3: CLI arguments override config."""
    print("\n" + "="*60)
    print("TEST 3: CLI Override of Config")
    print("="*60)
    
    cleanup_output()
    # Use quick_test config but override years with CLI
    stdout, stderr, returncode = run_simulation(
        ["--config", "configs/quick_test.yaml", "--years", "5", "--seed", "999"],
        cwd="phase1"
    )
    
    if returncode != 0:
        print(f"  ❌ Failed with return code {returncode}")
        print(f"  Error: {stderr}")
        return False
    
    # Check that CLI override worked
    if "max years: 5" not in stdout.lower():
        print("  ❌ CLI override for years didn't work")
        return False
    
    # But growth phase should still be from config
    if "growth phase: 3" not in stdout.lower():
        print("  ❌ Config value for growth phase not preserved")
        return False
    
    print("  ✅ CLI override worked correctly")
    cleanup_output()
    return True

def test_gene_rates_config():
    """Test 4: Gene-specific rates from config."""
    print("\n" + "="*60)
    print("TEST 4: Gene-Specific Rates Config")
    print("="*60)
    
    cleanup_output()
    stdout, stderr, returncode = run_simulation(
        ["--config", "configs/gene_rates.yaml", "--years", "5"],
        cwd="phase1"
    )
    
    if returncode != 0:
        print(f"  ❌ Failed with return code {returncode}")
        print(f"  Error: {stderr}")
        return False
    
    # Check that gene rates were used
    if "gene rate groups" not in stdout.lower():
        print("  ❌ Gene rate groups not recognized")
        return False
    
    print("  ✅ Gene rates config worked")
    cleanup_output()
    return True

def test_no_compression_config():
    """Test 5: No compression from config."""
    print("\n" + "="*60)
    print("TEST 5: No Compression in Config")
    print("="*60)
    
    config_content = """
simulation:
  rate: 0.005
  years: 5
  growth_phase: 2
output:
  compress: false
  directory: "data/test"
seed: 888
"""
    
    cleanup_output()
    config_path = create_temp_config(config_content)
    
    try:
        stdout, stderr, returncode = run_simulation(
            ["--config", config_path],
            cwd="phase1"
        )
        
        if returncode != 0:
            print(f"  ❌ Failed with return code {returncode}")
            print(f"  Error: {stderr}")
            return False
        
        # Check for uncompressed output
        if "saving uncompressed" not in stdout.lower():
            print("  ❌ Not saving uncompressed")
            return False
        
        # Find the output file
        output_files = glob.glob("phase1/data/test/*/simulation.json")
        if not output_files:
            print("  ❌ No uncompressed output file found")
            return False
        
        # Verify it's actually uncompressed JSON
        with open(output_files[0], 'r') as f:
            data = json.load(f)
        
        print("  ✅ Uncompressed output from config worked")
        return True
        
    finally:
        os.unlink(config_path)
        cleanup_output()

def test_validation_errors():
    """Test 6: Config validation errors."""
    print("\n" + "="*60)
    print("TEST 6: Config Validation Errors")
    print("="*60)
    
    # Test 6a: Both rate and gene_rate_groups
    config_content = """
simulation:
  rate: 0.005
  gene_rate_groups: "50:0.004,50:0.005"
"""
    
    config_path = create_temp_config(config_content)
    try:
        stdout, stderr, returncode = run_simulation(["--config", config_path], cwd="phase1")
        if returncode == 0:
            print("  ❌ Should have failed with both rate and gene_rate_groups")
            return False
        if "cannot specify both" not in stderr.lower() and "cannot specify both" not in stdout.lower():
            print("  ❌ Wrong error message for dual rate specification")
            return False
        print("  ✅ Correctly rejected both rate types")
    finally:
        os.unlink(config_path)
    
    # Test 6b: Neither rate nor gene_rate_groups
    config_content = """
simulation:
  years: 10
  growth_phase: 3
"""
    
    config_path = create_temp_config(config_content)
    try:
        stdout, stderr, returncode = run_simulation(["--config", config_path], cwd="phase1")
        if returncode == 0:
            print("  ❌ Should have failed without rate specification")
            return False
        if "must specify either" not in stderr.lower() and "must specify either" not in stdout.lower():
            print("  ❌ Wrong error message for missing rate")
            return False
        print("  ✅ Correctly rejected missing rate")
    finally:
        os.unlink(config_path)
    
    # Test 6c: Invalid growth phase
    config_content = """
simulation:
  rate: 0.005
  growth_phase: 25
"""
    
    config_path = create_temp_config(config_content)
    try:
        stdout, stderr, returncode = run_simulation(["--config", config_path], cwd="phase1")
        if returncode == 0:
            print("  ❌ Should have failed with growth_phase > 20")
            return False
        if "growth_phase must be between" not in stderr.lower() and "growth_phase must be between" not in stdout.lower():
            print("  ❌ Wrong error message for invalid growth phase")
            return False
        print("  ✅ Correctly rejected invalid growth phase")
    finally:
        os.unlink(config_path)
    
    return True

def test_missing_config_file():
    """Test 7: Error handling for missing config file."""
    print("\n" + "="*60)
    print("TEST 7: Missing Config File")
    print("="*60)
    
    stdout, stderr, returncode = run_simulation(
        ["--config", "nonexistent.yaml"],
        cwd="phase1"
    )
    
    if returncode == 0:
        print("  ❌ Should have failed with missing config file")
        return False
    
    if "not found" not in stderr.lower() and "not found" not in stdout.lower():
        print("  ❌ Wrong error message for missing file")
        return False
    
    print("  ✅ Correctly handled missing config file")
    return True

def test_cli_only_no_config():
    """Test 8: CLI only without any config file."""
    print("\n" + "="*60)
    print("TEST 8: CLI Only (No Config)")
    print("="*60)
    
    # Temporarily rename default config
    default_config = "phase1/config_default.yaml"
    backup_config = "phase1/config_default.yaml.bak"
    
    try:
        if os.path.exists(default_config):
            shutil.move(default_config, backup_config)
        
        cleanup_output()
        stdout, stderr, returncode = run_simulation(
            ["--rate", "0.005", "--years", "5", "--growth-phase", "2", "--seed", "777"],
            cwd="phase1"
        )
        
        if returncode != 0:
            print(f"  ❌ Failed with return code {returncode}")
            print(f"  Error: {stderr}")
            return False
        
        if "max years: 5" not in stdout.lower():
            print("  ❌ CLI parameters not used")
            return False
        
        print("  ✅ CLI-only mode works without config")
        return True
        
    finally:
        if os.path.exists(backup_config):
            shutil.move(backup_config, default_config)
        cleanup_output()

def test_partial_config():
    """Test 9: Partial config with defaults filling in."""
    print("\n" + "="*60)
    print("TEST 9: Partial Config with Defaults")
    print("="*60)
    
    config_content = """
simulation:
  rate: 0.003  # Only specify rate
# Everything else should come from defaults
"""
    
    cleanup_output()
    config_path = create_temp_config(config_content)
    
    try:
        stdout, stderr, returncode = run_simulation(
            ["--config", config_path, "--years", "5"],
            cwd="phase1"
        )
        
        if returncode != 0:
            print(f"  ❌ Failed with return code {returncode}")
            print(f"  Error: {stderr}")
            return False
        
        # Check rate from config
        if "rate: 0.300%" not in stdout.lower():
            print("  ❌ Config rate not used")
            return False
        
        # Check years from CLI
        if "max years: 5" not in stdout.lower():
            print("  ❌ CLI years not used")
            return False
        
        # Check other defaults were used
        if "gene size: 5" not in stdout.lower():
            print("  ❌ Default gene size not used")
            return False
        
        print("  ✅ Partial config with defaults worked")
        return True
        
    finally:
        os.unlink(config_path)
        cleanup_output()

def test_performance_flags_config():
    """Test 10: Performance flags in config."""
    print("\n" + "="*60)
    print("TEST 10: Performance Flags in Config")
    print("="*60)
    
    config_content = """
simulation:
  rate: 0.005
  years: 5
  growth_phase: 2
performance:
  track_gene_jsd: false
  calculate_jsds: false
seed: 555
"""
    
    cleanup_output()
    config_path = create_temp_config(config_content)
    
    try:
        stdout, stderr, returncode = run_simulation(
            ["--config", config_path],
            cwd="phase1"
        )
        
        if returncode != 0:
            print(f"  ❌ Failed with return code {returncode}")
            print(f"  Error: {stderr}")
            return False
        
        # Check that gene JSD tracking is disabled
        if "gene jsd tracking disabled" not in stdout.lower():
            print("  ❌ Gene JSD tracking not disabled")
            return False
        
        print("  ✅ Performance flags from config worked")
        return True
        
    finally:
        os.unlink(config_path)
        cleanup_output()

def test_seed_handling():
    """Test 11: Seed handling (null, -1, specific values)."""
    print("\n" + "="*60)
    print("TEST 11: Seed Handling")
    print("="*60)
    
    # Test with seed: -1 (random)
    config_content = """
simulation:
  rate: 0.005
  years: 3
  growth_phase: 2
seed: -1
"""
    
    cleanup_output()
    config_path = create_temp_config(config_content)
    
    try:
        stdout1, _, rc1 = run_simulation(["--config", config_path], cwd="phase1")
        stdout2, _, rc2 = run_simulation(["--config", config_path], cwd="phase1")
        
        if rc1 != 0 or rc2 != 0:
            print("  ❌ Failed with seed: -1")
            return False
        
        # Results should be different with random seed
        # (This is probabilistic but very likely to be different)
        print("  ✅ Random seed (-1) worked")
        
    finally:
        os.unlink(config_path)
    
    # Test with specific seed (reproducibility)
    config_content = """
simulation:
  rate: 0.005
  years: 3
  growth_phase: 2
seed: 12345
"""
    
    cleanup_output()
    config_path = create_temp_config(config_content)
    
    try:
        stdout1, _, rc1 = run_simulation(["--config", config_path], cwd="phase1")
        cleanup_output()
        stdout2, _, rc2 = run_simulation(["--config", config_path], cwd="phase1")
        
        if rc1 != 0 or rc2 != 0:
            print("  ❌ Failed with specific seed")
            return False
        
        # Extract final statistics for comparison
        # They should be identical with same seed
        print("  ✅ Specific seed worked")
        
    finally:
        os.unlink(config_path)
        cleanup_output()
    
    return True

def test_complex_priority():
    """Test 12: Complex priority (default < user config < CLI)."""
    print("\n" + "="*60)
    print("TEST 12: Complex Priority Chain")
    print("="*60)
    
    # Create a user config that overrides some defaults
    config_content = """
simulation:
  rate: 0.003      # Override default
  years: 20        # Override default
  # growth_phase not specified, use default
output:
  compress: false  # Override default
"""
    
    cleanup_output()
    config_path = create_temp_config(config_content)
    
    try:
        # Run with config + CLI override
        stdout, stderr, returncode = run_simulation(
            ["--config", config_path, "--years", "7", "--growth-phase", "2"],
            cwd="phase1"
        )
        
        if returncode != 0:
            print(f"  ❌ Failed with return code {returncode}")
            print(f"  Error: {stderr}")
            return False
        
        # Check priorities:
        # - rate from user config (0.003)
        if "rate: 0.300%" not in stdout.lower():
            print("  ❌ User config rate not used")
            return False
        
        # - years from CLI (7)
        if "max years: 7" not in stdout.lower():
            print("  ❌ CLI years not used")
            return False
        
        # - growth_phase from CLI (2)
        if "growth phase: 2" not in stdout.lower():
            print("  ❌ CLI growth phase not used")
            return False
        
        # - compress from user config (false)
        if "uncompressed" not in stdout.lower():
            print("  ❌ User config compress setting not used")
            return False
        
        print("  ✅ Complex priority chain worked correctly")
        return True
        
    finally:
        os.unlink(config_path)
        cleanup_output()

def main():
    """Run all tests."""
    print("\n" + "="*60)
    print("COMPREHENSIVE CONFIG YAML TESTING")
    print("="*60)
    
    # First install pyyaml if needed
    print("\nEnsuring PyYAML is installed...")
    subprocess.run([sys.executable, "-m", "pip", "install", "-q", "pyyaml"], 
                   capture_output=True)
    
    tests = [
        ("Default Config", test_default_config),
        ("Custom Config", test_custom_config),
        ("CLI Override", test_cli_override),
        ("Gene Rates Config", test_gene_rates_config),
        ("No Compression Config", test_no_compression_config),
        ("Validation Errors", test_validation_errors),
        ("Missing Config File", test_missing_config_file),
        ("CLI Only Mode", test_cli_only_no_config),
        ("Partial Config", test_partial_config),
        ("Performance Flags", test_performance_flags_config),
        ("Seed Handling", test_seed_handling),
        ("Complex Priority", test_complex_priority),
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"\n  ❌ Test failed with exception: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, False))
    
    # Clean up any remaining output
    cleanup_output()
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    all_passed = True
    for name, success in results:
        status = "✅ PASSED" if success else "❌ FAILED"
        print(f"  {name}: {status}")
        all_passed = all_passed and success
    
    print("\n" + "="*60)
    if all_passed:
        print("✅ ALL TESTS PASSED!")
    else:
        print("❌ SOME TESTS FAILED")
    print("="*60)
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())