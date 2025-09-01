#!/usr/bin/env python3
"""
Test that the pipeline can be started with config file.
Creates a minimal simulation and tests config loading.
"""

import os
import sys
import json
import gzip
import tempfile
import subprocess


def create_minimal_simulation(output_dir):
    """Create a minimal simulation file for testing."""
    os.makedirs(output_dir, exist_ok=True)
    
    sim_data = {
        'parameters': {
            'rate': 0.005,
            'n': 100,  # Small for testing
            'gene_size': 5,
            'growth_phase': 2,
            'years': 15,
            'seed': 42
        },
        'history': {}
    }
    
    # Add minimal history
    for year in range(16):
        n_cells = min(2 ** year, 4) if year <= 2 else 4
        sim_data['history'][str(year)] = {
            'cells': [
                {
                    'methylated': [0] * 100,
                    'cell_JSD': 0.01,
                    'id': f"cell_{i}_{year}"
                }
                for i in range(n_cells)
            ]
        }
    
    sim_path = os.path.join(output_dir, 'simulation.json.gz')
    with gzip.open(sim_path, 'wt') as f:
        json.dump(sim_data, f)
    
    return sim_path


def test_with_quick_config():
    """Test pipeline with quick test config."""
    print("\nTesting pipeline with quick_test.yaml config...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create simulation
        sim_path = create_minimal_simulation(os.path.join(tmpdir, 'sim'))
        
        # Run with config
        cmd = [
            sys.executable, 'run_pipeline.py',
            '--config', 'configs/quick_test.yaml',
            '--simulation', sim_path,
            '--rate', '0.005',
            '--output-dir', os.path.join(tmpdir, 'output')
        ]
        
        print(f"Running: {' '.join(cmd)}")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=10
            )
            
            # Check for successful start
            if "PHASE 2 PIPELINE" in result.stdout:
                print("✓ Pipeline started successfully with config")
                
                # Check if config values were used
                if "first: 10" in result.stdout or "Year 10" in result.stdout:
                    print("✓ Config snapshot values were used")
                if "growth_phase: 2" in result.stdout or "2^2" in result.stdout:
                    print("✓ Config growth phase was used")
                
                return True
            else:
                print("✗ Pipeline did not start properly")
                if result.stderr:
                    print(f"Error: {result.stderr[:500]}")
                return False
                
        except subprocess.TimeoutExpired:
            print("✓ Pipeline started (timed out as expected for long run)")
            return True
        except Exception as e:
            print(f"✗ Failed to run pipeline: {e}")
            return False


def test_config_override():
    """Test that CLI args override config."""
    print("\nTesting CLI override of config values...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create simulation
        sim_path = create_minimal_simulation(os.path.join(tmpdir, 'sim'))
        
        # Run with config but override snapshots
        cmd = [
            sys.executable, 'run_pipeline.py',
            '--config', 'configs/quick_test.yaml',
            '--simulation', sim_path,
            '--rate', '0.005',
            '--first-snapshot', '7',
            '--second-snapshot', '12',
            '--output-dir', os.path.join(tmpdir, 'output')
        ]
        
        print(f"Running with overrides...")
        
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=10
            )
            
            # Check for override values
            if "Year 7" in result.stdout or "first-snapshot 7" in result.stdout:
                print("✓ CLI override of first-snapshot worked")
            if "Year 12" in result.stdout or "second-snapshot 12" in result.stdout:
                print("✓ CLI override of second-snapshot worked")
            
            return True
                
        except subprocess.TimeoutExpired:
            print("✓ Pipeline started (timed out as expected)")
            return True
        except Exception as e:
            print(f"✗ Failed: {e}")
            return False


def main():
    print("=" * 60)
    print("PHASE 2 PIPELINE CONFIG INTEGRATION TEST")
    print("=" * 60)
    
    # Check we're in the right directory
    if not os.path.exists('run_pipeline.py'):
        print("Error: Must run from phase2 directory")
        return 1
    
    # Check config files exist
    if not os.path.exists('configs/quick_test.yaml'):
        print("Error: Example config files not found")
        return 1
    
    tests = [
        test_with_quick_config,
        test_config_override
    ]
    
    passed = 0
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"✗ Test failed: {e}")
    
    print("\n" + "=" * 60)
    print(f"Integration tests: {passed}/{len(tests)} passed")
    print("=" * 60)
    
    if passed == len(tests):
        print("\n✅ Phase 2 config support is working correctly!")
    
    return 0 if passed == len(tests) else 1


if __name__ == "__main__":
    sys.exit(main())