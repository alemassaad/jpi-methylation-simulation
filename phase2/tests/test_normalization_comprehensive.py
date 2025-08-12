#!/usr/bin/env python3
"""
Comprehensive test for normalization feature covering all edge cases.
Tests all combinations of:
- Normalization alone
- Normalization + uniform mixing
- Edge cases (all excluded, single individual, etc.)
"""

import os
import sys
import subprocess
import tempfile
import json
import gzip
import glob
import shutil

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

def run_phase1_simulation(years=20, growth_phase=4, rate=0.005, seed=42):
    """Run a quick phase1 simulation for testing."""
    print(f"\n{'='*60}")
    print("Running Phase 1 simulation...")
    print(f"{'='*60}")
    
    cmd = [
        "python", "../../phase1/run_simulation.py",
        "--rate", str(rate),
        "--years", str(years),
        "--growth-phase", str(growth_phase),
        "--seed", str(seed),
        "-n", "100"  # Small for speed
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print("ERROR running phase1:")
        print(result.stderr)
        return None
    
    # Find the simulation file
    pattern = f"../../phase1/data/rate_{rate:.5f}/grow{growth_phase}-sites100-years{years}-seed{seed}-*/simulation.json.gz"
    files = glob.glob(pattern)
    if not files:
        print(f"ERROR: No simulation file found matching {pattern}")
        return None
    
    return files[0]

def run_phase2_pipeline(simulation_file, test_name, **kwargs):
    """Run phase2 pipeline with given parameters."""
    print(f"\n{'='*60}")
    print(f"Test: {test_name}")
    print(f"{'='*60}")
    
    # Default parameters
    params = {
        "rate": 0.005,
        "first_snapshot": 10,
        "second_snapshot": 15,
        "individual_growth_phase": 3,  # Target ~8 cells
        "n_quantiles": 4,
        "cells_per_quantile": 2,  # 8 individuals total
        "mix_ratio": 80,
        "seed": 42
    }
    
    # Update with test-specific parameters
    params.update(kwargs)
    
    # Build command
    cmd = [
        "python", "../run_pipeline.py",
        "--rate", str(params["rate"]),
        "--simulation", simulation_file,
        "--first-snapshot", str(params["first_snapshot"]),
        "--second-snapshot", str(params["second_snapshot"]),
        "--individual-growth-phase", str(params["individual_growth_phase"]),
        "--n-quantiles", str(params["n_quantiles"]),
        "--cells-per-quantile", str(params["cells_per_quantile"]),
        "--mix-ratio", str(params["mix_ratio"]),
        "--seed", str(params["seed"])
    ]
    
    # Add optional flags
    if params.get("normalize_size"):
        cmd.append("--normalize-size")
    if params.get("uniform_mixing"):
        cmd.append("--uniform-mixing")
    if params.get("force_recreate"):
        cmd.append("--force-recreate")
    
    print(f"Running: {' '.join(cmd)}")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    success = result.returncode == 0
    
    # Extract key information from output
    output_lines = result.stdout.split('\n')
    
    # Look for normalization statistics
    norm_stats = {}
    for i, line in enumerate(output_lines):
        if "NORMALIZATION STATISTICS" in line:
            # Extract stats from following lines
            for j in range(i+1, min(i+10, len(output_lines))):
                if "Threshold (median - 0.5σ):" in output_lines[j]:
                    parts = output_lines[j].split(":")
                    if len(parts) > 1:
                        norm_stats['threshold'] = parts[1].strip().split()[0]
                elif "Overall retention:" in output_lines[j]:
                    parts = output_lines[j].split(":")
                    if len(parts) > 1:
                        norm_stats['retention'] = parts[1].strip()
    
    # Find output directory
    output_dir = None
    for line in output_lines:
        if "Output directory:" in line:
            output_dir = line.split("Output directory:")[1].strip()
            break
    
    # Check if files were created
    files_created = {}
    if output_dir and os.path.exists(output_dir):
        # Count individuals
        for group in ['mutant', 'control1', 'control2']:
            group_dir = os.path.join(output_dir, 'individuals', group)
            if os.path.exists(group_dir):
                files_created[group] = len(glob.glob(os.path.join(group_dir, '*.json.gz')))
            else:
                files_created[group] = 0
    
    return {
        'success': success,
        'test_name': test_name,
        'norm_stats': norm_stats,
        'files_created': files_created,
        'output_dir': output_dir,
        'stdout': result.stdout,
        'stderr': result.stderr
    }

def analyze_results(results):
    """Analyze and print test results."""
    print(f"\n{'='*80}")
    print("TEST RESULTS SUMMARY")
    print(f"{'='*80}")
    
    all_passed = True
    
    for r in results:
        status = "✓ PASSED" if r['success'] else "✗ FAILED"
        print(f"\n{status}: {r['test_name']}")
        
        if r['norm_stats']:
            print(f"  Normalization threshold: {r['norm_stats'].get('threshold', 'N/A')}")
            print(f"  Retention: {r['norm_stats'].get('retention', 'N/A')}")
        
        if r['files_created']:
            print(f"  Files created:")
            for group, count in r['files_created'].items():
                print(f"    {group}: {count} individuals")
        
        if not r['success']:
            all_passed = False
            print(f"  Error output:")
            if r['stderr']:
                print(f"    {r['stderr'][:500]}")
            # Look for specific error messages
            if "ERROR: All individuals were excluded" in r['stdout']:
                print("    → All individuals excluded by normalization (expected for edge case)")
            elif "WARNING:" in r['stdout']:
                for line in r['stdout'].split('\n'):
                    if "WARNING:" in line:
                        print(f"    → {line.strip()}")
    
    return all_passed

def main():
    """Run comprehensive normalization tests."""
    
    # Use existing simulation or create new one
    existing_sim = "../../phase1/data/rate_0.00500/grow5-sites100-years20-seed47-93ce/simulation.json.gz"
    if os.path.exists(existing_sim):
        print(f"Using existing simulation: {existing_sim}")
        simulation_file = existing_sim
    else:
        # Create a test simulation
        simulation_file = run_phase1_simulation(years=20, growth_phase=4, seed=42)
        if not simulation_file:
            print("ERROR: Could not create test simulation")
            sys.exit(1)
    
    print(f"\nUsing simulation: {simulation_file}")
    
    results = []
    
    # Test 1: Basic pipeline without normalization (baseline)
    results.append(run_phase2_pipeline(
        simulation_file, 
        "1. Baseline (no normalization, no uniform)",
        force_recreate=True
    ))
    
    # Test 2: Normalization only
    results.append(run_phase2_pipeline(
        simulation_file,
        "2. Normalization only",
        normalize_size=True,
        force_recreate=True
    ))
    
    # Test 3: Uniform mixing only
    results.append(run_phase2_pipeline(
        simulation_file,
        "3. Uniform mixing only", 
        uniform_mixing=True,
        force_recreate=True
    ))
    
    # Test 4: Both normalization and uniform mixing
    results.append(run_phase2_pipeline(
        simulation_file,
        "4. Normalization + Uniform mixing",
        normalize_size=True,
        uniform_mixing=True,
        force_recreate=True
    ))
    
    # Test 5: Edge case - very small growth (likely all excluded)
    results.append(run_phase2_pipeline(
        simulation_file,
        "5. Edge case: tiny growth phase",
        individual_growth_phase=1,  # Only 2 cells expected
        normalize_size=True,
        force_recreate=True
    ))
    
    # Test 6: Edge case - single individual per group
    results.append(run_phase2_pipeline(
        simulation_file,
        "6. Edge case: single individual",
        n_quantiles=1,
        cells_per_quantile=1,
        normalize_size=True,
        force_recreate=True
    ))
    
    # Test 7: Different mix ratios with normalization
    results.append(run_phase2_pipeline(
        simulation_file,
        "7. Mix ratio 50% with normalization",
        mix_ratio=50,
        normalize_size=True,
        uniform_mixing=True,
        force_recreate=True
    ))
    
    # Test 8: Edge case - 100% mixing (all snapshot cells)
    results.append(run_phase2_pipeline(
        simulation_file,
        "8. Edge case: 100% mix ratio",
        mix_ratio=100,
        normalize_size=True,
        uniform_mixing=True,
        force_recreate=True
    ))
    
    # Test 9: Reproducibility - run same config twice
    print(f"\n{'='*60}")
    print("Testing reproducibility...")
    print(f"{'='*60}")
    
    r1 = run_phase2_pipeline(
        simulation_file,
        "9a. Reproducibility run 1",
        normalize_size=True,
        uniform_mixing=True,
        seed=123,
        force_recreate=True
    )
    
    r2 = run_phase2_pipeline(
        simulation_file,
        "9b. Reproducibility run 2",
        normalize_size=True,
        uniform_mixing=True,
        seed=123,
        force_recreate=True
    )
    
    # Check if outputs are identical
    reproducible = False
    if r1['output_dir'] and r2['output_dir']:
        # Load and compare statistics
        stats1_file = os.path.join(r1['output_dir'], 'results', 'statistics.json')
        stats2_file = os.path.join(r2['output_dir'], 'results', 'statistics.json')
        
        if os.path.exists(stats1_file) and os.path.exists(stats2_file):
            with open(stats1_file) as f:
                stats1 = json.load(f)
            with open(stats2_file) as f:
                stats2 = json.load(f)
            
            # Compare key statistics
            reproducible = (
                stats1.get('mutant', {}).get('mean') == stats2.get('mutant', {}).get('mean') and
                stats1.get('control1', {}).get('mean') == stats2.get('control1', {}).get('mean')
            )
    
    r1['test_name'] = f"9. Reproducibility: {'✓ IDENTICAL' if reproducible else '✗ DIFFERENT'}"
    results.append(r1)
    
    # Analyze all results
    all_passed = analyze_results(results)
    
    # Clean up test data if requested
    if '--cleanup' in sys.argv:
        print(f"\n{'='*60}")
        print("Cleaning up test data...")
        for r in results:
            if r['output_dir'] and os.path.exists(r['output_dir']):
                shutil.rmtree(r['output_dir'])
                print(f"  Removed: {r['output_dir']}")
    
    # Final summary
    print(f"\n{'='*80}")
    if all_passed:
        print("✓ ALL TESTS PASSED")
    else:
        print("✗ SOME TESTS FAILED")
    print(f"{'='*80}")
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())