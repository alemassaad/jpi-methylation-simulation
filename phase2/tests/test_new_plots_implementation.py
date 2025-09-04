#!/usr/bin/env python3
"""
Test script for new plot implementations.
Tests each of the 5 new plots as they are implemented.
"""

import os
import sys
import json
import subprocess
import tempfile
import shutil

# Add paths
phase2_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
phase1_dir = os.path.join(os.path.dirname(phase2_dir), 'phase1')
sys.path.insert(0, phase2_dir)
sys.path.append(phase1_dir)

def run_small_simulation():
    """Run a small phase1 simulation for testing."""
    print("\nRunning small test simulation...")
    
    # Create a temporary directory for test data
    test_dir = tempfile.mkdtemp(prefix="test_plots_")
    
    # Run phase1 simulation
    cmd = [
        "python", os.path.join(phase1_dir, "run_simulation.py"),
        "--rate", "0.005",
        "--years", "60",
        "--growth-phase", "4",  # Small for quick test
        "--sites", "100",
        "--gene-size", "5",
        "--seed", "42",
        "--output", test_dir,
        "--no-compress"
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=phase1_dir)
    if result.returncode != 0:
        print(f"  ✗ Simulation failed: {result.stderr}")
        print(f"  Stdout: {result.stdout[:500]}")
        return None
    
    print(f"  ✓ Simulation completed")
    
    # Phase1 creates nested directory structure, search for simulation.json
    sim_files = []
    for root, dirs, files in os.walk(test_dir):
        for file in files:
            if file == "simulation.json":
                sim_files.append(os.path.join(root, file))
    
    if sim_files:
        sim_path = sim_files[0]  # Take the first one found
        print(f"  ✓ Created simulation: {sim_path}")
        return sim_path
    
    print(f"  ✗ Could not find simulation.json in {test_dir}")
    # List what we found
    for root, dirs, files in os.walk(test_dir):
        if files:
            print(f"    Directory: {root}")
            for file in files[:3]:  # Show first 3 files
                print(f"      - {file}")
    
    return None

def run_pipeline_test(simulation_path):
    """Run phase2 pipeline to test new plots."""
    print("\nRunning pipeline with new plots...")
    
    # Create output directory
    output_dir = tempfile.mkdtemp(prefix="pipeline_test_")
    
    cmd = [
        "python", os.path.join(phase2_dir, "run_pipeline.py"),
        "--simulation", simulation_path,
        "--first-snapshot", "30",
        "--second-snapshot", "50",
        "--individual-growth-phase", "3",
        "--n-quantiles", "3",
        "--cells-per-quantile", "2",
        "--seed", "42",
        "--output", output_dir,
        "--no-compress",
        "--verbose"
    ]
    
    print(f"  Running: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=phase2_dir)
    
    if result.returncode != 0:
        print(f"  ✗ Pipeline failed:")
        print(result.stderr)
        return None
        
    print(f"  ✓ Pipeline completed successfully")
    
    # Find results directory
    for root, dirs, files in os.walk(output_dir):
        if "results" in dirs:
            return os.path.join(root, "results")
    
    return None

def test_cell_methylation_histogram(results_dir):
    """Test that cell methylation histograms were created."""
    print("\n1. Testing Cell Methylation Snapshot Histogram...")
    
    expected_files = [
        "year30_methylation_histogram.png",
        "year50_methylation_histogram.png"
    ]
    
    all_found = True
    for file in expected_files:
        path = os.path.join(results_dir, file)
        if os.path.exists(path):
            size = os.path.getsize(path)
            print(f"  ✓ Found {file} ({size:,} bytes)")
        else:
            print(f"  ✗ Missing {file}")
            all_found = False
    
    return all_found

def test_cell_methylation_comparison(results_dir):
    """Test that cell methylation batch comparison was created."""
    print("\n2. Testing Cell Methylation Batch Comparison...")
    
    expected_files = [
        "cell_methylation_comparison.png",
        "cell_methylation_analysis.json"
    ]
    
    all_found = True
    for file in expected_files:
        path = os.path.join(results_dir, file)
        if os.path.exists(path):
            size = os.path.getsize(path)
            print(f"  ✓ Found {file} ({size:,} bytes)")
            
            # Check JSON content if it exists
            if file.endswith(".json"):
                with open(path, 'r') as f:
                    data = json.load(f)
                    if 'summary_statistics' in data:
                        print("    - Contains summary_statistics")
                    if 'statistical_tests' in data:
                        print("    - Contains statistical_tests")
                    if 'individual_means' in data:
                        print("    - Contains individual_means")
        else:
            print(f"  ✗ Missing {file}")
            all_found = False
    
    return all_found

def test_gene_jsd_histogram(results_dir):
    """Test that gene-level JSD histograms were created."""
    print("\n3. Testing Gene-level JSD Snapshot Histogram...")
    
    expected_files = [
        "year30_gene_jsd_histogram.png",
        "year50_gene_jsd_histogram.png"
    ]
    
    all_found = True
    for file in expected_files:
        path = os.path.join(results_dir, file)
        if os.path.exists(path):
            size = os.path.getsize(path)
            print(f"  ✓ Found {file} ({size:,} bytes)")
        else:
            print(f"  ✗ Missing {file}")
            all_found = False
    
    return all_found

def test_gene_jsd_trajectories(results_dir):
    """Test that gene-level JSD trajectories were created."""
    print("\n4. Testing Gene-level JSD Individual Trajectories...")
    
    trajectories_dir = os.path.join(results_dir, "individual_trajectories")
    if not os.path.exists(trajectories_dir):
        print(f"  ✗ Missing individual_trajectories directory")
        return False
    
    # Check for gene JSD trajectory files
    gene_jsd_files = [f for f in os.listdir(trajectories_dir) 
                      if "gene_jsd" in f and f.endswith(".png")]
    
    if gene_jsd_files:
        print(f"  ✓ Found {len(gene_jsd_files)} gene JSD trajectory files:")
        for f in gene_jsd_files[:3]:  # Show first 3
            size = os.path.getsize(os.path.join(trajectories_dir, f))
            print(f"    - {f} ({size:,} bytes)")
        return True
    else:
        print(f"  ✗ No gene JSD trajectory files found")
        return False

def test_gene_jsd_timeline(results_dir):
    """Test that gene-level JSD timeline was created."""
    print("\n5. Testing Gene-level JSD Original Timeline...")
    
    expected_file = "original_simulation_gene_jsd_timeline.png"
    path = os.path.join(results_dir, expected_file)
    
    if os.path.exists(path):
        size = os.path.getsize(path)
        print(f"  ✓ Found {expected_file} ({size:,} bytes)")
        return True
    else:
        print(f"  ✗ Missing {expected_file}")
        return False

def main():
    """Run all tests."""
    print("="*60)
    print("Testing New Plot Implementations")
    print("="*60)
    
    # Run small simulation
    sim_path = run_small_simulation()
    if not sim_path:
        print("\n❌ Failed to create test simulation")
        return 1
    
    # Run pipeline
    results_dir = run_pipeline_test(sim_path)
    if not results_dir:
        print("\n❌ Failed to run pipeline")
        return 1
    
    print(f"\n  Results directory: {results_dir}")
    
    # Test each new plot type
    print("\n" + "="*60)
    print("Testing New Plots")
    print("="*60)
    
    tests = [
        test_cell_methylation_histogram,  # Implemented
        test_cell_methylation_comparison, # Not yet
        test_gene_jsd_histogram,          # Not yet
        test_gene_jsd_trajectories,       # Not yet
        test_gene_jsd_timeline            # Not yet
    ]
    
    implemented = []
    not_implemented = []
    
    for test in tests:
        try:
            if test(results_dir):
                implemented.append(test.__name__)
            else:
                not_implemented.append(test.__name__)
        except Exception as e:
            print(f"  ⚠️ Test {test.__name__} raised exception: {e}")
            not_implemented.append(test.__name__)
    
    # Summary
    print("\n" + "="*60)
    print("Summary")
    print("="*60)
    
    if implemented:
        print(f"\n✅ Implemented ({len(implemented)}):")
        for name in implemented:
            print(f"  - {name.replace('test_', '').replace('_', ' ').title()}")
    
    if not_implemented:
        print(f"\n⏳ Not Yet Implemented ({len(not_implemented)}):")
        for name in not_implemented:
            print(f"  - {name.replace('test_', '').replace('_', ' ').title()}")
    
    # Cleanup temp directories (optional)
    print("\n  Keeping test files for inspection.")
    print(f"  Simulation: {os.path.dirname(sim_path)}")
    print(f"  Results: {os.path.dirname(results_dir)}")
    
    return 0 if len(implemented) > 0 else 1

if __name__ == "__main__":
    exit(main())