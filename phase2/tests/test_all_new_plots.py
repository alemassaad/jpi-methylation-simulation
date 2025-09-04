#!/usr/bin/env python3
"""
Comprehensive test for all 5 new plot implementations.
"""

import os
import sys
import subprocess

# Add paths
phase2_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
phase1_dir = os.path.join(os.path.dirname(phase2_dir), 'phase1')
sys.path.insert(0, phase2_dir)
sys.path.append(phase1_dir)

def run_full_pipeline_test():
    """Run a complete pipeline test with all new plots."""
    print("\n" + "="*60)
    print("Testing All New Plot Implementations")
    print("="*60)
    
    # Use existing simulation
    sim_path = "../phase1/data/gene_rates_20x0.00500/size16-sites100-genesize5-years60-seed123-20250904-180530/simulation.json"
    
    # Check if simulation exists
    full_sim_path = os.path.join(phase2_dir, sim_path)
    if not os.path.exists(full_sim_path):
        print(f"\n❌ Test simulation not found at {full_sim_path}")
        print("   Please run: python phase1/run_simulation.py --rate 0.005 --years 60 --growth-phase 4 --sites 100 --seed 123 --no-compress")
        return 1
    
    print(f"\n✓ Using existing simulation: {sim_path}")
    
    # Run pipeline with all features enabled
    print("\n" + "-"*60)
    print("Running pipeline with all new plots enabled...")
    print("-"*60)
    
    cmd = [
        "python", "run_pipeline.py",
        "--simulation", sim_path,
        "--first-snapshot", "30",
        "--second-snapshot", "50",
        "--individual-growth-phase", "3",
        "--n-quantiles", "3",
        "--cells-per-quantile", "2",
        "--seed", "888",
        "--no-compress",
        "--plot-individuals"  # Enable trajectory plots
    ]
    
    print(f"Command: {' '.join(cmd)}\n")
    
    # Run pipeline
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=phase2_dir)
    
    if result.returncode != 0:
        print(f"❌ Pipeline failed:")
        print(result.stderr[-2000:])  # Show last 2000 chars of error
        return 1
    
    # Extract output directory from stdout
    output_dir = None
    for line in result.stdout.split('\n'):
        if 'Output:' in line and 'data/' in line:
            output_dir = line.split('Output:')[1].strip()
            break
    
    if not output_dir:
        print("❌ Could not find output directory in pipeline output")
        return 1
    
    results_dir = os.path.join(phase2_dir, output_dir, "results")
    print(f"\n✓ Pipeline completed. Results in: {results_dir}")
    
    # Check for all expected files
    print("\n" + "-"*60)
    print("Checking for new plot files...")
    print("-"*60)
    
    expected_files = {
        # 1. Cell Methylation Snapshot Histograms
        "year30_methylation_histogram.png": "Cell Methylation Snapshot (Year 30)",
        "year50_methylation_histogram.png": "Cell Methylation Snapshot (Year 50)",
        
        # 2. Cell Methylation Batch Comparison
        "cell_methylation_comparison.png": "Cell Methylation Batch Comparison",
        "cell_methylation_analysis.json": "Cell Methylation Analysis Data",
        
        # 3. Gene JSD Snapshot Histograms
        "year30_gene_jsd_histogram.png": "Gene JSD Snapshot (Year 30)",
        "year50_gene_jsd_histogram.png": "Gene JSD Snapshot (Year 50)",
        
        # 5. Gene JSD Original Timeline
        "original_simulation_gene_jsd_timeline.png": "Gene JSD Original Timeline"
    }
    
    all_found = True
    found_count = 0
    
    for filename, description in expected_files.items():
        filepath = os.path.join(results_dir, filename)
        if os.path.exists(filepath):
            size = os.path.getsize(filepath)
            print(f"  ✓ {description:45} {size:>10,} bytes")
            found_count += 1
        else:
            print(f"  ✗ {description:45} MISSING")
            all_found = False
    
    # 4. Check for Gene JSD Individual Trajectories
    print("\nChecking individual trajectory plots...")
    traj_dir = os.path.join(results_dir, "individual_trajectories")
    if os.path.exists(traj_dir):
        gene_jsd_files = [f for f in os.listdir(traj_dir) if "gene_jsd" in f]
        if gene_jsd_files:
            print(f"  ✓ Gene JSD Trajectories:                      {len(gene_jsd_files)} files")
            for f in gene_jsd_files[:3]:  # Show first 3
                size = os.path.getsize(os.path.join(traj_dir, f))
                print(f"    - {f:40} {size:>10,} bytes")
            found_count += len(gene_jsd_files)
        else:
            print(f"  ✗ Gene JSD Trajectories:                      MISSING")
            all_found = False
    else:
        print(f"  ✗ Individual trajectories directory:          MISSING")
        all_found = False
    
    # Summary
    print("\n" + "="*60)
    print("Test Summary")
    print("="*60)
    
    if all_found:
        print(f"\n✅ SUCCESS: All new plots generated successfully!")
        print(f"   Total new plot files: {found_count}")
        print(f"\n   Results directory: {results_dir}")
        return 0
    else:
        print(f"\n❌ FAILURE: Some plots are missing")
        print(f"   Found {found_count} plot files")
        print(f"   Check {results_dir} for details")
        return 1

if __name__ == "__main__":
    exit(run_full_pipeline_test())