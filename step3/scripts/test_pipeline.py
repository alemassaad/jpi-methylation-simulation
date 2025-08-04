#!/usr/bin/env python3
"""
Test the complete step3 pipeline with small test data.
"""

import subprocess
import os
import sys

def run_command(cmd, description):
    """Run a command and report results."""
    print(f"\n{'='*60}")
    print(f"{description}")
    print(f"Command: {cmd}")
    print(f"{'='*60}")
    
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    
    if result.stdout:
        print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr)
    
    if result.returncode != 0:
        print(f"ERROR: Command failed with return code {result.returncode}")
        return False
    
    return True

def main():
    """Run the test pipeline."""
    
    print("TESTING STEP 3 PIPELINE WITH SMALL DATA")
    print("="*80)
    
    # Create test output directories
    os.makedirs("../data/test_individuals/mutant", exist_ok=True)
    os.makedirs("../data/test_individuals/control", exist_ok=True)
    os.makedirs("../test_plots", exist_ok=True)
    
    # Step 1: Create test individuals (both types)
    success = run_command(
        "python create_individuals.py "
        "--type both "
        "--test-mode "
        "--year60-snapshot ../data/test_snapshots/test_year60_snapshot.json.gz "
        "--lineage-base-dir ../../step2 "
        "--output-base-dir ../data/test_individuals",
        "Creating test individuals from test lineages"
    )
    
    if not success:
        print("\nFailed to create individuals!")
        return 1
    
    # Step 2: Plot distributions
    success = run_command(
        "python plot_distributions.py "
        "--mutant-dir ../data/test_individuals/mutant "
        "--control-dir ../data/test_individuals/control "
        "--output-dir ../test_plots",
        "Plotting JSD distributions"
    )
    
    if not success:
        print("\nFailed to plot distributions!")
        return 1
    
    print("\n" + "="*80)
    print("TEST PIPELINE COMPLETED SUCCESSFULLY!")
    print("="*80)
    print("\nCheck the following outputs:")
    print("  - Test individuals: step3/data/test_individuals/")
    print("  - Test plots: step3/test_plots/")
    print("  - Statistics: step3/data/results/jsd_distributions.json")
    
    return 0

if __name__ == "__main__":
    sys.exit(main())