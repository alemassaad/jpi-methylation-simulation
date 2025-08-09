#!/usr/bin/env python3
"""
Complete the analysis phase for already-processed individuals.
"""

import os
import sys
import json
import numpy as np

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'step1-prime'))

from pipeline_utils import load_all_petri_dishes
from pipeline_analysis import analyze_populations_from_dishes, plot_cell_level_distributions

def main():
    base_dir = "data/rate_0.005000"
    mutant_dir = os.path.join(base_dir, "individuals/mutant")
    control1_dir = os.path.join(base_dir, "individuals/control1")
    control2_dir = os.path.join(base_dir, "individuals/control2")
    results_dir = os.path.join(base_dir, "results")
    plots_dir = os.path.join(base_dir, "plots")
    
    print("Loading populations...")
    mutant_dishes = load_all_petri_dishes(mutant_dir)
    control1_dishes = load_all_petri_dishes(control1_dir)
    
    # Check if control2 exists, if not skip
    if os.path.exists(control2_dir) and len(os.listdir(control2_dir)) > 0:
        control2_dishes = load_all_petri_dishes(control2_dir)
    else:
        print("Control2 not found or empty, skipping analysis")
        return
    
    print(f"  Loaded {len(mutant_dishes)} mutant individuals")
    print(f"  Loaded {len(control1_dishes)} control1 individuals")
    print(f"  Loaded {len(control2_dishes)} control2 individuals")
    
    # Run analysis
    print("\nRunning analysis...")
    analysis_results = analyze_populations_from_dishes(
        mutant_dishes, control1_dishes, control2_dishes, results_dir
    )
    
    # Create cell-level plot
    cell_plot_path = os.path.join(plots_dir, "cell_level_jsd_distributions.png")
    plot_cell_level_distributions(mutant_dishes, control1_dishes, control2_dishes, cell_plot_path)
    
    # Print summary
    stats = analysis_results['statistics']
    print("\nResults:")
    print(f"  Mutant mean JSD: {stats['mutant']['mean']:.6f} ± {stats['mutant']['std']:.6f}")
    print(f"  Control1 mean JSD: {stats['control1']['mean']:.6f} ± {stats['control1']['std']:.6f}")
    print(f"  Control2 mean JSD: {stats['control2']['mean']:.6f} ± {stats['control2']['std']:.6f}")
    
    print("\nStatistical tests:")
    for comparison, values in stats['comparisons'].items():
        print(f"  {comparison}: p={values['p_value']:.6f}")

if __name__ == "__main__":
    main()