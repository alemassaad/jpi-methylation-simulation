#!/usr/bin/env python3
"""
Generate growth trajectory plots for each individual in phase2 pipeline.
Uses PetriDishPlotter to create consistent plots.
"""

import os
import sys
import glob
import argparse
import gzip
import json
from typing import List, Dict, Optional

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import PetriDish, PetriDishPlotter
from pipeline_utils import load_petri_dish


def plot_individual(filepath: str, output_dir: str, group_name: str, individual_id: int) -> bool:
    """
    Generate plots for a single individual with history.
    
    Args:
        filepath: Path to individual .json.gz file with history
        output_dir: Directory to save plots
        group_name: Group name (mutant/control1/control2)
        individual_id: Individual identifier
    
    Returns:
        bool: True if successful, False otherwise
    """
    try:
        # Load with history
        petri = load_petri_dish(filepath, include_cell_history=True)
        
        # Check if history exists
        if not hasattr(petri, 'history') or not petri.history:
            print(f"    No history found for {group_name} individual {individual_id:02d}")
            return False
        
        # Create plotter
        plotter = PetriDishPlotter(petri)
        
        # Generate title for plots
        title = f"{group_name.capitalize()} Individual {individual_id:02d}"
        
        # Create plot paths
        jsd_path = os.path.join(output_dir, f"{group_name}_{individual_id:02d}_jsd.png")
        meth_path = os.path.join(output_dir, f"{group_name}_{individual_id:02d}_methylation.png")
        combined_path = os.path.join(output_dir, f"{group_name}_{individual_id:02d}_combined.png")
        
        # Generate plots
        plotter.plot_jsd(title, jsd_path)
        plotter.plot_methylation(title, meth_path)
        plotter.plot_combined(title, combined_path)
        
        print(f"    Generated plots for {group_name} individual {individual_id:02d}")
        return True
        
    except Exception as e:
        print(f"    Error plotting {group_name} individual {individual_id:02d}: {e}")
        return False


def plot_all_individuals(pipeline_dir: str, plot_combined: bool = True) -> None:
    """
    Generate plots for all individuals in a pipeline run.
    
    Args:
        pipeline_dir: Path to pipeline results directory
        plot_combined: Whether to generate combined plots (default: True)
    """
    # Find group directories - check both old and new structures
    if os.path.exists(os.path.join(pipeline_dir, "individuals")):
        # New structure: individuals/mutant, individuals/control1, etc.
        mutant_dir = os.path.join(pipeline_dir, "individuals", "mutant")
        control1_dir = os.path.join(pipeline_dir, "individuals", "control1")
        control2_dir = os.path.join(pipeline_dir, "individuals", "control2")
    else:
        # Old structure: mutant/, control1/, control2/ directly
        mutant_dir = os.path.join(pipeline_dir, "mutant")
        control1_dir = os.path.join(pipeline_dir, "control1")
        control2_dir = os.path.join(pipeline_dir, "control2")
    
    # Create plots directory
    plots_dir = os.path.join(pipeline_dir, "individual_plots")
    os.makedirs(plots_dir, exist_ok=True)
    
    print(f"\nGenerating individual growth trajectory plots...")
    print(f"Output directory: {plots_dir}")
    
    success_count = 0
    total_count = 0
    
    # Process mutant individuals
    if os.path.exists(mutant_dir):
        print("\n  Processing mutant individuals...")
        mutant_files = sorted(glob.glob(os.path.join(mutant_dir, "individual_*.json.gz")))
        for i, filepath in enumerate(mutant_files):
            total_count += 1
            if plot_individual(filepath, plots_dir, "mutant", i):
                success_count += 1
    
    # Process control1 individuals
    if os.path.exists(control1_dir):
        print("\n  Processing control1 individuals...")
        control1_files = sorted(glob.glob(os.path.join(control1_dir, "individual_*.json.gz")))
        for i, filepath in enumerate(control1_files):
            total_count += 1
            if plot_individual(filepath, plots_dir, "control1", i):
                success_count += 1
    
    # Process control2 individuals (usually no history, but check anyway)
    if os.path.exists(control2_dir):
        print("\n  Processing control2 individuals...")
        control2_files = sorted(glob.glob(os.path.join(control2_dir, "individual_*.json.gz")))
        for i, filepath in enumerate(control2_files):
            total_count += 1
            if plot_individual(filepath, plots_dir, "control2", i):
                success_count += 1
    
    # Summary
    print(f"\n  Summary:")
    print(f"    Successfully plotted: {success_count}/{total_count} individuals")
    if success_count < total_count:
        print(f"    Note: {total_count - success_count} individuals had no history data")
    print(f"    Plots saved to: {plots_dir}")


def main():
    parser = argparse.ArgumentParser(
        description="Generate growth trajectory plots for phase2 individuals",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument("pipeline_dir", type=str,
                       help="Path to phase2 pipeline results directory")
    parser.add_argument("--no-combined", action='store_true',
                       help="Skip generation of combined plots")
    
    args = parser.parse_args()
    
    # Validate directory
    if not os.path.exists(args.pipeline_dir):
        print(f"Error: Directory not found: {args.pipeline_dir}")
        return 1
    
    # Check for group directories (support both old and new structures)
    has_groups = False
    if os.path.exists(os.path.join(args.pipeline_dir, "individuals")):
        # New structure
        for group in ["mutant", "control1", "control2"]:
            if os.path.exists(os.path.join(args.pipeline_dir, "individuals", group)):
                has_groups = True
                break
    else:
        # Old structure
        for group in ["mutant", "control1", "control2"]:
            if os.path.exists(os.path.join(args.pipeline_dir, group)):
                has_groups = True
                break
    
    if not has_groups:
        print(f"Error: No group directories found in {args.pipeline_dir}")
        print("Expected structure:")
        print("  - individuals/mutant/, individuals/control1/, individuals/control2/")
        print("  OR")
        print("  - mutant/, control1/, control2/ (legacy)")
        return 1
    
    # Generate plots
    plot_all_individuals(args.pipeline_dir, plot_combined=not args.no_combined)
    
    return 0


if __name__ == "__main__":
    exit(main())