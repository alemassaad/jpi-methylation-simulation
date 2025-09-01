#!/usr/bin/env python3
"""
Analyze the distribution of individual sizes after homeostasis.
This helps us decide on the best normalization strategy.
"""

import os
import sys
import glob
import json
import gzip
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

# Add parent directory to path
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish

def load_petri_dish(filepath):
    """Load a PetriDish from JSON file (simplified - just count cells)."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    # We only need the number of cells for size analysis
    num_cells = len(data.get('cells', []))
    
    # Create a simple object to hold the data
    class SimplePetri:
        def __init__(self, num_cells):
            self.cells = [None] * num_cells  # Just need the length
    
    return SimplePetri(num_cells)

def analyze_size_distribution(data_dir):
    """Analyze individual sizes in a phase2 output directory."""
    
    # Find all individual directories
    pattern = os.path.join(data_dir, "**", "individuals", "*", "individual_*.json.gz")
    files = glob.glob(pattern, recursive=True)
    
    if not files:
        print(f"No individual files found in {data_dir}")
        print(f"Searched pattern: {pattern}")
        return None
    
    print(f"Found {len(files)} individual files")
    
    # Group by type
    mutant_sizes = []
    control1_sizes = []
    control2_sizes = []
    
    for filepath in files:
        if 'mutant' in filepath:
            petri = load_petri_dish(filepath)
            mutant_sizes.append(len(petri.cells))
        elif 'control1' in filepath:
            petri = load_petri_dish(filepath)
            control1_sizes.append(len(petri.cells))
        elif 'control2' in filepath:
            petri = load_petri_dish(filepath)
            control2_sizes.append(len(petri.cells))
    
    return {
        'mutant': mutant_sizes,
        'control1': control1_sizes,
        'control2': control2_sizes
    }

def plot_distribution(sizes_dict, title="Individual Size Distribution"):
    """Create visualization of size distributions."""
    
    # Combine mutant and control1 for analysis (exclude control2 as it's pure snapshot)
    if sizes_dict['mutant'] and sizes_dict['control1']:
        combined = sizes_dict['mutant'] + sizes_dict['control1']
    elif sizes_dict['mutant']:
        combined = sizes_dict['mutant']
    elif sizes_dict['control1']:
        combined = sizes_dict['control1']
    else:
        print("No data to plot")
        return
    
    # Calculate statistics
    combined_array = np.array(combined)
    mean = np.mean(combined_array)
    median = np.median(combined_array)
    std = np.std(combined_array)
    
    # Calculate various thresholds
    threshold_20p = np.percentile(combined_array, 20)
    threshold_med_1std = median - std
    threshold_med_0_5std = median - 0.5 * std
    
    # Count how many would be kept with each threshold
    keep_20p = np.sum(combined_array >= threshold_20p)
    keep_med_1std = np.sum(combined_array >= threshold_med_1std)
    keep_med_0_5std = np.sum(combined_array >= threshold_med_0_5std)
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. Histogram with thresholds
    ax1 = axes[0, 0]
    ax1.hist(combined, bins=30, alpha=0.7, color='blue', edgecolor='black')
    ax1.axvline(threshold_20p, color='red', linestyle='--', alpha=0.5, label=f'20th percentile (comparison): {threshold_20p:.0f}')
    ax1.axvline(threshold_med_1std, color='green', linestyle='--', alpha=0.5, label=f'Median-1σ (comparison): {threshold_med_1std:.0f}')
    ax1.axvline(threshold_med_0_5std, color='orange', linestyle='-', linewidth=2, label=f'Median-0.5σ (USED): {threshold_med_0_5std:.0f}')
    ax1.axvline(median, color='black', linestyle='-', alpha=0.5, label=f'Median: {median:.0f}')
    ax1.set_xlabel('Individual Size (cells)')
    ax1.set_ylabel('Count')
    ax1.set_title('Size Distribution with Threshold Options')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # 2. Box plot by group
    ax2 = axes[0, 1]
    data_to_plot = []
    labels = []
    if sizes_dict['mutant']:
        data_to_plot.append(sizes_dict['mutant'])
        labels.append(f"Mutant (n={len(sizes_dict['mutant'])})")
    if sizes_dict['control1']:
        data_to_plot.append(sizes_dict['control1'])
        labels.append(f"Control1 (n={len(sizes_dict['control1'])})")
    
    bp = ax2.boxplot(data_to_plot, labels=labels, patch_artist=True)
    for patch in bp['boxes']:
        patch.set_facecolor('lightblue')
    ax2.set_ylabel('Individual Size (cells)')
    ax2.set_title('Size Distribution by Group')
    ax2.grid(True, alpha=0.3)
    
    # 3. Cumulative distribution
    ax3 = axes[1, 0]
    sorted_sizes = np.sort(combined_array)
    cumulative = np.arange(1, len(sorted_sizes) + 1) / len(sorted_sizes) * 100
    ax3.plot(sorted_sizes, cumulative, 'b-', linewidth=2)
    ax3.axhline(20, color='red', linestyle='--', alpha=0.5, label='20% excluded')
    ax3.axhline(80, color='red', linestyle='--', alpha=0.5, label='80% kept')
    ax3.axvline(threshold_20p, color='red', linestyle='--', alpha=0.5)
    ax3.axvline(threshold_med_1std, color='green', linestyle='--', alpha=0.5)
    ax3.set_xlabel('Individual Size (cells)')
    ax3.set_ylabel('Cumulative Percentage')
    ax3.set_title('Cumulative Distribution')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    
    # 4. Statistics summary
    ax4 = axes[1, 1]
    ax4.axis('off')
    
    stats_text = f"""
    POPULATION STATISTICS
    =====================
    Total individuals: {len(combined)}
    Mean: {mean:.1f} cells
    Median: {median:.1f} cells
    Std Dev: {std:.1f} cells
    Min: {np.min(combined_array):.0f} cells
    Max: {np.max(combined_array):.0f} cells
    CV: {(std/mean)*100:.1f}%
    
    THRESHOLD COMPARISON
    ====================
    Threshold options (for comparison):
    20th percentile:
      Threshold: {threshold_20p:.0f} cells
      Keep: {keep_20p}/{len(combined)} ({keep_20p/len(combined)*100:.1f}%)
    
    Median - 1 StDev:
      Threshold: {threshold_med_1std:.0f} cells
      Keep: {keep_med_1std}/{len(combined)} ({keep_med_1std/len(combined)*100:.1f}%)
    
    Median - 0.5 StDev (ACTUALLY USED):
      Threshold: {threshold_med_0_5std:.0f} cells
      Keep: {keep_med_0_5std}/{len(combined)} ({keep_med_0_5std/len(combined)*100:.1f}%)
    """
    
    ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes, fontsize=10,
             verticalalignment='top', fontfamily='monospace')
    
    plt.suptitle(title, fontsize=14, fontweight='bold')
    plt.tight_layout()
    
    # Save instead of showing
    output_file = 'size_distribution_analysis.png'
    plt.savefig(output_file, dpi=150, bbox_inches='tight')
    print(f"\nPlot saved to: {output_file}")
    plt.close()
    
    # Print summary to console
    print("\n" + "="*60)
    print("ANALYSIS SUMMARY")
    print("="*60)
    print(stats_text)
    
    return {
        'mean': mean,
        'median': median,
        'std': std,
        'threshold_20p': threshold_20p,
        'threshold_med_1std': threshold_med_1std,
        'threshold_med_0_5std': threshold_med_0_5std,
        'retention_20p': keep_20p/len(combined)*100,
        'retention_med_1std': keep_med_1std/len(combined)*100,
        'retention_med_0_5std': keep_med_0_5std/len(combined)*100
    }

def main():
    """Main function to analyze size distributions."""
    
    # Look for most recent phase2 output
    base_dir = "data"
    
    # Find directories with the new parameter structure
    pattern = os.path.join(base_dir, "**", "snap*to*-growth*-quant*-mix*")
    dirs = glob.glob(pattern, recursive=True)
    
    if not dirs:
        print("No phase2 output directories found")
        print(f"Searched pattern: {pattern}")
        return
    
    # Sort by modification time and use most recent
    dirs.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    latest_dir = dirs[0]
    
    print(f"Analyzing: {latest_dir}")
    print("="*60)
    
    # Analyze sizes
    sizes = analyze_size_distribution(latest_dir)
    
    if sizes:
        # Create visualization
        plot_distribution(sizes, title=f"Individual Sizes: {os.path.basename(latest_dir)}")
    
    # Also offer to analyze a specific directory
    print("\nTo analyze a specific directory, run:")
    print("  python analyze_individual_sizes.py /path/to/phase2/output")

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Analyze specific directory
        data_dir = sys.argv[1]
        print(f"Analyzing: {data_dir}")
        sizes = analyze_size_distribution(data_dir)
        if sizes:
            plot_distribution(sizes, title=f"Individual Sizes: {os.path.basename(data_dir)}")
    else:
        main()