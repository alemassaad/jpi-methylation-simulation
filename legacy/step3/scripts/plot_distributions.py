#!/usr/bin/env python3
"""
Load mutant and control individuals and plot their JSD distributions.
Compares individuals created from decile-based lineages vs uniform lineages.
"""

import json
import gzip
import numpy as np
import glob
import os
import plotly.graph_objects as go
import argparse


def load_individuals(pattern):
    """Load all individual files matching pattern and extract mean JSD values."""
    files = sorted(glob.glob(pattern))
    print(f"Found {len(files)} files matching {pattern}")
    
    mean_jsds = []
    for filepath in files:
        with gzip.open(filepath, 'rt') as f:
            data = json.load(f)
        mean_jsd = data['metadata']['mean_jsd']
        mean_jsds.append(mean_jsd)
    
    return np.array(mean_jsds)


def create_comparison_plot(mutant_jsds, control_jsds, output_path="jsd_distribution_comparison.png"):
    """Create a plot comparing mutant and control distributions."""
    
    # Create figure
    fig = go.Figure()
    
    # Add mutant points
    x_mutant = np.ones(len(mutant_jsds)) * 0 + np.random.normal(0, 0.02, len(mutant_jsds))
    fig.add_trace(go.Scatter(
        x=x_mutant,
        y=mutant_jsds,
        mode='markers',
        name='Mutant (Decile-based)',
        marker=dict(color='#1f77b4', size=10, opacity=0.6)
    ))
    
    # Add control points
    x_control = np.ones(len(control_jsds)) * 1 + np.random.normal(0, 0.02, len(control_jsds))
    fig.add_trace(go.Scatter(
        x=x_control,
        y=control_jsds,
        mode='markers',
        name='Control (Uniform)',
        marker=dict(color='#ff7f0e', size=10, opacity=0.6)
    ))
    
    # Add mean lines
    fig.add_hline(y=np.mean(mutant_jsds), line_dash="dash", line_color="#1f77b4", 
                  annotation_text=f"Mutant mean: {np.mean(mutant_jsds):.4f}")
    fig.add_hline(y=np.mean(control_jsds), line_dash="dash", line_color="#ff7f0e",
                  annotation_text=f"Control mean: {np.mean(control_jsds):.4f}")
    
    # Update layout
    fig.update_layout(
        title="JSD Distribution: Mutant vs Control Individuals",
        xaxis=dict(
            tickmode='array',
            tickvals=[0, 1],
            ticktext=['Mutant<br>(Decile-based)', 'Control<br>(Uniform)'],
            range=[-0.5, 1.5]
        ),
        yaxis=dict(title="Mean JSD per Individual"),
        showlegend=True,
        width=800,
        height=600,
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Update axes
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor='lightgray', zeroline=True, zerolinecolor='lightgray')
    
    # Save as PNG
    fig.write_image(output_path)
    print(f"Saved plot to {output_path}")
    
    # Also save as interactive HTML
    html_path = output_path.replace('.png', '.html')
    fig.write_html(html_path)
    print(f"Saved interactive plot to {html_path}")
    
    return fig


def save_statistics(mutant_jsds, control_jsds, output_path="jsd_distributions.json"):
    """Save statistical summary of the distributions."""
    
    stats = {
        "mutant": {
            "mean": float(np.mean(mutant_jsds)),
            "std": float(np.std(mutant_jsds)),
            "min": float(np.min(mutant_jsds)),
            "max": float(np.max(mutant_jsds)),
            "median": float(np.median(mutant_jsds)),
            "n": len(mutant_jsds),
            "values": mutant_jsds.tolist()
        },
        "control": {
            "mean": float(np.mean(control_jsds)),
            "std": float(np.std(control_jsds)),
            "min": float(np.min(control_jsds)),
            "max": float(np.max(control_jsds)),
            "median": float(np.median(control_jsds)),
            "n": len(control_jsds),
            "values": control_jsds.tolist()
        },
        "comparison": {
            "mean_difference": float(np.mean(mutant_jsds) - np.mean(control_jsds)),
            "t_test": None  # Could add statistical test here
        }
    }
    
    with open(output_path, 'w') as f:
        json.dump(stats, f, indent=2)
    
    print(f"Saved statistics to {output_path}")
    
    return stats


def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Plot JSD distributions for mutant and control individuals",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('--mutant-dir', 
                        default='../data/individuals/mutant',
                        help='Directory containing mutant individuals')
    parser.add_argument('--control-dir',
                        default='../data/individuals/control',
                        help='Directory containing control individuals')
    parser.add_argument('--output-dir',
                        default='../plots',
                        help='Output directory for plots and statistics')
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()
    
    print("=" * 60)
    print("JSD DISTRIBUTION COMPARISON")
    print("=" * 60)
    
    # Load mutant individuals
    print("\nLoading mutant individuals...")
    mutant_pattern = os.path.join(args.mutant_dir, "individual_*.json.gz")
    mutant_jsds = load_individuals(mutant_pattern)
    
    # Load control individuals
    print("\nLoading control individuals...")
    control_pattern = os.path.join(args.control_dir, "individual_*.json.gz")
    control_jsds = load_individuals(control_pattern)
    
    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create plot
    print("\nCreating comparison plot...")
    plot_path = os.path.join(args.output_dir, "jsd_distribution_comparison.png")
    create_comparison_plot(mutant_jsds, control_jsds, plot_path)
    
    # Save statistics
    print("\nCalculating statistics...")
    stats_path = os.path.join(args.output_dir, "../data/results/jsd_distributions.json")
    os.makedirs(os.path.dirname(stats_path), exist_ok=True)
    stats = save_statistics(mutant_jsds, control_jsds, stats_path)
    
    # Print summary
    print("\n" + "=" * 60)
    print("SUMMARY STATISTICS")
    print("=" * 60)
    print(f"\nMutant individuals (n={len(mutant_jsds)}):")
    print(f"  Mean JSD: {stats['mutant']['mean']:.4f} ± {stats['mutant']['std']:.4f}")
    print(f"  Range: [{stats['mutant']['min']:.4f}, {stats['mutant']['max']:.4f}]")
    
    print(f"\nControl individuals (n={len(control_jsds)}):")
    print(f"  Mean JSD: {stats['control']['mean']:.4f} ± {stats['control']['std']:.4f}")
    print(f"  Range: [{stats['control']['min']:.4f}, {stats['control']['max']:.4f}]")
    
    print(f"\nDifference (Mutant - Control):")
    print(f"  Mean difference: {stats['comparison']['mean_difference']:.4f}")
    
    print("\n" + "=" * 60)


if __name__ == "__main__":
    main()