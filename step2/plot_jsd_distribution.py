#!/usr/bin/env python3
"""
Plot the distribution of JSD scores across all cells from a snapshot.
Maintains the same visual style as the time-series plots.
"""

import json
import gzip
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from scipy import stats
import os
import sys


def load_snapshot(filepath):
    """Load snapshot data from json.gz file."""
    print(f"Loading snapshot from {filepath}...")
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    print(f"  Loaded {data['metadata']['num_cells']} cells from year {data['metadata']['extracted_year']}")
    return data


def create_jsd_combined_plot(snapshot_data, output_filename, n_bins=200):
    """Create step histogram plot of JSD distribution."""
    # Extract JSD values from all cells
    jsd_values = np.array([cell['jsd'] for cell in snapshot_data['cells']])
    
    # Calculate statistics
    mean_jsd = np.mean(jsd_values)
    std_jsd = np.std(jsd_values)
    median_jsd = np.median(jsd_values)
    p5_jsd = np.percentile(jsd_values, 5)
    p25_jsd = np.percentile(jsd_values, 25)
    p75_jsd = np.percentile(jsd_values, 75)
    p95_jsd = np.percentile(jsd_values, 95)
    
    # Create single plot
    fig = go.Figure()
    
    # Calculate bin size and histogram data
    data_range = np.max(jsd_values) - np.min(jsd_values)
    bin_size = data_range / n_bins
    counts, bin_edges = np.histogram(jsd_values, bins=n_bins)
    
    # Create step plot data
    x_step = []
    y_step = []
    for i in range(len(counts)):
        x_step.extend([bin_edges[i], bin_edges[i+1]])
        y_step.extend([counts[i], counts[i]])
    
    step_trace = go.Scatter(
        x=x_step,
        y=y_step,
        mode='lines',
        name='JSD Distribution',
        line=dict(color='#1f77b4', width=2),
        fill='tozeroy',
        fillcolor='rgba(31, 119, 180, 0.3)',
        showlegend=False,
        hovertemplate='JSD: %{x:.4f}<br>Count: %{y}<extra></extra>'
    )
    fig.add_trace(step_trace)
    
    # Add mean line
    fig.add_vline(
        x=mean_jsd,
        line_dash="solid",
        line_color="red",
        line_width=2,
        annotation_text=f"Mean: {mean_jsd:.4f}",
        annotation_position="top right"
    )
    
    
    # Add statistics annotation
    stats_text = (f"<b>Statistics</b><br>"
                  f"Mean: {mean_jsd:.4f}<br>"
                  f"SD: {std_jsd:.4f}<br>"
                  f"Median: {median_jsd:.4f}<br>"
                  f"5%: {p5_jsd:.4f}<br>"
                  f"95%: {p95_jsd:.4f}")
    
    fig.add_annotation(
        text=stats_text,
        xref="paper", yref="paper",
        x=0.98,
        y=0.97,
        showarrow=False,
        font=dict(size=11, family="Arial"),
        align="right",
        xanchor="right",
        yanchor="top",
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="#333333",
        borderwidth=1
    )
    
    # Update layout
    metadata = snapshot_data['metadata']
    fig.update_layout(
        title=dict(
            text=f"JSD Distribution at Year {metadata['extracted_year']}<br>"
                 f"<sub>{metadata['source_file']} | {metadata['num_cells']} cells</sub>",
            font=dict(size=16)
        ),
        xaxis_title="Jensen-Shannon Divergence",
        yaxis_title="Count",
        template="plotly_white",
        width=1200,
        height=600,
        showlegend=False
    )
    
    # Update axes
    fig.update_xaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray',
        zeroline=True,
        zerolinewidth=1,
        zerolinecolor='gray'
    )
    
    fig.update_yaxes(
        showgrid=True,
        gridwidth=1,
        gridcolor='lightgray'
    )
    
    return fig




def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print("Usage: python plot_jsd_distribution.py <snapshot.json.gz> [n_bins]")
        print("Example: python plot_jsd_distribution.py year50_snapshot.json.gz")
        print("         python plot_jsd_distribution.py year50_snapshot.json.gz 500")
        print("Default: n_bins=200")
        sys.exit(1)
    
    snapshot_file = sys.argv[1]
    n_bins = int(sys.argv[2]) if len(sys.argv) == 3 else 200
    
    # Load snapshot data
    snapshot_data = load_snapshot(snapshot_file)
    
    # Extract base filename for output
    base_filename = os.path.basename(snapshot_file).replace('.json.gz', '')
    
    # Create plots directory if it doesn't exist
    plots_dir = 'plots'
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir)
        print(f"Created directory: {plots_dir}")
    
    # Create and save combined plot
    print(f"\nCreating JSD distribution combined plot with {n_bins} bins...")
    fig_combined = create_jsd_combined_plot(snapshot_data, base_filename, n_bins)
    combined_output = os.path.join(plots_dir, f"{base_filename}_jsd_distribution_{n_bins}bins.png")
    fig_combined.write_image(combined_output, width=1200, height=600, scale=2)
    print(f"  Saved to: {combined_output}")
    
    # Print summary statistics
    jsd_values = [cell['jsd'] for cell in snapshot_data['cells']]
    print("\nJSD Distribution Summary:")
    print(f"  Mean: {np.mean(jsd_values):.4f}")
    print(f"  Std Dev: {np.std(jsd_values):.4f}")
    print(f"  Min: {np.min(jsd_values):.4f}")
    print(f"  Max: {np.max(jsd_values):.4f}")
    print(f"  Median: {np.median(jsd_values):.4f}")
    print(f"  5th percentile: {np.percentile(jsd_values, 5):.4f}")
    print(f"  95th percentile: {np.percentile(jsd_values, 95):.4f}")
    
    print("\nPlot generation complete!")


if __name__ == "__main__":
    main()