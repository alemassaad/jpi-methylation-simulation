#!/usr/bin/env python3
"""
Analysis and visualization functions for the step23 pipeline.
"""

import json
import gzip
import numpy as np
import glob
import os
import plotly.graph_objects as go
from scipy import stats


def plot_jsd_distribution(cells, bins, output_path, rate=None):
    """
    Create histogram of JSD distribution matching step2 style.
    
    Args:
        cells: List of cell dictionaries
        bins: Number of bins for histogram
        output_path: Path to save plot
        rate: Methylation rate (optional, for display)
    """
    print(f"\nPlotting JSD distribution...")
    
    # Extract JSD values
    jsd_values = np.array([cell['jsd'] for cell in cells])
    
    # Calculate statistics
    mean_jsd = np.mean(jsd_values)
    std_jsd = np.std(jsd_values)
    median_jsd = np.median(jsd_values)
    cv_jsd = std_jsd / mean_jsd if mean_jsd > 0 else 0  # Coefficient of variation
    mad_jsd = np.median(np.abs(jsd_values - median_jsd))  # Median Absolute Deviation
    p5_jsd = np.percentile(jsd_values, 5)
    p25_jsd = np.percentile(jsd_values, 25)
    p75_jsd = np.percentile(jsd_values, 75)
    p95_jsd = np.percentile(jsd_values, 95)
    
    # Create figure
    fig = go.Figure()
    
    # Calculate histogram data for step plot
    counts, bin_edges = np.histogram(jsd_values, bins=bins)
    
    # Create step plot data (like the original)
    x_step = []
    y_step = []
    for i in range(len(counts)):
        x_step.extend([bin_edges[i], bin_edges[i+1]])
        y_step.extend([counts[i], counts[i]])
    
    # Add step histogram with fill
    fig.add_trace(go.Scatter(
        x=x_step,
        y=y_step,
        mode='lines',
        name='JSD Distribution',
        line=dict(color='#1f77b4', width=2),
        fill='tozeroy',
        fillcolor='rgba(31, 119, 180, 0.3)',
        showlegend=False,
        hovertemplate='JSD: %{x:.4f}<br>Count: %{y}<extra></extra>'
    ))
    
    # Add mean line
    fig.add_vline(
        x=mean_jsd,
        line_dash="solid",
        line_color="red",
        line_width=2,
        annotation_text=f"Mean: {mean_jsd:.4f}",
        annotation_position="top right"
    )
    
    # Add statistics box in top right corner
    stats_text = (f"<b>Statistics</b><br>"
                  f"Mean: {mean_jsd:.4f}<br>"
                  f"Median: {median_jsd:.4f}<br>"
                  f"SD: {std_jsd:.4f}<br>"
                  f"CV: {cv_jsd:.3f}<br>"
                  f"MAD: {mad_jsd:.4f}<br>"
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
    year = cells[0]['age'] if cells else 'unknown'
    
    # Format rate as percentage if provided
    rate_text = ""
    if rate is not None:
        rate_percentage = rate * 100
        rate_text = f" | {rate_percentage:.1f}% methylation rate"
    
    fig.update_layout(
        title=dict(
            text=f"JSD Distribution at Year {year}<br>"
                 f"<sub>{len(cells)} cells{rate_text}</sub>",
            font=dict(size=16)
        ),
        xaxis_title="JSD Score",
        yaxis_title="Count",
        template="plotly_white",
        width=1200,
        height=600,
        showlegend=False
    )
    
    # Update axes with grid
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
    
    # Save
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_image(output_path, width=1200, height=600, scale=2)
    print(f"  Saved plot to {output_path}")


def load_individuals_from_dir(directory):
    """Load all individuals from a directory and extract mean JSD values."""
    pattern = os.path.join(directory, "individual_*.json.gz")
    files = sorted(glob.glob(pattern))
    
    if not files:
        raise ValueError(f"No individual files found in {directory}")
    
    print(f"  Loading {len(files)} individuals from {directory}")
    
    mean_jsds = []
    for filepath in files:
        with gzip.open(filepath, 'rt') as f:
            data = json.load(f)
        mean_jsd = data['metadata']['mean_jsd']
        mean_jsds.append(mean_jsd)
    
    return np.array(mean_jsds)


def analyze_populations(mutant_dir, control1_dir, control2_dir, output_dir):
    """
    Analyze and compare the three population groups.
    
    Returns:
        Dictionary with results and plot paths
    """
    print(f"\nAnalyzing populations...")
    
    # Load mean JSD values for each group
    mutant_jsds = load_individuals_from_dir(mutant_dir)
    control1_jsds = load_individuals_from_dir(control1_dir)
    control2_jsds = load_individuals_from_dir(control2_dir)
    
    print(f"  Mutant: {len(mutant_jsds)} individuals")
    print(f"  Control1: {len(control1_jsds)} individuals")
    print(f"  Control2: {len(control2_jsds)} individuals")
    
    # Calculate statistics
    stats_dict = {
        "mutant": {
            "mean": float(np.mean(mutant_jsds)),
            "std": float(np.std(mutant_jsds)),
            "median": float(np.median(mutant_jsds)),
            "min": float(np.min(mutant_jsds)),
            "max": float(np.max(mutant_jsds)),
            "n": len(mutant_jsds)
        },
        "control1": {
            "mean": float(np.mean(control1_jsds)),
            "std": float(np.std(control1_jsds)),
            "median": float(np.median(control1_jsds)),
            "min": float(np.min(control1_jsds)),
            "max": float(np.max(control1_jsds)),
            "n": len(control1_jsds)
        },
        "control2": {
            "mean": float(np.mean(control2_jsds)),
            "std": float(np.std(control2_jsds)),
            "median": float(np.median(control2_jsds)),
            "min": float(np.min(control2_jsds)),
            "max": float(np.max(control2_jsds)),
            "n": len(control2_jsds)
        }
    }
    
    # Statistical tests
    print("\n  Statistical comparisons:")
    
    # Mutant vs Control1
    t_stat_mc1, p_value_mc1 = stats.ttest_ind(mutant_jsds, control1_jsds)
    print(f"    Mutant vs Control1: t={t_stat_mc1:.3f}, p={p_value_mc1:.6f}")
    
    # Mutant vs Control2
    t_stat_mc2, p_value_mc2 = stats.ttest_ind(mutant_jsds, control2_jsds)
    print(f"    Mutant vs Control2: t={t_stat_mc2:.3f}, p={p_value_mc2:.6f}")
    
    # Control1 vs Control2
    t_stat_c1c2, p_value_c1c2 = stats.ttest_ind(control1_jsds, control2_jsds)
    print(f"    Control1 vs Control2: t={t_stat_c1c2:.3f}, p={p_value_c1c2:.6f}")
    
    stats_dict["comparisons"] = {
        "mutant_vs_control1": {
            "t_statistic": float(t_stat_mc1),
            "p_value": float(p_value_mc1)
        },
        "mutant_vs_control2": {
            "t_statistic": float(t_stat_mc2),
            "p_value": float(p_value_mc2)
        },
        "control1_vs_control2": {
            "t_statistic": float(t_stat_c1c2),
            "p_value": float(p_value_c1c2)
        }
    }
    
    # Create comparison plot
    plot_path = os.path.join(output_dir, "jsd_comparison.png")
    create_comparison_plot(mutant_jsds, control1_jsds, control2_jsds, plot_path)
    
    # Save statistics
    stats_path = os.path.join(output_dir, "statistics.json")
    os.makedirs(output_dir, exist_ok=True)
    with open(stats_path, 'w') as f:
        json.dump(stats_dict, f, indent=2)
    print(f"\n  Saved statistics to {stats_path}")
    
    # Save raw distributions
    distributions_path = os.path.join(output_dir, "jsd_distributions.json")
    distributions = {
        "mutant": mutant_jsds.tolist(),
        "control1": control1_jsds.tolist(),
        "control2": control2_jsds.tolist()
    }
    with open(distributions_path, 'w') as f:
        json.dump(distributions, f, indent=2)
    
    return {
        "statistics": stats_dict,
        "plot_path": plot_path,
        "stats_path": stats_path,
        "distributions_path": distributions_path
    }


def create_comparison_plot(mutant_jsds, control1_jsds, control2_jsds, output_path):
    """Create a simple scatter plot comparing all three distributions."""
    
    print(f"\n  Creating comparison plot...")
    
    # Create figure
    fig = go.Figure()
    
    # Add mutant points with jitter
    x_mutant = np.ones(len(mutant_jsds)) * 0 + np.random.normal(0, 0.02, len(mutant_jsds))
    fig.add_trace(go.Scatter(
        x=x_mutant,
        y=mutant_jsds,
        mode='markers',
        name='Mutant (Decile-based)',
        marker=dict(color='#1f77b4', size=10, opacity=0.6),
        showlegend=True
    ))
    
    # Add control1 points with jitter
    x_control1 = np.ones(len(control1_jsds)) * 1 + np.random.normal(0, 0.02, len(control1_jsds))
    fig.add_trace(go.Scatter(
        x=x_control1,
        y=control1_jsds,
        mode='markers',
        name='Control1 (Uniform)',
        marker=dict(color='#ff7f0e', size=10, opacity=0.6),
        showlegend=True
    ))
    
    # Add control2 points with jitter
    x_control2 = np.ones(len(control2_jsds)) * 2 + np.random.normal(0, 0.02, len(control2_jsds))
    fig.add_trace(go.Scatter(
        x=x_control2,
        y=control2_jsds,
        mode='markers',
        name='Control2 (Pure Year 60)',
        marker=dict(color='#2ca02c', size=10, opacity=0.6),
        showlegend=True
    ))
    
    # Add mean lines
    fig.add_hline(y=np.mean(mutant_jsds), line_dash="dash", line_color="#1f77b4", 
                  annotation_text=f"Mutant mean: {np.mean(mutant_jsds):.4f}",
                  annotation_position="left")
    fig.add_hline(y=np.mean(control1_jsds), line_dash="dash", line_color="#ff7f0e",
                  annotation_text=f"Control1 mean: {np.mean(control1_jsds):.4f}",
                  annotation_position="right")
    fig.add_hline(y=np.mean(control2_jsds), line_dash="dash", line_color="#2ca02c",
                  annotation_text=f"Control2 mean: {np.mean(control2_jsds):.4f}",
                  annotation_position="right")
    
    # Update layout
    fig.update_layout(
        title="JSD Distribution Comparison: Three Groups",
        xaxis=dict(
            tickmode='array',
            tickvals=[0, 1, 2],
            ticktext=['Mutant<br>(Decile-based)', 'Control1<br>(Uniform)', 'Control2<br>(Pure Year 60)'],
            range=[-0.5, 2.5]
        ),
        yaxis=dict(title="Mean JSD per Individual"),
        showlegend=True,
        width=1000,
        height=600,
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Update axes
    fig.update_xaxes(showgrid=False, zeroline=False)
    fig.update_yaxes(showgrid=True, gridcolor='lightgray', zeroline=True, zerolinecolor='lightgray')
    
    # Save PNG only (no HTML)
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_image(output_path, width=1000, height=600, scale=2)
    print(f"    Saved plot to {output_path}")