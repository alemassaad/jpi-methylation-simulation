#!/usr/bin/env python3
"""
Load mixed and control distributions and plot them together.
Simple visualization to compare the two sets of 30 mean JSD values.
"""

import json
import gzip
import numpy as np
import glob
import os
import plotly.graph_objects as go


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


def create_comparison_plot(mixed_jsds, control_jsds, output_path="jsd_distribution_comparison.png"):
    """Create a simple plot comparing mixed and control distributions."""
    
    # Create figure
    fig = go.Figure()
    
    # Add mixed points
    x_mixed = np.ones(len(mixed_jsds)) * 0 + np.random.normal(0, 0.02, len(mixed_jsds))
    fig.add_trace(go.Scatter(
        x=x_mixed,
        y=mixed_jsds,
        mode='markers',
        name='Mixed',
        marker=dict(color='#1f77b4', size=10, opacity=0.6)
    ))
    
    # Add control points
    x_control = np.ones(len(control_jsds)) * 1 + np.random.normal(0, 0.02, len(control_jsds))
    fig.add_trace(go.Scatter(
        x=x_control,
        y=control_jsds,
        mode='markers',
        name='Control',
        marker=dict(color='#ff7f0e', size=10, opacity=0.6)
    ))
    
    # Calculate statistics
    mixed_mean = np.mean(mixed_jsds)
    control_mean = np.mean(control_jsds)
    
    # Calculate quantiles
    mixed_q5, mixed_q25, mixed_q75, mixed_q95 = np.percentile(mixed_jsds, [5, 25, 75, 95])
    control_q5, control_q25, control_q75, control_q95 = np.percentile(control_jsds, [5, 25, 75, 95])
    
    # Add quantile lines for mixed (subtle)
    # 5-95 range (very light)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mixed_q5, y1=mixed_q5,
                  line=dict(color="#1f77b4", width=1, dash="dot"), opacity=0.3)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mixed_q95, y1=mixed_q95,
                  line=dict(color="#1f77b4", width=1, dash="dot"), opacity=0.3)
    
    # 25-75 range (slightly darker)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mixed_q25, y1=mixed_q25,
                  line=dict(color="#1f77b4", width=1.5, dash="dash"), opacity=0.5)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mixed_q75, y1=mixed_q75,
                  line=dict(color="#1f77b4", width=1.5, dash="dash"), opacity=0.5)
    
    # Add quantile lines for control (subtle)
    # 5-95 range (very light)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control_q5, y1=control_q5,
                  line=dict(color="#ff7f0e", width=1, dash="dot"), opacity=0.3)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control_q95, y1=control_q95,
                  line=dict(color="#ff7f0e", width=1, dash="dot"), opacity=0.3)
    
    # 25-75 range (slightly darker)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control_q25, y1=control_q25,
                  line=dict(color="#ff7f0e", width=1.5, dash="dash"), opacity=0.5)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control_q75, y1=control_q75,
                  line=dict(color="#ff7f0e", width=1.5, dash="dash"), opacity=0.5)
    
    # Add mean lines (solid, prominent)
    fig.add_shape(type="line", x0=-0.2, x1=0.2, y0=mixed_mean, y1=mixed_mean,
                  line=dict(color="#1f77b4", width=3))
    fig.add_shape(type="line", x0=0.8, x1=1.2, y0=control_mean, y1=control_mean,
                  line=dict(color="#ff7f0e", width=3))
    
    # Calculate statistics
    mixed_std = np.std(mixed_jsds)
    control_std = np.std(control_jsds)
    mixed_cv = mixed_std / mixed_mean
    control_cv = control_std / control_mean
    
    # Add statistics text
    fig.add_annotation(
        x=0, y=np.min(mixed_jsds) - 0.002,
        text=f"Mean: {mixed_mean:.4f}<br>Std: {mixed_std:.4f}<br>CV: {mixed_cv:.3f}",
        showarrow=False,
        font=dict(size=11, color="#1f77b4"),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    fig.add_annotation(
        x=1, y=np.min(control_jsds) - 0.002,
        text=f"Mean: {control_mean:.4f}<br>Std: {control_std:.4f}<br>CV: {control_cv:.3f}",
        showarrow=False,
        font=dict(size=11, color="#ff7f0e"),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    # Clean layout with more space
    fig.update_layout(
        title=dict(
            text="Mean JSD Distribution at Year 60<br><sub>30 individuals per group | Methylation rate: 0.5%/year</sub>",
            font=dict(size=16),
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            tickmode='array',
            tickvals=[0, 1],
            ticktext=['Mixed', 'Control'],
            range=[-0.5, 1.5],
            showgrid=False
        ),
        yaxis=dict(
            title='Mean JSD',
            showgrid=True,
            gridcolor='lightgray',
            range=[min(np.min(mixed_jsds), np.min(control_jsds)) - 0.01,
                   max(np.max(mixed_jsds), np.max(control_jsds)) + 0.002]
        ),
        plot_bgcolor='white',
        showlegend=False,
        height=650,
        width=500,
        margin=dict(l=70, r=40, t=80, b=120)
    )
    
    # Save plot
    os.makedirs("plots", exist_ok=True)
    output_path = os.path.join("plots", output_path)
    fig.write_image(output_path)
    print(f"\nPlot saved to: {output_path}")
    
    # Also save as HTML for interactive viewing
    html_path = output_path.replace('.png', '.html')
    fig.write_html(html_path)
    print(f"Interactive version saved to: {html_path}")
    
    return fig


def print_statistics(mixed_jsds, control_jsds):
    """Print summary statistics for both distributions."""
    print("\n" + "="*60)
    print("SUMMARY STATISTICS")
    print("="*60)
    
    print(f"\nMixed individuals (n={len(mixed_jsds)}):")
    print(f"  Mean: {np.mean(mixed_jsds):.4f}")
    print(f"  Std:  {np.std(mixed_jsds):.4f}")
    print(f"  Min:  {np.min(mixed_jsds):.4f}")
    print(f"  Max:  {np.max(mixed_jsds):.4f}")
    
    print(f"\nControl individuals (n={len(control_jsds)}):")
    print(f"  Mean: {np.mean(control_jsds):.4f}")
    print(f"  Std:  {np.std(control_jsds):.4f}")
    print(f"  Min:  {np.min(control_jsds):.4f}")
    print(f"  Max:  {np.max(control_jsds):.4f}")
    
    # Simple difference
    mean_diff = np.mean(mixed_jsds) - np.mean(control_jsds)
    print(f"\nMean difference (Mixed - Control): {mean_diff:.4f}")
    
    # Effect size (Cohen's d)
    pooled_std = np.sqrt((np.std(mixed_jsds)**2 + np.std(control_jsds)**2) / 2)
    cohens_d = mean_diff / pooled_std
    print(f"Cohen's d effect size: {cohens_d:.3f}")


def main():
    """Load distributions and create comparison plot."""
    print("Loading mixed and control distributions...")
    
    # Load mixed individuals
    mixed_pattern = "individuals/individual_*.json.gz"
    mixed_jsds = load_individuals(mixed_pattern)
    
    # Load control individuals
    control_pattern = "control_individuals/control_individual_*.json.gz"
    control_jsds = load_individuals(control_pattern)
    
    # Create plot
    create_comparison_plot(mixed_jsds, control_jsds)
    
    # Print statistics
    print_statistics(mixed_jsds, control_jsds)
    
    # Save raw data for later analysis
    output_data = {
        "mixed": {
            "values": mixed_jsds.tolist(),
            "n": len(mixed_jsds),
            "mean": float(np.mean(mixed_jsds)),
            "std": float(np.std(mixed_jsds))
        },
        "control": {
            "values": control_jsds.tolist(),
            "n": len(control_jsds),
            "mean": float(np.mean(control_jsds)),
            "std": float(np.std(control_jsds))
        }
    }
    
    with open("jsd_distributions.json", 'w') as f:
        json.dump(output_data, f, indent=2)
    print(f"\nRaw data saved to: jsd_distributions.json")


if __name__ == "__main__":
    main()