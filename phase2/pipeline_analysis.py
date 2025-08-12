#!/usr/bin/env python3
"""
Analysis and visualization functions for the phase2 pipeline.
Adapted to work directly with Cell and PetriDish objects.
"""

import json
import gzip
import numpy as np
import glob
import os
import sys
import plotly.graph_objects as go
from scipy import stats
from typing import List, Dict, Any, Optional, Tuple

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import load_petri_dish, load_all_petri_dishes, get_jsd_array


def plot_jsd_distribution_from_cells(cells: List[Cell], bins: int, output_path: str, 
                                    rate: Optional[float] = None) -> None:
    """
    Create histogram of JSD distribution from Cell objects.
    
    Args:
        cells: List of Cell objects
        bins: Number of bins for histogram
        output_path: Path to save plot
        rate: Methylation rate (optional, for display)
    """
    print(f"\nPlotting JSD distribution...")
    
    # Extract JSD values directly from Cell objects
    jsd_values = np.array([cell.JSD for cell in cells])
    
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
    year = cells[0].age if cells else 'unknown'
    
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


def get_mean_jsds_from_petri_dishes(dishes: List[PetriDish]) -> np.ndarray:
    """
    Extract mean JSD values from a list of PetriDish objects.
    
    Args:
        dishes: List of PetriDish objects
    
    Returns:
        Array of mean JSD values (one per dish)
    """
    mean_jsds = []
    for petri in dishes:
        jsd_values = [cell.JSD for cell in petri.cells]
        mean_jsd = np.mean(jsd_values) if jsd_values else 0.0
        mean_jsds.append(mean_jsd)
    
    return np.array(mean_jsds)


def analyze_populations_from_dishes(mutant_dishes: List[PetriDish], 
                                   control1_dishes: List[PetriDish],
                                   control2_dishes: List[PetriDish],
                                   output_dir: str) -> Dict[str, Any]:
    """
    Analyze and compare three population groups using PetriDish objects.
    
    Args:
        mutant_dishes: List of mutant PetriDish objects
        control1_dishes: List of control1 PetriDish objects
        control2_dishes: List of control2 PetriDish objects
        output_dir: Directory to save results
    
    Returns:
        Dictionary with results and plot paths
    """
    print(f"\nAnalyzing populations from PetriDish objects...")
    
    # Get mean JSD values for each group
    mutant_jsds = get_mean_jsds_from_petri_dishes(mutant_dishes)
    control1_jsds = get_mean_jsds_from_petri_dishes(control1_dishes)
    control2_jsds = get_mean_jsds_from_petri_dishes(control2_dishes)
    
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
    create_comparison_plot_from_jsds(mutant_jsds, control1_jsds, control2_jsds, plot_path)
    
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


def create_comparison_plot_from_jsds(mutant_jsds: np.ndarray, 
                                    control1_jsds: np.ndarray, 
                                    control2_jsds: np.ndarray, 
                                    output_path: str) -> None:
    """
    Create a clean scatter plot comparing all three distributions.
    
    Args:
        mutant_jsds: Array of mean JSD values for mutant group
        control1_jsds: Array of mean JSD values for control1 group
        control2_jsds: Array of mean JSD values for control2 group
        output_path: Path to save the plot
    """
    print(f"\n  Creating comparison plot...")
    
    # Set random seed for consistent jitter
    np.random.seed(42)
    
    # Create figure
    fig = go.Figure()
    
    # Add mutant points with jitter
    x_mutant = np.ones(len(mutant_jsds)) * 0 + np.random.normal(0, 0.02, len(mutant_jsds))
    fig.add_trace(go.Scatter(
        x=x_mutant,
        y=mutant_jsds,
        mode='markers',
        name='Mutant',
        marker=dict(color='#1f77b4', size=10, opacity=0.6)
    ))
    
    # Add control1 points with jitter
    x_control1 = np.ones(len(control1_jsds)) * 1 + np.random.normal(0, 0.02, len(control1_jsds))
    fig.add_trace(go.Scatter(
        x=x_control1,
        y=control1_jsds,
        mode='markers',
        name='Control1',
        marker=dict(color='#ff7f0e', size=10, opacity=0.6)
    ))
    
    # Add control2 points with jitter
    x_control2 = np.ones(len(control2_jsds)) * 2 + np.random.normal(0, 0.02, len(control2_jsds))
    fig.add_trace(go.Scatter(
        x=x_control2,
        y=control2_jsds,
        mode='markers',
        name='Control2',
        marker=dict(color='#2ca02c', size=10, opacity=0.6)
    ))
    
    # Calculate means and quantiles for horizontal lines
    mutant_mean = np.mean(mutant_jsds)
    control1_mean = np.mean(control1_jsds)
    control2_mean = np.mean(control2_jsds)
    
    # Calculate quantiles
    mutant_q5, mutant_q25, mutant_q75, mutant_q95 = np.percentile(mutant_jsds, [5, 25, 75, 95])
    control1_q5, control1_q25, control1_q75, control1_q95 = np.percentile(control1_jsds, [5, 25, 75, 95])
    control2_q5, control2_q25, control2_q75, control2_q95 = np.percentile(control2_jsds, [5, 25, 75, 95])
    
    # Add quantile lines for mutant
    # 5-95 range (very light, dotted)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mutant_q5, y1=mutant_q5,
                  line=dict(color="#1f77b4", width=1, dash="dot"), opacity=0.3)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mutant_q95, y1=mutant_q95,
                  line=dict(color="#1f77b4", width=1, dash="dot"), opacity=0.3)
    
    # 25-75 range (slightly darker, dashed)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mutant_q25, y1=mutant_q25,
                  line=dict(color="#1f77b4", width=1.5, dash="dash"), opacity=0.5)
    fig.add_shape(type="line", x0=-0.15, x1=0.15, y0=mutant_q75, y1=mutant_q75,
                  line=dict(color="#1f77b4", width=1.5, dash="dash"), opacity=0.5)
    
    # Add quantile lines for control1
    # 5-95 range (very light, dotted)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control1_q5, y1=control1_q5,
                  line=dict(color="#ff7f0e", width=1, dash="dot"), opacity=0.3)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control1_q95, y1=control1_q95,
                  line=dict(color="#ff7f0e", width=1, dash="dot"), opacity=0.3)
    
    # 25-75 range (slightly darker, dashed)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control1_q25, y1=control1_q25,
                  line=dict(color="#ff7f0e", width=1.5, dash="dash"), opacity=0.5)
    fig.add_shape(type="line", x0=0.85, x1=1.15, y0=control1_q75, y1=control1_q75,
                  line=dict(color="#ff7f0e", width=1.5, dash="dash"), opacity=0.5)
    
    # Add quantile lines for control2
    # 5-95 range (very light, dotted)
    fig.add_shape(type="line", x0=1.85, x1=2.15, y0=control2_q5, y1=control2_q5,
                  line=dict(color="#2ca02c", width=1, dash="dot"), opacity=0.3)
    fig.add_shape(type="line", x0=1.85, x1=2.15, y0=control2_q95, y1=control2_q95,
                  line=dict(color="#2ca02c", width=1, dash="dot"), opacity=0.3)
    
    # 25-75 range (slightly darker, dashed)
    fig.add_shape(type="line", x0=1.85, x1=2.15, y0=control2_q25, y1=control2_q25,
                  line=dict(color="#2ca02c", width=1.5, dash="dash"), opacity=0.5)
    fig.add_shape(type="line", x0=1.85, x1=2.15, y0=control2_q75, y1=control2_q75,
                  line=dict(color="#2ca02c", width=1.5, dash="dash"), opacity=0.5)
    
    # Add mean lines (solid, prominent) - these go last to be on top
    fig.add_shape(type="line", x0=-0.2, x1=0.2, y0=mutant_mean, y1=mutant_mean,
                  line=dict(color="#1f77b4", width=3))
    fig.add_shape(type="line", x0=0.8, x1=1.2, y0=control1_mean, y1=control1_mean,
                  line=dict(color="#ff7f0e", width=3))
    fig.add_shape(type="line", x0=1.8, x1=2.2, y0=control2_mean, y1=control2_mean,
                  line=dict(color="#2ca02c", width=3))
    
    # Calculate comprehensive statistics for annotations
    # Mutant stats
    mutant_median = np.median(mutant_jsds)
    mutant_std = np.std(mutant_jsds)
    mutant_cv = mutant_std / mutant_mean if mutant_mean != 0 else 0
    mutant_mad = np.median(np.abs(mutant_jsds - mutant_median))
    
    # Control1 stats
    control1_median = np.median(control1_jsds)
    control1_std = np.std(control1_jsds)
    control1_cv = control1_std / control1_mean if control1_mean != 0 else 0
    control1_mad = np.median(np.abs(control1_jsds - control1_median))
    
    # Control2 stats
    control2_median = np.median(control2_jsds)
    control2_std = np.std(control2_jsds)
    control2_cv = control2_std / control2_mean if control2_mean != 0 else 0
    control2_mad = np.median(np.abs(control2_jsds - control2_median))
    
    # Determine y position for annotations (below the lowest points)
    all_jsds = np.concatenate([mutant_jsds, control1_jsds, control2_jsds])
    y_min = np.min(all_jsds)
    y_max = np.max(all_jsds)
    y_range = y_max - y_min
    annotation_y = y_min - 0.08 * y_range  # Slightly more space for more lines
    
    # Add comprehensive statistics text below each group
    fig.add_annotation(
        x=0, y=annotation_y,
        text=(f"Mean: {mutant_mean:.4f}<br>"
              f"Median: {mutant_median:.4f}<br>"
              f"SD: {mutant_std:.4f}<br>"
              f"CV: {mutant_cv:.3f}<br>"
              f"MAD: {mutant_mad:.4f}<br>"
              f"5%: {mutant_q5:.4f}<br>"
              f"95%: {mutant_q95:.4f}"),
        showarrow=False,
        font=dict(size=9, color="#1f77b4"),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    fig.add_annotation(
        x=1, y=annotation_y,
        text=(f"Mean: {control1_mean:.4f}<br>"
              f"Median: {control1_median:.4f}<br>"
              f"SD: {control1_std:.4f}<br>"
              f"CV: {control1_cv:.3f}<br>"
              f"MAD: {control1_mad:.4f}<br>"
              f"5%: {control1_q5:.4f}<br>"
              f"95%: {control1_q95:.4f}"),
        showarrow=False,
        font=dict(size=9, color="#ff7f0e"),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    fig.add_annotation(
        x=2, y=annotation_y,
        text=(f"Mean: {control2_mean:.4f}<br>"
              f"Median: {control2_median:.4f}<br>"
              f"SD: {control2_std:.4f}<br>"
              f"CV: {control2_cv:.3f}<br>"
              f"MAD: {control2_mad:.4f}<br>"
              f"5%: {control2_q5:.4f}<br>"
              f"95%: {control2_q95:.4f}"),
        showarrow=False,
        font=dict(size=9, color="#2ca02c"),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    # Clean layout with proper spacing
    fig.update_layout(
        title=dict(
            text="Mean JSD Distribution Comparison",
            font=dict(size=16),
            x=0.5,
            xanchor='center'
        ),
        xaxis=dict(
            tickmode='array',
            tickvals=[0, 1, 2],
            ticktext=['Mutant', 'Control1', 'Control2'],
            range=[-0.5, 2.5],
            showgrid=False,
            title=""
        ),
        yaxis=dict(
            title='Mean JSD per Individual',
            showgrid=True,
            gridcolor='lightgray',
            range=[annotation_y - 0.02, y_max + 0.01]
        ),
        plot_bgcolor='white',
        showlegend=False,
        height=700,  # Increased height for stats
        width=700,
        margin=dict(l=80, r=40, t=60, b=150)  # Increased bottom margin for stats
    )
    
    # Save PNG only (no HTML)
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_image(output_path, width=700, height=700, scale=2)
    print(f"    Saved plot to {output_path}")


def plot_cell_level_distributions(mutant_dishes: List[PetriDish],
                                 control1_dishes: List[PetriDish],
                                 control2_dishes: List[PetriDish],
                                 output_path: str) -> None:
    """
    Create a plot showing cell-level JSD distributions (all cells from all individuals).
    
    Args:
        mutant_dishes: List of mutant PetriDish objects
        control1_dishes: List of control1 PetriDish objects
        control2_dishes: List of control2 PetriDish objects
        output_path: Path to save the plot
    """
    print(f"\n  Creating cell-level distribution plot...")
    
    # Collect all cell JSDs from each group
    mutant_cell_jsds = np.concatenate([get_jsd_array(p) for p in mutant_dishes])
    control1_cell_jsds = np.concatenate([get_jsd_array(p) for p in control1_dishes])
    control2_cell_jsds = np.concatenate([get_jsd_array(p) for p in control2_dishes])
    
    print(f"    Mutant: {len(mutant_cell_jsds)} cells total")
    print(f"    Control1: {len(control1_cell_jsds)} cells total")
    print(f"    Control2: {len(control2_cell_jsds)} cells total")
    
    # Create figure with overlapping histograms
    fig = go.Figure()
    
    # Add histograms for each group
    fig.add_trace(go.Histogram(
        x=mutant_cell_jsds,
        name='Mutant',
        opacity=0.5,
        marker_color='#1f77b4',
        nbinsx=100
    ))
    
    fig.add_trace(go.Histogram(
        x=control1_cell_jsds,
        name='Control1',
        opacity=0.5,
        marker_color='#ff7f0e',
        nbinsx=100
    ))
    
    fig.add_trace(go.Histogram(
        x=control2_cell_jsds,
        name='Control2',
        opacity=0.5,
        marker_color='#2ca02c',
        nbinsx=100
    ))
    
    # Update layout
    fig.update_layout(
        title=dict(
            text="Cell-Level JSD Distributions<br>"
                 "<sub>All cells from all individuals</sub>",
            font=dict(size=16)
        ),
        xaxis_title="JSD Score",
        yaxis_title="Count",
        barmode='overlay',
        template="plotly_white",
        width=1200,
        height=600,
        showlegend=True,
        legend=dict(
            x=0.7,
            y=0.95,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="#333333",
            borderwidth=1
        )
    )
    
    # Save
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_image(output_path, width=1200, height=600, scale=2)
    print(f"    Saved plot to {output_path}")