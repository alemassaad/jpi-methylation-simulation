#!/usr/bin/env python3
"""
Histogram plotting function that exactly matches the original plotting logic.
"""

import numpy as np
import plotly.graph_objects as go
from typing import Optional, Tuple, List, Dict, Any
import os


def plot_histogram_original(
    data: np.ndarray,
    # Basic plot configuration
    title: str,
    xlabel: str,
    ylabel: str = "Count",
    
    # Styling (matches original)
    histogram_color: str = "#d62728",  # Red for methylation, blue for JSD
    mean_line_color: str = "darkblue",  # Dark blue for mean
    
    # Histogram configuration
    bins: int = 200,
    
    # Statistics
    custom_stats: Optional[Dict[str, float]] = None,
    
    # Output
    output_path: str = None,
    
    # Additional metadata for subtitle
    year: Optional[int] = None,
    n_cells: Optional[int] = None,
    n_genes: Optional[int] = None,
    gene_rate_groups: Optional[List[Tuple[int, float]]] = None,
    
    # Plot type identifier
    plot_type: str = "methylation",  # "methylation" or "jsd"
) -> go.Figure:
    """
    Create a histogram plot matching the exact original style.
    
    Args:
        data: 1D numpy array of values to plot
        title: Main plot title (will be formatted based on plot_type)
        xlabel: X-axis label
        ylabel: Y-axis label
        histogram_color: Color for histogram
        mean_line_color: Color for mean vertical line
        bins: Number of histogram bins
        custom_stats: Pre-calculated statistics (if None, will calculate)
        output_path: Path to save the plot
        year: Year for display
        n_cells: Number of cells for subtitle (mutually exclusive with n_genes)
        n_genes: Number of genes for subtitle (mutually exclusive with n_cells)
        gene_rate_groups: Gene rate groups for subtitle
        plot_type: Type of plot ("methylation" or "jsd")
        
    Returns:
        Plotly figure object
    """

    # Validate that n_cells and n_genes are mutually exclusive
    if n_cells is not None and n_genes is not None:
        raise ValueError("Cannot specify both n_cells and n_genes - use only one")

    # Calculate statistics if not provided
    if custom_stats is None:
        mean_val = np.mean(data)
        std_val = np.std(data)
        median_val = np.median(data)
        cv_val = std_val / mean_val if mean_val > 0 else 0
        mad_val = np.median(np.abs(data - median_val))
        p5_val = np.percentile(data, 5)
        p95_val = np.percentile(data, 95)
    else:
        mean_val = custom_stats.get('mean', np.mean(data))
        std_val = custom_stats.get('std', np.std(data))
        median_val = custom_stats.get('median', np.median(data))
        cv_val = custom_stats.get('cv', 0)
        mad_val = custom_stats.get('mad', 0)
        p5_val = custom_stats.get('p5', np.percentile(data, 5))
        p95_val = custom_stats.get('p95', np.percentile(data, 95))
    
    # Create figure
    fig = go.Figure()
    
    # Calculate histogram data for step plot
    counts, bin_edges = np.histogram(data, bins=bins)
    
    # Create step plot data (exactly as original)
    x_step = []
    y_step = []
    for i in range(len(counts)):
        x_step.extend([bin_edges[i], bin_edges[i+1]])
        y_step.extend([counts[i], counts[i]])
    
    # Determine fill color based on plot type
    if plot_type == "methylation":
        fillcolor = 'rgba(214, 39, 40, 0.3)'  # Red with transparency
        hover_label = 'Methylation'
    else:  # jsd
        fillcolor = 'rgba(31, 119, 180, 0.3)'  # Blue with transparency
        hover_label = 'JSD'
    
    # Add step histogram with fill (exactly as original)
    fig.add_trace(go.Scatter(
        x=x_step,
        y=y_step,
        mode='lines',
        name=f'{title} Distribution',
        line=dict(color=histogram_color, width=2),
        fill='tozeroy',
        fillcolor=fillcolor,
        showlegend=False,
        hovertemplate=f'{hover_label}: %{{x:.4f}}<br>Count: %{{y}}<extra></extra>'
    ))
    
    # Add mean line (exactly as original)
    fig.add_vline(
        x=mean_val,
        line_dash="solid",
        line_color=mean_line_color,
        line_width=2,
        annotation_text=f"Mean: {mean_val:.4f}",
        annotation_position="top right"
    )
    
    # Add statistics box in top right corner (exactly as original)
    stats_text = (f"<b>Statistics</b><br>"
                  f"Mean: {mean_val:.4f}<br>"
                  f"Median: {median_val:.4f}<br>"
                  f"SD: {std_val:.4f}<br>"
                  f"CV: {cv_val:.3f}<br>"
                  f"MAD: {mad_val:.4f}<br>"
                  f"5%: {p5_val:.4f}<br>"
                  f"95%: {p95_val:.4f}")
    
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
    
    # Update layout (exactly as original)
    display_year = year if year is not None else 'unknown'
    
    # Format gene rate groups for display if provided
    rate_text = ""
    if gene_rate_groups is not None:
        rates = set(rate for _, rate in gene_rate_groups)
        if len(rates) == 1:
            rate_percentage = list(rates)[0] * 100
            rate_text = f" | {rate_percentage:.1f}% methylation rate"
        else:
            rate_text = f" | {len(gene_rate_groups)} gene rate groups"

    # Generate count text based on whether it's cells or genes
    if n_genes is not None:
        count_text = f"{n_genes} genes"
    elif n_cells is not None:
        count_text = f"{n_cells} cells"
    else:
        count_text = ""
    
    fig.update_layout(
        title=dict(
            text=f"{title} Distribution at Year {display_year}<br>"
                 f"<sub>{count_text}{rate_text}</sub>",
            font=dict(size=16)
        ),
        xaxis_title=xlabel,
        yaxis_title=ylabel,
        template="plotly_white",
        width=1200,
        height=600,
        showlegend=False
    )
    
    # Update axes with grid (exactly as original)
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
    if output_path:
        if os.path.dirname(output_path):
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig.write_image(output_path)
        print(f"  Saved plot to: {output_path}")
    
    return fig


def plot_cell_jsd_histogram_original(
    data: np.ndarray,
    output_path: str,
    bins: int = 200,
    year: Optional[int] = None,
    n_cells: Optional[int] = None,
    n_genes: Optional[int] = None,
    gene_rate_groups: Optional[List[Tuple[int, float]]] = None,
    custom_stats: Optional[Dict[str, float]] = None
) -> go.Figure:
    """
    Plot cell JSD histogram with original styling.
    Blue theme for JSD plots.

    Args:
        data: Array of JSD values to plot
        output_path: Path to save the plot
        bins: Number of histogram bins (default 200)
        year: Year for display
        n_cells: Number of cells for subtitle
        n_genes: Number of genes for subtitle
        gene_rate_groups: Gene rate groups for subtitle
        custom_stats: Pre-calculated statistics
    """
    return plot_histogram_original(
        data=data,
        title="Cell JSD",
        xlabel="Cell JSD",
        ylabel="Count",
        histogram_color="#1f77b4",  # Blue for JSD
        mean_line_color="red",  # Red mean line for JSD
        bins=bins,
        custom_stats=custom_stats,
        output_path=output_path,
        year=year,
        n_cells=n_cells,
        n_genes=n_genes,
        gene_rate_groups=gene_rate_groups,
        plot_type="jsd"
    )


def plot_cell_methylation_histogram_original(
    data: np.ndarray,
    output_path: str,
    bins: int = 200,
    year: Optional[int] = None,
    n_cells: Optional[int] = None,
    n_genes: Optional[int] = None,
    gene_rate_groups: Optional[List[Tuple[int, float]]] = None,
    custom_stats: Optional[Dict[str, float]] = None
) -> go.Figure:
    """
    Plot cell methylation proportion histogram with original styling.
    Red theme for methylation plots.

    Args:
        data: Array of methylation proportion values to plot
        output_path: Path to save the plot
        bins: Number of histogram bins (default 200)
        year: Year for display
        n_cells: Number of cells for subtitle
        n_genes: Number of genes for subtitle
        gene_rate_groups: Gene rate groups for subtitle
        custom_stats: Pre-calculated statistics
    """
    return plot_histogram_original(
        data=data,
        title="Cell Methylation Proportion",
        xlabel="Cell Methylation Proportion",
        ylabel="Count",
        histogram_color="#d62728",  # Red for methylation
        mean_line_color="darkblue",  # Dark blue mean line for methylation
        bins=bins,
        custom_stats=custom_stats,
        output_path=output_path,
        year=year,
        n_cells=n_cells,
        n_genes=n_genes,
        gene_rate_groups=gene_rate_groups,
        plot_type="methylation"
    )