#!/usr/bin/env python3
"""
Generic comparison plotting function for batch comparisons.
Creates scatter plots with jitter, mean lines, and percentile bands.
"""

import os
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from typing import Dict, Optional, Tuple


def calculate_adaptive_dtick(data_range: float) -> Optional[float]:
    """
    Calculate appropriate tick spacing based on data range.

    Args:
        data_range: The range of data values (max - min)

    Returns:
        Appropriate dtick value, or None to use Plotly auto mode
    """
    if data_range < 0.01:      # Very close values (e.g., 0.001 range)
        return 0.001
    elif data_range < 0.02:    # Very tight clustering
        return 0.002
    elif data_range < 0.05:    # Close values (e.g., methylation 0.15-0.18)
        return 0.005
    elif data_range < 0.1:     # Moderate range
        return 0.01
    elif data_range < 0.2:     # Medium range (common for methylation)
        return 0.02            # NEW: Better for methylation ranges like 0.1-0.25
    elif data_range < 0.3:     # Medium-wide range
        return 0.025           # NEW: Intermediate step
    elif data_range < 0.5:     # Wider range
        return 0.05
    elif data_range < 1.0:     # Full 0-1 range
        return 0.1
    else:                       # Large values
        return None  # Let Plotly auto-determine


def plot_comparison_generic(
    csv_path: str,
    value_column: str,
    title: str,
    ylabel: str,
    output_path: str,
    batch_colors: Optional[Dict[str, str]] = None,
    figsize: Tuple[int, int] = (1200, 600),
    scale: int = 2,
    seed: int = 42,
    verbose: bool = True,
    y_tick_spacing: Optional[float] = None,  # Manual override for tick spacing
    auto_tick_density: bool = True,          # Enable/disable adaptive spacing
    aggregation: str = 'mean',               # Aggregation method: 'mean' or 'std'
    y_range_padding: float = 0.1,            # Padding factor for y-axis range (0.1 = 10%)
    y_range_min: Optional[float] = None,     # Manual override for y-axis minimum
    y_range_max: Optional[float] = None      # Manual override for y-axis maximum
) -> None:
    """
    Create a scatter plot comparing values across batches.

    Args:
        csv_path: Path to CSV file with batch comparison data
        value_column: Column name to plot
        title: Plot title
        ylabel: Y-axis label
        output_path: Path to save the plot
        batch_colors: Optional custom colors for batches
        figsize: Figure size in pixels (width, height)
        scale: Scale factor for high DPI
        seed: Random seed for consistent jitter
        verbose: Whether to print progress
        y_tick_spacing: Manual override for y-axis tick spacing (e.g., 0.005)
        auto_tick_density: Whether to automatically adjust tick density based on data range
        aggregation: Aggregation method when gene_index present: 'mean' or 'std' (default='mean')
        y_range_padding: Padding factor for y-axis range (0.1 = 10% padding on each side)
        y_range_min: Manual override for y-axis minimum value
        y_range_max: Manual override for y-axis maximum value
    """
    if verbose:
        print(f"\n  Creating {title}...")

    # Default colors matching phase3 style
    if batch_colors is None:
        batch_colors = {
            'mutant': '#1f77b4',    # Blue
            'control1': '#ff7f0e',   # Orange
            'control2': '#2ca02c'    # Green
        }

    # Read CSV or accept DataFrame directly
    if isinstance(csv_path, pd.DataFrame):
        # Already a DataFrame, use directly
        df = csv_path.copy()
    elif isinstance(csv_path, str):
        # File path, read the CSV
        df = pd.read_csv(csv_path)
    else:
        raise ValueError(f"csv_path must be either a string file path or a pandas DataFrame, got {type(csv_path)}")

    # Check if this is the gene format (with gene_index column)
    if 'gene_index' in df.columns:
        # Gene format: check if we need to aggregate
        # If there are multiple gene_index values, aggregate by individual
        if df['gene_index'].nunique() > 1:
            # Multiple genes present: aggregate by individual
            # Check that value column exists
            if value_column not in df.columns:
                raise ValueError(f"Column '{value_column}' not found in CSV. Available columns: {df.columns.tolist()}")

            # Aggregate by individual using specified method
            if aggregation == 'std':
                df = df.groupby(['batch', 'individual_id'])[value_column].std().reset_index()
                if verbose:
                    print(f"    Aggregated gene-level data (20 genes per individual → std dev)")
            elif aggregation == 'mean':
                df = df.groupby(['batch', 'individual_id'])[value_column].mean().reset_index()
                if verbose:
                    print(f"    Aggregated gene-level data (20 genes per individual → mean)")
            else:
                raise ValueError(f"Unsupported aggregation method: {aggregation}. Use 'mean' or 'std'")
        else:
            # Single gene: use as-is (per-gene comparison mode)
            # Check that value column exists
            if value_column not in df.columns:
                raise ValueError(f"Column '{value_column}' not found in CSV. Available columns: {df.columns.tolist()}")
            # No aggregation needed, data is already per-individual for this specific gene
    else:
        # Cell format: use as-is
        # Check that value column exists
        if value_column not in df.columns:
            raise ValueError(f"Column '{value_column}' not found in CSV. Available columns: {df.columns.tolist()}")

    # Set random seed for consistent jitter
    np.random.seed(seed)

    # Create figure
    fig = go.Figure()

    # Process each batch
    batch_positions = {'mutant': 0, 'control1': 1, 'control2': 2}
    all_stats = {}
    all_values = []  # Collect all data values for range calculation

    for batch_name, x_pos in batch_positions.items():
        # Get data for this batch
        batch_data = df[df['batch'] == batch_name][value_column].values

        if len(batch_data) == 0:
            if verbose:
                print(f"    Warning: No data for {batch_name}")
            continue

        # Collect values for range calculation
        all_values.extend(batch_data)

        # Add points with jitter
        x_jittered = np.ones(len(batch_data)) * x_pos + np.random.normal(0, 0.02, len(batch_data))

        fig.add_trace(go.Scatter(
            x=x_jittered,
            y=batch_data,
            mode='markers',
            name=batch_name.capitalize(),
            marker=dict(
                color=batch_colors[batch_name],
                size=10,
                opacity=0.6
            )
        ))

        # Calculate statistics
        mean_val = np.mean(batch_data)
        median_val = np.median(batch_data)
        std_val = np.std(batch_data)
        cv_val = std_val / mean_val * 100 if mean_val != 0 else 0
        mad_val = np.median(np.abs(batch_data - median_val))
        q5, q25, q75, q95 = np.percentile(batch_data, [5, 25, 75, 95])

        # Store statistics for annotations
        all_stats[batch_name] = {
            'mean': mean_val,
            'median': median_val,
            'std': std_val,
            'cv': cv_val,
            'mad': mad_val,
            'q5': q5,
            'q25': q25,
            'q75': q75,
            'q95': q95,
            'n': len(batch_data)
        }

        # Add percentile lines
        # 5-95 range (very light, dotted)
        fig.add_shape(
            type="line",
            x0=x_pos - 0.15, x1=x_pos + 0.15,
            y0=q5, y1=q5,
            line=dict(color=batch_colors[batch_name], width=1, dash="dot"),
            opacity=0.3
        )
        fig.add_shape(
            type="line",
            x0=x_pos - 0.15, x1=x_pos + 0.15,
            y0=q95, y1=q95,
            line=dict(color=batch_colors[batch_name], width=1, dash="dot"),
            opacity=0.3
        )

        # 25-75 range (slightly darker, dashed)
        fig.add_shape(
            type="line",
            x0=x_pos - 0.15, x1=x_pos + 0.15,
            y0=q25, y1=q25,
            line=dict(color=batch_colors[batch_name], width=1.5, dash="dash"),
            opacity=0.5
        )
        fig.add_shape(
            type="line",
            x0=x_pos - 0.15, x1=x_pos + 0.15,
            y0=q75, y1=q75,
            line=dict(color=batch_colors[batch_name], width=1.5, dash="dash"),
            opacity=0.5
        )

        # Add mean line (solid, prominent) - goes last to be on top
        fig.add_shape(
            type="line",
            x0=x_pos - 0.2, x1=x_pos + 0.2,
            y0=mean_val, y1=mean_val,
            line=dict(color=batch_colors[batch_name], width=3)
        )

    # Add statistics annotations below the plot in three columns
    x_positions = {'mutant': 0.2, 'control1': 0.5, 'control2': 0.8}  # Three columns
    for batch_name, stats in all_stats.items():
        annotation_text = (
            f"<b>{batch_name.capitalize()}</b><br>"
            f"n={stats['n']}<br>"
            f"Mean: {stats['mean']:.4f}<br>"
            f"Median: {stats['median']:.4f}<br>"
            f"SD: {stats['std']:.4f}<br>"
            f"CV: {stats['cv']:.1f}%<br>"
            f"MAD: {stats['mad']:.4f}<br>"
            f"5%: {stats['q5']:.4f}<br>"
            f"95%: {stats['q95']:.4f}"
        )

        # Position annotations below the plot
        fig.add_annotation(
            x=x_positions[batch_name],
            y=-0.15,  # Below the plot
            xref="paper",  # Use paper coordinates for consistent positioning
            yref="paper",
            text=annotation_text,
            showarrow=False,
            font=dict(size=10, color=batch_colors[batch_name]),
            align="left",
            xanchor="center",
            yanchor="top"  # Anchor from top of text
        )

    # Calculate adaptive tick spacing for y-axis
    if all_values:
        y_min = np.min(all_values)
        y_max = np.max(all_values)
        data_range = y_max - y_min

        # Calculate padded y-axis range
        if y_range_padding > 0:
            padding = data_range * y_range_padding
            y_range = [y_min - padding, y_max + padding]
        else:
            y_range = [y_min, y_max]

        # Apply manual overrides if provided
        if y_range_min is not None:
            y_range[0] = y_range_min
        if y_range_max is not None:
            y_range[1] = y_range_max

        if verbose and y_range_padding > 0:
            print(f"    Y-axis range: [{y_range[0]:.4f}, {y_range[1]:.4f}] (padding={y_range_padding*100:.0f}%)")

        # Determine tick spacing
        if y_tick_spacing is not None:
            # Use manual override
            dtick = y_tick_spacing
            if verbose:
                print(f"    Using manual tick spacing: {dtick}")
        elif auto_tick_density:
            # Use adaptive calculation
            dtick = calculate_adaptive_dtick(data_range)
            if verbose and dtick:
                print(f"    Data range: {data_range:.4f}, using tick spacing: {dtick}")
        else:
            # No adaptive spacing
            dtick = None

        # Build y-axis configuration
        yaxis_config = dict(
            title=ylabel,
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False,
            range=y_range  # Add the calculated range
        )

        if dtick is not None:
            # Use linear mode with calculated spacing
            yaxis_config.update({
                'tickmode': 'linear',
                'dtick': dtick,
                'tick0': np.floor(y_min / dtick) * dtick  # Start at nice round number
            })
        else:
            # For large ranges or disabled adaptive mode, use auto mode with more ticks
            yaxis_config.update({
                'tickmode': 'auto',
                'nticks': 20  # Request more ticks
            })
    else:
        # Fallback if no data
        yaxis_config = dict(
            title=ylabel,
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False
        )

    # Update layout
    fig.update_layout(
        title=dict(
            text=title,
            font=dict(size=16, family='Arial, sans-serif')
        ),
        xaxis=dict(
            title='',  # No x-axis title
            tickmode='array',
            tickvals=[0, 1, 2],
            ticktext=['Mutant', 'Control1', 'Control2'],
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False,
            range=[-0.5, 2.5]
        ),
        yaxis=yaxis_config,
        plot_bgcolor='white',
        paper_bgcolor='white',
        width=figsize[0],
        height=figsize[1],
        showlegend=True,
        legend=dict(
            x=1.02,
            y=1,
            xanchor='left',
            yanchor='top'
        ),
        margin=dict(l=80, r=150, t=100, b=250),  # Increased bottom margin for statistics
        template='plotly_white'
    )

    # Create output directory if needed
    output_dir = os.path.dirname(output_path)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Save plot
    fig.write_image(output_path, scale=scale)
    if verbose:
        print(f"    Saved to {output_path}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Create batch comparison plot")
    parser.add_argument('--csv', required=True, help='Input CSV file path')
    parser.add_argument('--column', required=True, help='Column to plot')
    parser.add_argument('--title', required=True, help='Plot title')
    parser.add_argument('--ylabel', required=True, help='Y-axis label')
    parser.add_argument('--output', required=True, help='Output plot path')
    parser.add_argument('--quiet', action='store_true', help='Suppress progress output')

    args = parser.parse_args()

    plot_comparison_generic(
        csv_path=args.csv,
        value_column=args.column,
        title=args.title,
        ylabel=args.ylabel,
        output_path=args.output,
        verbose=not args.quiet
    )