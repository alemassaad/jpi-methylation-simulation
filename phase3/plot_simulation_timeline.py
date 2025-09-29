#!/usr/bin/env python3
"""
Plot simulation timeline from extracted CSV files.
Generic plotting logic that handles all metric types (cell JSD, cell methylation, gene JSD, gene methylation).

This script provides a unified approach to timeline plotting, using a single metric-agnostic function
that can be configured for different data types through a configuration dictionary. All plots use
consistent percentile bands (5-95%, 25-75%) and show mean lines. Gene plots automatically load
cell counts from corresponding cell CSV files for secondary y-axis display.

Key features:
- Single generic function handles all 4 metric types
- Automatic cell count loading for gene metrics from cell CSVs
- Consistent percentile bands across all plots
- No median lines for gene metrics
- No combined plot generation

Usage:
    python plot_simulation_timeline.py tables/
    python plot_simulation_timeline.py tables/ --plots cell_jsd gene_jsd
    python plot_simulation_timeline.py tables/ --output-dir plots/ --title "Experiment 1"
"""

import argparse
import os
import sys
from typing import Dict, List, Optional, Tuple
import pandas as pd
import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    print("Warning: Plotly not installed. Install with: pip install plotly kaleido", file=sys.stderr)


# Metric configuration dictionary
METRIC_CONFIGS = {
    'cell_jsd': {
        'title': 'Cell JSD Score vs Time',
        'ylabel': 'Cell JSD Score',
        'color_rgb': 'rgb(99, 110, 250)',  # Blue
        'color_rgba_15': 'rgba(99, 110, 250, 0.15)',
        'color_rgba_25': 'rgba(99, 110, 250, 0.25)',
        'is_percentage': False,
        'format': '.4f',
        'has_cell_count': True,
        'show_median': False,
        'show_minmax': False,
        'metric_name': 'Cell JSD'
    },
    'cell_methylation': {
        'title': 'Cell Methylation Proportion vs Time',
        'ylabel': 'Cell Methylation Proportion (%)',
        'color_rgb': 'rgb(239, 85, 59)',  # Red
        'color_rgba_15': 'rgba(239, 85, 59, 0.15)',
        'color_rgba_25': 'rgba(239, 85, 59, 0.25)',
        'is_percentage': True,
        'format': '.2f',
        'has_cell_count': True,
        'show_median': False,
        'show_minmax': False,
        'metric_name': 'Cell Methylation'
    },
    'gene_jsd': {
        'title': 'Gene JSD Score vs Time',
        'ylabel': 'Gene JSD Score',
        'color_rgb': 'rgb(31, 119, 180)',  # Light blue
        'color_rgba_15': 'rgba(31, 119, 180, 0.1)',
        'color_rgba_25': 'rgba(31, 119, 180, 0.25)',
        'is_percentage': False,
        'format': '.4f',
        'has_cell_count': False,
        'show_median': False,  # No median for gene metrics
        'show_minmax': False,
        'metric_name': 'Gene JSD'
    },
    'gene_methylation': {
        'title': 'Gene Methylation Proportion vs Time',
        'ylabel': 'Gene Methylation Proportion (%)',
        'color_rgb': 'rgb(44, 160, 44)',  # Green
        'color_rgba_15': 'rgba(44, 160, 44, 0.1)',
        'color_rgba_25': 'rgba(44, 160, 44, 0.25)',
        'is_percentage': True,
        'format': '.2f',
        'has_cell_count': False,
        'show_median': False,  # No median for gene metrics
        'show_minmax': False,  # Use percentile bands like other plots
        'metric_name': 'Gene Methylation'
    }
}


def load_timeline_csv(csv_path: str) -> pd.DataFrame:
    """
    Load timeline CSV file.

    Args:
        csv_path: Path to CSV file

    Returns:
        DataFrame with years as index and values/genes as columns
    """
    if not os.path.exists(csv_path):
        raise FileNotFoundError(f"CSV file not found: {csv_path}")

    df = pd.read_csv(csv_path, index_col=0)
    df.index.name = 'year'
    return df


def calculate_timeline_statistics(df: pd.DataFrame) -> Dict:
    """
    Calculate statistics for timeline data.

    Args:
        df: DataFrame with years as rows and individual values as columns

    Returns:
        Dictionary with years and statistics arrays
    """
    stats = {
        'years': [],
        'mean': [],
        'median': [],
        'p5': [],
        'p25': [],
        'p75': [],
        'p95': [],
        'min': [],
        'max': [],
        'count': []
    }

    for year in df.index:
        row = df.loc[year]
        # Filter out NaN values (cells may vary during homeostasis)
        values = row.dropna().values

        if len(values) == 0:
            continue

        stats['years'].append(year)
        stats['mean'].append(np.mean(values))
        stats['median'].append(np.median(values))
        stats['p5'].append(np.percentile(values, 5))
        stats['p25'].append(np.percentile(values, 25))
        stats['p75'].append(np.percentile(values, 75))
        stats['p95'].append(np.percentile(values, 95))
        stats['min'].append(np.min(values))
        stats['max'].append(np.max(values))
        stats['count'].append(len(values))

    return stats


def load_cell_counts_from_csv(cell_csv_path: str) -> Dict[int, int]:
    """
    Load cell counts by year from a cell metric CSV.

    Args:
        cell_csv_path: Path to cell metric CSV file

    Returns:
        Dictionary mapping year to cell count
    """
    if not os.path.exists(cell_csv_path):
        return {}

    df = pd.read_csv(cell_csv_path, index_col=0)
    cell_counts = {}
    for year in df.index:
        row = df.loc[year]
        cell_counts[year] = row.dropna().shape[0]
    return cell_counts


def detect_growth_phase(cell_counts: List[int]) -> Optional[int]:
    """
    Detect the end of growth phase from population sizes.

    Args:
        cell_counts: List of cell counts per year

    Returns:
        Year where growth phase ends, or None if not detected
    """
    if len(cell_counts) < 2:
        return None

    # Look for transition from exponential growth (2^n) to homeostasis
    for i in range(1, min(len(cell_counts), 20)):  # Check first 20 years max
        if i > 1:
            expected = 2 ** i
            actual = cell_counts[i]
            # If actual deviates significantly from expected exponential
            if actual < expected * 0.9 or actual > expected * 1.1:
                return i - 1

    # If no clear transition found, check for plateau
    for i in range(2, len(cell_counts)):
        if cell_counts[i] == cell_counts[i-1]:
            return i - 1

    return None


def plot_timeline_generic(
    csv_path: str,
    metric_type: str,
    output_path: str = None,
    cell_count_csv: str = None,
    title: str = None,
    show_cell_count: bool = True,
    width: int = 1200,
    height: int = 500
) -> go.Figure:
    """
    Generic timeline plotting function for any metric type.

    Args:
        csv_path: Path to metric CSV file
        metric_type: Type of metric ('cell_jsd', 'cell_methylation', 'gene_jsd', 'gene_methylation')
        output_path: Optional path to save plot
        cell_count_csv: Optional path to cell CSV for loading cell counts (for gene metrics)
        title: Optional custom title
        show_cell_count: Whether to show cell count on secondary y-axis
        width: Plot width in pixels
        height: Plot height in pixels

    Returns:
        Plotly figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly is required for plotting")

    if metric_type not in METRIC_CONFIGS:
        raise ValueError(f"Unknown metric type: {metric_type}")

    # Get configuration for this metric
    config = METRIC_CONFIGS[metric_type]

    # Load data
    df = load_timeline_csv(csv_path)

    # Convert to percentage if needed
    if config['is_percentage']:
        df = df * 100

    # Calculate statistics
    stats = calculate_timeline_statistics(df)
    years = stats['years']

    # Load or extract cell counts
    cell_counts = None
    growth_phase = None

    if show_cell_count:
        if config['has_cell_count']:
            # Cell metrics have counts in their own stats
            cell_counts = stats['count']
        elif cell_count_csv:
            # Gene metrics need to load from cell CSV
            cell_count_dict = load_cell_counts_from_csv(cell_count_csv)
            cell_counts = [cell_count_dict.get(year, 0) for year in years]

        # Detect growth phase if we have cell counts
        if cell_counts:
            growth_phase = detect_growth_phase(cell_counts)

    # Create figure with optional secondary y-axis
    if show_cell_count and cell_counts:
        fig = make_subplots(specs=[[{"secondary_y": True}]])
    else:
        fig = go.Figure()

    # Add vertical line for growth phase (only if we have cell counts)
    if growth_phase and growth_phase > 0 and growth_phase < years[-1]:
        fig.add_vline(
            x=growth_phase,
            line_dash="dash",
            line_color="gray",
            line_width=1,
            annotation_text="End of growth phase",
            annotation_position="top"
        )

    # Add percentile bands or min-max range
    if config.get('show_minmax'):
        # Min-max range (for gene methylation)
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=stats['max'] + stats['min'][::-1],
                fill='toself',
                fillcolor=config['color_rgba_15'],
                line=dict(width=0),
                name='Min-Max range',
                showlegend=True,
                hoverinfo='skip'
            )
        )
    else:
        # Standard percentile bands
        # Add 5-95 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=stats['p95'] + stats['p5'][::-1],
                fill='toself',
                fillcolor=config['color_rgba_15'],
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            )
        )

        # Add 25-75 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=stats['p75'] + stats['p25'][::-1],
                fill='toself',
                fillcolor=config['color_rgba_25'],
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            )
        )

    # Add mean line
    unit = '%' if config['is_percentage'] else ''
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['mean'],
            mode='lines',
            name=f'Mean {config["metric_name"]}',
            line=dict(color=config['color_rgb'], width=2.5),
            hovertemplate=f'Year: %{{x}}<br>Mean: %{{y:{config["format"]}}}{unit}<extra></extra>'
        )
    )

    # Add median line if configured
    if config.get('show_median'):
        fig.add_trace(
            go.Scatter(
                x=years,
                y=stats['median'],
                mode='lines',
                name=f'Median {config["metric_name"]}',
                line=dict(color=config['color_rgb'], width=2, dash='dash'),
                hovertemplate=f'Year: %{{x}}<br>Median: %{{y:{config["format"]}}}{unit}<extra></extra>'
            )
        )

    # Add cell count on secondary y-axis if available
    if show_cell_count and cell_counts:
        fig.add_trace(
            go.Scatter(
                x=years,
                y=cell_counts,
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
                yaxis='y2' if isinstance(fig, go.FigureWidget) else None
            ),
            secondary_y=True if show_cell_count and cell_counts else False
        )

    # Update layout
    if title is None:
        title = config['title']

    fig.update_layout(
        title=dict(text=title, font=dict(size=16)),
        xaxis_title="Age (years)",
        height=height,
        margin=dict(t=120),
        hovermode='x unified',
        showlegend=True,
        legend=dict(
            yanchor="top", y=0.99, xanchor="left", x=0.01,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1
        ),
        template='plotly_white'
    )

    # Set axis titles
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')

    if show_cell_count and cell_counts:
        fig.update_yaxes(title_text=config['ylabel'], secondary_y=False,
                        showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
    else:
        fig.update_yaxes(title_text=config['ylabel'], showgrid=True, gridcolor='rgba(0,0,0,0.1)')

    # Add annotation with statistics
    final_idx = -1
    final_year = years[final_idx]

    # Build annotation text based on metric type
    annotation_lines = [f"<b>Final Statistics (Year {final_year}):</b>"]

    # Add population/gene count
    if config['has_cell_count'] and cell_counts:
        annotation_lines.append(f"Population: {cell_counts[final_idx]} cells")
    elif metric_type.startswith('gene'):
        n_values = stats['count'][0] if stats['count'] else 20
        annotation_lines.append(f"Number of Genes: {n_values}")

    # Add mean
    annotation_lines.append(
        f"{config['metric_name']} Mean: {stats['mean'][final_idx]:{config['format']}}{unit}"
    )

    # Add median if shown
    if config.get('show_median'):
        annotation_lines.append(
            f"{config['metric_name']} Median: {stats['median'][final_idx]:{config['format']}}{unit}"
        )

    # Add ranges
    if config.get('show_minmax'):
        annotation_lines.append(
            f"{config['metric_name']} Range: [{stats['min'][final_idx]:{config['format']}}{unit}, "
            f"{stats['max'][final_idx]:{config['format']}}{unit}]"
        )
    else:
        annotation_lines.append(
            f"{config['metric_name']} 25-75%: [{stats['p25'][final_idx]:{config['format']}}{unit}, "
            f"{stats['p75'][final_idx]:{config['format']}}{unit}]"
        )
        annotation_lines.append(
            f"{config['metric_name']} 5-95%: [{stats['p5'][final_idx]:{config['format']}}{unit}, "
            f"{stats['p95'][final_idx]:{config['format']}}{unit}]"
        )

    # Add growth phase if detected
    if growth_phase:
        annotation_lines.append(f"Growth phase: Years 0-{growth_phase}")

    annotation_text = '<br>'.join(annotation_lines)

    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.5, y=1.15,
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="rgba(0, 0, 0, 0.2)",
        borderwidth=1,
        font=dict(size=10),
        align="center",
        xanchor="center",
        yanchor="top"
    )

    if output_path:
        fig.write_image(output_path, width=width, height=height, scale=2)

    return fig



def plot_all_timelines(tables_dir: str, output_dir: str = None,
                      title_prefix: str = None, plots_to_generate: List[str] = None) -> None:
    """
    Plot all timeline types from CSV files using the generic function.

    Args:
        tables_dir: Directory containing CSV files
        output_dir: Optional output directory for plots
        title_prefix: Optional prefix for plot titles
        plots_to_generate: List of plots to generate (default: all)
    """
    # Check directory exists
    if not os.path.exists(tables_dir):
        raise FileNotFoundError(f"Tables directory not found: {tables_dir}")

    # Default output directory
    if output_dir is None:
        output_dir = os.path.join(tables_dir, '..', 'plots')
    os.makedirs(output_dir, exist_ok=True)

    # Default plots
    if plots_to_generate is None:
        plots_to_generate = ['cell_jsd', 'cell_methylation', 'gene_jsd', 'gene_methylation']

    print(f"\nGenerating timeline plots from: {tables_dir}")
    print(f"Output directory: {output_dir}")

    # Define metrics and their CSV files
    metrics = [
        ('cell_jsd', 'simulation_cell_jsd.csv', None),
        ('cell_methylation', 'simulation_cell_methylation.csv', None),
        ('gene_jsd', 'simulation_gene_jsd.csv', 'simulation_cell_jsd.csv'),  # Use cell CSV for counts
        ('gene_methylation', 'simulation_gene_methylation.csv', 'simulation_cell_methylation.csv')
    ]

    # Generate individual timeline plots
    for metric_type, csv_file, cell_count_csv in metrics:
        if metric_type not in plots_to_generate:
            continue

        csv_path = os.path.join(tables_dir, csv_file)
        if not os.path.exists(csv_path):
            continue

        print(f"  - Creating {metric_type} timeline...")

        cell_count_path = os.path.join(tables_dir, cell_count_csv) if cell_count_csv else None
        output_path = os.path.join(output_dir, f'{csv_file.replace(".csv", "_timeline.png")}')

        custom_title = None
        if title_prefix:
            config = METRIC_CONFIGS[metric_type]
            custom_title = f"{title_prefix} - {config['title']}"

        try:
            plot_timeline_generic(
                csv_path=csv_path,
                metric_type=metric_type,
                output_path=output_path,
                cell_count_csv=cell_count_path,
                title=custom_title
            )
        except Exception as e:
            print(f"    Warning: Failed to create {metric_type} plot: {e}")

    print("\nTimeline plots generated successfully!")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description='Plot simulation timeline from extracted CSV files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot all timelines from CSV directory
  python plot_simulation_timeline.py results/tables/

  # Plot specific metrics only
  python plot_simulation_timeline.py results/tables/ --plots cell_jsd gene_jsd

  # Save to specific directory with custom title
  python plot_simulation_timeline.py results/tables/ --output-dir timeline_plots/ --title "Experiment 1"

  # Plot without cell count on secondary axis
  python plot_simulation_timeline.py results/tables/ --skip-cell-count
        """
    )

    parser.add_argument('tables_dir',
                       help='Directory containing CSV files from extract_simulation_timeline.py')

    parser.add_argument('--output-dir', '-o',
                       help='Output directory for plots (default: tables_dir/../plots)')

    parser.add_argument('--plots', '-p', nargs='+',
                       choices=['cell_jsd', 'cell_methylation', 'gene_jsd', 'gene_methylation'],
                       help='Specific plots to generate (default: all)')

    parser.add_argument('--title', '-t',
                       help='Title prefix for plots')

    parser.add_argument('--skip-cell-count', action='store_true',
                       help='Do not show cell count on secondary y-axis')

    parser.add_argument('--width', type=int, default=1200,
                       help='Plot width in pixels (default: 1200)')

    parser.add_argument('--height', type=int, default=500,
                       help='Plot height in pixels for single plots (default: 500)')

    args = parser.parse_args()

    # Check for plotly
    if not HAS_PLOTLY:
        print("Error: Plotly is required for plotting.", file=sys.stderr)
        print("Install with: pip install plotly kaleido", file=sys.stderr)
        return 1

    try:
        # Generate plots
        plot_all_timelines(
            tables_dir=args.tables_dir,
            output_dir=args.output_dir,
            title_prefix=args.title,
            plots_to_generate=args.plots
        )
        return 0
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error generating plots: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())