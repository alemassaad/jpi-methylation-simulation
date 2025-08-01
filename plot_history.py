#!/usr/bin/env python3
"""
Plot JSD vs time from simulation history using Plotly.
Shows median with 5-95 percentile bands.
"""

import json
import argparse
import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import gzip

def load_history(filepath):
    """Load history from JSON file."""
    print(f"Loading history from {filepath}...")
    
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            return json.load(f)
    else:
        with open(filepath, 'r') as f:
            return json.load(f)

def calculate_statistics(history):
    """Calculate statistics for JSD and methylation at each time point."""
    stats = {
        'years': [],
        'jsd': {
            'median': [],
            'mean': [],
            'p5': [],
            'p25': [],
            'p75': [],
            'p95': [],
            'min': [],
            'max': []
        },
        'methylation': {
            'median': [],
            'mean': [],
            'p5': [],
            'p25': [],
            'p75': [],
            'p95': []
        }
    }
    
    # Sort years numerically
    sorted_years = sorted([int(year) for year in history.keys()])
    
    for year in sorted_years:
        year_data = history[str(year)]
        
        # Extract values
        jsd_values = np.array([cell['jsd'] for cell in year_data])
        meth_values = np.array([cell['methylation_proportion'] for cell in year_data]) * 100  # Convert to percentage
        
        stats['years'].append(year)
        
        # JSD statistics
        stats['jsd']['median'].append(np.median(jsd_values))
        stats['jsd']['mean'].append(np.mean(jsd_values))
        stats['jsd']['p5'].append(np.percentile(jsd_values, 5))
        stats['jsd']['p25'].append(np.percentile(jsd_values, 25))
        stats['jsd']['p75'].append(np.percentile(jsd_values, 75))
        stats['jsd']['p95'].append(np.percentile(jsd_values, 95))
        stats['jsd']['min'].append(np.min(jsd_values))
        stats['jsd']['max'].append(np.max(jsd_values))
        
        # Methylation statistics
        stats['methylation']['median'].append(np.median(meth_values))
        stats['methylation']['mean'].append(np.mean(meth_values))
        stats['methylation']['p5'].append(np.percentile(meth_values, 5))
        stats['methylation']['p25'].append(np.percentile(meth_values, 25))
        stats['methylation']['p75'].append(np.percentile(meth_values, 75))
        stats['methylation']['p95'].append(np.percentile(meth_values, 95))
    
    return stats

def create_plot(stats, filepath, show_methylation=True):
    """Create interactive Plotly visualization."""
    
    # Extract filename for title
    filename = os.path.basename(filepath).replace('.json', '').replace('.gz', '')
    
    # Determine subplot configuration
    if show_methylation:
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('Jensen-Shannon Divergence vs Time', 
                          'Methylation Proportion vs Time'),
            vertical_spacing=0.12,
            row_heights=[0.5, 0.5]
        )
        row_jsd, row_meth = 1, 2
    else:
        fig = go.Figure()
        row_jsd, row_meth = None, None
    
    years = stats['years']
    
    # ===== JSD Plot =====
    # Add 5-95 percentile band (lightest)
    fig.add_trace(
        go.Scatter(
            x=years + years[::-1],
            y=stats['jsd']['p95'] + stats['jsd']['p5'][::-1],
            fill='toself',
            fillcolor='rgba(99, 110, 250, 0.15)',
            line=dict(color='rgba(255,255,255,0)'),
            name='5-95 percentile',
            showlegend=True,
            hoverinfo='skip'
        ),
        row=row_jsd, col=1 if show_methylation else None
    )
    
    # Add 25-75 percentile band (medium)
    fig.add_trace(
        go.Scatter(
            x=years + years[::-1],
            y=stats['jsd']['p75'] + stats['jsd']['p25'][::-1],
            fill='toself',
            fillcolor='rgba(99, 110, 250, 0.25)',
            line=dict(color='rgba(255,255,255,0)'),
            name='25-75 percentile',
            showlegend=True,
            hoverinfo='skip'
        ),
        row=row_jsd, col=1 if show_methylation else None
    )
    
    # Add median line
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['jsd']['median'],
            mode='lines',
            name='Median JSD',
            line=dict(color='rgb(99, 110, 250)', width=2.5),
            hovertemplate='Year: %{x}<br>Median JSD: %{y:.4f}<extra></extra>'
        ),
        row=row_jsd, col=1 if show_methylation else None
    )
    
    # Add mean line (dashed)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['jsd']['mean'],
            mode='lines',
            name='Mean JSD',
            line=dict(color='rgb(99, 110, 250)', width=1.5, dash='dash'),
            hovertemplate='Year: %{x}<br>Mean JSD: %{y:.4f}<extra></extra>'
        ),
        row=row_jsd, col=1 if show_methylation else None
    )
    
    # ===== Methylation Plot =====
    if show_methylation:
        # Add 5-95 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=stats['methylation']['p95'] + stats['methylation']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=False,
                hoverinfo='skip'
            ),
            row=row_meth, col=1
        )
        
        # Add 25-75 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=stats['methylation']['p75'] + stats['methylation']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=False,
                hoverinfo='skip'
            ),
            row=row_meth, col=1
        )
        
        # Add median line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=stats['methylation']['median'],
                mode='lines',
                name='Median Methylation',
                line=dict(color='rgb(239, 85, 59)', width=2.5),
                hovertemplate='Year: %{x}<br>Median Methylation: %{y:.2f}%<extra></extra>'
            ),
            row=row_meth, col=1
        )
        
        # Add mean line (dashed)
        fig.add_trace(
            go.Scatter(
                x=years,
                y=stats['methylation']['mean'],
                mode='lines',
                name='Mean Methylation',
                line=dict(color='rgb(239, 85, 59)', width=1.5, dash='dash'),
                hovertemplate='Year: %{x}<br>Mean Methylation: %{y:.2f}%<extra></extra>'
            ),
            row=row_meth, col=1
        )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"Simulation Results: {filename}",
            font=dict(size=16)
        ),
        height=800 if show_methylation else 500,
        hovermode='x unified',
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1
        ),
        template='plotly_white'
    )
    
    # Update axes
    if show_methylation:
        fig.update_xaxes(title_text="Age (years)", row=2, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_xaxes(title_text="", row=1, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="JSD", row=1, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Methylation (%)", row=2, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    else:
        fig.update_xaxes(title_text="Age (years)", showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Jensen-Shannon Divergence", showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    
    # Add annotations with final statistics
    final_idx = -1
    final_year = years[final_idx]
    
    annotation_text = (
        f"<b>Final Statistics (Year {final_year}):</b><br>"
        f"JSD Median: {stats['jsd']['median'][final_idx]:.4f}<br>"
        f"JSD Mean: {stats['jsd']['mean'][final_idx]:.4f}<br>"
        f"JSD Range: [{stats['jsd']['p5'][final_idx]:.4f}, {stats['jsd']['p95'][final_idx]:.4f}]"
    )
    
    if show_methylation:
        annotation_text += (
            f"<br><br>Methylation Median: {stats['methylation']['median'][final_idx]:.2f}%<br>"
            f"Methylation Mean: {stats['methylation']['mean'][final_idx]:.2f}%<br>"
            f"Methylation Range: [{stats['methylation']['p5'][final_idx]:.2f}%, "
            f"{stats['methylation']['p95'][final_idx]:.2f}%]"
        )
    
    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.98, y=0.98,
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="rgba(0, 0, 0, 0.2)",
        borderwidth=1,
        font=dict(size=10),
        align="left",
        xanchor="right",
        yanchor="top"
    )
    
    return fig

def create_jsd_plot(stats, filename):
    """Create JSD-only plot."""
    fig = go.Figure()
    
    years = stats['years']
    
    # Add 5-95 percentile band (lightest)
    fig.add_trace(
        go.Scatter(
            x=years + years[::-1],
            y=stats['jsd']['p95'] + stats['jsd']['p5'][::-1],
            fill='toself',
            fillcolor='rgba(99, 110, 250, 0.15)',
            line=dict(color='rgba(255,255,255,0)'),
            name='5-95 percentile',
            showlegend=True,
            hoverinfo='skip'
        )
    )
    
    # Add 25-75 percentile band (medium)
    fig.add_trace(
        go.Scatter(
            x=years + years[::-1],
            y=stats['jsd']['p75'] + stats['jsd']['p25'][::-1],
            fill='toself',
            fillcolor='rgba(99, 110, 250, 0.25)',
            line=dict(color='rgba(255,255,255,0)'),
            name='25-75 percentile',
            showlegend=True,
            hoverinfo='skip'
        )
    )
    
    # Add mean line (solid)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['jsd']['mean'],
            mode='lines',
            name='Mean JSD',
            line=dict(color='rgb(99, 110, 250)', width=2.5),
            hovertemplate='Year: %{x}<br>Mean JSD: %{y:.4f}<extra></extra>'
        )
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"Jensen-Shannon Divergence vs Time<br><sub>{filename}</sub>",
            font=dict(size=16)
        ),
        xaxis_title="Age (years)",
        yaxis_title="Jensen-Shannon Divergence",
        height=500,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1
        ),
        template='plotly_white'
    )
    
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    
    # Add annotation with final statistics
    final_idx = -1
    final_year = years[final_idx]
    
    annotation_text = (
        f"<b>Final Statistics (Year {final_year}):</b><br>"
        f"Mean: {stats['jsd']['mean'][final_idx]:.4f}<br>"
        f"5-95%: [{stats['jsd']['p5'][final_idx]:.4f}, {stats['jsd']['p95'][final_idx]:.4f}]"
    )
    
    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.98, y=0.98,
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="rgba(0, 0, 0, 0.2)",
        borderwidth=1,
        font=dict(size=10),
        align="left",
        xanchor="right",
        yanchor="top"
    )
    
    return fig

def create_methylation_plot(stats, filename):
    """Create methylation-only plot."""
    fig = go.Figure()
    
    years = stats['years']
    
    # Add 5-95 percentile band
    fig.add_trace(
        go.Scatter(
            x=years + years[::-1],
            y=stats['methylation']['p95'] + stats['methylation']['p5'][::-1],
            fill='toself',
            fillcolor='rgba(239, 85, 59, 0.15)',
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
            y=stats['methylation']['p75'] + stats['methylation']['p25'][::-1],
            fill='toself',
            fillcolor='rgba(239, 85, 59, 0.25)',
            line=dict(color='rgba(255,255,255,0)'),
            name='25-75 percentile',
            showlegend=True,
            hoverinfo='skip'
        )
    )
    
    # Add mean line (solid)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['methylation']['mean'],
            mode='lines',
            name='Mean Methylation',
            line=dict(color='rgb(239, 85, 59)', width=2.5),
            hovertemplate='Year: %{x}<br>Mean Methylation: %{y:.2f}%<extra></extra>'
        )
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"Methylation Proportion vs Time<br><sub>{filename}</sub>",
            font=dict(size=16)
        ),
        xaxis_title="Age (years)",
        yaxis_title="Methylation (%)",
        height=500,
        showlegend=True,
        legend=dict(
            yanchor="top",
            y=0.99,
            xanchor="left",
            x=0.01,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1
        ),
        template='plotly_white'
    )
    
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    
    # Add annotation with final statistics
    final_idx = -1
    final_year = years[final_idx]
    
    annotation_text = (
        f"<b>Final Statistics (Year {final_year}):</b><br>"
        f"Mean: {stats['methylation']['mean'][final_idx]:.2f}%<br>"
        f"5-95%: [{stats['methylation']['p5'][final_idx]:.2f}%, "
        f"{stats['methylation']['p95'][final_idx]:.2f}%]"
    )
    
    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.98, y=0.98,
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="rgba(0, 0, 0, 0.2)",
        borderwidth=1,
        font=dict(size=10),
        align="left",
        xanchor="right",
        yanchor="top"
    )
    
    return fig

def main():
    parser = argparse.ArgumentParser(description='Plot JSD and methylation from simulation history using Plotly')
    parser.add_argument('history_file', type=str, 
                        help='Path to history JSON file (can be .json or .json.gz)')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output base name (default: same as input, creates _jsd.png and _methylation.png)')
    parser.add_argument('--jsd-only', action='store_true',
                        help='Only create JSD plot')
    parser.add_argument('--methylation-only', action='store_true',
                        help='Only create methylation plot')
    
    args = parser.parse_args()
    
    # Check if file exists
    if not os.path.exists(args.history_file):
        # Try in history directory
        history_path = os.path.join('history', args.history_file)
        if os.path.exists(history_path):
            args.history_file = history_path
        else:
            print(f"Error: File not found: {args.history_file}")
            return 1
    
    # Determine output base name
    if args.output is None:
        base_name = os.path.basename(args.history_file)
        base_name = base_name.replace('.json.gz', '').replace('.json', '')
    else:
        base_name = args.output
    
    # Load and process data
    history = load_history(args.history_file)
    
    print(f"Processing {len(history)} time points...")
    stats = calculate_statistics(history)
    
    # Extract filename for plot titles
    filename = os.path.basename(args.history_file).replace('.json', '').replace('.gz', '')
    
    # Create and save JSD plot
    if not args.methylation_only:
        try:
            print("Creating JSD plot...")
            fig_jsd = create_jsd_plot(stats, filename)
            jsd_output = f"{base_name}_jsd.png"
            fig_jsd.write_image(jsd_output, width=1200, height=500, scale=2)
            print(f"  JSD plot saved to {jsd_output}")
        except Exception as e:
            print(f"Error creating JSD plot: {e}")
            return 1
    
    # Create and save methylation plot
    if not args.jsd_only:
        try:
            print("Creating methylation plot...")
            fig_meth = create_methylation_plot(stats, filename)
            meth_output = f"{base_name}_methylation.png"
            fig_meth.write_image(meth_output, width=1200, height=500, scale=2)
            print(f"  Methylation plot saved to {meth_output}")
        except Exception as e:
            print(f"Error creating methylation plot: {e}")
            return 1
    
    print("Done!")
    return 0

if __name__ == "__main__":
    exit(main())