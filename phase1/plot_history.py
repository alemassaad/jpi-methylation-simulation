#!/usr/bin/env python3
"""
Plot JSD vs time from phase1 simulation history using Plotly.
Shows mean with 5-95 percentile bands.
Adapted for the hierarchical directory structure and PetriDish simulations.
"""

import json
import argparse
import os
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import gzip
import glob

def find_simulation_file(path_pattern):
    """Find simulation file from pattern (supports wildcards)."""
    matches = glob.glob(path_pattern)
    if not matches:
        # Try in data directory
        data_pattern = os.path.join('data', path_pattern)
        matches = glob.glob(data_pattern)
    
    if not matches:
        return None
    
    # If multiple matches, pick the first .json.gz file
    for match in matches:
        if match.endswith('simulation.json.gz'):
            return match
    
    # Otherwise return first match
    return matches[0]

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
        'population_size': [],
        'jsd': {
            'mean': [],
            'median': [],
            'p5': [],
            'p25': [],
            'p75': [],
            'p95': [],
            'min': [],
            'max': []
        },
        'methylation': {
            'mean': [],
            'median': [],
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
        stats['population_size'].append(len(year_data))
        
        # JSD statistics
        stats['jsd']['mean'].append(np.mean(jsd_values))
        stats['jsd']['median'].append(np.median(jsd_values))
        stats['jsd']['p5'].append(np.percentile(jsd_values, 5))
        stats['jsd']['p25'].append(np.percentile(jsd_values, 25))
        stats['jsd']['p75'].append(np.percentile(jsd_values, 75))
        stats['jsd']['p95'].append(np.percentile(jsd_values, 95))
        stats['jsd']['min'].append(np.min(jsd_values))
        stats['jsd']['max'].append(np.max(jsd_values))
        
        # Methylation statistics
        stats['methylation']['mean'].append(np.mean(meth_values))
        stats['methylation']['median'].append(np.median(meth_values))
        stats['methylation']['p5'].append(np.percentile(meth_values, 5))
        stats['methylation']['p25'].append(np.percentile(meth_values, 25))
        stats['methylation']['p75'].append(np.percentile(meth_values, 75))
        stats['methylation']['p95'].append(np.percentile(meth_values, 95))
    
    return stats

def detect_growth_phase(stats):
    """Detect the growth phase duration by looking for population stabilization."""
    pop_sizes = stats['population_size']
    
    # Look for when population stops doubling predictably
    for i in range(1, len(pop_sizes)):
        if i > 1:  # Skip year 0 and 1
            expected = 2 ** i
            actual = pop_sizes[i]
            # If population is not exactly 2^i, we've left growth phase
            if actual != expected:
                return i - 1
    
    # Default to 13 if not detected
    return 13

def create_jsd_plot(stats, filename):
    """Create JSD-only plot with growth phase indicator and cell count."""
    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    years = stats['years']
    growth_phase = detect_growth_phase(stats)
    
    # Add vertical line at end of growth phase
    if growth_phase > 0 and growth_phase < years[-1]:
        fig.add_vline(
            x=growth_phase,
            line_dash="dash",
            line_color="gray",
            line_width=1,
            annotation_text=f"End of growth phase",
            annotation_position="top"
        )
    
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
        secondary_y=False
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
        secondary_y=False
    )
    
    # Add mean line (solid)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['jsd']['mean'],
            mode='lines',
            name='Mean JSD',
            line=dict(color='rgb(99, 110, 250)', width=2.5),
            hovertemplate='Year: %{x}<br>Mean JSD: %{y:.4f}<br>Population: %{customdata}<extra></extra>',
            customdata=stats['population_size']
        ),
        secondary_y=False
    )
    
    # Add cell count on secondary y-axis
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['population_size'],
            mode='lines',
            name='Cell Count',
            line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
            hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>'
        ),
        secondary_y=True
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"JSD Score vs Time<br><sub>{filename}</sub>",
            font=dict(size=16)
        ),
        xaxis_title="Age (years)",
        height=500,
        margin=dict(t=120),  # Add top margin for annotation
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
    
    # Set axis titles and formatting
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(title_text="JSD Score", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
    
    # Add annotation with final statistics
    final_idx = -1
    final_year = years[final_idx]
    final_pop = stats['population_size'][final_idx]
    
    annotation_text = (
        f"<b>Final Statistics (Year {final_year}):</b><br>"
        f"Population: {final_pop} cells<br>"
        f"Mean JSD: {stats['jsd']['mean'][final_idx]:.4f}<br>"
        f"25-75%: [{stats['jsd']['p25'][final_idx]:.4f}, {stats['jsd']['p75'][final_idx]:.4f}]<br>"
        f"5-95%: [{stats['jsd']['p5'][final_idx]:.4f}, {stats['jsd']['p95'][final_idx]:.4f}]"
    )
    
    if growth_phase > 0:
        annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
    
    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.5, y=1.15,  # Centered above the plot
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="rgba(0, 0, 0, 0.2)",
        borderwidth=1,
        font=dict(size=10),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    return fig

def create_methylation_plot(stats, filename):
    """Create methylation-only plot with growth phase indicator and cell count."""
    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    years = stats['years']
    growth_phase = detect_growth_phase(stats)
    
    # Add vertical line at end of growth phase
    if growth_phase > 0 and growth_phase < years[-1]:
        fig.add_vline(
            x=growth_phase,
            line_dash="dash",
            line_color="gray",
            line_width=1,
            annotation_text=f"End of growth phase",
            annotation_position="top"
        )
    
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
        ),
        secondary_y=False
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
        ),
        secondary_y=False
    )
    
    # Add mean line (solid)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['methylation']['mean'],
            mode='lines',
            name='Mean Methylation',
            line=dict(color='rgb(239, 85, 59)', width=2.5),
            hovertemplate='Year: %{x}<br>Mean Methylation: %{y:.2f}%<br>Population: %{customdata}<extra></extra>',
            customdata=stats['population_size']
        ),
        secondary_y=False
    )
    
    # Add cell count on secondary y-axis
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['population_size'],
            mode='lines',
            name='Cell Count',
            line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
            hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>'
        ),
        secondary_y=True
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"Methylation Proportion vs Time<br><sub>{filename}</sub>",
            font=dict(size=16)
        ),
        xaxis_title="Age (years)",
        height=500,
        margin=dict(t=120),  # Add top margin for annotation
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
    
    # Set axis titles and formatting
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(title_text="Methylation (%)", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
    
    # Add annotation with final statistics
    final_idx = -1
    final_year = years[final_idx]
    final_pop = stats['population_size'][final_idx]
    
    annotation_text = (
        f"<b>Final Statistics (Year {final_year}):</b><br>"
        f"Population: {final_pop} cells<br>"
        f"Mean: {stats['methylation']['mean'][final_idx]:.2f}%<br>"
        f"25-75%: [{stats['methylation']['p25'][final_idx]:.2f}%, "
        f"{stats['methylation']['p75'][final_idx]:.2f}%]<br>"
        f"5-95%: [{stats['methylation']['p5'][final_idx]:.2f}%, "
        f"{stats['methylation']['p95'][final_idx]:.2f}%]"
    )
    
    if growth_phase > 0:
        annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
    
    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.5, y=1.15,  # Centered above the plot
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="rgba(0, 0, 0, 0.2)",
        borderwidth=1,
        font=dict(size=10),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    return fig

def create_combined_plot(stats, filename):
    """Create combined JSD and methylation plot."""
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=('JSD Score vs Time', 
                      'Methylation Proportion vs Time'),
        vertical_spacing=0.12,
        row_heights=[0.5, 0.5],
        specs=[[{"secondary_y": True}],
               [{"secondary_y": True}]]
    )
    
    years = stats['years']
    growth_phase = detect_growth_phase(stats)
    
    # Add vertical lines at end of growth phase
    if growth_phase > 0 and growth_phase < years[-1]:
        fig.add_vline(
            x=growth_phase,
            line_dash="dash",
            line_color="gray",
            line_width=1,
            annotation_text=f"End of growth phase",
            annotation_position="top",
            row=1, col=1
        )
        fig.add_vline(
            x=growth_phase,
            line_dash="dash",
            line_color="gray",
            line_width=1,
            row=2, col=1
        )
    
    # ===== JSD Plot (Row 1) =====
    # Add 5-95 percentile band
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
        row=1, col=1, secondary_y=False
    )
    
    # Add 25-75 percentile band
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
        row=1, col=1, secondary_y=False
    )
    
    # Add mean JSD line
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['jsd']['mean'],
            mode='lines',
            name='Mean JSD',
            line=dict(color='rgb(99, 110, 250)', width=2.5),
            hovertemplate='Year: %{x}<br>Mean JSD: %{y:.4f}<extra></extra>'
        ),
        row=1, col=1, secondary_y=False
    )
    
    # ===== Methylation Plot (Row 2) =====
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
        row=2, col=1, secondary_y=False
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
        row=2, col=1, secondary_y=False
    )
    
    # Add mean methylation line
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['methylation']['mean'],
            mode='lines',
            name='Mean Methylation',
            line=dict(color='rgb(239, 85, 59)', width=2.5),
            hovertemplate='Year: %{x}<br>Mean Methylation: %{y:.2f}%<extra></extra>'
        ),
        row=2, col=1, secondary_y=False
    )
    
    # Add cell count traces on secondary y-axes
    # Cell count for JSD plot (row 1)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['population_size'],
            mode='lines',
            name='Cell Count',
            line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
            hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
            showlegend=True
        ),
        row=1, col=1, secondary_y=True
    )
    
    # Cell count for methylation plot (row 2) - don't show in legend to avoid duplication
    fig.add_trace(
        go.Scatter(
            x=years,
            y=stats['population_size'],
            mode='lines',
            name='Cell Count',
            line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
            hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
            showlegend=False
        ),
        row=2, col=1, secondary_y=True
    )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text=f"Simulation Results: {filename}",
            font=dict(size=16)
        ),
        height=800,
        margin=dict(t=120),  # Add top margin for annotation
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
    fig.update_xaxes(title_text="", row=1, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_xaxes(title_text="Age (years)", row=2, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(title_text="JSD Score", row=1, col=1, secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(title_text="Cell Count", row=1, col=1, secondary_y=True, showgrid=False)
    fig.update_yaxes(title_text="Methylation (%)", row=2, col=1, secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(title_text="Cell Count", row=2, col=1, secondary_y=True, showgrid=False)
    
    # Add annotation with final statistics
    final_idx = -1
    final_year = years[final_idx]
    final_pop = stats['population_size'][final_idx]
    
    annotation_text = (
        f"<b>Final Statistics (Year {final_year}):</b><br>"
        f"Population: {final_pop} cells<br>"
        f"JSD Mean: {stats['jsd']['mean'][final_idx]:.4f}<br>"
        f"JSD 25-75%: [{stats['jsd']['p25'][final_idx]:.4f}, {stats['jsd']['p75'][final_idx]:.4f}]<br>"
        f"JSD 5-95%: [{stats['jsd']['p5'][final_idx]:.4f}, {stats['jsd']['p95'][final_idx]:.4f}]<br>"
        f"Methylation Mean: {stats['methylation']['mean'][final_idx]:.2f}%<br>"
        f"Methylation 25-75%: [{stats['methylation']['p25'][final_idx]:.2f}%, "
        f"{stats['methylation']['p75'][final_idx]:.2f}%]<br>"
        f"Methylation 5-95%: [{stats['methylation']['p5'][final_idx]:.2f}%, "
        f"{stats['methylation']['p95'][final_idx]:.2f}%]"
    )
    
    if growth_phase > 0:
        annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
    
    fig.add_annotation(
        text=annotation_text,
        xref="paper", yref="paper",
        x=0.5, y=1.15,  # Centered above the plot
        showarrow=False,
        bgcolor="rgba(255, 255, 255, 0.9)",
        bordercolor="rgba(0, 0, 0, 0.2)",
        borderwidth=1,
        font=dict(size=10),
        align="center",
        xanchor="center",
        yanchor="top"
    )
    
    return fig

def main():
    parser = argparse.ArgumentParser(
        description='Plot JSD and methylation from phase1 simulation history',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument('history_file', type=str, 
                        help='Path to simulation.json.gz file (supports wildcards)')
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output directory (default: plots/ in simulation directory)')
    parser.add_argument('--jsd-only', action='store_true',
                        help='Only create JSD plot')
    parser.add_argument('--methylation-only', action='store_true',
                        help='Only create methylation plot')
    parser.add_argument('--combined', action='store_true',
                        help='Create combined plot (both JSD and methylation in one figure)')
    
    args = parser.parse_args()
    
    # Find simulation file
    if '*' in args.history_file or '?' in args.history_file:
        filepath = find_simulation_file(args.history_file)
        if not filepath:
            print(f"Error: No simulation file found matching pattern: {args.history_file}")
            return 1
    else:
        filepath = args.history_file
        if not os.path.exists(filepath):
            print(f"Error: File not found: {filepath}")
            return 1
    
    # Determine output directory
    if args.output:
        plots_dir = args.output
    else:
        # Create plots directory in the same location as the simulation
        sim_dir = os.path.dirname(filepath)
        plots_dir = os.path.join(sim_dir, 'plots')
    
    if not os.path.exists(plots_dir):
        os.makedirs(plots_dir, exist_ok=True)
        print(f"Created {plots_dir}/ directory")
    
    # Extract base name for output files
    # Handle hierarchical structure: extract parameters from path
    if 'grow' in filepath and 'sites' in filepath:
        # Extract parameters from directory name
        dir_name = os.path.basename(os.path.dirname(filepath))
        base_name = dir_name
    else:
        base_name = os.path.basename(filepath).replace('.json.gz', '').replace('.json', '')
    
    # Load and process data
    history = load_history(filepath)
    
    print(f"Processing {len(history)} time points...")
    stats = calculate_statistics(history)
    
    # Extract filename for plot titles
    filename = base_name
    
    # Create requested plots
    if args.combined:
        try:
            print("Creating combined plot...")
            fig_combined = create_combined_plot(stats, filename)
            combined_output = os.path.join(plots_dir, f"{base_name}_combined.png")
            fig_combined.write_image(combined_output, width=1200, height=800, scale=2)
            print(f"  Combined plot saved to {combined_output}")
        except Exception as e:
            print(f"Error creating combined plot: {e}")
            return 1
    else:
        # Create separate plots
        if not args.methylation_only:
            try:
                print("Creating JSD plot...")
                fig_jsd = create_jsd_plot(stats, filename)
                jsd_output = os.path.join(plots_dir, f"{base_name}_jsd.png")
                fig_jsd.write_image(jsd_output, width=1200, height=500, scale=2)
                print(f"  JSD plot saved to {jsd_output}")
            except Exception as e:
                print(f"Error creating JSD plot: {e}")
                return 1
        
        if not args.jsd_only:
            try:
                print("Creating methylation plot...")
                fig_meth = create_methylation_plot(stats, filename)
                meth_output = os.path.join(plots_dir, f"{base_name}_methylation.png")
                fig_meth.write_image(meth_output, width=1200, height=500, scale=2)
                print(f"  Methylation plot saved to {meth_output}")
            except Exception as e:
                print(f"Error creating methylation plot: {e}")
                return 1
    
    print("Done!")
    return 0

if __name__ == "__main__":
    exit(main())