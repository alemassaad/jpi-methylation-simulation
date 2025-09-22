"""
Gene-level JSD plotting methods.
Extended plotting functionality for gene-level JSD analysis.
Moved from phase1/cell.py to consolidate visualization in phase3.
"""
import os
import numpy as np
import colorsys
from typing import Optional

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

# Import phase1 modules
import sys
phase1_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1')
sys.path.insert(0, phase1_path)
from cell import PetriDish, Cell


def plot_gene_jsd_heatmap(petri: PetriDish, title: str = None, output_path: str = None,
                          width: int = 1200, height: int = 800):
    """
    Create heatmap showing each gene's JSD evolution over time.
    
    Args:
        petri: PetriDish object with gene_jsd_history
        title: Optional title for the plot
        output_path: Optional path to save the plot
        width: Plot width in pixels
        height: Plot height in pixels
        
    Returns:
        Plotly figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
    
    if not hasattr(petri, 'gene_jsd_history') or not petri.gene_jsd_history:
        raise ValueError("No gene JSD history found. Enable gene JSD tracking with --track-gene-jsd flag.")
    
    # Extract data from gene_jsd_history
    years = sorted([int(y) for y in petri.gene_jsd_history.keys()])
    n_genes = len(petri.gene_jsd_history[str(years[0])])
    
    # Create 2D array: rows = genes, columns = years
    heatmap_data = []
    for gene_idx in range(n_genes):
        gene_row = []
        for year in years:
            gene_jsd_value = petri.gene_jsd_history[str(year)][gene_idx]
            gene_row.append(gene_jsd_value)
        heatmap_data.append(gene_row)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=heatmap_data,
        x=years,
        y=list(range(n_genes)),
        colorscale='Blues',
        colorbar=dict(title="Gene JSD Score"),
        hovertemplate='Year: %{x}<br>Gene: %{y}<br>JSD: %{z:.4f}<extra></extra>'
    ))
    
    # Add horizontal lines every 50 genes if using gene-rate-groups
    if hasattr(petri.cells[0], 'gene_rate_groups') and petri.cells[0].gene_rate_groups:
        gene_group_boundaries = []
        cumulative = 0
        for n_genes_in_group, _ in petri.cells[0].gene_rate_groups:
            cumulative += n_genes_in_group
            if cumulative < n_genes:  # Don't add line after last group
                gene_group_boundaries.append(cumulative - 0.5)
        
        # Add horizontal lines
        for boundary in gene_group_boundaries:
            fig.add_hline(
                y=boundary,
                line_dash="solid",
                line_color="white",
                line_width=2,
                opacity=0.8
            )
    
    # Calculate max JSD for annotation
    max_jsd = max(max(row) for row in heatmap_data)
    
    # Add annotation box
    stats_text = (f"<b>Heatmap Info</b><br>"
                 f"Genes: {n_genes}<br>"
                 f"Years: {len(years)}<br>"
                 f"Max JSD: {max_jsd:.3f}")
    
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
    if title is None:
        title = "Gene JSD Evolution Heatmap"
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=16)),
        xaxis_title="Age (years)",
        yaxis_title="Gene Index",
        width=width,
        height=height,
        template='plotly_white',
        showlegend=False
    )
    
    # Update axes with grid
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    
    # Save if path provided
    if output_path:
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        fig.write_image(output_path, width=width, height=height, scale=2)
        print(f"Gene JSD heatmap saved to {output_path}")
    else:
        fig.show()
    
    return fig


def plot_gene_jsd_by_rate_group(petri: PetriDish, title: str = None, output_path: str = None,
                                width: int = 1200, height: int = 600):
    """
    Line plot showing mean JSD for each gene rate group over time.
    
    Args:
        petri: PetriDish object with gene_jsd_history
        title: Optional title for the plot
        output_path: Optional path to save the plot
        width: Plot width in pixels
        height: Plot height in pixels
        
    Returns:
        Plotly figure object
    """
    if not HAS_PLOTLY:
        raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
    
    if not hasattr(petri, 'gene_jsd_history') or not petri.gene_jsd_history:
        raise ValueError("No gene JSD history found. Enable gene JSD tracking with --track-gene-jsd flag.")
    
    # Check if we have gene rate groups
    if not hasattr(petri.cells[0], 'gene_rate_groups') or not petri.cells[0].gene_rate_groups:
        raise ValueError("No gene rate groups found. This plot requires --gene-rate-groups parameter.")
    
    # Extract gene rate groups information
    gene_rate_groups = petri.cells[0].gene_rate_groups
    
    # Extract data from gene_jsd_history
    years = sorted([int(y) for y in petri.gene_jsd_history.keys()])
    
    # Calculate mean JSD for each rate group at each time point
    fig = go.Figure()
    
    # Define colors for each group (blue to red gradient)
    n_groups = len(gene_rate_groups)
    if n_groups == 4:
        colors = ['#0066cc', '#3399ff', '#ff6600', '#cc0000']
    else:
        # Generate color gradient
        colors = []
        for i in range(n_groups):
            hue = 240 - (240 * i / max(1, n_groups - 1)) / 360  # Blue to red
            rgb = colorsys.hsv_to_rgb(hue, 0.8, 0.9)
            colors.append(f'rgb({int(rgb[0]*255)}, {int(rgb[1]*255)}, {int(rgb[2]*255)})')
    
    # Calculate statistics for each group over time
    group_stats = []
    gene_start_idx = 0
    
    for group_idx, (n_genes, rate) in enumerate(gene_rate_groups):
        gene_end_idx = gene_start_idx + n_genes
        group_means = []
        group_p25 = []
        group_p75 = []
        
        for year in years:
            gene_jsds_for_group = petri.gene_jsd_history[str(year)][gene_start_idx:gene_end_idx]
            group_means.append(np.mean(gene_jsds_for_group))
            group_p25.append(np.percentile(gene_jsds_for_group, 25))
            group_p75.append(np.percentile(gene_jsds_for_group, 75))
        
        group_stats.append({
            'mean': group_means,
            'p25': group_p25,
            'p75': group_p75,
            'rate': rate,
            'n_genes': n_genes
        })
        gene_start_idx = gene_end_idx
    
    # Add confidence bands (25-75 percentile) for each group
    for group_idx, stats in enumerate(group_stats):
        # Add confidence band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=stats['p75'] + stats['p25'][::-1],
                fill='toself',
                fillcolor=colors[group_idx].replace('rgb', 'rgba').replace(')', ', 0.2)'),
                line=dict(color='rgba(255,255,255,0)'),
                name=f'Group {group_idx+1} (25-75%)',
                showlegend=False,
                hoverinfo='skip'
            )
        )
    
    # Add mean lines for each group
    for group_idx, stats in enumerate(group_stats):
        rate_pct = stats['rate'] * 100
        fig.add_trace(
            go.Scatter(
                x=years,
                y=stats['mean'],
                mode='lines',
                name=f'Group {group_idx+1}: {rate_pct:.1f}% rate',
                line=dict(color=colors[group_idx], width=2),
                hovertemplate=f'Year: %{{x}}<br>Mean Gene JSD: %{{y:.4f}}<br>Group {group_idx+1} ({stats["n_genes"]} genes @ {rate_pct:.1f}%)<extra></extra>'
            )
        )
    
    # Add annotation box with group information
    group_info_lines = ["<b>Gene Groups:</b>"]
    for i, (n_genes, rate) in enumerate(gene_rate_groups):
        rate_pct = rate * 100
        group_info_lines.append(f"Group {i+1}: {n_genes} genes @ {rate_pct:.1f}%")
    
    group_info_text = "<br>".join(group_info_lines)
    
    fig.add_annotation(
        text=group_info_text,
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
    if title is None:
        title = "Gene JSD by Methylation Rate Group"
    
    fig.update_layout(
        title=dict(text=title, font=dict(size=16)),
        xaxis_title="Age (years)",
        yaxis_title="Mean Gene JSD Score",
        width=width,
        height=height,
        template='plotly_white',
        showlegend=True,
        legend=dict(
            yanchor="top", y=0.99, xanchor="left", x=0.01,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1
        ),
        hovermode='x unified'
    )
    
    # Update axes with grid
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    
    # Save if path provided
    if output_path:
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        fig.write_image(output_path, width=width, height=height, scale=2)
        print(f"Gene rate group comparison saved to {output_path}")
    else:
        fig.show()
    
    return fig


def plot_gene_jsd_trajectory(petri: PetriDish, title: str = None, output_path: str = None,
                            width: int = 1200, height: int = 600):
    """
    Plot gene-level JSD trajectory with percentile bands matching cell-level plot style.
    
    For each timepoint in cell_history, calculates gene-level JSD values
    (one per gene) and plots mean, median, and percentile bands.
    
    Args:
        petri: PetriDish object with cell_history
        title: Optional title for the plot
        output_path: Optional path to save the plot
        width: Plot width in pixels
        height: Plot height in pixels
    """
    if not HAS_PLOTLY:
        print("Warning: Plotly not installed. Skipping gene JSD trajectory plot.")
        return
    
    # Check if cell history exists
    if not hasattr(petri, 'cell_history') or not petri.cell_history:
        print("No cell history available for gene JSD trajectory.")
        return
    
    print("Calculating gene-level JSD trajectory...")
    
    # Initialize statistics dictionary (matching cell-level structure)
    gene_stats = {
        'years': [],
        'gene_jsd': {
            'mean': [], 'median': [], 'p5': [], 'p25': [],
            'p75': [], 'p95': [], 'min': [], 'max': []
        },
        'population_size': []  # Track cell count for secondary axis
    }
    
    # Extract years and calculate gene JSDs for each year
    years = sorted([int(y) for y in petri.cell_history.keys()])
    
    for year in years:
        # Get cells for this year
        year_data = petri.cell_history[str(year)]
        
        # Handle both formats: list of cells or dict with 'cells' key
        if isinstance(year_data, dict) and 'cells' in year_data:
            cells = year_data['cells']
        else:
            cells = year_data
        
        # Create temporary PetriDish to calculate gene JSDs
        temp_petri = PetriDish(
            gene_rate_groups=petri.gene_rate_groups,
            n=petri.n,  # Pass n to match actual cell sizes
            gene_size=petri.gene_size,  # Pass gene_size to match configuration
            growth_phase=1  # Not used for calculation
        )
        
        # Convert dictionary cells to Cell objects if needed
        if cells and len(cells) > 0 and isinstance(cells[0], dict):
            from_dict_cells = []
            for cell_dict in cells:
                cell = Cell.from_dict(cell_dict, gene_rate_groups=petri.gene_rate_groups)
                from_dict_cells.append(cell)
            temp_petri.cells = from_dict_cells
        else:
            temp_petri.cells = cells
        
        # Calculate gene-level JSD values (one per gene)
        gene_jsds = temp_petri.calculate_gene_jsd()
        
        # Store year and population size
        gene_stats['years'].append(year)
        gene_stats['population_size'].append(len(cells))
        
        # Calculate all percentiles from the 20 gene values
        gene_stats['gene_jsd']['mean'].append(np.mean(gene_jsds))
        gene_stats['gene_jsd']['median'].append(np.median(gene_jsds))
        gene_stats['gene_jsd']['p5'].append(np.percentile(gene_jsds, 5))
        gene_stats['gene_jsd']['p25'].append(np.percentile(gene_jsds, 25))
        gene_stats['gene_jsd']['p75'].append(np.percentile(gene_jsds, 75))
        gene_stats['gene_jsd']['p95'].append(np.percentile(gene_jsds, 95))
        gene_stats['gene_jsd']['min'].append(np.min(gene_jsds))
        gene_stats['gene_jsd']['max'].append(np.max(gene_jsds))
    
    # Create figure with secondary y-axis
    fig = make_subplots(specs=[[{"secondary_y": True}]])
    
    years = gene_stats['years']
    
    # Add percentile bands (matching cell-level style)
    # 5-95 percentile band (lighter)
    fig.add_trace(
        go.Scatter(
            x=years + years[::-1],
            y=gene_stats['gene_jsd']['p95'] + gene_stats['gene_jsd']['p5'][::-1],
            fill='toself',
            fillcolor='rgba(99, 110, 250, 0.15)',
            line=dict(color='rgba(255,255,255,0)'),
            name='5-95 percentile',
            showlegend=True,
            hoverinfo='skip'
        ),
        secondary_y=False
    )
    
    # 25-75 percentile band (darker)
    fig.add_trace(
        go.Scatter(
            x=years + years[::-1],
            y=gene_stats['gene_jsd']['p75'] + gene_stats['gene_jsd']['p25'][::-1],
            fill='toself',
            fillcolor='rgba(99, 110, 250, 0.25)',
            line=dict(color='rgba(255,255,255,0)'),
            name='25-75 percentile',
            showlegend=True,
            hoverinfo='skip'
        ),
        secondary_y=False
    )
    
    # Add mean line (blue, solid)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=gene_stats['gene_jsd']['mean'],
            mode='lines',
            name='Mean Gene JSD',
            line=dict(color='rgb(99, 110, 250)', width=2.5),
            hovertemplate='Year: %{x}<br>Mean Gene JSD: %{y:.4f}<extra></extra>'
        ),
        secondary_y=False
    )
    
    # Add median line (orange, dashed)
    fig.add_trace(
        go.Scatter(
            x=years,
            y=gene_stats['gene_jsd']['median'],
            mode='lines',
            name='Median Gene JSD',
            line=dict(color='rgb(255, 127, 14)', width=2, dash='dash'),
            hovertemplate='Year: %{x}<br>Median Gene JSD: %{y:.4f}<extra></extra>'
        ),
        secondary_y=False
    )
    
    # Add cell count on secondary y-axis
    fig.add_trace(
        go.Scatter(
            x=years,
            y=gene_stats['population_size'],
            mode='lines',
            name='Cell Count',
            line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
            hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
            showlegend=True
        ),
        secondary_y=True
    )
    
    # Add statistics annotation box (matching cell-level style)
    final_idx = -1
    annotation_text = (
        f"<b>Final Statistics (Year {gene_stats['years'][final_idx]})</b><br>"
        f"Population Size: {gene_stats['population_size'][final_idx]} cells<br>"
        f"Gene JSD Mean: {gene_stats['gene_jsd']['mean'][final_idx]:.4f}<br>"
        f"Gene JSD Median: {gene_stats['gene_jsd']['median'][final_idx]:.4f}<br>"
        f"Gene JSD 25-75%: [{gene_stats['gene_jsd']['p25'][final_idx]:.4f}, "
        f"{gene_stats['gene_jsd']['p75'][final_idx]:.4f}]<br>"
        f"Gene JSD 5-95%: [{gene_stats['gene_jsd']['p5'][final_idx]:.4f}, "
        f"{gene_stats['gene_jsd']['p95'][final_idx]:.4f}]"
    )
    
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
    
    # Update layout (matching cell-level style)
    default_title = "Gene JSD Trajectory"
    fig.update_layout(
        title=dict(
            text=title or default_title,
            font=dict(size=16)
        ),
        xaxis_title="Year",
        yaxis_title="Gene JSD Score",
        template="plotly_white",
        width=width,
        height=height,
        showlegend=True,
        legend=dict(
            yanchor="bottom",
            y=0.01,
            xanchor="right",
            x=0.99,
            bgcolor="rgba(255, 255, 255, 0.8)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1
        ),
        margin=dict(t=100),  # Add top margin for annotation
        hovermode='x unified'
    )
    
    # Update secondary y-axis
    fig.update_yaxes(title_text="Cell Count", secondary_y=True)
    
    # Save if path provided
    if output_path:
        os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
        fig.write_image(output_path, width=width, height=height, scale=2)
        print(f"Gene JSD trajectory saved to {output_path}")
    else:
        fig.show()
    
    return fig