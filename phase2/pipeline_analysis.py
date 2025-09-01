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
from pipeline_utils import load_petri_dish, load_all_petri_dishes, get_cell_jsd_array


def plot_cell_jsd_distribution(cells: List[Cell], bins: int, output_path: str, 
                              rate: Optional[float] = None, year: Optional[int] = None) -> None:
    """
    Create histogram of cell JSD distribution from Cell objects.
    
    Args:
        cells: List of Cell objects
        bins: Number of bins for histogram
        output_path: Path to save plot
        rate: Methylation rate (optional, for display)
        year: Year to display in title (optional, otherwise uses cell.age)
    """
    print(f"\nPlotting Cell JSD distribution...")
    
    # Extract cell JSD values directly from Cell objects
    cell_jsd_values = np.array([cell.cell_JSD for cell in cells])
    
    # Calculate statistics
    mean_jsd = np.mean(cell_jsd_values)
    std_jsd = np.std(cell_jsd_values)
    median_jsd = np.median(cell_jsd_values)
    cv_jsd = std_jsd / mean_jsd if mean_jsd > 0 else 0  # Coefficient of variation
    mad_jsd = np.median(np.abs(cell_jsd_values - median_jsd))  # Median Absolute Deviation
    p5_jsd = np.percentile(cell_jsd_values, 5)
    p25_jsd = np.percentile(cell_jsd_values, 25)
    p75_jsd = np.percentile(cell_jsd_values, 75)
    p95_jsd = np.percentile(cell_jsd_values, 95)
    
    # Create figure
    fig = go.Figure()
    
    # Calculate histogram data for step plot
    counts, bin_edges = np.histogram(cell_jsd_values, bins=bins)
    
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
        name='Cell JSD Distribution',
        line=dict(color='#1f77b4', width=2),
        fill='tozeroy',
        fillcolor='rgba(31, 119, 180, 0.3)',
        showlegend=False,
        hovertemplate='Cell JSD: %{x:.4f}<br>Count: %{y}<extra></extra>'
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
    # Use provided year if given, otherwise try cell.age
    if year is not None:
        display_year = year
    else:
        display_year = cells[0].age if cells else 'unknown'
    
    # Format rate as percentage if provided
    rate_text = ""
    if rate is not None:
        rate_percentage = rate * 100
        rate_text = f" | {rate_percentage:.1f}% methylation rate"
    
    fig.update_layout(
        title=dict(
            text=f"Cell JSD Distribution at Year {display_year}<br>"
                 f"<sub>{len(cells)} cells{rate_text}</sub>",
            font=dict(size=16)
        ),
        xaxis_title="Cell JSD Score",
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


def get_mean_cell_jsds_from_petri_dishes(dishes: List[PetriDish]) -> np.ndarray:
    """
    Extract mean cell JSD values from a list of PetriDish objects.
    
    Args:
        dishes: List of PetriDish objects
    
    Returns:
        Array of mean cell JSD values (one per dish)
    """
    mean_cell_jsds = []
    for petri in dishes:
        cell_jsd_values = [cell.cell_JSD for cell in petri.cells]
        mean_cell_jsd = np.mean(cell_jsd_values) if cell_jsd_values else 0.0
        mean_cell_jsds.append(mean_cell_jsd)
    
    return np.array(mean_cell_jsds)


def generate_gene_jsd_analysis(mutant_dishes: List[PetriDish],
                               control1_dishes: List[PetriDish], 
                               control2_dishes: List[PetriDish],
                               output_dir: str) -> Dict[str, Any]:
    """
    Generate gene-level JSD analysis for all genes across all batches.
    
    Args:
        mutant_dishes: List of mutant PetriDish objects
        control1_dishes: List of control1 PetriDish objects
        control2_dishes: List of control2 PetriDish objects
        output_dir: Directory to save results
        
    Returns:
        Dictionary with gene JSD analysis
    """
    print(f"\nGenerating gene JSD analysis...")
    
    # Calculate gene JSDs for each individual in each batch
    def calculate_batch_gene_jsds(dishes: List[PetriDish]) -> Dict[int, List[float]]:
        """Calculate gene JSDs for all genes across all individuals in a batch."""
        if not dishes or not dishes[0].cells:
            return {}
        
        # Get number of genes from first cell
        n_sites = dishes[0].cells[0].n
        gene_size = dishes[0].cells[0].gene_size
        n_genes = n_sites // gene_size
        
        # Initialize storage for each gene
        gene_jsds = {i: [] for i in range(n_genes)}
        
        # Calculate gene JSD for each individual
        for petri in dishes:
            individual_gene_jsds = petri.calculate_gene_jsd()
            for gene_idx, jsd_value in enumerate(individual_gene_jsds):
                gene_jsds[gene_idx].append(float(jsd_value))
        
        return gene_jsds
    
    # Calculate for each batch
    print("  Calculating gene JSDs for mutant batch...")
    mutant_gene_jsds = calculate_batch_gene_jsds(mutant_dishes)
    
    print("  Calculating gene JSDs for control1 batch...")
    control1_gene_jsds = calculate_batch_gene_jsds(control1_dishes)
    
    print("  Calculating gene JSDs for control2 batch...")
    control2_gene_jsds = calculate_batch_gene_jsds(control2_dishes)
    
    # Get number of genes
    n_genes = len(mutant_gene_jsds) if mutant_gene_jsds else 0
    
    # Calculate individual averages across all genes
    individual_averages = {
        "mutant": [],
        "control1": [],
        "control2": []
    }
    
    # Calculate average for each individual across all genes
    for batch_name, batch_jsds, n_individuals in [
        ("mutant", mutant_gene_jsds, len(mutant_dishes)),
        ("control1", control1_gene_jsds, len(control1_dishes)),
        ("control2", control2_gene_jsds, len(control2_dishes))
    ]:
        if batch_jsds and n_individuals > 0:
            # For each individual
            for ind_idx in range(n_individuals):
                # Collect all gene JSD values for this individual
                individual_gene_values = []
                for gene_idx in range(n_genes):
                    if gene_idx in batch_jsds and ind_idx < len(batch_jsds[gene_idx]):
                        individual_gene_values.append(batch_jsds[gene_idx][ind_idx])
                
                # Calculate mean across all genes for this individual
                if individual_gene_values:
                    individual_mean = float(np.mean(individual_gene_values))
                    individual_averages[batch_name].append(individual_mean)
    
    # Build gene JSD distributions structure (keep raw data)
    gene_jsd_distributions = {}
    for gene_idx in range(n_genes):
        gene_jsd_distributions[f"gene_{gene_idx}"] = {
            "mutant": mutant_gene_jsds.get(gene_idx, []),
            "control1": control1_gene_jsds.get(gene_idx, []),
            "control2": control2_gene_jsds.get(gene_idx, [])
        }
    
    # Build metadata
    gene_metadata = {
        "n_genes": n_genes,
        "n_individuals_per_batch": {
            "mutant": len(mutant_dishes),
            "control1": len(control1_dishes),
            "control2": len(control2_dishes)
        }
    }
    
    # Check if gene-specific rates were used
    if mutant_dishes and mutant_dishes[0].cells and hasattr(mutant_dishes[0].cells[0], 'gene_rate_groups'):
        gene_rate_groups = mutant_dishes[0].cells[0].gene_rate_groups
        if gene_rate_groups:
            # Convert to metadata format
            rate_groups = []
            start_idx = 0
            for n_genes_in_group, rate in gene_rate_groups:
                rate_groups.append({
                    "start": start_idx,
                    "end": start_idx + n_genes_in_group - 1,
                    "rate": float(rate)
                })
                start_idx += n_genes_in_group
            gene_metadata["gene_rate_groups"] = rate_groups
    
    # Create final structure with individual averages at the top
    gene_jsd_analysis = {
        "individual_averages": individual_averages,  # NEW: Individual-level averages across all genes
        "gene_jsd_distributions": gene_jsd_distributions,
        "gene_metadata": gene_metadata
    }
    
    # Save to file
    gene_jsd_path = os.path.join(output_dir, "gene_jsd_analysis.json")
    os.makedirs(output_dir, exist_ok=True)
    with open(gene_jsd_path, 'w') as f:
        json.dump(gene_jsd_analysis, f, indent=2)
    
    print(f"  Saved gene JSD analysis to {gene_jsd_path}")
    print(f"    Analyzed {n_genes} genes across {sum(gene_metadata['n_individuals_per_batch'].values())} individuals")
    
    # Print individual averages summary
    for batch_name in ["mutant", "control1", "control2"]:
        if individual_averages[batch_name]:
            mean_of_averages = np.mean(individual_averages[batch_name])
            print(f"    {batch_name.capitalize()}: {len(individual_averages[batch_name])} individuals, mean gene JSD = {mean_of_averages:.4f}")
    
    return {
        "gene_jsd_analysis": gene_jsd_analysis,
        "gene_jsd_path": gene_jsd_path
    }


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
    
    # Get mean cell JSD values for each group
    mutant_jsds = get_mean_cell_jsds_from_petri_dishes(mutant_dishes)
    control1_jsds = get_mean_cell_jsds_from_petri_dishes(control1_dishes)
    control2_jsds = get_mean_cell_jsds_from_petri_dishes(control2_dishes)
    
    print(f"  Mutant: {len(mutant_jsds)} individuals")
    print(f"  Control1: {len(control1_jsds)} individuals")
    print(f"  Control2: {len(control2_jsds)} individuals")
    
    # Create consolidated cell JSD analysis structure
    cell_jsd_analysis = {
        "summary_statistics": {
            "mutant": {
                "mean": float(np.mean(mutant_jsds)),
                "std": float(np.std(mutant_jsds)),
                "median": float(np.median(mutant_jsds)),
                "min": float(np.min(mutant_jsds)),
                "max": float(np.max(mutant_jsds)),
                "n_individuals": len(mutant_jsds)
            },
            "control1": {
                "mean": float(np.mean(control1_jsds)),
                "std": float(np.std(control1_jsds)),
                "median": float(np.median(control1_jsds)),
                "min": float(np.min(control1_jsds)),
                "max": float(np.max(control1_jsds)),
                "n_individuals": len(control1_jsds)
            },
            "control2": {
                "mean": float(np.mean(control2_jsds)),
                "std": float(np.std(control2_jsds)),
                "median": float(np.median(control2_jsds)),
                "min": float(np.min(control2_jsds)),
                "max": float(np.max(control2_jsds)),
                "n_individuals": len(control2_jsds)
            }
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
    
    cell_jsd_analysis["statistical_tests"] = {
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
    
    # Add individual means (renamed from distributions for clarity)
    cell_jsd_analysis["individual_means"] = {
        "mutant": mutant_jsds.tolist(),
        "control1": control1_jsds.tolist(),
        "control2": control2_jsds.tolist()
    }
    
    # Create comparison plot for cell JSDs
    plot_path = os.path.join(output_dir, "cell_jsd_comparison.png")
    create_comparison_plot_from_jsds(mutant_jsds, control1_jsds, control2_jsds, plot_path)
    
    # Save consolidated cell JSD analysis
    cell_jsd_path = os.path.join(output_dir, "cell_jsd_analysis.json")
    os.makedirs(output_dir, exist_ok=True)
    with open(cell_jsd_path, 'w') as f:
        json.dump(cell_jsd_analysis, f, indent=2)
    print(f"\n  Saved cell JSD analysis to {cell_jsd_path}")
    
    return {
        "cell_jsd_analysis": cell_jsd_analysis,
        "plot_path": plot_path,
        "cell_jsd_path": cell_jsd_path
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
    
    # Calculate y range for plot
    all_jsds = np.concatenate([mutant_jsds, control1_jsds, control2_jsds])
    y_min = np.min(all_jsds)
    y_max = np.max(all_jsds)
    
    # Add comprehensive statistics text ABOVE the plot using paper coordinates
    # Position them in three columns at the top
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.2, y=1.02,  # Left column for Mutant
        text=(f"<b>Mutant</b><br>"
              f"Mean: {mutant_mean:.4f}<br>"
              f"Median: {mutant_median:.4f}<br>"
              f"SD: {mutant_std:.4f}<br>"
              f"CV: {mutant_cv:.3f}<br>"
              f"MAD: {mutant_mad:.4f}<br>"
              f"5%: {mutant_q5:.4f}<br>"
              f"95%: {mutant_q95:.4f}"),
        showarrow=False,
        font=dict(size=9, color="#1f77b4"),
        align="left",
        xanchor="center",
        yanchor="bottom"
    )
    
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.5, y=1.02,  # Middle column for Control1
        text=(f"<b>Control1</b><br>"
              f"Mean: {control1_mean:.4f}<br>"
              f"Median: {control1_median:.4f}<br>"
              f"SD: {control1_std:.4f}<br>"
              f"CV: {control1_cv:.3f}<br>"
              f"MAD: {control1_mad:.4f}<br>"
              f"5%: {control1_q5:.4f}<br>"
              f"95%: {control1_q95:.4f}"),
        showarrow=False,
        font=dict(size=9, color="#ff7f0e"),
        align="left",
        xanchor="center",
        yanchor="bottom"
    )
    
    fig.add_annotation(
        xref="paper", yref="paper",
        x=0.8, y=1.02,  # Right column for Control2
        text=(f"<b>Control2</b><br>"
              f"Mean: {control2_mean:.4f}<br>"
              f"Median: {control2_median:.4f}<br>"
              f"SD: {control2_std:.4f}<br>"
              f"CV: {control2_cv:.3f}<br>"
              f"MAD: {control2_mad:.4f}<br>"
              f"5%: {control2_q5:.4f}<br>"
              f"95%: {control2_q95:.4f}"),
        showarrow=False,
        font=dict(size=9, color="#2ca02c"),
        align="left",
        xanchor="center",
        yanchor="bottom"
    )
    
    # Clean layout with proper spacing
    fig.update_layout(
        title=dict(
            text="Mean Cell JSD per Individual",
            font=dict(size=16),
            x=0.5,
            xanchor='center',
            y=0.98,  # Move title up (default is ~0.9)
            yanchor='top'
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
            title='Mean Cell JSD per Individual',
            showgrid=True,
            gridcolor='lightgray',
            range=[y_min - 0.01, y_max + 0.01]  # Simple range without annotation adjustment
        ),
        plot_bgcolor='white',
        showlegend=False,
        height=800,  # Keep height for good proportions
        width=700,
        margin=dict(l=80, r=40, t=180, b=60)  # Increased TOP margin for stats above
    )
    
    # Save PNG only (no HTML)
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_image(output_path, width=700, height=800, scale=2)
    print(f"    Saved plot to {output_path}")


def plot_gene_jsd_distribution_comparison(snapshot1_cells, snapshot2_cells, 
                                          year1, year2, output_path):
    """
    Overlapping histograms comparing gene JSD distributions at two time points.
    
    Args:
        snapshot1_cells: List of Cell objects from first snapshot
        snapshot2_cells: List of Cell objects from second snapshot
        year1: Year of first snapshot
        year2: Year of second snapshot
        output_path: Path to save the plot
    """
    # Import plotly here to avoid dependency issues
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
    
    print(f"\n  Creating gene JSD distribution comparison plot...")
    
    # Calculate gene JSD for each snapshot
    # Create temporary PetriDish objects to calculate gene JSD
    from cell import PetriDish
    
    petri1 = PetriDish()
    petri1.cells = snapshot1_cells
    gene_jsds1 = petri1.calculate_gene_jsd()
    
    petri2 = PetriDish()
    petri2.cells = snapshot2_cells
    gene_jsds2 = petri2.calculate_gene_jsd()
    
    print(f"    Year {year1}: {len(gene_jsds1)} genes")
    print(f"    Year {year2}: {len(gene_jsds2)} genes")
    
    # Calculate statistics
    mean1 = np.mean(gene_jsds1)
    mean2 = np.mean(gene_jsds2)
    median1 = np.median(gene_jsds1)
    median2 = np.median(gene_jsds2)
    max1 = np.max(gene_jsds1)
    max2 = np.max(gene_jsds2)
    
    # Create figure
    fig = go.Figure()
    
    # Add histogram for year 1
    fig.add_trace(go.Histogram(
        x=gene_jsds1,
        name=f'Year {year1}',
        opacity=0.5,
        marker_color='#1f77b4',
        nbinsx=50,
        histnorm='',
        hovertemplate=f'Year {year1}<br>Gene JSD: %{{x:.4f}}<br>Count: %{{y}}<extra></extra>'
    ))
    
    # Add histogram for year 2
    fig.add_trace(go.Histogram(
        x=gene_jsds2,
        name=f'Year {year2}',
        opacity=0.5,
        marker_color='#ff7f0e',
        nbinsx=50,
        histnorm='',
        hovertemplate=f'Year {year2}<br>Gene JSD: %{{x:.4f}}<br>Count: %{{y}}<extra></extra>'
    ))
    
    # Add statistics annotation
    stats_text = (f"<b>Year {year1} Stats</b><br>"
                  f"Mean: {mean1:.3f}<br>"
                  f"Median: {median1:.3f}<br>"
                  f"Max: {max1:.3f}<br><br>"
                  f"<b>Year {year2} Stats</b><br>"
                  f"Mean: {mean2:.3f}<br>"
                  f"Median: {median2:.3f}<br>"
                  f"Max: {max2:.3f}")
    
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
    fig.update_layout(
        title=dict(
            text=f"Gene JSD Distribution: Year {year1} vs Year {year2}",
            font=dict(size=16)
        ),
        xaxis_title="Gene JSD Score",
        yaxis_title="Number of Genes",
        barmode='overlay',
        template="plotly_white",
        width=1200,
        height=600,
        showlegend=True,
        legend=dict(
            x=0.7,
            y=0.3,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="#333333",
            borderwidth=1
        )
    )
    
    # Update axes with grid
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    
    # Save
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_image(output_path, width=1200, height=600, scale=2)
    print(f"    Saved gene JSD distribution comparison to {output_path}")


def plot_top_variable_genes(petri_dishes, n_top=20, output_path=None):
    """
    Bar chart showing the most heterogeneous genes.
    
    Args:
        petri_dishes: List of PetriDish objects to analyze
        n_top: Number of top genes to show (default: 20)
        output_path: Path to save the plot (optional)
    """
    # Import plotly here to avoid dependency issues
    try:
        import plotly.graph_objects as go
    except ImportError:
        raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
    
    print(f"\n  Creating top {n_top} variable genes plot...")
    
    if not petri_dishes:
        raise ValueError("No PetriDish objects provided")
    
    # Calculate gene JSD for all dishes and average across individuals
    all_gene_jsds = []
    for petri in petri_dishes:
        gene_jsds = petri.calculate_gene_jsd()
        all_gene_jsds.append(gene_jsds)
    
    # Convert to numpy array and calculate mean across individuals
    all_gene_jsds = np.array(all_gene_jsds)  # Shape: (n_individuals, n_genes)
    mean_gene_jsds = np.mean(all_gene_jsds, axis=0)  # Mean across individuals
    
    n_genes = len(mean_gene_jsds)
    print(f"    Analyzing {n_genes} genes across {len(petri_dishes)} individuals")
    
    # Get top variable genes
    top_indices = np.argsort(mean_gene_jsds)[-n_top:][::-1]  # Top n_top, highest first
    top_jsds = mean_gene_jsds[top_indices]
    top_gene_names = [f"Gene {i}" for i in top_indices]
    
    # Check if we have gene rate groups to add rate information
    rate_info = []
    if (petri_dishes and petri_dishes[0].cells and 
        hasattr(petri_dishes[0].cells[0], 'gene_rate_groups') and 
        petri_dishes[0].cells[0].gene_rate_groups):
        
        gene_rate_groups = petri_dishes[0].cells[0].gene_rate_groups
        gene_to_group = {}
        gene_start_idx = 0
        
        for group_idx, (n_genes_in_group, rate) in enumerate(gene_rate_groups):
            for gene_idx in range(gene_start_idx, gene_start_idx + n_genes_in_group):
                gene_to_group[gene_idx] = (group_idx + 1, rate * 100)  # Convert to percentage
            gene_start_idx += n_genes_in_group
        
        # Add rate info to gene names
        for i, gene_idx in enumerate(top_indices):
            if gene_idx in gene_to_group:
                group_num, rate_pct = gene_to_group[gene_idx]
                top_gene_names[i] = f"Gene {gene_idx} (Group {group_num}: {rate_pct:.1f}%)"
                rate_info.append(f"Group {group_num}")
            else:
                rate_info.append("Unknown")
    else:
        rate_info = ["N/A"] * n_top
    
    # Create color gradient based on JSD values (light to dark blue)
    # Normalize JSD values to 0-1 range for color mapping
    if len(top_jsds) > 0:
        jsd_min = np.min(top_jsds)
        jsd_max = np.max(top_jsds)
        jsd_normalized = (top_jsds - jsd_min) / (jsd_max - jsd_min) if jsd_max > jsd_min else np.ones_like(top_jsds)
    else:
        jsd_normalized = np.array([])
    
    # Generate colors (light blue to dark blue)
    colors = []
    for norm_val in jsd_normalized:
        # Interpolate between light blue and dark blue
        light_blue = [173, 216, 230]  # Light blue RGB
        dark_blue = [25, 25, 112]     # Dark blue RGB
        color_rgb = [
            int(light_blue[i] + (dark_blue[i] - light_blue[i]) * norm_val)
            for i in range(3)
        ]
        colors.append(f'rgb({color_rgb[0]}, {color_rgb[1]}, {color_rgb[2]})')
    
    # Create figure
    fig = go.Figure()
    
    # Add horizontal bar chart
    fig.add_trace(go.Bar(
        x=top_jsds,
        y=top_gene_names,
        orientation='h',
        marker=dict(
            color=colors,
            line=dict(color='white', width=0.5)
        ),
        hovertemplate='%{y}<br>Gene JSD: %{x:.4f}<extra></extra>'
    ))
    
    # Add statistics annotation
    stats_text = (f"<b>Top {n_top} Variable Genes</b><br>"
                  f"Total Genes: {n_genes}<br>"
                  f"Individuals: {len(petri_dishes)}<br>"
                  f"Highest JSD: {np.max(top_jsds):.4f}<br>"
                  f"Mean of Top {n_top}: {np.mean(top_jsds):.4f}")
    
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
    fig.update_layout(
        title=dict(
            text=f"Top {n_top} Most Heterogeneous Genes",
            font=dict(size=16)
        ),
        xaxis_title="Gene JSD Score",
        yaxis_title="",
        template="plotly_white",
        width=1200,
        height=max(600, 25 * n_top),  # Adjust height based on number of genes
        showlegend=False,
        margin=dict(l=250, r=40, t=60, b=60)  # Extra left margin for gene names
    )
    
    # Update axes with grid
    fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
    fig.update_yaxes(showgrid=False)
    
    # Save
    if output_path:
        if os.path.dirname(output_path):
            os.makedirs(os.path.dirname(output_path), exist_ok=True)
        fig.write_image(output_path, width=1200, height=max(600, 25 * n_top), scale=2)
        print(f"    Saved top variable genes plot to {output_path}")
    else:
        fig.show()
    
    return fig


def plot_gene_jsd_distributions(gene_jsd_data: Dict[str, Any], output_dir: str, max_genes: int = None) -> None:
    """
    Create individual gene JSD plots with exact same design as cell-level JSD plot.
    
    Args:
        gene_jsd_data: The gene_jsd_analysis dictionary with distributions and metadata
        output_dir: Directory to save plots (results/)
        max_genes: Maximum number of genes to plot (None for all)
    """
    print(f"\nGenerating gene JSD distribution plots...")
    
    # Extract data
    gene_distributions = gene_jsd_data['gene_jsd_distributions']
    metadata = gene_jsd_data['gene_metadata']
    n_genes = metadata['n_genes']
    
    # Determine how many genes to plot
    genes_to_plot = min(n_genes, max_genes) if max_genes else n_genes
    
    # Create output directory
    gene_plots_dir = os.path.join(output_dir, 'gene_jsd_plots')
    os.makedirs(gene_plots_dir, exist_ok=True)
    
    print(f"  Creating plots for {genes_to_plot} genes...")
    
    # Process each gene
    for gene_idx in range(genes_to_plot):
        gene_key = f"gene_{gene_idx}"
        if gene_key not in gene_distributions:
            continue
        
        gene_data = gene_distributions[gene_key]
        
        # Extract arrays of JSD values for this gene
        mutant_jsds = np.array(gene_data.get('mutant', []))
        control1_jsds = np.array(gene_data.get('control1', []))
        control2_jsds = np.array(gene_data.get('control2', []))
        
        # Skip if no data
        if len(mutant_jsds) == 0 and len(control1_jsds) == 0 and len(control2_jsds) == 0:
            continue
        
        # Create plot using EXACT same design as cell-level
        plot_path = os.path.join(gene_plots_dir, f'gene_{gene_idx:03d}_jsd_comparison.png')
        create_gene_comparison_plot(mutant_jsds, control1_jsds, control2_jsds, 
                                   plot_path, gene_idx, metadata)
        
        # Progress indicator
        if (gene_idx + 1) % 10 == 0:
            print(f"    Processed {gene_idx + 1}/{genes_to_plot} genes...")
    
    print(f"  ✓ Generated {genes_to_plot} gene JSD plots in {gene_plots_dir}")


def create_gene_comparison_plot(mutant_jsds: np.ndarray, 
                               control1_jsds: np.ndarray, 
                               control2_jsds: np.ndarray, 
                               output_path: str,
                               gene_idx: int,
                               metadata: Dict[str, Any]) -> None:
    """
    Create a scatter plot for a single gene's JSD distribution.
    Exact copy of create_comparison_plot_from_jsds but for gene-level data.
    
    Args:
        mutant_jsds: Array of JSD values for mutant group
        control1_jsds: Array of JSD values for control1 group
        control2_jsds: Array of JSD values for control2 group
        output_path: Path to save the plot
        gene_idx: Index of the gene being plotted
        metadata: Gene metadata including rate groups if present
    """
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
    if len(mutant_jsds) > 0:
        mutant_mean = np.mean(mutant_jsds)
        mutant_q5, mutant_q25, mutant_q75, mutant_q95 = np.percentile(mutant_jsds, [5, 25, 75, 95])
        mutant_median = np.median(mutant_jsds)
        mutant_std = np.std(mutant_jsds)
        mutant_cv = mutant_std / mutant_mean if mutant_mean != 0 else 0
        mutant_mad = np.median(np.abs(mutant_jsds - mutant_median))
    else:
        mutant_mean = mutant_median = mutant_std = mutant_cv = mutant_mad = 0
        mutant_q5 = mutant_q25 = mutant_q75 = mutant_q95 = 0
    
    if len(control1_jsds) > 0:
        control1_mean = np.mean(control1_jsds)
        control1_q5, control1_q25, control1_q75, control1_q95 = np.percentile(control1_jsds, [5, 25, 75, 95])
        control1_median = np.median(control1_jsds)
        control1_std = np.std(control1_jsds)
        control1_cv = control1_std / control1_mean if control1_mean != 0 else 0
        control1_mad = np.median(np.abs(control1_jsds - control1_median))
    else:
        control1_mean = control1_median = control1_std = control1_cv = control1_mad = 0
        control1_q5 = control1_q25 = control1_q75 = control1_q95 = 0
    
    if len(control2_jsds) > 0:
        control2_mean = np.mean(control2_jsds)
        control2_q5, control2_q25, control2_q75, control2_q95 = np.percentile(control2_jsds, [5, 25, 75, 95])
        control2_median = np.median(control2_jsds)
        control2_std = np.std(control2_jsds)
        control2_cv = control2_std / control2_mean if control2_mean != 0 else 0
        control2_mad = np.median(np.abs(control2_jsds - control2_median))
    else:
        control2_mean = control2_median = control2_std = control2_cv = control2_mad = 0
        control2_q5 = control2_q25 = control2_q75 = control2_q95 = 0
    
    # Add quantile lines for mutant (if data exists)
    if len(mutant_jsds) > 0:
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
        
        # Mean line (solid, prominent)
        fig.add_shape(type="line", x0=-0.2, x1=0.2, y0=mutant_mean, y1=mutant_mean,
                      line=dict(color="#1f77b4", width=3))
    
    # Add quantile lines for control1 (if data exists)
    if len(control1_jsds) > 0:
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
        
        # Mean line (solid, prominent)
        fig.add_shape(type="line", x0=0.8, x1=1.2, y0=control1_mean, y1=control1_mean,
                      line=dict(color="#ff7f0e", width=3))
    
    # Add quantile lines for control2 (if data exists)
    if len(control2_jsds) > 0:
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
        
        # Mean line (solid, prominent)
        fig.add_shape(type="line", x0=1.8, x1=2.2, y0=control2_mean, y1=control2_mean,
                      line=dict(color="#2ca02c", width=3))
    
    # Calculate y range for plot
    all_jsds = np.concatenate([jsds for jsds in [mutant_jsds, control1_jsds, control2_jsds] if len(jsds) > 0])
    if len(all_jsds) > 0:
        y_min = np.min(all_jsds)
        y_max = np.max(all_jsds)
    else:
        y_min = 0
        y_max = 1
    
    # Add comprehensive statistics text ABOVE the plot using paper coordinates
    # Position them in three columns at the top
    if len(mutant_jsds) > 0:
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.2, y=1.02,  # Left column for Mutant
            text=(f"<b>Mutant</b><br>"
                  f"Mean: {mutant_mean:.4f}<br>"
                  f"Median: {mutant_median:.4f}<br>"
                  f"SD: {mutant_std:.4f}<br>"
                  f"CV: {mutant_cv:.3f}<br>"
                  f"MAD: {mutant_mad:.4f}<br>"
                  f"5%: {mutant_q5:.4f}<br>"
                  f"95%: {mutant_q95:.4f}"),
            showarrow=False,
            font=dict(size=9, color="#1f77b4"),
            align="left",
            xanchor="center",
            yanchor="bottom"
        )
    
    if len(control1_jsds) > 0:
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.5, y=1.02,  # Middle column for Control1
            text=(f"<b>Control1</b><br>"
                  f"Mean: {control1_mean:.4f}<br>"
                  f"Median: {control1_median:.4f}<br>"
                  f"SD: {control1_std:.4f}<br>"
                  f"CV: {control1_cv:.3f}<br>"
                  f"MAD: {control1_mad:.4f}<br>"
                  f"5%: {control1_q5:.4f}<br>"
                  f"95%: {control1_q95:.4f}"),
            showarrow=False,
            font=dict(size=9, color="#ff7f0e"),
            align="left",
            xanchor="center",
            yanchor="bottom"
        )
    
    if len(control2_jsds) > 0:
        fig.add_annotation(
            xref="paper", yref="paper",
            x=0.8, y=1.02,  # Right column for Control2
            text=(f"<b>Control2</b><br>"
                  f"Mean: {control2_mean:.4f}<br>"
                  f"Median: {control2_median:.4f}<br>"
                  f"SD: {control2_std:.4f}<br>"
                  f"CV: {control2_cv:.3f}<br>"
                  f"MAD: {control2_mad:.4f}<br>"
                  f"5%: {control2_q5:.4f}<br>"
                  f"95%: {control2_q95:.4f}"),
            showarrow=False,
            font=dict(size=9, color="#2ca02c"),
            align="left",
            xanchor="center",
            yanchor="bottom"
        )
    
    # Add gene rate annotation if applicable
    gene_rate_text = ""
    if 'gene_rate_groups' in metadata and metadata['gene_rate_groups']:
        for group in metadata['gene_rate_groups']:
            if group['start'] <= gene_idx <= group['end']:
                gene_rate_text = f" | Rate: {group['rate']*100:.2f}%"
                break
    
    # Clean layout with proper spacing (EXACT same as cell-level)
    fig.update_layout(
        title=dict(
            text=f"Gene {gene_idx} JSD Distribution{gene_rate_text}",
            font=dict(size=16),
            x=0.5,
            xanchor='center',
            y=0.98,  # Move title up (default is ~0.9)
            yanchor='top'
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
            title='Gene JSD Score',
            showgrid=True,
            gridcolor='lightgray',
            range=[y_min - 0.01, y_max + 0.01] if len(all_jsds) > 0 else [0, 1]  # Simple range
        ),
        plot_bgcolor='white',
        showlegend=False,
        height=800,  # Keep height for good proportions
        width=700,
        margin=dict(l=80, r=40, t=180, b=60)  # Increased TOP margin for stats above
    )
    
    # Save PNG only
    if os.path.dirname(output_path):
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
    fig.write_image(output_path, width=700, height=800, scale=2)


def plot_gene_jsd_individual_comparison(gene_jsd_path, output_dir, verbose=False):
    """
    Create scatter plot comparing individual-averaged gene JSDs across batches.
    
    Similar to plot_cell_jsd_comparison but uses individual_averages from gene_jsd_analysis.json
    
    Args:
        gene_jsd_path: Path to gene_jsd_analysis.json
        output_dir: Directory to save plot
        verbose: Print progress messages
    """
    if verbose:
        print(f"\n📊 Creating individual-averaged gene JSD comparison plot...")
    
    # Load gene JSD analysis data
    with open(gene_jsd_path, 'r') as f:
        data = json.load(f)
    
    # Extract individual averages
    if 'individual_averages' not in data:
        print("  ⚠️  No individual_averages found in gene JSD analysis")
        return
    
    individual_avgs = data['individual_averages']
    
    # Create figure with plotly
    fig = go.Figure()
    
    # Color mapping (same as cell plot)
    colors = {
        'mutant': 'blue',
        'control1': 'green', 
        'control2': 'red'
    }
    
    # Track all values for y-axis range
    all_values = []
    
    # Plot each batch with jittered points
    for i, (batch_name, batch_values) in enumerate(individual_avgs.items()):
        if not batch_values:
            continue
            
        # Create jittered x positions
        x_base = i
        x_positions = np.random.normal(x_base, 0.05, len(batch_values))
        
        # Add scatter trace
        fig.add_trace(go.Scatter(
            x=x_positions,
            y=batch_values,
            mode='markers',
            name=batch_name.capitalize(),
            marker=dict(
                color=colors.get(batch_name, 'gray'),
                size=6,
                opacity=0.6
            ),
            hovertemplate='%{y:.4f}<extra></extra>'
        ))
        
        all_values.extend(batch_values)
    
    # Calculate statistics for each batch
    stats_y = 1.02  # Position above plot using paper coordinates
    
    # Mutant stats
    if 'mutant' in individual_avgs and individual_avgs['mutant']:
        mutant_mean = np.mean(individual_avgs['mutant'])
        mutant_std = np.std(individual_avgs['mutant'])
        fig.add_annotation(
            x=0,
            y=stats_y,
            text=f"Mean: {mutant_mean:.4f} ± {mutant_std:.4f}",
            showarrow=False,
            font=dict(size=9, color="blue"),
            align="left",
            xanchor="center",
            yanchor="bottom",
            xref="x",
            yref="paper"
        )
    
    # Control1 stats
    if 'control1' in individual_avgs and individual_avgs['control1']:
        control1_mean = np.mean(individual_avgs['control1'])
        control1_std = np.std(individual_avgs['control1'])
        fig.add_annotation(
            x=1,
            y=stats_y,
            text=f"Mean: {control1_mean:.4f} ± {control1_std:.4f}",
            showarrow=False,
            font=dict(size=9, color="green"),
            align="left",
            xanchor="center",
            yanchor="bottom",
            xref="x",
            yref="paper"
        )
    
    # Control2 stats
    if 'control2' in individual_avgs and individual_avgs['control2']:
        control2_mean = np.mean(individual_avgs['control2'])
        control2_std = np.std(individual_avgs['control2'])
        fig.add_annotation(
            x=2,
            y=stats_y,
            text=f"Mean: {control2_mean:.4f} ± {control2_std:.4f}",
            showarrow=False,
            font=dict(size=9, color="red"),
            align="left",
            xanchor="center",
            yanchor="bottom",
            xref="x",
            yref="paper"
        )
    
    # Calculate y-axis range
    if all_values:
        y_min = min(all_values)
        y_max = max(all_values)
        y_range = [y_min - 0.01, y_max + 0.01]
    else:
        y_range = [0, 1]
    
    # Update layout (matching cell plot style)
    fig.update_layout(
        title=dict(
            text="Individual-Averaged Gene JSD Comparison",
            font=dict(size=16),
            x=0.5,
            xanchor='center',
            y=0.98,  # High position like cell plot
            yanchor='top'
        ),
        xaxis=dict(
            tickmode='array',
            tickvals=[0, 1, 2],
            ticktext=['Mutant', 'Control 1', 'Control 2'],
            range=[-0.5, 2.5],
            showgrid=False,
            title=""
        ),
        yaxis=dict(
            title='Individual Average Gene JSD',
            showgrid=True,
            gridcolor='lightgray',
            range=y_range
        ),
        plot_bgcolor='white',
        showlegend=False,
        height=800,
        width=1200,
        margin=dict(l=80, r=40, t=180, b=60)  # Increased top margin for stats
    )
    
    # Save plot
    output_path = os.path.join(output_dir, 'gene_jsd_individual_comparison.png')
    os.makedirs(output_dir, exist_ok=True)
    fig.write_image(output_path, format='png', scale=1.5)
    
    if verbose:
        print(f"  ✓ Saved individual-averaged gene JSD comparison to {output_path}")