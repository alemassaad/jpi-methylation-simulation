#!/usr/bin/env python3
"""
Analyze populations and generate all plots.
This script performs all analysis and visualization tasks.
"""

import argparse
import os
import sys
import json
import gzip
from typing import List, Dict

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import PetriDish, Cell, PetriDishPlotter
from core.pipeline_utils import (
    load_snapshot_cells, load_all_petri_dishes, load_petri_dish,
    extract_gene_jsd_from_history, smart_open
)
from core.pipeline_analysis import (
    plot_cell_jsd_distribution,
    plot_cell_methylation_proportion_histogram,
    plot_gene_jsd_distribution,
    plot_gene_methylation_proportion_histogram,
    analyze_cell_methylation_proportion_comparison,
    analyze_populations_from_dishes,
    plot_gene_jsd_distributions,
    plot_gene_jsd_individual_comparison,
    generate_gene_methylation_analysis,
    plot_gene_methylation_individual_comparison,
    generate_gene_jsd_analysis
)
from core.plot_paths import PlotPaths


def load_metadata(base_dir: str) -> Dict:
    """Load metadata from extraction stage."""
    metadata_path = os.path.join(base_dir, "snapshots", "metadata.json")
    with open(metadata_path, 'r') as f:
        return json.load(f)


def plot_snapshot_distributions(
    cells: List[Cell],
    year: int,
    bins: int,
    plot_paths: PlotPaths,
    gene_rate_groups: List
) -> None:
    """Plot all distributions for a snapshot."""
    print(f"\nPlotting distributions for year {year} snapshot...")
    
    # Cell JSD distribution
    plot_path = plot_paths.get_cell_jsd_distribution_path(year)
    plot_cell_jsd_distribution(cells, bins, plot_path, 
                              gene_rate_groups=gene_rate_groups, year=year)
    print(f"  ✓ Cell JSD distribution")
    
    # Cell methylation proportion
    plot_path = plot_paths.get_cell_methylation_distribution_path(year)
    plot_cell_methylation_proportion_histogram(cells, bins, plot_path,
                                              gene_rate_groups=gene_rate_groups, year=year)
    print(f"  ✓ Cell methylation proportion")
    
    # Gene JSD distribution
    plot_path = plot_paths.get_gene_jsd_distribution_path(year)
    plot_gene_jsd_distribution(cells, bins=20, output_path=plot_path,
                              gene_rate_groups=gene_rate_groups, year=year)
    print(f"  ✓ Gene JSD distribution")
    
    # Gene methylation proportion
    plot_path = plot_paths.get_gene_methylation_distribution_path(year)
    plot_gene_methylation_proportion_histogram(cells, bins=20,
                                              output_path=plot_path,
                                              gene_rate_groups=gene_rate_groups,
                                              year=year)
    print(f"  ✓ Gene methylation proportion")


def plot_individual_trajectories(base_dir: str, plot_paths: PlotPaths) -> None:
    """Plot growth trajectories for mutant and control1 individuals."""
    print("\nGenerating Individual Growth Trajectories...")
    
    individuals_dir = os.path.join(base_dir, "individuals")
    
    # Plot mutant individuals
    mutant_dir = os.path.join(individuals_dir, "mutant")
    print("\n  Creating plots for mutant individuals...")
    mutant_plot_count = 0
    
    for i, filename in enumerate(sorted(os.listdir(mutant_dir)), 1):
        if filename.endswith(('.json', '.json.gz')):
            filepath = os.path.join(mutant_dir, filename)
            try:
                # Load with history
                petri = load_petri_dish(filepath, include_cell_history=True)
                
                # Check if history exists
                if hasattr(petri, 'cell_history') and petri.cell_history:
                    plotter = PetriDishPlotter(petri)
                    
                    # Generate plots
                    jsd_path = plot_paths.get_individual_cell_jsd_path('mutant', i)
                    meth_path = plot_paths.get_individual_cell_methylation_path('mutant', i)
                    gene_jsd_path = plot_paths.get_individual_gene_jsd_path('mutant', i)
                    
                    plotter.plot_jsd(f"Mutant Individual {i:02d}", jsd_path)
                    plotter.plot_cell_methylation_proportion(
                        f"Mutant Individual {i:02d} Cell Methylation Proportion", meth_path
                    )
                    plotter.plot_gene_jsd_trajectory(
                        f"Mutant Individual {i:02d} Gene JSD", gene_jsd_path
                    )
                    
                    mutant_plot_count += 1
                    print(f"    ✓ Plotted mutant_{i:02d}")
            except Exception as e:
                print(f"    ⚠️  Could not plot mutant_{i:02d}: {str(e)}")
    
    # Plot control1 individuals
    control1_dir = os.path.join(individuals_dir, "control1")
    print("\n  Creating plots for control1 individuals...")
    control1_plot_count = 0
    
    for i, filename in enumerate(sorted(os.listdir(control1_dir)), 1):
        if filename.endswith(('.json', '.json.gz')):
            filepath = os.path.join(control1_dir, filename)
            try:
                # Load with history
                petri = load_petri_dish(filepath, include_cell_history=True)
                
                # Check if history exists
                if hasattr(petri, 'cell_history') and petri.cell_history:
                    plotter = PetriDishPlotter(petri)
                    
                    # Generate plots
                    jsd_path = plot_paths.get_individual_cell_jsd_path('control1', i)
                    meth_path = plot_paths.get_individual_cell_methylation_path('control1', i)
                    gene_jsd_path = plot_paths.get_individual_gene_jsd_path('control1', i)
                    
                    plotter.plot_jsd(f"Control1 Individual {i:02d}", jsd_path)
                    plotter.plot_cell_methylation_proportion(
                        f"Control1 Individual {i:02d} Cell Methylation Proportion", meth_path
                    )
                    plotter.plot_gene_jsd_trajectory(
                        f"Control1 Individual {i:02d} Gene JSD", gene_jsd_path
                    )
                    
                    control1_plot_count += 1
                    print(f"    ✓ Plotted control1_{i:02d}")
            except Exception as e:
                print(f"    ⚠️  Could not plot control1_{i:02d}: {str(e)}")
    
    print(f"\n  Skipping control2 (pure snapshots, no growth trajectories)")
    total_plots = (mutant_plot_count + control1_plot_count) * 3
    print(f"  ✓ Generated {total_plots} trajectory plots")


def generate_timeline_plots(simulation_path: str, results_dir: str) -> None:
    """Generate timeline plots from original simulation."""
    print("\nGenerating Original Simulation Timeline Plots...")
    
    try:
        print(f"  Loading simulation with full history...")
        
        # Load simulation data
        with smart_open(simulation_path, 'r') as f:
            sim_data = json.load(f)
        
        # Convert history to PetriDish format
        if 'history' in sim_data and sim_data['history']:
            from core.pipeline_utils import dict_to_cell
            
            # Get parameters from simulation
            params = sim_data.get('parameters', {})
            rate = params.get('rate')
            gene_rate_groups = params.get('gene_rate_groups')
            n_sites = params.get('n', 1000)
            gene_size = params.get('gene_size', 5)
            
            # Create PetriDish with appropriate configuration
            if gene_rate_groups:
                rate_groups = [(g[0], g[1]) for g in gene_rate_groups]
                original_petri = PetriDish(gene_rate_groups=rate_groups, n=n_sites, gene_size=gene_size)
            else:
                original_petri = PetriDish(rate=rate, n=n_sites, gene_size=gene_size)
            
            # Build cell history from simulation data
            original_petri.cell_history = {}
            original_petri.years_simulated = len(sim_data['history']) - 1
            
            for year_str, year_data in sim_data['history'].items():
                if 'cells' in year_data:
                    cells = []
                    for cell_dict in year_data['cells']:
                        # Handle new lean format
                        if 'cpg_sites' in cell_dict:
                            cell_data = {
                                'methylated': cell_dict['cpg_sites'],
                                'cell_jsd': cell_dict.get('cell_jsd', cell_dict.get('cell_JSD', 0.0)),
                                'cell_methylation_proportion': sum(cell_dict['cpg_sites']) / len(cell_dict['cpg_sites']),
                                'age': cell_dict.get('age', int(year_str))
                            }
                        else:
                            cell_data = {
                                'methylated': cell_dict.get('methylated', []),
                                'cell_jsd': cell_dict.get('cell_jsd', 0.0),
                                'cell_methylation_proportion': cell_dict.get('cell_methylation_proportion', 0.0),
                                'age': cell_dict.get('age', int(year_str))
                            }
                        cells.append(cell_data)
                    original_petri.cell_history[year_str] = cells
            
            # Extract gene_jsd_history
            gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(sim_data)
            if gene_jsd_hist:
                original_petri.gene_jsd_history = gene_jsd_hist
                original_petri.mean_gene_jsd_history = mean_hist
                original_petri.median_gene_jsd_history = median_hist
            
            # Create plotter and generate timeline plots
            plot_paths = PlotPaths(results_dir)
            plotter = PetriDishPlotter(original_petri)
            
            # JSD timeline
            jsd_timeline_path = plot_paths.get_cell_jsd_timeline_path()
            plotter.plot_jsd("Original Simulation JSD Timeline", jsd_timeline_path)
            print(f"  ✓ Cell JSD timeline")
            
            # Methylation timeline
            meth_timeline_path = plot_paths.get_cell_methylation_timeline_path()
            plotter.plot_cell_methylation_proportion(
                "Original Simulation Cell Methylation Proportion Timeline", meth_timeline_path
            )
            print(f"  ✓ Cell methylation timeline")
            
            # Gene JSD timeline
            if gene_jsd_hist:
                gene_jsd_timeline_path = plot_paths.get_gene_jsd_timeline_path()
                plotter.plot_gene_jsd_timeline(
                    "Original Simulation Gene JSD Timeline", gene_jsd_timeline_path
                )
                print(f"  ✓ Gene JSD timeline")
                
                # Gene JSD heatmap
                heatmap_path = os.path.join(results_dir, "simulation_gene_jsd_heatmap.png")
                plotter.plot_gene_jsd_heatmap(
                    title="Gene JSD Evolution Heatmap",
                    output_path=heatmap_path
                )
                print(f"  ✓ Gene JSD heatmap")
            else:
                print(f"  ⚠️  No gene JSD data in simulation")
        
    except Exception as e:
        print(f"  ⚠️  Could not generate timeline plots: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Analyze populations and generate all plots",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--base-dir", required=True, type=str,
                       help="Base directory with all data")
    parser.add_argument("--simulation", required=True, type=str,
                       help="Original simulation path for timeline plots")
    
    # Optional arguments
    parser.add_argument("--bins", type=int, default=200,
                       help="Number of bins for histograms")
    parser.add_argument("--max-gene-plots", type=int, default=None,
                       help="Maximum number of gene plots to generate")
    
    args = parser.parse_args()
    
    # Load metadata
    metadata = load_metadata(args.base_dir)
    gene_rate_groups = [tuple(g) for g in metadata['gene_rate_groups']]
    first_year = metadata['first_snapshot_year']
    second_year = metadata['second_snapshot_year']
    
    print("=" * 60)
    print("ANALYZE AND PLOT")
    print("=" * 60)
    print(f"Base directory: {args.base_dir}")
    print(f"Simulation: {args.simulation}")
    print(f"Gene rate groups: {gene_rate_groups}")
    print(f"Snapshots: year {first_year}, year {second_year}")
    print(f"Bins: {args.bins}")
    print("=" * 60)
    
    # Create results directory
    results_dir = os.path.join(args.base_dir, "results")
    os.makedirs(results_dir, exist_ok=True)
    
    # Initialize plot paths
    plot_paths = PlotPaths(results_dir)
    plot_paths.create_all_directories()
    
    # Determine compression from existing files
    snapshots_dir = os.path.join(args.base_dir, "snapshots")
    first_snapshot_gz = os.path.join(snapshots_dir, f"year{first_year}_snapshot.json.gz")
    use_compression = os.path.exists(first_snapshot_gz)
    ext = ".json.gz" if use_compression else ".json"
    
    # ========================================================================
    # Plot Snapshot Distributions
    # ========================================================================
    print("\n" + "=" * 60)
    print("STAGE 7a: Plot Snapshot Distributions")
    print("=" * 60)
    
    # Load and plot first snapshot
    first_snapshot_path = os.path.join(snapshots_dir, f"year{first_year}_snapshot{ext}")
    first_snapshot_cells = load_snapshot_cells(first_snapshot_path)
    plot_snapshot_distributions(
        first_snapshot_cells, first_year, args.bins, plot_paths, gene_rate_groups
    )
    
    # Load and plot second snapshot
    second_snapshot_path = os.path.join(snapshots_dir, f"year{second_year}_snapshot{ext}")
    second_snapshot_cells = load_snapshot_cells(second_snapshot_path)
    plot_snapshot_distributions(
        second_snapshot_cells, second_year, args.bins, plot_paths, gene_rate_groups
    )
    
    # ========================================================================
    # Plot Individual Trajectories
    # ========================================================================
    print("\n" + "=" * 60)
    print("STAGE 7b: Plot Individual Trajectories")
    print("=" * 60)
    
    plot_individual_trajectories(args.base_dir, plot_paths)
    
    # ========================================================================
    # Analyze Populations
    # ========================================================================
    print("\n" + "=" * 60)
    print("STAGE 7c: Analyze and Compare Populations")
    print("=" * 60)
    
    # Load all individuals
    individuals_dir = os.path.join(args.base_dir, "individuals")
    
    print("\n  Loading all individuals...")
    mutant_dishes = load_all_petri_dishes(os.path.join(individuals_dir, "mutant"))
    control1_dishes = load_all_petri_dishes(os.path.join(individuals_dir, "control1"))
    control2_dishes = load_all_petri_dishes(os.path.join(individuals_dir, "control2"))
    
    print(f"    Loaded {len(mutant_dishes)} mutant individuals")
    print(f"    Loaded {len(control1_dishes)} control1 individuals")
    print(f"    Loaded {len(control2_dishes)} control2 individuals")
    
    # Cell-level JSD analysis
    print("\n  Analyzing cell-level JSD...")
    analysis_results = analyze_populations_from_dishes(
        mutant_dishes, control1_dishes, control2_dishes, plot_paths
    )
    
    # Cell methylation analysis
    print("\n  Analyzing cell methylation proportions...")
    methylation_results = analyze_cell_methylation_proportion_comparison(
        mutant_dishes, control1_dishes, control2_dishes, plot_paths
    )
    
    # Gene-level JSD analysis
    print("\n  Analyzing gene-level JSD...")
    gene_analysis_results = generate_gene_jsd_analysis(
        mutant_dishes, control1_dishes, control2_dishes, plot_paths
    )
    
    # Generate gene JSD distribution plots
    if 'gene_jsd_analysis' in gene_analysis_results:
        print("\n  Generating gene JSD distribution plots...")
        plot_gene_jsd_distributions(
            gene_analysis_results['gene_jsd_analysis'], 
            plot_paths,
            max_genes=args.max_gene_plots
        )
        
        # Generate individual-averaged gene JSD comparison plot
        gene_jsd_path = plot_paths.get_gene_jsd_analysis_path()
        if os.path.exists(gene_jsd_path):
            plot_gene_jsd_individual_comparison(
                gene_jsd_path=gene_jsd_path,
                plot_paths=plot_paths,
                verbose=True
            )
    
    # Gene methylation analysis
    print("\n  Analyzing gene methylation proportions...")
    gene_methylation_results = generate_gene_methylation_analysis(
        mutant_dishes, control1_dishes, control2_dishes, plot_paths
    )
    
    # Generate gene methylation comparison plot
    gene_methylation_path = plot_paths.get_gene_methylation_analysis_path()
    if os.path.exists(gene_methylation_path):
        plot_gene_methylation_individual_comparison(
            gene_methylation_path=gene_methylation_path,
            plot_paths=plot_paths,
            verbose=True
        )
    
    # ========================================================================
    # Generate Timeline Plots
    # ========================================================================
    print("\n" + "=" * 60)
    print("STAGE 7d: Generate Timeline Plots")
    print("=" * 60)
    
    generate_timeline_plots(args.simulation, results_dir)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"All results saved to: {results_dir}")


if __name__ == "__main__":
    main()