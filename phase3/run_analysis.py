#!/usr/bin/env python3
"""
Phase 3: Analysis and Visualization Pipeline
Reads data from phase2 output and generates all plots and analysis.
"""

import argparse
import os
import sys
import json
import time
from datetime import datetime
from typing import Dict, List, Optional, Any

# Try to import YAML support
try:
    import yaml
except ImportError:
    yaml = None

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from core import (
    load_snapshot_cells, load_all_petri_dishes, load_petri_dish,
    extract_gene_jsd_from_history, load_phase2_metadata, smart_open,
    PlotPaths,
    plot_cell_jsd_distribution, plot_cell_methylation_proportion_histogram,
    plot_gene_jsd_distribution, plot_gene_methylation_proportion_histogram,
    analyze_cell_methylation_proportion_comparison, analyze_populations_from_dishes,
    plot_gene_jsd_distributions, plot_gene_jsd_individual_comparison,
    generate_gene_methylation_analysis, plot_gene_methylation_individual_comparison,
    generate_gene_jsd_analysis
)
from core.petri_dish_plotter import PetriDishPlotter


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    config = {}
    
    if yaml is None:
        return config
    
    # Load default config if it exists
    default_path = os.path.join(os.path.dirname(__file__), "configs", "default.yaml")
    if os.path.exists(default_path):
        try:
            with open(default_path, 'r') as f:
                default_config = yaml.safe_load(f) or {}
                config.update(default_config)
        except Exception as e:
            print(f"Warning: Could not load default config: {e}")
    
    # Load user config if provided
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f) or {}
                config.update(user_config)
        except Exception as e:
            print(f"Error loading config file {config_path}: {e}")
    
    return config


def generate_run_id(phase2_dir: str, bins: int) -> str:
    """Generate a unique run ID for the analysis."""
    # Extract info from phase2 directory name
    phase2_basename = os.path.basename(phase2_dir.rstrip('/'))
    
    # Add analysis-specific params
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    
    return f"analysis_bins{bins}_{timestamp}"


def plot_snapshot_distributions(
    cells: List, year: int, bins: int,
    plot_paths: PlotPaths, gene_rate_groups: List
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


def plot_individual_trajectories(individuals_dir: str, plot_paths: PlotPaths) -> None:
    """Plot growth trajectories for mutant and control1 individuals."""
    print("\nGenerating Individual Growth Trajectories...")
    
    # Plot mutant individuals
    mutant_dir = os.path.join(individuals_dir, "mutant")
    if os.path.exists(mutant_dir):
        print("\n  Creating plots for mutant individuals...")
        mutant_plot_count = 0
        
        for i, filename in enumerate(sorted(os.listdir(mutant_dir)), 1):
            if filename.endswith(('.json', '.json.gz')):
                filepath = os.path.join(mutant_dir, filename)
                try:
                    petri = load_petri_dish(filepath, include_cell_history=True)
                    
                    if hasattr(petri, 'cell_history') and petri.cell_history:
                        plotter = PetriDishPlotter(petri)
                        
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
    if os.path.exists(control1_dir):
        print("\n  Creating plots for control1 individuals...")
        control1_plot_count = 0
        
        for i, filename in enumerate(sorted(os.listdir(control1_dir)), 1):
            if filename.endswith(('.json', '.json.gz')):
                filepath = os.path.join(control1_dir, filename)
                try:
                    petri = load_petri_dish(filepath, include_cell_history=True)
                    
                    if hasattr(petri, 'cell_history') and petri.cell_history:
                        plotter = PetriDishPlotter(petri)
                        
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


def generate_timeline_plots(simulation_path: str, results_dir: str) -> None:
    """Generate timeline plots from original simulation."""
    print("\nGenerating Original Simulation Timeline Plots...")
    
    try:
        print(f"  Loading simulation with full history...")
        
        # Load simulation data
        with smart_open(simulation_path, 'r') as f:
            sim_data = json.load(f)
        
        if 'history' in sim_data and sim_data['history']:
            # Create PetriDish from simulation
            from cell import PetriDish
            
            params = sim_data.get('config', sim_data.get('parameters', {}))  # Support both old and new format
            rate = params.get('rate')
            gene_rate_groups = params.get('gene_rate_groups')
            n_sites = params.get('n', 1000)
            gene_size = params.get('gene_size', 5)
            
            if gene_rate_groups:
                rate_groups = [(g[0], g[1]) for g in gene_rate_groups]
                original_petri = PetriDish(gene_rate_groups=rate_groups, n=n_sites, gene_size=gene_size)
            else:
                original_petri = PetriDish(rate=rate, n=n_sites, gene_size=gene_size)
            
            # Build cell history
            original_petri.cell_history = {}
            original_petri.years_simulated = len(sim_data['history']) - 1
            
            for year_str, year_data in sim_data['history'].items():
                if 'cells' in year_data:
                    cells = []
                    for cell_dict in year_data['cells']:
                        if 'cpg_sites' in cell_dict:
                            cell_data = {
                                'methylated': cell_dict['cpg_sites'],
                                'cell_jsd': cell_dict.get('cell_jsd', 0.0),
                                'cell_methylation_proportion': sum(cell_dict['cpg_sites']) / len(cell_dict['cpg_sites']),
                                'age': int(year_str)
                            }
                        else:
                            cell_data = cell_dict
                        cells.append(cell_data)
                    original_petri.cell_history[year_str] = cells
            
            # Extract gene JSD history
            gene_jsd_hist = extract_gene_jsd_from_history(sim_data)
            if gene_jsd_hist:
                original_petri.gene_jsd_history = gene_jsd_hist
            
            # Generate plots
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
    
    except Exception as e:
        print(f"  ⚠️  Could not generate timeline plots: {e}")


def main():
    parser = argparse.ArgumentParser(
        description="Phase 3: Analyze populations and generate plots",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--phase2-dir", required=True, type=str,
                       help="Path to phase2 output directory")
    parser.add_argument("--simulation", required=True, type=str,
                       help="Path to original phase1 simulation for timeline plots")
    
    # Optional arguments
    parser.add_argument("--bins", type=int, default=200,
                       help="Number of bins for histograms")
    parser.add_argument("--max-gene-plots", type=int, default=None,
                       help="Maximum number of per-gene plots to generate")
    parser.add_argument("--output-dir", type=str, default=None,
                       help="Output directory (default: auto-generate)")
    parser.add_argument("--config", type=str,
                       help="Path to configuration file")
    
    args = parser.parse_args()
    
    # Load config
    config = load_config(args.config)
    for key, value in config.items():
        if not hasattr(args, key) or getattr(args, key) is None:
            setattr(args, key, value)
    
    # Validate inputs
    if not os.path.exists(args.phase2_dir):
        print(f"Error: Phase2 directory not found: {args.phase2_dir}")
        sys.exit(1)
    
    if not os.path.exists(args.simulation):
        print(f"Error: Simulation file not found: {args.simulation}")
        sys.exit(1)
    
    # Determine output directory
    if args.output_dir:
        results_base = args.output_dir
    else:
        run_id = generate_run_id(args.phase2_dir, args.bins)
        results_base = os.path.join("data", run_id)
    
    results_dir = os.path.join(results_base, "results")
    
    # Initialize plot paths
    plot_paths = PlotPaths(results_dir)
    plot_paths.create_all_directories()
    
    print("=" * 60)
    print("PHASE 3: ANALYSIS AND VISUALIZATION")
    print("=" * 60)
    print(f"Phase2 data: {args.phase2_dir}")
    print(f"Simulation: {args.simulation}")
    print(f"Output: {results_dir}")
    print(f"Bins: {args.bins}")
    if args.max_gene_plots:
        print(f"Max gene plots: {args.max_gene_plots}")
    print("=" * 60)
    
    start_time = time.time()
    
    # Load metadata
    metadata = load_phase2_metadata(args.phase2_dir)
    
    if 'snapshots' in metadata:
        gene_rate_groups = [tuple(g) for g in metadata['snapshots']['gene_rate_groups']]
        first_year = metadata['snapshots']['first_snapshot_year']
        second_year = metadata['snapshots']['second_snapshot_year']
        print(f"\nDetected parameters:")
        print(f"  Gene rate groups: {gene_rate_groups}")
        print(f"  Snapshots: year {first_year}, year {second_year}")
    else:
        print("Warning: Could not load snapshots metadata")
        gene_rate_groups = None
        first_year = 30
        second_year = 50
    
    # Determine file extension
    snapshots_dir = os.path.join(args.phase2_dir, "snapshots")
    use_gz = os.path.exists(os.path.join(snapshots_dir, f"year{first_year}_snapshot.json.gz"))
    ext = ".json.gz" if use_gz else ".json"
    
    # ========================================================================
    # Stage 1: Plot Snapshot Distributions
    # ========================================================================
    print("\n" + "=" * 60)
    print("Stage 1: Plot Snapshot Distributions")
    print("=" * 60)
    
    # First snapshot
    first_snapshot_path = os.path.join(snapshots_dir, f"year{first_year}_snapshot{ext}")
    if os.path.exists(first_snapshot_path):
        first_snapshot_cells = load_snapshot_cells(first_snapshot_path)
        plot_snapshot_distributions(
            first_snapshot_cells, first_year, args.bins, plot_paths, gene_rate_groups
        )
    
    # Second snapshot
    second_snapshot_path = os.path.join(snapshots_dir, f"year{second_year}_snapshot{ext}")
    if os.path.exists(second_snapshot_path):
        second_snapshot_cells = load_snapshot_cells(second_snapshot_path)
        plot_snapshot_distributions(
            second_snapshot_cells, second_year, args.bins, plot_paths, gene_rate_groups
        )
    
    # ========================================================================
    # Stage 2: Plot Individual Trajectories
    # ========================================================================
    print("\n" + "=" * 60)
    print("Stage 2: Plot Individual Trajectories")
    print("=" * 60)
    
    individuals_dir = os.path.join(args.phase2_dir, "individuals")
    plot_individual_trajectories(individuals_dir, plot_paths)
    
    # ========================================================================
    # Stage 3: Analyze and Compare Populations
    # ========================================================================
    print("\n" + "=" * 60)
    print("Stage 3: Analyze and Compare Populations")
    print("=" * 60)
    
    # Load all individuals
    print("\n  Loading all individuals...")
    mutant_dishes = load_all_petri_dishes(os.path.join(individuals_dir, "mutant"))
    control1_dishes = load_all_petri_dishes(os.path.join(individuals_dir, "control1"))
    control2_dishes = load_all_petri_dishes(os.path.join(individuals_dir, "control2"))
    
    print(f"    Loaded {len(mutant_dishes)} mutant individuals")
    print(f"    Loaded {len(control1_dishes)} control1 individuals")
    print(f"    Loaded {len(control2_dishes)} control2 individuals")
    
    # Cell-level analysis
    print("\n  Analyzing cell-level JSD...")
    analyze_populations_from_dishes(
        mutant_dishes, control1_dishes, control2_dishes, plot_paths
    )
    
    # Cell methylation analysis
    print("\n  Analyzing cell methylation proportions...")
    analyze_cell_methylation_proportion_comparison(
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
        
        # Generate comparison plot
        gene_jsd_path = plot_paths.get_gene_jsd_analysis_path()
        if os.path.exists(gene_jsd_path):
            plot_gene_jsd_individual_comparison(
                gene_jsd_path=gene_jsd_path,
                plot_paths=plot_paths,
                verbose=True
            )
    
    # Gene methylation analysis
    print("\n  Analyzing gene methylation proportions...")
    generate_gene_methylation_analysis(
        mutant_dishes, control1_dishes, control2_dishes, plot_paths
    )
    
    # Generate comparison plot
    gene_methylation_path = plot_paths.get_gene_methylation_analysis_path()
    if os.path.exists(gene_methylation_path):
        plot_gene_methylation_individual_comparison(
            gene_methylation_path=gene_methylation_path,
            plot_paths=plot_paths,
            verbose=True
        )
    
    # ========================================================================
    # Stage 4: Generate Timeline Plots
    # ========================================================================
    print("\n" + "=" * 60)
    print("Stage 4: Generate Timeline Plots")
    print("=" * 60)
    
    generate_timeline_plots(args.simulation, results_dir)
    
    # Summary
    elapsed_time = time.time() - start_time
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Total time: {elapsed_time:.2f} seconds")
    print(f"Results saved to: {results_dir}")


if __name__ == "__main__":
    main()