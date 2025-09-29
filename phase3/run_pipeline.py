#!/usr/bin/env python3
"""
Simplified pipeline driver that directly orchestrates data extraction and plotting.
Uses timeline CSVs as the single source of truth for histogram generation.
"""

import os
import sys
import glob
import json
import gzip
import csv
import argparse
import subprocess
import numpy as np
import pandas as pd
from datetime import datetime
from typing import Dict, List, Optional, Tuple, Any

# Import plotting functions
from plot_histogram_original import (
    plot_cell_methylation_histogram_original,
    plot_cell_jsd_histogram_original,
    plot_histogram_original
)


def load_snapshot_metadata(phase2_dir: str) -> Dict[str, Any]:
    """
    Load snapshot metadata to get snapshot years and gene configuration.

    Returns dict with keys:
    - first_snapshot_year: int
    - second_snapshot_year: int
    - gene_rate_groups: list
    - n_sites: int
    - gene_size: int
    """
    metadata_path = os.path.join(phase2_dir, 'snapshots', 'metadata.json')
    if not os.path.exists(metadata_path):
        raise FileNotFoundError(f"Snapshot metadata not found: {metadata_path}")

    with open(metadata_path, 'r') as f:
        return json.load(f)


def extract_snapshot_from_timeline(timeline_csv: str, target_year: int, column_prefix: str) -> np.ndarray:
    """
    Extract values for a specific year from timeline CSV.

    Args:
        timeline_csv: Path to timeline CSV file
        target_year: Year to extract
        column_prefix: 'value' for cells, 'gene' for genes

    Returns:
        Array of values for that year
    """
    df = pd.read_csv(timeline_csv)
    year_row = df[df['year'] == target_year]

    if year_row.empty:
        raise ValueError(f"Year {target_year} not found in {os.path.basename(timeline_csv)}")

    # Extract columns matching pattern (value_0, value_1, ... or gene_0, gene_1, ...)
    value_columns = [col for col in df.columns if col.startswith(f'{column_prefix}_')]
    if not value_columns:
        raise ValueError(f"No columns with prefix '{column_prefix}_' found in {os.path.basename(timeline_csv)}")

    values = year_row[value_columns].iloc[0].dropna().values

    return values


def parse_gene_rates(gene_rates_str: str) -> Optional[List[Tuple[int, float]]]:
    """
    Parse gene rate groups from string format.
    """
    if not gene_rates_str:
        return None

    try:
        parts = gene_rates_str.split(',')
        gene_rate_groups = []
        for part in parts:
            count, rate = part.split(':')
            gene_rate_groups.append((int(count), float(rate)))
        return gene_rate_groups
    except:
        return None


def main():
    """
    Main pipeline execution.
    """
    parser = argparse.ArgumentParser(
        description='Run analysis pipeline using timeline CSVs for histogram generation',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument('--phase2-dir', required=True,
                       help='Path to phase2 output directory containing snapshots/')

    # Optional arguments
    parser.add_argument('--skip-timeline', action='store_true',
                       help='Skip timeline extraction from Phase 1 simulation')
    parser.add_argument('--results-dir', default=None,
                       help='Base directory for results (default: phase2-dir/results)')
    parser.add_argument('--gene-rates',
                       default="5:0.004,5:0.005,5:0.006,5:0.007",
                       help='Gene rate groups (e.g., "5:0.004,5:0.005,5:0.006,5:0.007")')
    parser.add_argument('--skip-plots', action='store_true',
                       help='Only extract CSVs, skip plotting')
    parser.add_argument('--skip-comparison', action='store_true',
                       help='Skip batch comparison extraction and plotting')
    parser.add_argument('--skip-gene-comparison', action='store_true',
                       help='Skip per-gene comparison plots (40 plots for individual genes)')
    parser.add_argument('--verbose', action='store_true', default=True,
                       help='Print detailed progress messages')

    args = parser.parse_args()

    # Setup directories
    if args.results_dir is None:
        results_dir = os.path.join(os.path.abspath(args.phase2_dir), 'results')
    else:
        results_dir = args.results_dir

    tables_dir = os.path.join(results_dir, 'tables')
    plots_dir = os.path.join(results_dir, 'plots')
    os.makedirs(tables_dir, exist_ok=True)
    os.makedirs(plots_dir, exist_ok=True)

    # Parse gene rates from command line (may be overridden by metadata)
    gene_rate_groups_cli = parse_gene_rates(args.gene_rates)

    print("=" * 60)
    print("ANALYSIS PIPELINE (TIMELINE-BASED)")
    print("=" * 60)
    print(f"Phase2 directory: {args.phase2_dir}")
    print(f"Tables directory: {tables_dir}")
    print(f"Plots directory: {plots_dir}")
    print("=" * 60)

    try:
        # ====================================================================
        # Stage 1: Extract Timeline from Simulation (automatic)
        # ====================================================================
        # Auto-detect Phase 1 simulation from phase2-dir structure
        # Phase2 dir pattern: .../phase1_dir/phase2_subdir/
        # Simulation is at: .../phase1_dir/simulation.json[.gz]

        simulation_path = None
        if not args.skip_timeline:
            # Try to find the Phase 1 simulation
            phase2_abs = os.path.abspath(args.phase2_dir)
            phase1_dir = os.path.dirname(phase2_abs)

            # Check for simulation.json or simulation.json.gz
            for sim_file in ['simulation.json', 'simulation.json.gz']:
                potential_path = os.path.join(phase1_dir, sim_file)
                if os.path.exists(potential_path):
                    simulation_path = potential_path
                    break

            if simulation_path:
                print("\nStage 1: Extracting timeline from Phase 1 simulation...")
                print("-" * 60)
                print(f"  Auto-detected simulation: {os.path.basename(phase1_dir)}/{os.path.basename(simulation_path)}")

                # Run extract_simulation_timeline.py as subprocess
                timeline_script = os.path.join(os.path.dirname(__file__), 'extract_simulation_timeline.py')
                timeline_cmd = [
                    sys.executable, timeline_script,
                    '--simulation', simulation_path,
                    '--output-dir', tables_dir
                ]

                if args.verbose:
                    print(f"  Extracting timeline data to {tables_dir}/")

                try:
                    result = subprocess.run(timeline_cmd, capture_output=True, text=True)
                    if result.returncode == 0:
                        print("  ✓ Timeline extraction completed")
                        # List created files
                        timeline_files = ['simulation_cell_jsd.csv', 'simulation_cell_methylation.csv',
                                        'simulation_gene_jsd.csv', 'simulation_gene_methylation.csv']
                        for f in timeline_files:
                            fpath = os.path.join(tables_dir, f)
                            if os.path.exists(fpath):
                                size_kb = os.path.getsize(fpath) / 1024
                                print(f"    - {f} ({size_kb:.1f} KB)")
                    else:
                        print(f"  ⚠️  Timeline extraction failed: {result.stderr}")
                except Exception as e:
                    print(f"  ⚠️  Error running timeline extraction: {str(e)}")
            else:
                print("\nStage 1: No Phase 1 simulation found (skipping timeline extraction)")
                print("  Looked in:", phase1_dir)
        else:
            print("\nStage 1: Skipping timeline extraction (--skip-timeline)")

        # ====================================================================
        # Stage 2: Generate histogram plots from timeline CSVs
        # ====================================================================
        if not args.skip_plots:
            print("\nStage 2: Generating histogram plots from timeline data...")
            print("-" * 60)

            # Load snapshot metadata to get years and gene configuration
            try:
                metadata = load_snapshot_metadata(args.phase2_dir)
                snapshot_years = [metadata['first_snapshot_year'], metadata['second_snapshot_year']]
                gene_rate_groups = metadata.get('gene_rate_groups', gene_rate_groups_cli)
                print(f"  Snapshot years from metadata: {snapshot_years}")
                if gene_rate_groups:
                    print(f"  Gene rate groups: {gene_rate_groups}")
            except FileNotFoundError as e:
                print(f"  ⚠️  {str(e)}")
                print("  Cannot generate histograms without snapshot metadata")
                snapshot_years = []
                gene_rate_groups = gene_rate_groups_cli

            if snapshot_years:
                # Define timeline CSV paths
                timeline_files = {
                    'cell_jsd': os.path.join(tables_dir, 'simulation_cell_jsd.csv'),
                    'cell_methylation': os.path.join(tables_dir, 'simulation_cell_methylation.csv'),
                    'gene_jsd': os.path.join(tables_dir, 'simulation_gene_jsd.csv'),
                    'gene_methylation': os.path.join(tables_dir, 'simulation_gene_methylation.csv')
                }

                # Check if timeline CSVs exist
                missing_files = [name for name, path in timeline_files.items() if not os.path.exists(path)]
                if missing_files:
                    print(f"  ⚠️  Missing timeline CSVs: {missing_files}")
                    print("  Run Stage 1 (timeline extraction) first or use --skip-plots")
                else:
                    plot_count = 0

                    for year in snapshot_years:
                        print(f"\n  Processing year {year}...")

                        try:
                            # Extract cell data for this year
                            cell_jsd_values = extract_snapshot_from_timeline(
                                timeline_files['cell_jsd'], year, 'value'
                            )
                            cell_meth_values = extract_snapshot_from_timeline(
                                timeline_files['cell_methylation'], year, 'value'
                            )

                            # Verify both arrays have same length (same cells)
                            if len(cell_jsd_values) != len(cell_meth_values):
                                print(f"    ⚠️  Mismatch: {len(cell_jsd_values)} JSD values vs {len(cell_meth_values)} methylation values")
                                continue

                            print(f"    Found {len(cell_jsd_values)} cells")

                            # Plot cell JSD histogram
                            jsd_output = os.path.join(plots_dir, f'year{year}_cell_jsd_histogram.png')
                            plot_cell_jsd_histogram_original(
                                data=cell_jsd_values,
                                output_path=jsd_output,
                                year=year,
                                n_cells=len(cell_jsd_values),
                                gene_rate_groups=gene_rate_groups
                            )
                            plot_count += 1

                            # Plot cell methylation histogram
                            meth_output = os.path.join(plots_dir, f'year{year}_cell_methylation_histogram.png')
                            plot_cell_methylation_histogram_original(
                                data=cell_meth_values,
                                output_path=meth_output,
                                year=year,
                                n_cells=len(cell_meth_values),
                                gene_rate_groups=gene_rate_groups
                            )
                            plot_count += 1

                            # Extract and plot gene data if available
                            try:
                                gene_jsd_values = extract_snapshot_from_timeline(
                                    timeline_files['gene_jsd'], year, 'gene'
                                )
                                gene_meth_values = extract_snapshot_from_timeline(
                                    timeline_files['gene_methylation'], year, 'gene'
                                )

                                print(f"    Found {len(gene_jsd_values)} genes")

                                # Plot gene JSD histogram
                                gene_jsd_output = os.path.join(plots_dir, f'year{year}_gene_jsd_histogram.png')
                                plot_histogram_original(
                                    data=gene_jsd_values,
                                    title="Gene JSD",
                                    xlabel="Gene JSD",
                                    ylabel="Count",
                                    histogram_color="#1f77b4",  # Blue for JSD
                                    mean_line_color="red",
                                    bins=20,  # Only 20 genes
                                    output_path=gene_jsd_output,
                                    year=year,
                                    n_cells=len(gene_jsd_values),
                                    gene_rate_groups=gene_rate_groups,
                                    plot_type="jsd"
                                )
                                plot_count += 1

                                # Plot gene methylation histogram
                                gene_meth_output = os.path.join(plots_dir, f'year{year}_gene_methylation_histogram.png')
                                plot_histogram_original(
                                    data=gene_meth_values,
                                    title="Gene Methylation Proportion",
                                    xlabel="Gene Methylation Proportion",
                                    ylabel="Count",
                                    histogram_color="#d62728",  # Red for methylation
                                    mean_line_color="darkblue",
                                    bins=20,  # Only 20 genes
                                    output_path=gene_meth_output,
                                    year=year,
                                    n_cells=len(gene_meth_values),
                                    gene_rate_groups=gene_rate_groups,
                                    plot_type="methylation"
                                )
                                plot_count += 1

                            except Exception as e:
                                print(f"    Note: Gene data not available for year {year}: {str(e)}")

                            print(f"    ✓ Generated plots for year {year}")

                        except Exception as e:
                            print(f"    ⚠️  Error processing year {year}: {str(e)}")
                            continue

                    if plot_count > 0:
                        print(f"\n✓ Generated {plot_count} histogram plots total")
        else:
            print("\nStage 2: Skipping histogram plots (--skip-plots)")

        # ====================================================================
        # Stage 3: Generate Timeline Plots from Stage 1 CSVs
        # ====================================================================
        if not args.skip_plots and not args.skip_timeline:
            # Check if timeline CSVs exist from Stage 1
            timeline_csv_files = [
                os.path.join(tables_dir, 'simulation_cell_jsd.csv'),
                os.path.join(tables_dir, 'simulation_cell_methylation.csv'),
                os.path.join(tables_dir, 'simulation_gene_jsd.csv'),
                os.path.join(tables_dir, 'simulation_gene_methylation.csv')
            ]

            if any(os.path.exists(f) for f in timeline_csv_files):
                print("\nStage 3: Generating timeline plots from simulation CSVs...")
                print("-" * 60)

                from plot_simulation_timeline import plot_all_timelines

                try:
                    # Generate timeline plots
                    plot_all_timelines(
                        tables_dir=tables_dir,
                        output_dir=plots_dir,
                        title_prefix=None,  # Could be customized
                        plots_to_generate=None  # Generate all available
                    )
                    print("  ✓ Timeline plots generated successfully")
                except Exception as e:
                    print(f"  ⚠️  Error generating timeline plots: {str(e)}")
            else:
                if args.verbose:
                    print("\nStage 3: No timeline CSVs found (Stage 1 may have been skipped)")

        # ====================================================================
        # Stage 4: Extract Batch Comparison Data
        # ====================================================================
        comparison_cell_csv_path = None
        comparison_gene_csv_path = None
        if not args.skip_comparison:
            print("\n" + "=" * 60)
            print("Stage 4: Extract Batch Comparison Data")
            print("=" * 60)

            # Check if individuals directory exists
            individuals_dir = os.path.join(args.phase2_dir, "individuals")

            if os.path.exists(individuals_dir):
                from extract_batch_comparison import extract_cell_comparison, extract_gene_comparison

                comparison_cell_csv_path = os.path.join(tables_dir, "batch_comparison_cell.csv")
                comparison_gene_csv_path = os.path.join(tables_dir, "batch_comparison_gene.csv")

                # Extract cell-level metrics
                try:
                    extract_cell_comparison(
                        individuals_dir=individuals_dir,
                        output_csv=comparison_cell_csv_path,
                        verbose=args.verbose
                    )
                    print(f"✓ Saved cell comparison data to {os.path.basename(comparison_cell_csv_path)}")
                except Exception as e:
                    print(f"⚠️ Error extracting cell comparison: {str(e)}")
                    comparison_cell_csv_path = None

                # Extract gene-level metrics
                try:
                    extract_gene_comparison(
                        individuals_dir=individuals_dir,
                        output_csv=comparison_gene_csv_path,
                        verbose=args.verbose
                    )
                    print(f"✓ Saved gene comparison data to {os.path.basename(comparison_gene_csv_path)}")
                except Exception as e:
                    print(f"⚠️ Error extracting gene comparison: {str(e)}")
                    comparison_gene_csv_path = None
            else:
                print(f"⚠️ Individuals directory not found: {individuals_dir}")
                print("   Skipping batch comparison extraction")

        # ====================================================================
        # Stage 5: Generate Comparison Plots
        # ====================================================================
        comparison_plot_count = 0
        if not args.skip_comparison and not args.skip_plots and (comparison_cell_csv_path or comparison_gene_csv_path):
            print("\n" + "=" * 60)
            print("Stage 5: Generate Comparison Plots")
            print("=" * 60)

            from plot_comparison_generic import plot_comparison_generic

            # Define cell metrics to plot
            cell_metrics = [
                {
                    'column': 'cell_jsd_mean',
                    'title': 'Mean Cell JSD per Individual Across Batches',
                    'ylabel': 'Mean Cell JSD (averaged over ~288 cells)',
                    'filename': 'cell_jsd_mean_comparison.png'
                },
                {
                    'column': 'cell_methylation_mean',
                    'title': 'Mean Cell Methylation per Individual Across Batches',
                    'ylabel': 'Mean Cell Methylation (averaged over ~288 cells)',
                    'filename': 'cell_methylation_mean_comparison.png'
                }
            ]

            # Define gene metrics to plot
            gene_metrics = [
                {
                    'column': 'gene_jsd',
                    'title': 'Mean Gene JSD per Individual Across Batches',
                    'ylabel': 'Mean Gene JSD (averaged over 20 genes)',
                    'filename': 'gene_jsd_mean_comparison.png'
                },
                {
                    'column': 'gene_methylation',
                    'title': 'Mean Gene Methylation per Individual Across Batches',
                    'ylabel': 'Mean Gene Methylation (averaged over 20 genes)',
                    'filename': 'gene_methylation_mean_comparison.png'
                }
            ]

            # Generate cell comparison plots
            if comparison_cell_csv_path:
                for metric in cell_metrics:
                    try:
                        plot_comparison_generic(
                            csv_path=comparison_cell_csv_path,
                            value_column=metric['column'],
                            title=metric['title'],
                            ylabel=metric['ylabel'],
                            output_path=os.path.join(plots_dir, metric['filename']),
                            verbose=args.verbose,
                            y_range_padding=0.1  # 10% padding to prevent edge clipping
                        )
                        comparison_plot_count += 1
                    except Exception as e:
                        print(f"  ⚠️ Error creating {metric['title']}: {str(e)}")

            # Generate gene comparison plots (mean)
            if comparison_gene_csv_path:
                for metric in gene_metrics:
                    try:
                        plot_comparison_generic(
                            csv_path=comparison_gene_csv_path,
                            value_column=metric['column'],
                            title=metric['title'],
                            ylabel=metric['ylabel'],
                            output_path=os.path.join(plots_dir, metric['filename']),
                            verbose=args.verbose,
                            aggregation='mean',  # Explicitly specify mean aggregation
                            y_range_padding=0.1  # 10% padding to prevent edge clipping
                        )
                        comparison_plot_count += 1
                    except Exception as e:
                        print(f"  ⚠️ Error creating {metric['title']}: {str(e)}")

                # NEW: Generate gene standard deviation plots
                gene_std_metrics = [
                    {
                        'column': 'gene_jsd',
                        'title': 'Gene JSD Standard Deviation per Individual Across Batches',
                        'ylabel': 'Std Dev of Gene JSD (across 20 genes)',
                        'filename': 'gene_jsd_std_comparison.png'
                    },
                    {
                        'column': 'gene_methylation',
                        'title': 'Gene Methylation Standard Deviation per Individual Across Batches',
                        'ylabel': 'Std Dev of Gene Methylation (across 20 genes)',
                        'filename': 'gene_methylation_std_comparison.png'
                    }
                ]

                for metric in gene_std_metrics:
                    try:
                        plot_comparison_generic(
                            csv_path=comparison_gene_csv_path,
                            value_column=metric['column'],
                            title=metric['title'],
                            ylabel=metric['ylabel'],
                            output_path=os.path.join(plots_dir, metric['filename']),
                            verbose=args.verbose,
                            aggregation='std',  # Use std aggregation
                            y_range_padding=0.1  # 10% padding to prevent edge clipping
                        )
                        comparison_plot_count += 1
                    except Exception as e:
                        print(f"  ⚠️ Error creating {metric['title']}: {str(e)}")

            print(f"\n✓ Generated {comparison_plot_count} comparison plots")

        # ====================================================================
        # Stage 6: Generate Per-Gene Comparison Plots
        # ====================================================================
        per_gene_plot_count = 0
        if not args.skip_comparison and not args.skip_gene_comparison and not args.skip_plots and comparison_gene_csv_path:
            print("\n" + "=" * 60)
            print("Stage 6: Generate Per-Gene Comparison Plots")
            print("=" * 60)

            from plot_comparison_by_gene import plot_comparison_by_gene

            # Create output directory for per-gene plots
            per_gene_dir = os.path.join(plots_dir, 'comparison_by_gene')

            try:
                per_gene_plot_count = plot_comparison_by_gene(
                    gene_csv_path=comparison_gene_csv_path,
                    output_dir=per_gene_dir,
                    verbose=args.verbose
                )
                print(f"✓ Generated {per_gene_plot_count} per-gene comparison plots")
            except Exception as e:
                print(f"⚠️ Error creating per-gene comparison plots: {str(e)}")

        # ====================================================================
        # Complete
        # ====================================================================
        print("\n" + "=" * 60)
        print("PIPELINE COMPLETE")
        print("=" * 60)
        print(f"Results directory: {results_dir}")
        print(f"  Tables: {tables_dir}")
        print(f"  Plots: {plots_dir}")

        if not args.skip_plots:
            total_plots = plot_count + comparison_plot_count + per_gene_plot_count
            if total_plots > 0:
                print(f"\nGenerated {total_plots} total plots:")
                if plot_count > 0:
                    print(f"  - Histogram plots: {plot_count}")
                if comparison_plot_count > 0:
                    print(f"  - Comparison plots: {comparison_plot_count}")
                if per_gene_plot_count > 0:
                    print(f"  - Per-gene plots: {per_gene_plot_count}")
        print("=" * 60)

    except Exception as e:
        print(f"\n✗ Pipeline failed: {str(e)}")
        import traceback
        if args.verbose:
            traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()