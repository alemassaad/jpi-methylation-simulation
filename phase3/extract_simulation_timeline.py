#!/usr/bin/env python3
"""
Extract full timeline data from Phase 1 simulation into separate CSV files.
Creates 4 CSVs: cell JSD, cell methylation, gene JSD, and gene methylation.
Each row is a year, columns are individual values (value_* for cells, gene_* for genes).
"""

import os
import sys
import json
import csv
import argparse
import numpy as np
from typing import Dict, List, Any, Optional, Tuple

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from data_extractors import smart_open, dict_to_cell


def extract_cell_jsd_timeline(history: Dict[str, Any], start_year: int = 0,
                              end_year: Optional[int] = None) -> Tuple[Dict[int, List[float]], int]:
    """
    Extract cell JSD values for each year.

    Args:
        history: Dictionary with year keys (as strings) containing cell data
        start_year: First year to extract (inclusive)
        end_year: Last year to extract (inclusive), None for all

    Returns:
        Tuple of (data_dict, max_cells) where:
        - data_dict: {year: [jsd_values]}
        - max_cells: Maximum number of cells across all years
    """
    data = {}
    max_cells = 0

    # Determine year range
    available_years = sorted([int(y) for y in history.keys()])
    if end_year is None:
        end_year = max(available_years)

    for year in range(start_year, end_year + 1):
        year_str = str(year)
        if year_str not in history:
            continue

        year_data = history[year_str]
        cells = year_data.get('cells', [])

        # Extract JSD values
        jsd_values = []
        for cell_dict in cells:
            if 'cell_jsd' in cell_dict:
                jsd_values.append(cell_dict['cell_jsd'])
            else:
                # Calculate JSD if not present (shouldn't happen in modern format)
                cell = dict_to_cell(cell_dict)
                if hasattr(cell, 'cell_jsd'):
                    jsd_values.append(cell.cell_jsd)
                else:
                    jsd_values.append(0.0)

        data[year] = jsd_values
        max_cells = max(max_cells, len(jsd_values))

    return data, max_cells


def extract_cell_methylation_timeline(history: Dict[str, Any], start_year: int = 0,
                                     end_year: Optional[int] = None) -> Tuple[Dict[int, List[float]], int]:
    """
    Extract cell methylation proportions for each year.

    Args:
        history: Dictionary with year keys (as strings) containing cell data
        start_year: First year to extract (inclusive)
        end_year: Last year to extract (inclusive), None for all

    Returns:
        Tuple of (data_dict, max_cells) where:
        - data_dict: {year: [methylation_proportions]}
        - max_cells: Maximum number of cells across all years
    """
    data = {}
    max_cells = 0

    # Determine year range
    available_years = sorted([int(y) for y in history.keys()])
    if end_year is None:
        end_year = max(available_years)

    for year in range(start_year, end_year + 1):
        year_str = str(year)
        if year_str not in history:
            continue

        year_data = history[year_str]
        cells = year_data.get('cells', [])

        # Extract methylation proportions
        meth_values = []
        for cell_dict in cells:
            if 'cell_methylation_proportion' in cell_dict:
                meth_values.append(cell_dict['cell_methylation_proportion'])
            elif 'cpg_sites' in cell_dict:
                # Calculate from cpg_sites
                sites = cell_dict['cpg_sites']
                proportion = sum(sites) / len(sites) if sites else 0.0
                meth_values.append(proportion)
            else:
                meth_values.append(0.0)

        data[year] = meth_values
        max_cells = max(max_cells, len(meth_values))

    return data, max_cells


def extract_gene_jsd_timeline(history: Dict[str, Any], start_year: int = 0,
                             end_year: Optional[int] = None) -> Tuple[Dict[int, List[float]], int]:
    """
    Extract gene JSD values for each year.

    Args:
        history: Dictionary with year keys (as strings) containing gene JSD data
        start_year: First year to extract (inclusive)
        end_year: Last year to extract (inclusive), None for all

    Returns:
        Tuple of (data_dict, n_genes) where:
        - data_dict: {year: [gene_jsd_values]}
        - n_genes: Number of genes (should be consistent, typically 20)
    """
    data = {}
    n_genes = 0

    # Determine year range
    available_years = sorted([int(y) for y in history.keys()])
    if end_year is None:
        end_year = max(available_years)

    for year in range(start_year, end_year + 1):
        year_str = str(year)
        if year_str not in history:
            continue

        year_data = history[year_str]

        # Gene JSD should be directly available
        if 'gene_jsd' in year_data:
            gene_jsds = year_data['gene_jsd']
            data[year] = gene_jsds
            n_genes = max(n_genes, len(gene_jsds))
        else:
            # If not available, we'd need to calculate from cells
            # This requires creating a PetriDish and calculating
            cells = year_data.get('cells', [])
            if cells:
                # Create Cell objects
                cell_objects = [dict_to_cell(cell_dict) for cell_dict in cells]
                # Get parameters from first cell
                if cell_objects:
                    first_cell = cell_objects[0]
                    n = first_cell.n
                    gene_size = first_cell.gene_size if hasattr(first_cell, 'gene_size') else 5
                    gene_rate_groups = first_cell.gene_rate_groups if hasattr(first_cell, 'gene_rate_groups') else None

                    # Create PetriDish and calculate gene JSD
                    petri = PetriDish(n=n, gene_size=gene_size, gene_rate_groups=gene_rate_groups, cells=cell_objects)
                    gene_jsds = petri.calculate_gene_jsd()
                    data[year] = gene_jsds
                    n_genes = max(n_genes, len(gene_jsds))

    return data, n_genes


def extract_gene_methylation_timeline(history: Dict[str, Any], config: Dict[str, Any],
                                     start_year: int = 0, end_year: Optional[int] = None) -> Tuple[Dict[int, List[float]], int]:
    """
    Extract gene methylation proportions for each year.

    Args:
        history: Dictionary with year keys (as strings) containing cell data
        config: Configuration dictionary with gene_size info
        start_year: First year to extract (inclusive)
        end_year: Last year to extract (inclusive), None for all

    Returns:
        Tuple of (data_dict, n_genes) where:
        - data_dict: {year: [gene_methylation_proportions]}
        - n_genes: Number of genes (should be consistent, typically 20)
    """
    data = {}
    n_genes = 0
    gene_size = config.get('gene_size', 5)

    # Determine year range
    available_years = sorted([int(y) for y in history.keys()])
    if end_year is None:
        end_year = max(available_years)

    for year in range(start_year, end_year + 1):
        year_str = str(year)
        if year_str not in history:
            continue

        year_data = history[year_str]
        cells = year_data.get('cells', [])

        if not cells:
            continue

        # Need to calculate gene methylation from cells
        cell_objects = [dict_to_cell(cell_dict) for cell_dict in cells]

        # Get parameters from first cell
        first_cell = cell_objects[0]
        n = first_cell.n
        gene_size = first_cell.gene_size if hasattr(first_cell, 'gene_size') else gene_size
        gene_rate_groups = first_cell.gene_rate_groups if hasattr(first_cell, 'gene_rate_groups') else None

        # Create PetriDish to access gene_methylation_proportions
        petri = PetriDish(n=n, gene_size=gene_size, gene_rate_groups=gene_rate_groups, cells=cell_objects)

        # Get gene methylation proportions
        gene_meth_props = petri.gene_methylation_proportions
        data[year] = gene_meth_props
        n_genes = max(n_genes, len(gene_meth_props))

    return data, n_genes


def save_to_csv(data_dict: Dict[int, List[float]], output_path: str,
                max_values: int, column_prefix: str = "value") -> None:
    """
    Save timeline data to CSV file.

    Args:
        data_dict: {year: [values]} dictionary
        output_path: Path for output CSV file
        max_values: Maximum number of values (determines number of columns)
        column_prefix: Prefix for column names ("value" for cells, "gene" for genes)
    """
    # Create header
    header = ['year'] + [f"{column_prefix}_{i}" for i in range(max_values)]

    # Sort years
    years = sorted(data_dict.keys())

    # Write CSV
    with open(output_path, 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(header)

        for year in years:
            values = data_dict[year]
            # Pad with empty strings if needed
            row_values = values + [''] * (max_values - len(values))
            row = [year] + row_values
            writer.writerow(row)

    print(f"  Saved {len(years)} years × {max_values} {column_prefix}s to {output_path}")


def main():
    """Main entry point for timeline extraction."""
    parser = argparse.ArgumentParser(
        description="Extract full timeline data from Phase 1 simulation to CSV files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--simulation', type=str, required=True,
                       help='Path to Phase 1 simulation.json[.gz] file')
    parser.add_argument('--output-dir', type=str, default='tables',
                       help='Directory for output CSV files')
    parser.add_argument('--start-year', type=int, default=0,
                       help='First year to extract (inclusive)')
    parser.add_argument('--end-year', type=int, default=None,
                       help='Last year to extract (inclusive), default is all')
    parser.add_argument('--metrics', type=str, nargs='+',
                       choices=['cell_jsd', 'cell_methylation', 'gene_jsd', 'gene_methylation'],
                       default=['cell_jsd', 'cell_methylation', 'gene_jsd', 'gene_methylation'],
                       help='Which metrics to extract')
    parser.add_argument('--verbose', action='store_true',
                       help='Enable verbose output')

    args = parser.parse_args()

    # Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # Load simulation data
    print(f"Loading simulation from {args.simulation}...")
    with smart_open(args.simulation, 'r') as f:
        sim_data = json.load(f)

    # Extract history and config
    history = sim_data.get('history', {})
    config = sim_data.get('config', {})

    if not history:
        print("Error: No history found in simulation file")
        return 1

    # Determine year range
    available_years = sorted([int(y) for y in history.keys()])
    start_year = args.start_year
    end_year = args.end_year if args.end_year is not None else max(available_years)

    print(f"Extracting years {start_year} to {end_year}")
    print(f"Extracting metrics: {', '.join(args.metrics)}")

    # Extract each metric
    if 'cell_jsd' in args.metrics:
        print("\nExtracting cell JSD timeline...")
        data, max_cells = extract_cell_jsd_timeline(history, start_year, end_year)
        output_path = os.path.join(args.output_dir, 'simulation_cell_jsd.csv')
        save_to_csv(data, output_path, max_cells, column_prefix="value")

    if 'cell_methylation' in args.metrics:
        print("\nExtracting cell methylation timeline...")
        data, max_cells = extract_cell_methylation_timeline(history, start_year, end_year)
        output_path = os.path.join(args.output_dir, 'simulation_cell_methylation.csv')
        save_to_csv(data, output_path, max_cells, column_prefix="value")

    if 'gene_jsd' in args.metrics:
        print("\nExtracting gene JSD timeline...")
        data, n_genes = extract_gene_jsd_timeline(history, start_year, end_year)
        output_path = os.path.join(args.output_dir, 'simulation_gene_jsd.csv')
        save_to_csv(data, output_path, n_genes, column_prefix="gene")

    if 'gene_methylation' in args.metrics:
        print("\nExtracting gene methylation timeline...")
        data, n_genes = extract_gene_methylation_timeline(history, config, start_year, end_year)
        output_path = os.path.join(args.output_dir, 'simulation_gene_methylation.csv')
        save_to_csv(data, output_path, n_genes, column_prefix="gene")

    print(f"\nExtraction complete. Files saved to {args.output_dir}/")
    return 0


if __name__ == "__main__":
    sys.exit(main())