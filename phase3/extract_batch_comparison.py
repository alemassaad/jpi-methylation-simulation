#!/usr/bin/env python3
"""
Extract batch comparison data from individual PetriDish JSON files.
Calculates mean metrics for each individual and saves to CSV.
"""

import os
import sys
import json
import gzip
import csv
import numpy as np
from typing import Dict, List, Any, Optional, Tuple
import glob

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from data_extractors import smart_open, dict_to_cell


def load_petri_dish(filepath: str) -> PetriDish:
    """
    Load a PetriDish from a JSON file.

    Phase 2 individual files have this structure:
    - 'config': Configuration and metadata
    - 'history': Year-by-year growth history (years 23-38, ~68 cells at end)
    - 'individual_final': MIXED population after combining grown cells with snapshot cells (~288 cells)

    We MUST use 'individual_final' because it contains the final mixed population.
    The 'history' only contains the growth phase BEFORE mixing.

    Args:
        filepath: Path to individual JSON file

    Returns:
        PetriDish object with cells from the final mixed population
    """
    with smart_open(filepath, 'r') as f:
        data = json.load(f)

    # Handle Phase 2 format with history and individual_final
    if 'individual_final' in data:
        # CORRECT: Use the final mixed population
        individual_data = data['individual_final']
        cells = [dict_to_cell(cell_dict) for cell_dict in individual_data['cells']]

        # Log cell count for validation
        if len(cells) < 150:
            print(f"    Warning: Only {len(cells)} cells in individual_final (expected ~210 with default config)")

    elif 'history' in data:
        # FALLBACK: Use the last year in history (WARNING: this is pre-mixing!)
        print(f"    WARNING: Using history instead of individual_final - this may be pre-mixing data!")
        history = data['history']
        # Get the last year (highest number)
        years = [int(y) for y in history.keys()]
        last_year = str(max(years))
        cells = [dict_to_cell(cell_dict) for cell_dict in history[last_year]['cells']]
        print(f"    Using year {last_year} with {len(cells)} cells (likely pre-mixing)")

    elif 'cells' in data:
        # Direct cells format (for backward compatibility or different format)
        cells = [dict_to_cell(cell_dict) for cell_dict in data['cells']]
    else:
        raise ValueError(
            f"Cannot find cells in {filepath}. "
            f"Expected 'individual_final' (mixed population) or 'cells'. "
            f"Found keys: {list(data.keys())}"
        )

    if cells:
        # Get configuration from first cell
        first_cell = cells[0]
        n = first_cell.n
        gene_size = first_cell.gene_size if hasattr(first_cell, 'gene_size') else 5
        gene_rate_groups = first_cell.gene_rate_groups if hasattr(first_cell, 'gene_rate_groups') else None

        # Create PetriDish with matching parameters
        petri = PetriDish(
            n=n,
            gene_size=gene_size,
            gene_rate_groups=gene_rate_groups,
            cells=cells
        )
    else:
        raise ValueError(f"No cells found in {filepath}")

    # Load metadata if present
    if 'config' in data and 'metadata' in data['config']:
        petri.metadata = data['config']['metadata']
    elif 'metadata' in data:
        petri.metadata = data['metadata']

    return petri


def calculate_cell_metrics(petri: PetriDish) -> Dict[str, float]:
    """
    Calculate cell-level comparison metrics for one individual.

    Args:
        petri: PetriDish object

    Returns:
        Dictionary with cell metric names and values
    """
    metrics = {}

    # Cell-level JSD: mean of all cells' JSD values
    cell_jsds = [cell.cell_jsd for cell in petri.cells]
    metrics['cell_jsd_mean'] = float(np.mean(cell_jsds)) if cell_jsds else 0.0

    # Cell-level methylation: mean of all cells' methylation proportions
    cell_methylations = []
    for cell in petri.cells:
        if hasattr(cell, 'cpg_sites') and cell.cpg_sites:
            prop = np.sum(cell.cpg_sites) / len(cell.cpg_sites)
        else:
            prop = 0.0
        cell_methylations.append(prop)
    metrics['cell_methylation_mean'] = float(np.mean(cell_methylations)) if cell_methylations else 0.0

    # Add cell count
    metrics['cell_count'] = len(petri.cells)

    return metrics




def _process_individuals(individuals_dir: str, verbose: bool = True) -> List[Tuple[str, int, Any]]:
    """
    Helper function to load all individuals from batches.

    Args:
        individuals_dir: Directory containing test2/, test1/, control/ subdirectories
        verbose: Whether to print progress

    Returns:
        List of tuples: (batch_name, individual_id, petri_dish)
    """
    if verbose:
        print(f"\nProcessing individuals from {individuals_dir}")

    # Check that individuals directory exists
    if not os.path.exists(individuals_dir):
        raise FileNotFoundError(f"Individuals directory not found: {individuals_dir}")

    all_individuals = []

    # Process each batch (Control, Test 1, Test 2 order)
    for batch_name in ['control', 'test1', 'test2']:
        batch_dir = os.path.join(individuals_dir, batch_name)

        if not os.path.exists(batch_dir):
            if verbose:
                print(f"  Warning: {batch_name} directory not found, skipping")
            continue

        # Find all individual JSON files
        json_files = sorted(glob.glob(os.path.join(batch_dir, "individual_*.json*")))

        if verbose:
            print(f"  Processing {batch_name}: {len(json_files)} individuals")

        for idx, filepath in enumerate(json_files):
            try:
                # Load PetriDish
                petri = load_petri_dish(filepath)
                all_individuals.append((batch_name, idx, petri))

            except Exception as e:
                if verbose:
                    print(f"    Warning: Error processing {filepath}: {e}")
                continue

    return all_individuals


def extract_cell_comparison(individuals_dir: str, output_csv: str, verbose: bool = True) -> None:
    """
    Extract cell-level comparison data from all individuals and save to CSV.

    Args:
        individuals_dir: Directory containing test2/, test1/, control/ subdirectories
        output_csv: Path to output CSV file for cell metrics
        verbose: Whether to print progress
    """
    # Load all individuals
    all_individuals = _process_individuals(individuals_dir, verbose)

    if verbose:
        print(f"\nExtracting cell-level metrics...")

    # Prepare data for CSV
    csv_rows = []

    # Process each individual
    for batch_name, individual_id, petri in all_individuals:
        # Calculate cell metrics
        metrics = calculate_cell_metrics(petri)

        # Add to CSV rows
        row = {
            'batch': batch_name,
            'individual_id': individual_id,
            'cell_jsd_mean': metrics['cell_jsd_mean'],
            'cell_methylation_mean': metrics['cell_methylation_mean'],
            'cell_count': metrics['cell_count']
        }
        csv_rows.append(row)

    # Create output directory if needed
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Write CSV
    if csv_rows:
        fieldnames = ['batch', 'individual_id', 'cell_jsd_mean', 'cell_methylation_mean', 'cell_count']

        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)

        if verbose:
            print(f"  Saved {len(csv_rows)} rows to {output_csv}")
    else:
        raise ValueError("No data extracted from any individuals")


def extract_gene_comparison(individuals_dir: str, output_csv: str, verbose: bool = True) -> None:
    """
    Extract gene-level comparison data from all individuals and save to CSV.

    Outputs the full distribution: one row per gene per individual (20 rows per individual).

    Args:
        individuals_dir: Directory containing test2/, test1/, control/ subdirectories
        output_csv: Path to output CSV file for gene metrics
        verbose: Whether to print progress
    """
    # Load all individuals
    all_individuals = _process_individuals(individuals_dir, verbose)

    if verbose:
        print(f"\nExtracting gene-level metrics (full distribution)...")

    # Prepare data for CSV
    csv_rows = []

    # Process each individual
    for batch_name, individual_id, petri in all_individuals:
        # Calculate gene-level JSDs
        try:
            gene_jsds = petri.calculate_gene_jsd()
        except:
            gene_jsds = []

        # Calculate gene-level methylation proportions
        if petri.cells and hasattr(petri.cells[0], 'gene_size'):
            gene_size = petri.cells[0].gene_size
            n_sites = petri.cells[0].n
            n_genes = n_sites // gene_size

            # Calculate methylation for each gene across all cells
            for gene_idx in range(n_genes):
                start_idx = gene_idx * gene_size
                end_idx = (gene_idx + 1) * gene_size

                # Average methylation for this gene across all cells
                gene_meth_values = []
                for cell in petri.cells:
                    if hasattr(cell, 'cpg_sites'):
                        gene_sites = cell.cpg_sites[start_idx:end_idx]
                        gene_meth = np.sum(gene_sites) / len(gene_sites) if gene_sites else 0.0
                        gene_meth_values.append(gene_meth)

                # Calculate mean methylation for this gene
                gene_methylation = float(np.mean(gene_meth_values)) if gene_meth_values else 0.0

                # Get JSD for this gene (if available)
                gene_jsd = float(gene_jsds[gene_idx]) if gene_idx < len(gene_jsds) else 0.0

                # Add row for this gene
                row = {
                    'batch': batch_name,
                    'individual_id': individual_id,
                    'gene_index': gene_idx,
                    'gene_jsd': gene_jsd,
                    'gene_methylation': gene_methylation
                }
                csv_rows.append(row)

    # Create output directory if needed
    output_dir = os.path.dirname(output_csv)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    # Write CSV
    if csv_rows:
        fieldnames = ['batch', 'individual_id', 'gene_index', 'gene_jsd', 'gene_methylation']

        with open(output_csv, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)

        if verbose:
            n_individuals = len(all_individuals)
            print(f"  Saved {len(csv_rows)} rows ({n_individuals} individuals × 20 genes) to {output_csv}")
    else:
        raise ValueError("No data extracted from any individuals")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Extract batch comparison data to separate cell and gene CSV files")
    parser.add_argument('--individuals-dir', required=True,
                       help='Directory containing test2/test1/control subdirectories')
    parser.add_argument('--cell-csv', required=False,
                       help='Output CSV file path for cell metrics only')
    parser.add_argument('--gene-csv', required=False,
                       help='Output CSV file path for gene metrics only')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress progress output')

    args = parser.parse_args()

    # Determine output paths
    if args.cell_csv or args.gene_csv:
        # Use provided paths
        if args.cell_csv:
            extract_cell_comparison(
                individuals_dir=args.individuals_dir,
                output_csv=args.cell_csv,
                verbose=not args.quiet
            )
        if args.gene_csv:
            extract_gene_comparison(
                individuals_dir=args.individuals_dir,
                output_csv=args.gene_csv,
                verbose=not args.quiet
            )
    else:
        # Default: create both cell and gene CSVs with standard names
        base_dir = os.path.dirname(args.individuals_dir)
        tables_dir = os.path.join(base_dir, 'results', 'tables')
        os.makedirs(tables_dir, exist_ok=True)

        extract_cell_comparison(
            individuals_dir=args.individuals_dir,
            output_csv=os.path.join(tables_dir, 'batch_comparison_cell.csv'),
            verbose=not args.quiet
        )
        extract_gene_comparison(
            individuals_dir=args.individuals_dir,
            output_csv=os.path.join(tables_dir, 'batch_comparison_gene.csv'),
            verbose=not args.quiet
        )