#!/usr/bin/env python3
"""
Generate per-gene comparison plots for all 20 genes.
Creates separate comparison plots for each gene's JSD and methylation proportion.
"""

import os
import pandas as pd
from typing import Optional
from plot_comparison_generic import plot_comparison_generic


def plot_comparison_by_gene(
    gene_csv_path: str,
    output_dir: str,
    verbose: bool = True,
    figsize: tuple = (1200, 600),
    scale: int = 2
) -> int:
    """
    Generate per-gene comparison plots for all 20 genes.

    Args:
        gene_csv_path: Path to batch_comparison_gene.csv file
        output_dir: Directory to save the per-gene plots
        verbose: Whether to print progress messages
        figsize: Figure size in pixels (width, height)
        scale: Scale factor for high DPI

    Returns:
        Number of plots successfully created
    """
    if verbose:
        print(f"\nGenerating per-gene comparison plots...")
        print(f"  Reading data from: {gene_csv_path}")
        print(f"  Output directory: {output_dir}")

    # Read the gene comparison CSV
    try:
        df = pd.read_csv(gene_csv_path)
    except Exception as e:
        print(f"Error reading CSV file: {str(e)}")
        return 0

    # Verify required columns exist - handle both possible column names for methylation
    methylation_col = None
    if 'gene_methylation_proportion' in df.columns:
        methylation_col = 'gene_methylation_proportion'
    elif 'gene_methylation' in df.columns:
        methylation_col = 'gene_methylation'
    else:
        print(f"Error: No gene methylation column found")
        print(f"Available columns: {df.columns.tolist()}")
        return 0

    required_columns = ['gene_index', 'gene_jsd', 'batch', 'individual_id']
    missing_columns = [col for col in required_columns if col not in df.columns]
    if missing_columns:
        print(f"Error: Missing required columns: {missing_columns}")
        print(f"Available columns: {df.columns.tolist()}")
        return 0

    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)

    # Get unique gene indices (should be 0-19)
    unique_genes = sorted(df['gene_index'].unique())
    n_genes = len(unique_genes)

    if verbose:
        print(f"  Found {n_genes} unique genes: {unique_genes[0]}-{unique_genes[-1]}")
        print(f"  Generating 2 plots per gene ({n_genes * 2} total)...")

    plot_count = 0

    # Loop through each gene
    for gene_idx in unique_genes:
        if verbose and gene_idx % 5 == 0:  # Progress indicator every 5 genes
            print(f"    Processing gene {gene_idx}...")

        # Filter data for this specific gene
        gene_df = df[df['gene_index'] == gene_idx].copy()

        # Check if we have data for all batches
        batches_present = gene_df['batch'].unique()
        if len(batches_present) < 3:
            if verbose:
                print(f"    Warning: Gene {gene_idx} missing some batches (found: {batches_present})")

        # Generate JSD comparison plot for this gene
        jsd_output_path = os.path.join(output_dir, f"gene_{gene_idx}_jsd_comparison.png")
        try:
            # Pass the filtered dataframe directly
            plot_comparison_generic(
                csv_path=gene_df,  # Pass DataFrame directly
                value_column='gene_jsd',
                title=f'Gene {gene_idx} JSD Comparison Across Batches',
                ylabel='Gene JSD',
                output_path=jsd_output_path,
                figsize=figsize,
                scale=scale,
                verbose=False,  # Suppress individual plot messages
                y_range_padding=0.1  # 10% padding to prevent edge clipping
            )

            plot_count += 1

        except Exception as e:
            if verbose:
                print(f"    Error creating JSD plot for gene {gene_idx}: {str(e)}")

        # Generate methylation comparison plot for this gene
        methylation_output_path = os.path.join(output_dir, f"gene_{gene_idx}_methylation_comparison.png")
        try:
            plot_comparison_generic(
                csv_path=gene_df,  # Pass DataFrame directly
                value_column=methylation_col,  # Use the detected column name
                title=f'Gene {gene_idx} Methylation Comparison Across Batches',
                ylabel='Gene Methylation Proportion',
                output_path=methylation_output_path,
                figsize=figsize,
                scale=scale,
                verbose=False,  # Suppress individual plot messages
                y_range_padding=0.1  # 10% padding to prevent edge clipping
            )

            plot_count += 1

        except Exception as e:
            if verbose:
                print(f"    Error creating methylation plot for gene {gene_idx}: {str(e)}")

    if verbose:
        print(f"  ✓ Successfully created {plot_count}/{n_genes * 2} plots")
        if plot_count < n_genes * 2:
            print(f"  ⚠ {n_genes * 2 - plot_count} plots failed to generate")

    return plot_count


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description="Generate per-gene comparison plots from batch comparison data"
    )
    parser.add_argument(
        '--gene-csv',
        required=True,
        help='Path to batch_comparison_gene.csv file'
    )
    parser.add_argument(
        '--output-dir',
        default='plots/comparison_by_gene',
        help='Directory to save per-gene plots (default: plots/comparison_by_gene)'
    )
    parser.add_argument(
        '--quiet',
        action='store_true',
        help='Suppress progress messages'
    )

    args = parser.parse_args()

    # Run the plotting
    n_plots = plot_comparison_by_gene(
        gene_csv_path=args.gene_csv,
        output_dir=args.output_dir,
        verbose=not args.quiet
    )

    if not args.quiet:
        print(f"\nComplete: Generated {n_plots} plots in {args.output_dir}")