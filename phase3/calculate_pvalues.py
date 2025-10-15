#!/usr/bin/env python
"""
Calculate p-values for batch comparisons using pairwise t-tests and one-way ANOVA.
Pairwise tests: Student's t-test and Welch's t-test
ANOVA tests: Fisher's ANOVA and Welch's ANOVA
Generates tables for all 6 comparison metrics.
"""

import os
import sys
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from typing import Dict, Tuple


def calculate_pvalues_for_metric(
    df: pd.DataFrame,
    value_column: str,
    aggregation: str = None
) -> Dict[str, Tuple[float, float]]:
    """
    Calculate p-values for all pairwise batch comparisons.

    Args:
        df: DataFrame with batch comparison data
        value_column: Column name containing the values to compare
        aggregation: 'mean' or 'std' for gene data aggregation (None for cell data)

    Returns:
        Dictionary mapping comparison pairs to (student_pvalue, welch_pvalue)
    """
    # If gene data with multiple genes, aggregate first
    if 'gene_index' in df.columns and aggregation is not None:
        if df['gene_index'].nunique() > 1:
            if aggregation == 'mean':
                df = df.groupby(['batch', 'individual_id'])[value_column].mean().reset_index()
            elif aggregation == 'std':
                df = df.groupby(['batch', 'individual_id'])[value_column].std().reset_index()
            else:
                raise ValueError(f"Unsupported aggregation: {aggregation}")

    # Extract data for each batch
    batches = {}
    for batch_name in ['control', 'test1', 'test2']:
        batch_data = df[df['batch'] == batch_name][value_column].values
        batches[batch_name] = batch_data

    # Define pairwise comparisons
    comparisons = [
        ('control', 'test1'),
        ('control', 'test2'),
        ('test1', 'test2')
    ]

    # Calculate p-values for each comparison
    results = {}
    for batch1, batch2 in comparisons:
        data1 = batches[batch1]
        data2 = batches[batch2]

        # Student's t-test (equal variance assumed)
        student_stat, student_pval = stats.ttest_ind(data1, data2, equal_var=True)

        # Welch's t-test (equal variance NOT assumed)
        welch_stat, welch_pval = stats.ttest_ind(data1, data2, equal_var=False)

        pair_key = f"{batch1} vs {batch2}"
        results[pair_key] = (student_pval, welch_pval)

    return results


def format_pvalue(pval: float) -> str:
    """Format p-value with appropriate precision."""
    if pval < 0.001:
        return f"{pval:.2e}"
    else:
        return f"{pval:.4f}"


def format_fstat(fstat: float) -> str:
    """Format F-statistic with appropriate precision."""
    return f"{fstat:.4f}"


def format_df(df_tuple: Tuple[float, float]) -> str:
    """Format degrees of freedom."""
    df1, df2 = df_tuple
    # If df2 is an integer, format as integer
    if abs(df2 - round(df2)) < 0.01:
        return f"({int(df1)}, {int(df2)})"
    else:
        # Welch's ANOVA has non-integer df2
        return f"({int(df1)}, {df2:.2f})"


def calculate_anova_for_metric(
    df: pd.DataFrame,
    value_column: str,
    aggregation: str = None
) -> Dict[str, Tuple[float, Tuple[float, float], float]]:
    """
    Calculate ANOVA for comparing all three batches.

    Args:
        df: DataFrame with batch comparison data
        value_column: Column name containing the values to compare
        aggregation: 'mean' or 'std' for gene data aggregation (None for cell data)

    Returns:
        Dictionary with keys 'fisher' and 'welch',
        values are tuples of (f_statistic, (df1, df2), p_value)
    """
    # If gene data with multiple genes, aggregate first
    if 'gene_index' in df.columns and aggregation is not None:
        if df['gene_index'].nunique() > 1:
            if aggregation == 'mean':
                df = df.groupby(['batch', 'individual_id'])[value_column].mean().reset_index()
            elif aggregation == 'std':
                df = df.groupby(['batch', 'individual_id'])[value_column].std().reset_index()
            else:
                raise ValueError(f"Unsupported aggregation: {aggregation}")

    # Extract data for each batch
    data_control = df[df['batch'] == 'control'][value_column].values
    data_test1 = df[df['batch'] == 'test1'][value_column].values
    data_test2 = df[df['batch'] == 'test2'][value_column].values

    # Calculate sample sizes
    n_control = len(data_control)
    n_test1 = len(data_test1)
    n_test2 = len(data_test2)
    n_total = n_control + n_test1 + n_test2
    k_groups = 3  # number of groups

    results = {}

    # Fisher's ANOVA (assumes equal variances)
    f_stat_fisher, p_val_fisher = stats.f_oneway(data_control, data_test1, data_test2)
    df1 = k_groups - 1  # 2
    df2_fisher = n_total - k_groups  # typically 6 for 9 individuals
    results['fisher'] = (f_stat_fisher, (df1, df2_fisher), p_val_fisher)

    # Welch's ANOVA (does not assume equal variances)
    # Use Alexander-Govern test or Welch's ANOVA
    # For scipy >= 1.15.0, f_oneway supports equal_var parameter
    try:
        # Try using equal_var parameter (scipy >= 1.15.0)
        f_stat_welch_result = stats.f_oneway(data_control, data_test1, data_test2, equal_var=False)
        # Extract F-statistic and p-value
        if hasattr(f_stat_welch_result, 'statistic'):
            f_stat_welch = f_stat_welch_result.statistic
            p_val_welch = f_stat_welch_result.pvalue
        else:
            f_stat_welch, p_val_welch = f_stat_welch_result

        # Calculate Welch-Satterthwaite degrees of freedom
        # For Welch's ANOVA, df2 is adjusted based on variances
        vars = [np.var(data_control, ddof=1), np.var(data_test1, ddof=1), np.var(data_test2, ddof=1)]
        ns = [n_control, n_test1, n_test2]

        # Simplified Welch df calculation
        numerator = sum([(1 - ni/n_total) * vi / ni for ni, vi in zip(ns, vars)])**2
        denominator = sum([((1 - ni/n_total)**2 * vi**2) / (ni**2 * (ni - 1)) for ni, vi in zip(ns, vars)])
        df2_welch = numerator / denominator if denominator > 0 else df2_fisher

        results['welch'] = (f_stat_welch, (df1, df2_welch), p_val_welch)
    except TypeError:
        # Fallback for older scipy versions without equal_var parameter
        # Use regular f_oneway and note limitation
        f_stat_welch, p_val_welch = stats.f_oneway(data_control, data_test1, data_test2)
        results['welch'] = (f_stat_welch, (df1, df2_fisher), p_val_welch)

    return results


def print_comparison_table(
    metric_name: str,
    pvalue_results: Dict[str, Tuple[float, float]]
):
    """
    Print a formatted table of p-values for one metric.

    Args:
        metric_name: Name of the metric being compared
        pvalue_results: Dictionary mapping pairs to (student_pval, welch_pval)
    """
    print(f"\n{'='*70}")
    print(f"{metric_name}")
    print(f"{'='*70}")

    # Format header without backslash in f-string
    header = f"{'Comparison':<25} {'Student t-test':<20} {'Welch t-test':<20}"
    print(header)
    print(f"{'-'*70}")

    for pair, (student_pval, welch_pval) in pvalue_results.items():
        print(f"{pair:<25} {format_pvalue(student_pval):<20} {format_pvalue(welch_pval):<20}")


def save_pvalues_to_csv(
    all_results: Dict[str, Dict[str, Tuple[float, float]]],
    output_path: str
):
    """
    Save all p-value results to a CSV file.

    Args:
        all_results: Dictionary mapping metric names to their p-value results
        output_path: Path to save CSV file
    """
    rows = []
    for metric_name, pvalue_results in all_results.items():
        for pair, (student_pval, welch_pval) in pvalue_results.items():
            rows.append({
                'metric': metric_name,
                'comparison': pair,
                'students_t_pvalue': student_pval,
                'welchs_t_pvalue': welch_pval
            })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"\n✓ Saved all p-values to {output_path}")


def print_anova_table(
    all_anova_results: Dict[str, Dict[str, Tuple[float, Tuple[float, float], float]]]
):
    """
    Print formatted ANOVA results table.

    Args:
        all_anova_results: Dictionary mapping metric names to their ANOVA results
    """
    print("\n" + "="*90)
    print("ONE-WAY ANOVA RESULTS")
    print("="*90)

    # Header
    header = f"{'Metric':<28} {'Test Type':<18} {'F-statistic':<15} {'df':<15} {'p-value':<12}"
    print(header)
    print("-"*90)

    # Print results for each metric
    for metric_name, anova_results in all_anova_results.items():
        # Fisher's ANOVA
        fisher_f, fisher_df, fisher_p = anova_results['fisher']
        print(f"{metric_name:<28} {'Fisher ANOVA':<18} {format_fstat(fisher_f):<15} {format_df(fisher_df):<15} {format_pvalue(fisher_p):<12}")

        # Welch's ANOVA
        welch_f, welch_df, welch_p = anova_results['welch']
        print(f"{'':<28} {'Welch ANOVA':<18} {format_fstat(welch_f):<15} {format_df(welch_df):<15} {format_pvalue(welch_p):<12}")

        # Separator between metrics
        print("-"*90)

    print("="*90)


def save_anova_to_csv(
    all_anova_results: Dict[str, Dict[str, Tuple[float, Tuple[float, float], float]]],
    output_path: str
):
    """
    Save ANOVA results to CSV file.

    Args:
        all_anova_results: Dictionary mapping metric names to their ANOVA results
        output_path: Path to save CSV file
    """
    rows = []
    for metric_name, anova_results in all_anova_results.items():
        # Fisher's ANOVA
        fisher_f, fisher_df, fisher_p = anova_results['fisher']
        rows.append({
            'metric': metric_name,
            'test_type': 'Fisher ANOVA',
            'f_statistic': fisher_f,
            'df1': fisher_df[0],
            'df2': fisher_df[1],
            'p_value': fisher_p
        })

        # Welch's ANOVA
        welch_f, welch_df, welch_p = anova_results['welch']
        rows.append({
            'metric': metric_name,
            'test_type': 'Welch ANOVA',
            'f_statistic': welch_f,
            'df1': welch_df[0],
            'df2': welch_df[1],
            'p_value': welch_p
        })

    df = pd.DataFrame(rows)
    df.to_csv(output_path, index=False)
    print(f"\n✓ Saved ANOVA results to {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description='Calculate p-values for batch comparisons',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument('--cell-csv', required=True,
                       help='Path to cell comparison CSV')
    parser.add_argument('--gene-csv', required=True,
                       help='Path to gene comparison CSV')
    parser.add_argument('--output-csv', default=None,
                       help='Optional: Save pairwise p-values to CSV file')
    parser.add_argument('--anova-csv', default=None,
                       help='Optional: Save ANOVA results to CSV file')
    parser.add_argument('--quiet', action='store_true',
                       help='Suppress table output, only save CSV')

    args = parser.parse_args()

    # Read CSV files
    try:
        df_cell = pd.read_csv(args.cell_csv)
        df_gene = pd.read_csv(args.gene_csv)
    except FileNotFoundError as e:
        print(f"Error: {e}")
        sys.exit(1)

    print("\n" + "="*90)
    print("STATISTICAL TESTS FOR BATCH COMPARISONS")
    print("="*90)
    print(f"Cell CSV: {args.cell_csv}")
    print(f"Gene CSV: {args.gene_csv}")
    print("\nPairwise Tests:")
    print("  - Student's t-test (assumes equal variance)")
    print("  - Welch's t-test (does not assume equal variance)")
    print("\nOmnibus Tests:")
    print("  - Fisher's ANOVA (assumes equal variance)")
    print("  - Welch's ANOVA (does not assume equal variance)")
    print("\nNote: Bonferroni correction should be applied during interpretation")
    print("="*90)

    # Store all results for CSV export
    all_results = {}

    # 1. Cell JSD Mean
    metric_name = "Cell JSD Mean"
    pvalue_results = calculate_pvalues_for_metric(df_cell, 'cell_jsd_mean')
    all_results[metric_name] = pvalue_results
    if not args.quiet:
        print_comparison_table(metric_name, pvalue_results)

    # 2. Cell Methylation Mean
    metric_name = "Cell Methylation Mean"
    pvalue_results = calculate_pvalues_for_metric(df_cell, 'cell_methylation_mean')
    all_results[metric_name] = pvalue_results
    if not args.quiet:
        print_comparison_table(metric_name, pvalue_results)

    # 3. Gene JSD Mean
    metric_name = "Gene JSD Mean"
    pvalue_results = calculate_pvalues_for_metric(df_gene, 'gene_jsd', aggregation='mean')
    all_results[metric_name] = pvalue_results
    if not args.quiet:
        print_comparison_table(metric_name, pvalue_results)

    # 4. Gene Methylation Mean
    metric_name = "Gene Methylation Mean"
    pvalue_results = calculate_pvalues_for_metric(df_gene, 'gene_methylation', aggregation='mean')
    all_results[metric_name] = pvalue_results
    if not args.quiet:
        print_comparison_table(metric_name, pvalue_results)

    # 5. Gene JSD Standard Deviation
    metric_name = "Gene JSD Std Dev"
    pvalue_results = calculate_pvalues_for_metric(df_gene, 'gene_jsd', aggregation='std')
    all_results[metric_name] = pvalue_results
    if not args.quiet:
        print_comparison_table(metric_name, pvalue_results)

    # 6. Gene Methylation Standard Deviation
    metric_name = "Gene Methylation Std Dev"
    pvalue_results = calculate_pvalues_for_metric(df_gene, 'gene_methylation', aggregation='std')
    all_results[metric_name] = pvalue_results
    if not args.quiet:
        print_comparison_table(metric_name, pvalue_results)

    # Save to CSV if requested
    if args.output_csv:
        save_pvalues_to_csv(all_results, args.output_csv)

    # ========================================================================
    # ANOVA Calculations
    # ========================================================================
    print("\n" + "="*90)
    print("CALCULATING ANOVA (Fisher's and Welch's)")
    print("="*90)

    all_anova_results = {}

    # 1. Cell JSD Mean
    metric_name = "Cell JSD Mean"
    anova_results = calculate_anova_for_metric(df_cell, 'cell_jsd_mean')
    all_anova_results[metric_name] = anova_results

    # 2. Cell Methylation Mean
    metric_name = "Cell Methylation Mean"
    anova_results = calculate_anova_for_metric(df_cell, 'cell_methylation_mean')
    all_anova_results[metric_name] = anova_results

    # 3. Gene JSD Mean
    metric_name = "Gene JSD Mean"
    anova_results = calculate_anova_for_metric(df_gene, 'gene_jsd', aggregation='mean')
    all_anova_results[metric_name] = anova_results

    # 4. Gene Methylation Mean
    metric_name = "Gene Methylation Mean"
    anova_results = calculate_anova_for_metric(df_gene, 'gene_methylation', aggregation='mean')
    all_anova_results[metric_name] = anova_results

    # 5. Gene JSD Standard Deviation
    metric_name = "Gene JSD Std Dev"
    anova_results = calculate_anova_for_metric(df_gene, 'gene_jsd', aggregation='std')
    all_anova_results[metric_name] = anova_results

    # 6. Gene Methylation Standard Deviation
    metric_name = "Gene Methylation Std Dev"
    anova_results = calculate_anova_for_metric(df_gene, 'gene_methylation', aggregation='std')
    all_anova_results[metric_name] = anova_results

    # Print ANOVA table
    if not args.quiet:
        print_anova_table(all_anova_results)

    # Save ANOVA to CSV if requested
    if args.anova_csv:
        save_anova_to_csv(all_anova_results, args.anova_csv)

    print("\n" + "="*90)
    print("ALL CALCULATIONS COMPLETE")
    print("="*90)


if __name__ == "__main__":
    main()
