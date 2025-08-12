#!/usr/bin/env python3
"""
Detailed comparison of individuals between step23 and step23-prime.
Tests statistical equivalence without requiring exact matches.
"""

import json
import gzip
import numpy as np
import glob
import os
from scipy import stats
from typing import Dict, List, Tuple
import warnings
warnings.filterwarnings('ignore')

def load_individual(filepath: str) -> Dict:
    """Load an individual file and extract cells."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    cells = data['cells'] if 'cells' in data else data
    return cells

def get_individual_statistics(cells: List[Dict]) -> Dict:
    """Calculate comprehensive statistics for an individual."""
    jsds = [c['jsd'] for c in cells]
    return {
        'n_cells': len(cells),
        'mean': np.mean(jsds),
        'std': np.std(jsds),
        'median': np.median(jsds),
        'min': np.min(jsds),
        'max': np.max(jsds),
        'q25': np.percentile(jsds, 25),
        'q75': np.percentile(jsds, 75),
        'cv': np.std(jsds) / np.mean(jsds) if np.mean(jsds) > 0 else 0
    }

def compare_group(group_name: str, step23_dir: str, prime_dir: str) -> Dict:
    """Compare all individuals in a group between pipelines."""
    
    print(f"\n{'='*60}")
    print(f"{group_name.upper()} INDIVIDUALS")
    print(f"{'='*60}")
    
    # Get file lists
    step23_files = sorted(glob.glob(f'{step23_dir}/individuals/{group_name}/individual_*.json.gz'))
    prime_files = sorted(glob.glob(f'{prime_dir}/individuals/{group_name}/individual_*.json.gz'))
    
    print(f"Files found: step23={len(step23_files)}, step23-prime={len(prime_files)}")
    
    if len(step23_files) != len(prime_files):
        print(f"❌ Different number of individuals!")
        return None
    
    # Collect statistics for all individuals
    step23_stats = []
    prime_stats = []
    cell_count_matches = 0
    
    for i, (s23_file, pr_file) in enumerate(zip(step23_files, prime_files)):
        try:
            # Load and analyze
            s23_cells = load_individual(s23_file)
            pr_cells = load_individual(pr_file)
            
            s23_stat = get_individual_statistics(s23_cells)
            pr_stat = get_individual_statistics(pr_cells)
            
            step23_stats.append(s23_stat)
            prime_stats.append(pr_stat)
            
            # Check cell counts
            if s23_stat['n_cells'] == pr_stat['n_cells']:
                cell_count_matches += 1
            else:
                print(f"  Individual {i:02d}: Cell count mismatch! {s23_stat['n_cells']} vs {pr_stat['n_cells']}")
        
        except Exception as e:
            print(f"  Error loading individual {i:02d}: {e}")
    
    # Calculate group-level statistics
    s23_means = [s['mean'] for s in step23_stats]
    pr_means = [s['mean'] for s in prime_stats]
    
    results = {
        'n_individuals': len(step23_stats),
        'cell_count_matches': cell_count_matches,
        'all_cells_correct': cell_count_matches == len(step23_stats),
        
        # Group statistics
        'step23': {
            'mean_of_means': np.mean(s23_means),
            'std_of_means': np.std(s23_means),
            'median_of_means': np.median(s23_means),
            'min_mean': np.min(s23_means),
            'max_mean': np.max(s23_means)
        },
        'prime': {
            'mean_of_means': np.mean(pr_means),
            'std_of_means': np.std(pr_means),
            'median_of_means': np.median(pr_means),
            'min_mean': np.min(pr_means),
            'max_mean': np.max(pr_means)
        },
        
        # Differences
        'mean_difference': abs(np.mean(s23_means) - np.mean(pr_means)),
        'std_difference': abs(np.std(s23_means) - np.std(pr_means))
    }
    
    # Statistical tests
    if len(s23_means) > 1:
        # T-test to see if means are significantly different
        t_stat, p_value = stats.ttest_ind(s23_means, pr_means)
        results['t_test'] = {'t_stat': t_stat, 'p_value': p_value}
        
        # KS test to compare distributions
        ks_stat, ks_p = stats.ks_2samp(s23_means, pr_means)
        results['ks_test'] = {'ks_stat': ks_stat, 'p_value': ks_p}
        
        # Correlation between the two sets
        correlation = np.corrcoef(s23_means, pr_means)[0, 1]
        results['correlation'] = correlation
    
    return results

def print_group_results(group_name: str, results: Dict):
    """Print formatted results for a group."""
    
    print(f"\n📊 {group_name.upper()} STATISTICS:")
    print(f"{'='*50}")
    
    # Basic counts
    print(f"✓ Individuals: {results['n_individuals']}")
    print(f"✓ Cell counts match: {results['cell_count_matches']}/{results['n_individuals']}")
    
    # Statistics comparison
    print(f"\nGroup Mean JSD:")
    print(f"  step23:       {results['step23']['mean_of_means']:.6f} ± {results['step23']['std_of_means']:.6f}")
    print(f"  step23-prime: {results['prime']['mean_of_means']:.6f} ± {results['prime']['std_of_means']:.6f}")
    print(f"  Difference:   {results['mean_difference']:.6f}")
    
    # Statistical tests
    if 't_test' in results:
        p_val = results['t_test']['p_value']
        sig = "✅ Not significant" if p_val > 0.05 else "⚠️ Significant"
        print(f"\nT-test: p={p_val:.4f} {sig}")
        
        ks_p = results['ks_test']['p_value']
        ks_sig = "✅ Same distribution" if ks_p > 0.05 else "⚠️ Different distribution"
        print(f"KS test: p={ks_p:.4f} {ks_sig}")
        
        corr = results['correlation']
        print(f"Correlation: {corr:.4f}")
    
    # Verdict
    if results['mean_difference'] < 0.01:
        print(f"\n✅ EQUIVALENT: Mean difference < 0.01")
    elif results['mean_difference'] < 0.05:
        print(f"⚠️ SIMILAR: Mean difference < 0.05")
    else:
        print(f"❌ DIFFERENT: Mean difference > 0.05")

def compare_cross_group_statistics(step23_dir: str, prime_dir: str):
    """Compare the statistical test results between groups."""
    
    print(f"\n{'='*60}")
    print("CROSS-GROUP STATISTICAL TESTS")
    print(f"{'='*60}")
    
    # Load statistics files if they exist
    step23_stats_file = f"{step23_dir}/results/statistics.json"
    prime_stats_file = f"{prime_dir}/results/statistics.json"
    
    if os.path.exists(step23_stats_file) and os.path.exists(prime_stats_file):
        with open(step23_stats_file) as f:
            s23_stats = json.load(f)
        with open(prime_stats_file) as f:
            pr_stats = json.load(f)
        
        print("\nP-values from t-tests:")
        print(f"{'Comparison':<25} {'step23':>12} {'phase2':>12} {'Difference':>12}")
        print("-" * 62)
        
        for comp in ['mutant_vs_control1', 'mutant_vs_control2', 'control1_vs_control2']:
            if 'comparisons' in s23_stats and 'comparisons' in pr_stats:
                s23_p = s23_stats['comparisons'][comp]['p_value']
                pr_p = pr_stats['comparisons'][comp]['p_value']
                diff = abs(s23_p - pr_p)
                
                print(f"{comp:<25} {s23_p:>12.6f} {pr_p:>12.6f} {diff:>12.6f}")
                
                # Check if conclusions would be the same
                s23_sig = s23_p < 0.05
                pr_sig = pr_p < 0.05
                if s23_sig == pr_sig:
                    print(f"  ✅ Same conclusion (both {'significant' if s23_sig else 'not significant'})")
                else:
                    print(f"  ⚠️ Different conclusions!")
    else:
        print("Statistics files not found, skipping cross-group tests")

def main():
    """Run comprehensive comparison of individuals."""
    
    print("="*60)
    print("COMPREHENSIVE INDIVIDUAL COMPARISON")
    print("step23 vs step23-prime")
    print("="*60)
    
    # Directories
    step23_dir = "../step23/data/rate_0.005000"
    prime_dir = "data/rate_0.005000"
    
    # Check directories exist
    if not os.path.exists(step23_dir):
        print(f"❌ step23 directory not found: {step23_dir}")
        return
    if not os.path.exists(prime_dir):
        print(f"❌ step23-prime directory not found: {prime_dir}")
        return
    
    # Compare each group
    all_results = {}
    for group in ['mutant', 'control1', 'control2']:
        results = compare_group(group, step23_dir, prime_dir)
        if results:
            all_results[group] = results
            print_group_results(group, results)
    
    # Cross-group statistics
    compare_cross_group_statistics(step23_dir, prime_dir)
    
    # Final summary
    print(f"\n{'='*60}")
    print("FINAL VERDICT")
    print(f"{'='*60}")
    
    if len(all_results) == 3:
        # Check overall equivalence
        all_equivalent = all(r['mean_difference'] < 0.01 for r in all_results.values())
        all_similar = all(r['mean_difference'] < 0.05 for r in all_results.values())
        all_cells_match = all(r['all_cells_correct'] for r in all_results.values())
        
        print("\n📋 Summary:")
        print(f"  Cell counts: {'✅ All match' if all_cells_match else '❌ Some mismatches'}")
        
        for group, res in all_results.items():
            diff = res['mean_difference']
            symbol = "✅" if diff < 0.01 else "⚠️" if diff < 0.05 else "❌"
            print(f"  {group.capitalize()}: Δ={diff:.6f} {symbol}")
        
        if all_equivalent:
            print("\n✅ PIPELINES ARE STATISTICALLY EQUIVALENT")
            print("All group means within 0.01 JSD - excellent agreement!")
        elif all_similar:
            print("\n⚠️ PIPELINES ARE SIMILAR")
            print("All group means within 0.05 JSD - acceptable agreement")
        else:
            print("\n❌ PIPELINES SHOW SIGNIFICANT DIFFERENCES")
            print("Some groups differ by >0.05 JSD - investigate further")
        
        # Correlation summary
        print("\nCorrelations between pipelines:")
        for group, res in all_results.items():
            if 'correlation' in res:
                corr = res['correlation']
                print(f"  {group}: {corr:.4f}")
    else:
        print("❌ Could not complete comparison")

if __name__ == "__main__":
    main()