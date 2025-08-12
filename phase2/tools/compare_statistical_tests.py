#!/usr/bin/env python3
"""
Compare statistical test results between step23 and step23-prime.
"""

import json
import os

def load_stats(filepath):
    """Load statistics file."""
    if os.path.exists(filepath):
        with open(filepath) as f:
            return json.load(f)
    return None

def main():
    print("STATISTICAL TEST COMPARISON")
    print("=" * 50)
    
    # Load statistics files
    s23_stats = load_stats("../step23/data/rate_0.005000/results/statistics.json")
    pr_stats = load_stats("data/rate_0.005000/results/statistics.json")
    
    if not s23_stats:
        print("❌ step23 statistics file not found")
        return
    
    if not pr_stats:
        print("❌ step23-prime statistics file not found")
        return
    
    # Compare group statistics
    print("\nGROUP MEAN JSD:")
    print(f"{'Group':<15} {'step23':>12} {'phase2':>12} {'Difference':>12}")
    print("-" * 52)
    
    for group in ['mutant', 'control1', 'control2']:
        if group in s23_stats and group in pr_stats:
            s23_mean = s23_stats[group]['mean']
            pr_mean = pr_stats[group]['mean']
            diff = abs(s23_mean - pr_mean)
            
            symbol = "✅" if diff < 0.01 else "⚠️" if diff < 0.05 else "❌"
            print(f"{group:<15} {s23_mean:>12.6f} {pr_mean:>12.6f} {diff:>12.6f} {symbol}")
    
    # Compare p-values
    print("\nSTATISTICAL TEST P-VALUES:")
    print(f"{'Comparison':<25} {'step23':>12} {'phase2':>12} {'Same?':>8}")
    print("-" * 58)
    
    if 'comparisons' in s23_stats and 'comparisons' in pr_stats:
        for comp in s23_stats['comparisons'].keys():
            if comp in pr_stats['comparisons']:
                s23_p = s23_stats['comparisons'][comp]['p_value']
                pr_p = pr_stats['comparisons'][comp]['p_value']
                
                # Check if both lead to same conclusion (sig at 0.05)
                s23_sig = s23_p < 0.05
                pr_sig = pr_p < 0.05
                same = "✅ Yes" if s23_sig == pr_sig else "❌ No"
                
                print(f"{comp:<25} {s23_p:>12.6f} {pr_p:>12.6f} {same:>8}")
    
    print("\n" + "=" * 50)
    
    # Summary
    all_close = True
    for group in ['mutant', 'control1', 'control2']:
        if group in s23_stats and group in pr_stats:
            diff = abs(s23_stats[group]['mean'] - pr_stats[group]['mean'])
            if diff > 0.01:
                all_close = False
                break
    
    if all_close:
        print("✅ Group means are statistically equivalent (<0.01 difference)")
    else:
        print("⚠️ Some group means show differences (>0.01)")

if __name__ == "__main__":
    main()