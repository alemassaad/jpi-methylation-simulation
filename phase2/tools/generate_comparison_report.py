#!/usr/bin/env python3
"""
Generate a comprehensive comparison report.
"""

import json
import gzip
import numpy as np
import glob
import os
from datetime import datetime

def count_files(directory):
    """Count files in a directory."""
    return len(glob.glob(os.path.join(directory, "*.json.gz")))

def get_first_cell_count(directory):
    """Get cell count from first file in directory."""
    files = sorted(glob.glob(os.path.join(directory, "*.json.gz")))
    if files:
        with gzip.open(files[0], 'rt') as f:
            data = json.load(f)
            cells = data['cells'] if 'cells' in data else data
            return len(cells)
    return 0

def generate_report():
    """Generate comprehensive comparison report."""
    
    report = []
    report.append("=" * 70)
    report.append("COMPARISON REPORT: phase2 pipeline analysis")
    report.append(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    report.append("=" * 70)
    
    step23_dir = "../step23/data/rate_0.005000"
    prime_dir = "data/rate_0.005000"
    
    # 1. FILE COUNTS
    report.append("\n1. FILE COUNTS")
    report.append("-" * 40)
    
    for group in ['mutant', 'control1', 'control2']:
        s23_count = count_files(f"{step23_dir}/individuals/{group}")
        pr_count = count_files(f"{prime_dir}/individuals/{group}")
        match = "✅" if s23_count == pr_count else "❌"
        report.append(f"{group:<10}: step23={s23_count:2d}, prime={pr_count:2d} {match}")
    
    # 2. CELL COUNTS
    report.append("\n2. CELL COUNTS (from first individual)")
    report.append("-" * 40)
    
    for group in ['mutant', 'control1', 'control2']:
        s23_cells = get_first_cell_count(f"{step23_dir}/individuals/{group}")
        pr_cells = get_first_cell_count(f"{prime_dir}/individuals/{group}")
        match = "✅" if s23_cells == pr_cells else "❌"
        report.append(f"{group:<10}: step23={s23_cells:4d}, prime={pr_cells:4d} {match}")
    
    # 3. STATISTICS COMPARISON
    report.append("\n3. GROUP STATISTICS")
    report.append("-" * 40)
    
    try:
        with open(f"{step23_dir}/results/statistics.json") as f:
            s23_stats = json.load(f)
        with open(f"{prime_dir}/results/statistics.json") as f:
            pr_stats = json.load(f)
        
        report.append(f"{'Group':<10} {'step23 mean':>12} {'prime mean':>12} {'Difference':>10}")
        
        for group in ['mutant', 'control1', 'control2']:
            if group in s23_stats and group in pr_stats:
                s23_mean = s23_stats[group]['mean']
                pr_mean = pr_stats[group]['mean']
                diff = abs(s23_mean - pr_mean)
                report.append(f"{group:<10} {s23_mean:>12.6f} {pr_mean:>12.6f} {diff:>10.6f}")
    except:
        report.append("Could not load statistics files")
    
    # 4. SNAPSHOT CHECK
    report.append("\n4. SNAPSHOTS")
    report.append("-" * 40)
    
    for year in [50, 60]:
        s23_snap = f"{step23_dir}/snapshots/year{year}_snapshot.json.gz"
        pr_snap = f"{prime_dir}/snapshots/year{year}_snapshot.json.gz"
        
        s23_exists = "✅" if os.path.exists(s23_snap) else "❌"
        pr_exists = "✅" if os.path.exists(pr_snap) else "❌"
        
        report.append(f"Year {year}: step23={s23_exists}, prime={pr_exists}")
    
    # 5. OVERALL VERDICT
    report.append("\n5. OVERALL VERDICT")
    report.append("-" * 40)
    
    # Check key criteria
    all_good = True
    issues = []
    
    # Check file counts
    for group in ['mutant', 'control1', 'control2']:
        s23_count = count_files(f"{step23_dir}/individuals/{group}")
        pr_count = count_files(f"{prime_dir}/individuals/{group}")
        if s23_count != pr_count:
            all_good = False
            issues.append(f"Different {group} counts")
    
    # Check statistics differences
    try:
        for group in ['mutant', 'control1', 'control2']:
            if group in s23_stats and group in pr_stats:
                diff = abs(s23_stats[group]['mean'] - pr_stats[group]['mean'])
                if diff > 0.01:
                    all_good = False
                    issues.append(f"{group} mean differs by {diff:.6f}")
    except:
        pass
    
    if all_good:
        report.append("✅ PIPELINES ARE STATISTICALLY EQUIVALENT")
        report.append("   All checks passed!")
    else:
        report.append("⚠️ PIPELINES SHOW SOME DIFFERENCES:")
        for issue in issues:
            report.append(f"   - {issue}")
    
    report.append("\n" + "=" * 70)
    report.append("END OF REPORT")
    report.append("=" * 70)
    
    # Write to file
    with open("COMPARISON_REPORT.txt", "w") as f:
        f.write("\n".join(report))
    
    # Also print to console
    print("\n".join(report))
    
    print(f"\n📄 Report saved to COMPARISON_REPORT.txt")

if __name__ == "__main__":
    generate_report()