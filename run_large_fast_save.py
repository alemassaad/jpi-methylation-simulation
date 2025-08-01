#!/usr/bin/env python3
"""
Same as baseline but with fast JSON saving.
Parameters: rate=0.5%, n=1000, m=10000, t_max=100
"""

import time
import statistics
import json
import gzip
from cell import History, GENE_SIZE

def save_history_fast(history, filename):
    """Fast save using native JSON serialization and compression."""
    import os
    
    history_dir = "history"
    if not os.path.exists(history_dir):
        os.makedirs(history_dir)
    
    filepath = os.path.join(history_dir, filename + ".json.gz")
    print(f"Saving compressed history to {filepath}")
    
    save_start = time.time()
    
    # Use gzip compression and native JSON serialization
    with gzip.open(filepath, 'wt', encoding='utf-8', compresslevel=1) as f:
        json.dump(history.history, f, separators=(',', ':'))
    
    save_time = time.time() - save_start
    file_size_mb = os.path.getsize(filepath) / (1024 * 1024)
    
    print(f"  Save time: {save_time:.2f} seconds")
    print(f"  File size: {file_size_mb:.2f} MB (compressed)")
    
    return save_time

def main():
    # Simulation parameters
    rate = 0.005  # 0.5%
    n = 1000      # 1000 CpG sites per cell
    m = 10000     # 10,000 cells
    t_max = 100   # 100 years
    gene_size = GENE_SIZE  # Use default (5)
    
    print("=" * 60)
    print("LARGE-SCALE METHYLATION SIMULATION (FAST SAVE)")
    print("=" * 60)
    print(f"Parameters:")
    print(f"  Methylation rate: {rate:.1%}")
    print(f"  CpG sites per cell (n): {n:,}")
    print(f"  Number of cells (m): {m:,}")
    print(f"  Simulation years (t_max): {t_max}")
    print(f"  Gene size: {gene_size}")
    print("=" * 60)
    
    # Start timing
    print("\nStarting simulation...")
    total_start_time = time.time()
    
    # Create History with timing
    print("Creating cell population...")
    creation_start = time.time()
    history = History(m=m, n=n, rate=rate, gene_size=gene_size, t_max=t_max)
    creation_time = time.time() - creation_start
    print(f"  Cell creation time: {creation_time:.2f} seconds")
    
    # Run simulation with timing
    print("\nAging cells...")
    aging_start = time.time()
    
    # Age year by year with progress updates
    for year in range(1, t_max + 1):
        year_start = time.time()
        
        # Age all cells by 1 year
        for cell in history.cells:
            cell.age_1_year()
        
        # Store current state
        history.history[year] = [cell.to_dict() for cell in history.cells]
        
        year_time = time.time() - year_start
        elapsed = time.time() - aging_start
        
        # Progress output every year
        if year <= 10 or year % 5 == 0:  # More frequent updates at start, then every 5 years
            avg_rate = year / elapsed if elapsed > 0 else 0
            remaining = (t_max - year) / avg_rate if avg_rate > 0 else 0
            print(f"  Year {year:3}/{t_max}: {year_time:5.2f}s | Total elapsed: {elapsed:6.1f}s | Est. remaining: {remaining:6.1f}s")
    
    aging_time = time.time() - aging_start
    print(f"\n  Total aging time: {aging_time:.2f} seconds ({aging_time/60:.2f} minutes)")
    
    # Extract and display final statistics with timing
    print("\nCalculating final statistics...")
    stats_start = time.time()
    final_year_data = history.history[t_max]
    methylation_props = [cell['methylation_proportion'] for cell in final_year_data]
    jsd_scores = [cell['jsd'] for cell in final_year_data]
    
    avg_methylation = statistics.mean(methylation_props)
    std_methylation = statistics.stdev(methylation_props)
    min_methylation = min(methylation_props)
    max_methylation = max(methylation_props)
    avg_jsd = statistics.mean(jsd_scores)
    std_jsd = statistics.stdev(jsd_scores)
    stats_time = time.time() - stats_start
    print(f"  Statistics calculation time: {stats_time:.2f} seconds")
    
    print("\nFinal Statistics:")
    print(f"  Average methylation: {avg_methylation:.2%}")
    print(f"  Std dev methylation: {std_methylation:.4f}")
    print(f"  Min methylation: {min_methylation:.2%}")
    print(f"  Max methylation: {max_methylation:.2%}")
    print(f"  Average JSD: {avg_jsd:.4f}")
    print(f"  Std dev JSD: {std_jsd:.4f}")
    
    # Save history with timing (FAST VERSION)
    print("\nSaving results (compressed)...")
    filename = f"baseline_rate_{rate:.1%}_m{m}_n{n}_t{t_max}"
    save_time = save_history_fast(history, filename)
    
    # Total time
    total_time = time.time() - total_start_time
    print("\n" + "=" * 60)
    print("TIMING SUMMARY:")
    print(f"  Cell creation: {creation_time:.2f} seconds")
    print(f"  Aging simulation: {aging_time:.2f} seconds ({aging_time/60:.2f} minutes)")
    print(f"  Statistics: {stats_time:.2f} seconds")
    print(f"  Saving to JSON.gz: {save_time:.2f} seconds")
    print(f"  TOTAL TIME: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    print("=" * 60)
    
    print("\nSimulation complete!")
    print(f"\nTo plot the results, run:")
    print(f"  python plot_history.py {filename}.json.gz")


if __name__ == "__main__":
    main()