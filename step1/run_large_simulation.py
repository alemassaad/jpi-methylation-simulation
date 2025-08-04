#!/usr/bin/env python3
"""
Large-scale methylation simulation with configurable parameters and fast compressed saving.
"""

import argparse
import time
import statistics
import json
import gzip
from cell import History, GENE_SIZE

def save_history_fast(history, filename):
    """Fast save using native JSON serialization and compression."""
    import os
    
    history_dir = "data"
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

def parse_arguments():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="Run large-scale methylation simulation with configurable parameters",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Simulation parameters
    parser.add_argument('-r', '--rate', type=float, default=0.005,
                        help='Methylation rate per site per year (0.01 = 1%%, 0.005 = 0.5%%)')
    parser.add_argument('-m', '--cells', type=int, default=10000,
                        help='Number of cells in the population')
    parser.add_argument('-n', '--sites', type=int, default=1000,
                        help='Number of CpG sites per cell')
    parser.add_argument('-t', '--years', type=int, default=100,
                        help='Number of years to simulate')
    parser.add_argument('-g', '--gene-size', type=int, default=GENE_SIZE,
                        help='Size of genes (sites per gene)')
    
    # Optional parameters
    parser.add_argument('-o', '--output', type=str, default=None,
                        help='Output filename (auto-generated if not specified)')
    parser.add_argument('--seed', type=int, default=None,
                        help='Random seed for reproducibility')
    parser.add_argument('--no-compress', action='store_true',
                        help='Save uncompressed JSON instead of compressed')
    
    return parser.parse_args()

def main():
    # Parse command-line arguments
    args = parse_arguments()
    
    # Set random seed if provided
    if args.seed is not None:
        import random
        random.seed(args.seed)
        print(f"Random seed set to: {args.seed}")
    
    # Extract parameters from args
    rate = args.rate
    n = args.sites
    m = args.cells
    t_max = args.years
    gene_size = args.gene_size
    
    # Validate parameters
    if n % gene_size != 0:
        print(f"Error: Number of sites ({n}) must be divisible by gene size ({gene_size})")
        return 1
    
    print("=" * 60)
    print("LARGE-SCALE METHYLATION SIMULATION (FAST SAVE)")
    print("=" * 60)
    print(f"Parameters:")
    print(f"  Methylation rate: {rate:.3%}")
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
    
    # Save history with timing
    if args.output:
        filename = args.output
    else:
        filename = f"simulation_rate_{rate:.6f}_m{m}_n{n}_t{t_max}"
    
    if args.no_compress:
        print("\nSaving results (uncompressed)...")
        save_start = time.time()
        history.save_history(filename, directory="data")
        save_time = time.time() - save_start
        extension = ".json"
    else:
        print("\nSaving results (compressed)...")
        save_time = save_history_fast(history, filename)
        extension = ".json.gz"
    
    # Total time
    total_time = time.time() - total_start_time
    print("\n" + "=" * 60)
    print("TIMING SUMMARY:")
    print(f"  Cell creation: {creation_time:.2f} seconds")
    print(f"  Aging simulation: {aging_time:.2f} seconds ({aging_time/60:.2f} minutes)")
    print(f"  Statistics: {stats_time:.2f} seconds")
    print(f"  Saving to {extension}: {save_time:.2f} seconds")
    print(f"  TOTAL TIME: {total_time:.2f} seconds ({total_time/60:.2f} minutes)")
    print("=" * 60)
    
    print("\nSimulation complete!")
    print(f"\nTo plot the results, run:")
    print(f"  python plot_history.py data/{filename}{extension}")


if __name__ == "__main__":
    main()