#!/usr/bin/env python3
"""
Extract a snapshot of cells at a specific year from simulation history.
"""

import json
import gzip
import sys
import os


def extract_year_snapshot(input_file, year, output_file):
    """
    Extract all cells from a specific year and save as snapshot.
    
    Args:
        input_file: Path to input history json.gz file
        year: Year to extract (integer)
        output_file: Path to output snapshot json.gz file
    """
    print(f"Extracting year {year} snapshot from {input_file}...")
    
    # Load the full history
    print("  Loading history data...")
    with gzip.open(input_file, 'rt') as f:
        history = json.load(f)
    
    # Check if year exists
    year_str = str(year)
    if year_str not in history:
        available_years = sorted([int(y) for y in history.keys()])
        print(f"  Error: Year {year} not found in history!")
        print(f"  Available years: {available_years[:5]}...{available_years[-5:]}")
        return False
    
    # Extract the cells for the specified year
    cells_data = history[year_str]
    print(f"  Found {len(cells_data)} cells at year {year}")
    
    # Create snapshot data structure
    # The actual cells might not have n, rate, gene_size stored individually
    # Let's check what fields are available
    metadata = {
        "source_file": os.path.basename(input_file),
        "extracted_year": year,
        "num_cells": len(cells_data)
    }
    
    # Try to extract additional metadata if available
    if cells_data:
        sample_cell = cells_data[0]
        if "n" in sample_cell:
            metadata["n"] = sample_cell["n"]
        if "rate" in sample_cell:
            metadata["rate"] = sample_cell["rate"]
        if "gene_size" in sample_cell:
            metadata["gene_size"] = sample_cell["gene_size"]
        
        # Print available fields for debugging
        print(f"  Cell fields: {list(sample_cell.keys())}")
    
    snapshot = {
        "metadata": metadata,
        "cells": cells_data
    }
    
    # Save snapshot
    print(f"  Saving snapshot to {output_file}...")
    with gzip.open(output_file, 'wt') as f:
        json.dump(snapshot, f, indent=2)
    
    # Report file size
    file_size_mb = os.path.getsize(output_file) / (1024 * 1024)
    print(f"  Saved successfully! File size: {file_size_mb:.2f} MB")
    
    return True


def main():
    if len(sys.argv) != 4:
        print("Usage: python extract_snapshot.py <input_history.json.gz> <year> <output_snapshot.json.gz>")
        print("Example: python extract_snapshot.py ../history/simulation_rate_0.005_m10000_n1000_t100.json.gz 50 snapshots/year50_snapshot.json.gz")
        sys.exit(1)
    
    input_file = sys.argv[1]
    year = int(sys.argv[2])
    output_file = sys.argv[3]
    
    # Create output directory if needed
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")
    
    # Extract snapshot
    success = extract_year_snapshot(input_file, year, output_file)
    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()