#!/usr/bin/env python3
"""
Extract year 60 cells from the original simulation.
This creates a snapshot similar to step2's extract_snapshot.py but specifically for year 60.
"""

import json
import gzip
import sys
import os

def extract_year_60(input_filepath, output_filepath):
    """Extract year 60 data from simulation history file."""
    print(f"Loading simulation data from {input_filepath}...")
    
    with gzip.open(input_filepath, 'rt') as f:
        history_data = json.load(f)
    
    # Extract year 60 data
    year_60_key = "60"
    if year_60_key not in history_data:
        print(f"Error: Year 60 not found in simulation data!")
        print(f"Available years: {sorted([int(k) for k in history_data.keys()])}")
        return False
    
    year_60_cells = history_data[year_60_key]
    print(f"Found {len(year_60_cells)} cells at year 60")
    
    # Create snapshot data structure
    snapshot_data = {
        "metadata": {
            "source_file": os.path.basename(input_filepath),
            "extracted_year": 60,
            "num_cells": len(year_60_cells),
            "description": "Year 60 snapshot from original simulation"
        },
        "cells": year_60_cells
    }
    
    # Save to compressed JSON
    print(f"Saving snapshot to {output_filepath}...")
    with gzip.open(output_filepath, 'wt') as f:
        json.dump(snapshot_data, f)
    
    print(f"Successfully saved {len(year_60_cells)} cells to {output_filepath}")
    return True


if __name__ == "__main__":
    # Default paths
    input_file = "../history/simulation_rate_0.005_m10000_n1000_t100.json.gz"
    output_file = "year60_original_snapshot.json.gz"
    
    # Allow command line override
    if len(sys.argv) > 1:
        input_file = sys.argv[1]
    if len(sys.argv) > 2:
        output_file = sys.argv[2]
    
    extract_year_60(input_file, output_file)