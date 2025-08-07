#!/usr/bin/env python3
"""
Create test data for step3 pipeline testing.
This creates small test control lineages to match the test mutant lineages from step2.
"""

import json
import gzip
import os
import sys
sys.path.append('../..')

from step1.cell import Cell
import random

def create_test_control_lineages():
    """Create 10 test control lineages with 8 cells each."""
    
    # Create test control lineages directory
    os.makedirs('../../step2/test_lineages_control', exist_ok=True)
    
    # Set seed for reproducibility
    random.seed(123)
    
    print("Creating test control lineages...")
    
    for i in range(10):
        # Create 8 cells with varying JSD
        cells = []
        for j in range(8):
            # Create a cell with some methylation
            cell = Cell(n=15, rate=0.01, gene_size=5)
            
            # Age it a bit to get some methylation
            for _ in range(random.randint(30, 70)):
                cell.age_1_year()
            
            cells.append(cell.to_dict())
        
        # Create lineage data structure
        lineage_data = {
            "metadata": {
                "lineage_id": i,
                "lineage_type": "control",
                "source_decile": 0,  # 0 for uniform sampling
                "within_decile_index": i,
                "source_year": 50,
                "final_year": 60,
                "generations": 3,  # Test uses 3 years
                "final_cell_count": 8,
                "expected_cell_count": 8,
                "source_cell": cells[0]  # Just use first cell as source
            },
            "cells": cells
        }
        
        # Save lineage
        filename = f"control_{i:02d}.json.gz"
        filepath = os.path.join('../../step2/test_lineages_control', filename)
        
        with gzip.open(filepath, 'wt') as f:
            json.dump(lineage_data, f, indent=2)
        
        print(f"  Created {filepath}")
    
    print("\nTest control lineages created!")


def create_test_year60_snapshot():
    """Create a small test year 60 snapshot."""
    
    print("\nCreating test year 60 snapshot...")
    
    # Create 100 cells at year 60
    cells = []
    random.seed(456)
    
    for i in range(100):
        cell = Cell(n=15, rate=0.01, gene_size=5)
        
        # Age to year 60
        for _ in range(60):
            cell.age_1_year()
        
        cells.append(cell.to_dict())
    
    # Create snapshot
    snapshot_data = {
        "metadata": {
            "source_file": "test_simulation.json.gz",
            "extracted_year": 60,
            "num_cells": len(cells)
        },
        "cells": cells
    }
    
    # Save
    os.makedirs('../data/test_snapshots', exist_ok=True)
    filepath = '../data/test_snapshots/test_year60_snapshot.json.gz'
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(snapshot_data, f, indent=2)
    
    print(f"  Created {filepath}")
    print(f"  Contains {len(cells)} cells")


if __name__ == "__main__":
    create_test_control_lineages()
    create_test_year60_snapshot()
    print("\nTest data creation complete!")