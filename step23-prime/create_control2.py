#!/usr/bin/env python3
"""
Create control2 individuals (pure year 60 cells).
"""

import os
import sys

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'step1-prime'))

from pipeline_utils import load_snapshot_cells, create_pure_snapshot_petri, save_petri_dish

def main():
    # Parameters matching the original run
    rate = 0.005
    seed = 42
    expected_individuals = 30
    expected_final_cells = 5120  # Same as mixed individuals
    
    # Paths
    base_dir = "data/rate_0.005000"
    year60_path = os.path.join(base_dir, "snapshots/year60_snapshot.json.gz")
    control2_dir = os.path.join(base_dir, "individuals/control2")
    
    print("Creating Control2 individuals...")
    print(f"  Loading year 60 snapshot from {year60_path}...")
    year60_cells = load_snapshot_cells(year60_path)
    print(f"  Loaded {len(year60_cells)} cells")
    
    print(f"  Creating {expected_individuals} control2 individuals...")
    print(f"    Each with {expected_final_cells} pure year 60 cells")
    
    for i in range(expected_individuals):
        print(f"    Creating individual {i+1}/{expected_individuals}")
        
        # Create PetriDish with pure year 60 cells
        petri = create_pure_snapshot_petri(year60_cells, n_cells=expected_final_cells,
                                          rate=rate, seed=seed + 300 + i)
        
        # Add metadata
        if not hasattr(petri, 'metadata'):
            petri.metadata = {}
        petri.metadata.update({
            'individual_id': i,
            'individual_type': 'control2',
            'source': 'pure_year60',
            'year': 60
        })
        
        # Save
        filepath = os.path.join(control2_dir, f"individual_{i:02d}.json.gz")
        save_petri_dish(petri, filepath)
    
    print(f"  Created and saved {expected_individuals} control2 individuals")

if __name__ == "__main__":
    main()