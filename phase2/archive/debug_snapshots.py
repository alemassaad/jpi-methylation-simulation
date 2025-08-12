#!/usr/bin/env python3
"""
Debug why snapshots are different between runs.
"""

import os
import json
import gzip
import hashlib

def analyze_snapshot(filepath):
    """Analyze a snapshot file."""
    if not os.path.exists(filepath):
        print(f"  File not found: {filepath}")
        return None
    
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    # Get metadata
    metadata = data.get('metadata', {})
    cells = data.get('cells', data if isinstance(data, list) else [])
    
    # Get first few cells' JSDs for comparison
    first_jsds = [c['jsd'] for c in cells[:5]]
    
    # Hash the actual cell data (not metadata)
    cells_str = json.dumps(cells, sort_keys=True)
    cells_hash = hashlib.md5(cells_str.encode()).hexdigest()[:8]
    
    return {
        'metadata': metadata,
        'n_cells': len(cells),
        'first_jsds': first_jsds,
        'cells_hash': cells_hash,
        'file_size': os.path.getsize(filepath)
    }

def main():
    print("DEBUGGING SNAPSHOT DIFFERENCES")
    print("="*50)
    
    # Check snapshots from test runs
    for run in ['run1', 'run2']:
        print(f"\n{run.upper()}:")
        snap_path = f"data_test_{run}/rate_0.010000/snapshots/year50_snapshot.json.gz"
        
        info = analyze_snapshot(snap_path)
        if info:
            print(f"  Metadata: {info['metadata']}")
            print(f"  Cells: {info['n_cells']}")
            print(f"  First JSDs: {info['first_jsds'][:3]}")
            print(f"  Cell data hash: {info['cells_hash']}")
            print(f"  File size: {info['file_size']} bytes")
    
    # Compare the actual extraction process
    print("\n" + "="*50)
    print("TESTING EXTRACTION DETERMINISM")
    print("="*50)
    
    # Extract twice from same source
    simulation = "../phase1/data/simulation_rate_0.010000_g3_m22_n100_t60_seed42.json.gz"
    if os.path.exists(simulation):
        print(f"\nExtracting year 50 twice from same simulation...")
        
        import sys
        sys.path.append('..')
        sys.path.append('../phase1')
        from pipeline_utils import load_snapshot_as_cells, save_snapshot_cells
        
        # Extract 1
        cells1 = load_snapshot_as_cells(simulation, 50)
        save_snapshot_cells(cells1, "./test_snap1.json.gz")
        
        # Extract 2  
        cells2 = load_snapshot_as_cells(simulation, 50)
        save_snapshot_cells(cells2, "./test_snap2.json.gz")
        
        # Compare
        info1 = analyze_snapshot("test_snap1.json.gz")
        info2 = analyze_snapshot("test_snap2.json.gz")
        
        if info1 and info2:
            if info1['cells_hash'] == info2['cells_hash']:
                print("✅ Extractions are identical!")
            else:
                print("❌ Extractions differ!")
                print(f"  Hash 1: {info1['cells_hash']}")
                print(f"  Hash 2: {info2['cells_hash']}")
        
        # Clean up
        os.remove("test_snap1.json.gz")
        os.remove("test_snap2.json.gz")
    else:
        print(f"Simulation not found: {simulation}")

if __name__ == "__main__":
    main()