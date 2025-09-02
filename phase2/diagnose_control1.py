#!/usr/bin/env python3
"""
Diagnose why control1 has uniform rate instead of gene rate groups.
"""

import sys
import os
import json
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish
from core.pipeline_utils import load_snapshot_cells

def diagnose():
    """Diagnose the control1 issue."""
    print("="*60)
    print("Diagnosing Control1 Rate Issue")
    print("="*60)
    print()
    
    # Path to problematic data
    base_path = "data/gene_rates_5x0.00400_5x0.00450_5x0.00500_5x0.00550-grow9-sites100-years50/snap30to50-growth6-quant4x2-mix70un-seed42-20250902030844"
    
    # Check snapshot
    print("1. Checking Year 30 Snapshot:")
    print("-" * 40)
    snapshot_path = os.path.join(base_path, "snapshots/year30_snapshot.json")
    
    with open(snapshot_path, 'r') as f:
        snapshot_data = json.load(f)
    
    print(f"  Metadata rate: {snapshot_data['metadata'].get('rate')}")
    print(f"  Metadata gene_rate_groups: {snapshot_data['metadata'].get('gene_rate_groups')}")
    
    # Load cells from snapshot
    print("\n2. Loading Cells from Snapshot:")
    print("-" * 40)
    cells = load_snapshot_cells(snapshot_path)
    
    if cells:
        first_cell = cells[0]
        print(f"  First cell rate: {first_cell.rate}")
        print(f"  First cell gene_rate_groups: {first_cell.gene_rate_groups}")
        print(f"  First cell n: {first_cell.n}")
        print(f"  First cell gene_size: {first_cell.gene_size}")
    
    # Test PetriDish creation
    print("\n3. Testing PetriDish Creation from Cell:")
    print("-" * 40)
    if cells:
        petri = PetriDish.from_cells(cells[0], growth_phase=7)
        print(f"  PetriDish.rate: {petri.rate}")
        print(f"  PetriDish.gene_rate_groups: {petri.gene_rate_groups}")
        
        # Test prepare_for_save
        print("\n4. Testing prepare_for_save:")
        print("-" * 40)
        save_data = petri.prepare_for_save(include_gene_metrics=False)
        print(f"  Saved rate: {save_data['metadata'].get('rate')}")
        print(f"  Saved gene_rate_groups: {save_data['metadata'].get('gene_rate_groups')}")
    
    # Check actual control1 file
    print("\n5. Checking Actual Control1 Individual:")
    print("-" * 40)
    control1_path = os.path.join(base_path, "individuals/control1/individual_01.json")
    
    with open(control1_path, 'r') as f:
        # Just read the metadata section
        content = f.read()
        # Find metadata section
        import re
        match = re.search(r'"metadata":\s*{[^}]*"rate":\s*([^,}]+)', content)
        if match:
            print(f"  Found rate in file: {match.group(1)}")
        match = re.search(r'"metadata":\s*{[^}]*"gene_rate_groups":\s*([^,}]+)', content)
        if match:
            print(f"  Found gene_rate_groups in file: {match.group(1)}")
    
    print("\n" + "="*60)
    print("Diagnosis Complete")
    print("="*60)

if __name__ == "__main__":
    diagnose()