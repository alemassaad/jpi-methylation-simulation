#!/usr/bin/env python3
"""
Test that gene rate groups are preserved through the pipeline.
"""

import sys
import os
import json
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish

def test_cell_with_gene_rates():
    """Test Cell creation with gene rate groups."""
    print("Testing Cell with gene rate groups...")
    
    gene_rate_groups = [(5, 0.004), (5, 0.005), (5, 0.006), (5, 0.007)]
    cell = Cell(n=100, gene_rate_groups=gene_rate_groups, gene_size=5)
    
    print(f"  Cell.rate: {cell.rate}")
    print(f"  Cell.gene_rate_groups: {cell.gene_rate_groups}")
    
    assert cell.rate is None, f"Cell.rate should be None, got {cell.rate}"
    assert cell.gene_rate_groups == gene_rate_groups, "Gene rate groups not preserved"
    
    print("  ✓ Cell correctly has gene rate groups\n")
    return cell

def test_cell_from_dict():
    """Test Cell.from_dict with gene rate groups."""
    print("Testing Cell.from_dict with gene rate groups...")
    
    gene_rate_groups = [(5, 0.004), (5, 0.005), (5, 0.006), (5, 0.007)]
    data = {
        'methylated': [0] * 100,
        'cell_JSD': 0.0,
        'age': 30
    }
    
    cell = Cell.from_dict(data, rate=None, gene_rate_groups=gene_rate_groups, gene_size=5)
    
    print(f"  Cell.rate: {cell.rate}")
    print(f"  Cell.gene_rate_groups: {cell.gene_rate_groups}")
    
    assert cell.rate is None, f"Cell.rate should be None, got {cell.rate}"
    assert cell.gene_rate_groups == gene_rate_groups, "Gene rate groups not preserved"
    
    print("  ✓ Cell.from_dict correctly preserves gene rate groups\n")
    return cell

def test_petri_from_cells():
    """Test PetriDish.from_cells with gene rate groups."""
    print("Testing PetriDish.from_cells with gene rate groups...")
    
    # Create cell with gene rates
    gene_rate_groups = [(5, 0.004), (5, 0.005), (5, 0.006), (5, 0.007)]
    cell = Cell(n=100, gene_rate_groups=gene_rate_groups, gene_size=5)
    
    # Create PetriDish from cell
    petri = PetriDish.from_cells(cell, growth_phase=7, calculate_cell_jsds=True)
    
    print(f"  PetriDish.rate: {petri.rate}")
    print(f"  PetriDish.gene_rate_groups: {petri.gene_rate_groups}")
    
    assert petri.rate is None, f"PetriDish.rate should be None, got {petri.rate}"
    assert petri.gene_rate_groups == gene_rate_groups, f"Gene rate groups not preserved, got {petri.gene_rate_groups}"
    
    print("  ✓ PetriDish.from_cells correctly preserves gene rate groups\n")
    return petri

def test_petri_prepare_for_save():
    """Test that prepare_for_save correctly saves gene rate config."""
    print("Testing PetriDish.prepare_for_save with gene rate groups...")
    
    # Create PetriDish with gene rates
    gene_rate_groups = [(5, 0.004), (5, 0.005), (5, 0.006), (5, 0.007)]
    cell = Cell(n=100, gene_rate_groups=gene_rate_groups, gene_size=5)
    petri = PetriDish.from_cells(cell, growth_phase=7, calculate_cell_jsds=True)
    
    # Prepare for save
    data = petri.prepare_for_save(include_gene_metrics=False)
    
    print(f"  Saved metadata['rate']: {data['metadata'].get('rate')}")
    print(f"  Saved metadata['gene_rate_groups']: {data['metadata'].get('gene_rate_groups')}")
    
    assert data['metadata']['rate'] is None, f"Saved rate should be None, got {data['metadata']['rate']}"
    assert data['metadata']['gene_rate_groups'] == gene_rate_groups, "Gene rate groups not saved correctly"
    
    print("  ✓ prepare_for_save correctly saves gene rate configuration\n")
    return data

def main():
    print("="*60)
    print("Testing Gene Rate Preservation Through Pipeline")
    print("="*60)
    print()
    
    try:
        cell = test_cell_with_gene_rates()
        cell2 = test_cell_from_dict()
        petri = test_petri_from_cells()
        data = test_petri_prepare_for_save()
        
        print("="*60)
        print("All tests PASSED! ✓")
        print("Gene rate groups are correctly preserved at each step")
        print("="*60)
        
        # Print final check
        print("\nFinal verification:")
        print(f"  If rate is None: {data['metadata']['rate'] is None}")
        print(f"  If gene_rate_groups present: {data['metadata']['gene_rate_groups'] is not None}")
        print(f"  Actual groups: {data['metadata']['gene_rate_groups']}")
        
        return 0
    except AssertionError as e:
        print(f"\n✗ Test failed: {e}")
        return 1
    except Exception as e:
        print(f"\n✗ Unexpected error: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())