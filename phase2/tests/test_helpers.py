#!/usr/bin/env python3
"""
Helper functions for testing uniform mixing implementation.
Creates lightweight mock data for fast tests.
"""

import sys
import os
import random
import json
import gzip
import copy

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish


def create_mock_cells(n_cells=100, n_sites=50, year=60, seed=42):
    """Create lightweight mock cells for testing."""
    random.seed(seed)
    
    cells = []
    for i in range(n_cells):
        cell = Cell(n=n_sites, rate=0.005, gene_size=5, seed=seed + i)
        
        # Create unique methylation pattern for each cell
        # Use a gradient to make cells distinguishable
        methylation_rate = 0.1 + (0.3 * i / n_cells)
        for j in range(n_sites):
            if random.random() < methylation_rate:
                cell.cpg_sites[j] = 1
                cell.methylated[j] = True
        
        cell.age = year
        # cell.cell_JSD is already calculated in Cell.__init__
        cells.append(cell)
    
    return cells


def create_mock_petri_dish(n_cells=10, n_sites=50, rate=0.005, seed=42):
    """Create simple PetriDish for testing."""
    petri = PetriDish(rate=rate, n=n_sites, gene_size=5, seed=seed)
    petri.cells = create_mock_cells(n_cells, n_sites, seed=seed)
    petri.year = 60
    return petri


def create_mock_simulation_file(filepath, years=10, growth_phase=3, n_sites=50, seed=42):
    """Create minimal simulation.json.gz for testing."""
    random.seed(seed)
    
    simulation_data = {
        'parameters': {
            'rate': 0.005,
            'n': n_sites,
            'gene_size': 5,
            'growth_phase': growth_phase,
            'years': years,
            'seed': seed
        },
        'history': {}
    }
    
    # Create simple growth pattern
    n_cells = 1
    for year in range(years + 1):
        if year <= growth_phase:
            n_cells = min(2 ** year, 2 ** growth_phase)
        else:
            # Homeostasis - vary around target
            n_cells = 2 ** growth_phase + random.randint(-2, 2)
        
        # Create cells for this year
        year_cells = []
        for i in range(n_cells):
            cell_data = {
                'methylated': [random.random() < 0.2 for _ in range(n_sites)],
                'age': year,
                'cell_jsd': random.random() * 0.1
            }
            year_cells.append(cell_data)
        
        simulation_data['history'][str(year)] = year_cells
    
    # Save as compressed JSON
    with gzip.open(filepath, 'wt') as f:
        json.dump(simulation_data, f)
    
    return filepath


def cells_are_identical(cell1, cell2):
    """Check if two cells have identical methylation patterns."""
    if len(cell1.methylated) != len(cell2.methylated):
        return False
    
    for i in range(len(cell1.methylated)):
        if cell1.methylated[i] != cell2.methylated[i]:
            return False
    
    return True


def get_cell_fingerprint(cell):
    """Create a unique fingerprint for a cell based on its methylation pattern."""
    return tuple(cell.methylated)


def count_matching_cells(cells1, cells2):
    """Count how many cells match between two lists at corresponding positions."""
    matches = 0
    for i in range(min(len(cells1), len(cells2))):
        if cells_are_identical(cells1[i], cells2[i]):
            matches += 1
    return matches


def verify_no_duplicates(cells):
    """Check that no cell appears twice in the list."""
    fingerprints = [get_cell_fingerprint(cell) for cell in cells]
    return len(fingerprints) == len(set(fingerprints))