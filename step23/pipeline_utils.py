#!/usr/bin/env python3
"""
Utility functions for the step23 unified pipeline.
Refactored from step2 and step3 utilities.
"""

import json
import gzip
import numpy as np
import random
import copy
import os
import sys
sys.path.append('..')  # To import from parent directory

from step1.cell import Cell, GENE_SIZE, BASELINE_METHYLATION_DISTRIBUTION


def load_simulation(filepath):
    """Load full simulation history from json.gz file."""
    print(f"Loading simulation from {filepath}...")
    with gzip.open(filepath, 'rt') as f:
        history = json.load(f)
    years = sorted([int(y) for y in history.keys()])
    print(f"  Loaded simulation with years {years[0]} to {years[-1]}")
    return history


def extract_snapshot(simulation_file, year):
    """
    Extract cells from a specific year of the simulation.
    
    Returns:
        List of cell dictionaries
    """
    print(f"\nExtracting year {year} snapshot...")
    
    # Load simulation
    with gzip.open(simulation_file, 'rt') as f:
        history = json.load(f)
    
    # Extract specified year
    year_str = str(year)
    if year_str not in history:
        available_years = sorted([int(y) for y in history.keys()])
        raise ValueError(f"Year {year} not found. Available: {available_years}")
    
    cells = history[year_str]
    print(f"  Extracted {len(cells)} cells from year {year}")
    
    return cells


def save_snapshot(cells, filepath):
    """Save snapshot to compressed JSON file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    snapshot_data = {
        "metadata": {
            "num_cells": len(cells),
            "extracted_year": cells[0]['age'] if cells else 0
        },
        "cells": cells
    }
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(snapshot_data, f, indent=2)
    
    print(f"  Saved snapshot to {filepath}")


def load_individual(filepath):
    """Load individual cells from json.gz file."""
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    # Handle both formats: direct cell list or with metadata
    if isinstance(data, list):
        return data
    elif isinstance(data, dict) and 'cells' in data:
        return data['cells']
    else:
        raise ValueError(f"Unknown individual file format in {filepath}")


def save_individual(cells, filepath, metadata=None):
    """Save individual cells to compressed JSON file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # Calculate statistics
    jsd_values = [cell['jsd'] for cell in cells]
    
    data = {
        "metadata": {
            "num_cells": len(cells),
            "mean_jsd": float(np.mean(jsd_values)),
            "std_jsd": float(np.std(jsd_values)),
            "min_jsd": float(np.min(jsd_values)),
            "max_jsd": float(np.max(jsd_values))
        },
        "cells": cells
    }
    
    # Add custom metadata if provided
    if metadata:
        data["metadata"].update(metadata)
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(data, f, indent=2)


def sample_by_deciles(cells, cells_per_decile=3, seed=None):
    """
    Sample cells from each JSD decile.
    
    Returns:
        List of cell dictionaries (30 cells total for default)
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    print(f"  Sampling {cells_per_decile} cells from each JSD decile...")
    
    # Sort cells by JSD
    cells_with_jsd = [(cell, cell['jsd']) for cell in cells]
    cells_with_jsd.sort(key=lambda x: x[1])
    
    # Calculate decile boundaries
    n_cells = len(cells_with_jsd)
    decile_size = n_cells // 10
    
    sampled_cells = []
    
    for decile in range(10):
        # Define range for this decile
        start_idx = decile * decile_size
        end_idx = n_cells if decile == 9 else (decile + 1) * decile_size
        
        # Get cells in this decile
        decile_cells = cells_with_jsd[start_idx:end_idx]
        
        # Sample from this decile
        if len(decile_cells) >= cells_per_decile:
            sampled_indices = random.sample(range(len(decile_cells)), cells_per_decile)
            for idx in sampled_indices:
                cell, _ = decile_cells[idx]
                # Add decile info
                cell_copy = cell.copy()
                cell_copy['source_decile'] = decile + 1
                sampled_cells.append(cell_copy)
    
    print(f"    Sampled {len(sampled_cells)} cells total")
    return sampled_cells


def sample_uniform(cells, n_cells=30, seed=None):
    """
    Sample cells uniformly from the population.
    
    Returns:
        List of cell dictionaries
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    print(f"  Sampling {n_cells} cells uniformly...")
    
    if n_cells > len(cells):
        raise ValueError(f"Cannot sample {n_cells} from {len(cells)} cells")
    
    sampled_indices = random.sample(range(len(cells)), n_cells)
    sampled_cells = []
    
    for idx in sampled_indices:
        cell_copy = cells[idx].copy()
        cell_copy['source_decile'] = 0  # 0 indicates uniform sampling
        sampled_cells.append(cell_copy)
    
    print(f"    Sampled {len(sampled_cells)} cells")
    return sampled_cells


def dict_to_cell(cell_dict):
    """Convert dictionary representation to Cell object."""
    cell = Cell(
        n=len(cell_dict['cpg_sites']),
        rate=cell_dict['rate'],
        gene_size=cell_dict.get('gene_size', GENE_SIZE)
    )
    cell.cpg_sites = cell_dict['cpg_sites'].copy()
    cell.age_val = cell_dict['age']
    
    # Set computed values
    cell.methylation_proportion = cell_dict['methylation_proportion']
    cell.methylation_distribution = cell_dict['methylation_distribution'].copy()
    cell.JSD = cell_dict['jsd']
    
    return cell


def grow_individual_inplace(filepath, years, rate=None):
    """
    Grow an individual by cell division and aging, modifying the file in place.
    Each year: cells divide (duplicate), then all cells age.
    
    Args:
        filepath: Path to individual file
        years: Number of years to grow
        rate: Methylation rate (if None, uses cell's rate)
    
    Returns:
        Final number of cells
    """
    print(f"    Growing individual for {years} years...")
    
    for year in range(years):
        # Load current cells
        cells_dict = load_individual(filepath)
        
        # Convert to Cell objects
        cells = [dict_to_cell(cd) for cd in cells_dict]
        
        # Division phase: each cell creates a copy
        new_cells = []
        for cell in cells:
            daughter1 = copy.deepcopy(cell)
            daughter2 = copy.deepcopy(cell)
            new_cells.extend([daughter1, daughter2])
        
        # Aging phase: all cells age by 1 year
        for cell in new_cells:
            cell.age_1_year()
        
        # Convert back to dictionaries
        cells_dict = [cell.to_dict() for cell in new_cells]
        
        # Save back to same file
        save_individual(cells_dict, filepath)
        
        if year == 0 or year == years - 1:
            print(f"      Year {year + 1}: {len(cells_dict)} cells")
    
    return len(cells_dict)


def mix_with_year60_inplace(filepath, year60_cells, mix_ratio=80, seed=None):
    """
    Add year 60 cells to an existing individual to create mixture.
    
    Args:
        filepath: Path to individual file
        year60_cells: List of all year 60 cells to sample from
        mix_ratio: Percentage of year 60 cells in final mixture (default 80%)
        seed: Random seed for sampling
    
    Returns:
        Total number of cells after mixing
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Load grown cells
    grown_cells = load_individual(filepath)
    n_grown = len(grown_cells)
    
    # Calculate total cells needed
    # If grown cells are (100-mix_ratio)%, calculate total
    percentage_grown = 100 - mix_ratio
    n_total = int(n_grown * 100 / percentage_grown)
    n_to_add = n_total - n_grown
    
    print(f"    Mixing: {n_grown} grown + {n_to_add} from year 60 = {n_total} total")
    
    # Sample cells from year 60
    if n_to_add > len(year60_cells):
        raise ValueError(f"Need {n_to_add} year 60 cells but only {len(year60_cells)} available")
    
    sampled_indices = random.sample(range(len(year60_cells)), n_to_add)
    added_cells = [year60_cells[i].copy() for i in sampled_indices]
    
    # Mark sources
    for cell in grown_cells:
        cell['source'] = 'grown'
    for cell in added_cells:
        cell['source'] = 'year60'
    
    # Combine all cells
    all_cells = grown_cells + added_cells
    
    # Shuffle to mix
    random.shuffle(all_cells)
    
    # Save back with updated metadata
    metadata = {
        "n_grown": n_grown,
        "n_added": n_to_add,
        "mix_ratio": mix_ratio
    }
    save_individual(all_cells, filepath, metadata)
    
    return len(all_cells)


def create_control2_individual(year60_cells, n_cells, seed=None):
    """
    Create a control2 individual by sampling from year 60 snapshot.
    
    Args:
        year60_cells: List of all year 60 cells
        n_cells: Number of cells to sample
        seed: Random seed
    
    Returns:
        List of sampled cells
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    if n_cells > len(year60_cells):
        raise ValueError(f"Cannot sample {n_cells} from {len(year60_cells)} cells")
    
    sampled_indices = random.sample(range(len(year60_cells)), n_cells)
    sampled_cells = []
    
    for idx in sampled_indices:
        cell = year60_cells[idx].copy()
        cell['source'] = 'year60_control2'
        sampled_cells.append(cell)
    
    return sampled_cells