#!/usr/bin/env python3
"""
Utility functions for phase2 pipeline using PetriDish and Cell classes.
"""

import json
import gzip
import numpy as np
import random
import os
import sys
import glob
import copy
from typing import List, Dict, Tuple, Optional, Any

# Add parent directory to path to import from phase1
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import Cell, PetriDish, GENE_SIZE, BASELINE_METHYLATION_DISTRIBUTION


def dict_to_cell(cell_dict: Dict[str, Any]) -> Cell:
    """
    Convert dictionary from JSON to Cell object.
    
    Args:
        cell_dict: Dictionary with cell data
    
    Returns:
        Cell: Initialized Cell object
    """
    # Create cell with proper parameters
    cell = Cell(
        n=len(cell_dict['cpg_sites']),
        rate=cell_dict['rate'],
        gene_size=cell_dict.get('gene_size', GENE_SIZE)
    )
    
    # Set attributes directly
    cell.cpg_sites = cell_dict['cpg_sites']
    cell.age = cell_dict['age']
    cell.methylation_proportion = cell_dict['methylation_proportion']
    cell.methylation_distribution = cell_dict['methylation_distribution']
    cell.JSD = cell_dict['jsd']
    
    return cell


def load_snapshot_as_cells(simulation_file: str, year: int) -> List[Cell]:
    """
    Load a year from simulation and convert to Cell objects.
    
    Args:
        simulation_file: Path to phase1 simulation
        year: Year to extract (50 or 60)
    
    Returns:
        List[Cell]: List of Cell objects
    """
    print(f"Loading year {year} snapshot from {simulation_file}...")
    
    with gzip.open(simulation_file, 'rt') as f:
        history = json.load(f)
    
    year_str = str(year)
    if year_str not in history:
        available_years = sorted([int(y) for y in history.keys()])
        raise ValueError(f"Year {year} not found. Available: {available_years}")
    
    # Convert all cell dicts to Cell objects
    cell_dicts = history[year_str]
    cells = [dict_to_cell(cd) for cd in cell_dicts]
    
    print(f"  Loaded {len(cells)} cells from year {year}")
    return cells


def save_snapshot_cells(cells: List[Cell], filepath: str) -> None:
    """
    Save a list of Cell objects as a snapshot.
    
    Args:
        cells: List of Cell objects
        filepath: Output path for compressed JSON
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    snapshot_data = {
        "metadata": {
            "num_cells": len(cells),
            "year": cells[0].age if cells else 0
        },
        "cells": [cell.to_dict() for cell in cells]
    }
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(snapshot_data, f, indent=2)
    
    print(f"  Saved {len(cells)} cells to {filepath}")


def load_snapshot_cells(filepath: str) -> List[Cell]:
    """
    Load a snapshot file and convert to Cell objects.
    
    Args:
        filepath: Path to snapshot file
    
    Returns:
        List[Cell]: List of Cell objects
    """
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    cell_dicts = data['cells']
    cells = [dict_to_cell(cd) for cd in cell_dicts]
    
    return cells


def save_petri_dish(petri: PetriDish, filepath: str, metadata: Optional[Dict] = None) -> None:
    """
    Save a PetriDish to compressed JSON file.
    
    Args:
        petri: PetriDish object
        filepath: Output path
        metadata: Extra metadata dict
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # Calculate statistics
    jsd_values = [cell.JSD for cell in petri.cells]
    
    data = {
        "metadata": {
            "num_cells": len(petri.cells),
            "mean_jsd": float(np.mean(jsd_values)) if jsd_values else 0.0,
            "std_jsd": float(np.std(jsd_values)) if jsd_values else 0.0,
            "min_jsd": float(np.min(jsd_values)) if jsd_values else 0.0,
            "max_jsd": float(np.max(jsd_values)) if jsd_values else 0.0,
            "median_jsd": float(np.median(jsd_values)) if jsd_values else 0.0,
            "year": petri.year,
            "growth_phase": petri.growth_phase,
            "rate": petri.rate
        },
        "cells": [cell.to_dict() for cell in petri.cells]
    }
    
    # Add custom metadata if provided
    if metadata:
        data["metadata"].update(metadata)
    
    # Add petri metadata if it exists
    if hasattr(petri, 'metadata'):
        data["metadata"].update(petri.metadata)
    
    with gzip.open(filepath, 'wt') as f:
        json.dump(data, f, indent=2)


def load_petri_dish(filepath: str) -> PetriDish:
    """
    Load a PetriDish from file.
    
    Returns:
        PetriDish: Reconstructed object with cells and metadata
    """
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    # Get first cell to extract parameters
    if data['cells']:
        first_cell = data['cells'][0]
        rate = first_cell.get('rate', 0.005)
        n = len(first_cell['cpg_sites'])
        gene_size = first_cell.get('gene_size', GENE_SIZE)
    else:
        # Defaults if no cells
        rate = 0.005
        n = 1000
        gene_size = GENE_SIZE
    
    # Create PetriDish (starts with 1 cell, we'll replace)
    petri = PetriDish(rate=rate, n=n, gene_size=gene_size, seed=None)
    
    # Replace cells with loaded cells
    petri.cells = [dict_to_cell(cd) for cd in data['cells']]
    
    # Restore year if available
    if 'year' in data['metadata']:
        petri.year = data['metadata']['year']
    
    # Add metadata attribute if not exists
    if not hasattr(petri, 'metadata'):
        petri.metadata = {}
    
    # Restore metadata
    petri.metadata = data.get('metadata', {})
    
    return petri


def load_all_petri_dishes(directory: str) -> List[PetriDish]:
    """
    Load all PetriDish objects from a directory.
    
    Args:
        directory: Directory containing .json.gz files
    
    Returns:
        List of PetriDish objects
    """
    files = sorted(glob.glob(os.path.join(directory, "*.json.gz")))
    dishes = []
    
    for filepath in files:
        try:
            petri = load_petri_dish(filepath)
            dishes.append(petri)
        except Exception as e:
            print(f"  Warning: Could not load {filepath}: {e}")
    
    return dishes


def sample_by_quantiles(cells: List[Cell], n_quantiles: int = 10, 
                        cells_per_quantile: int = 3, seed: int = 42) -> List[Tuple[Cell, int]]:
    """
    Sample cells from quantiles based on JSD.
    
    Args:
        cells: List[Cell] - all cells to sample from
        n_quantiles: Number of quantiles (4=quartiles, 10=deciles)
        cells_per_quantile: Cells to sample per quantile
        seed: Random seed
    
    Returns:
        List[Tuple[Cell, int]]: List of (cell, quantile_index) tuples
    """
    random.seed(seed)
    np.random.seed(seed)  # Also seed numpy for full reproducibility
    
    print(f"  Sampling {cells_per_quantile} cells from each of {n_quantiles} quantiles...")
    
    # Sort cells by JSD
    sorted_cells = sorted(cells, key=lambda c: c.JSD)
    
    # Calculate quantile boundaries
    n_cells = len(sorted_cells)
    quantile_size = n_cells // n_quantiles
    
    sampled = []
    
    for q in range(n_quantiles):
        # Define quantile boundaries
        start_idx = q * quantile_size
        if q == n_quantiles - 1:
            # Last quantile includes any remainder
            end_idx = n_cells
        else:
            end_idx = (q + 1) * quantile_size
        
        # Get cells in this quantile
        quantile_cells = sorted_cells[start_idx:end_idx]
        
        # Sample from this quantile
        if len(quantile_cells) >= cells_per_quantile:
            sample_indices = random.sample(range(len(quantile_cells)), cells_per_quantile)
            for idx in sample_indices:
                # Deep copy the cell to avoid reference issues
                cell_copy = copy.deepcopy(quantile_cells[idx])
                sampled.append((cell_copy, q))
        else:
            # If not enough cells, take all
            for cell in quantile_cells:
                cell_copy = copy.deepcopy(cell)
                sampled.append((cell_copy, q))
    
    print(f"    Sampled {len(sampled)} cells total")
    return sampled


def sample_uniform(cells: List[Cell], n_samples: int = 30, seed: int = 42) -> List[Cell]:
    """
    Sample cells uniformly.
    
    Args:
        cells: List of Cell objects to sample from
        n_samples: Number of cells to sample
        seed: Random seed
    
    Returns:
        List[Cell]: Sampled cells
    """
    random.seed(seed)
    np.random.seed(seed)  # Also seed numpy for full reproducibility
    
    print(f"  Sampling {n_samples} cells uniformly...")
    
    if n_samples > len(cells):
        raise ValueError(f"Cannot sample {n_samples} from {len(cells)} cells")
    
    sample_indices = random.sample(range(len(cells)), n_samples)
    sampled = [copy.deepcopy(cells[idx]) for idx in sample_indices]
    
    print(f"    Sampled {len(sampled)} cells")
    return sampled


def mix_petri_with_snapshot(petri: PetriDish, snapshot_cells: List[Cell], 
                           mix_ratio: float = 0.8, seed: Optional[int] = None) -> int:
    """
    Add cells from snapshot to existing PetriDish.
    
    Args:
        petri: PetriDish to modify (in place)
        snapshot_cells: List[Cell] to sample from
        mix_ratio: Target percentage of snapshot cells (0.8 = 80%)
        seed: Random seed for sampling
    
    Returns:
        int: Total cells after mixing
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)  # Also seed numpy for full reproducibility
    
    current_count = len(petri.cells)
    
    # Calculate how many cells to add
    # If current cells are (1-mix_ratio), calculate total
    percentage_current = 1.0 - mix_ratio
    target_total = int(current_count / percentage_current)
    n_to_add = target_total - current_count
    
    print(f"    Mixing: {current_count} current + {n_to_add} from snapshot = {target_total} total")
    
    if n_to_add > len(snapshot_cells):
        # If we need more cells than available, sample with replacement
        print(f"    Warning: Need {n_to_add} cells but only {len(snapshot_cells)} available - sampling with replacement")
        sample_indices = [random.randint(0, len(snapshot_cells)-1) for _ in range(n_to_add)]
    else:
        # Sample without replacement
        sample_indices = random.sample(range(len(snapshot_cells)), n_to_add)
    added_cells = [copy.deepcopy(snapshot_cells[idx]) for idx in sample_indices]
    
    # Add to petri
    petri.cells.extend(added_cells)
    
    # Shuffle to mix thoroughly
    random.shuffle(petri.cells)
    
    return len(petri.cells)


def create_pure_snapshot_petri(snapshot_cells: List[Cell], n_cells: int = 5120, 
                              rate: float = 0.005, seed: Optional[int] = None) -> PetriDish:
    """
    Create a PetriDish with pure snapshot cells (for control2).
    
    Args:
        snapshot_cells: Cells to sample from
        n_cells: Number of cells to sample
        rate: Methylation rate for the PetriDish
        seed: Random seed
    
    Returns:
        PetriDish: New PetriDish with sampled cells
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)  # Also seed numpy for full reproducibility
    
    if n_cells > len(snapshot_cells):
        raise ValueError(f"Cannot sample {n_cells} from {len(snapshot_cells)} cells")
    
    # Sample cells
    sample_indices = random.sample(range(len(snapshot_cells)), n_cells)
    sampled_cells = [copy.deepcopy(snapshot_cells[idx]) for idx in sample_indices]
    
    # Get parameters from first cell
    if sampled_cells:
        n = len(sampled_cells[0].cpg_sites)
        gene_size = sampled_cells[0].gene_size
    else:
        n = 1000
        gene_size = GENE_SIZE
    
    # Create PetriDish
    petri = PetriDish(rate=rate, n=n, gene_size=gene_size, seed=None)
    petri.cells = sampled_cells
    petri.year = sampled_cells[0].age if sampled_cells else 60
    
    return petri


def grow_petri_for_years(petri: PetriDish, years: int, growth_phase: Optional[int] = None, verbose: bool = True) -> None:
    """
    Grow a PetriDish for specified years with optional homeostasis after growth phase.
    
    Args:
        petri: PetriDish to grow (modified in place)
        years: Total number of years to simulate
        growth_phase: Years of exponential growth before homeostasis (None = pure exponential)
        verbose: Print progress
    """
    if growth_phase is not None and growth_phase > years:
        raise ValueError(f"growth_phase ({growth_phase}) cannot exceed total years ({years})")
    
    for year in range(years):
        initial_count = len(petri.cells)
        current_year = year + 1  # 1-indexed for clarity
        
        if growth_phase is None or current_year <= growth_phase:
            # Growth phase: divide + methylate
            if not verbose:
                petri.divide_cells()
                petri.methylate_cells()
            else:
                # Temporarily redirect print output
                import io
                from contextlib import redirect_stdout
                
                f = io.StringIO()
                with redirect_stdout(f):
                    petri.divide_cells()
                    petri.methylate_cells()
                
                final_count = len(petri.cells)
                print(f"      Year {current_year} (growth): {initial_count} → {final_count} cells")
        else:
            # Homeostasis phase: divide + cull + methylate
            if not verbose:
                petri.divide_cells()
                petri.random_cull_cells()
                petri.methylate_cells()
            else:
                # Temporarily redirect print output
                import io
                from contextlib import redirect_stdout
                
                f = io.StringIO()
                with redirect_stdout(f):
                    petri.divide_cells()
                    intermediate_count = len(petri.cells)
                    petri.random_cull_cells()
                    petri.methylate_cells()
                
                final_count = len(petri.cells)
                print(f"      Year {current_year} (homeostasis): {initial_count} → {final_count} cells")
        
        petri.year += 1


def get_petri_statistics(petri: PetriDish) -> Dict[str, float]:
    """
    Calculate statistics for a PetriDish.
    
    Returns:
        dict: Statistics including mean/std/min/max JSD
    """
    if not petri.cells:
        return {
            'n_cells': 0,
            'mean_jsd': 0.0,
            'std_jsd': 0.0,
            'min_jsd': 0.0,
            'max_jsd': 0.0,
            'median_jsd': 0.0
        }
    
    jsd_values = [cell.JSD for cell in petri.cells]
    
    return {
        'n_cells': len(petri.cells),
        'mean_jsd': float(np.mean(jsd_values)),
        'std_jsd': float(np.std(jsd_values)),
        'min_jsd': float(np.min(jsd_values)),
        'max_jsd': float(np.max(jsd_values)),
        'median_jsd': float(np.median(jsd_values))
    }


def check_petri_files_state(directory: str, expected_cells: int = 1024) -> Dict:
    """
    Check state of all PetriDish files in directory.
    
    Returns:
        dict: State information (counts, corrupted, etc.)
    """
    files = sorted(glob.glob(os.path.join(directory, "*.json.gz")))
    
    result = {
        'total_files': len(files),
        'cell_counts': {},
        'corrupted': [],
        'all_expected': False,
        'all_above': False,
        'expected_count': 0,
        'above_count': 0,
        'valid_files': []
    }
    
    if not files:
        return result
    
    for filepath in files:
        fname = os.path.basename(filepath)
        try:
            # Quick load to check cell count
            with gzip.open(filepath, 'rt') as f:
                data = json.load(f)
                cell_count = len(data.get('cells', []))
                result['cell_counts'][fname] = cell_count
                result['valid_files'].append(filepath)
                
                if cell_count == expected_cells:
                    result['expected_count'] += 1
                elif cell_count > expected_cells:
                    result['above_count'] += 1
                    
        except Exception as e:
            result['corrupted'].append((fname, str(e)))
            result['cell_counts'][fname] = None
    
    # Check if all valid files meet expectations
    valid_counts = [c for c in result['cell_counts'].values() if c is not None]
    if valid_counts:
        result['all_expected'] = all(c == expected_cells for c in valid_counts)
        result['all_above'] = all(c > expected_cells for c in valid_counts)
    
    return result


def get_jsd_array(petri: PetriDish) -> np.ndarray:
    """
    Extract JSD values as numpy array.
    
    Args:
        petri: PetriDish object
    
    Returns:
        numpy array of JSD values
    """
    return np.array([cell.JSD for cell in petri.cells])