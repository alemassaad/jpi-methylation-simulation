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
from typing import List, Dict, Tuple, Optional, Any, Union, TextIO

# Add parent directory to path to import from phase1
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import Cell, PetriDish, GENE_SIZE, BASELINE_METHYLATION_DISTRIBUTION, RATE


def smart_open(filepath: str, mode: str = 'r') -> Union[TextIO, gzip.GzipFile]:
    """
    Automatically open file based on extension.
    
    Args:
        filepath: Path to file (can end with .json or .json.gz)
        mode: File mode ('r' for read, 'w' for write)
    
    Returns:
        File handle (text mode)
    
    Raises:
        ValueError: If file extension is not .json or .json.gz
    """
    if filepath.endswith('.json.gz') or filepath.endswith('.gz'):
        # Text mode for gzip (rt/wt)
        return gzip.open(filepath, mode + 't' if 't' not in mode else mode)
    elif filepath.endswith('.json'):
        # Regular text mode
        return open(filepath, mode)
    else:
        raise ValueError(f"Unsupported file extension: {filepath}. Expected .json or .json.gz")


def dict_to_cell(cell_dict: Dict[str, Any], rate: Optional[float] = None) -> Cell:
    """
    Convert dictionary from JSON to Cell object.
    
    Args:
        cell_dict: Dictionary with cell data
        rate: Optional rate override (deprecated, use cell's own rate)
    
    Returns:
        Cell: Initialized Cell object
    """
    # Determine rate configuration
    if 'gene_rate_groups' in cell_dict:
        # Gene-specific rates
        gene_rate_groups = [tuple(group) for group in cell_dict['gene_rate_groups']]
        cell = Cell(
            n=len(cell_dict['cpg_sites']),
            gene_rate_groups=gene_rate_groups,
            gene_size=cell_dict.get('gene_size', GENE_SIZE)
        )
    else:
        # Uniform rate (backward compatible)
        cell_rate = cell_dict.get('rate', rate or RATE)
        cell = Cell(
            n=len(cell_dict['cpg_sites']),
            rate=cell_rate,
            gene_size=cell_dict.get('gene_size', GENE_SIZE)
        )
    
    # Set attributes directly
    cell.cpg_sites = cell_dict['cpg_sites']
    cell.age = cell_dict['age']
    cell.methylation_proportion = cell_dict['methylation_proportion']
    cell.methylation_distribution = cell_dict['methylation_distribution']
    cell.cell_JSD = cell_dict.get('cell_jsd', cell_dict.get('jsd', 0.0))  # Support old format for compatibility
    
    # Restore site_rates if available (faster) or rebuild
    if 'site_rates' in cell_dict:
        # Try to import numpy for efficiency
        try:
            import numpy as np
            cell.site_rates = np.array(cell_dict['site_rates'], dtype=np.float64)
        except ImportError:
            cell.site_rates = cell_dict['site_rates']
    else:
        cell._build_site_rates()
    
    return cell


def load_snapshot_as_cells(simulation_file: str, year: int) -> List[Cell]:
    """
    Load a year from simulation and convert to Cell objects.
    Expects new format with 'parameters' and 'history' sections.
    
    Args:
        simulation_file: Path to phase1 simulation
        year: Year to extract (50 or 60)
    
    Returns:
        List[Cell]: List of Cell objects
    """
    print(f"Loading year {year} snapshot from {simulation_file}...")
    
    # Handle both compressed and uncompressed files
    if simulation_file.endswith('.gz'):
        with smart_open(simulation_file, 'r') as f:
            data = json.load(f)
    else:
        with open(simulation_file, 'r') as f:
            data = json.load(f)
    
    year_str = str(year)
    
    # Extract parameters and history (new format)
    params = data['parameters']
    history = data['history']
    
    if year_str not in history:
        available_years = sorted([int(y) for y in history.keys()])
        raise ValueError(f"Year {year} not found. Available: {available_years}")
    
    year_data = history[year_str]
    cell_dicts = year_data['cells']
    
    # Extract rate configuration from parameters
    rate = params.get('rate')
    gene_rate_groups = params.get('gene_rate_groups')
    gene_size = params.get('gene_size', GENE_SIZE)
    
    # Use Cell.from_dict with parameters
    cells = [Cell.from_dict(cd, rate=rate, gene_rate_groups=gene_rate_groups, 
                           gene_size=gene_size) for cd in cell_dicts]
    
    print(f"  Loaded {len(cells)} cells from year {year}")
    return cells


def save_snapshot_cells(cells: List[Cell], filepath: str, compress: bool = True) -> None:
    """
    Save a list of Cell objects as a snapshot.
    
    Args:
        cells: List of Cell objects
        filepath: Output path for JSON file
        compress: If True, save as .json.gz; if False, save as .json
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # Adjust filepath based on compress flag
    if compress and not filepath.endswith('.gz'):
        if filepath.endswith('.json'):
            filepath = filepath + '.gz'
        else:
            filepath = filepath + '.json.gz'
    elif not compress and filepath.endswith('.gz'):
        filepath = filepath[:-3]  # Remove .gz extension
    
    # Extract rate configuration from first cell
    metadata = {
        "num_cells": len(cells),
        "year": cells[0].age if cells else 0
    }
    
    if cells:
        first_cell = cells[0]
        if hasattr(first_cell, 'rate') and first_cell.rate is not None:
            metadata['rate'] = first_cell.rate
        if hasattr(first_cell, 'gene_rate_groups') and first_cell.gene_rate_groups:
            metadata['gene_rate_groups'] = first_cell.gene_rate_groups
        metadata['gene_size'] = getattr(first_cell, 'gene_size', GENE_SIZE)
    
    snapshot_data = {
        "metadata": metadata,
        "cells": [cell.to_dict() for cell in cells]
    }
    
    # Use smart_open which handles compression based on file extension
    with smart_open(filepath, 'w') as f:
        json.dump(snapshot_data, f, indent=2)
    
    print(f"  Saved {len(cells)} cells to {filepath}")


def load_snapshot_cells(filepath: str) -> List[Cell]:
    """
    Load a snapshot file and convert to Cell objects.
    Auto-detects whether file is compressed or not.
    
    Args:
        filepath: Path to snapshot file (.json or .json.gz)
    
    Returns:
        List[Cell]: List of Cell objects
    """
    # Use smart_open for automatic format detection
    with smart_open(filepath, 'r') as f:
        data = json.load(f)
    
    cell_dicts = data['cells']
    
    # Detect format based on first cell
    if cell_dicts and 'methylated' in cell_dicts[0] and 'cpg_sites' not in cell_dicts[0]:
        # New lean format - use Cell.from_dict
        metadata = data.get('metadata', {})
        rate = metadata.get('rate')
        gene_rate_groups = metadata.get('gene_rate_groups')
        gene_size = metadata.get('gene_size', GENE_SIZE)
        cells = [Cell.from_dict(cd, rate=rate, gene_rate_groups=gene_rate_groups, 
                               gene_size=gene_size) for cd in cell_dicts]
    else:
        # Old format - use dict_to_cell
        cells = [dict_to_cell(cd) for cd in cell_dicts]
    
    return cells


def save_petri_dish(petri: PetriDish, filepath: str, metadata: Optional[Dict] = None, 
                   include_cell_history: bool = False, include_gene_jsd: bool = False,
                   include_gene_metrics: bool = False,
                   compress: bool = True) -> None:
    """
    Save a PetriDish to JSON file (compressed or uncompressed), optionally with history.
    
    Args:
        petri: PetriDish object
        filepath: Output path
        metadata: Extra metadata dict
        include_cell_history: Whether to save cell history
        include_gene_jsd: Whether to save gene JSD history
        include_gene_metrics: Whether to calculate and save gene-level JSD and mean methylation
        compress: If True, save as .json.gz; if False, save as .json
    """
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    
    # Calculate cell JSD statistics
    cell_jsd_values = [cell.cell_JSD for cell in petri.cells]
    
    data = {
        "metadata": {
            "num_cells": len(petri.cells),
            "mean_cell_jsd": float(np.mean(cell_jsd_values)) if cell_jsd_values else 0.0,
            "std_cell_jsd": float(np.std(cell_jsd_values)) if cell_jsd_values else 0.0,
            "min_cell_jsd": float(np.min(cell_jsd_values)) if cell_jsd_values else 0.0,
            "max_cell_jsd": float(np.max(cell_jsd_values)) if cell_jsd_values else 0.0,
            "median_cell_jsd": float(np.median(cell_jsd_values)) if cell_jsd_values else 0.0,
            "year": petri.year,
            "growth_phase": petri.growth_phase,
            "rate": petri.rate,
            "has_cell_history": include_cell_history and hasattr(petri, 'cell_history') and bool(petri.cell_history),
            "has_gene_jsd": include_gene_jsd and hasattr(petri, 'gene_jsd_history') and bool(petri.gene_jsd_history),
            "gene_size": getattr(petri, 'gene_size', GENE_SIZE),
            "n_genes": getattr(petri, 'n_genes', petri.n // GENE_SIZE)
        },
        "cells": [cell.to_dict() for cell in petri.cells]
    }
    
    # Add cell history if requested and available
    if include_cell_history and hasattr(petri, 'cell_history') and petri.cell_history:
        data["cell_history"] = petri.cell_history
    
    # Add gene JSD history if requested and available
    if include_gene_jsd and hasattr(petri, 'gene_jsd_history') and petri.gene_jsd_history:
        data["gene_jsd_history"] = petri.gene_jsd_history
    
    # Add custom metadata if provided
    if metadata:
        data["metadata"].update(metadata)
    
    # Add petri metadata if it exists
    if hasattr(petri, 'metadata'):
        data["metadata"].update(petri.metadata)
    
    # Add gene-level metrics if requested
    if include_gene_metrics and len(petri.cells) > 1:
        try:
            # Calculate gene JSDs and mean methylation
            gene_jsds = petri.calculate_gene_jsds()
            gene_means = petri.calculate_gene_mean_methylation()
            
            # Add to metadata
            data["metadata"]["gene_jsds"] = [float(jsd) for jsd in gene_jsds]
            data["metadata"]["gene_mean_methylation"] = [float(mean) for mean in gene_means]
            data["metadata"]["n_genes"] = len(gene_jsds)
        except Exception as e:
            print(f"  Warning: Could not calculate gene metrics: {e}")
    
    # Adjust filepath based on compress flag
    if compress and not filepath.endswith('.gz'):
        if filepath.endswith('.json'):
            filepath = filepath + '.gz'
        else:
            filepath = filepath + '.json.gz'
    elif not compress and filepath.endswith('.gz'):
        filepath = filepath[:-3]  # Remove .gz extension
    
    # Use smart_open which handles compression based on file extension
    with smart_open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


def load_petri_dish(filepath: str, include_cell_history: bool = False, include_gene_jsd: bool = False) -> PetriDish:
    """
    Load a PetriDish from file, optionally with history.
    Auto-detects whether file is compressed or not.
    
    Args:
        filepath: Path to .json or .json.gz file
        include_cell_history: Whether to load cell history if available
        include_gene_jsd: Whether to load gene JSD history if available
        
    Returns:
        PetriDish: Reconstructed object with cells, metadata, and optionally history
    """
    # Use smart_open for automatic format detection
    with smart_open(filepath, 'r') as f:
        data = json.load(f)
    
    # Get first cell to extract parameters and determine rate configuration
    if data['cells']:
        first_cell_dict = data['cells'][0]
        # Handle both old format (cpg_sites) and new format (methylated)
        if 'cpg_sites' in first_cell_dict:
            n = len(first_cell_dict['cpg_sites'])
        elif 'methylated' in first_cell_dict:
            n = len(first_cell_dict['methylated'])
        else:
            raise ValueError("Cell dict has neither 'cpg_sites' nor 'methylated' field")
        
        # Get gene_size from metadata or first cell
        gene_size = data.get('metadata', {}).get('gene_size', first_cell_dict.get('gene_size', GENE_SIZE))
        
        # Determine rate configuration from first cell
        if 'gene_rate_groups' in first_cell_dict:
            gene_rate_groups = [tuple(group) for group in first_cell_dict['gene_rate_groups']]
            petri = PetriDish(
                gene_rate_groups=gene_rate_groups,
                n=n,
                gene_size=data.get('metadata', {}).get('gene_size', gene_size),
                seed=None
            )
        else:
            rate = first_cell_dict.get('rate', data.get('metadata', {}).get('rate', 0.005))
            petri = PetriDish(
                rate=rate,
                n=n,
                gene_size=data.get('metadata', {}).get('gene_size', gene_size),
                seed=None
            )
    else:
        # Defaults if no cells
        petri = PetriDish(
            rate=0.005,
            n=1000,
            gene_size=GENE_SIZE,
            seed=None
        )
    
    # Replace cells with loaded cells
    # Detect format based on first cell
    if data['cells'] and 'methylated' in data['cells'][0] and 'cpg_sites' not in data['cells'][0]:
        # New lean format - use Cell.from_dict
        rate = data.get('metadata', {}).get('rate')
        gene_rate_groups = data.get('metadata', {}).get('gene_rate_groups')
        gene_size = data.get('metadata', {}).get('gene_size', GENE_SIZE)
        petri.cells = [Cell.from_dict(cd, rate=rate, gene_rate_groups=gene_rate_groups, 
                                     gene_size=gene_size) for cd in data['cells']]
    else:
        # Old format or phase2 format - use dict_to_cell
        petri.cells = [dict_to_cell(cd) for cd in data['cells']]
    
    # Validate consistency
    petri._validate_cell_rate_consistency()
    
    # Restore year if available
    if 'year' in data['metadata']:
        petri.year = data['metadata']['year']
    
    # Add metadata attribute if not exists
    if not hasattr(petri, 'metadata'):
        petri.metadata = {}
    
    # Restore metadata
    petri.metadata = data.get('metadata', {})
    
    # Restore cell history if requested and available
    if include_cell_history and 'cell_history' in data:
        petri.cell_history = data['cell_history']
        petri.track_cell_history = True
    elif include_cell_history and 'history' in data:  # Backward compatibility
        petri.cell_history = data['history']
        petri.track_cell_history = True
    else:
        # Clear any auto-initialized history if not loading it
        petri.cell_history = {}
        petri.track_cell_history = False
    
    # Restore gene JSD history if requested and available
    if include_gene_jsd and 'gene_jsd_history' in data:
        petri.gene_jsd_history = data['gene_jsd_history']
        petri.track_gene_jsd = True
    else:
        petri.gene_jsd_history = {}
        petri.track_gene_jsd = False
        petri.history_enabled = False
    
    return petri


def load_all_petri_dishes(directory: str, include_cell_history: bool = False) -> List[PetriDish]:
    """
    Load all PetriDish objects from a directory.
    
    Args:
        directory: Directory containing .json or .json.gz files
        include_cell_history: Whether to load cell history if available
    
    Returns:
        List of PetriDish objects
    """
    # Look for both compressed and uncompressed files
    gz_files = glob.glob(os.path.join(directory, "*.json.gz"))
    json_files = glob.glob(os.path.join(directory, "*.json"))
    files = sorted(gz_files + json_files)
    dishes = []
    
    for filepath in files:
        try:
            petri = load_petri_dish(filepath, include_cell_history=include_cell_history)
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
    sorted_cells = sorted(cells, key=lambda c: c.cell_JSD)
    
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


def create_control2_with_uniform_base(
    snapshot_cells: List[Cell],
    uniform_pool: List[Cell],
    uniform_indices: List[int],
    target_size: int,
    rate: float = 0.005,
    seed: Optional[int] = None
) -> PetriDish:
    """
    Create Control2 individual with uniform base + additional snapshot cells.
    
    When using --uniform-mixing, Control2 shares the same 80% base cells as
    mutant and control1, plus additional random snapshot cells for the remaining 20%.
    
    Args:
        snapshot_cells: Full second snapshot cells to sample additional from
        uniform_pool: The shared cells used in all individuals (e.g., 80%)
        uniform_indices: Indices of cells already used in uniform pool
        target_size: Total number of cells needed
        rate: Methylation rate for PetriDish
        seed: Random seed for additional sampling
    
    Returns:
        PetriDish with uniform base + additional unique snapshot cells
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Start with the uniform pool (deep copy to avoid reference issues)
    combined_cells = [copy.deepcopy(cell) for cell in uniform_pool]
    
    # Calculate how many additional cells needed
    n_additional = target_size - len(uniform_pool)
    
    if n_additional > 0:
        # Create set of available indices (excluding those already used)
        all_indices = set(range(len(snapshot_cells)))
        used_indices = set(uniform_indices)
        available_indices = list(all_indices - used_indices)
        
        if len(available_indices) >= n_additional:
            # Sample from available indices
            additional_indices = random.sample(available_indices, n_additional)
        else:
            # Not enough unique cells, need to sample with replacement
            print(f"      ⚠️ Not enough unique cells for Control2 additional sampling")
            print(f"         Need {n_additional} but only {len(available_indices)} available")
            print(f"         Sampling with replacement from available cells")
            additional_indices = [random.choice(available_indices) for _ in range(n_additional)]
        
        # Add the additional cells
        additional_cells = [copy.deepcopy(snapshot_cells[idx]) for idx in additional_indices]
        combined_cells.extend(additional_cells)
    elif n_additional < 0:
        # Edge case: uniform pool is larger than target (e.g., 100% mix ratio)
        # Just use subset of uniform pool
        combined_cells = combined_cells[:target_size]
    
    # Get parameters from first cell
    if combined_cells:
        n = len(combined_cells[0].cpg_sites)
        gene_size = combined_cells[0].gene_size
    else:
        n = 1000
        gene_size = GENE_SIZE
    
    # Create PetriDish
    petri = PetriDish(rate=rate, n=n, gene_size=gene_size, seed=None)
    petri.cells = combined_cells
    petri.year = combined_cells[0].age if combined_cells else 60
    
    return petri


def grow_petri_for_years(petri: PetriDish, years: int, growth_phase: Optional[int] = None, 
                        verbose: bool = True, track_history: bool = False, 
                        start_year: Optional[int] = None, track_gene_jsd: bool = False) -> None:
    """
    Grow a PetriDish for specified years with optional homeostasis after growth phase.
    Now uses PetriDish's built-in grow_with_homeostasis method when history tracking is enabled.
    
    Args:
        petri: PetriDish to grow (modified in place)
        years: Total number of years to simulate
        growth_phase: Years of exponential growth before homeostasis (None = pure exponential)
        verbose: Print progress
        track_history: Enable cell history tracking
        start_year: Starting year for history tracking
        track_gene_jsd: Enable gene JSD tracking
    """
    if growth_phase is not None and growth_phase > years:
        raise ValueError(f"growth_phase ({growth_phase}) cannot exceed total years ({years})")
    
    # If history tracking is enabled, use PetriDish's new method
    if track_history and start_year is not None:
        petri.enable_history_tracking(start_year, track_gene_jsd=track_gene_jsd)
        petri.grow_with_homeostasis(years, growth_phase, verbose, record_history=True)
        return
    
    # Otherwise, use the original implementation (for backward compatibility)
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
        dict: Statistics including mean/std/min/max cell JSD
    """
    if not petri.cells:
        return {
            'n_cells': 0,
            'mean_cell_jsd': 0.0,
            'std_cell_jsd': 0.0,
            'min_cell_jsd': 0.0,
            'max_cell_jsd': 0.0,
            'median_cell_jsd': 0.0
        }
    
    cell_jsd_values = [cell.cell_JSD for cell in petri.cells]
    
    return {
        'n_cells': len(petri.cells),
        'mean_cell_jsd': float(np.mean(cell_jsd_values)),
        'std_cell_jsd': float(np.std(cell_jsd_values)),
        'min_cell_jsd': float(np.min(cell_jsd_values)),
        'max_cell_jsd': float(np.max(cell_jsd_values)),
        'median_cell_jsd': float(np.median(cell_jsd_values))
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
            with smart_open(filepath, 'r') as f:
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


def get_cell_jsd_array(petri: PetriDish) -> np.ndarray:
    """
    Extract cell JSD values as numpy array.
    
    Args:
        petri: PetriDish object
    
    Returns:
        numpy array of cell JSD values
    """
    return np.array([cell.cell_JSD for cell in petri.cells])


def calculate_population_statistics(dishes: List[PetriDish], group_name: str) -> Dict:
    """
    Calculate size statistics for a group of PetriDish objects.
    
    Returns:
        dict with mean, std, median, min, max, cv, iqr, outliers
    """
    sizes = [len(dish.cells) for dish in dishes]
    
    if not sizes:
        return {}
    
    q1, median, q3 = np.percentile(sizes, [25, 50, 75])
    iqr = q3 - q1
    
    # Identify outliers (> 1.5 * IQR from quartiles)
    outliers = []
    for i, size in enumerate(sizes):
        if size < (q1 - 1.5 * iqr) or size > (q3 + 1.5 * iqr):
            outliers.append({'index': i, 'size': size})
    
    stats = {
        'group': group_name,
        'n': len(sizes),
        'sizes': sizes,
        'mean': float(np.mean(sizes)),
        'std': float(np.std(sizes)),
        'median': float(median),
        'min': int(min(sizes)),
        'max': int(max(sizes)),
        'q1': float(q1),
        'q3': float(q3),
        'iqr': float(iqr),
        'cv': float(np.std(sizes) / np.mean(sizes)) if np.mean(sizes) > 0 else 0,
        'outliers': outliers
    }
    
    return stats


def print_mixing_statistics(mutant_stats: Dict, control1_stats: Dict, combined_median: int):
    """
    Print formatted statistics before mixing.
    """
    print("\n  Population Size Distribution (before mixing):")
    print("  " + "="*50)
    
    # Mutant stats
    print(f"  Mutant individuals (n={mutant_stats['n']}):")
    print(f"    Mean: {mutant_stats['mean']:.1f} ± {mutant_stats['std']:.1f} cells")
    print(f"    Median: {mutant_stats['median']:.0f} cells")
    print(f"    Range: {mutant_stats['min']}-{mutant_stats['max']} cells")
    print(f"    CV: {mutant_stats['cv']*100:.1f}%")
    if mutant_stats['outliers']:
        print(f"    ⚠️ Outliers: {len(mutant_stats['outliers'])} individuals")
    
    # Control1 stats
    print(f"\n  Control1 individuals (n={control1_stats['n']}):")
    print(f"    Mean: {control1_stats['mean']:.1f} ± {control1_stats['std']:.1f} cells")
    print(f"    Median: {control1_stats['median']:.0f} cells")
    print(f"    Range: {control1_stats['min']}-{control1_stats['max']} cells")
    print(f"    CV: {control1_stats['cv']*100:.1f}%")
    if control1_stats['outliers']:
        print(f"    ⚠️ Outliers: {len(control1_stats['outliers'])} individuals")
    
    # Combined stats
    print(f"\n  Combined (n={mutant_stats['n'] + control1_stats['n']}):")
    print(f"    Median for mixing: {combined_median} cells")
    
    # Warnings
    max_cv = max(mutant_stats['cv'], control1_stats['cv'])
    if max_cv > 0.20:
        print(f"\n  ⚠️ WARNING: High variation in individual sizes (max CV={max_cv*100:.1f}%)")
        print(f"     Uniform mixing may not be appropriate")
        print(f"     Consider using default independent mixing")
    
    print("  " + "="*50)


def normalize_individuals_for_uniform_mixing(
    mutant_dishes: List[PetriDish],
    control1_dishes: List[PetriDish], 
    seed: Optional[int] = None
) -> Tuple[List[PetriDish], List[PetriDish], int]:
    """
    Normalize all individuals to the same size by sampling down to minimum.
    This ensures fair uniform mixing by eliminating size variation between individuals.
    
    Args:
        mutant_dishes: List of mutant PetriDish objects
        control1_dishes: List of control1 PetriDish objects  
        seed: Random seed for reproducibility
        
    Returns:
        Tuple of (normalized_mutant_dishes, normalized_control1_dishes, normalized_size)
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Collect all sizes
    all_dishes = mutant_dishes + control1_dishes
    all_sizes = [len(dish.cells) for dish in all_dishes]
    
    if not all_sizes:
        print("  WARNING: No individuals to normalize for uniform mixing!")
        return [], [], 0
    
    # Find minimum size to normalize to
    min_size = min(all_sizes)
    max_size = max(all_sizes)
    
    print(f"\n  === NORMALIZING INDIVIDUALS FOR UNIFORM MIXING ===")
    print(f"  Individual sizes: min={min_size}, max={max_size}")
    print(f"  Normalizing all individuals to {min_size} cells")
    
    def normalize_dish(dish: PetriDish) -> PetriDish:
        """Normalize a single dish by random sampling."""
        current_size = len(dish.cells)
        if current_size == min_size:
            # Already correct size, just deep copy
            normalized_cells = [copy.deepcopy(cell) for cell in dish.cells]
        else:
            # Randomly sample down to min_size
            sampled_indices = random.sample(range(current_size), min_size)
            normalized_cells = [copy.deepcopy(dish.cells[idx]) for idx in sampled_indices]
        
        # Create new PetriDish with normalized cells
        normalized_dish = PetriDish(
            rate=dish.methylation_rate if hasattr(dish, 'methylation_rate') else 0.005,
            n=len(normalized_cells[0].cpg_sites) if normalized_cells else 1000,
            gene_size=normalized_cells[0].gene_size if normalized_cells else 5,
            seed=None
        )
        normalized_dish.cells = normalized_cells
        normalized_dish.year = dish.year
        
        # Preserve history if it exists
        if hasattr(dish, 'cell_history'):
            normalized_dish.cell_history = dish.cell_history
        
        return normalized_dish
    
    # Normalize all dishes
    normalized_mutant = [normalize_dish(dish) for dish in mutant_dishes]
    normalized_control1 = [normalize_dish(dish) for dish in control1_dishes]
    
    # Verify normalization
    new_sizes = [len(dish.cells) for dish in normalized_mutant + normalized_control1]
    if any(size != min_size for size in new_sizes):
        print(f"  ⚠️ ERROR: Normalization failed! Sizes: {set(new_sizes)}")
    else:
        print(f"  ✓ Successfully normalized {len(new_sizes)} individuals to {min_size} cells each")
    
    return normalized_mutant, normalized_control1, min_size


def create_uniform_mixing_pool(snapshot_cells: List[Cell], 
                               normalized_size: int, 
                               mix_ratio: float,
                               seed: Optional[int] = None) -> Tuple[List[Cell], List[int]]:
    """
    Create a shared pool of snapshot cells for uniform mixing.
    
    Args:
        snapshot_cells: Available cells to sample from
        normalized_size: Normalized size of all individuals (after normalization)
        mix_ratio: Percentage from snapshot (0-1)
        seed: Random seed
        
    Returns:
        Tuple of (cells to be shared across all individuals, indices used from snapshot)
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Calculate how many cells needed
    # Handle edge case: 100% mix ratio means all cells from snapshot
    if mix_ratio >= 1.0:
        target_total = normalized_size
        n_snapshot_cells = normalized_size
    else:
        target_total = int(normalized_size / (1 - mix_ratio))
        n_snapshot_cells = int(target_total * mix_ratio)
    
    print(f"    Creating uniform mixing pool:")
    print(f"      Normalized individual size: {normalized_size} cells")
    print(f"      Target total size: {target_total} cells")
    print(f"      Snapshot cells needed: {n_snapshot_cells} cells")
    
    # Sample the cells
    if n_snapshot_cells > len(snapshot_cells):
        print(f"      ⚠️ Sampling with replacement ({n_snapshot_cells} > {len(snapshot_cells)})")
        indices = [random.randint(0, len(snapshot_cells)-1) for _ in range(n_snapshot_cells)]
    else:
        indices = random.sample(range(len(snapshot_cells)), n_snapshot_cells)
    
    # Create deep copies
    pool = [copy.deepcopy(snapshot_cells[idx]) for idx in indices]
    
    print(f"      Created pool of {len(pool)} cells")
    
    return pool, indices


def mix_petri_uniform(petri: PetriDish, 
                      uniform_pool: List[Cell],
                      target_ratio: float) -> int:
    """
    Mix PetriDish with cells from uniform pool.
    Now simplified: all individuals are already normalized, just add the entire pool.
    
    Args:
        petri: PetriDish to mix into (already normalized)
        uniform_pool: Pre-sampled shared cells (sized exactly for the normalized individuals)
        target_ratio: What fraction should be from pool (used for 100% edge case only)
        
    Returns:
        Total cells after mixing
    """
    # Handle edge case: 100% mix ratio means replace all cells
    if target_ratio >= 1.0:
        # Replace all cells with snapshot cells
        current_size = len(petri.cells)
        petri.cells = [copy.deepcopy(cell) for cell in uniform_pool[:current_size]]
        random.shuffle(petri.cells)
        return len(petri.cells)
    
    # Normal case: individuals are already normalized, just add the entire pool
    added_cells = [copy.deepcopy(cell) for cell in uniform_pool]
    petri.cells.extend(added_cells)
    
    # Shuffle to mix
    random.shuffle(petri.cells)
    
    return len(petri.cells)


def normalize_populations(mutant_dishes: List[PetriDish], 
                         control1_dishes: List[PetriDish],
                         seed: Optional[int] = None) -> Tuple[List[PetriDish], List[PetriDish], int]:
    """
    Normalize populations to the same size using median - 0.5 stdev threshold.
    
    This threshold was chosen to:
    - Exclude individuals that are somewhat below average (likely due to stochastic effects)
    - Preserve most of the population (less aggressive than median - 1σ)
    - Adapt to actual variance in the population (unlike fixed percentiles)
    
    Algorithm:
    1. Calculate median and standard deviation of all individual sizes
    2. Set threshold = median - 0.5 * standard_deviation
    3. Exclude individuals below this threshold
    4. Randomly sample cells from larger individuals to match threshold size
    
    Args:
        mutant_dishes: List of mutant PetriDish objects
        control1_dishes: List of control1 PetriDish objects
        seed: Random seed for reproducibility
        
    Returns:
        Tuple of (normalized_mutant_dishes, normalized_control1_dishes, threshold_size)
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Collect all sizes
    all_dishes = mutant_dishes + control1_dishes
    all_sizes = [len(dish.cells) for dish in all_dishes]
    
    # Handle edge case: no individuals
    if not all_sizes:
        print("  WARNING: No individuals to normalize!")
        return [], [], 0
    
    # Handle edge case: only one individual
    if len(all_sizes) == 1:
        print(f"  WARNING: Only 1 individual found with {all_sizes[0]} cells")
        return mutant_dishes, control1_dishes, all_sizes[0]
    
    # Calculate threshold using median - 0.5 * stdev
    median_size = np.median(all_sizes)
    std_size = np.std(all_sizes)
    threshold_size = int(median_size - 0.5 * std_size)
    
    # Ensure threshold is positive and reasonable
    if threshold_size < 1:
        print(f"  WARNING: Calculated threshold ({threshold_size}) < 1, using minimum size instead")
        threshold_size = min(all_sizes)
    
    print(f"\n  === NORMALIZATION STATISTICS (median - 0.5σ method) ===")
    print(f"  All individual sizes: min={min(all_sizes)}, max={max(all_sizes)}")
    print(f"  Median: {median_size:.1f} cells")
    print(f"  Std Dev: {std_size:.1f} cells")
    print(f"  Threshold calculation: {median_size:.1f} - 0.5 * {std_size:.1f} = {threshold_size} cells")
    
    # Process mutant dishes
    normalized_mutant = []
    mutant_excluded = 0
    mutant_trimmed = 0
    mutant_kept_as_is = 0
    
    for i, dish in enumerate(mutant_dishes):
        current_size = len(dish.cells)
        
        if current_size < threshold_size:
            # Exclude this individual
            mutant_excluded += 1
            print(f"    Mutant {i:02d}: {current_size} cells - EXCLUDED (below threshold)")
        elif current_size == threshold_size:
            # Keep as is
            mutant_kept_as_is += 1
            normalized_mutant.append(dish)
            print(f"    Mutant {i:02d}: {current_size} cells - kept as is")
        else:
            # Trim to threshold size
            mutant_trimmed += 1
            # Create a deep copy to avoid modifying original
            new_dish = copy.deepcopy(dish)
            # Randomly sample cells to match threshold
            new_dish.cells = random.sample(new_dish.cells, threshold_size)
            normalized_mutant.append(new_dish)
            print(f"    Mutant {i:02d}: {current_size} → {threshold_size} cells (trimmed)")
    
    # Process control1 dishes
    normalized_control1 = []
    control1_excluded = 0
    control1_trimmed = 0
    control1_kept_as_is = 0
    
    for i, dish in enumerate(control1_dishes):
        current_size = len(dish.cells)
        
        if current_size < threshold_size:
            # Exclude this individual
            control1_excluded += 1
            print(f"    Control1 {i:02d}: {current_size} cells - EXCLUDED (below threshold)")
        elif current_size == threshold_size:
            # Keep as is
            control1_kept_as_is += 1
            normalized_control1.append(dish)
            print(f"    Control1 {i:02d}: {current_size} cells - kept as is")
        else:
            # Trim to threshold size
            control1_trimmed += 1
            # Create a deep copy to avoid modifying original
            new_dish = copy.deepcopy(dish)
            # Randomly sample cells to match threshold
            new_dish.cells = random.sample(new_dish.cells, threshold_size)
            normalized_control1.append(new_dish)
            print(f"    Control1 {i:02d}: {current_size} → {threshold_size} cells (trimmed)")
    
    # Print summary
    print(f"\n  Mutant summary:")
    print(f"    Kept: {len(normalized_mutant)} ({mutant_kept_as_is} as-is, {mutant_trimmed} trimmed)")
    print(f"    Excluded: {mutant_excluded}")
    print(f"  Control1 summary:")
    print(f"    Kept: {len(normalized_control1)} ({control1_kept_as_is} as-is, {control1_trimmed} trimmed)")
    print(f"    Excluded: {control1_excluded}")
    
    total_kept = len(normalized_mutant) + len(normalized_control1)
    total_original = len(mutant_dishes) + len(control1_dishes)
    retention_rate = (total_kept / total_original) * 100 if total_original > 0 else 0
    print(f"\n  Overall retention: {total_kept}/{total_original} ({retention_rate:.1f}%)")
    
    # Handle edge case: all individuals excluded
    if total_kept == 0:
        print("\n  ERROR: All individuals were excluded by normalization!")
        print(f"  Threshold was {threshold_size} but all individuals were smaller")
        print("  Consider adjusting growth parameters or using a different threshold")
        # Return empty lists but with a valid threshold for metadata
        return [], [], threshold_size
    
    # Handle edge case: only one group has individuals
    if len(normalized_mutant) == 0:
        print("\n  WARNING: All mutant individuals were excluded!")
    if len(normalized_control1) == 0:
        print("\n  WARNING: All control1 individuals were excluded!")
    
    return normalized_mutant, normalized_control1, threshold_size