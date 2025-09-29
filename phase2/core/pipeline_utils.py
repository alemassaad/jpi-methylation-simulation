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
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any, Union, TextIO

# Add parent directory to path to import from phase1
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import Cell, PetriDish, GENE_SIZE, BASELINE_METHYLATION_DISTRIBUTION, RATE, rate_to_gene_rate_groups


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


def dict_to_cell(cell_dict: Dict[str, Any]) -> Cell:
    """
    Convert dictionary from JSON to Cell object.
    Now expects gene_rate_groups in each cell (new format from Phase 1).
    
    Args:
        cell_dict: Dictionary with cell data
    
    Returns:
        Cell: Initialized Cell object
    """
    # Extract gene_rate_groups from the cell itself (required)
    gene_rate_groups = cell_dict.get('gene_rate_groups')
    if not gene_rate_groups:
        raise ValueError(
            "Cell data missing gene_rate_groups! "
            "This cell appears to be from an old simulation. "
            "Please re-run Phase 1 with the latest code."
        )
    
    # Convert to tuples
    gene_rate_groups = [tuple(group) for group in gene_rate_groups]
    
    # Get cpg_sites from cell dictionary
    sites = cell_dict['cpg_sites']
    
    # Create cell with gene_rate_groups only
    cell = Cell(
        n=len(sites),
        gene_rate_groups=gene_rate_groups,
        gene_size=cell_dict.get('gene_size', GENE_SIZE)
    )
    
    # Set the sites
    cell.cpg_sites = sites
    
    # Set all enhanced attributes from Phase 1
    cell.age = cell_dict.get('age', 0)
    # Handle both old and new naming for backward compatibility
    cell.cell_jsd = cell_dict.get('cell_jsd', cell_dict.get('cell_JSD', 0.0))
    
    # Use stored methylation stats if available (more efficient)
    if 'n_methylated' in cell_dict:
        cell.methylation_proportion = cell_dict.get('methylation_proportion', 0.0)
    else:
        # Recalculate if missing
        n_methylated = sum(1 for site in sites if site == 1)
        cell.methylation_proportion = n_methylated / len(sites) if sites else 0.0
    
    # Restore site_rates if available (faster) or rebuild
    if 'site_rates' in cell_dict:
        try:
            import numpy as np
            cell.site_rates = np.array(cell_dict['site_rates'], dtype=np.float64)
        except ImportError:
            cell.site_rates = cell_dict['site_rates']
    else:
        cell._build_site_rates()
    
    return cell


def load_snapshot_as_cells(simulation_file: str, year: int, 
                          expected_gene_rate_groups: List[Tuple[int, float]] = None) -> List[Cell]:
    """
    Load a year from simulation and convert to Cell objects.
    Now validates gene_rate_groups consistency.
    
    Args:
        simulation_file: Path to phase1 simulation
        year: Year to extract (50 or 60)
        expected_gene_rate_groups: Expected gene rate groups for validation
    
    Returns:
        List[Cell]: List of Cell objects
    """
    print(f"Loading year {year} snapshot from {simulation_file}...")
    
    # Handle both compressed and uncompressed files
    with smart_open(simulation_file, 'r') as f:
        data = json.load(f)
    
    year_str = str(year)
    
    # Extract config/parameters and history (support both old and new format)
    params = data.get('config', data.get('parameters', {}))
    history = data['history']
    
    if year_str not in history:
        available_years = sorted([int(y) for y in history.keys()])
        raise ValueError(f"Year {year} not found. Available: {available_years}")
    
    # Get gene_rate_groups from parameters
    gene_rate_groups = params.get('gene_rate_groups')
    if not gene_rate_groups:
        # Handle old format files (convert rate to gene_rate_groups)
        rate = params.get('rate')
        if rate:
            n = params.get('n', 1000)
            gene_size = params.get('gene_size', GENE_SIZE)
            gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
        else:
            raise ValueError("Simulation has neither gene_rate_groups nor rate!")
    
    # Convert to tuples
    gene_rate_groups = [tuple(group) for group in gene_rate_groups]
    
    # Validate against expected if provided
    if expected_gene_rate_groups and gene_rate_groups != expected_gene_rate_groups:
        raise ValueError(
            f"Gene rate groups mismatch in snapshot!\n"
            f"  Expected: {expected_gene_rate_groups}\n"
            f"  Found in simulation: {gene_rate_groups}"
        )
    
    # Load cells
    year_data = history[year_str]
    cell_dicts = year_data['cells']
    cells = []
    
    for i, cd in enumerate(cell_dicts):
        # Check if new format (has gene_rate_groups in cell)
        if 'gene_rate_groups' in cd:
            cell = dict_to_cell(cd)
        else:
            # Old format - use Cell.from_dict with parameters
            gene_size = params.get('gene_size', GENE_SIZE)
            cell = Cell.from_dict(cd, gene_rate_groups=gene_rate_groups, gene_size=gene_size)
        
        # Validate each cell's gene_rate_groups
        if cell.gene_rate_groups != gene_rate_groups:
            raise ValueError(
                f"Cell {i} from year {year} has inconsistent gene_rate_groups!\n"
                f"  Expected: {gene_rate_groups}\n"
                f"  Cell has: {cell.gene_rate_groups}"
            )
        
        cells.append(cell)
    
    print(f"  ✓ Loaded {len(cells)} cells with gene_rate_groups: {gene_rate_groups}")
    return cells


# Note: save_snapshot_cells() has been removed - snapshots are now direct copies
# from phase1 simulation years, saved by extract_snapshots.py


def load_snapshot_cells(filepath: str) -> List[Cell]:
    """
    Load a snapshot file and convert to Cell objects.
    
    Handles snapshot format with year key: {"30": {"cells": [...], "gene_jsd": [...]}}
    The year key wrapper is stripped and the cells are extracted.
    
    Args:
        filepath: Path to snapshot file (.json or .json.gz)
    
    Returns:
        List[Cell]: List of Cell objects
    """
    # Use smart_open for automatic format detection
    with smart_open(filepath, 'r') as f:
        data = json.load(f)
    
    # Handle year key wrapper (e.g., {"30": {...}})
    if len(data) == 1 and isinstance(list(data.keys())[0], str) and list(data.keys())[0].isdigit():
        # Has year key wrapper, extract the year data
        year_key = list(data.keys())[0]
        year_data = data[year_key]
    else:
        # Assume it's direct year data (backward compatibility just in case)
        year_data = data
    
    # Direct access to cells
    cell_dicts = year_data['cells']
    cells = []
    
    # Need to get gene_rate_groups from metadata.json if cells don't have it
    gene_rate_groups_fallback = None
    
    # Try to load metadata.json for fallback gene_rate_groups
    snapshot_dir = os.path.dirname(filepath)
    metadata_file = os.path.join(snapshot_dir, 'metadata.json')
    if os.path.exists(metadata_file):
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
            gene_rate_groups_fallback = metadata.get('gene_rate_groups')
            if gene_rate_groups_fallback:
                gene_rate_groups_fallback = [tuple(g) for g in gene_rate_groups_fallback]
            gene_size_fallback = metadata.get('gene_size', GENE_SIZE)
    
    for cd in cell_dicts:
        if 'gene_rate_groups' in cd:
            # Cell has gene_rate_groups embedded (current phase1 format)
            cells.append(dict_to_cell(cd))
        else:
            # Old format - use fallback from metadata.json
            if not gene_rate_groups_fallback:
                raise ValueError(
                    "Cells don't have gene_rate_groups and no metadata.json found! "
                    "Please re-run phase1 with latest code."
                )
            cell = Cell.from_dict(cd, gene_rate_groups=gene_rate_groups_fallback, gene_size=gene_size_fallback)
            cells.append(cell)
    
    # Validate all cells have same gene_rate_groups
    if cells:
        PetriDish.validate_cells_compatible(cells)
    
    # Log what was preserved (if not already logged elsewhere)
    if 'gene_jsd' in year_data and len(cells) < 10:  # Only log for small test runs
        print(f"  Note: Snapshot contains gene_jsd data ({len(year_data['gene_jsd'])} values)")
    
    return cells


def save_petri_dish(petri: PetriDish, filepath: str, metadata: Optional[Dict] = None, 
                   include_cell_history: bool = False, include_gene_jsd: bool = False,
                   include_gene_metrics: bool = False,
                   compress: bool = True) -> None:
    """
    Save a PetriDish using Phase 1's save_history method.
    
    Args:
        petri: PetriDish object
        filepath: Output path
        metadata: Extra metadata dict (will be added to config)
        include_cell_history: Ignored - always included
        include_gene_jsd: Ignored - always included
        include_gene_metrics: Deprecated, ignored
        compress: If True, save as .json.gz; if False, save as .json
    """
    dir_path = os.path.dirname(filepath)
    if dir_path:  # Only create directory if there is one
        os.makedirs(dir_path, exist_ok=True)
    
    # Ensure filepath has correct extension
    if compress and not filepath.endswith('.gz'):
        if filepath.endswith('.json'):
            filepath = filepath + '.gz'
        else:
            filepath = filepath + '.json.gz'
    elif not compress and filepath.endswith('.gz'):
        filepath = filepath[:-3]  # Remove .gz
    
    # Convert to absolute path if not already
    abs_filepath = os.path.abspath(filepath)
    
    # Use Phase 1's save_history method
    # This saves config + complete history with cells and gene_jsd
    actual_path = petri.save_history(
        filename=abs_filepath,  # Pass the absolute filepath
        directory="",  # Empty since filepath is absolute
        compress=compress
    )
    
    # If we have Phase 2 specific metadata, we need to add it to the saved file
    if hasattr(petri, 'metadata') and petri.metadata:
        # Load, add metadata, save again
        with smart_open(actual_path, 'r') as f:
            data = json.load(f)
        
        # Add Phase 2 specific metadata to config
        data['config']['phase2_metadata'] = petri.metadata
        
        # Add individual_final with current state (post-mixing)
        data['individual_final'] = {
            'cells': [cell.to_dict() for cell in petri.cells],
            'gene_jsd': petri.calculate_gene_jsd()
        }
        
        with smart_open(actual_path, 'w') as f:
            if compress:
                json.dump(data, f, separators=(',', ':'))
            else:
                json.dump(data, f, indent=2)


def load_petri_dish(filepath: str, include_cell_history: bool = False, include_gene_jsd: bool = False) -> PetriDish:
    """
    Load a PetriDish from Phase 1 format file.
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
    
    # Handle control2 format (config + individual_final only, no history)
    if 'config' in data and 'individual_final' in data and 'history' not in data:
        # This is the simplified control2 format
        config = data['config']
        
        # Convert cells from individual_final
        cells = [dict_to_cell(cell_dict) for cell_dict in data['individual_final']['cells']]
        
        # Create PetriDish from cells
        petri = PetriDish.from_cells(
            cells,
            growth_phase=config.get('growth_phase'),
            metadata=config.get('phase2_metadata')
        )
        
        # Set the year from config
        petri.year = config.get('years', 0)
        
        # Restore Phase 2 metadata if present
        if 'phase2_metadata' in config:
            petri.metadata = config['phase2_metadata']
        
        return petri
    
    # Handle Phase 1 format (config + history)
    elif 'config' in data and 'history' in data:
        # This is Phase 1 format
        config = data['config']
        
        # Get the last year's cells to recreate PetriDish
        years = sorted([int(y) for y in data['history'].keys()])
        last_year = str(max(years))
        last_year_data = data['history'][last_year]
        
        # Convert cells from dict to Cell objects
        cells = [dict_to_cell(cell_dict) for cell_dict in last_year_data['cells']]
        
        # Create PetriDish from cells
        petri = PetriDish.from_cells(
            cells,
            growth_phase=config.get('growth_phase'),
            gene_rate_groups=config.get('gene_rate_groups'),
            seed=config.get('seed')
        )
        
        # Set the correct year
        petri.year = max(years)
        
        # Restore history if requested
        if include_cell_history:
            petri.cell_history = data['history']
        
        # Restore gene JSD history if requested
        if include_gene_jsd:
            petri.gene_jsd_history = {}
            for year_str, year_data in data['history'].items():
                if 'gene_jsd' in year_data:
                    petri.gene_jsd_history[year_str] = year_data['gene_jsd']
        
        # Restore Phase 2 metadata if present
        if 'phase2_metadata' in config:
            petri.metadata = config['phase2_metadata']
        
        return petri
    
    # Handle old Phase 2 format for backward compatibility
    # Load cells first
    cells = []
    if 'cells' in data and data['cells']:
        for cd in data['cells']:
            if 'gene_rate_groups' in cd:
                # New format with gene_rate_groups in each cell
                cells.append(dict_to_cell(cd))
            else:
                # Old format - get from metadata
                metadata = data.get('metadata', {})
                gene_rate_groups = metadata.get('gene_rate_groups')
                if not gene_rate_groups:
                    # Convert from rate if available
                    rate = metadata.get('rate')
                    if rate:
                        n = len(cd['cpg_sites'])
                        gene_size = metadata.get('gene_size', GENE_SIZE)
                        gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
                    else:
                        raise ValueError("PetriDish data has no rate configuration!")
                
                gene_size = metadata.get('gene_size', GENE_SIZE)
                cell = Cell.from_dict(cd, gene_rate_groups=gene_rate_groups, gene_size=gene_size)
                cells.append(cell)
    
    # Validate all cells have same gene_rate_groups
    if cells:
        PetriDish.validate_cells_compatible(cells)
        
        # Create PetriDish with gene_rate_groups from cells
        gene_rate_groups = cells[0].gene_rate_groups
        n = cells[0].n
        gene_size = cells[0].gene_size
        
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=gene_rate_groups,
            n=n,
            gene_size=gene_size,
            seed=None
        )
    else:
        # No cells - create empty PetriDish with metadata info
        metadata = data.get('metadata', {})
        gene_rate_groups = metadata.get('gene_rate_groups')
        if not gene_rate_groups:
            # Convert from rate if available
            rate = metadata.get('rate', 0.005)
            n = metadata.get('n', 1000)
            gene_size = metadata.get('gene_size', GENE_SIZE)
            gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
        else:
            gene_rate_groups = [tuple(group) for group in gene_rate_groups]
            n = metadata.get('n', 1000)
            gene_size = metadata.get('gene_size', GENE_SIZE)
        
        petri = PetriDish(
            gene_rate_groups=gene_rate_groups,
            n=n,
            gene_size=gene_size,
            seed=None
        )
    
    # Restore year if available
    metadata = data.get('metadata', {})
    if 'year' in metadata:
        petri.year = metadata['year']
    
    # Restore metadata
    petri.metadata = metadata
    
    # Restore cell history if requested and available
    if include_cell_history and 'cell_history' in data:
        petri.cell_history = data['cell_history']
    elif include_cell_history and 'history' in data:  # Backward compatibility
        petri.cell_history = data['history']
    else:
        # Clear any auto-initialized history if not loading it
        petri.cell_history = {}
    
    # Restore gene JSD history if requested and available
    if include_gene_jsd:
        # Use the helper to extract from either location
        gene_jsd_hist = extract_gene_jsd_from_history(data)
        petri.gene_jsd_history = gene_jsd_hist
    else:
        petri.gene_jsd_history = {}
    
    return petri


def extract_gene_jsd_from_history(sim_data: dict) -> dict:
    """
    Extract gene JSD data from simulation history structure.
    Reads from the correct location where phase1 stores the data (nested in history).
    
    Args:
        sim_data: Simulation data dictionary
        
    Returns:
        Dict with string years as keys, or empty dict if not found
    """
    gene_jsd_history = {}
    
    # Read from the correct location (nested in history)
    if 'history' in sim_data:
        for year_key, year_data in sim_data['history'].items():
            # Force to string to ensure consistency
            year_str = str(year_key)
            
            if 'gene_jsd' in year_data:
                gene_jsd_history[year_str] = year_data['gene_jsd']
    
    return gene_jsd_history


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
    sorted_cells = sorted(cells, key=lambda c: c.cell_jsd)
    
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
    
    # Update metadata to reflect new cell count
    if hasattr(petri, 'update_metadata'):
        petri.update_metadata({'num_cells': len(petri.cells)})
    
    return len(petri.cells)


def create_pure_snapshot_petri(snapshot_cells: List[Cell], n_cells: int = 5120, 
                              seed: Optional[int] = None) -> PetriDish:
    """
    Create a PetriDish with pure snapshot cells (for control2).
    Uses gene_rate_groups from the cells themselves.
    
    Args:
        snapshot_cells: Cells to sample from
        n_cells: Number of cells to sample
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
    
    if not sampled_cells:
        raise ValueError("No cells sampled")
    
    # Validate all cells have same gene_rate_groups
    PetriDish.validate_cells_compatible(sampled_cells)
    
    # Get configuration from first cell (all should be the same)
    first_cell = sampled_cells[0]
    
    # Create PetriDish with gene_rate_groups from cells
    petri = PetriDish(
        cells=sampled_cells,  # Pass ALL cells directly - no bogus history!
        gene_rate_groups=first_cell.gene_rate_groups,
        n=first_cell.n,
        gene_size=first_cell.gene_size,
        growth_phase=None  # Static population - won't be aged
    )
    
    # Set metadata
    petri.metadata = {
        'creation_method': 'pure_snapshot',
        'n_cells_sampled': n_cells,
        'num_cells': len(sampled_cells)
    }
    
    return petri


def create_control2_with_common_base(
    snapshot_cells: List[Cell],
    base_dir: str,
    target_size: int,
    seed: Optional[int] = None
) -> PetriDish:
    """
    Create Control2 individual with common base + additional snapshot cells.
    
    When using common mixing, Control2 shares the same base cells as
    mutant and control1, plus additional random snapshot cells if needed.
    Uses gene_rate_groups from the cells themselves.
    
    Args:
        snapshot_cells: Full second snapshot cells to sample additional from
        base_dir: Base directory to load common pool from
        target_size: Total number of cells needed
        seed: Random seed for additional sampling
    
    Returns:
        PetriDish with common base + additional unique snapshot cells
    """
    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)
    
    # Load the common pool from file
    common_pool = load_common_pool(base_dir)
    
    # Start with the common pool (deep copy to avoid reference issues)
    combined_cells = [copy.deepcopy(cell) for cell in common_pool]
    
    # Calculate how many additional cells needed
    n_additional = target_size - len(common_pool)
    
    if n_additional > 0:
        # Sample additional cells randomly from snapshot
        # Since we're not tracking indices anymore, just sample randomly
        additional_indices = random.sample(range(len(snapshot_cells)), 
                                         min(n_additional, len(snapshot_cells)))
        
        if n_additional > len(snapshot_cells):
            # Need more cells than available - error instead of sampling with replacement
            from .validation import ValidationError
            raise ValidationError(
                f"Insufficient snapshot cells for Control2: "
                f"need {n_additional} additional cells but only {len(snapshot_cells)} available. "
                f"This shouldn't happen with proper common mixing setup."
            )
        
        # Add the additional cells
        additional_cells = [copy.deepcopy(snapshot_cells[idx]) for idx in additional_indices]
        combined_cells.extend(additional_cells)
    elif n_additional < 0:
        # Edge case: common pool is larger than target (e.g., 100% mix ratio)
        # Just use subset of common pool
        combined_cells = combined_cells[:target_size]
    
    # Validate all cells have same gene_rate_groups
    if not combined_cells:
        raise ValueError("No cells in combined pool")
    
    PetriDish.validate_cells_compatible(combined_cells)
    
    # Get configuration from first cell (all should be the same)
    first_cell = combined_cells[0]
    
    # Create PetriDish with gene_rate_groups from cells
    petri = PetriDish(
        cells=combined_cells,  # Pass ALL cells directly - no bogus history!
        gene_rate_groups=first_cell.gene_rate_groups,
        n=first_cell.n,
        gene_size=first_cell.gene_size,
        growth_phase=None  # Static population - won't be aged
    )
    
    # Set metadata
    petri.metadata = {
        'creation_method': 'common_base_plus_snapshot',
        'common_base_size': len(common_pool),
        'additional_cells': n_additional if n_additional > 0 else 0,
        'target_size': target_size,
        'num_cells': len(combined_cells)
    }
    
    return petri


def grow_petri_for_years(petri: PetriDish, years: int, growth_phase: Optional[int] = None, 
                        track_history: bool = False, 
                        start_year: Optional[int] = None) -> None:
    """
    Grow a PetriDish for specified years using Phase 1's growth methods.
    
    Args:
        petri: PetriDish to grow (modified in place)
        years: Total number of years to simulate
        growth_phase: Years of exponential growth before homeostasis
        track_history: Ignored - history is always tracked now
        start_year: Ignored - PetriDish knows its starting year
    """
    if growth_phase is None:
        # Pure exponential growth
        petri.grow_exponentially(years, verbose=True)
    elif growth_phase >= years:
        # Only growth, no homeostasis
        petri.grow_exponentially(years, verbose=True)
    else:
        # Growth followed by homeostasis
        homeostasis_years = years - growth_phase
        target_population = 2 ** growth_phase
        
        petri.grow_exponentially(growth_phase, verbose=True)
        petri.maintain_homeostasis(homeostasis_years, target_population=target_population, verbose=True)


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
    
    cell_jsd_values = [cell.cell_jsd for cell in petri.cells]
    
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
    return np.array([cell.cell_jsd for cell in petri.cells])


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



def create_common_mixing_pool(snapshot_cells: List[Cell], 
                              normalized_size: int, 
                              mix_ratio: float,
                              seed: Optional[int] = None) -> List[Cell]:
    """
    Create a shared pool of snapshot cells for common mixing.
    
    Args:
        snapshot_cells: Available cells to sample from
        normalized_size: Normalized size of all individuals (after normalization)
        mix_ratio: Percentage from snapshot (0-1)
        seed: Random seed
        
    Returns:
        List of cells to be shared across all individuals
    
    Raises:
        ValidationError: If insufficient snapshot cells available
    """
    from .validation import ValidationError
    
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
    
    print(f"    Creating common mixing pool:")
    print(f"      Normalized individual size: {normalized_size} cells")
    print(f"      Target total size: {target_total} cells")
    print(f"      Snapshot cells needed: {n_snapshot_cells} cells")
    
    # Check if we have enough cells
    if n_snapshot_cells > len(snapshot_cells):
        raise ValidationError(
            f"Insufficient snapshot cells for uniform mixing: "
            f"need {n_snapshot_cells}, have {len(snapshot_cells)}. "
            f"Reduce mix_ratio (currently {mix_ratio:.1%}) or use a larger snapshot."
        )
    
    # Sample the cells
    indices = random.sample(range(len(snapshot_cells)), n_snapshot_cells)
    
    # Create deep copies
    pool = [copy.deepcopy(snapshot_cells[idx]) for idx in indices]
    
    print(f"      Created pool of {len(pool)} cells")
    
    return pool


def save_common_pool(pool: List[Cell], base_dir: str, compress: bool = False):
    """
    Save the common mixing pool to a file for reuse by Control2.
    
    Args:
        pool: List of Cell objects to save
        base_dir: Base directory for the phase2 run
        compress: Whether to compress the output file
    """
    import os
    
    # Create individuals directory if it doesn't exist
    individuals_dir = os.path.join(base_dir, 'individuals')
    os.makedirs(individuals_dir, exist_ok=True)
    
    # Determine filename
    ext = '.json.gz' if compress else '.json'
    pool_path = os.path.join(individuals_dir, f'common_pool{ext}')
    
    # Convert cells to dicts - just save the list directly
    pool_data = [cell.to_dict() for cell in pool]
    
    # Save
    print(f"    Saving common pool to: {pool_path}")
    with smart_open(pool_path, 'w') as f:
        json.dump(pool_data, f, indent=2)
    print(f"      Saved {len(pool_data)} cells to common pool")


def load_common_pool(base_dir: str) -> List[Cell]:
    """
    Load the common mixing pool from file.
    
    Args:
        base_dir: Base directory for the phase2 run
        
    Returns:
        List of Cell objects from the saved pool
        
    Raises:
        FileNotFoundError: If pool file doesn't exist
    """
    import os
    import glob
    
    # Look for pool file (could be compressed or not)
    individuals_dir = os.path.join(base_dir, 'individuals')
    pool_files = glob.glob(os.path.join(individuals_dir, 'common_pool.json*'))
    
    if not pool_files:
        raise FileNotFoundError(
            f"No common pool file found in {individuals_dir}. "
            f"Expected common_pool.json or common_pool.json.gz"
        )
    
    pool_path = pool_files[0]  # Take first match
    print(f"    Loading common pool from: {pool_path}")
    
    with smart_open(pool_path, 'r') as f:
        pool_data = json.load(f)
    
    # Convert dicts back to Cell objects (now just a list)
    pool = [dict_to_cell(cell_dict) for cell_dict in pool_data]
    print(f"      Loaded {len(pool)} cells from common pool")
    
    return pool


def validate_pool_compatibility(pool: List[Cell], target_cells: List[Cell]):
    """
    Validate that pool cells are compatible with target cells for mixing.
    
    Args:
        pool: Uniform pool cells
        target_cells: Cells to be mixed with (e.g., grown individuals)
        
    Raises:
        ValidationError: If ages or gene_rate_groups don't match
    """
    from .validation import ValidationError
    
    if not pool or not target_cells:
        return  # Nothing to validate if either is empty
    
    # Check pool internal consistency - all pool cells should have same age
    pool_ages = set(cell.age for cell in pool)
    if len(pool_ages) > 1:
        raise ValidationError(
            f"Pool cells have inconsistent ages: {sorted(pool_ages)}. "
            f"All pool cells should have the same age."
        )
    
    # Check gene_rate_groups consistency
    # Get first cell's gene_rate_groups as reference
    ref_groups = pool[0].gene_rate_groups if hasattr(pool[0], 'gene_rate_groups') else None
    
    if ref_groups is not None:
        # Check all pool cells have same gene_rate_groups
        for i, cell in enumerate(pool[1:], 1):
            if hasattr(cell, 'gene_rate_groups'):
                if cell.gene_rate_groups != ref_groups:
                    raise ValidationError(
                        f"Pool cell {i} has different gene_rate_groups than cell 0. "
                        f"All pool cells must have identical gene_rate_groups."
                    )
        
        # Check target cells have compatible gene_rate_groups
        for cell in target_cells:
            if hasattr(cell, 'gene_rate_groups'):
                if cell.gene_rate_groups != ref_groups:
                    raise ValidationError(
                        f"Target cells have incompatible gene_rate_groups with pool. "
                        f"Pool: {ref_groups}, Target: {cell.gene_rate_groups}"
                    )
    
    # Note: We don't validate ages between pool and target as they may differ
    # (e.g., grown cells have different ages than snapshot cells)


def mix_petri_common(petri: PetriDish, 
                     common_pool: List[Cell],
                     target_ratio: float) -> int:
    """
    Mix PetriDish with cells from common pool.
    All individuals have been normalized to same size, just add the entire pool.
    
    Args:
        petri: PetriDish to mix into (normalized to threshold size)
        common_pool: Pre-sampled shared cells (sized for the threshold-normalized individuals)
        target_ratio: What fraction should be from pool (used for 100% edge case only)
        
    Returns:
        Total cells after mixing
        
    Raises:
        ValidationError: If pool and petri cells are incompatible
    """
    # Validate compatibility before mixing
    validate_pool_compatibility(common_pool, petri.cells)
    
    # Handle edge case: 100% mix ratio means replace all cells
    if target_ratio >= 1.0:
        # Replace all cells with snapshot cells
        current_size = len(petri.cells)
        new_cells = [copy.deepcopy(cell) for cell in common_pool[:current_size]]
        random.shuffle(new_cells)
        petri.cells = new_cells
        # Update metadata to reflect change
        petri.update_metadata({'num_cells': len(new_cells)})
        return len(petri.cells)
    
    # Normal case: individuals are normalized to threshold, just add the entire pool
    added_cells = [copy.deepcopy(cell) for cell in common_pool]
    petri.cells.extend(added_cells)
    
    # Shuffle to mix
    random.shuffle(petri.cells)
    
    # Update metadata to reflect new cell count
    petri.update_metadata({'num_cells': len(petri.cells)})
    
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
            print(f"    Mutant {i+1:02d}: {current_size} cells - EXCLUDED (below threshold)")
        elif current_size == threshold_size:
            # Keep as is
            mutant_kept_as_is += 1
            normalized_mutant.append(dish)
            print(f"    Mutant {i+1:02d}: {current_size} cells - kept as is")
        else:
            # Trim to threshold size
            mutant_trimmed += 1
            # Create a deep copy to avoid modifying original
            new_dish = copy.deepcopy(dish)
            # Randomly sample cells to match threshold
            new_dish.cells = random.sample(new_dish.cells, threshold_size)
            normalized_mutant.append(new_dish)
            print(f"    Mutant {i+1:02d}: {current_size} → {threshold_size} cells (trimmed)")
    
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
            print(f"    Control1 {i+1:02d}: {current_size} cells - EXCLUDED (below threshold)")
        elif current_size == threshold_size:
            # Keep as is
            control1_kept_as_is += 1
            normalized_control1.append(dish)
            print(f"    Control1 {i+1:02d}: {current_size} cells - kept as is")
        else:
            # Trim to threshold size
            control1_trimmed += 1
            # Create a deep copy to avoid modifying original
            new_dish = copy.deepcopy(dish)
            # Randomly sample cells to match threshold
            new_dish.cells = random.sample(new_dish.cells, threshold_size)
            normalized_control1.append(new_dish)
            print(f"    Control1 {i+1:02d}: {current_size} → {threshold_size} cells (trimmed)")
    
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


def save_all_individuals(
    dishes: List[PetriDish],
    output_dir: str,
    batch_name: str,
    compress: bool = True,
    include_cell_history: bool = True,
    include_gene_metrics: bool = True
) -> None:
    """
    Save all individuals to disk after all processing complete.
    
    Args:
        dishes: List of PetriDish objects to save
        output_dir: Directory to save individuals
        batch_name: Name for logging ('mutant', 'control1', 'control2')
        compress: Whether to compress output files
        include_cell_history: Whether to include cell history in save
        include_gene_metrics: Whether to include gene metrics
    """
    if not dishes:
        print(f"  No {batch_name} individuals to save")
        return
        
    os.makedirs(output_dir, exist_ok=True)
    
    for petri in dishes:
        individual_id = petri.metadata.get('individual_id', 1)
        filepath = os.path.join(output_dir, f"individual_{individual_id:02d}.json")
        if compress:
            filepath += '.gz'
        
        save_petri_dish(
            petri, filepath,
            include_cell_history=include_cell_history,
            include_gene_metrics=include_gene_metrics,
            compress=compress
        )
    
    print(f"  Saved {len(dishes)} {batch_name} individuals to disk")