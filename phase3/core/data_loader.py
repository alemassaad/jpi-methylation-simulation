#!/usr/bin/env python3
"""
Data loading functions for phase3 analysis.
Minimal subset of functions from phase2's pipeline_utils.py needed for reading data.
"""

import json
import gzip
import os
import glob
from typing import List, Dict, Any, Optional, Tuple, Union, TextIO

import sys
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish


def smart_open(filepath: str, mode: str = 'r') -> Union[TextIO, gzip.GzipFile]:
    """
    Open a file, automatically handling gzip compression based on extension.
    
    Args:
        filepath: Path to file
        mode: Open mode ('r', 'w', etc.)
    
    Returns:
        File handle (either regular or gzip)
    """
    if filepath.endswith('.gz'):
        if 'b' not in mode and 't' not in mode:
            mode = mode + 't'  # Default to text mode for gzip
        return gzip.open(filepath, mode)
    else:
        return open(filepath, mode)


def dict_to_cell(cell_dict: Dict[str, Any]) -> Cell:
    """
    Convert a dictionary to a Cell object.
    Handles both old and new JSON formats.
    
    Args:
        cell_dict: Dictionary containing cell data
        
    Returns:
        Cell object
    """
    # Handle new lean format
    if 'cpg_sites' in cell_dict:
        methylated = cell_dict['cpg_sites']
    else:
        methylated = cell_dict.get('methylated', [])
    
    # Get age
    age = cell_dict.get('age', 0)
    
    # Get gene_rate_groups if available
    gene_rate_groups = None
    if 'gene_rate_groups' in cell_dict:
        gene_rate_groups = [tuple(g) for g in cell_dict['gene_rate_groups']]
    
    # Create cell with appropriate configuration
    n = len(methylated)
    if gene_rate_groups:
        # Calculate gene_size from gene_rate_groups and n
        total_genes_in_groups = sum(g[0] for g in gene_rate_groups)
        gene_size = n // total_genes_in_groups if total_genes_in_groups > 0 else 5
        cell = Cell(gene_rate_groups=gene_rate_groups, n=n, gene_size=gene_size)
    else:
        # Default to a reasonable rate if not specified
        cell = Cell(rate=0.005, n=n)
    
    cell.methylated = methylated
    cell.age = age
    
    # Set cell_jsd if available
    if 'cell_jsd' in cell_dict:
        cell.cell_jsd = cell_dict['cell_jsd']
    elif 'cell_JSD' in cell_dict:  # Handle old format
        cell.cell_jsd = cell_dict['cell_JSD']
    
    return cell


def load_snapshot_cells(filepath: str) -> List[Cell]:
    """
    Load snapshot cells from a JSON file.
    
    Args:
        filepath: Path to snapshot JSON file
        
    Returns:
        List of Cell objects
    """
    with smart_open(filepath, 'r') as f:
        data = json.load(f)
    
    cells = []
    for cell_dict in data['cells']:
        cell = dict_to_cell(cell_dict)
        cells.append(cell)
    
    return cells


def load_petri_dish(filepath: str, include_cell_history: bool = False) -> PetriDish:
    """
    Load a PetriDish from a JSON file.
    
    Args:
        filepath: Path to PetriDish JSON file
        include_cell_history: Whether to load cell history
        
    Returns:
        PetriDish object
    """
    with smart_open(filepath, 'r') as f:
        data = json.load(f)
    
    # Create PetriDish with appropriate parameters
    cells = [dict_to_cell(cell_dict) for cell_dict in data['cells']]
    
    if cells:
        # Use first cell's configuration
        petri = PetriDish()
        petri.cells = cells
        
        # Set gene_rate_groups from first cell if available
        if hasattr(cells[0], 'gene_rate_groups') and cells[0].gene_rate_groups:
            petri.gene_rate_groups = cells[0].gene_rate_groups
    else:
        petri = PetriDish()
    
    # Load metadata if present
    if 'metadata' in data:
        petri.metadata = data['metadata']
    
    # Load history if requested and available
    if include_cell_history and 'cell_history' in data:
        petri.cell_history = data['cell_history']
        
    # Load gene JSD histories if available
    if 'gene_jsd_history' in data:
        petri.gene_jsd_history = data['gene_jsd_history']
    if 'mean_gene_jsd_history' in data:
        petri.mean_gene_jsd_history = data['mean_gene_jsd_history']
    if 'median_gene_jsd_history' in data:
        petri.median_gene_jsd_history = data['median_gene_jsd_history']
    
    return petri


def load_all_petri_dishes(directory: str, include_cell_history: bool = False) -> List[PetriDish]:
    """
    Load all PetriDish objects from a directory.
    
    Args:
        directory: Directory containing individual JSON files
        include_cell_history: Whether to load cell history
        
    Returns:
        List of PetriDish objects
    """
    dishes = []
    
    # Find all individual files
    json_files = sorted(glob.glob(os.path.join(directory, "individual_*.json*")))
    
    for filepath in json_files:
        try:
            petri = load_petri_dish(filepath, include_cell_history)
            dishes.append(petri)
        except Exception as e:
            print(f"Warning: Could not load {filepath}: {e}")
    
    return dishes


def extract_gene_jsd_from_history(sim_data: dict) -> Tuple[Dict, Dict, Dict]:
    """
    Extract gene JSD history from simulation data.
    
    Args:
        sim_data: Simulation data dictionary
        
    Returns:
        Tuple of (gene_jsd_history, mean_history, median_history)
    """
    import numpy as np
    
    gene_jsd_history = {}
    mean_gene_jsd_history = {}
    median_gene_jsd_history = {}
    
    history = sim_data.get('history', {})
    
    for year_str, year_data in history.items():
        if 'gene_jsd' in year_data:
            gene_jsds = year_data['gene_jsd']
            # Use string keys for consistency
            gene_jsd_history[year_str] = gene_jsds
            mean_gene_jsd_history[year_str] = np.mean(gene_jsds)
            median_gene_jsd_history[year_str] = np.median(gene_jsds)
    
    return gene_jsd_history, mean_gene_jsd_history, median_gene_jsd_history


def load_phase2_metadata(base_dir: str) -> Dict[str, Any]:
    """
    Load metadata files from phase2 output.
    
    Args:
        base_dir: Base directory of phase2 output
        
    Returns:
        Dictionary with metadata
    """
    metadata = {}
    
    # Load snapshots metadata
    snapshots_meta_path = os.path.join(base_dir, "snapshots", "metadata.json")
    if os.path.exists(snapshots_meta_path):
        with open(snapshots_meta_path, 'r') as f:
            metadata['snapshots'] = json.load(f)
    
    # Load mixing metadata
    mixing_meta_path = os.path.join(base_dir, "individuals", "mixing_metadata.json")
    if os.path.exists(mixing_meta_path):
        with open(mixing_meta_path, 'r') as f:
            metadata['mixing'] = json.load(f)
    
    return metadata