#!/usr/bin/env python3
"""
Minimal data extraction utilities for phase3.
Contains only the essential functions used by other modules.
"""

import gzip
import sys
import os
from typing import Dict, Any

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell


def smart_open(filepath: str, mode: str = 'r'):
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

    Args:
        cell_dict: Dictionary containing cell data

    Returns:
        Cell object
    """
    # Extract cpg_sites from cell dictionary
    if 'cpg_sites' not in cell_dict:
        raise ValueError("Cell dictionary must contain 'cpg_sites'")
    cpg_sites_data = cell_dict['cpg_sites']

    # Extract gene_rate_groups if present
    gene_rate_groups = cell_dict.get('gene_rate_groups', None)

    # Create Cell object
    n = len(cpg_sites_data)

    if gene_rate_groups:
        # Convert to list of tuples if needed
        if gene_rate_groups and isinstance(gene_rate_groups[0], list):
            gene_rate_groups = [tuple(g) for g in gene_rate_groups]
        cell = Cell(gene_rate_groups=gene_rate_groups, n=n)
    else:
        # Default cell with uniform rate
        cell = Cell(rate=0.005, n=n)

    # Set methylation state
    cell.cpg_sites = cpg_sites_data  # This is the actual methylation array

    # Set other attributes if present
    if 'cell_jsd' in cell_dict:
        cell.cell_jsd = cell_dict['cell_jsd']
    if 'age' in cell_dict:
        cell.age = cell_dict['age']

    return cell