#!/usr/bin/env python3
"""
Extract snapshots from phase1 simulation.
This script extracts cells from specific years and saves them as snapshots.
"""

import argparse
import os
import sys
import json
import gzip
from typing import Optional, Dict, List, Tuple

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, rate_to_gene_rate_groups
from core.pipeline_utils import (
    load_snapshot_as_cells, save_snapshot_cells, smart_open
)
from core.validation import PipelineValidator, ValidationError


def extract_and_save_snapshot(
    simulation_path: str,
    year: int,
    output_path: str,
    gene_rate_groups: List[Tuple[int, float]],
    compress: bool = True
) -> List[Cell]:
    """
    Extract snapshot from simulation and save to disk.
    
    Args:
        simulation_path: Path to phase1 simulation file
        year: Year to extract
        output_path: Where to save snapshot
        gene_rate_groups: Expected gene rate groups for validation
        compress: Whether to compress output
    
    Returns:
        List of Cell objects
    """
    print(f"\nExtracting year {year} snapshot...")
    
    # Extract cells
    cells = load_snapshot_as_cells(
        simulation_path, year,
        expected_gene_rate_groups=gene_rate_groups
    )
    print(f"  Extracted {len(cells)} cells")
    
    # Save to disk
    save_snapshot_cells(cells, output_path, compress=compress)
    print(f"  Saved to {output_path}")
    
    return cells


def save_metadata(
    output_dir: str,
    gene_rate_groups: List[Tuple[int, float]],
    n_sites: int,
    gene_size: int,
    first_year: int,
    second_year: int,
    simulation_path: str
) -> None:
    """
    Save metadata about the extraction for other scripts to use.
    """
    metadata = {
        'gene_rate_groups': gene_rate_groups,
        'n_sites': n_sites,
        'gene_size': gene_size,
        'first_snapshot_year': first_year,
        'second_snapshot_year': second_year,
        'source_simulation': simulation_path,
        'extraction_timestamp': str(os.path.getmtime(simulation_path))
    }
    
    metadata_path = os.path.join(output_dir, 'metadata.json')
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"\nSaved extraction metadata to {metadata_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Extract snapshots from phase1 simulation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--simulation", required=True, type=str,
                       help="Path to phase1 simulation file")
    parser.add_argument("--first-snapshot", required=True, type=int,
                       help="First year to extract (e.g., 30)")
    parser.add_argument("--second-snapshot", required=True, type=int,
                       help="Second year to extract (e.g., 50)")
    parser.add_argument("--output-dir", required=True, type=str,
                       help="Base output directory")
    
    # Optional arguments
    parser.add_argument("--force-reload", action='store_true',
                       help="Force extraction even if snapshots exist")
    parser.add_argument("--no-compress", action='store_true',
                       help="Save uncompressed JSON files")
    
    args = parser.parse_args()
    
    # Determine compression
    use_compression = not args.no_compress
    ext = ".json.gz" if use_compression else ".json"
    
    # Create snapshots directory
    snapshots_dir = os.path.join(args.output_dir, "snapshots")
    os.makedirs(snapshots_dir, exist_ok=True)
    
    # Define output paths
    first_snapshot_path = os.path.join(snapshots_dir, f"year{args.first_snapshot}_snapshot{ext}")
    second_snapshot_path = os.path.join(snapshots_dir, f"year{args.second_snapshot}_snapshot{ext}")
    
    print("=" * 60)
    print("EXTRACT SNAPSHOTS")
    print("=" * 60)
    print(f"Simulation: {args.simulation}")
    print(f"Output directory: {snapshots_dir}")
    print(f"Years to extract: {args.first_snapshot}, {args.second_snapshot}")
    print(f"Compression: {'enabled' if use_compression else 'disabled'}")
    print("=" * 60)
    
    # Load simulation to get parameters
    print("\nLoading simulation parameters...")
    with smart_open(args.simulation, 'r') as f:
        sim_data = json.load(f)
    
    sim_params = sim_data.get('config', sim_data.get('parameters', {}))  # Support both old and new format
    n_sites = sim_params.get('n', 1000)
    gene_size = sim_params.get('gene_size', 5)
    
    # Get gene_rate_groups from simulation
    sim_gene_rate_groups = sim_params.get('gene_rate_groups')
    if not sim_gene_rate_groups:
        # Convert rate to gene_rate_groups for old simulations
        if 'rate' in sim_params:
            sim_gene_rate_groups = rate_to_gene_rate_groups(
                sim_params['rate'], n_sites, gene_size
            )
        else:
            raise ValueError("Simulation has no rate configuration!")
    
    # Convert to tuples
    gene_rate_groups = [tuple(group) for group in sim_gene_rate_groups]
    print(f"Gene rate groups: {gene_rate_groups}")
    print(f"Sites: {n_sites}, Gene size: {gene_size}")
    
    # Initialize validator
    validator = PipelineValidator(verbose=True)
    
    # Extract first snapshot
    if os.path.exists(first_snapshot_path) and not args.force_reload:
        print(f"\n✓ First snapshot already exists: {first_snapshot_path}")
        print("  Use --force-reload to re-extract")
    else:
        first_cells = extract_and_save_snapshot(
            args.simulation,
            args.first_snapshot,
            first_snapshot_path,
            gene_rate_groups,
            use_compression
        )
        
        # Validate
        try:
            validator.validate_snapshot(
                cells=first_cells,
                expected_year=args.first_snapshot,
                expected_gene_rate_groups=gene_rate_groups,
                min_cells=100  # Minimum needed for sampling
            )
            print("  ✓ Validation passed")
        except ValidationError as e:
            print(f"\n❌ Snapshot validation failed: {e}")
            sys.exit(1)
    
    # Extract second snapshot
    if os.path.exists(second_snapshot_path) and not args.force_reload:
        print(f"\n✓ Second snapshot already exists: {second_snapshot_path}")
        print("  Use --force-reload to re-extract")
    else:
        second_cells = extract_and_save_snapshot(
            args.simulation,
            args.second_snapshot,
            second_snapshot_path,
            gene_rate_groups,
            use_compression
        )
        
        # Validate
        try:
            validator.validate_snapshot(
                cells=second_cells,
                expected_year=args.second_snapshot,
                expected_gene_rate_groups=gene_rate_groups,
                min_cells=100
            )
            print("  ✓ Validation passed")
        except ValidationError as e:
            print(f"\n❌ Snapshot validation failed: {e}")
            sys.exit(1)
    
    # Save metadata
    save_metadata(
        snapshots_dir,
        gene_rate_groups,
        n_sites,
        gene_size,
        args.first_snapshot,
        args.second_snapshot,
        args.simulation
    )
    
    print("\n" + "=" * 60)
    print("EXTRACTION COMPLETE")
    print("=" * 60)
    print(f"Snapshots saved in: {snapshots_dir}")
    print("Files created:")
    print(f"  - year{args.first_snapshot}_snapshot{ext}")
    print(f"  - year{args.second_snapshot}_snapshot{ext}")
    print(f"  - metadata.json")


if __name__ == "__main__":
    main()