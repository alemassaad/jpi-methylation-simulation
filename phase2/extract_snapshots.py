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
from typing import Optional, Dict, List, Tuple, Any

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, rate_to_gene_rate_groups
from core.pipeline_utils import smart_open
from core.validation import PipelineValidator, ValidationError


def extract_and_save_snapshot(
    simulation_path: str,
    year: int,
    output_path: str,
    compress: bool = True
) -> Dict:
    """
    Extract snapshot from simulation and save as direct copy with year key.
    
    The snapshot format is {"year": {data}} where data is an exact copy
    from the phase1 simulation's history[year].
    
    Args:
        simulation_path: Path to phase1 simulation file
        year: Year to extract
        output_path: Where to save snapshot
        compress: Whether to compress output
    
    Returns:
        Dictionary with year data (cells and optionally gene_jsd)
    """
    print(f"\nExtracting year {year} snapshot...")
    
    # Load simulation and extract year data directly
    with smart_open(simulation_path, 'r') as f:
        data = json.load(f)
    
    year_str = str(year)
    history = data['history']
    
    if year_str not in history:
        available_years = sorted([int(y) for y in history.keys()])
        raise ValueError(f"Year {year} not found. Available: {available_years}")
    
    # Get the year data directly (no transformation)
    year_data = history[year_str]
    print(f"  Extracted {len(year_data.get('cells', []))} cells")
    
    # Save year data with the year key included (like in phase1)
    snapshot_with_key = {year_str: year_data}
    with smart_open(output_path, 'w') as f:
        json.dump(snapshot_with_key, f, indent=2)
    print(f"  Saved to {output_path}")
    
    return year_data


def save_metadata(
    output_dir: str,
    simulation_path: str,
    first_year: int,
    second_year: int
) -> None:
    """
    Save metadata about the extraction for other scripts to use.
    Gets configuration from simulation file.
    """
    # Load simulation config
    with smart_open(simulation_path, 'r') as f:
        data = json.load(f)
    
    config = data.get('config', data.get('parameters', {}))
    
    # Extract gene_rate_groups
    gene_rate_groups = config.get('gene_rate_groups')
    if not gene_rate_groups:
        # Convert from rate for old format
        if 'rate' in config:
            gene_rate_groups = rate_to_gene_rate_groups(
                config['rate'], 
                config.get('n', 1000), 
                config.get('gene_size', 5)
            )
        else:
            raise ValueError("Simulation has no rate configuration!")
    
    metadata = {
        'gene_rate_groups': gene_rate_groups,
        'n_sites': config.get('n', 1000),
        'gene_size': config.get('gene_size', 5),
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
    # Compression options (mutually exclusive)
    compress_group = parser.add_mutually_exclusive_group()
    compress_group.add_argument("--compress", action='store_true', dest='compress',
                               help="Compress output files (.json.gz)")
    compress_group.add_argument("--no-compress", action='store_false', dest='compress',
                               help="Don't compress output files (.json)")
    parser.set_defaults(compress=False)  # Default to no compression
    
    args = parser.parse_args()
    
    # Determine compression
    use_compression = args.compress
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
    
    # Extract first snapshot
    if os.path.exists(first_snapshot_path):
        print(f"\n✓ First snapshot already exists: {first_snapshot_path}")
        print("  Skipping extraction (using cached file)")
    else:
        year_data = extract_and_save_snapshot(
            args.simulation,
            args.first_snapshot,
            first_snapshot_path,
            use_compression
        )
        print("  ✓ Extraction complete")
        if 'gene_jsd' in year_data:
            print(f"  ✓ Preserved gene_jsd data ({len(year_data['gene_jsd'])} values)")
    
    # Extract second snapshot
    if os.path.exists(second_snapshot_path):
        print(f"\n✓ Second snapshot already exists: {second_snapshot_path}")
        print("  Skipping extraction (using cached file)")
    else:
        year_data = extract_and_save_snapshot(
            args.simulation,
            args.second_snapshot,
            second_snapshot_path,
            use_compression
        )
        print("  ✓ Extraction complete")
        if 'gene_jsd' in year_data:
            print(f"  ✓ Preserved gene_jsd data ({len(year_data['gene_jsd'])} values)")
    
    # Save metadata
    save_metadata(
        snapshots_dir,
        args.simulation,
        args.first_snapshot,
        args.second_snapshot
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