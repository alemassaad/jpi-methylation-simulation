#!/usr/bin/env python3
"""
Create control2 individuals (pure second snapshot).
This script creates control individuals from the second snapshot without growth.
"""

import argparse
import os
import sys
import json
import glob
import numpy as np
from typing import List, Dict, Optional

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import PetriDish, Cell
from core.pipeline_utils import (
    load_snapshot_cells, create_pure_snapshot_petri,
    create_control2_with_common_base, save_petri_dish,
    check_petri_files_state
)
from core.validation import PipelineValidator, ValidationError


def load_snapshots_metadata(base_dir: str) -> Dict:
    """Load metadata from extraction stage."""
    metadata_path = os.path.join(base_dir, "snapshots", "metadata.json")
    with open(metadata_path, 'r') as f:
        return json.load(f)


def load_mixing_metadata(base_dir: str) -> Optional[Dict]:
    """Load mixing metadata if it exists."""
    metadata_path = os.path.join(base_dir, "individuals", "mixing_metadata.json")
    if os.path.exists(metadata_path):
        with open(metadata_path, 'r') as f:
            return json.load(f)
    return None


def determine_control2_count(base_dir: str) -> int:
    """Determine how many control2 individuals to create."""
    individuals_dir = os.path.join(base_dir, "individuals")
    
    # Count existing individuals
    mutant_count = len(glob.glob(os.path.join(individuals_dir, "mutant", "individual_*.json*")))
    control1_count = len(glob.glob(os.path.join(individuals_dir, "control1", "individual_*.json*")))
    
    # After normalization, use average
    num_control2 = (mutant_count + control1_count) // 2
    print(f"  Control2 count: {num_control2}")
    print(f"    (Average of {mutant_count} mutant + {control1_count} control1)")
    
    return num_control2


def calculate_expected_final_cells(mixing_metadata: Optional[Dict]) -> int:
    """Calculate expected final cell count based on mixing parameters."""
    if not mixing_metadata:
        # Default if no mixing metadata
        return 640
    
    # Common mixing is always true now
    normalized_size = mixing_metadata.get('normalized_size', 128)
    mix_ratio = mixing_metadata.get('mix_ratio', 80) / 100
    expected_final_cells = int(normalized_size / (1 - mix_ratio))
    
    return expected_final_cells


def main():
    parser = argparse.ArgumentParser(
        description="Create control2 individuals from second snapshot",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--base-dir", required=True, type=str,
                       help="Base directory with snapshots and individuals")
    
    # Optional arguments
    parser.add_argument("--seed", type=int, default=342,
                       help="Random seed for reproducibility (default: 342)")
    # Removed force-recreate since we always use timestamped directories
    
    args = parser.parse_args()
    
    # Set random seed
    np.random.seed(args.seed)
    
    # Load metadata
    snapshots_metadata = load_snapshots_metadata(args.base_dir)
    mixing_metadata = load_mixing_metadata(args.base_dir)
    
    gene_rate_groups = [tuple(g) for g in snapshots_metadata['gene_rate_groups']]
    second_year = snapshots_metadata['second_snapshot_year']
    
    print("=" * 60)
    print("CREATE CONTROL2 INDIVIDUALS")
    print("=" * 60)
    print(f"Base directory: {args.base_dir}")
    print(f"Gene rate groups: {gene_rate_groups}")
    print(f"Second snapshot year: {second_year}")
    print(f"Seed: {args.seed}")
    
    if mixing_metadata:
        print(f"Mixing mode: Common (always enabled)")
        print(f"Mix ratio: {mixing_metadata['mix_ratio']}%")
        if mixing_metadata.get('normalized_size'):
            print(f"Normalized size: {mixing_metadata['normalized_size']}")
    print("=" * 60)
    
    # Create control2 directory
    control2_dir = os.path.join(args.base_dir, "individuals", "control2")
    os.makedirs(control2_dir, exist_ok=True)
    
    # Determine compression from existing files
    snapshots_dir = os.path.join(args.base_dir, "snapshots")
    second_snapshot_gz = os.path.join(snapshots_dir, f"year{second_year}_snapshot.json.gz")
    use_compression = os.path.exists(second_snapshot_gz)
    ext = ".json.gz" if use_compression else ".json"
    
    # Determine how many control2 individuals to create
    num_control2 = determine_control2_count(args.base_dir)
    
    # Calculate expected final cells
    expected_final_cells = calculate_expected_final_cells(mixing_metadata)
    print(f"  Expected final cells: {expected_final_cells}")
    
    # Check if control2 already exists
    control2_state = check_petri_files_state(control2_dir, expected_final_cells)
    
    if control2_state['total_files'] >= num_control2:
        print(f"\n✓ Control2 individuals already exist ({control2_state['total_files']} files)")
        print("  Skipping creation")
        return
    
    print(f"\n✓ Creating {num_control2} control2 individuals...")
    
    # Clean up any existing control2 files
    for old_file in glob.glob(os.path.join(control2_dir, f"individual_*{ext}")):
        os.remove(old_file)
    
    # Load second snapshot
    second_snapshot_path = os.path.join(snapshots_dir, f"year{second_year}_snapshot{ext}")
    print(f"\nLoading second snapshot from {second_snapshot_path}...")
    second_snapshot_cells = load_snapshot_cells(second_snapshot_path)
    print(f"  Loaded {len(second_snapshot_cells)} cells")
    
    # Adjust target size if snapshot has fewer cells
    actual_control2_size = min(expected_final_cells, len(second_snapshot_cells))
    if actual_control2_size < expected_final_cells:
        print(f"  Warning: Snapshot has only {len(second_snapshot_cells)} cells")
        print(f"  Adjusting control2 size to {actual_control2_size}")
    
    # Create control2 individuals
    control2_dishes = []
    
    # Note: common pool is now loaded from file inside create_control2_with_common_base
    print(f"  Using common base from saved pool file")
    
    for i in range(num_control2):
        print(f"  Creating individual {i+1}/{num_control2}")
        
        file_index = i + 1
        
        # Create PetriDish using saved common pool
        # The function will load the pool from file internally
        petri = create_control2_with_common_base(
            second_snapshot_cells,
            args.base_dir,  # Pass base_dir instead of pool and indices
            actual_control2_size,
            seed=args.seed + i
        )
        
        # Add minimal metadata
        if not hasattr(petri, 'metadata'):
            petri.metadata = {}
        petri.metadata.update({
            'individual_id': file_index,
            'individual_type': 'control2',
            'mix_ratio': mixing_metadata.get('mix_ratio') if mixing_metadata else 0,
            'normalized_size': mixing_metadata.get('normalized_size') if mixing_metadata else 0
        })
        
        control2_dishes.append(petri)
        
        # Save without history (control2 doesn't need it)
        filepath = os.path.join(control2_dir, f"individual_{file_index:02d}{ext}")
        save_petri_dish(
            petri, filepath,
            include_cell_history=False,  # No history needed
            include_gene_metrics=False,  # No longer saving gene metrics
            compress=use_compression
        )
    
    print(f"\n✓ Created and saved {len(control2_dishes)} control2 individuals")
    
    # Validate control2
    validator = PipelineValidator()
    try:
        validator.validate_control2(
            control2_dishes=control2_dishes,
            second_snapshot_year=second_year,
            expected_count=num_control2
        )
        print("  ✓ Validation passed")
    except ValidationError as e:
        print(f"\n❌ Control2 validation failed: {e}")
        sys.exit(1)
    
    print("\n" + "=" * 60)
    print("CONTROL2 CREATION COMPLETE")
    print("=" * 60)
    print(f"Created {num_control2} control2 individuals")
    print(f"Saved to: {control2_dir}")


if __name__ == "__main__":
    main()