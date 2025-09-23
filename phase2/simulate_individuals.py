#!/usr/bin/env python3
"""
Create, grow, and mix individuals (mutant and control1).
This script handles the simulation of individual populations.
"""

import argparse
import os
import sys
import json
import numpy as np
import random
from typing import List, Dict, Optional, Tuple

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import PetriDish, Cell
from core.pipeline_utils import (
    load_snapshot_cells, sample_by_quantiles, sample_uniform,
    normalize_populations,
    create_uniform_mixing_pool, mix_petri_uniform, save_all_individuals
)
from core.individual_helpers import (
    create_individual, process_batch_growth, process_batch_mixing
)
from core.validation import PipelineValidator, ValidationError


def load_metadata(base_dir: str) -> Dict:
    """Load metadata from extraction stage."""
    metadata_path = os.path.join(base_dir, "snapshots", "metadata.json")
    with open(metadata_path, 'r') as f:
        return json.load(f)


def save_mixing_metadata(
    base_dir: str,
    mix_ratio: int,
    normalized_size: int  # No longer optional
) -> None:
    """Save mixing metadata for control2 creation."""
    metadata = {
        'uniform_mixing': True,  # Always true now
        'mix_ratio': mix_ratio,
        'normalized_size': normalized_size,
        'normalization_threshold': normalized_size  # Keep for backwards compat in Control2
    }
    
    metadata_path = os.path.join(base_dir, "individuals", "mixing_metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    print(f"\nSaved mixing metadata to {metadata_path}")


def create_mutant_individuals(
    snapshot_cells: List[Cell],
    n_quantiles: int,
    cells_per_quantile: int,
    growth_phase: int,
    seed: int
) -> List[PetriDish]:
    """Create mutant individuals using quantile sampling."""
    print(f"\n✓ Creating mutant individuals...")
    print(f"  Sampling by quantiles ({n_quantiles} quantiles, {cells_per_quantile} each)...")
    
    dishes = []
    sampled = sample_by_quantiles(
        snapshot_cells, 
        n_quantiles=n_quantiles,
        cells_per_quantile=cells_per_quantile,
        seed=seed
    )
    
    for i, (cell, quantile) in enumerate(sampled):
        file_index = i + 1
        additional_metadata = {
            'source_quantile': quantile,
            'n_quantiles': n_quantiles
        }
        
        petri = create_individual(
            cell=cell,
            individual_type='mutant',
            individual_id=file_index,
            growth_phase=growth_phase,
            additional_metadata=additional_metadata
        )
        dishes.append(petri)
    
    print(f"  Created {len(dishes)} mutant individuals")
    return dishes


def create_control1_individuals(
    snapshot_cells: List[Cell],
    n_individuals: int,
    growth_phase: int,
    seed: int
) -> List[PetriDish]:
    """Create control1 individuals using uniform sampling."""
    print(f"\n✓ Creating control1 individuals...")
    print(f"  Sampling uniformly...")
    
    dishes = []
    sampled_cells = sample_uniform(
        snapshot_cells,
        n_samples=n_individuals,
        seed=seed + 1000
    )
    
    for i, cell in enumerate(sampled_cells):
        file_index = i + 1
        additional_metadata = {'source': 'uniform'}
        
        petri = create_individual(
            cell=cell,
            individual_type='control1',
            individual_id=file_index,
            growth_phase=growth_phase,
            additional_metadata=additional_metadata
        )
        dishes.append(petri)
    
    print(f"  Created {len(dishes)} control1 individuals")
    return dishes


def calculate_timeline_metrics(first_snapshot: int, second_snapshot: int, growth_phase: int) -> Dict:
    """Calculate timeline metrics."""
    timeline_duration = second_snapshot - first_snapshot
    homeostasis_years = timeline_duration - growth_phase
    target_cells = 2 ** growth_phase
    
    return {
        'timeline_duration': timeline_duration,
        'homeostasis_years': homeostasis_years,
        'target_cells': target_cells
    }


def main():
    parser = argparse.ArgumentParser(
        description="Create, grow, and mix individuals",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Required arguments
    parser.add_argument("--base-dir", required=True, type=str,
                       help="Base directory with snapshots")
    parser.add_argument("--n-quantiles", required=True, type=int,
                       help="Number of quantiles for sampling")
    parser.add_argument("--cells-per-quantile", required=True, type=int,
                       help="Cells to sample per quantile")
    parser.add_argument("--growth-phase", required=True, type=int,
                       help="Years of exponential growth")
    parser.add_argument("--mix-ratio", required=True, type=int,
                       help="Percentage of second snapshot in mix (0-100)")
    
    # Optional arguments
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed for reproducibility")
    # Removed force-recreate since we always use timestamped directories
    
    args = parser.parse_args()
    
    # Set random seeds
    random.seed(args.seed)
    np.random.seed(args.seed)
    
    # Load metadata from extraction stage
    metadata = load_metadata(args.base_dir)
    gene_rate_groups = [tuple(g) for g in metadata['gene_rate_groups']]
    first_year = metadata['first_snapshot_year']
    second_year = metadata['second_snapshot_year']
    
    # Calculate timeline metrics
    metrics = calculate_timeline_metrics(first_year, second_year, args.growth_phase)
    timeline_duration = metrics['timeline_duration']
    homeostasis_years = metrics['homeostasis_years']
    expected_population = metrics['target_cells']
    
    # Calculate expected individuals
    expected_individuals = args.n_quantiles * args.cells_per_quantile
    
    print("=" * 60)
    print("SIMULATE INDIVIDUALS")
    print("=" * 60)
    print(f"Base directory: {args.base_dir}")
    print(f"Gene rate groups: {gene_rate_groups}")
    print(f"Timeline: year {first_year} → year {second_year}")
    print(f"Growth: {timeline_duration} years ({args.growth_phase} growth + {homeostasis_years} homeostasis)")
    print(f"Target population: {expected_population} cells")
    print(f"Mix ratio: {args.mix_ratio}% year {second_year}, {100-args.mix_ratio}% grown")
    print(f"Uniform mixing: Always enabled")
    print(f"Normalization: Always enabled (median - 0.5σ)")
    print(f"Seed: {args.seed}")
    print("=" * 60)
    
    # Create directories
    individuals_dir = os.path.join(args.base_dir, "individuals")
    mutant_dir = os.path.join(individuals_dir, "mutant")
    control1_dir = os.path.join(individuals_dir, "control1")
    
    for dir_path in [mutant_dir, control1_dir]:
        os.makedirs(dir_path, exist_ok=True)
    
    # Initialize validator
    validator = PipelineValidator()
    
    # Determine compression from snapshot files
    snapshots_dir = os.path.join(args.base_dir, "snapshots")
    first_snapshot_gz = os.path.join(snapshots_dir, f"year{first_year}_snapshot.json.gz")
    use_compression = os.path.exists(first_snapshot_gz)
    ext = ".json.gz" if use_compression else ".json"
    
    # ========================================================================
    # STAGE 3: Create Initial Individuals
    # ========================================================================
    print("\n" + "=" * 60)
    print("STAGE 3: Create Initial Individuals")
    print("=" * 60)
    
    # Load first snapshot
    first_snapshot_path = os.path.join(snapshots_dir, f"year{first_year}_snapshot{ext}")
    print(f"\nLoading first snapshot from {first_snapshot_path}...")
    first_snapshot_cells = load_snapshot_cells(first_snapshot_path)
    print(f"  Loaded {len(first_snapshot_cells)} cells")
    
    # Create mutant individuals
    mutant_dishes = create_mutant_individuals(
        first_snapshot_cells,
        args.n_quantiles,
        args.cells_per_quantile,
        args.growth_phase,
        args.seed
    )
    
    # Create control1 individuals
    control1_dishes = create_control1_individuals(
        first_snapshot_cells,
        expected_individuals,
        args.growth_phase,
        args.seed
    )
    
    # Validate initial individuals
    try:
        validator.validate_initial_individuals(
            mutant_dishes=mutant_dishes,
            control1_dishes=control1_dishes,
            expected_count=expected_individuals,
            expected_gene_rate_groups=gene_rate_groups,
            snapshot_year=first_year
        )
        print("  ✓ Validation passed")
    except ValidationError as e:
        print(f"\n❌ Initial individuals validation failed: {e}")
        sys.exit(1)
    
    # ========================================================================
    # STAGE 4: Grow Individuals
    # ========================================================================
    print("\n" + "=" * 60)
    print(f"STAGE 4: Grow Individuals ({timeline_duration} years)")
    print("=" * 60)
    
    # Grow mutant individuals
    process_batch_growth(
        dishes=mutant_dishes,
        batch_name='mutant',
        years=timeline_duration,
        growth_phase=args.growth_phase,
        expected_population=expected_population,
        start_year=first_year
    )
    
    # Grow control1 individuals
    process_batch_growth(
        dishes=control1_dishes,
        batch_name='control1',
        years=timeline_duration,
        growth_phase=args.growth_phase,
        expected_population=expected_population,
        start_year=first_year
    )
    
    # Validate grown individuals
    extinction_report = validator.validate_grown_individuals(
        mutant_dishes=mutant_dishes,
        control1_dishes=control1_dishes,
        expected_population=expected_population,
        growth_years=timeline_duration,
        allow_extinction=True
    )
    
    if extinction_report['total_extinct'] > 0:
        print(f"\n  ⚠️  Warning: {extinction_report['total_extinct']} individuals went extinct")
        if extinction_report['complete_batch_extinction']:
            print("\n❌ Fatal: Entire batch extinct")
            sys.exit(1)
    
    # ========================================================================
    # STAGE 5: Mix Populations
    # ========================================================================
    print("\n" + "=" * 60)
    print("STAGE 5: Mix Populations")
    print("=" * 60)
    
    # Load second snapshot
    second_snapshot_path = os.path.join(snapshots_dir, f"year{second_year}_snapshot{ext}")
    print(f"\nLoading second snapshot from {second_snapshot_path}...")
    second_snapshot_cells = load_snapshot_cells(second_snapshot_path)
    print(f"  Loaded {len(second_snapshot_cells)} cells")
    
    # Apply normalization (always mandatory now)
    print("\n  === APPLYING SIZE NORMALIZATION ===")
    print("  Using median - 0.5σ threshold")
    
    mutant_dishes, control1_dishes, normalized_size = normalize_populations(
        mutant_dishes, control1_dishes,
        seed=args.seed + 5000
    )
    
    # Update metadata
    for new_id, dish in enumerate(mutant_dishes, 1):
        if not hasattr(dish, 'metadata'):
            dish.metadata = {}
        dish.metadata['individual_id'] = new_id
        dish.metadata['normalized'] = True
        dish.metadata['normalization_threshold'] = normalized_size
    
    for new_id, dish in enumerate(control1_dishes, 1):
        if not hasattr(dish, 'metadata'):
            dish.metadata = {}
        dish.metadata['individual_id'] = new_id
        dish.metadata['normalized'] = True
        dish.metadata['normalization_threshold'] = normalized_size
    
    print(f"  Normalized all individuals to {normalized_size} cells")
    
    # Create uniform pool
    print("\n  === UNIFORM MIXING ===")
    print(f"  Creating uniform mixing pool...")
    uniform_pool = create_uniform_mixing_pool(
        second_snapshot_cells,
        normalized_size,
        args.mix_ratio / 100,
        seed=args.seed + 1000
    )
    
    # Save the uniform pool to file
    from phase2.core.pipeline_utils import save_uniform_pool
    save_uniform_pool(uniform_pool, args.base_dir, compress=use_compression)
    
    # Validate the pool
    validator.validate_uniform_pool(uniform_pool, len(uniform_pool), second_year)
    
    # Mix all individuals
    print(f"\n  Mixing with uniform pool...")
    
    for i, petri in enumerate(mutant_dishes):
        individual_id = i + 1
        initial_size = len(petri.cells)
        final_size = mix_petri_uniform(petri, uniform_pool, args.mix_ratio / 100)
        print(f"    Mutant {individual_id:02d}: {initial_size} → {final_size} cells")
        
        if not hasattr(petri, 'metadata'):
            petri.metadata = {}
        petri.metadata['mixed'] = True
        petri.metadata['mix_mode'] = 'uniform'
        petri.metadata['mix_ratio'] = args.mix_ratio
    
    for i, petri in enumerate(control1_dishes):
        individual_id = i + 1
        initial_size = len(petri.cells)
        final_size = mix_petri_uniform(petri, uniform_pool, args.mix_ratio / 100)
        print(f"    Control1 {individual_id:02d}: {initial_size} → {final_size} cells")
        
        if not hasattr(petri, 'metadata'):
            petri.metadata = {}
        petri.metadata['mixed'] = True
        petri.metadata['mix_mode'] = 'uniform'
        petri.metadata['mix_ratio'] = args.mix_ratio
    
    # Validate mixed populations
    try:
        validator.validate_mixed_populations(
            mutant_dishes=mutant_dishes,
            control1_dishes=control1_dishes,
            mix_ratio=args.mix_ratio,
            uniform_mixing=True,  # Always true now
            second_snapshot_year=second_year
        )
        print("  ✓ Mixing validation passed")
    except ValidationError as e:
        print(f"\n❌ Mixing validation failed: {e}")
        sys.exit(1)
    
    # ========================================================================
    # Save All Individuals
    # ========================================================================
    print("\n  === SAVING ALL INDIVIDUALS ===")
    
    # Save mutant individuals
    save_all_individuals(
        dishes=mutant_dishes,
        output_dir=mutant_dir,
        batch_name='mutant',
        compress=use_compression,
        include_cell_history=True,
        include_gene_metrics=True
    )
    
    # Save control1 individuals
    save_all_individuals(
        dishes=control1_dishes,
        output_dir=control1_dir,
        batch_name='control1',
        compress=use_compression,
        include_cell_history=True,
        include_gene_metrics=True
    )
    
    # Save mixing metadata
    save_mixing_metadata(
        args.base_dir,
        args.mix_ratio,
        normalized_size
    )
    
    print("\n" + "=" * 60)
    print("SIMULATION COMPLETE")
    print("=" * 60)
    print(f"Created {len(mutant_dishes)} mutant individuals")
    print(f"Created {len(control1_dishes)} control1 individuals")
    print(f"All individuals saved to: {individuals_dir}")


if __name__ == "__main__":
    main()