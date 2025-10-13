#!/usr/bin/env python
"""
Phase 2 Pipeline - Single-file data generation pipeline.

This consolidated pipeline replaces the previous multi-script architecture
with direct function calls for better performance and maintainability.
"""

import argparse
import os
import sys
import json
import gzip
import glob
import random
import copy
import time
from datetime import datetime
from typing import List, Dict, Tuple, Optional, Any, Union, TextIO

# Import dependencies
import numpy as np

# Try to import YAML support
try:
    import yaml
    HAS_YAML = True
except ImportError:
    HAS_YAML = False
    print("Warning: PyYAML not installed. Config file support disabled.")
    print("Install with: pip install pyyaml")

# Add parent directory to path for phase1 imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))

from cell import Cell, PetriDish, rate_to_gene_rate_groups, GENE_SIZE


# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

def smart_open(filepath: str, mode: str = 'r') -> Union[TextIO, gzip.GzipFile]:
    """Open file based on extension (.json or .json.gz)."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, mode + 't' if 't' not in mode else mode)
    elif filepath.endswith('.json'):
        return open(filepath, mode)
    else:
        raise ValueError(f"Unsupported extension: {filepath}")


def load_snapshot_cells(filepath: str) -> List[Cell]:
    """Load cells from snapshot file with year key wrapper."""
    with smart_open(filepath, 'r') as f:
        data = json.load(f)

    # Handle year key wrapper (e.g., {"30": {...}})
    if len(data) == 1 and isinstance(list(data.keys())[0], str):
        year_key = list(data.keys())[0]
        if year_key.isdigit():
            year_data = data[year_key]
        else:
            year_data = data
    else:
        year_data = data

    cell_dicts = year_data['cells']

    # Try to load metadata for fallback gene_rate_groups
    snapshot_dir = os.path.dirname(filepath)
    metadata_file = os.path.join(snapshot_dir, 'metadata.json')
    gene_rate_groups_fallback = None
    gene_size_fallback = GENE_SIZE

    if os.path.exists(metadata_file):
        with open(metadata_file, 'r') as f:
            metadata = json.load(f)
            if 'gene_rate_groups' in metadata:
                gene_rate_groups_fallback = [tuple(g) for g in metadata['gene_rate_groups']]
            gene_size_fallback = metadata.get('gene_size', GENE_SIZE)

    # Convert to Cell objects
    cells = []
    for cd in cell_dicts:
        if 'gene_rate_groups' in cd:
            # New format
            gene_rate_groups = [tuple(g) for g in cd['gene_rate_groups']]
            cell = Cell(
                n=len(cd['cpg_sites']),
                gene_rate_groups=gene_rate_groups,
                gene_size=cd.get('gene_size', GENE_SIZE)
            )
            cell.cpg_sites = cd['cpg_sites']
            cell.age = cd.get('age', 0)
            cell.cell_jsd = cd.get('cell_jsd', cd.get('cell_JSD', 0.0))
        else:
            # Old format - use fallback
            if not gene_rate_groups_fallback:
                raise ValueError("Cells missing gene_rate_groups and no metadata.json found!")
            cell = Cell.from_dict(cd, gene_rate_groups=gene_rate_groups_fallback,
                                gene_size=gene_size_fallback)
        cells.append(cell)

    # Validate compatibility
    if cells:
        PetriDish.validate_cells_compatible(cells)

    return cells


def sample_by_quantiles(cells: List[Cell], n_quantiles: int,
                       cells_per_quantile: int, seed: int) -> List[Tuple[Cell, int]]:
    """Sample cells from quantiles based on JSD."""
    random.seed(seed)
    np.random.seed(seed)

    # Sort by JSD
    sorted_cells = sorted(cells, key=lambda c: c.cell_jsd)
    n_cells = len(sorted_cells)
    quantile_size = n_cells // n_quantiles

    sampled = []
    for q in range(n_quantiles):
        start_idx = q * quantile_size
        end_idx = n_cells if q == n_quantiles - 1 else (q + 1) * quantile_size
        quantile_cells = sorted_cells[start_idx:end_idx]

        # Sample from this quantile
        if len(quantile_cells) >= cells_per_quantile:
            sample_indices = random.sample(range(len(quantile_cells)), cells_per_quantile)
            for idx in sample_indices:
                sampled.append((copy.deepcopy(quantile_cells[idx]), q))
        else:
            for cell in quantile_cells:
                sampled.append((copy.deepcopy(cell), q))

    return sampled


def sample_uniform(cells: List[Cell], n_samples: int, seed: int) -> List[Cell]:
    """Sample cells uniformly."""
    random.seed(seed)
    np.random.seed(seed)

    if n_samples > len(cells):
        raise ValueError(f"Cannot sample {n_samples} from {len(cells)} cells")

    sample_indices = random.sample(range(len(cells)), n_samples)
    return [copy.deepcopy(cells[idx]) for idx in sample_indices]


def normalize_populations(mutant_dishes: List[PetriDish],
                         control1_dishes: List[PetriDish],
                         seed: int) -> Tuple[List[PetriDish], List[PetriDish], int]:
    """Normalize populations using median - 0.5σ threshold."""
    random.seed(seed)
    np.random.seed(seed)

    all_dishes = mutant_dishes + control1_dishes
    all_sizes = [len(dish.cells) for dish in all_dishes]

    if not all_sizes:
        return [], [], 0

    if len(all_sizes) == 1:
        return mutant_dishes, control1_dishes, all_sizes[0]

    # Calculate threshold
    median_size = np.median(all_sizes)
    std_size = np.std(all_sizes)
    threshold_size = int(median_size - 0.5 * std_size)

    if threshold_size < 1:
        threshold_size = min(all_sizes)

    print(f"\n  Normalization: median={median_size:.1f}, σ={std_size:.1f}, threshold={threshold_size}")

    # Process mutant
    normalized_mutant = []
    for i, dish in enumerate(mutant_dishes):
        size = len(dish.cells)
        if size < threshold_size:
            print(f"    Mutant {i+1:02d}: {size} cells - EXCLUDED")
        elif size == threshold_size:
            normalized_mutant.append(dish)
            print(f"    Mutant {i+1:02d}: {size} cells - kept")
        else:
            new_dish = copy.deepcopy(dish)
            new_dish.cells = random.sample(new_dish.cells, threshold_size)
            normalized_mutant.append(new_dish)
            print(f"    Mutant {i+1:02d}: {size} → {threshold_size} cells")

    # Process control1
    normalized_control1 = []
    for i, dish in enumerate(control1_dishes):
        size = len(dish.cells)
        if size < threshold_size:
            print(f"    Control1 {i+1:02d}: {size} cells - EXCLUDED")
        elif size == threshold_size:
            normalized_control1.append(dish)
            print(f"    Control1 {i+1:02d}: {size} cells - kept")
        else:
            new_dish = copy.deepcopy(dish)
            new_dish.cells = random.sample(new_dish.cells, threshold_size)
            normalized_control1.append(new_dish)
            print(f"    Control1 {i+1:02d}: {size} → {threshold_size} cells")

    total_kept = len(normalized_mutant) + len(normalized_control1)
    total_orig = len(mutant_dishes) + len(control1_dishes)
    print(f"  Retained: {total_kept}/{total_orig} ({100*total_kept/total_orig:.1f}%)")

    return normalized_mutant, normalized_control1, threshold_size


def save_petri_dish(petri: PetriDish, filepath: str, compress: bool = True) -> None:
    """Save PetriDish to file."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    # Ensure correct extension
    if compress and not filepath.endswith('.gz'):
        filepath = filepath + '.gz' if filepath.endswith('.json') else filepath + '.json.gz'
    elif not compress and filepath.endswith('.gz'):
        filepath = filepath[:-3]

    # Use Phase 1's save method
    abs_filepath = os.path.abspath(filepath)
    actual_path = petri.save_history(filename=abs_filepath, directory="", compress=compress)

    # Add Phase 2 metadata and individual_final
    if hasattr(petri, 'metadata') and petri.metadata:
        with smart_open(actual_path, 'r') as f:
            data = json.load(f)

        data['config']['phase2_metadata'] = petri.metadata
        data['individual_final'] = {
            'cells': [cell.to_dict() for cell in petri.cells],
            'gene_jsd': petri.calculate_gene_jsd()
        }

        with smart_open(actual_path, 'w') as f:
            if compress:
                json.dump(data, f, separators=(',', ':'))
            else:
                json.dump(data, f, indent=2)


# ============================================================================
# MAIN PIPELINE CLASS
# ============================================================================

class Phase2Pipeline:
    """Consolidated Phase 2 pipeline."""

    def __init__(self, args):
        self.args = args
        self.simulation_path = args.simulation
        self.first_year = args.first_snapshot
        self.second_year = args.second_year
        self.n_quantiles = args.n_quantiles
        self.cells_per_quantile = args.cells_per_quantile
        self.growth_phase = args.individual_growth_phase
        self.mix_ratio = args.mix_ratio
        self.seed = args.seed
        self.compress = args.compress

        # Set random seeds globally
        random.seed(self.seed)
        np.random.seed(self.seed)

        # Parse simulation path and generate output directory
        self.output_dir = self._generate_output_dir()

        # Create directories
        self.snapshots_dir = os.path.join(self.output_dir, "snapshots")
        self.individuals_dir = os.path.join(self.output_dir, "individuals")
        os.makedirs(self.snapshots_dir, exist_ok=True)
        os.makedirs(os.path.join(self.individuals_dir, "mutant"), exist_ok=True)
        os.makedirs(os.path.join(self.individuals_dir, "control1"), exist_ok=True)
        os.makedirs(os.path.join(self.individuals_dir, "control2"), exist_ok=True)

    def _generate_output_dir(self) -> str:
        """Generate output directory path."""
        # Parse phase1 simulation path
        sim_dir = os.path.dirname(self.simulation_path)

        # Generate phase2 subdirectory name
        timestamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        subdir_name = (f"snap{self.first_year}to{self.second_year}-"
                      f"growth{self.growth_phase}-"
                      f"quant{self.n_quantiles}x{self.cells_per_quantile}-"
                      f"mix{self.mix_ratio}-"
                      f"seed{self.seed}-{timestamp}")

        return os.path.join(sim_dir, subdir_name)

    def extract_snapshots(self) -> Tuple[List[Cell], List[Cell]]:
        """Extract snapshots from simulation."""
        print("\n" + "="*60)
        print("STAGE 1-2: Extract Snapshots")
        print("="*60)

        # Load simulation
        print(f"Loading simulation: {self.simulation_path}")
        with smart_open(self.simulation_path, 'r') as f:
            sim_data = json.load(f)

        # Get gene_rate_groups from config
        config = sim_data.get('config', sim_data.get('parameters', {}))
        gene_rate_groups = config.get('gene_rate_groups')
        if not gene_rate_groups:
            rate = config.get('rate')
            if rate:
                n = config.get('n', 1000)
                gene_size = config.get('gene_size', GENE_SIZE)
                gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
            else:
                raise ValueError("Simulation has no rate configuration!")

        # Save metadata
        metadata = {
            'gene_rate_groups': gene_rate_groups,
            'n_sites': config.get('n', 1000),
            'gene_size': config.get('gene_size', GENE_SIZE),
            'first_snapshot_year': self.first_year,
            'second_snapshot_year': self.second_year,
            'source_simulation': self.simulation_path
        }
        with open(os.path.join(self.snapshots_dir, 'metadata.json'), 'w') as f:
            json.dump(metadata, f, indent=2)

        # Extract snapshots
        history = sim_data['history']
        ext = ".json.gz" if self.compress else ".json"

        for year in [self.first_year, self.second_year]:
            year_str = str(year)
            if year_str not in history:
                available = sorted([int(y) for y in history.keys()])
                raise ValueError(f"Year {year} not found. Available: {available}")

            snapshot_path = os.path.join(self.snapshots_dir, f"year{year}_snapshot{ext}")
            year_data = history[year_str]
            snapshot_with_key = {year_str: year_data}

            with smart_open(snapshot_path, 'w') as f:
                json.dump(snapshot_with_key, f, indent=2)

            print(f"  Extracted year {year}: {len(year_data['cells'])} cells")

        # Load and return cells
        snapshot1_path = os.path.join(self.snapshots_dir, f"year{self.first_year}_snapshot{ext}")
        snapshot2_path = os.path.join(self.snapshots_dir, f"year{self.second_year}_snapshot{ext}")

        snapshot1_cells = load_snapshot_cells(snapshot1_path)
        snapshot2_cells = load_snapshot_cells(snapshot2_path)

        print(f"\n  Loaded {len(snapshot1_cells)} cells from year {self.first_year}")
        print(f"  Loaded {len(snapshot2_cells)} cells from year {self.second_year}")

        return snapshot1_cells, snapshot2_cells

    def create_individuals(self, snapshot1_cells: List[Cell]) -> Tuple[List[PetriDish], List[PetriDish]]:
        """Create initial individuals (mutant and control1)."""
        print("\n" + "="*60)
        print("STAGE 3: Create Initial Individuals")
        print("="*60)

        expected_count = self.n_quantiles * self.cells_per_quantile

        # Create mutant individuals (quantile sampling)
        print(f"\n  Creating mutant individuals (quantile sampling)...")
        mutant_dishes = []
        sampled = sample_by_quantiles(snapshot1_cells, self.n_quantiles,
                                     self.cells_per_quantile, self.seed)

        for i, (cell, quantile) in enumerate(sampled):
            petri = PetriDish.from_cells(
                [cell],
                growth_phase=self.growth_phase,
                metadata={
                    'individual_id': i + 1,
                    'individual_type': 'mutant',
                    'source_quantile': quantile,
                    'n_quantiles': self.n_quantiles,
                    'initial_year': self.first_year,
                    'growth_phase': self.growth_phase
                }
            )
            mutant_dishes.append(petri)

        print(f"    Created {len(mutant_dishes)} mutant individuals")

        # Create control1 individuals (uniform sampling)
        print(f"\n  Creating control1 individuals (uniform sampling)...")
        control1_dishes = []
        sampled_cells = sample_uniform(snapshot1_cells, expected_count, self.seed + 1000)

        for i, cell in enumerate(sampled_cells):
            petri = PetriDish.from_cells(
                [cell],
                growth_phase=self.growth_phase,
                metadata={
                    'individual_id': i + 1,
                    'individual_type': 'control1',
                    'initial_year': self.first_year,
                    'growth_phase': self.growth_phase
                }
            )
            control1_dishes.append(petri)

        print(f"    Created {len(control1_dishes)} control1 individuals")

        return mutant_dishes, control1_dishes

    def grow_populations(self, mutant_dishes: List[PetriDish],
                        control1_dishes: List[PetriDish]) -> None:
        """Grow individuals over time."""
        print("\n" + "="*60)
        print(f"STAGE 4: Grow Individuals ({self.second_year - self.first_year} years)")
        print("="*60)

        timeline_duration = self.second_year - self.first_year

        # Grow mutant
        print("\n  Growing mutant individuals...")
        for i, petri in enumerate(mutant_dishes):
            print(f"    Individual {i+1}/{len(mutant_dishes)}")
            petri.grow_exponentially(min(self.growth_phase, timeline_duration), verbose=False)
            if timeline_duration > self.growth_phase:
                homeostasis_years = timeline_duration - self.growth_phase
                petri.maintain_homeostasis(homeostasis_years, verbose=False)

        # Grow control1
        print("\n  Growing control1 individuals...")
        for i, petri in enumerate(control1_dishes):
            print(f"    Individual {i+1}/{len(control1_dishes)}")
            petri.grow_exponentially(min(self.growth_phase, timeline_duration), verbose=False)
            if timeline_duration > self.growth_phase:
                homeostasis_years = timeline_duration - self.growth_phase
                petri.maintain_homeostasis(homeostasis_years, verbose=False)

        print(f"\n  Growth complete")

    def normalize_and_mix(self, mutant_dishes: List[PetriDish],
                         control1_dishes: List[PetriDish],
                         snapshot2_cells: List[Cell]) -> Tuple[List[PetriDish], List[PetriDish], List[Cell], int]:
        """Normalize populations and mix with snapshot."""
        print("\n" + "="*60)
        print("STAGE 5: Normalize and Mix Populations")
        print("="*60)

        # Normalize
        print("\n  Applying size normalization (median - 0.5σ)...")
        mutant_dishes, control1_dishes, normalized_size = normalize_populations(
            mutant_dishes, control1_dishes, self.seed + 5000
        )

        # Update metadata with individual IDs after normalization
        for i, dish in enumerate(mutant_dishes, 1):
            dish.metadata['individual_id'] = i
            dish.metadata['normalized_size'] = normalized_size

        for i, dish in enumerate(control1_dishes, 1):
            dish.metadata['individual_id'] = i
            dish.metadata['normalized_size'] = normalized_size

        # Create common mixing pool
        print(f"\n  Creating common mixing pool...")
        mix_ratio_fraction = self.mix_ratio / 100

        if mix_ratio_fraction >= 1.0:
            target_total = normalized_size
            n_snapshot_cells = normalized_size
        else:
            target_total = int(normalized_size / (1 - mix_ratio_fraction))
            n_snapshot_cells = int(target_total * mix_ratio_fraction)

        print(f"    Normalized size: {normalized_size}")
        print(f"    Target total: {target_total}")
        print(f"    Snapshot cells needed: {n_snapshot_cells}")

        if n_snapshot_cells > len(snapshot2_cells):
            raise ValueError(f"Insufficient snapshot cells: need {n_snapshot_cells}, have {len(snapshot2_cells)}")

        random.seed(self.seed + 1000)
        indices = random.sample(range(len(snapshot2_cells)), n_snapshot_cells)
        common_pool = [copy.deepcopy(snapshot2_cells[idx]) for idx in indices]

        print(f"    Created pool of {len(common_pool)} cells")

        # Save common pool
        ext = '.json.gz' if self.compress else '.json'
        pool_path = os.path.join(self.individuals_dir, f'common_pool{ext}')
        pool_data = [cell.to_dict() for cell in common_pool]
        with smart_open(pool_path, 'w') as f:
            json.dump(pool_data, f, indent=2)

        # Mix all individuals
        print(f"\n  Mixing individuals with common pool...")

        for i, petri in enumerate(mutant_dishes):
            if mix_ratio_fraction >= 1.0:
                size = len(petri.cells)
                petri.cells = [copy.deepcopy(cell) for cell in common_pool[:size]]
                random.shuffle(petri.cells)
            else:
                added = [copy.deepcopy(cell) for cell in common_pool]
                petri.cells.extend(added)
                random.shuffle(petri.cells)

            petri.metadata['mix_ratio'] = self.mix_ratio
            print(f"    Mutant {i+1:02d}: {len(petri.cells)} cells")

        for i, petri in enumerate(control1_dishes):
            if mix_ratio_fraction >= 1.0:
                size = len(petri.cells)
                petri.cells = [copy.deepcopy(cell) for cell in common_pool[:size]]
                random.shuffle(petri.cells)
            else:
                added = [copy.deepcopy(cell) for cell in common_pool]
                petri.cells.extend(added)
                random.shuffle(petri.cells)

            petri.metadata['mix_ratio'] = self.mix_ratio
            print(f"    Control1 {i+1:02d}: {len(petri.cells)} cells")

        # Save mixing metadata
        mixing_metadata = {
            'mix_ratio': self.mix_ratio,
            'normalized_size': normalized_size
        }
        with open(os.path.join(self.individuals_dir, 'mixing_metadata.json'), 'w') as f:
            json.dump(mixing_metadata, f, indent=2)

        return mutant_dishes, control1_dishes, common_pool, normalized_size

    def create_control2(self, snapshot2_cells: List[Cell], common_pool: List[Cell],
                       target_size: int, num_individuals: int) -> List[PetriDish]:
        """Create control2 individuals."""
        print("\n" + "="*60)
        print("STAGE 6: Create Control2")
        print("="*60)

        print(f"\n  Creating {num_individuals} control2 individuals...")

        control2_dishes = []
        for i in range(num_individuals):
            # Start with common pool
            combined_cells = [copy.deepcopy(cell) for cell in common_pool]

            # Add additional cells if needed
            n_additional = target_size - len(common_pool)
            if n_additional > 0:
                random.seed(self.seed + 300 + i)
                additional_indices = random.sample(range(len(snapshot2_cells)),
                                                  min(n_additional, len(snapshot2_cells)))
                additional = [copy.deepcopy(snapshot2_cells[idx]) for idx in additional_indices]
                combined_cells.extend(additional)
            elif n_additional < 0:
                combined_cells = combined_cells[:target_size]

            # Create PetriDish
            petri = PetriDish(
                cells=combined_cells,
                gene_rate_groups=combined_cells[0].gene_rate_groups,
                n=combined_cells[0].n,
                gene_size=combined_cells[0].gene_size,
                growth_phase=None
            )

            petri.metadata = {
                'individual_id': i + 1,
                'individual_type': 'control2',
                'mix_ratio': self.mix_ratio,
                'normalized_size': target_size
            }

            control2_dishes.append(petri)
            print(f"    Control2 {i+1:02d}: {len(petri.cells)} cells")

        return control2_dishes

    def save_outputs(self, mutant_dishes: List[PetriDish],
                    control1_dishes: List[PetriDish],
                    control2_dishes: List[PetriDish]) -> None:
        """Save all individuals to disk."""
        print("\n" + "="*60)
        print("Saving Outputs")
        print("="*60)

        # Save mutant
        mutant_dir = os.path.join(self.individuals_dir, "mutant")
        for petri in mutant_dishes:
            individual_id = petri.metadata['individual_id']
            filepath = os.path.join(mutant_dir, f"individual_{individual_id:02d}.json")
            save_petri_dish(petri, filepath, self.compress)
        print(f"  Saved {len(mutant_dishes)} mutant individuals")

        # Save control1
        control1_dir = os.path.join(self.individuals_dir, "control1")
        for petri in control1_dishes:
            individual_id = petri.metadata['individual_id']
            filepath = os.path.join(control1_dir, f"individual_{individual_id:02d}.json")
            save_petri_dish(petri, filepath, self.compress)
        print(f"  Saved {len(control1_dishes)} control1 individuals")

        # Save control2 (special format)
        control2_dir = os.path.join(self.individuals_dir, "control2")
        ext = ".json.gz" if self.compress else ".json"

        for petri in control2_dishes:
            individual_id = petri.metadata['individual_id']
            filepath = os.path.join(control2_dir, f"individual_{individual_id:02d}{ext}")

            # Create simplified structure for control2
            control2_data = {
                "config": {
                    "gene_rate_groups": petri.cells[0].gene_rate_groups,
                    "n": petri.cells[0].n,
                    "gene_size": petri.cells[0].gene_size,
                    "growth_phase": None,
                    "years": self.second_year,
                    "seed": None,
                    "phase2_metadata": petri.metadata
                },
                "individual_final": {
                    "cells": [cell.to_dict() for cell in petri.cells],
                    "gene_jsd": petri.calculate_gene_jsd()
                }
            }

            with smart_open(filepath, 'w') as f:
                if self.compress:
                    json.dump(control2_data, f, separators=(',', ':'))
                else:
                    json.dump(control2_data, f, indent=2)

        print(f"  Saved {len(control2_dishes)} control2 individuals")

    def run(self) -> None:
        """Execute complete pipeline."""
        start_time = time.time()

        print("="*80)
        print("PHASE 2 PIPELINE")
        print("="*80)
        print(f"Simulation: {self.simulation_path}")
        print(f"Output: {self.output_dir}")
        print(f"Snapshots: years {self.first_year}, {self.second_year}")
        print(f"Quantiles: {self.n_quantiles} × {self.cells_per_quantile} cells")
        print(f"Growth: {self.growth_phase} years")
        print(f"Mix ratio: {self.mix_ratio}%")
        print(f"Seed: {self.seed}")
        print("="*80)

        # Execute stages
        snapshot1_cells, snapshot2_cells = self.extract_snapshots()
        mutant_dishes, control1_dishes = self.create_individuals(snapshot1_cells)
        self.grow_populations(mutant_dishes, control1_dishes)
        mutant_dishes, control1_dishes, common_pool, normalized_size = \
            self.normalize_and_mix(mutant_dishes, control1_dishes, snapshot2_cells)

        # Determine control2 count
        num_control2 = (len(mutant_dishes) + len(control1_dishes)) // 2
        control2_dishes = self.create_control2(snapshot2_cells, common_pool,
                                               normalized_size, num_control2)

        self.save_outputs(mutant_dishes, control1_dishes, control2_dishes)

        # Summary
        elapsed = time.time() - start_time
        print("\n" + "="*80)
        print("PIPELINE COMPLETE")
        print("="*80)
        print(f"Total time: {elapsed:.2f} seconds")
        print(f"Output directory: {self.output_dir}")
        print(f"Created:")
        print(f"  - {len(mutant_dishes)} mutant individuals")
        print(f"  - {len(control1_dishes)} control1 individuals")
        print(f"  - {len(control2_dishes)} control2 individuals")
        print("\nNext step: Run phase3 for analysis")
        print(f"  cd ../phase3")
        print(f"  python run_pipeline.py --phase2-dir {self.output_dir}")
        print("="*80)


# ============================================================================
# MAIN ENTRY POINT
# ============================================================================

def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """Load configuration from YAML file."""
    config = {}

    if not HAS_YAML:
        return config

    # Load default config
    default_path = os.path.join(os.path.dirname(__file__), "config_default.yaml")
    if os.path.exists(default_path):
        with open(default_path, 'r') as f:
            config.update(yaml.safe_load(f) or {})

    # Load custom config
    if config_path and os.path.exists(config_path):
        with open(config_path, 'r') as f:
            config.update(yaml.safe_load(f) or {})

    return config


def handle_glob_patterns(simulation_path: str) -> str:
    """Handle glob patterns in simulation path."""
    if '*' not in simulation_path:
        return simulation_path

    matching_files = glob.glob(simulation_path)
    if not matching_files:
        print(f"Error: No files matching pattern: {simulation_path}")
        sys.exit(1)

    if len(matching_files) == 1:
        return matching_files[0]

    # Multiple files - let user choose
    print(f"\nFound {len(matching_files)} simulations:")
    for i, f in enumerate(sorted(matching_files), 1):
        display = f.split("phase1/data/")[1] if "phase1/data/" in f else os.path.basename(f)
        print(f"  [{i}] {display}")

    while True:
        try:
            sel = input("\nSelect simulation (number or 'q' to quit): ").strip()
            if sel.lower() == 'q':
                sys.exit(0)
            idx = int(sel) - 1
            if 0 <= idx < len(matching_files):
                return matching_files[idx]
            print(f"Invalid selection. Enter 1-{len(matching_files)}")
        except (ValueError, KeyboardInterrupt):
            sys.exit(0)


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 Pipeline - Data Generation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Config file
    parser.add_argument("--config", type=str, help="YAML configuration file")

    # Required
    parser.add_argument("--simulation", type=str, required=True,
                       help="Path to phase1 simulation file")

    # Snapshot parameters
    parser.add_argument("--first-snapshot", type=int, default=None,
                       help="First snapshot year")
    parser.add_argument("--second-snapshot", type=int, default=None,
                       help="Second snapshot year (renamed from --second-year)")

    # Sampling parameters
    parser.add_argument("--n-quantiles", type=int, default=None,
                       help="Number of quantiles")
    parser.add_argument("--cells-per-quantile", type=int, default=None,
                       help="Cells per quantile")

    # Growth parameters
    parser.add_argument("--individual-growth-phase", type=int, default=None,
                       help="Years of exponential growth")
    parser.add_argument("--mix-ratio", type=int, default=None,
                       help="Percentage from second snapshot")

    # Other
    parser.add_argument("--seed", type=int, default=None,
                       help="Random seed")

    compress_group = parser.add_mutually_exclusive_group()
    compress_group.add_argument("--compress", action='store_true', dest='compress')
    compress_group.add_argument("--no-compress", action='store_false', dest='compress')
    parser.set_defaults(compress=None)

    args = parser.parse_args()

    # Load config and apply defaults
    config = load_config(args.config)

    if args.first_snapshot is None:
        args.first_snapshot = config.get('first_snapshot', 30)
    if args.second_snapshot is None:
        args.second_snapshot = config.get('second_snapshot', 50)
    if args.n_quantiles is None:
        args.n_quantiles = config.get('n_quantiles', 4)
    if args.cells_per_quantile is None:
        args.cells_per_quantile = config.get('cells_per_quantile', 2)
    if args.individual_growth_phase is None:
        args.individual_growth_phase = config.get('individual_growth_phase', 6)
    if args.mix_ratio is None:
        args.mix_ratio = config.get('mix_ratio', 70)
    if args.seed is None:
        args.seed = config.get('seed', 42)
    if args.compress is None:
        args.compress = config.get('compress', False)

    # Add second_year alias for internal use
    args.second_year = args.second_snapshot

    # Handle glob patterns
    args.simulation = handle_glob_patterns(args.simulation)

    # Validate simulation file
    if not os.path.exists(args.simulation):
        print(f"Error: Simulation file not found: {args.simulation}")
        sys.exit(1)

    # Run pipeline
    pipeline = Phase2Pipeline(args)
    pipeline.run()


if __name__ == "__main__":
    main()
