"""
Helper functions for creating and managing individuals in the pipeline.
Reduces code duplication between mutant and control1 processing.
"""

import os
from typing import List, Dict, Optional, Tuple
from cell import PetriDish, Cell
from .pipeline_utils import save_petri_dish, grow_petri_for_years, mix_petri_with_snapshot


def create_individual(
    cell: Cell,
    individual_type: str,
    individual_id: int,
    growth_phase: int,
    additional_metadata: Optional[Dict] = None
) -> PetriDish:
    """
    Create an individual PetriDish from a cell (in-memory only).
    
    Args:
        cell: The founding cell for this individual
        individual_type: 'mutant' or 'control1'
        individual_id: Unique ID for this individual (matches filename)
        growth_phase: Growth phase for the PetriDish
        additional_metadata: Extra metadata to include
        
    Returns:
        The created PetriDish object
    """
    # Build metadata
    metadata = {
        'individual_id': individual_id,
        'individual_type': individual_type,
        'initial_year': 50  # Standard for phase2
    }
    if additional_metadata:
        metadata.update(additional_metadata)
    
    # Create PetriDish using new cleaner method
    petri = PetriDish.from_cells(
        cell,
        growth_phase=growth_phase,
        calculate_cell_jsds=True,
        metadata=metadata
    )
    
    return petri


def grow_individual(
    petri: PetriDish,
    years: int,
    growth_phase: int,
    verbose: bool = True,
    start_year: Optional[int] = None
) -> None:
    """
    Grow an individual for specified years (in-memory only).
    
    Args:
        petri: PetriDish to grow
        years: Number of years to grow
        growth_phase: Years of exponential growth
        verbose: Print progress
        start_year: Starting year for history tracking
    """
    # Grow with history tracking
    grow_petri_for_years(
        petri, years, growth_phase,
        verbose=verbose,
        track_history=True,
        start_year=start_year
    )


def mix_individual(
    petri: PetriDish,
    snapshot_cells: List[Cell],
    mix_ratio: float,
    seed: int
) -> int:
    """
    Mix an individual with snapshot cells (in-memory only).
    
    Args:
        petri: PetriDish to mix
        snapshot_cells: Cells from second snapshot
        mix_ratio: Fraction from snapshot (0.0 to 1.0)
        seed: Random seed for mixing
        
    Returns:
        Total number of cells after mixing
    """
    # Mix with snapshot
    total_cells = mix_petri_with_snapshot(
        petri, snapshot_cells,
        mix_ratio=mix_ratio,
        seed=seed
    )
    
    # Update metadata
    if not hasattr(petri, 'metadata'):
        petri.metadata = {}
    petri.metadata.update({
        'mixed': True,
        'mix_ratio': int(mix_ratio * 100),
        'final_cells': total_cells
    })
    
    return total_cells


def process_batch_growth(
    dishes: List[PetriDish],
    batch_name: str,
    years: int,
    growth_phase: int,
    expected_population: int,
    start_year: int,
    verbose: bool = True
) -> None:
    """
    Process growth for a batch of individuals (mutant or control1) in memory.
    
    Args:
        dishes: List of PetriDish objects to grow
        batch_name: Name of batch ('mutant' or 'control1')
        years: Number of years to grow
        growth_phase: Years of exponential growth
        expected_population: Expected cell count after growth
        start_year: Starting year for history
        verbose: Print progress
    """
    print(f"\n  ✓ Growing {batch_name} individuals to ~{expected_population} cells...")
    
    for i, petri in enumerate(dishes):
        current_cells = len(petri.cells)
        individual_id = petri.metadata.get('individual_id', i + 1)
        
        if current_cells == 1:
            # Fresh individual, needs full growth
            print(f"    Individual {individual_id:02d}: {current_cells} → ~{expected_population} cells")
            grow_individual(
                petri, years, growth_phase,
                verbose, start_year
            )
        elif expected_population * 0.5 <= current_cells <= expected_population * 1.5:
            # Already grown (with homeostasis variation)
            print(f"    Individual {individual_id:02d}: Already at {current_cells} cells")
        else:
            # Already mixed or in unexpected state
            print(f"    Individual {individual_id:02d}: Already mixed ({current_cells} cells)")


def process_batch_mixing(
    dishes: List[PetriDish],
    batch_name: str,
    snapshot_cells: List[Cell],
    mix_ratio: int,
    expected_population: int,
    expected_final_cells: int,
    base_seed: int
) -> None:
    """
    Process mixing for a batch of individuals (mutant or control1) in memory.
    
    Args:
        dishes: List of PetriDish objects to mix
        batch_name: Name of batch ('mutant' or 'control1')
        snapshot_cells: Cells from second snapshot
        mix_ratio: Percentage from snapshot (0-100)
        expected_population: Expected cells before mixing
        expected_final_cells: Expected cells after mixing
        base_seed: Base random seed (will add index)
    """
    print(f"\n  ✓ Mixing {batch_name} individuals with year snapshot cells...")
    
    for i, petri in enumerate(dishes):
        individual_id = petri.metadata.get('individual_id', i + 1)
        current_cells = len(petri.cells)
        
        # Check if within expected range (homeostasis causes variation)
        if expected_population * 0.5 <= current_cells <= expected_population * 1.5:
            print(f"    Individual {individual_id:02d}: Mixing {current_cells} → {expected_final_cells} cells")
            
            mix_individual(
                petri, snapshot_cells,
                mix_ratio / 100.0,  # Convert percentage to fraction
                base_seed + i
            )
        elif current_cells > expected_population * 1.5:
            print(f"    Individual {individual_id:02d}: Already mixed ({current_cells} cells)")