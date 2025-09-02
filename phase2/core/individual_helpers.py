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
    output_dir: str,
    compress: bool = True,
    additional_metadata: Optional[Dict] = None
) -> PetriDish:
    """
    Create an individual PetriDish from a cell and save it.
    
    Args:
        cell: The founding cell for this individual
        individual_type: 'mutant' or 'control1'
        individual_id: Unique ID for this individual (matches filename)
        growth_phase: Growth phase for the PetriDish
        output_dir: Directory to save the individual
        compress: Whether to compress the output file
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
    
    # Save immediately
    filepath = os.path.join(output_dir, f"individual_{individual_id:02d}.json")
    if compress:
        filepath += '.gz'
    save_petri_dish(petri, filepath, compress=compress)
    
    return petri


def grow_individual(
    petri: PetriDish,
    years: int,
    growth_phase: int,
    output_dir: str,
    verbose: bool = True,
    compress: bool = True,
    start_year: Optional[int] = None
) -> None:
    """
    Grow an individual for specified years and save.
    
    Args:
        petri: PetriDish to grow
        years: Number of years to grow
        growth_phase: Years of exponential growth
        output_dir: Directory to save the grown individual
        verbose: Print progress
        compress: Whether to compress the output file
        start_year: Starting year for history tracking
    """
    # Grow with history tracking
    grow_petri_for_years(
        petri, years, growth_phase,
        verbose=verbose,
        track_history=True,
        start_year=start_year
    )
    
    # Save with history
    individual_id = petri.metadata.get('individual_id', 1)
    filepath = os.path.join(output_dir, f"individual_{individual_id:02d}.json")
    if compress:
        filepath += '.gz'
    save_petri_dish(
        petri, filepath,
        include_cell_history=True,
        include_gene_metrics=True,
        compress=compress
    )


def mix_individual(
    petri: PetriDish,
    snapshot_cells: List[Cell],
    mix_ratio: float,
    output_dir: str,
    seed: int,
    compress: bool = True
) -> int:
    """
    Mix an individual with snapshot cells and save.
    
    Args:
        petri: PetriDish to mix
        snapshot_cells: Cells from second snapshot
        mix_ratio: Fraction from snapshot (0.0 to 1.0)
        output_dir: Directory to save the mixed individual
        seed: Random seed for mixing
        compress: Whether to compress the output file
        
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
    
    # Save with history
    individual_id = petri.metadata.get('individual_id', 1)
    filepath = os.path.join(output_dir, f"individual_{individual_id:02d}.json")
    if compress:
        filepath += '.gz'
    save_petri_dish(
        petri, filepath,
        include_cell_history=True,
        include_gene_metrics=True,
        compress=compress
    )
    
    return total_cells


def process_batch_growth(
    dishes: List[PetriDish],
    batch_name: str,
    output_dir: str,
    years: int,
    growth_phase: int,
    expected_population: int,
    start_year: int,
    compress: bool = True,
    verbose: bool = True
) -> None:
    """
    Process growth for a batch of individuals (mutant or control1).
    
    Args:
        dishes: List of PetriDish objects to grow
        batch_name: Name of batch ('mutant' or 'control1')
        output_dir: Directory to save grown individuals
        years: Number of years to grow
        growth_phase: Years of exponential growth
        expected_population: Expected cell count after growth
        start_year: Starting year for history
        compress: Whether to compress output files
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
                output_dir, verbose, compress, start_year
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
    output_dir: str,
    snapshot_cells: List[Cell],
    mix_ratio: int,
    expected_population: int,
    expected_final_cells: int,
    base_seed: int,
    compress: bool = True
) -> None:
    """
    Process mixing for a batch of individuals (mutant or control1).
    
    Args:
        dishes: List of PetriDish objects to mix
        batch_name: Name of batch ('mutant' or 'control1')
        output_dir: Directory to save mixed individuals
        snapshot_cells: Cells from second snapshot
        mix_ratio: Percentage from snapshot (0-100)
        expected_population: Expected cells before mixing
        expected_final_cells: Expected cells after mixing
        base_seed: Base random seed (will add index)
        compress: Whether to compress output files
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
                output_dir,
                base_seed + i,
                compress
            )
        elif current_cells > expected_population * 1.5:
            print(f"    Individual {individual_id:02d}: Already mixed ({current_cells} cells)")