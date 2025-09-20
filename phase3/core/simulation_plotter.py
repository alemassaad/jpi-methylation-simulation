"""
Simulation history plotting utilities.
Replaces phase1/plot_history.py functionality.
"""
import json
import gzip
import glob
import os
import sys
from typing import Optional, Dict, Any

# Add phase1 to path for PetriDish access
phase1_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1')
sys.path.insert(0, phase1_path)
from cell import PetriDish, Cell

# Import our plotting modules
from .petri_dish_plotter import PetriDishPlotter


def find_simulation_file(path_pattern: str) -> Optional[str]:
    """
    Find simulation file from pattern (supports wildcards).
    
    Args:
        path_pattern: Path pattern with optional wildcards
        
    Returns:
        Path to simulation file or None if not found
    """
    matches = glob.glob(path_pattern)
    if not matches:
        # Try in data directory
        data_pattern = os.path.join('data', path_pattern)
        matches = glob.glob(data_pattern)
    
    if not matches:
        # Try common patterns
        if not path_pattern.endswith('.json.gz'):
            gz_pattern = path_pattern + '/simulation.json.gz'
            matches = glob.glob(gz_pattern)
    
    if not matches:
        return None
    
    # If multiple matches, pick the first .json.gz file
    for match in matches:
        if match.endswith('simulation.json.gz'):
            return match
    
    # Otherwise return first match
    return matches[0]


def load_simulation_history(filepath: str) -> Dict[str, Any]:
    """
    Load history from phase1 simulation file.
    
    Args:
        filepath: Path to simulation.json or simulation.json.gz
        
    Returns:
        Dictionary containing simulation data
    """
    print(f"Loading history from {filepath}...")
    
    if filepath.endswith('.gz'):
        with gzip.open(filepath, 'rt') as f:
            return json.load(f)
    else:
        with open(filepath, 'r') as f:
            return json.load(f)


def create_petri_from_history(data: Dict[str, Any]) -> PetriDish:
    """
    Create a PetriDish object from simulation history data.
    
    Args:
        data: Dictionary from load_simulation_history
        
    Returns:
        PetriDish object with history loaded
    """
    params = data['parameters']
    history = data['history']
    
    # Extract parameters
    gene_rate_groups = params.get('gene_rate_groups')
    n = params.get('n', 1000)
    gene_size = params.get('gene_size', 5)
    growth_phase = params.get('growth_phase', 13)
    seed = params.get('seed')
    
    # Create PetriDish with proper parameters
    petri = PetriDish(
        gene_rate_groups=gene_rate_groups,
        n=n,
        gene_size=gene_size,
        growth_phase=growth_phase,
        seed=seed,
        calculate_cell_jsds=False
    )
    
    # Load cell history
    petri.cell_history = {}
    for year_str, year_data in history.items():
        if 'cells' in year_data:
            petri.cell_history[year_str] = year_data['cells']
        else:
            # Handle old format where history[year] is directly the cells list
            petri.cell_history[year_str] = year_data
    
    # Load gene JSD history if available
    if any('gene_jsd' in year_data for year_data in history.values() if isinstance(year_data, dict)):
        petri.gene_jsd_history = {}
        for year_str, year_data in history.items():
            if isinstance(year_data, dict) and 'gene_jsd' in year_data:
                petri.gene_jsd_history[year_str] = year_data['gene_jsd']
    
    # Load mean and median gene JSD histories if available
    if any('mean_gene_jsd' in year_data for year_data in history.values() if isinstance(year_data, dict)):
        petri.mean_gene_jsd_history = {}
        petri.median_gene_jsd_history = {}
        for year_str, year_data in history.items():
            if isinstance(year_data, dict):
                if 'mean_gene_jsd' in year_data:
                    petri.mean_gene_jsd_history[year_str] = year_data['mean_gene_jsd']
                if 'median_gene_jsd' in year_data:
                    petri.median_gene_jsd_history[year_str] = year_data['median_gene_jsd']
    
    # Set the year to the last year in history
    years = [int(y) for y in petri.cell_history.keys()]
    petri.year = max(years) if years else 0
    
    # Set attributes from parameters
    petri.n = n
    petri.gene_size = gene_size
    petri.gene_rate_groups = gene_rate_groups
    
    return petri


def plot_simulation_history(simulation_path: str, output_dir: str = None, 
                           title_prefix: str = None, plots_to_generate: list = None):
    """
    Main plotting function for phase1 simulations.
    
    Args:
        simulation_path: Path to simulation.json.gz file
        output_dir: Optional output directory for plots
        title_prefix: Optional prefix for plot titles
        plots_to_generate: List of plot types to generate (default: all)
                         Options: ['jsd', 'methylation', 'combined', 'gene_jsd']
    """
    # Find and load simulation
    filepath = find_simulation_file(simulation_path)
    if not filepath:
        raise FileNotFoundError(f"Could not find simulation file: {simulation_path}")
    
    # Load data
    data = load_simulation_history(filepath)
    
    # Create PetriDish from history
    petri = create_petri_from_history(data)
    
    # Create plotter
    plotter = PetriDishPlotter(petri)
    
    # Determine output directory
    if output_dir is None:
        # Use the same directory as the simulation file
        output_dir = os.path.dirname(filepath)
    
    # Create output directory if needed
    os.makedirs(output_dir, exist_ok=True)
    
    # Base path for output files
    base_name = os.path.splitext(os.path.basename(filepath))[0]
    if base_name == 'simulation':
        base_name = 'plot'
    base_path = os.path.join(output_dir, base_name)
    
    # Default to all plots
    if plots_to_generate is None:
        plots_to_generate = ['jsd', 'methylation', 'combined']
        if hasattr(petri, 'gene_jsd_history') and petri.gene_jsd_history:
            plots_to_generate.append('gene_jsd')
    
    # Generate requested plots
    title = title_prefix or "Simulation Results"
    
    print(f"\nGenerating plots in: {output_dir}")
    
    if 'jsd' in plots_to_generate:
        print("  - Creating JSD plot...")
        plotter.plot_jsd(f"{title} - Cell JSD", f"{base_path}_jsd.png")
    
    if 'methylation' in plots_to_generate:
        print("  - Creating methylation plot...")
        plotter.plot_cell_methylation_proportion(f"{title} - Methylation", 
                                                f"{base_path}_methylation.png")
    
    if 'combined' in plots_to_generate:
        print("  - Creating combined plot...")
        plotter.plot_combined(title, f"{base_path}_combined.png")
    
    if 'gene_jsd' in plots_to_generate and hasattr(petri, 'gene_jsd_history') and petri.gene_jsd_history:
        print("  - Creating gene JSD plot...")
        plotter.plot_gene_jsd(f"{title} - Gene JSD", f"{base_path}_gene_jsd.png")
    
    print("\nPlots generated successfully!")
    return plotter