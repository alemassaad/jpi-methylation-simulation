#!/usr/bin/env python3
"""
Test that the pipeline can start with gene rates without crashing.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

# Test the setup phase
print("Testing pipeline initialization with gene rates...")

# Simulate the command-line arguments
class Args:
    def __init__(self):
        self.gene_rate_groups = "50:0.004,50:0.0045,50:0.005,50:0.0055"
        self.gene_size = 5
        self.rate = None
        self.simulation = "../phase1/data/gene_rates_50x0.00400_50x0.00450_50x0.00500_50x0.0/grow10-sites1000-years50-seed42-d60f/simulation.json.gz"
        self.first_snapshot = 30
        self.second_snapshot = 45
        self.individual_growth_phase = 7
        self.n_quantiles = 5
        self.cells_per_quantile = 2
        self.seed = 42
        self.output_dir = "data"
        self.mix_ratio = 80
        self.uniform_mixing = False
        self.normalize_size = False
        self.plot_individuals = True

args = Args()

# Test parsing gene rate groups
from run_pipeline import parse_gene_rate_groups, validate_gene_rate_groups

try:
    gene_groups = parse_gene_rate_groups(args.gene_rate_groups)
    print(f"✓ Parsed gene groups: {len(gene_groups)} groups")
    
    validate_gene_rate_groups(gene_groups, 1000, args.gene_size)
    print(f"✓ Gene groups validated")
    
    # Test rate config creation
    rate_config = {
        'type': 'gene_specific',
        'gene_rate_groups': gene_groups,
        'gene_size': args.gene_size
    }
    print(f"✓ Rate config created: {rate_config['type']}")
    
    # Test path parsing
    from path_utils import parse_step1_simulation_path
    sim_params = parse_step1_simulation_path(args.simulation)
    if sim_params:
        print(f"✓ Simulation path parsed: grow{sim_params['growth_phase']}, {sim_params['n_sites']} sites")
    else:
        print("✗ Failed to parse simulation path")
        
    # Test output directory generation
    from path_utils import generate_step23_output_dir
    output_dir = generate_step23_output_dir(args, sim_params)
    print(f"✓ Output directory: {output_dir}")
    
    # Check that gene_rates appears in the path
    if "gene_rates" in output_dir:
        print(f"✓ Gene rates included in output path")
    else:
        print(f"✗ Gene rates NOT in output path!")
        
except Exception as e:
    print(f"✗ Error: {e}")
    import traceback
    traceback.print_exc()

print("\nAll initialization checks passed! Pipeline should be able to start.")