#!/usr/bin/env python3
"""
Utilities for parsing and generating hierarchical paths for phase1 and phase2.
"""

import os
import re
import numpy as np
from datetime import datetime
from typing import Dict, Optional


def parse_phase1_simulation_path(filepath: str) -> Optional[Dict]:
    """
    Parse parameters from phase1 simulation path.
    
    Handles multiple formats:
    - Old: simulation_rate_0.005000_g13_n100_m957_t100_seed42.json.gz
    - Newer: data/rate_0.00500/grow13-sites100-years100-seed42-a3f2/simulation.json.gz
    - Newest: data/gene_rates_200x0.00500/size8192-sites1000-genesize5-years100-seed42-YYYYMMDD-HHMMSS/simulation.json.gz
    
    Returns:
        Dictionary with rate, growth_phase, n_sites, sim_years, sim_seed
    """
    # Try new hierarchical format first (both compressed and uncompressed)
    if 'simulation.json' in filepath:
        # Extract from directory structure
        
        # Pattern for newest format with size instead of grow
        # .../gene_rates_.../size8192-sites1000-genesize5-years100-seed42-YYYYMMDD-HHMMSS/simulation.json.gz
        pattern = r"gene_rates_[^/]+/size(\d+)-sites(\d+)-genesize(\d+)-years(\d+)-(seed\d+|noseed)-"
        match = re.search(pattern, filepath)
        if match:
            size = int(match.group(1))
            growth_phase = int(np.log2(size)) if size > 0 else 0
            seed_str = match.group(5)
            seed = int(seed_str.replace('seed', '')) if seed_str != 'noseed' else None
            return {
                'rate': None,  # Gene-specific rates, not a single rate
                'growth_phase': growth_phase,
                'n_sites': int(match.group(2)),
                'gene_size': int(match.group(3)),
                'sim_years': int(match.group(4)),
                'sim_seed': seed
            }
        
        # Pattern for older format with grow
        # Pattern for uniform rate: .../rate_X.XXXXX/growG-sitesN-yearsT-seedS-HASH/simulation.json.gz
        pattern = r"rate_([\d.]+)/grow(\d+)-sites(\d+)-years(\d+)-(seed\d+|noseed)-"
        match = re.search(pattern, filepath)
        if match:
            seed_str = match.group(5)
            seed = int(seed_str.replace('seed', '')) if seed_str != 'noseed' else None
            return {
                'rate': float(match.group(1)),
                'growth_phase': int(match.group(2)),
                'n_sites': int(match.group(3)),
                'sim_years': int(match.group(4)),
                'sim_seed': seed
            }
        
        # Pattern for gene-specific rates with grow: .../gene_rates_.../growG-sitesN-yearsT-seedS-HASH/simulation.json.gz
        pattern = r"gene_rates_[^/]+/grow(\d+)-sites(\d+)-years(\d+)-(seed\d+|noseed)-"
        match = re.search(pattern, filepath)
        if match:
            seed_str = match.group(4)
            seed = int(seed_str.replace('seed', '')) if seed_str != 'noseed' else None
            return {
                'rate': None,  # Gene-specific rates, not a single rate
                'growth_phase': int(match.group(1)),
                'n_sites': int(match.group(2)),
                'sim_years': int(match.group(3)),
                'sim_seed': seed
            }
    
    # Try old flat format
    basename = os.path.basename(filepath)
    # Old format with m: simulation_rate_X.XXXXXX_gG_mM_nN_tT_seedS.json.gz (m before n)
    pattern = r"simulation_rate_([\d.]+)_g(\d+)_m\d+_n(\d+)_t(\d+)(?:_seed(\d+)|_noseed)?"
    match = re.match(pattern, basename)
    if match:
        return {
            'rate': float(match.group(1)),
            'growth_phase': int(match.group(2)),
            'n_sites': int(match.group(3)),
            'sim_years': int(match.group(4)),
            'sim_seed': int(match.group(5)) if match.group(5) else None
        }
    
    # Try new flat format without m: simulation_rate_X.XXXXXX_gG_nN_tT_seedS.json.gz
    pattern = r"simulation_rate_([\d.]+)_g(\d+)_n(\d+)_t(\d+)(?:_seed(\d+)|_noseed)?"
    match = re.match(pattern, basename)
    if match:
        return {
            'rate': float(match.group(1)),
            'growth_phase': int(match.group(2)),
            'n_sites': int(match.group(3)),
            'sim_years': int(match.group(4)),
            'sim_seed': int(match.group(5)) if match.group(5) else None
        }
    
    return None


def generate_phase2_output_dir(args, sim_params: Dict) -> str:
    """
    Generate output directory for phase2 under the phase1 simulation directory.
    
    Structure:
    {phase1_simulation_dir}/snap30to50-growth7-quant10x3-mix80-seed42-TIMESTAMP/
    
    Args:
        args: Command line arguments with pipeline parameters (must include simulation path)
        sim_params: Dictionary from parse_phase1_simulation_path (not used anymore but kept for compatibility)
        
    Returns:
        Full path to output directory under the phase1 simulation directory
    """
    # Extract the parent directory of the simulation file
    # This is where phase1 stored its output
    simulation_dir = os.path.dirname(os.path.abspath(args.simulation))
    
    # Build phase2-specific directory name with parameters
    params_str = (f"snap{args.first_snapshot}to{args.second_snapshot}-"
                  f"growth{args.individual_growth_phase}-"
                  f"quant{args.n_quantiles}x{args.cells_per_quantile}-"
                  f"mix{args.mix_ratio}-"
                  f"seed{args.seed}")
    
    # Add timestamp for uniqueness (YYYYMMDDHHMMSS format)
    # This allows multiple phase2 runs with same parameters
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    phase2_subdir = f"{params_str}-{timestamp}"
    
    # Return path under the simulation directory
    return os.path.join(simulation_dir, phase2_subdir)


