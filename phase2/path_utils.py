#!/usr/bin/env python3
"""
Utilities for parsing and generating hierarchical paths for phase1 and phase2.
"""

import os
import re
import hashlib
from typing import Dict, Optional


def parse_step1_simulation_path(filepath: str) -> Optional[Dict]:
    """
    Parse parameters from phase1 simulation path.
    
    Handles both old and new formats:
    - Old: simulation_rate_0.005000_g13_n100_m957_t100_seed42.json.gz
    - New: data/rate_0.00500/grow13-sites100-years100-seed42-a3f2/simulation.json.gz
    
    Returns:
        Dictionary with rate, growth_phase, n_sites, sim_years, sim_seed
    """
    # Try new hierarchical format first
    if 'simulation.json.gz' in filepath:
        # Extract from directory structure
        # Pattern: .../rate_X.XXXXX/growG-sitesN-yearsT-seedS-HASH/simulation.json.gz
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


def generate_step23_output_dir(args, sim_params: Dict) -> str:
    """
    Generate hierarchical output directory for phase2.
    
    Structure:
    data/rate_0.00500-grow13-sites100-years100/snap50-quant10x3-grow10-mix80-seed42-HASH/
    
    Args:
        args: Command line arguments with pipeline parameters
        sim_params: Dictionary from parse_step1_simulation_path
        
    Returns:
        Full path to output directory
    """
    # Level 1: Rate and source simulation params
    level1 = (f"rate_{args.rate:.5f}-"
              f"grow{sim_params['growth_phase']}-"
              f"sites{sim_params['n_sites']}-"
              f"years{sim_params['sim_years']}")
    
    # Level 2: Pipeline params in logical flow order
    params_str = (f"snap{args.snapshot_year}-"
                  f"quant{args.n_quantiles}x{args.cells_per_quantile}-"
                  f"grow{args.growth_years}-"
                  f"mix{args.mix_ratio}-"
                  f"seed{args.seed}")
    
    # Add 4-char hash for uniqueness
    # Include all parameters in hash to ensure uniqueness
    hash_input = (f"{args.rate:.5f}-{sim_params['growth_phase']}-"
                  f"{sim_params['n_sites']}-{sim_params['sim_years']}-"
                  f"{args.snapshot_year}-{args.n_quantiles}-"
                  f"{args.cells_per_quantile}-{args.growth_years}-"
                  f"{args.mix_ratio}-{args.seed}")
    hash_str = hashlib.md5(hash_input.encode()).hexdigest()[:4]
    level2 = f"{params_str}-{hash_str}"
    
    return os.path.join(args.output_dir, level1, level2)


def compute_hash(params_str: str, length: int = 4) -> str:
    """
    Compute a hash of the given parameters string.
    
    Args:
        params_str: String containing all parameters
        length: Number of hex characters to return (default 4)
        
    Returns:
        Hash string of specified length
    """
    return hashlib.md5(params_str.encode()).hexdigest()[:length]