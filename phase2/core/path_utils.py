#!/usr/bin/env python3
"""
Utilities for parsing and generating hierarchical paths for phase1 and phase2.
"""

import os
import re
import numpy as np
from datetime import datetime
from typing import Dict, Optional


def parse_step1_simulation_path(filepath: str) -> Optional[Dict]:
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


def generate_step23_output_dir(args, sim_params: Dict) -> str:
    """
    Generate hierarchical output directory for phase2.
    
    Structure:
    data/rate_0.00500-grow13-sites100-years100/snap50to60-growth7-quant10x3-mix80-seed42-HASH/
    
    Args:
        args: Command line arguments with pipeline parameters
        sim_params: Dictionary from parse_step1_simulation_path
        
    Returns:
        Full path to output directory
    """
    # Level 1: Rate and source simulation params
    # Generate rate string based on configuration
    if hasattr(args, 'gene_rate_groups') and args.gene_rate_groups:
        # Parse and format gene rate groups
        groups = []
        for group in args.gene_rate_groups.split(','):
            n, rate = group.split(':')
            groups.append(f"{n}x{float(rate):.5f}")
        rate_str = "gene_rates_" + "_".join(groups)
        # Truncate if too long to avoid filesystem issues
        if len(rate_str) > 50:
            rate_str = rate_str[:47] + "..."
    elif hasattr(args, 'rate') and args.rate:
        rate_str = f"rate_{args.rate:.5f}"
    else:
        # No rate specified - must have been inferred from simulation
        # Extract from simulation path
        import re
        match = re.search(r'(gene_rates_[^/]+|rate_[\d.]+)', args.simulation)
        if match:
            rate_str = match.group(1)
        else:
            rate_str = "rate_inferred"
    
    level1 = (f"{rate_str}-"
              f"grow{sim_params['growth_phase']}-"
              f"sites{sim_params['n_sites']}-"
              f"years{sim_params['sim_years']}")
    
    # Level 2: Pipeline params in logical flow order
    # Build mix suffix: 'u' for uniform, 'n' for normalized, 'un' for both
    mix_suffix = ""
    if hasattr(args, 'uniform_mixing') and args.uniform_mixing:
        mix_suffix += "u"
    if hasattr(args, 'normalize_size') and args.normalize_size:
        mix_suffix += "n"
    
    params_str = (f"snap{args.first_snapshot}to{args.second_snapshot}-"
                  f"growth{args.individual_growth_phase}-"
                  f"quant{args.n_quantiles}x{args.cells_per_quantile}-"
                  f"mix{args.mix_ratio}{mix_suffix}-"
                  f"seed{args.seed}")
    
    # Add timestamp for uniqueness (YYYYMMDDHHMMSS format)
    timestamp = datetime.now().strftime("%Y%m%d%H%M%S")
    level2 = f"{params_str}-{timestamp}"
    
    return os.path.join(args.output_dir, level1, level2)


