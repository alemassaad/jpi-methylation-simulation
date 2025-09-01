#!/usr/bin/env python3
"""
Test with the actual path the user is trying to use.
"""

import sys
import os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from path_utils import parse_step1_simulation_path

# Test with the actual path
test_path = "../phase1/data/gene_rates_50x0.00400_50x0.00450_50x0.00500_50x0.0/grow10-sites1000-years50-seed42-d60f/simulation.json.gz"

print(f"Testing path: {test_path}")
result = parse_step1_simulation_path(test_path)

if result:
    print("✓ Path parsed successfully!")
    print(f"  Growth phase: {result['growth_phase']}")
    print(f"  Sites: {result['n_sites']}")
    print(f"  Years: {result['sim_years']}")
    print(f"  Seed: {result['sim_seed']}")
    print(f"  Rate: {result['rate']}")
else:
    print("✗ Failed to parse path")
    print("\nDEBUGGING:")
    
    # Check what pattern would match
    import re
    
    # Check uniform rate pattern
    pattern1 = r"rate_([\d.]+)/grow(\d+)-sites(\d+)-years(\d+)-(seed\d+|noseed)-"
    match1 = re.search(pattern1, test_path)
    print(f"  Uniform rate pattern match: {bool(match1)}")
    
    # Check gene rates pattern  
    pattern2 = r"gene_rates_[^/]+/grow(\d+)-sites(\d+)-years(\d+)-(seed\d+|noseed)-"
    match2 = re.search(pattern2, test_path)
    print(f"  Gene rates pattern match: {bool(match2)}")
    
    if match2:
        print(f"    Groups found: {match2.groups()}")