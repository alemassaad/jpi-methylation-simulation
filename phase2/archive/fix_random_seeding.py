#!/usr/bin/env python3
"""
Fix script to add proper random seeding to phase2 pipeline.
This will modify the files to ensure reproducibility.
"""

import os
import sys

def apply_fixes():
    print("Applying random seed fixes to phase2...")
    
    # Fix 1: Add global seeding to run_pipeline.py
    print("\n1. Adding global seeding to run_pipeline.py...")
    
    run_pipeline_path = "run_pipeline.py"
    with open(run_pipeline_path, 'r') as f:
        lines = f.readlines()
    
    # Find where to insert (after argument parsing)
    insert_index = None
    for i, line in enumerate(lines):
        if "run_pipeline(args)" in line:
            # Go back to find the run_pipeline function definition
            for j in range(i, 0, -1):
                if "def run_pipeline(args):" in line[j]:
                    # Insert after the docstring and initial setup
                    for k in range(j, len(lines)):
                        if "start_time = time.time()" in lines[k]:
                            insert_index = k + 1
                            break
                    break
            break
    
    if insert_index:
        global_seed_code = """    
    # Set global random seeds for reproducibility (matching step23)
    import random
    import numpy as np
    random.seed(args.seed)
    np.random.seed(args.seed)
    print(f"Random seeds set to {args.seed}")
    
"""
        lines.insert(insert_index, global_seed_code)
        
        with open(run_pipeline_path, 'w') as f:
            f.writelines(lines)
        print(f"  ✓ Added global seeding after line {insert_index}")
    else:
        print(f"  ❌ Could not find insertion point")
    
    # Fix 2: Add np.random seeding to pipeline_utils.py
    print("\n2. Adding np.random.seed() to pipeline_utils.py functions...")
    
    utils_path = "pipeline_utils.py"
    with open(utils_path, 'r') as f:
        content = f.read()
    
    # Fix sample_uniform
    old_uniform = """def sample_uniform(cells: List[Cell], n_samples: int = 30, seed: int = 42) -> List[Cell]:
    \"\"\"
    Sample cells uniformly.
    
    Args:
        cells: List of Cell objects to sample from
        n_samples: Number of cells to sample
        seed: Random seed
    
    Returns:
        List[Cell]: Sampled cells
    \"\"\"
    random.seed(seed)"""
    
    new_uniform = """def sample_uniform(cells: List[Cell], n_samples: int = 30, seed: int = 42) -> List[Cell]:
    \"\"\"
    Sample cells uniformly.
    
    Args:
        cells: List of Cell objects to sample from
        n_samples: Number of cells to sample
        seed: Random seed
    
    Returns:
        List[Cell]: Sampled cells
    \"\"\"
    random.seed(seed)
    np.random.seed(seed)  # Also seed numpy for consistency with step23"""
    
    if old_uniform in content:
        content = content.replace(old_uniform, new_uniform)
        print("  ✓ Fixed sample_uniform")
    
    # Fix sample_by_quantiles
    old_quantiles = """    random.seed(seed)
    
    print(f\"  Sampling {cells_per_quantile} cells from each of {n_quantiles} quantiles...\")"""
    
    new_quantiles = """    random.seed(seed)
    np.random.seed(seed)  # Also seed numpy for consistency with step23
    
    print(f\"  Sampling {cells_per_quantile} cells from each of {n_quantiles} quantiles...\")"""
    
    if old_quantiles in content:
        content = content.replace(old_quantiles, new_quantiles)
        print("  ✓ Fixed sample_by_quantiles")
    
    # Fix mix_petri_with_snapshot
    old_mix = """    if seed is not None:
        random.seed(seed)"""
    
    new_mix = """    if seed is not None:
        random.seed(seed)
        np.random.seed(seed)  # Also seed numpy for consistency with step23"""
    
    if old_mix in content:
        content = content.replace(old_mix, new_mix)
        print("  ✓ Fixed mix_petri_with_snapshot")
    
    # Fix create_pure_snapshot_petri
    old_pure = """    if seed is not None:
        random.seed(seed)"""
    
    # Make sure we don't replace the same pattern twice
    content = content.replace("def create_pure_snapshot_petri", "MARKER_PURE")
    content = content.replace(old_pure, new_mix, 1)  # Replace first occurrence after marker
    content = content.replace("MARKER_PURE", "def create_pure_snapshot_petri")
    print("  ✓ Fixed create_pure_snapshot_petri")
    
    with open(utils_path, 'w') as f:
        f.write(content)
    
    print("\n✅ Fixes applied successfully!")
    print("\nTo use:")
    print("1. Clean individuals: ./clean_individuals.sh 0.005")
    print("2. Re-run pipeline: python run_pipeline.py --rate 0.005 --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz")
    print("3. Compare results with step23")

if __name__ == "__main__":
    response = input("This will modify run_pipeline.py and pipeline_utils.py. Continue? (y/n): ")
    if response.lower() == 'y':
        apply_fixes()
    else:
        print("Aborted.")