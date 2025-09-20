#!/usr/bin/env python3
"""
Main driver script for the reorganized phase2 pipeline.
Orchestrates the execution of all pipeline stages.
"""

import argparse
import os
import sys
import subprocess
import time
import json
import glob
from typing import Optional, Dict, Any

# Try to import YAML support
try:
    import yaml
except ImportError:
    yaml = None

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from phase2.core.path_utils import parse_phase1_simulation_path, generate_phase2_output_dir


def load_config(config_path: Optional[str] = None) -> Dict[str, Any]:
    """Load configuration from YAML file(s)."""
    config = {}
    
    if yaml is None:
        return config
    
    # Load default config if it exists
    default_path = os.path.join(os.path.dirname(__file__), "configs", "config_default.yaml")
    if os.path.exists(default_path):
        try:
            with open(default_path, 'r') as f:
                default_config = yaml.safe_load(f) or {}
                config.update(default_config)
        except Exception as e:
            print(f"Warning: Could not load default config: {e}")
    
    # Load user config if provided
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f) or {}
                # Simple merge (no deep merge for simplicity)
                config.update(user_config)
        except Exception as e:
            print(f"Error loading config file {config_path}: {e}")
            sys.exit(1)
    
    return config


def merge_config_and_args(config: Dict, args: argparse.Namespace) -> argparse.Namespace:
    """Merge config values with command-line arguments."""
    for key, value in config.items():
        if not hasattr(args, key) or getattr(args, key) is None:
            setattr(args, key, value)
    return args


def run_command(cmd: list, description: str) -> bool:
    """Run a command and handle errors."""
    print(f"\n{'='*60}")
    print(f"Running: {description}")
    print(f"{'='*60}")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=False, text=True)
        print(f"✓ {description} completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ {description} failed with exit code {e.returncode}")
        return False
    except FileNotFoundError:
        print(f"❌ Could not find command: {cmd[0]}")
        return False


def handle_glob_patterns(simulation_path: str) -> str:
    """Handle glob patterns in simulation file path."""
    if '*' in simulation_path:
        matching_files = glob.glob(simulation_path)
        if not matching_files:
            print(f"Error: No files matching pattern: {simulation_path}")
            sys.exit(1)
        
        if len(matching_files) == 1:
            simulation_path = matching_files[0]
            print(f"Found simulation: {simulation_path}")
        else:
            # Multiple files found
            print(f"\n{'='*60}")
            print("MULTIPLE SIMULATIONS FOUND")
            print(f"{'='*60}")
            print(f"Found {len(matching_files)} simulation files")
            print("\nAvailable simulations:")
            
            matching_files.sort()
            for i, filepath in enumerate(matching_files, 1):
                if "phase1/data/" in filepath:
                    display_path = filepath.split("phase1/data/")[1]
                else:
                    display_path = os.path.basename(filepath)
                print(f"  [{i}] {display_path}")
            
            print("\nPlease select a simulation by number (or 'q' to quit):")
            while True:
                try:
                    selection = input("Selection: ").strip()
                    if selection.lower() == 'q':
                        print("Aborted by user.")
                        sys.exit(0)
                    
                    idx = int(selection) - 1
                    if 0 <= idx < len(matching_files):
                        simulation_path = matching_files[idx]
                        print(f"\nSelected: {simulation_path}")
                        break
                    else:
                        print(f"Invalid selection. Please enter 1-{len(matching_files)}")
                except ValueError:
                    print("Invalid input. Please enter a number or 'q'")
                except KeyboardInterrupt:
                    print("\n\nAborted by user.")
                    sys.exit(0)
    
    return simulation_path


def main():
    parser = argparse.ArgumentParser(
        description="Phase 2 Pipeline - Main Driver",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # Config file
    parser.add_argument("--config", type=str,
                       help="Path to YAML configuration file")
    
    # Required arguments
    parser.add_argument("--simulation", type=str, required=True,
                       help="Path to phase1 simulation file")
    
    # Snapshot parameters
    parser.add_argument("--first-snapshot", type=int, default=30,
                       help="First year to extract")
    parser.add_argument("--second-snapshot", type=int, default=50,
                       help="Second year to extract")
    
    # Sampling parameters
    parser.add_argument("--n-quantiles", type=int, default=10,
                       help="Number of quantiles for sampling")
    parser.add_argument("--cells-per-quantile", type=int, default=3,
                       help="Cells to sample per quantile")
    
    # Growth parameters
    parser.add_argument("--individual-growth-phase", type=int, default=7,
                       help="Years of exponential growth")
    parser.add_argument("--mix-ratio", type=int, default=80,
                       help="Percentage of second snapshot in mix")
    
    # Mixing options
    parser.add_argument("--uniform-mixing", action='store_true',
                       help="Use same snapshot cells for all individuals")
    parser.add_argument("--normalize-size", action='store_true',
                       help="Normalize all individuals to same size")
    
    # Visualization
    parser.add_argument("--bins", type=int, default=200,
                       help="Number of bins for histograms")
    parser.add_argument("--max-gene-plots", type=int, default=None,
                       help="Maximum number of gene plots")
    
    # Other options
    parser.add_argument("--seed", type=int, default=42,
                       help="Random seed")
    parser.add_argument("--output-dir", type=str, default="data",
                       help="Output directory")
    parser.add_argument("--no-compress", action='store_true',
                       help="Save uncompressed JSON files")
    parser.add_argument("--force-reload", action='store_true',
                       help="Force reload of snapshots")
    parser.add_argument("--force-recreate", action='store_true',
                       help="Force recreation of individuals")
    
    args = parser.parse_args()
    
    # Load config and merge
    config = load_config(args.config)
    args = merge_config_and_args(config, args)
    
    # Handle glob patterns
    args.simulation = handle_glob_patterns(args.simulation)
    
    # Validate simulation file
    if not os.path.exists(args.simulation):
        print(f"Error: Simulation file not found: {args.simulation}")
        sys.exit(1)
    
    if not (args.simulation.endswith('.json') or args.simulation.endswith('.json.gz')):
        print(f"Error: Simulation file must be .json or .json.gz")
        sys.exit(1)
    
    # Parse simulation parameters and generate output directory
    sim_params = parse_phase1_simulation_path(args.simulation)
    if not sim_params:
        print(f"Error: Could not parse simulation parameters")
        sys.exit(1)
    
    base_dir = generate_phase2_output_dir(args, sim_params)
    
    print("=" * 80)
    print("PHASE 2 PIPELINE - NEW ARCHITECTURE")
    print("=" * 80)
    print(f"Simulation: {args.simulation}")
    print(f"Output directory: {base_dir}")
    print(f"Snapshots: year {args.first_snapshot}, year {args.second_snapshot}")
    print(f"Growth phase: {args.individual_growth_phase} years")
    print(f"Mix ratio: {args.mix_ratio}%")
    print(f"Uniform mixing: {args.uniform_mixing}")
    print(f"Normalize size: {args.normalize_size}")
    print(f"Seed: {args.seed}")
    print("=" * 80)
    
    start_time = time.time()
    
    # Stage 1-2: Extract snapshots
    extract_cmd = [
        'python', 'extract_snapshots.py',
        '--simulation', args.simulation,
        '--first-snapshot', str(args.first_snapshot),
        '--second-snapshot', str(args.second_snapshot),
        '--output-dir', base_dir,
    ]
    if args.force_reload:
        extract_cmd.append('--force-reload')
    if args.no_compress:
        extract_cmd.append('--no-compress')
    
    if not run_command(extract_cmd, "Extract Snapshots (Stages 1-2)"):
        sys.exit(1)
    
    # Stage 3-5: Simulate individuals
    simulate_cmd = [
        'python', 'simulate_individuals.py',
        '--base-dir', base_dir,
        '--n-quantiles', str(args.n_quantiles),
        '--cells-per-quantile', str(args.cells_per_quantile),
        '--growth-phase', str(args.individual_growth_phase),
        '--mix-ratio', str(args.mix_ratio),
        '--seed', str(args.seed),
    ]
    if args.uniform_mixing:
        simulate_cmd.append('--uniform-mixing')
    if args.normalize_size:
        simulate_cmd.append('--normalize-size')
    if args.force_recreate:
        simulate_cmd.append('--force-recreate')
    
    if not run_command(simulate_cmd, "Simulate Individuals (Stages 3-5)"):
        sys.exit(1)
    
    # Stage 6: Create control2
    control2_cmd = [
        'python', 'create_control2.py',
        '--base-dir', base_dir,
        '--seed', str(args.seed + 300),
    ]
    if args.force_recreate:
        control2_cmd.append('--force-recreate')
    
    if not run_command(control2_cmd, "Create Control2 (Stage 6)"):
        sys.exit(1)
    
    # Analysis moved to phase3 - no longer part of phase2 pipeline
    
    # Summary
    elapsed_time = time.time() - start_time
    print("\n" + "=" * 80)
    print("PHASE 2 COMPLETE - Data Generation")
    print("=" * 80)
    print(f"Total time: {elapsed_time:.2f} seconds")
    print(f"Output directory: {base_dir}")
    print("\nGenerated directories:")
    print(f"  - snapshots/    : Extracted cell snapshots")
    print(f"  - individuals/  : Simulated populations")
    print("\nNext step: Run phase3 for analysis")
    print(f"  cd ../phase3")
    print(f"  python run_analysis.py --phase2-dir ../{base_dir} --simulation {args.simulation}")
    print("=" * 80)


if __name__ == "__main__":
    main()