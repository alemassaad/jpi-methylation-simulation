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
    """Load configuration from YAML file(s). Always loads defaults first."""
    config = {}
    
    if yaml is None:
        print("Warning: PyYAML not installed. Config file support disabled.")
        print("Install with: pip install pyyaml")
        return config
    
    # Always load default config first
    default_path = os.path.join(os.path.dirname(__file__), "config_default.yaml")
    if os.path.exists(default_path):
        try:
            with open(default_path, 'r') as f:
                default_config = yaml.safe_load(f) or {}
                config.update(default_config)
        except Exception as e:
            print(f"Warning: Could not load default config: {e}")
    else:
        print(f"Warning: Default config not found at {default_path}")
    
    # Load user config if provided (overrides defaults)
    if config_path and os.path.exists(config_path):
        try:
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f) or {}
                config.update(user_config)
        except Exception as e:
            print(f"Error loading config file {config_path}: {e}")
            sys.exit(1)
    elif config_path:
        print(f"Warning: Specified config file not found: {config_path}")
    
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
    parser.add_argument("--first-snapshot", type=int, default=None,
                       help="First year to extract (default from config: 30)")
    parser.add_argument("--second-snapshot", type=int, default=None,
                       help="Second year to extract (default from config: 50)")
    
    # Sampling parameters
    parser.add_argument("--n-quantiles", type=int, default=None,
                       help="Number of quantiles for sampling (default from config: 4)")
    parser.add_argument("--cells-per-quantile", type=int, default=None,
                       help="Cells to sample per quantile (default from config: 2)")
    
    # Growth parameters
    parser.add_argument("--individual-growth-phase", type=int, default=None,
                       help="Years of exponential growth (default from config: 6)")
    parser.add_argument("--mix-ratio", type=int, default=None,
                       help="Percentage of second snapshot in mix (default from config: 70)")
    
    # Mixing options (uniform mixing is now always enabled)
    
    # Other options
    parser.add_argument("--seed", type=int, default=None,
                       help="Random seed (default from config: 42)")
    
    # Compression options (mutually exclusive)
    compress_group = parser.add_mutually_exclusive_group()
    compress_group.add_argument("--compress", action='store_true', dest='compress',
                               help="Compress output files (.json.gz)")
    compress_group.add_argument("--no-compress", action='store_false', dest='compress',
                               help="Don't compress output files (.json)")
    parser.set_defaults(compress=None)  # Let config decide
    
    parser.add_argument("--force-recreate", action='store_true',
                       help="Force recreation of individuals")
    
    args = parser.parse_args()
    
    # Always load default config, then merge with CLI args
    config = load_config(args.config)  # Now always loads defaults
    
    # Map config keys to CLI argument names
    config_mapping = {
        'first_snapshot': 'first_snapshot',
        'second_snapshot': 'second_snapshot',
        'n_quantiles': 'n_quantiles',
        'cells_per_quantile': 'cells_per_quantile',
        'individual_growth_phase': 'individual_growth_phase',
        'mix_ratio': 'mix_ratio',
        'seed': 'seed',
        'compress': 'compress',
        'verbose': 'verbose'
    }
    
    # Apply config values where CLI args are None
    for config_key, arg_name in config_mapping.items():
        if config_key in config and hasattr(args, arg_name) and getattr(args, arg_name) is None:
            setattr(args, arg_name, config[config_key])
    
    # Ensure required values have defaults if not specified
    if args.first_snapshot is None:
        args.first_snapshot = 30
    if args.second_snapshot is None:
        args.second_snapshot = 50
    if args.n_quantiles is None:
        args.n_quantiles = 4
    if args.cells_per_quantile is None:
        args.cells_per_quantile = 2
    if args.individual_growth_phase is None:
        args.individual_growth_phase = 6
    if args.mix_ratio is None:
        args.mix_ratio = 70
    if args.seed is None:
        args.seed = 42
    if args.compress is None:
        args.compress = False  # Default if not in config
    
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
    print(f"Uniform mixing: Always enabled")
    print(f"Normalization: Always enabled (median - 0.5σ)")
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
    if args.compress:
        extract_cmd.append('--compress')
    else:
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
    # Uniform mixing is now always enabled
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
    print(f"  python run_analysis.py --phase2-dir {base_dir} --simulation {args.simulation}")
    print("=" * 80)


if __name__ == "__main__":
    main()