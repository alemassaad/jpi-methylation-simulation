#!/usr/bin/env python3
"""
Plot phase1 simulation history.
Replacement for phase1/plot_history.py with all plotting functionality.

Usage:
    python plot_simulation.py ../phase1/data/*/simulation.json.gz
    python plot_simulation.py ../phase1/data/*/simulation.json.gz --output-dir plots/
    python plot_simulation.py path/to/simulation.json.gz --plots jsd methylation
"""
import argparse
import sys
import os

# Add current directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from core.simulation_plotter import plot_simulation_history


def main():
    parser = argparse.ArgumentParser(
        description='Plot simulation history from phase1 simulations',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Plot all available plots
  python plot_simulation.py ../phase1/data/*/simulation.json.gz
  
  # Plot specific types only
  python plot_simulation.py simulation.json.gz --plots jsd methylation
  
  # Save to specific directory
  python plot_simulation.py simulation.json.gz --output-dir plots/
  
  # Add custom title prefix
  python plot_simulation.py simulation.json.gz --title "Experiment 1"
        """
    )
    
    parser.add_argument('simulation', 
                       help='Path to simulation.json.gz file (supports wildcards)')
    
    parser.add_argument('--output-dir', '-o',
                       help='Output directory for plots (default: same as simulation)')
    
    parser.add_argument('--title', '-t',
                       help='Title prefix for plots')
    
    parser.add_argument('--plots', '-p', nargs='+',
                       choices=['jsd', 'methylation', 'combined', 'gene_jsd'],
                       help='Specific plots to generate (default: all available)')
    
    args = parser.parse_args()
    
    try:
        # Generate plots
        plot_simulation_history(
            simulation_path=args.simulation,
            output_dir=args.output_dir,
            title_prefix=args.title,
            plots_to_generate=args.plots
        )
        return 0
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        return 1
    except Exception as e:
        print(f"Error generating plots: {e}", file=sys.stderr)
        import traceback
        traceback.print_exc()
        return 1


if __name__ == '__main__':
    sys.exit(main())