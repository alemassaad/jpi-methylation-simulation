#!/usr/bin/env python3
"""
DEPRECATED: Plotting functionality has moved to phase3.

Use the new plotting interface:
  cd phase3
  python plot_simulation.py ../phase1/data/*/simulation.json.gz

Or for programmatic use:
  from phase3.core.simulation_plotter import plot_simulation_history
  plot_simulation_history('path/to/simulation.json.gz')

This placeholder script is kept for backward compatibility.
"""
import sys

print(__doc__)
print("\nError: plot_history.py has been moved to phase3.")
print("Please use the new interface as shown above.")
sys.exit(1)