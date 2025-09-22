# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Global User Instructions (All Projects)

### General Rules
1. **NEVER make any code changes without explicit consent** - always explain what changes would fix an issue first, then wait for approval
2. When encountering errors or problems, first diagnose and explain the solution, showing exactly what would be changed
3. Only take destructive actions (delete, remove, overwrite) when explicitly instructed
4. When asked about files or code, explain/analyze but don't modify unless specifically asked
5. Don't commit or push unless explicitly requested
6. Keep responses concise and to the point
7. Always show: 1) What the problem is, 2) What files/functions would be modified, 3) What the changes would be - BEFORE making any edits

### File Operations
- **Require explicit permission for**: edit, modify, write, create, delete, remove, overwrite, move
- **Safe operations (no permission needed)**: read, analyze, explain, list, search, grep, find

### Preferences
- **Verbosity**: concise
- **Proactive actions**: minimal

## Project-Specific Instructions

### Critical Documentation
**PLOT_DOCUMENTATION.md in the root directory is a critical document that should NEVER be deleted. It explains all output plots from the simulation and must be updated whenever new plots are added.**

### Dictation Note
**When the user says "jean" or "gin" via dictation, they mean "gene".**

### Python Command Convention
**ALWAYS use `python` instead of `python3` for all commands in this repository.**

### No Backward Compatibility Policy
**DO NOT implement backward compatibility unless explicitly requested.**
- Clean breaks are preferred over compatibility layers
- When making breaking changes, document them clearly but don't maintain old behavior

### Testing Commands
**Run tests with `python` (not python3) from the test directory:**
```bash
# Phase 1 Tests
cd phase1/tests
python test_small.py           # Quick validation
python test_comprehensive.py   # Full feature tests

# Phase 2 Tests  
cd phase2/tests
python test_validation.py      # Data validation functions
python test_pipeline_with_config.py  # Pipeline with config
```

## Installation
```bash
pip install -r requirements.txt  # Installs numpy, scipy, plotly, kaleido, pyyaml
```

Required dependencies:
- `numpy` (>=1.19.0): Optional but recommended for performance
- `scipy` (>=1.7.0): Required for statistical analysis
- `pyyaml` (>=6.0): For YAML config file support
- `plotly` (>=5.0.0, <6.0.0): For interactive visualizations
- `kaleido` (0.2.1): For PNG export from plotly

## Commands

### Run Simulation (Phase 1)
```bash
cd phase1

# Quick test (100 sites, 50 years, 512 cells)
python run_simulation.py --config config_default.yaml

# Full simulation with all options
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000 --seed 42

# Gene-specific methylation rates
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006" --gene-size 5

# JSD calculations are now always enabled - no performance flags needed
```

### Run Data Generation (Phase 2)
```bash
cd phase2

# Complete pipeline with defaults
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz

# Run individual stages (for debugging/custom workflows)
python extract_snapshots.py --simulation ../phase1/data/*/simulation.json.gz --output-dir data/my_run
python simulate_individuals.py --base-dir data/my_run --n-quantiles 10 --cells-per-quantile 3
python create_control2.py --base-dir data/my_run
```

### Run Analysis (Phase 3)
```bash
cd phase3

# Standard analysis on phase2 data
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{simulation_directory}/simulation.json.gz

# Using config file
python run_analysis.py --config configs/quick_analysis.yaml \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{simulation_directory}/simulation.json.gz
```

### Plot Phase 1 Simulation
```bash
# Use phase3 for plotting phase1 results:
cd phase3
python plot_simulation.py ../phase1/data/*/simulation.json.gz
```

## High-Level Architecture

### Three-Phase Architecture
- **Phase 1**: Core simulation engine - generates cell populations over time
- **Phase 2**: Data generation pipeline - creates structured datasets from phase1
- **Phase 3**: Analysis pipeline - generates all plots and statistics from phase2 data

### Core Classes (phase1/cell.py)

#### Cell Class
Individual cell with CpG sites and methylation state.

**Key Attributes (saved to file)**:
- `cpg_sites`: List of 0s and 1s representing methylation state
- `cell_jsd`: Jensen-Shannon divergence from baseline (0-1 range)
- `age`: Cell age in years
- `gene_rate_groups`: Rate configuration as [(n_genes, rate), ...]
- `cell_methylation_proportion`: Proportion of methylated sites (0-1)
- `n_methylated`: Count of methylated sites

**Key Methods**:
- `methylate()`: Apply stochastic methylation
- `create_daughter_cell()`: Mitosis (identical copy)
- `to_dict()`: Serialize for saving
- `from_dict()`: Deserialize from saved data

#### PetriDish Class
Population manager for cells.

**Life Cycle**:
- Years 0-growth_phase: Exponential growth (1 → 2^growth_phase cells)
- Years growth_phase+: Homeostasis (~2^growth_phase cells maintained)

**Key Methods**:
- `divide_cells()`: Population doubling
- `random_cull_cells()`: Homeostatic ~50% survival
- `calculate_gene_jsd()`: Calculate JSD for each gene across population
- `from_cells()`: Factory method - creates PetriDish from cell(s)
- `save_history()`: Save simulation to JSON/JSON.gz

### Import Structure
```python
# Phase 1 imports
from cell import Cell, PetriDish

# Phase 2 imports (must add phase1 to path first)
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish, Cell

# Phase 3 imports (must add both phase1 and phase2 to path)
sys.path.append(os.path.join(project_root, 'phase1'))
sys.path.append(os.path.join(project_root, 'phase2'))
```

## Current Repository Structure (After Major Cleanup)

```
jpi-methylation-simulation/
├── phase1/ (5 files)
│   ├── README.md
│   ├── cell.py                 # Core simulation engine
│   ├── config_default.yaml     # Default configuration
│   ├── plot_history.py         # Plot simulation history (deprecated - use phase3)
│   └── run_simulation.py       # Main simulation script
│
├── phase2/ (10 files)
│   ├── README.md
│   ├── __init__.py
│   ├── configs/
│   │   └── config_default.yaml # Default configuration
│   ├── core/ (5 files)
│   │   ├── __init__.py
│   │   ├── individual_helpers.py  # Individual creation helpers
│   │   ├── path_utils.py          # Path generation utilities
│   │   ├── pipeline_utils.py      # Data handling utilities
│   │   └── validation.py          # Data validation
│   ├── run_pipeline.py         # Main pipeline driver
│   ├── extract_snapshots.py    # Extract snapshots (stages 1-2)
│   ├── simulate_individuals.py # Create/grow populations (stages 3-5)
│   └── create_control2.py      # Create controls (stage 6)
│
└── phase3/ (9 files)
    ├── README.md
    ├── requirements.txt
    ├── configs/
    │   ├── default.yaml         # Default analysis config
    │   └── quick_analysis.yaml  # Quick analysis config
    ├── core/ (4 files)
    │   ├── __init__.py
    │   ├── analysis_functions.py  # All plotting/analysis functions
    │   ├── data_loader.py         # Data loading utilities
    │   └── plot_paths.py          # Plot path management
    ├── run_analysis.py          # Main analysis script
    └── plot_simulation.py        # Plot phase1 simulation results
```

## JSON Format (Lean Format)
```json
{
  "config": {
    "gene_rate_groups": [[50, 0.004], [50, 0.006]],  // Gene-specific rates
    "n": 1000,                  // CpG sites per cell
    "gene_size": 5,             // Sites per gene
    "growth_phase": 13,         // Growth duration
    "years": 100,               // Total simulation time
    "seed": 42
  },
  "history": {
    "0": {
      "cells": [
        {
          "cpg_sites": [0, 0, 1, ...],        // Methylation state
          "cell_jsd": 0.0,                    // Cell's divergence
          "age": 0,                            // Cell age
          "gene_rate_groups": [[50, 0.004], [50, 0.006]],
          "cell_methylation_proportion": 0.0,  // Proportion methylated
          "n_methylated": 0                    // Count methylated
        }
      ],
      "gene_jsd": [0.0, 0.0, ...]  // Optional: per-gene JSDs (if calculate_cell_jsds=True)
    }
  }
}
```

## Performance Note
JSD calculations are always enabled. Both cell-level and gene-level JSDs are computed for all simulations.

## Dictionary Key Convention
All history dictionaries use STRING keys for years to ensure JSON compatibility:
- Always use `str(year)` when accessing dictionary values
- Convert to int for sorting, then back to string for access
- Pattern: `years = sorted([int(y) for y in history.keys()])`

## Biological Concepts

### CpG Sites
Cytosine-guanine dinucleotides where methylation occurs in DNA:
- Initially unmethylated (0)
- Can become methylated (1) with probability `rate` per year
- Once methylated, remain methylated (irreversible)

### Jensen-Shannon Divergence (JSD)
- **Cell-level JSD**: Measures individual cell divergence from baseline (0-1 range)
- **Gene-level JSD**: Measures population heterogeneity per gene

## Common Issues and Solutions
- **"No module named 'yaml'"**: Run `pip install pyyaml`
- **Plots not generating**: Install with `pip install plotly kaleido`
- **Out of memory**: Use `--no-jsds` flag to skip JSD calculations
- **Slow simulation**: Reduce sites with `--sites 100` or years with `--years 50`
- **ImportError in phase2/phase3**: Ensure sys.path additions are correct (see Import Structure section)

## Cleanup Summary (2025-01-20)
Repository simplified from ~100+ files to 24 essential files:
- **Phase 1**: Reduced from 27 to 5 files (81% reduction)
- **Phase 2**: Reduced from 70+ to 10 files (86% reduction)  
- **Phase 3**: Already minimal with 9 files

All test files and redundant configs removed. Core functionality preserved.

## Major Refactoring (2025-01-20)

### Plotting Separation
- **Moved all plotting from phase1 to phase3**: phase1/cell.py reduced from 2,831 to 1,278 lines (55% reduction)
- **PetriDishPlotter**: Now in phase3/core/petri_dish_plotter.py
- **plot_history.py**: Functionality moved to phase3/plot_simulation.py
- **Clear separation**: phase1 = simulation, phase3 = visualization

### Method Consolidation
- **Gene JSD**: Single method `calculate_gene_jsd()` - removed duplicate `calculate_gene_jsds()`
- **Factory method**: Single `from_cells()` method handles both single cell and list inputs - removed `from_snapshot_cell()`
- **Cleaner API**: No redundant properties or duplicate methods