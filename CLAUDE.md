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

### Linting and Testing Commands
**This repository doesn't have configured linting or type checking commands. If you need to check code quality:**
```bash
# To add linting (if requested by user):
pip install ruff mypy
ruff check .
mypy phase1 phase2
```

## Installation
```bash
pip install -r requirements.txt  # Installs numpy, scipy, plotly, kaleido, pyyaml
```

Required dependencies:
- `numpy` (>=1.19.0): Optional but recommended for performance
- `scipy` (>=1.7.0): Required for statistical analysis in phase2
- `pyyaml` (>=6.0): For YAML config file support
- `plotly` (>=5.0.0, <6.0.0): For interactive visualizations
- `kaleido` (0.2.1): For PNG export from plotly

## Commands

### Run Simulation (Phase 1)
```bash
cd phase1

# Quick test with new defaults (100 sites, 50 years, 64 cells)
python run_simulation.py --rate 0.005

# Using config file (recommended)
python run_simulation.py --config configs/production.yaml

# Full command with all options
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000 --seed 42 --no-compress

# Gene-specific methylation rates
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006" --gene-size 5

# Performance options
python run_simulation.py --rate 0.005 --no-cell-jsds  # Skip cell JSD calculations
python run_simulation.py --rate 0.005 --no-gene-jsd   # Skip gene JSD tracking
python run_simulation.py --rate 0.005 --no-jsds       # Skip ALL JSDs
```

### Run Pipeline (Phase 2 - NEW MODULAR ARCHITECTURE)
```bash
cd phase2

# Complete pipeline with defaults
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz

# Standard analysis
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz \
  --first-snapshot 30 --second-snapshot 50 \
  --n-quantiles 10 --cells-per-quantile 3 \
  --individual-growth-phase 7 --mix-ratio 80

# Using config file (recommended)
python run_pipeline.py --config configs/quick_test.yaml --simulation ../phase1/data/*/simulation.json.gz

# With advanced options
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz \
  --first-snapshot 30 --second-snapshot 50 \
  --uniform-mixing --normalize-size --no-compress

# Run individual stages (for debugging/custom workflows)
python extract_snapshots.py --simulation ../phase1/data/*/simulation.json.gz --output-dir data/my_run
python simulate_individuals.py --base-dir data/my_run --n-quantiles 10 --cells-per-quantile 3
python create_control2.py --base-dir data/my_run
python analyze_and_plot.py --base-dir data/my_run --simulation ../phase1/data/*/simulation.json.gz
```

### Run Tests
```bash
# Phase 1 tests
cd phase1/tests
python test_small.py
python test_comprehensive.py
python test_edge_cases.py
python test_gene_jsd.py
python test_key_consistency.py
python test_rate_consistency.py

# Phase 2 tests  
cd phase2/tests
python test_gene_jsd_extraction.py
python test_validation.py
python test_cleanup_simple.py
python test_static_petridish.py
python test_rate_consistency_phase2.py

# Test configs
cd phase2/tests/config
python test_config_simple.py
python test_config_phase2.py
python test_pipeline_with_config.py
```

### Visualization Scripts
```bash
# Plot simulation history
cd phase1
python plot_history.py ../data/rate_0.00500/.../simulation.json.gz

# Generate individual plots for existing run
cd phase2
python plot_individuals.py data/.../results

# Analyze individual population sizes
cd phase2/visualization
python analyze_individual_sizes.py ../../data/.../individuals
```

## High-Level Architecture

### Core Classes (phase1/cell.py)
- **Cell**: Individual cell with CpG sites and methylation state
  - `methylated`: List of 0s and 1s representing methylation state
  - `cell_JSD`: Jensen-Shannon divergence from baseline (0-1 range)
  - `methylate()`: Apply stochastic methylation
  - `create_daughter_cell()`: Mitosis (identical copy)
  - `to_dict()`/`from_dict()`: Lean serialization

- **PetriDish**: Population manager for cells
  - Manages growth phase (exponential) and homeostasis (steady state)
  - `divide_cells()`: Population doubling
  - `random_cull_cells()`: Homeostatic ~50% survival
  - `calculate_gene_jsd()`: Calculate JSD for each gene across population
  - `save_history()`: Save in lean JSON format

### Pipeline Structure (phase2/) - NEW MODULAR ARCHITECTURE (2025-01-20)
Phase 2 has been reorganized into 4 independent scripts + 1 driver:

**Main Scripts**:
- `run_pipeline.py`: Main driver that orchestrates all stages
- `extract_snapshots.py`: Stages 1-2 (extract snapshots from phase1)
- `simulate_individuals.py`: Stages 3-5 (create, grow, mix populations)
- `create_control2.py`: Stage 6 (create control populations)
- `analyze_and_plot.py`: Stage 7 (all analysis and visualization)

**Core Modules**:
- `core/pipeline_utils.py`: Cell/PetriDish utilities, data I/O
- `core/pipeline_analysis.py`: All plotting and analysis functions
- `core/individual_helpers.py`: Individual creation and growth
- `core/validation.py`: Data validation functions
- `core/path_utils.py`: Path generation and parsing
- `core/plot_paths.py`: Plot organization and naming

**7-Stage Pipeline**:
1. **Stage 1-2**: Extract both snapshots (years 30 & 50)
2. **Stage 3**: Create initial individuals (quantile/uniform sampling)
3. **Stage 4**: Grow individuals (exponential + homeostasis)
4. **Stage 5**: Mix with snapshot cells
5. **Stage 6**: Create control2 individuals
6. **Stage 7**: Analysis and visualization

**Benefits**:
- Modular: Each script has single responsibility
- Restartable: Can resume from any stage
- Debuggable: Run stages individually
- Cacheable: Snapshots saved for reuse
- Flexible: Custom workflows possible

### Import Structure
```python
# Phase 1 imports
from cell import Cell, PetriDish

# Phase 2 imports (must add phase1 to path first)
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish, Cell
from pipeline_utils import (...)
```

## Key Implementation Details

### JSON Format (New Lean Format)
```json
{
  "parameters": {
    "rate": 0.005,              // Or null if gene-specific
    "gene_rate_groups": null,   // Or [[50, 0.004], [50, 0.006]]
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
          "cpg_sites": [0, 0, 1, ...],  // Methylation state
          "cell_JSD": 0.0                // Cell's divergence
        }
      ],
      "gene_jsd": [0.0, 0.0, ...]  // Optional: per-gene JSDs
    }
  }
}
```
- ~90% file size reduction compared to legacy format
- Parameters stored once at top level (not per cell)
- Compressed (.json.gz) by default for production

### Directory Structure
```
# Phase 1 output
phase1/data/gene_rates_200x0.00500/size8192-sites1000-genesize5-years100-seed42-YYYYMMDDHHMMSS/
  ├── simulation.json.gz      # Full history (lean format)
  ├── jsd_history.png         # Cell JSD trajectory plot
  └── methylation_history.png # Methylation trajectory plot

# Phase 2 output (NEW STRUCTURE)
phase2/data/{rate_info}/snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}-seed{S}-{timestamp}/
  ├── snapshots/              # Extracted snapshots (cached for reuse)
  │   ├── year30_snapshot.json.gz
  │   ├── year50_snapshot.json.gz
  │   └── metadata.json      # Gene rate groups, parameters
  ├── individuals/            # Simulated populations
  │   ├── mutant/            # Quantile-sampled individuals
  │   ├── control1/          # Uniformly-sampled individuals
  │   ├── control2/          # Pure snapshot individuals
  │   └── mixing_metadata.json  # Mixing configuration
  └── results/               # All analysis outputs
      ├── cell_metrics/      # Cell-level analysis
      │   ├── distributions/
      │   ├── comparisons/
      │   ├── individual_trajectories/
      │   └── timeline/
      └── gene_metrics/      # Gene-level analysis
          ├── distributions/
          ├── comparisons/
          ├── per_gene/
          └── timeline/
```

### Config File System
Both phases support YAML configuration with CLI override capability:
- Command-line arguments > User config > Default config > Hardcoded defaults
- Phase 1: `configs/production.yaml`, `configs/quick_test.yaml`
- Phase 2: `configs/config_default.yaml`, `configs/quick_test.yaml`, `configs/uniform_mixing.yaml`

### Performance Flags
- `calculate_cell_jsds=True`: Calculate individual cell JSDs (renamed from calculate_jsds)
- `track_gene_jsd=True`: Track population-level gene JSDs
- CLI: `--no-cell-jsds`, `--no-gene-jsd`, `--no-jsds` (disable ALL)

### Rate Consistency
All cells in a PetriDish must have identical `gene_rate_groups` configuration. Validation occurs:
- When creating PetriDish with existing cells
- Using `validate_cell_consistency()` instance method
- Using `validate_cells_compatible(cells)` static method

### Phase 2 Auto-Inference
Phase 2 automatically infers gene_rate_groups from the simulation file, eliminating the need to specify rates if using the same configuration as the simulation.

## Recently Fixed Issues

### Phase 2 Simulation Selection with Wildcards (FIXED)
**Previous Problem**: When using wildcards (`*`) in the `--simulation` flag, phase2 would find ALL matching simulations and automatically select the first one without confirmation.

**Solution Implemented**: Modified `phase2/run_pipeline.py` (lines 1793-1845) to:
- When a single file matches the pattern: Use it automatically
- When multiple files match: Display a numbered list of all matches and prompt the user to select one
- User can enter a number to select a specific simulation or 'q' to quit

**Current Behavior**: 
```bash
# Using wildcards is now safe - it will prompt for selection if multiple matches
python run_pipeline.py --simulation "../phase1/data/*/simulation.json.gz"

# Output when multiple matches:
# ========================================================
# MULTIPLE SIMULATIONS FOUND
# ========================================================
# Found 3 simulation files matching pattern: ../phase1/data/*/simulation.json.gz
# 
# Available simulations:
#   [1] rate_0.00400/size8192-sites1000-genesize5-years100-seed42-20250117123456/simulation.json.gz
#   [2] rate_0.00500/size8192-sites1000-genesize5-years100-seed42-20250117134521/simulation.json.gz
#   [3] rate_0.00600/size8192-sites1000-genesize5-years100-seed42-20250117145632/simulation.json.gz
#
# Please select a simulation by number (or 'q' to quit):
# Selection: 2
```

## Recent Breaking Changes

### Lean JSON Format (90% size reduction)
- New format stores parameters once at top level, not per cell
- No backward compatibility - old simulations must be re-run
- Compressed (.json.gz) by default for production

### Config File System
- Both phases support YAML configuration with CLI override capability
- Command-line arguments > User config > Default config > Hardcoded defaults
- Phase 1: `configs/production.yaml`, `configs/quick_test.yaml`
- Phase 2: `configs/config_default.yaml`, `configs/quick_test.yaml`

### JSD Naming Convention
- Cell-level JSD metrics prefixed with `cell_` (e.g., `cell_JSD`)
- Gene-level metrics prefixed with `gene_` (e.g., `gene_jsds`)
- Performance flag renamed: `calculate_jsds` → `calculate_cell_jsds`

### Data Structure Changes
- Individual IDs start at 1 (not 0) and match filenames throughout pipeline
- `gene_mean_methylation` returns proportions (0.0-1.0) not counts
- Timestamp-based directory names (YYYYMMDDHHMMSS) replace MD5 hashes

## Plot Architecture

### Consistent Plot Design (As of commit dfaa763)

All cell-level and gene-level plots now use **identical, consistent styling**:
- **Distribution plots**: Step plots with filled areas, mean lines, and comprehensive statistics boxes
- **Comparison plots**: Scatter points with jitter, quantile lines, and three-column statistics above plot
- **Individual trajectories**: Time series tracking metrics over years with consistent styling
- **Perfect symmetry**: Cell and gene metrics have identical plot types and styling

### Color Scheme Convention
**Metric-based coloring for consistency**:
- **JSD plots** (cell & gene): Blue (#1f77b4) with blue fill, red mean line
- **Methylation plots** (cell & gene): Red (#d62728) with red fill, dark blue mean line
- **Batch comparisons**: Mutant (#1f77b4), Control1 (#ff7f0e), Control2 (#2ca02c)

### Key Plotting Functions in phase2/core/pipeline_analysis.py

**Distribution plots (identical structure for cell and gene)**:
- `plot_cell_jsd_distribution()`: Line 26 - Blue step plot, 200 bins, red mean line
- `plot_gene_jsd_distribution()`: Line 324 - Blue step plot, 20 bins, red mean line  
- `plot_cell_methylation_proportion_histogram()`: Line 172 - Red step plot, 200 bins, blue mean line
- `plot_gene_methylation_proportion_histogram()`: Line 493 - Red step plot, 20 bins, blue mean line

**Comparison plots (all using same template design)**:
- `create_comparison_plot_from_jsds()`: Line 1234 - Template for batch comparisons
- `create_gene_individual_comparison_plot()`: Line 1464 - For gene-level individual averages
- `create_gene_methylation_comparison_plot()`: Line 2232 - For gene methylation averages
- `create_gene_comparison_plot()`: Line 1895 - For per-gene distributions

**Plot generation flow**:
- `plot_gene_jsd_distributions()`: Line 1842 - Creates individual plots for each gene
- Individual gene plots stored in: `results/gene_metrics/per_gene/gene_XXX_jsd.png`

**OOP Properties for clean access**:
- `cell.cell_jsd`: Individual cell's JSD value
- `petri.gene_jsds`: List of 20 gene JSD values (calculated property)
- `petri.gene_methylation_proportions`: List of 20 gene methylation proportions (calculated property)

### Statistical Annotations
All distribution plots include:
- Mean, Median, SD, CV (coefficient of variation)
- MAD (median absolute deviation)  
- Percentiles: 5%, 25%, 75%, 95%
- Statistics box in top-right corner with consistent formatting

All comparison plots include:
- Three-column statistics positioned above plot at y=1.02 (one column per batch)
- Mean lines (solid, thick)
- Quantile lines (25-75% dashed, 5-95% dotted)
- Jittered scatter points with consistent opacity

### Plot Organization (via PlotPaths class)
Plots are organized into metric-based subdirectories:
- `results/cell_metrics/`: All cell-level plots
  - `distributions/`: Snapshot histograms
  - `comparisons/`: Batch comparisons
  - `individual_trajectories/`: Per-individual time series
- `results/gene_metrics/`: All gene-level plots
  - `distributions/`: Snapshot histograms
  - `comparisons/`: Batch comparisons  
  - `per_gene/`: Individual gene distributions (gene_000_jsd.png through gene_019_jsd.png)
  - `individual_trajectories/`: Per-individual gene trajectories

## Biological Concepts

### CpG Sites
Cytosine-guanine dinucleotides where methylation occurs in DNA:
- Initially unmethylated (0)
- Can become methylated (1) with probability `rate` per year
- Once methylated, remain methylated (irreversible)
- Grouped into genes for distribution analysis

### Jensen-Shannon Divergence (JSD)
- **Cell-level JSD**: Measures individual cell divergence from baseline (0-1 range)
- **Gene-level JSD**: Measures population heterogeneity per gene

### Epigenetic Drift
Gradual accumulation of methylation changes over time:
- Random methylation events occur stochastically
- Methylation patterns inherited during cell division
- Population-level patterns emerge from single-cell dynamics

## Troubleshooting

### Common Issues and Solutions
- **"No module named 'yaml'"**: Run `pip install pyyaml`
- **Plots not generating**: Install plotly and kaleido with `pip install plotly kaleido`
- **Out of memory with large simulations**: Use `--no-jsds` flag to skip JSD calculations
- **Simulation running slow**: Reduce sites with `--sites 100` or years with `--years 50`

### Linting and Type Checking
This repository doesn't have configured linting or type checking commands. To add them:
```bash
pip install ruff mypy  # Install tools
ruff check .  # Run linter
mypy phase1 phase2  # Run type checker
```

## Dictionary Key Convention
All history dictionaries use STRING keys for years to ensure JSON compatibility:
- `cell_history` uses string keys: "0", "1", "2", etc.
- `gene_jsd_history` uses string keys: "0", "1", "2", etc.
- `mean_gene_jsd_history` and `median_gene_jsd_history` follow same convention
- **When sorting years**: Convert to int for sorting, then back to string for access
- **Example pattern**: 
  ```python
  years = sorted([int(y) for y in history.keys()])  # Sort as integers
  for year in years:
      data = history[str(year)]  # Access with string key
  ```
- **Important**: Always use `str(year)` when accessing dictionary values

## Recent Fixes

### Dictionary Key Type Consistency (Fixed 2025-01-18)
- **Issue**: Type comparison error in phase1 plots: `'<' not supported between instances of 'int' and 'str'`
- **Root cause**: `plot_gene_jsd_heatmap()` and `plot_gene_jsd_by_rate_group()` sorted string keys without conversion
- **Fix**: Convert keys to int for sorting, use string keys for access
- **Files changed**: `phase1/cell.py` (lines 2175-2184, 2299-2330), `phase2/core/pipeline_utils.py`

### Phase 2 Simulation Selection with Wildcards (Fixed)
- When using wildcards (`*`) in `--simulation` flag, prompts for selection if multiple files match
- Single match: Uses file automatically
- Multiple matches: Shows numbered list for user selection

### Gene JSD Data Key Types (Fixed 2025-01-18)
- **Issue**: Phase2 KeyError when generating gene JSD timeline plots
- **Root cause**: `extract_gene_jsd_from_history()` returned integer keys, but phase1's `plot_gene_jsd_timeline()` expected string keys
- **Fix**: Modified `extract_gene_jsd_from_history()` to use string keys consistently
- **Files changed**: `phase2/core/pipeline_utils.py`, `phase2/tests/test_gene_jsd_extraction.py`