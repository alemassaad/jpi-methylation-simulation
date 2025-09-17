# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important User Instructions

### Permission Requirements
**NEVER make any code changes without explicit consent**
- Always explain what changes would fix an issue first, then wait for approval
- When encountering errors, first diagnose and explain the solution, showing exactly what would be changed  
- Only take destructive actions (delete, remove, overwrite) when explicitly instructed
- When asked about files or code, explain/analyze but don't modify unless specifically asked
- Don't commit or push unless explicitly requested
- Show: 1) What the problem is, 2) What files/functions would be modified, 3) What the changes would be - BEFORE making any edits

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

### Run Analysis Pipeline (Phase 2)
```bash
cd phase2

# Auto-infer rates from simulation (no rate needed)
python run_pipeline.py --simulation ../phase1/data/.../simulation.json.gz --first-snapshot 50 --second-snapshot 60

# Using config file (recommended)
python run_pipeline.py --config configs/quick_test.yaml --simulation ../phase1/data/.../simulation.json.gz

# With options
python run_pipeline.py --simulation ../phase1/data/.../simulation.json.gz --first-snapshot 50 --second-snapshot 60 --uniform-mixing --normalize-size --plot-individuals --no-compress
```

### Run Tests
```bash
# Phase 1 tests
cd phase1/tests
python test_small.py
python test_comprehensive.py
python test_edge_cases.py
python test_gene_jsd.py
python test_new_format.py
python test_config.py

# Phase 2 tests  
cd phase2/tests
python test_reproducibility_robust.py
python test_gene_rate_support.py
python test_final_integration.py
python verify_compression_fix.py
python verify_determinism.py

# Test configs
cd phase2/tests/config
python test_config_simple.py
python test_config_phase2.py
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

### Pipeline Structure (phase2/)
- **8-stage pipeline**: Snapshot extraction → sampling → growth → mixing → analysis
- **Core modules**:
  - `run_pipeline.py`: Main orchestrator
  - `core/pipeline_utils.py`: Cell/PetriDish utilities
  - `core/pipeline_analysis.py`: Visualization functions
  - `core/path_utils.py`: Path generation and parsing
  - `core/validation.py`: Data validation functions
  - `core/individual_helpers.py`: Individual creation utilities
  - `core/plot_paths.py`: Plot organization and naming

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

# Phase 2 output
phase2/data/rate_0.00500-grow13-sites1000-years100/snap50to60-growth7-quant10x3-mix80-seed42-XXXX/
  ├── individuals/            # Individual PetriDish objects (mutant/control1/control2)
  ├── snapshots/             # Extracted snapshot data
  ├── results/               # Analysis outputs and plots
  └── individual_plots/      # Growth trajectories (if --plot-individuals used)
```

### Config File System
Both phases support YAML configuration with CLI override capability:
- Command-line arguments > User config > Default config > Hardcoded defaults
- Phase 1: `configs/production.yaml`, `configs/quick_test.yaml`
- Phase 2: `configs/standard.yaml`, `configs/quick_test.yaml`

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

## Recent Breaking Changes

### Lean JSON Format
- New format reduces file sizes by ~90%
- No backward compatibility with old format
- Old simulations must be re-run

### Config System  
- Both phases now support YAML configuration files
- Preferred method for complex runs

### JSD Naming Convention
- All cell-level JSD metrics prefixed with `cell_`
- Gene-level metrics prefixed with `gene_`
- Renamed `calculate_jsds` → `calculate_cell_jsds`

### Individual ID Convention
- All individual IDs start at 1 (not 0)
- IDs always match filenames through all pipeline stages

### Gene Metrics as Proportions
- `gene_mean_methylation` returns proportions (0.0-1.0) instead of counts

### Timestamp-based Directory Names
- Both phases use YYYYMMDDHHMMSS timestamps instead of hashes
- Enables chronological sorting and easy identification

### Gene Methylation Proportion Plots
- Added gene-level methylation proportion comparison plots
- `gene_methylation_comparison.png` shows batch comparison with exact same layout as gene JSD plot
- `gene_methylation_analysis.json` contains per-individual statistics
- Individual trajectory plots include `_methylation.png` suffix for gene methylation

### Removed Plots
- `gene_jsd_snapshot_comparison.png` removed as not useful with only 20 data points
- Function `plot_gene_jsd_distribution_comparison` removed from pipeline

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