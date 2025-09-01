# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Important Guidelines

### Python Command Convention
**ALWAYS use `python` instead of `python3` for all commands in this repository.**

### No Backward Compatibility Policy
**DO NOT implement backward compatibility unless explicitly requested.**
- Backward compatibility creates multiple code paths for single features
- It increases complexity and introduces bugs
- Makes maintenance harder
- Clean breaks are preferred over compatibility layers
- When making breaking changes, document them clearly but don't maintain old behavior

## Project Overview

Biologically realistic DNA methylation simulation modeling epigenetic drift through cell growth, division, and homeostasis. Uses object-oriented design with Cell and PetriDish classes to simulate how cells accumulate methylation patterns over time.

### Repository Status
- Main branch: main
- Git repository: Yes
- **Latest Major Change**: Lean JSON format (90% file size reduction)
- **Config System**: Both phases now support YAML configuration files
- **Active Development**: Gene JSD analysis and visualization

## Recent Breaking Changes

### 🚨 New Lean JSON Format (No Backward Compatibility)
The simulation now uses a lean JSON format that reduces file sizes by ~90%:
- **Old format**: Stored redundant data (rate, gene_size, site_rates) in every cell
- **New format**: Stores parameters once, cells contain only essential data
- **Impact**: Phase 2 requires new format; old simulations must be re-run

### ✨ Config File System
Both phases now support YAML configuration files to reduce command-line complexity:
```bash
# Phase 1
python run_simulation.py --config configs/production.yaml

# Phase 2
python run_pipeline.py --config configs/quick_test.yaml --simulation PATH --rate 0.005
```

### 🎯 Clarified JSD Flags
- Renamed `calculate_jsds` → `calculate_cell_jsds` for clarity
- New independent control over different JSD types

### 🔧 Compression Consistency Fix
- Phase 2 now correctly respects compression settings for all batches (mutant, control1, control2)
- Automatically matches input file format (.json vs .json.gz) unless overridden with `--no-compress`
- All output files within a run now have consistent compression

### 📊 JSD Naming Convention Clarification
- All cell-level JSD metrics now prefixed with `cell_` for clarity
- Metadata keys: `mean_cell_jsd`, `std_cell_jsd`, `min_cell_jsd`, `max_cell_jsd`, `median_cell_jsd`
- Functions: `get_cell_jsd_array()`, `plot_cell_jsd_distribution()`
- Output file: `cell_jsd_comparison.png` (renamed from `jsd_comparison.png`)

### 🔢 Individual ID Convention
- All individual IDs now start at 1 (not 0) for intuitive numbering
- First individual: `individual_01.json` with `individual_id: 1`
- Consistent across all batches (mutant, control1, control2)
- IDs always match filenames through all pipeline stages

### 📈 Gene Metrics as Proportions
- `gene_mean_methylation` now returns proportions (0.0-1.0) instead of counts
- Makes metrics more interpretable (e.g., 0.5 = 50% methylated)
- Consistent with other proportion-based metrics
- Gene-level JSD already clearly named with `gene_` prefix
- **Breaking change**: Old JSON files will have old keys (`mean_jsd` etc.)

### 🧬 Gene-Level Metrics for Individuals
- Each individual now includes gene-level statistics after mixing
- New metadata fields:
  - `gene_jsds`: Array of JSD values, one per gene (heterogeneity measure)
  - `gene_mean_methylation`: Array of mean methylation levels (0-5), one per gene
  - `n_genes`: Number of genes analyzed
- Calculated only for final mixed populations (Stage 7)
- Enables analysis of gene-specific methylation patterns and rate effects

### 📦 Consolidated JSON Output
- **NEW**: `cell_jsd_analysis.json` consolidates all cell JSD data
- **Breaking change**: Removes `statistics.json` and `jsd_distributions.json` files
- Structure:
  - `summary_statistics`: Mean, std, median, min, max per batch
  - `statistical_tests`: T-test results between batches
  - `individual_means`: Mean cell JSD per individual
- Clearer semantic naming: `n_individuals` instead of `n`
- No backward compatibility - clean break from old format

## Biological Background

### CpG Sites
Cytosine-guanine dinucleotides where methylation occurs in DNA. In this simulation:
- Initially unmethylated (0)
- Can become methylated (1) with probability `rate` per year
- Once methylated, remain methylated (irreversible in this model)
- Grouped into genes for distribution analysis

### Epigenetic Drift
The gradual accumulation of methylation changes over time:
- Random methylation events occur stochastically
- Methylation patterns are inherited during cell division
- Population-level patterns emerge from single-cell dynamics
- Modeled as a Poisson-like process with fixed rate

### Cell-level Jensen-Shannon Divergence (cell_JSD)
A symmetric measure of difference between probability distributions:
- Stored as `cell_JSD` attribute on Cell objects
- Quantifies how far a cell's methylation pattern has diverged from baseline
- Calculated as average of KL divergences to midpoint distribution
- Range: 0 (identical) to 1 (maximally different)
- Used to track epigenetic age and cellular heterogeneity

### Gene-level Jensen-Shannon Divergence (gene_JSD)
Population-level measure of methylation heterogeneity for each gene:
- Calculated per gene across all cells in PetriDish
- Counts methylation levels (0-5 sites) across population
- Compares distribution to baseline (all unmethylated)
- Stored as list of JSDs, one per gene
- Tracks evolution of gene-specific methylation patterns
- Different from cell_JSD: gene_JSD measures population heterogeneity, cell_JSD measures individual divergence

## Architecture

### Import Structure
```python
# Phase 1 imports
from cell import Cell, PetriDish

# Phase 2 imports (must add phase1 to path first)
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish, Cell
from pipeline_utils import (...)
from pipeline_analysis import (...)
from path_utils import (...)
```

### Core Classes (phase1/cell.py)
```python
# Uniform methylation rate
Cell(n=1000, rate=0.005)        # Individual cell with CpG sites
PetriDish(rate=0.005, growth_phase=13, calculate_cell_jsds=True)  # Population manager

# Gene-specific methylation rates
Cell(n=1000, gene_rate_groups=[(50, 0.004), (50, 0.006)])
PetriDish(gene_rate_groups=[(50, 0.004), (50, 0.006)], growth_phase=13)
```

Key methods:
- `Cell.methylate()`: Apply stochastic methylation
- `Cell.to_dict()`: Lean serialization (only methylated array + cell_JSD)
- `Cell.from_dict()`: Deserialize from lean format
- `Cell.create_daughter_cell()`: Mitosis (identical copy)
- `PetriDish.divide_cells()`: Population doubling
- `PetriDish.random_cull_cells()`: Homeostatic ~50% survival
- `PetriDish.calculate_gene_jsd()`: Calculate JSD for each gene across population
- `PetriDish.calculate_gene_jsds()`: Alias for calculate_gene_jsd (clearer name)
- `PetriDish.calculate_gene_mean_methylation()`: Mean methylation level per gene
- `PetriDish.save_history()`: Save in new lean format

### JSON Format Structure (New Lean Format)

```json
{
  "parameters": {
    "rate": 0.005,                    // Or null if gene-specific
    "gene_rate_groups": null,         // Or [[50, 0.004], [50, 0.006]]
    "n": 1000,                        // CpG sites per cell
    "gene_size": 5,                   // Sites per gene
    "growth_phase": 13,               // Growth duration
    "years": 100,                     // Total simulation time
    "seed": 42,                       // Random seed
    "track_cell_history": true,       // Cell tracking enabled
    "track_gene_jsd": true           // Gene JSD tracking enabled
  },
  "history": {
    "0": {
      "cells": [
        {
          "methylated": [0, 0, 1, 0, 1, ...],  // Only essential data
          "cell_JSD": 0.0                      // Cell's divergence
        }
      ],
      "gene_jsd": [0.0, 0.0, ...]      // Optional: per-gene JSDs
    },
    "1": { ... },
    // ... more years
  }
}
```

**Key improvements over legacy format:**
- No redundant data per cell (rate, gene_size, site_rates removed)
- Parameters stored once at top level
- ~90% file size reduction (8KB → 0.8KB per cell)
- Cleaner structure with clear separation

### Main Pipeline Structure
- **phase1/**: Simulation engine (single cell → population growth → homeostasis)
  - Growth phase (years 0-growth_phase): Deterministic 2^year cells
  - Steady state (years growth_phase+1-T): Stochastic ~2^growth_phase cells
  - Key files: `cell.py` (core classes), `run_simulation.py` (CLI), `plot_history.py` (visualization)
- **phase2/**: Analysis pipeline (quantile sampling → growth → mixing → statistics)
  - 8-stage pipeline from snapshot extraction to statistical analysis
  - Key files: `run_pipeline.py` (main), `pipeline_utils.py` (helpers), `pipeline_analysis.py` (analysis)
  - `plot_individuals.py` - Standalone script to generate plots for existing runs
  - `path_utils.py` - Path generation and parsing utilities

### Directory Structure
```
# Phase 1 output:
phase1/data/rate_0.00500/grow13-sites1000-years100-seed42-XXXX/
  ├── simulation.json.gz      # Full history (lean format)
  ├── simulation.json         # If --no-compress used
  ├── jsd_history.png         # Cell JSD trajectory plot
  └── methylation_history.png # Methylation trajectory plot

# Phase 2 output:
phase2/data/rate_0.00500-grow13-sites1000-years100/snap50to60-growth7-quant10x3-mix80-seed42-XXXX/
  ├── individuals/            # Individual PetriDish objects
  │   ├── mutant/            # Quantile-sampled individuals
  │   ├── control1/          # Uniformly-sampled individuals
  │   └── control2/          # Pure second snapshot individuals
  ├── snapshots/             # Extracted snapshot data
  ├── results/               # Analysis outputs and plots
  └── individual_plots/      # Growth trajectories (if --plot-individuals used)
```

## Commands

### Configuration Files

Both phases support YAML configuration files with CLI override capability:

**Phase 1 Config Structure:**
```yaml
simulation:
  rate: 0.005           # Or use gene_rate_groups
  sites: 100           # CpG sites (default changed from 1000)
  years: 50            # Simulation duration (default changed from 100)
  gene_size: 5
  growth_phase: 6      # 2^6 = 64 cells (default changed from 13)

output:
  directory: "data"
  compress: false      # Default changed for easier inspection

performance:
  track_gene_jsd: true
  calculate_cell_jsds: true   # Renamed from calculate_jsds

seed: 42
```

**Phase 2 Config Structure:**
```yaml
input:
  simulation: null     # Required
  rate: 0.005         # Or gene_rate_groups

snapshots:
  first: 50
  second: 60

individuals:
  growth_phase: 7
  n_quantiles: 10
  cells_per_quantile: 3

mixing:
  ratio: 80
  uniform: false
  normalize_size: false

visualization:
  bins: 200
  plot_individuals: false

output:
  directory: "data"
  compress: true

seed: 42
verbose: false
```

### Run Simulation (Phase 1)

```bash
# Using config file (recommended)
python run_simulation.py --config configs/production.yaml

# Quick test with new defaults (faster)
python run_simulation.py --rate 0.005

# Full command with all options
python run_simulation.py \
    --rate 0.005 \
    --years 100 \
    --growth-phase 13 \
    --sites 1000 \
    --gene-size 5 \
    --seed 42 \
    --no-compress      # Save as .json instead of .json.gz

# Gene-specific methylation rates
python run_simulation.py \
    --gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055" \
    --gene-size 5

# Performance tuning with clarified flags
python run_simulation.py --rate 0.005 --no-cell-jsds  # Skip cell JSD calculations
python run_simulation.py --rate 0.005 --no-gene-jsd   # Skip gene JSD tracking
python run_simulation.py --rate 0.005 --no-jsds       # Skip ALL JSDs (maximum speed)
```

### Run Analysis Pipeline (Phase 2)

```bash
# Using config file
python run_pipeline.py \
    --config configs/quick_test.yaml \
    --simulation ../phase1/data/.../simulation.json.gz \
    --rate 0.005

# Standard run
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../phase1/data/.../simulation.json.gz \
    --first-snapshot 50 \
    --second-snapshot 60 \
    --seed 42

# With options
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../phase1/data/.../simulation.json.gz \
    --first-snapshot 50 \
    --second-snapshot 60 \
    --uniform-mixing \          # Use same snapshot cells for all
    --normalize-size \          # Normalize individual sizes
    --plot-individuals \        # Generate growth trajectories
    --no-compress              # Uncompressed output
```

### Run Tests
```bash
# Phase 1 tests
cd phase1/tests
python test_small.py
python test_comprehensive.py
python test_edge_cases.py
python test_gene_jsd.py

# Test new format
cd phase1
python test_new_format.py    # Tests lean JSON format
python test_config.py         # Tests config system

# Phase 2 tests  
cd phase2/tests
python test_reproducibility_robust.py
python test_gene_rate_support.py
python test_final_integration.py

# Test phase2 config
cd phase2
python test_config_simple.py
python test_config_phase2.py
```

## Key Implementation Details

### Performance Flags (Clarified)
```python
# Phase 1 PetriDish parameters
calculate_cell_jsds=True   # Calculate individual cell JSDs (renamed from calculate_jsds)
track_gene_jsd=True        # Track population-level gene JSDs

# CLI flags
--no-cell-jsds   # Disable cell JSD calculations only
--no-gene-jsd    # Disable gene JSD tracking only  
--no-jsds        # Disable ALL JSD calculations
--no-compress    # Save as .json instead of .json.gz
```

### File Size Comparison
**Legacy format (old):**
- Each cell: ~8KB of redundant data
- 8,192 cells: ~70MB uncompressed
- Compression ratio: ~10x

**Lean format (new):**
- Each cell: ~0.8KB (only essential data)
- 8,192 cells: ~7MB uncompressed  
- Compression ratio: ~70-75x
- **Overall: ~90% reduction in file size**

### Default Parameter Changes
Phase 1 defaults optimized for faster testing:
- `sites`: 1000 → 100 (10x faster)
- `years`: 100 → 50 (2x faster)
- `growth_phase`: 13 → 6 (128x smaller population)
- `compress`: true → false (easier inspection)

Production runs should override these with appropriate values.

### Compression Options
```bash
# Compressed (default for production)
--compress or compress: true in config
Output: simulation.json.gz (~70x smaller)

# Uncompressed (default for testing)
--no-compress or compress: false in config  
Output: simulation.json (human-readable)
```

## Development Guidelines

### Testing Philosophy
- Tests are standalone Python scripts, not using pytest framework
- Each test file can be run directly: `python test_name.py`
- Tests include comprehensive output with ✓ marks for passed checks
- Return exit code 0 for success, 1 for failure

### Config Priority
Configuration follows this priority (highest to lowest):
1. Command-line arguments
2. User config file (--config)
3. Default config file (config_default.yaml)
4. Hardcoded defaults in code

### Important Changes from Legacy
1. **JSON Format**: No backward compatibility - old simulations must be re-run
2. **Parameter Names**: `calculate_jsds` → `calculate_cell_jsds`
3. **File Sizes**: Expect ~90% smaller files with new format
4. **Config Files**: Preferred method for complex runs
5. **Defaults**: Optimized for testing, not production

## Common Troubleshooting

### "Year X not found" in Phase 2
- Simulation uses new lean format
- Re-run phase 1 simulation with latest code

### Large file sizes
- Ensure using compressed output (default in production)
- New format is ~90% smaller than legacy

### JSD values all zero
- Check if `--no-cell-jsds` or `--no-jsds` flag was used
- Verify `calculate_cell_jsds: true` in config

### Config not working
- Ensure PyYAML installed: `pip install pyyaml`
- Check YAML syntax (proper indentation)
- Verify config file path

## Performance Considerations

### Typical Runtimes (with new defaults)
- Small test: ~2 seconds (growth=2, years=10, sites=100)
- Standard test: ~5 seconds (growth=6, years=50, sites=100)
- Production: ~2-5 minutes (growth=13, years=100, sites=1000)

### Memory Usage (with lean format)
- Simulation history: ~10MB for 100 years, 8192 cells (was ~100MB)
- Cached snapshots: ~1MB per year (was ~10MB)
- Pipeline intermediates: ~5MB per stage (was ~50MB)

### Optimization Tips
- Use `--no-jsds` for maximum speed if JSDs not needed
- Use compressed output for production runs
- New lean format loads/saves much faster
- Config files reduce startup parsing time