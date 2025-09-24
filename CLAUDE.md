# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Phase 1: Run Simulation
```bash
cd phase1

# Quick test (100 sites, 50 years, 512 cells)
python run_simulation.py --config config_default.yaml

# Full simulation with all options
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000 --seed 42

# Gene-specific methylation rates
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006" --gene-size 5

# Output saved to data/gene_rates_*/size*-seed*-{timestamp}/simulation.json[.gz]
```

### Phase 2: Data Generation Pipeline
```bash
cd phase2

# Phase 2 outputs directly to Phase 1 directory - no output-dir needed!
# Complete pipeline (common mixing always enabled)
python run_pipeline.py --simulation ../phase1/data/gene_rates_*/size*-seed*-*/simulation.json

# With custom parameters (override config defaults)
python run_pipeline.py --simulation ../phase1/data/gene_rates_*/size*-seed*-*/simulation.json \
    --n-quantiles 10 --cells-per-quantile 3 --mix-ratio 80

# Run individual stages (for debugging/custom workflows)
python extract_snapshots.py --simulation ../phase1/data/*/simulation.json.gz --output-dir data/my_run
python simulate_individuals.py --base-dir data/my_run --n-quantiles 10 --cells-per-quantile 3
python create_control2.py --base-dir data/my_run
```

### Phase 3: Analysis and Visualization
```bash
cd phase3

# Standard analysis on phase2 data
python run_analysis.py \
    --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/ \
    --simulation ../phase1/data/{phase1_dir}/simulation.json.gz

# Using config file
python run_analysis.py --config configs/default.yaml \
    --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/ \
    --simulation ../phase1/data/{phase1_dir}/simulation.json.gz

# Plot phase1 simulation directly
python plot_simulation.py ../phase1/data/*/simulation.json.gz
```

## High-Level Architecture

### Three-Phase Pipeline
1. **Phase 1**: Core simulation engine - generates cell populations over time
2. **Phase 2**: Data generation pipeline - creates structured datasets from phase1
3. **Phase 3**: Analysis pipeline - generates all plots and statistics from phase2 data

### Key Classes
- `Cell`: Individual cell with methylation state and JSD calculations
- `PetriDish`: Population of cells with growth/homeostasis dynamics

### Directory Structure
```
phase1/data/gene_rates_*/size*-seed*-{phase1_timestamp}/      # Phase 1 directory
├── simulation.json                                           # Phase 1 simulation
└── snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}u-seed{S}-{phase2_timestamp}/
    ├── snapshots/
    │   ├── year{N}_snapshot.json[.gz]     # Snapshots with year key wrapper
    │   └── metadata.json                   # Extraction metadata
    └── individuals/
        ├── mutant/                         # Quantile-sampled populations
        ├── control1/                       # Random-sampled populations
        ├── control2/                       # Pure snapshot populations
        ├── common_pool.json[.gz]          # Shared snapshot cells (array)
        └── mixing_metadata.json           # Mixing parameters
```

## Phase 2 Data Flow - Key Implementation Details

### Snapshot Extraction (Stage 1-2)
- Loads Phase 1 simulation using `smart_open()` (handles .json/.json.gz)
- Accesses history with string keys: `history["30"]`
- Preserves year key wrapper in snapshots: `{"30": {...}}`
- Direct copy - NO transformation of cell data

### Individual Creation (Stage 3-5)
- **Mutant**: Quantile-based sampling (sorts by cell_jsd)
- **Control1**: Random uniform sampling
- Growth phase: exponential for N years, then homeostasis
- Common mixing: all individuals get identical snapshot cells
- Normalization: median - 0.5σ threshold, excludes/trims populations

### Control2 Creation (Stage 6)
- Pure second snapshot populations
- Uses common pool from Stage 5
- No growth phase (no cell_history)

### Critical JSON Structures

**Phase 1 History (string year keys):**
```json
{
    "history": {
        "0": {"cells": [...], "gene_jsd": [...]},
        "30": {...},
        "50": {...}
    }
}
```

**Snapshot Format (preserves year wrapper):**
```json
{
    "30": {
        "cells": [...],
        "gene_jsd": [...]
    }
}
```

**Common Pool (direct array):**
```json
[
    {"cpg_sites": [...], "gene_rate_groups": [...], "cell_jsd": 0.456, ...}
]
```

### Import Structure
```python
# Phase 2/3 must add phase1 to path
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish, Cell
```

## Project-Specific Instructions

### Critical Rules
- **ALWAYS use `python` instead of `python3`** for all commands
- **NEVER delete PLOT_DOCUMENTATION.md** - critical documentation for all plots
- **NO backward compatibility** unless explicitly requested - prefer clean breaks
- **Dictation note**: "jean" or "gin" means "gene"

### Configuration System
1. Default config loaded first (e.g., `config_default.yaml`)
2. Custom config file overrides defaults
3. Command-line arguments override everything

### File I/O Patterns
- Always use `smart_open()` for .json/.json.gz handling
- Text mode for JSON operations (`'rt'`/`'wt'` for gzip)
- Deep copies when necessary to avoid reference issues
- Year keys always strings: `year_str = str(year)`

### Validation Points
- Cell compatibility: check `gene_rate_groups` consistency
- Population size: normalization handles homeostasis variation
- Mixing validation: ensure pool compatibility before mixing

## Installation & Dependencies

```bash
pip install -r requirements.txt
```

Required packages:
- `numpy` (>=1.19.0): Performance optimization
- `scipy` (>=1.7.0): Statistical analysis
- `pyyaml` (>=6.0): Configuration files
- `plotly` (>=5.0.0, <6.0.0): Interactive visualization
- `kaleido` (0.2.1): PNG export from plotly

## Common Issues & Solutions

### Phase 2 Issues
- **Memory**: Reduce `n_quantiles` or `cells_per_quantile`
- **"Gene rate groups mismatch"**: Inconsistent simulation
- **"Insufficient cells"**: Snapshot too small for sampling
- **High variation**: Normalization automatically applied

### General Issues
- **ImportError**: Check sys.path additions
- **No plots**: Install `plotly kaleido`
- **Slow simulation**: Reduce sites or years
- **Out of memory**: Use compression or smaller parameters

## Performance Tips
- Snapshots cached within run directory
- Generate data once, analyze multiple times
- Use compressed files for large simulations
- Run stages individually for debugging

## Testing
**Note: Test files were removed in the 2025-01-20 cleanup. Tests are currently not available.**

## Recent Updates

### September 2025
- Phase 2 outputs directly to Phase 1 simulation directory
- Removed `--output-dir` argument from Phase 2
- Common mixing improvements
- Pool saved to file instead of indices

### January 2025
- Major refactoring: ~100 files reduced to 24
- Moved plotting from phase1 to phase3
- Consolidated duplicate methods
- Clear phase separation