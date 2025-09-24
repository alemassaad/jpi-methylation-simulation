# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands

### Run Simulation (Phase 1)
```bash
cd phase1

# Quick test (100 sites, 50 years, 64 cells)
python run_simulation.py --config config_default.yaml

# Full simulation with all options
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000 --seed 42

# Gene-specific methylation rates
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006" --gene-size 5

# Save output to phase1/data/
```

### Run Data Generation (Phase 2)
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

# Debug mixing
python simulate_individuals.py --base-dir data/my_run \
    --n-quantiles 4 --cells-per-quantile 3 \
    --growth-phase 6 --mix-ratio 80
```

### Run Analysis (Phase 3)
```bash
cd phase3

# Standard analysis on phase2 data
python run_analysis.py \
    --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/ \
    --simulation ../phase1/data/{phase1_dir}/simulation.json.gz

# Using config file
python run_analysis.py --config configs/quick_analysis.yaml \
    --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/ \
    --simulation ../phase1/data/{phase1_dir}/simulation.json.gz

# Plot phase1 simulation directly
python plot_simulation.py ../phase1/data/*/simulation.json.gz
```

## Phase 2 Data Flow Architecture - Detailed Mechanisms

### Overview
Phase 2 transforms Phase 1 simulations into structured datasets through a 6-stage pipeline. The pipeline uses common mixing where all individuals receive identical snapshot cells for reproducible analysis.

### Stage-by-Stage Data Flow

#### Stage 1-2: Extract Snapshots (`extract_snapshots.py`)

**Data Loading:**
1. Loads Phase 1 simulation file using `smart_open()` (handles .json/.json.gz)
2. Reads entire simulation into memory: `data = json.load(f)`
3. Accesses history with string keys: `history = data['history']`
4. Extracts specific year: `year_data = history[str(year)]` (e.g., `history["30"]`)

**Data Processing:**
- Direct copy of year data - NO transformation
- Preserves exact structure from Phase 1 including `cells` array and optional `gene_jsd`

**Data Writing:**
```python
# Snapshot format with year key wrapper
snapshot_with_key = {
    "30": {  # String year key
        "cells": [...],      # Direct copy from Phase 1
        "gene_jsd": [...]    # Preserved if present
    }
}
```

**Files Created:**
- `snapshots/year30_snapshot.json[.gz]` - First snapshot with year key
- `snapshots/year50_snapshot.json[.gz]` - Second snapshot with year key  
- `snapshots/metadata.json` - Configuration from simulation

**Metadata Structure:**
```json
{
    "gene_rate_groups": [[50, 0.004], [50, 0.006]],
    "n_sites": 1000,
    "gene_size": 5,
    "first_snapshot_year": 30,
    "second_snapshot_year": 50,
    "source_simulation": "../phase1/data/.../simulation.json",
    "extraction_timestamp": "1234567890.123"
}
```

#### Stage 3-5: Simulate Individuals (`simulate_individuals.py`)

##### Stage 3: Create Initial Individuals

**Data Loading:**
1. Loads metadata: `metadata = json.load(open('snapshots/metadata.json'))`
2. Loads first snapshot using `load_snapshot_cells()`:
   - Opens file with `smart_open()`
   - Detects year key wrapper: `if len(data) == 1 and list(data.keys())[0].isdigit()`
   - Strips wrapper: `year_data = data[year_key]`
   - Converts to Cell objects: `cells = [dict_to_cell(cd) for cd in year_data['cells']]`

**Cell Object Reconstruction (`dict_to_cell`):**
```python
# Each cell dictionary contains:
{
    "cpg_sites": [0, 1, 0, ...],           # Methylation state
    "gene_rate_groups": [[50, 0.004], ...], # Embedded in each cell
    "cell_jsd": 0.123,                      # Cell's JSD
    "age": 30,                              # Cell age
    "cell_methylation_proportion": 0.15,    # Cached calculation
    "n_methylated": 150                     # Cached count
}
```

**Sampling Process:**
- **Mutant**: `sample_by_quantiles()` - sorts by cell_jsd, divides into quantiles
- **Control1**: `sample_uniform()` - random sampling

##### Stage 4: Grow Individuals

**Processing:**
1. Each PetriDish grows for `timeline_duration` years
2. Growth tracked in `petri.cell_history` dictionary
3. Exponential phase: doubles each year for `growth_phase` years
4. Homeostasis: divide → cull ~50% → methylate

**Memory Management:**
- Full history maintained during growth
- Deep copies of cells at each year

##### Stage 5: Mix Populations  

**Normalization (`normalize_populations`):**
1. Calculate all population sizes
2. Compute threshold: `median - 0.5 * std_dev`
3. Process each individual:
   - Below threshold: EXCLUDE (remove from dataset)
   - At threshold: KEEP as is
   - Above threshold: TRIM to threshold (random sampling)
4. Returns normalized copies and threshold size

**Common Pool Creation (`create_common_mixing_pool`):**
```python
# Pool size calculation
if mix_ratio >= 1.0:
    n_snapshot_cells = normalized_size
else:
    target_total = int(normalized_size / (1 - mix_ratio))
    n_snapshot_cells = int(target_total * mix_ratio)
```

**Pool Saving (`save_common_pool`):**
```python
# Saves as simple array of cell dictionaries
pool_data = [cell.to_dict() for cell in pool]
# Written to: individuals/common_pool.json[.gz]
```

**Mixing Process (`mix_petri_common`):**
1. Validates compatibility with `validate_pool_compatibility()`
2. For 100% ratio: replaces all cells
3. Otherwise: adds entire pool to normalized individual
4. Shuffles combined cells

**Data Writing:**

Each individual saved as PetriDish with `save_petri_dish()`:
```json
{
    "cells": [...],           # Final mixed population
    "metadata": {
        // Essential identification
        "individual_id": 1,
        "individual_type": "mutant",  // or "control1", "control2"
        
        // Mutant-specific
        "source_quantile": 3,         // Only for mutant
        "n_quantiles": 10,            // Only for mutant
        
        // Growth parameters (mutant/control1 only)
        "initial_year": 30,           // Starting snapshot year
        "growth_phase": 6,            // Years of exponential growth
        
        // Mixing parameters (all types)
        "mix_ratio": 80,              // Percentage from snapshot
        "normalized_size": 128        // Size after normalization
    },
    "cell_history": {         # Optional, mutant/control1 only
        "30": [...],
        "31": [...],
        ...
    }
}
```

**Files Created:**
- `individuals/mutant/individual_01.json[.gz]` through `individual_N.json[.gz]`
- `individuals/control1/individual_01.json[.gz]` through `individual_N.json[.gz]`
- `individuals/common_pool.json[.gz]` - Shared snapshot cells (array format)
- `individuals/mixing_metadata.json` - Minimal mixing parameters: `{"mix_ratio": 80, "normalized_size": 128}`

#### Stage 6: Create Control2 (`create_control2.py`)

**Data Loading:**
1. Loads mixing metadata: `mixing_metadata = json.load(open('individuals/mixing_metadata.json'))`
2. Loads second snapshot cells
3. Loads common pool: `common_pool = load_common_pool(base_dir)`
   - Looks for `individuals/common_pool.json[.gz]`
   - Converts back to Cell objects

**Control2 Creation (`create_control2_with_common_base`):**
1. Starts with entire common pool (deep copy)
2. If need more cells: samples additional from snapshot
3. Creates PetriDish with combined cells

**Data Writing:**
Same PetriDish format but without `cell_history` (no growth phase)

### Key JSON Structures and Access Patterns

#### Phase 1 Simulation Format
```json
{
    "config": {
        "gene_rate_groups": [[50, 0.004], [50, 0.006]],
        "n": 1000,
        "gene_size": 5,
        "years": 100
    },
    "history": {
        "0": {                    // String keys for years
            "cells": [...],
            "gene_jsd": [...]     // Optional
        },
        "30": {...},
        "50": {...}
    }
}
```

#### Snapshot Format (Phase 2)
```json
{
    "30": {                      // Year key wrapper (string)
        "cells": [
            {
                "cpg_sites": [0, 1, ...],
                "gene_rate_groups": [[50, 0.004], ...],
                "cell_jsd": 0.123,
                "age": 30
            }
        ],
        "gene_jsd": [...]        // Optional, preserved from Phase 1
    }
}
```

#### Common Pool Format
```json
[                                // Direct array, no wrapper
    {
        "cpg_sites": [0, 1, ...],
        "gene_rate_groups": [[50, 0.004], ...],
        "cell_jsd": 0.456,
        "age": 50
    },
    ...
]
```

### Critical Implementation Details

#### File I/O Pattern
- Always use `smart_open()` for automatic .json/.json.gz handling
- Text mode for all JSON operations (`'rt'`/`'wt'` for gzip)
- Deep copies when loading/saving to avoid reference issues

#### Year Key Convention
- Phase 1 history uses string keys: `history["30"]`
- Snapshots preserve year key wrapper: `{"30": {...}}`
- Always convert: `year_str = str(year)`

#### Cell Compatibility Validation
- Every cell stores its own `gene_rate_groups`
- `PetriDish.validate_cells_compatible()` ensures consistency
- Validation at load, mix, and save points

#### Memory Optimization
- Snapshots cached within run (skip if exists)
- Compression determined from first snapshot file
- Deep copies only when necessary (mixing, normalization)

#### Common Pool Lifecycle
1. Created once during Stage 5
2. Saved to `common_pool.json[.gz]`
3. Loaded by Control2 in Stage 6
4. Ensures identical base across all individual types

### Configuration System

#### Default Configuration (`config_default.yaml`)
Always loaded automatically. Key parameters:
- `first_snapshot`: 20 (year)
- `second_snapshot`: 30 (year)
- `n_quantiles`: 4 (quartiles)
- `cells_per_quantile`: 2
- `individual_growth_phase`: 5 (32 cells)
- `mix_ratio`: 85 (percentage)
- `compress`: false
- `seed`: 24

#### Configuration Priority
1. Default config loaded first
2. Custom config file overrides defaults
3. Command-line arguments override everything

### Output Directory Structure
```
phase1/data/gene_rates_*/size*-seed*-{phase1_timestamp}/      # Phase 1 directory
├── simulation.json                                           # Phase 1 simulation
└── snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}u-seed{S}-{phase2_timestamp}/
    ├── snapshots/
    │   ├── year30_snapshot.json[.gz]      # First snapshot with year key
    │   ├── year50_snapshot.json[.gz]      # Second snapshot with year key
    │   └── metadata.json                   # Extraction metadata
    └── individuals/
        ├── mutant/
        │   └── individual_*.json[.gz]      # Mutant populations
        ├── control1/
        │   └── individual_*.json[.gz]      # Control1 populations
        ├── control2/
        │   └── individual_*.json[.gz]      # Control2 populations
        ├── common_pool.json[.gz]           # Shared snapshot cells (array)
        └── mixing_metadata.json            # Mixing parameters
```

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

## High-Level Architecture

### Three-Phase Architecture
- **Phase 1**: Core simulation engine - generates cell populations over time
- **Phase 2**: Data generation pipeline - creates structured datasets from phase1
- **Phase 3**: Analysis pipeline - generates all plots and statistics from phase2 data

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
**Note: Test files were removed in the 2025-01-20 cleanup. Tests are currently not available.**

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

## Recent Updates (2025)

### September 2025 Changes
- **Phase 2 Output Structure**: Phase 2 now outputs directly to Phase 1 simulation directory
- Removed `--output-dir` argument from Phase 2 pipeline
- Common mixing system improvements in Phase 2
- Pool saved to file instead of indices
- Added validation for pool compatibility
- Improved error handling for insufficient snapshot cells

### January 2025 Major Refactoring
- Removed ~100+ files, reduced to 24 essential files
- Moved all plotting from phase1 to phase3 (55% code reduction)
- Consolidated duplicate methods (single `calculate_gene_jsd()`, single `from_cells()`)
- Removed all test files as part of cleanup
- Clear separation: phase1 = simulation, phase3 = visualization

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

### Phase 2 Specific Issues

#### Memory Management
- Large simulations may require batch processing
- Reduce n_quantiles or cells_per_quantile if needed
- Use compression for large datasets

#### Validation Failures
- "Gene rate groups mismatch": Inconsistent simulation
- "Insufficient cells": Snapshot too small for sampling
- "Complete batch extinction": Growth parameters too aggressive

#### High Population Variation
- Normalization automatically applied
- Check growth parameters if variation persists

### General Issues
- **"No module named 'yaml'"**: Run `pip install pyyaml`
- **Plots not generating**: Install with `pip install plotly kaleido`
- **Out of memory**: Use smaller parameters or compression
- **Slow simulation**: Reduce sites with `--sites 100` or years with `--years 50`
- **ImportError in phase2/phase3**: Ensure sys.path additions are correct (see Import Structure section)

## Performance Tips
- Snapshots cached within run directory
- Use compressed files for large simulations
- Run stages individually for debugging
- Generate data once, analyze multiple times with phase3
- JSD calculations are always enabled (no performance flags needed)