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

## High-Level Architecture

### Three-Phase Architecture
- **Phase 1**: Core simulation engine - generates cell populations over time
- **Phase 2**: Data generation pipeline - creates structured datasets from phase1
- **Phase 3**: Analysis pipeline - generates all plots and statistics from phase2 data

### Phase 2 Pipeline Architecture

#### Overview
Phase 2 generates structured datasets from Phase 1 simulations through a 6-stage pipeline. All individuals receive common mixing (same snapshot cells) for reproducible analysis.

#### Pipeline Stages

##### Stage 1-2: Extract Snapshots (`extract_snapshots.py`)
**Purpose**: Extract cells from two time points in the simulation as direct copies

**Process**:
1. Load simulation file
2. Extract year data from history[year_str] (e.g., history["30"])
3. Save as direct copy with year key wrapper
4. Create metadata.json from simulation config

**Snapshot Format**:
```json
{
  "30": {
    "cells": [...],
    "gene_jsd": [...]  // Preserved if present
  }
}
```

##### Stage 3-5: Simulate Individuals (`simulate_individuals.py`)

###### Stage 3: Create Initial Individuals
**Mutant Creation**:
- Uses quantile-based sampling for systematic coverage
- Cells sorted by cell JSD, divided into quantiles
- Sample specified cells per quantile
- Creates `n_quantiles × cells_per_quantile` individuals

**Control1 Creation**:
- Uses common random sampling
- Creates same number as mutant individuals
- No quantile structure, pure random selection

###### Stage 4: Grow Individuals
1. Each individual grows for `timeline_duration` years
2. Exponential growth for `growth_phase` years
3. Homeostasis via random culling afterward
4. Target population: 2^growth_phase cells
5. Tracks full history during growth

###### Stage 5: Mix Populations
**Common Mixing Process**:
1. **Normalize individuals**: All populations set to same size
2. **Create common pool**: Select cells from second snapshot
3. **Apply mixing**: Each individual gets identical snapshot cells
4. **Final composition**: `mix_ratio%` from snapshot, rest from grown

##### Stage 6: Create Control2 (`create_control2.py`)
1. Determine count based on mutant/control1 populations
2. Load common pool from `common_pool.json`
3. Apply same common pool to all individuals
4. Fill remaining slots with random sampling if needed
5. No growth phase - pure snapshot cells

### Key Data Structures

#### Cell Object (from phase1)
```python
{
  "cpg_sites": [0, 0, 1, ...],           # Methylation state
  "cell_jsd": 0.123,                      # Jensen-Shannon divergence
  "age": 30,                              # Cell age in years
  "gene_rate_groups": [[50, 0.004], ...], # Gene-specific rates
  "cell_methylation_proportion": 0.15,    # Proportion methylated
  "n_methylated": 150                     # Count methylated
}
```

#### PetriDish Metadata
```python
{
  "individual_id": 1,
  "individual_type": "mutant",
  "source_quantile": 3,           # For mutant only
  "n_quantiles": 10,              # For mutant only
  "initial_year": 30,
  "mixed": true,
  "mix_mode": "common",
  "mix_ratio": 80,
  "normalized": true,              # If size normalization applied
  "normalization_threshold": 128
}
```

### Phase 2 Common Mixing System

#### Overview
Phase 2 implements a sophisticated common mixing system where all individuals receive identical snapshot cells. This ensures reproducibility and eliminates sampling variation between individuals. The common pool is saved to a file for reuse by Control2.

#### Mixing Implementation Flow

##### Step 1: Mandatory Size Normalization
**Location**: `phase2/core/pipeline_utils.py:1217-1351` (`normalize_populations()`)

**Algorithm**:
```python
threshold_size = median(all_sizes) - 0.5 * std_dev(all_sizes)
```

**Process**:
1. Calculate median and standard deviation of all population sizes
2. Set threshold = median - 0.5σ (excludes mild outliers)
3. Process each individual:
   - Size < threshold: EXCLUDE from dataset
   - Size = threshold: KEEP AS IS
   - Size > threshold: TRIM to threshold size (random sampling)
4. Create deep copies to preserve originals
5. Update metadata with normalization info

##### Step 2: Create and Save Common Mixing Pool
**Location**: `phase2/core/pipeline_utils.py:991-1046` (`create_common_mixing_pool()`)

**Key Features (2025-09-23)**:
- No longer returns indices (only pool cells)
- Throws `ValidationError` if insufficient cells (no sampling with replacement)
- Pool is saved to `common_pool.json` file

**Pool Size Calculation**:
```python
if mix_ratio >= 1.0:
    n_snapshot_cells = normalized_size  # Complete replacement
else:
    target_total = int(normalized_size / (1 - mix_ratio))
    n_snapshot_cells = int(target_total * mix_ratio)
```

##### Step 3: Apply Common Pool to All Individuals
**Location**: `phase2/core/pipeline_utils.py:1172-1214` (`mix_petri_common()`)

**Process**:
1. Validate pool compatibility with target cells
2. For 100% mix ratio: Replace all cells with snapshot
3. Otherwise: Add entire common pool to normalized individual
4. Shuffle cells for thorough mixing
5. Update metadata with mixing parameters

### Core Module Functions

#### pipeline_utils.py (Key Functions)
- `smart_open()`: Handle .json and .json.gz files
- `dict_to_cell()`: Convert JSON to Cell objects with gene_rate_groups
- `load_snapshot_cells()`: Load snapshots with year key handling
- `sample_by_quantiles()`: Quantile-based sampling for mutants
- `sample_common()`: Random common sampling for control1
- `normalize_populations()`: Normalize to consistent size
- `create_common_mixing_pool()`: Create shared pool of snapshot cells
- `save_common_pool()`: Save common pool to file for Control2 reuse
- `load_common_pool()`: Load common pool from saved file
- `validate_pool_compatibility()`: Validate pool/target cell compatibility
- `mix_petri_common()`: Apply common pool to individual
- `create_control2_with_common_base()`: Control2 with common base

#### individual_helpers.py
- `create_individual()`: Initialize PetriDish
- `grow_individual()`: Simulate growth with homeostasis
- `process_batch_growth()`: Batch growth processing

#### path_utils.py
- `parse_phase1_simulation_path()`: Extract simulation parameters
- `generate_phase2_output_dir()`: Create output path under Phase 1 directory with timestamp

#### validation.py
- `ValidationError`: Critical failures
- `ValidationWarning`: Non-fatal issues
- `ValidationConfig`: Configure strictness
- `PipelineValidator`: Main validation orchestrator

### Output Directory Structure (Sept 2025)
```
phase1/data/gene_rates_*/size*-seed*-{phase1_timestamp}/      # Phase 1 directory
├── simulation.json                                           # Phase 1 simulation
└── snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}u-seed{S}-{phase2_timestamp}/
    ├── snapshots/
    │   ├── year30_snapshot.json[.gz]      # First snapshot cells
    │   ├── year50_snapshot.json[.gz]      # Second snapshot cells
    │   └── metadata.json                  # Gene rates, parameters
    └── individuals/
        ├── mutant/
        │   └── individual_*.json[.gz]     # Mutant populations
        ├── control1/
        │   └── individual_*.json[.gz]     # Control1 populations
        ├── control2/
        │   └── individual_*.json[.gz]     # Control2 populations
        ├── common_pool.json[.gz]          # Shared snapshot cells
        └── mixing_metadata.json            # Mixing parameters
```

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
- **Note**: Verbose output is always enabled in phase2 (no configuration option)

#### Configuration Priority
1. Default config loaded first
2. Custom config file overrides defaults
3. Command-line arguments override everything

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
      "gene_jsd": [0.0, 0.0, ...]  // Optional: per-gene JSDs
    }
  }
}
```

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

## Critical Implementation Details

### Common Mixing Implementation
- Mandatory size normalization always applied
- Uses median - 0.5σ threshold method
- Excludes individuals below threshold
- Trims larger individuals to threshold size
- Ensures all individuals have exactly the same size

### Extinction Handling
- Natural extinction allowed during growth
- Extinct individuals removed during normalization
- Validation tracks extinction rate
- Complete batch extinction triggers failure

### Gene Rate Groups
- Every cell stores its own gene_rate_groups
- Validated at load time for consistency
- Passed through entire pipeline
- Critical for multi-rate simulations