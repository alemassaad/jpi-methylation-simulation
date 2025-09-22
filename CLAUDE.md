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

# Complete pipeline (config loaded automatically, uniform mixing always enabled)
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz

# With custom parameters (override config defaults)
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz \
    --n-quantiles 10 --cells-per-quantile 3 --mix-ratio 80

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

## Phase 2 Pipeline Architecture

### Overview
Phase 2 generates structured datasets from Phase 1 simulations through a 6-stage pipeline. All individuals receive uniform mixing (same snapshot cells) for reproducible analysis.

### Pipeline Stages

#### Stage 1-2: Extract Snapshots (`extract_snapshots.py`)
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

**Output**:
- `snapshots/year{N}_snapshot.json[.gz]` - Direct copy of year data with year key
- `snapshots/metadata.json` - Configuration extracted from simulation

**Key Features**:
- No transformation - exact copy from phase1
- Year key included for self-contained files
- Preserves all fields (gene_jsd, etc.)

#### Stage 3-5: Simulate Individuals (`simulate_individuals.py`)

##### Stage 3: Create Initial Individuals
**Purpose**: Initialize populations from first snapshot

**Mutant Creation**:
- Uses quantile-based sampling for systematic coverage
- Cells sorted by methylation, divided into quantiles
- Sample specified cells per quantile
- Creates `n_quantiles × cells_per_quantile` individuals

**Control1 Creation**:
- Uses uniform random sampling
- Creates same number as mutant individuals
- No quantile structure, pure random selection

##### Stage 4: Grow Individuals
**Purpose**: Simulate population growth over time

**Process**:
1. Each individual grows for `timeline_duration` years
2. Exponential growth for `growth_phase` years
3. Homeostasis via random culling afterward
4. Target population: 2^growth_phase cells
5. Tracks full history during growth

**Timeline Calculation**:
```
timeline_duration = second_snapshot - first_snapshot
homeostasis_years = timeline_duration - growth_phase
```

##### Stage 5: Mix Populations
**Purpose**: Combine grown populations with second snapshot

**Uniform Mixing Process**:
1. **Normalize individuals**: All populations set to same size
2. **Create uniform pool**: Select cells from second snapshot
3. **Apply mixing**: Each individual gets identical snapshot cells
4. **Final composition**: `mix_ratio%` from snapshot, rest from grown

**Key Feature**: Uniform mixing ensures all individuals receive exactly the same snapshot cells, eliminating sampling variation

#### Stage 6: Create Control2 (`create_control2.py`)
**Purpose**: Pure snapshot populations without growth

**Process**:
1. Determine count based on mutant/control1 populations
2. If uniform mixing was used, apply same uniform base
3. Fill remaining slots with random sampling
4. No growth phase - pure snapshot cells
5. Matches experimental population sizes

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
  "mix_mode": "uniform",
  "mix_ratio": 80,
  "normalized": true,              # If size normalization applied
  "normalization_threshold": 128
}
```

### Configuration System

#### Default Configuration (`config_default.yaml`)
Always loaded automatically. Key parameters:
- `first_snapshot`: 30 (year)
- `second_snapshot`: 50 (year)
- `n_quantiles`: 4 (quartiles)
- `cells_per_quantile`: 3
- `individual_growth_phase`: 6 (64 cells)
- `mix_ratio`: 70 (percentage)
- `normalize_size`: true
- `compress`: false
- `seed`: 42

#### Configuration Priority
1. Default config loaded first
2. Custom config file overrides defaults
3. Command-line arguments override everything

### Output Directory Structure
```
data/{rate_info}/snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}u-seed{S}-{timestamp}/
├── snapshots/
│   ├── year30_snapshot.json[.gz]      # First snapshot cells
│   ├── year50_snapshot.json[.gz]      # Second snapshot cells
│   └── metadata.json                  # Gene rates, parameters
├── individuals/
│   ├── mutant/
│   │   └── individual_*.json[.gz]     # Mutant populations
│   ├── control1/
│   │   └── individual_*.json[.gz]     # Control1 populations
│   ├── control2/
│   │   └── individual_*.json[.gz]     # Control2 populations
│   └── mixing_metadata.json           # Mixing parameters, uniform pool
```

### Validation System

#### PipelineValidator Class
Comprehensive validation at each stage:

**Snapshot Validation**:
- Cell count sufficient for sampling
- Age consistency
- Gene_rate_groups match expectations
- Methylation biologically valid

**Individual Validation**:
- Correct count created
- Metadata properly set
- Gene rates consistent

**Growth Validation**:
- Population within expected range
- Extinction handling (configurable)
- History tracking complete

**Mixing Validation**:
- Uniform pool correctly applied
- Size normalization successful
- Final populations valid

### Core Module Functions

#### pipeline_utils.py
- `smart_open()`: Handle .json and .json.gz files
- `dict_to_cell()`: Convert JSON to Cell objects
- `load_snapshot_cells()`: Load snapshots with year key handling
- `sample_by_quantiles()`: Quantile-based sampling
- `sample_uniform()`: Random uniform sampling
- `normalize_populations()`: Size normalization
- `create_uniform_mixing_pool()`: Generate uniform cells
- `mix_petri_uniform()`: Apply uniform mixing
Note: `save_snapshot_cells()` removed - snapshots are now direct copies saved by extract_snapshots.py

#### individual_helpers.py
- `create_individual()`: Initialize PetriDish
- `grow_individual()`: Simulate growth
- `mix_individual()`: Mix with snapshot
- `process_batch_growth()`: Batch growth processing
- `process_batch_mixing()`: Batch mixing processing

#### path_utils.py
- `parse_phase1_simulation_path()`: Extract simulation parameters
- `generate_phase2_output_dir()`: Create output path with timestamp

#### validation.py
- `ValidationError`: Critical failures
- `ValidationWarning`: Non-fatal issues
- `ValidationConfig`: Configure strictness
- `PipelineValidator`: Main validation orchestrator

### Critical Implementation Details

#### Uniform Mixing Implementation
1. All individuals normalized to median size
2. Uniform pool created from second snapshot
3. Pool size = `normalized_size × mix_ratio / (1 - mix_ratio)`
4. Same pool cells added to every individual
5. Indices tracked in `mixing_metadata.json`

#### Extinction Handling
- Natural extinction allowed during growth
- Extinct individuals removed during normalization
- Validation tracks extinction rate
- Complete batch extinction triggers failure

#### Gene Rate Groups
- Every cell stores its own gene_rate_groups
- Validated at load time for consistency
- Passed through entire pipeline
- Critical for multi-rate simulations

### Common Issues and Solutions

#### Memory Management
- Large simulations may require batch processing
- Reduce n_quantiles or cells_per_quantile if needed
- Use compression for large datasets

#### Validation Failures
- "Gene rate groups mismatch": Inconsistent simulation
- "Insufficient cells": Snapshot too small for sampling
- "Complete batch extinction": Growth parameters too aggressive

#### Performance Tips
- Snapshots cached within run directory
- Use compressed files for large simulations
- Run stages individually for debugging
- Generate data once, analyze multiple times with phase3

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
│   ├── config_default.yaml     # Default configuration (loaded automatically)
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