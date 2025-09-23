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

## Phase 2 Mixing System - Complete Architecture

### Overview
Phase 2 implements a sophisticated uniform mixing system where all individuals receive identical snapshot cells. This ensures reproducibility and eliminates sampling variation between individuals. The uniform pool is saved to a file for reuse by Control2.

### Mixing Implementation Flow

The mixing system is implemented in `phase2/simulate_individuals.py` (lines 319-433) and consists of:

#### Step 1: Mandatory Size Normalization
**Location**: `phase2/core/pipeline_utils.py:1217-1351` (`normalize_populations()`)

**Purpose**: Normalize all individuals to consistent size using median - 0.5σ threshold method

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

**Implementation** (`simulate_individuals.py:323-343`):
```python
mutant_dishes, control1_dishes, normalized_size = normalize_populations(
    mutant_dishes, control1_dishes,
    seed=args.seed + 5000
)
# Update metadata for all normalized individuals
```

#### Step 2: Create and Save Uniform Mixing Pool
**Location**: `phase2/core/pipeline_utils.py:991-1046` (`create_uniform_mixing_pool()`)

**Purpose**: Generate a shared pool of snapshot cells that all individuals will receive

**Key Changes (2025-09-23)**:
- No longer returns indices (only pool cells)
- Throws `ValidationError` if insufficient cells (no sampling with replacement)
- Pool is saved to `uniform_pool.json` file

**Pool Size Calculation**:
```python
if mix_ratio >= 1.0:
    n_snapshot_cells = normalized_size  # Complete replacement
else:
    target_total = int(normalized_size / (1 - mix_ratio))
    n_snapshot_cells = int(target_total * mix_ratio)
```

**Error Handling**:
```python
if n_snapshot_cells > len(snapshot_cells):
    raise ValidationError(
        f"Insufficient snapshot cells for uniform mixing: "
        f"need {n_snapshot_cells}, have {len(snapshot_cells)}. "
        f"Reduce mix_ratio or use a larger snapshot."
    )
```

**Implementation** (`simulate_individuals.py:348-360`):
```python
uniform_pool = create_uniform_mixing_pool(
    second_snapshot_cells,
    normalized_size,
    args.mix_ratio / 100,
    seed=args.seed + 1000
)

# Save the uniform pool to file
save_uniform_pool(uniform_pool, args.base_dir, compress=use_compression)

# Validate the pool
validator.validate_uniform_pool(uniform_pool, len(uniform_pool), second_year)
```

#### Step 3: Apply Uniform Pool to All Individuals
**Location**: `phase2/core/pipeline_utils.py:1172-1214` (`mix_petri_uniform()`)

**Purpose**: Mix each individual with the exact same pool of snapshot cells

**Key Changes (2025-09-23)**:
- Added validation before mixing via `validate_pool_compatibility()`
- Checks age consistency and gene_rate_groups compatibility

**Process**:
1. Validate pool compatibility with target cells
2. For 100% mix ratio: Replace all cells with snapshot
3. Otherwise: Add entire uniform pool to normalized individual
4. Shuffle cells for thorough mixing
5. Update metadata with mixing parameters

**Implementation** (`simulate_individuals.py:363-387`):
```python
for petri in mutant_dishes:
    final_size = mix_petri_uniform(petri, uniform_pool, args.mix_ratio / 100)
    petri.metadata['mixed'] = True
    petri.metadata['mix_mode'] = 'uniform'
    petri.metadata['mix_ratio'] = args.mix_ratio
```

### Key Mixing Functions

#### `normalize_populations()`
- **Purpose**: Normalize populations to consistent size
- **Returns**: `(normalized_mutant, normalized_control1, threshold_size)`
- **Key feature**: Deep copies preserve original data

#### `create_uniform_mixing_pool()` *(Updated 2025-09-23)*
- **Purpose**: Create shared pool of snapshot cells
- **Returns**: `List[Cell]` (pool cells only, no indices)
- **Key features**: 
  - Throws error if insufficient cells (no sampling with replacement)
  - Returns cells for direct saving to file

#### `save_uniform_pool()` *(New 2025-09-23)*
- **Purpose**: Save uniform pool to file for Control2 reuse
- **Location**: `pipeline_utils.py:1049-1079`
- **Saves to**: `individuals/uniform_pool.json[.gz]`
- **Includes**: Pool cells, cell count, timestamp

#### `load_uniform_pool()` *(New 2025-09-23)*
- **Purpose**: Load uniform pool from saved file
- **Location**: `pipeline_utils.py:1082-1118`
- **Returns**: `List[Cell]` from saved pool

#### `validate_pool_compatibility()` *(New 2025-09-23)*
- **Purpose**: Validate pool/target cell compatibility
- **Location**: `pipeline_utils.py:1121-1169`
- **Validates**: Age consistency, gene_rate_groups matching

#### `mix_petri_uniform()` *(Updated 2025-09-23)*
- **Purpose**: Apply uniform pool to individual
- **Returns**: Total cells after mixing
- **Key features**: 
  - Validates compatibility before mixing
  - Special handling for 100% mix ratio

### Metadata Tracking *(Updated 2025-09-23)*

**mixing_metadata.json Structure**:
```json
{
    "uniform_mixing": true,
    "mix_ratio": 80,
    "normalized_size": 100,
    "normalization_threshold": 133
}
```
*Note: `uniform_pool_indices` removed - pool saved directly to file*

**uniform_pool.json Structure** *(New)*:
```json
{
    "pool": [/* array of cell objects */],
    "n_cells": 400,
    "timestamp": "2025-09-23T00:18:52.345366"
}
```

**Individual PetriDish Metadata**:
```python
{
    'individual_id': 1,
    'individual_type': 'mutant',
    'source_quantile': 3,        # For mutants only
    'normalized': True,
    'normalization_threshold': 133,
    'mixed': True,
    'mix_mode': 'uniform',
    'mix_ratio': 80
}
```

### Edge Cases and Validation

1. **No individuals**: Return empty lists
2. **Single individual**: Skip normalization
3. **All excluded**: Trigger validation failure
4. **100% mix ratio**: Complete replacement logic
5. **Insufficient snapshot cells**: Sample with replacement
6. **High CV (>20%)**: Warning displayed

### Common Issues and Solutions

1. **High population variation**:
   - Normalization automatically applied
   - Check growth parameters if variation persists
   
2. **Memory issues**:
   - Reduce cells_per_quantile
   - Use compression (--compress)
   
3. **Validation failures**:
   - Check extinction rates
   - Verify snapshot has sufficient cells
   - Ensure growth parameters are reasonable

## Recent Updates (2025)

### September 2025 Changes
- **Phase 2 Output Structure**: Phase 2 now outputs directly to Phase 1 simulation directory
- Removed `--output-dir` argument from Phase 2 pipeline
- Uniform mixing system improvements in Phase 2
- Pool saved to file instead of indices
- Added validation for pool compatibility
- Improved error handling for insufficient snapshot cells

### January 2025 Major Refactoring
- Removed ~100+ files, reduced to 24 essential files
- Moved all plotting from phase1 to phase3 (55% code reduction)
- Consolidated duplicate methods (single `calculate_gene_jsd()`, single `from_cells()`)
- Removed all test files as part of cleanup
- Clear separation: phase1 = simulation, phase3 = visualization

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

# Phase 2 now outputs directly to Phase 1 directory - no output-dir needed!
# Complete pipeline (uniform mixing always enabled)
python run_pipeline.py --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json

# With custom parameters (override config defaults)
python run_pipeline.py --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json \
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
2. Apply same uniform pool (from mixing_metadata.json)
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
- `first_snapshot`: 20 (year)
- `second_snapshot`: 30 (year)
- `n_quantiles`: 3 (tertiles)
- `cells_per_quantile`: 2
- `individual_growth_phase`: 5 (32 cells)
- `mix_ratio`: 60 (percentage)
- `compress`: false
- `seed`: 42

#### Configuration Priority
1. Default config loaded first
2. Custom config file overrides defaults
3. Command-line arguments override everything

### Output Directory Structure (Updated Sept 2025)
```
data/gene_rates_*/size*-seed*-{phase1_timestamp}/      # Phase 1 directory
├── simulation.json                                    # Phase 1 simulation
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
        ├── uniform_pool.json[.gz]          # Shared snapshot cells
        └── mixing_metadata.json            # Mixing parameters
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

#### pipeline_utils.py (Mixing & Normalization Focus - Updated 2025-09-23)

**Normalization Functions**:
- `normalize_populations(mutant_dishes, control1_dishes, seed)`:
  - Calculates threshold: median - 0.5 * stdev
  - Excludes individuals below threshold (extinction handling)
  - Trims larger individuals to threshold size
  - Returns: (normalized_mutant, normalized_control1, threshold_size)

**Mixing Functions**:
- `create_uniform_mixing_pool(snapshot_cells, normalized_size, mix_ratio, seed)`:
  - Calculates required pool size from normalized size and ratio
  - **Updated**: Throws ValidationError if insufficient cells (no sampling with replacement)
  - **Updated**: Returns only pool cells (no indices)
  
- `save_uniform_pool(pool, base_dir, compress)` *(New)*:
  - Saves pool cells to `individuals/uniform_pool.json[.gz]`
  - Includes timestamp and cell count metadata
  
- `load_uniform_pool(base_dir)` *(New)*:
  - Loads pool from saved file
  - Converts dicts back to Cell objects
  
- `validate_pool_compatibility(pool, target_cells)` *(New)*:
  - Validates age consistency within pool
  - Validates gene_rate_groups compatibility
  
- `mix_petri_uniform(petri, uniform_pool, target_ratio)`:
  - **Updated**: Validates compatibility before mixing
  - Special handling for 100% mix ratio (complete replacement)
  - Shuffles cells for thorough mixing
  - Returns: total_cells_after_mixing

- `create_control2_with_uniform_base(snapshot_cells, base_dir, target_size, seed)`:
  - **Updated**: Loads pool from file instead of using indices
  - Creates Control2 with uniform base plus additional cells

**Statistics Functions**:
- `calculate_population_statistics(dishes, group_name)`: Calculate size statistics
- `print_mixing_statistics(mutant_stats, control1_stats, combined_median)`: Display pre-mixing stats

**Other Key Functions**:
- `smart_open()`: Handle .json and .json.gz files
- `dict_to_cell()`: Convert JSON to Cell objects
- `load_snapshot_cells()`: Load snapshots with year key handling
- `sample_by_quantiles()`: Quantile-based sampling for mutants
- `sample_uniform()`: Random uniform sampling for control1
- `save_all_individuals()`: Batch save after all processing complete

Note: `save_snapshot_cells()` removed - snapshots are now direct copies saved by extract_snapshots.py

#### individual_helpers.py
- `create_individual()`: Initialize PetriDish
- `grow_individual()`: Simulate growth
- `mix_individual()`: Mix with snapshot (deprecated - use uniform mixing)
- `process_batch_growth()`: Batch growth processing
- `process_batch_mixing()`: Batch mixing processing (deprecated)

#### path_utils.py
- `parse_phase1_simulation_path()`: Extract simulation parameters
- `generate_phase2_output_dir()`: Create output path under Phase 1 directory with timestamp

#### validation.py
- `ValidationError`: Critical failures
- `ValidationWarning`: Non-fatal issues
- `ValidationConfig`: Configure strictness
- `PipelineValidator`: Main validation orchestrator

### Critical Implementation Details

#### Uniform Mixing Implementation

**Normalization Process**:
- Mandatory size normalization always applied
- Uses median - 0.5σ threshold method
- Excludes individuals below threshold
- Trims larger individuals to threshold size
- Ensures all individuals have exactly the same size

**Uniform Pool Creation**:
- Calculates pool size based on normalized size and mix ratio
- Edge case handling for 100% mix ratio
- Samples from second snapshot with or without replacement
- Returns both cells and indices for Control2 creation

**Mixing Application**:
- Adds entire uniform pool to each normalized individual
- Shuffles cells for thorough mixing
- Updates metadata with mixing parameters
- Final composition: `mix_ratio%` from snapshot, rest from grown cells

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
- **Out of memory**: Use smaller parameters or compression
- **Slow simulation**: Reduce sites with `--sites 100` or years with `--years 50`
- **ImportError in phase2/phase3**: Ensure sys.path additions are correct (see Import Structure section)

