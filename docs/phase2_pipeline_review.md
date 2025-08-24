# Phase 2 Pipeline Comprehensive Review

## Command-Line Arguments Analysis

### 1. REQUIRED ARGUMENTS

#### `--rate` (float, required)
**Purpose**: Methylation rate for creating new individuals
**Why Required**: Needed to create PetriDish objects with proper methylation parameters
**Usage in Pipeline**:
- Used when creating PetriDish objects from sampled cells
- Should match the original simulation's rate
**Potential Issue**: Could be auto-detected from simulation file to avoid mismatch

#### `--simulation` (str, required)
**Purpose**: Path to phase1 simulation.json.gz file
**Why Required**: Source of all cell data for analysis
**Usage in Pipeline**:
- Loaded to extract snapshots at specified years
- Contains complete cell history

### 2. QUANTILE SAMPLING PARAMETERS

#### `--n-quantiles` (int, default=10)
**Purpose**: Divides JSD distribution into N equal parts
**Why It Exists**: Allows sampling from extremes of JSD distribution
**Usage**:
- Creates "mutant" individuals by sampling from each quantile
- 10 = deciles, 4 = quartiles, etc.
**Logic**: More quantiles = better coverage of JSD spectrum

#### `--cells-per-quantile` (int, default=3)
**Purpose**: How many cells to sample from each quantile
**Why It Exists**: Controls sample size per quantile
**Usage**:
- Total mutant individuals = n_quantiles × cells_per_quantile
- Default: 10 × 3 = 30 mutant individuals
**Logic**: Balance between statistical power and computational cost

### 3. SNAPSHOT AND GROWTH PARAMETERS

#### `--first-snapshot` (int, default=50)
**Purpose**: Year to extract initial cell population
**Why It Exists**: Starting point for individual creation
**Usage**:
- Extract cells from this year
- Sample cells for creating individuals
- Starting age for growth simulations

#### `--second-snapshot` (int, default=60)
**Purpose**: Year to extract second cell population
**Why It Exists**: Comparison point after aging
**Usage**:
- Extract cells for mixing
- Defines total aging time (second - first)
- Creates control2 (pure snapshot) individuals

#### `--individual-growth-phase` (int, default=7)
**Purpose**: Years of exponential growth before homeostasis
**Why It Exists**: Controls population size (2^N cells)
**Usage**:
- Determines target population size
- Must be ≤ (second_snapshot - first_snapshot)
**Formula**: 7 → 128 cells, 8 → 256 cells

#### `--mix-ratio` (int, default=80)
**Purpose**: Percentage of snapshot cells in final mix
**Why Default 80%**: Balances snapshot vs grown cells
**Usage**:
- 80% from second snapshot
- 20% from grown individuals
- Applied during Stage 6

### 4. VISUALIZATION PARAMETERS

#### `--bins` (int, default=200)
**Purpose**: Number of bins for JSD distribution histogram
**Why 200**: Fine-grained visualization
**Usage**:
- Stage 2: Plot JSD distribution
- More bins = higher resolution

#### `--plot-individuals` (action='store_true', default=False)
**Purpose**: Generate growth trajectory plots
**Why It Exists**: Track individual growth dynamics
**Usage**:
- Enables history tracking during growth
- Creates plots in individual_plots/
- Adds computational overhead

### 5. MIXING PARAMETERS

#### `--uniform-mixing` (action='store_true', default=False)
**Purpose**: All groups share same 80% base cells
**Why It Exists**: Makes groups more comparable
**Usage**:
- TRUE: Same random 80% for all groups
- FALSE: Different random 80% per individual
**Impact**: Reduces inter-group variance

#### `--normalize-size` (action='store_true', default=False)
**Purpose**: Normalize population sizes before mixing
**Why It Exists**: Account for stochastic variation
**Usage**:
- Uses median - 0.5σ as threshold
- Equalizes population sizes
**Logic**: Adaptive threshold that excludes slightly below-average individuals while preserving most of the population

### 6. OTHER PARAMETERS

#### `--seed` (int, default=42)
**Purpose**: Random seed for reproducibility
**Why 42**: Hitchhiker's Guide reference
**Usage**:
- Set globally for all random operations
- Ensures reproducible sampling

#### `--output-dir` (str, default="data")
**Purpose**: Base directory for output
**Why "data"**: Standard convention
**Usage**:
- Creates hierarchical structure underneath

#### `--force-reload` (action='store_true', default=False)
**Purpose**: Bypass snapshot cache
**Why It Exists**: Force fresh load from simulation
**Usage**:
- Ignores cached snapshots
- Useful after simulation changes

#### `--force-recreate` (action='store_true', default=False)
**Purpose**: Recreate individuals even if exist
**Why It Exists**: Force regeneration
**Usage**:
- Overwrites existing individual files
- Useful for parameter changes

## Pipeline Flow Analysis

### STAGE 1: Extract First Snapshot
**What**: Load cells from first_snapshot year
**Files Created**: `snapshots/year{N}_snapshot.json.gz`
**Caching**: Skips if exists (unless --force-reload)

### STAGE 2: Plot JSD Distribution
**What**: Create histogram of cell JSDs
**Parameters Used**: --bins
**Files Created**: `results/year{N}_jsd_distribution_{bins}bins.png`

### STAGE 3: Create Initial Individuals
**What**: Sample cells to create individuals
**Mutant**: Quantile-based sampling (n_quantiles × cells_per_quantile)
**Control1**: Uniform random sampling (same count as mutant)
**Files Created**: `individuals/mutant/individual_XX.json.gz`

### STAGE 4: Grow Individuals
**What**: Grow each individual for timeline_duration years
**Growth**: individual_growth_phase years exponential
**Homeostasis**: Remaining years steady state
**History**: Tracked if --plot-individuals

### STAGE 5: Extract Second Snapshot
**What**: Load cells from second_snapshot year
**Files Created**: `snapshots/year{N}_snapshot.json.gz`
**Purpose**: Source for mixing

### STAGE 6: Mix Populations
**What**: Combine grown + snapshot cells
**Ratio**: --mix-ratio (default 80% snapshot)
**Modes**:
  - Independent: Each individual gets different random cells
  - Uniform: All share same random cells (--uniform-mixing)

### STAGE 7: Create Control2
**What**: Pure second snapshot individuals
**Purpose**: Baseline comparison (no growth)
**Uses**: Same base cells as others if --uniform-mixing

### STAGE 8: Statistical Analysis
**What**: Compare populations
**Outputs**:
  - T-tests between groups
  - Scatter plots
  - Summary statistics

## Issues Found

### 1. Parameter Name Inconsistencies
- `include_history` vs `include_cell_history` (FIXED)
- `total_growth_years` vs `timeline_duration` (FIXED)

### 2. Potential Improvements

#### Auto-detect Rate
The --rate parameter could be extracted from simulation file to avoid mismatches.

#### Normalize-size Logic
Uses median - 0.5σ threshold to exclude slightly below-average individuals while preserving most of the population. This adaptive threshold accounts for actual variance in the data.

#### Missing Parameters
- No way to specify control count independently
- No parameter for minimum cells in plot (hardcoded to 100)

### 3. Variable Naming Issues

#### In run_pipeline.py:
- `timeline_duration` - GOOD (replaced total_growth_years)
- `homeostasis_years` - GOOD
- `expected_population` - GOOD (clear meaning)
- `expected_individuals` - Could be `individuals_per_group`

#### In pipeline_utils.py:
Need to check function signatures for consistency

## Function Signature Review

### Key Functions to Check:
1. `load_petri_dish()` - uses `include_cell_history`
2. `load_all_petri_dishes()` - uses `include_history` 
3. `save_petri_dish()` - uses `include_cell_history`
4. `grow_petri_for_years()` - parameters?

## Recommendations

### High Priority Fixes:
1. ✅ Standardize parameter names across all functions
2. ✅ Fix timeline validation messages
3. ⬜ Check all function calls match signatures

### Medium Priority:
1. Auto-detect rate from simulation
2. Document normalization logic
3. Add --n-control parameter

### Low Priority:
1. Rename expected_individuals
2. Add --min-cells parameter for plots
3. Consider renaming mix-ratio to mix-percent

## Next Steps
1. Check all function signatures in pipeline_utils.py
2. Verify all function calls match signatures
3. Test complete pipeline end-to-end
4. Update documentation