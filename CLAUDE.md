# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Biologically realistic DNA methylation simulation modeling epigenetic drift through cell growth, division, and homeostasis. Uses object-oriented design with Cell and PetriDish classes to simulate how cells accumulate methylation patterns over time.

### Repository Status
- 8 untracked test files in phase2/tests/ directory (new tests not yet committed)
- Main branch: main
- Git repository: Yes
- Active development: Gene JSD plotting implementation (see GENE_JSD_IMPLEMENTATION_PLAN.md)

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

### Jensen-Shannon Divergence (cell_JSD)
A symmetric measure of difference between probability distributions:
- Stored as `cell_JSD` attribute on Cell objects (renamed from `JSD` for clarity)
- Quantifies how far a cell's methylation pattern has diverged from baseline
- Calculated as average of KL divergences to midpoint distribution
- Range: 0 (identical) to 1 (maximally different)
- Used to track epigenetic age and cellular heterogeneity
- Named `cell_JSD` to distinguish from `gene_JSD` calculations

### Gene-level JSD (gene_JSD)
Population-level measure of methylation heterogeneity for each gene:
- Calculated per gene across all cells in PetriDish
- Counts methylation levels (0-5 sites) across population
- Compares distribution to baseline (all unmethylated)
- Stored as list of JSDs, one per gene
- Tracks evolution of gene-specific methylation patterns
- Different from cell_JSD: gene_JSD measures population heterogeneity per gene, cell_JSD measures individual cell divergence
- Implementation in progress: See GENE_JSD_IMPLEMENTATION_PLAN.md for roadmap

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
PetriDish(rate=0.005, growth_phase=13)  # Population manager

# Gene-specific methylation rates
Cell(n=1000, gene_rate_groups=[(50, 0.004), (50, 0.006)])
PetriDish(gene_rate_groups=[(50, 0.004), (50, 0.006)], growth_phase=13)
```

Key methods:
- `Cell.methylate()`: Apply stochastic methylation (formerly age_1_year)
- `Cell.create_daughter_cell()`: Mitosis (identical copy)
- `PetriDish.divide_cells()`: Population doubling
- `PetriDish.random_cull_cells()`: Homeostatic ~50% survival
- `PetriDish.calculate_gene_jsd()`: Calculate JSD for each gene across population
- `PetriDish.enable_history_tracking()`: Enable cell/gene history tracking

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
  ├── simulation.json.gz      # Full history (all years)
  ├── jsd_history.png         # cell_JSD trajectory plot
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
      ├── mutant_00_jsd.png
      ├── mutant_00_methylation.png
      └── mutant_00_combined.png
```

## Commands

### Linting and Type Checking
```bash
# No automated linting or type checking is configured
# Manual code quality checklist:
# - Follow PEP 8 conventions (snake_case for functions/variables)
# - Use descriptive names (e.g., cell_JSD not JSD)
# - Add type hints for complex functions (Optional[Dict], List[float])
# - Ensure all tests pass before committing
```

### Run Simulation
```bash
cd phase1
# Uniform methylation rate (all genes methylate at same rate)
python3 run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Gene-specific methylation rates (different gene groups methylate at different rates)
python3 run_simulation.py --gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055" \
    --years 100 --growth-phase 13 --seed 42

# Gene JSD tracking enabled by default:
python3 run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Disable gene JSD tracking for performance:
python3 run_simulation.py --rate 0.005 --years 100 --no-gene-jsd --seed 42

# Quick test: --years 20 --growth-phase 4
# Performance mode: --no-jsds (disables all JSD calculations)
```

Expected output:
```
Running simulation with parameters:
  Methylation rate: 0.005
  CpG sites: 1000
  Gene size: 5
  Growth phase: 13 years (target: 8192 cells)
  Total years: 100
  Random seed: 42

Year 0: 1 cells
Year 1: 2 cells (growth phase)
Year 2: 4 cells (growth phase)
...
Year 13: 8192 cells (reached target)
Year 14: 8145 cells (steady state)
...
Simulation complete. Output saved to:
data/rate_0.00500/grow13-sites1000-years100-seed42-a3f7/
```

### Run Analysis Pipeline

**Uniform Rate (from Phase 1 --rate):**
```bash
cd phase2
python3 run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 --seed 42
```

**Gene-Specific Rates (NEW - from Phase 1 --gene-rate-groups):**
```bash
cd phase2
python3 run_pipeline.py --gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055" \
    --simulation ../phase1/data/gene_rates_*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 --seed 42
```

**With Uniform Mixing (FIXED - no more warnings):**
```bash
python3 run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/.../simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --uniform-mixing --seed 42
```

**With Individual Growth Plots:**
```bash
python3 run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/.../simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --plot-individuals --seed 42
```

**Generate plots for existing run:**
```bash
python3 plot_individuals.py data/.../pipeline_output_dir/
```

**Quick test parameters:** `--n-quantiles 4 --cells-per-quantile 1 --individual-growth-phase 2`
```

Expected pipeline output:
```
============================================================
TIMELINE BREAKDOWN
============================================================
📊 Individual Simulation Timeline:
   Year 50: Sample individual cells
   Years 50→57: Exponential growth (1 → 128 cells)
   Years 57→60: Homeostasis (~128 cells)
   Year 60: Extract for mixing

📈 Growth Summary:
   Total aging time: 10 years
   Exponential phase: 7 years (70%)
   Homeostasis phase: 3 years (30%)
   Target population: 128 cells (2^7)

================================================================================
PHASE 2 PIPELINE
================================================================================
Rate: 0.005 (or gene_rates_50x0.00400_50x0.00450_... for gene-specific)
Simulation: ../phase1/data/...
Output: data/...
Quantiles: 10 (deciles)
Cells per quantile: 3
Total individuals: 30 per group

============================================================
STAGE 1: Extract Year 50 Snapshot
============================================================
  Found 8192 cells (or loading from cache)

============================================================
STAGE 2: Plot JSD Distribution
============================================================
  Plotting Cell JSD distribution...
  Mean: 0.0456, Median: 0.0448, SD: 0.0089

[Stages 3-8 continue with similar formatting...]

Pipeline complete. Results in:
data/rate_0.00500-grow13-sites1000-years100/snap50to60-growth7-quant10x3-mix80-seed42-b8c9/
```

Note: Pipeline output is logged to `pipeline_output.log` in phase2 directory

### Run Tests
```bash
# Phase 1 tests (all are standalone scripts)
cd phase1/tests
python3 test_small.py           # Quick validation (10 years, 16 cells)
python3 test_comprehensive.py   # Full features
python3 test_edge_cases.py      # Edge cases
python3 test_gene_jsd.py        # Gene JSD tracking functionality

# Phase 2 tests  
cd phase2/tests
# Core functionality tests
python3 test_reproducibility_robust.py   # Reproducibility
python3 test_dynamic_mix_year.py        # Dynamic year calculations
python3 test_history_tracking.py        # History tracking and plotting
python3 test_normalization_comprehensive.py  # Size normalization

# Gene-specific rate tests
python3 test_gene_rate_basic.py         # Gene-specific rate support
python3 test_gene_rate_support.py       # Gene rate pipeline integration
python3 test_petri_gene_rates.py        # PetriDish gene rate functions

# Integration and pipeline tests
python3 test_final_integration.py       # Full pipeline integration test
python3 test_pipeline_start.py          # Pipeline initialization
python3 test_helper_functions.py        # Utility function tests
python3 test_plot_generation.py         # Plot generation tests
python3 test_real_path.py               # Path utility tests

# Run all uniform mixing tests
python3 run_all_uniform_tests.py        # Complete uniform mixing test suite
```

### Compare Pipeline Runs
```bash
cd phase2/tools
python3 compare_two_runs.py --dir1 path1 --dir2 path2
```

## Key Implementation Details

### Biological Model
- **Methylation**: Random with rate per site per year (default 0.005 = 0.5%)
- **Inheritance**: Daughter cells inherit parent's methylation pattern exactly
- **Population dynamics**: Growth phase → steady state via division and culling
- **Gene-level distribution**: Tracks methylation per gene (5 sites/gene by default)

### Development Roadmap
Active development is tracked in `GENE_JSD_IMPLEMENTATION_PLAN.md`:
- **Part 1**: Fix missing cell JSD labels in plots
- **Part 2**: Implement gene JSD visualizations (heatmaps, distributions, correlations)
- **Part 3**: Add gene JSD statistics and analysis tools
- Follows existing phase2 design patterns for consistency

### Variable Methylation Rates
The simulation supports two methylation modes:

**Uniform Rate (default):**
- All genes methylate at the same rate
- Specified via `--rate` parameter
- Example: `--rate 0.005` (all sites methylate at 0.5% per year)

**Gene-Specific Rates:**
- Different gene groups can have different methylation rates
- Specified via `--gene-rate-groups` parameter
- Format: `"n1:rate1,n2:rate2,..."` where n is number of genes
- Example: `"50:0.004,50:0.0045,50:0.005,50:0.0055"` 
  - First 50 genes: 0.4% per year
  - Next 50 genes: 0.45% per year
  - Next 50 genes: 0.5% per year
  - Last 50 genes: 0.55% per year
- Total genes must equal n_sites/gene_size
- Cannot specify both `--rate` and `--gene-rate-groups`
- All cells in a PetriDish share the same rate configuration

### Growth Phase vs Steady State

**Growth Phase (Years 0 to growth_phase):**
- Population doubles each year through synchronous division
- Deterministic population size: exactly 2^year cells
- No cell death or culling
- Models embryonic/developmental expansion
- Example: growth_phase=13 → 8192 cells at year 13

**Steady State (Years growth_phase+1 to T):**
- Population maintained around target size through homeostasis
- Each year: divide (double) → cull (~50% survival) → methylate
- Stochastic population varies around 2^growth_phase
- Models adult tissue maintenance
- Natural variation: ±10% typical

### Pipeline Stages (phase2)

1. **Extract First Snapshot** (e.g., year 50)
   - Load simulation data
   - Extract all cells at specified year
   - Cache in `snapshots/year{N}_snapshot.json.gz`

2. **Plot cell_JSD Distribution** 
   - Calculate cell_JSD for each cell
   - Create step histogram with 200 bins
   - Overlay statistics: Mean, Median, SD, CV, MAD, 5%, 95%
   - Save as `year{N}_jsd_distribution_200bins.png`

3. **Create Initial Individuals** (mutant: quantile-based, control1: uniform)
   - Sort cells by cell_JSD
   - Mutant: Sample from quantiles (e.g., 3 cells from each decile)
   - Control1: Random uniform sampling
   - Each individual starts as single sampled cell

4. **Grow N Years** (PetriDish simulation)
   - Individual growth phase: exponential to 2^phase cells
   - Then homeostasis: divide → cull → methylate
   - Track full history if --plot-individuals flag set

5. **Extract Second Snapshot** (first + growth_years)
   - Load year-matched snapshot
   - Cache for reuse

6. **Mix Populations** (default 80% snapshot, 20% grown)
   - Add snapshot cells to grown individuals
   - Options: uniform mixing (same cells), size normalization

7. **Create Control2** (pure second snapshot)
   - Same number of individuals
   - Same total cells per individual
   - Pure snapshot cells, no grown component

8. **Analysis** (t-tests, scatter plots)
   - Individual-level cell_JSD means
   - Cell-level cell_JSD distributions
   - Pairwise t-tests between groups
   - Scatter plots with group comparisons

### Reproducibility
- Set seeds globally AND in each function
- Use deep copy for Cell objects (avoid reference bugs)
- MD5 hashing for unique directory names
- Hierarchical caching of intermediate results
- Deterministic file I/O order

### Uniform Mixing Implementation
When using `--uniform-mixing` flag, uses a **3-step normalization approach** (fixed in latest version):

**Step 1: Normalize Individuals**
- All individuals normalized to same size via `normalize_individuals_for_uniform_mixing()`
- Uses minimum size across all individuals as target
- Random sampling ensures diversity while eliminating size variation

**Step 2: Create Uniform Pool**
- `create_uniform_mixing_pool()` creates shared snapshot cells
- Pool sized exactly for normalized individuals + mix ratio
- All individuals receive identical pool cells

**Step 3: Mix with Simplified Logic**
- `mix_petri_uniform()` adds entire pool to each individual
- All individuals end up with identical final size
- Eliminates "pool size mismatch" warnings

**Key Benefits:**
- **Perfect size consistency**: All final individuals identical size
- **Fair comparison**: No size-related sampling bias
- **No warnings**: Pool always matches individual needs
- **Preserved diversity**: Random sampling maintains cell variation

**Control2 Integration:**
- Uses same normalization-based size calculation
- `create_control2_with_uniform_base()` ensures size matching

### History Tracking and Plotting
- **Cell history**: Renamed from `history` to `cell_history` for clarity
- **Gene JSD history**: New `gene_jsd_history` tracks population-level gene JSDs
- **--plot-individuals flag**: Enables growth trajectory tracking for phase2 individuals
- **Gene JSD tracking**: Enabled by default in phase1 (use --no-gene-jsd to disable)
- **PetriDishPlotter class**: Unified plotting for both phases (cell_JSD, methylation, combined)
- **History format**: Year-indexed dictionary with cell states at each time point
- **Automatic plotting**: Generates plots at pipeline completion when flag is used
- **Standalone plotting**: Use `plot_individuals.py` for existing runs
- **Plot types**: cell_JSD trajectory, methylation trajectory, combined view

## Important Constants
```python
N = 1000                    # CpG sites per cell
RATE = 0.005               # Methylation rate (0.5%)
GENE_SIZE = 5              # Sites per gene
DEFAULT_GROWTH_PHASE = 13  # → 8192 cells (2^13)
DEFAULT_MIX_RATIO = 80     # 80% snapshot, 20% grown
DEFAULT_QUANTILES = 10     # Deciles for sampling
DEFAULT_CELLS_PER_QUANTILE = 3  # 30 total individuals
```

## Development Guidelines

### Testing Philosophy
- Tests are standalone Python scripts, not using pytest framework
- Each test file can be run directly: `python3 test_name.py`
- Tests include comprehensive output with ✓ marks for passed checks
- Return exit code 0 for success, 1 for failure
- Group related tests in single files with a `run_all_tests()` function

### Important Parameter Names (phase2)
- Use `--first-snapshot` and `--second-snapshot` (NOT --snapshot-year or --growth-years)
- Growth years are calculated as: second_snapshot - first_snapshot
- Use `--individual-growth-phase` for growth phase duration (exponential growth before homeostasis)
- Use `--plot-individuals` to enable growth trajectory tracking

### Object-Oriented Usage
```python
# Direct manipulation - no dict conversions needed
from pipeline_utils import load_snapshot_as_cells, grow_petri_for_years
cells = load_snapshot_as_cells(simulation_file, year=50)
petri = PetriDish(rate=0.005)
petri.cells = cells
grow_petri_for_years(petri, years=10)

# Access cell data
for cell in petri.cells:
    methylation_count = sum(cell.methylated)
    cell_jsd = cell.cell_JSD  # Direct attribute access
```

### Type Consistency
- Use float literals in distributions: `[1.0, 0.0, 0.0]` not `[1, 0, 0]`
- Deep copy cells when creating populations to avoid reference issues
- Dynamic year calculations, not hardcoded values
- Ensure numpy arrays use consistent dtypes

### File I/O Patterns
- Simulation data: `simulation.json.gz` (compressed JSON)
- Snapshots: `year{N}_snapshot.json.gz`
- PetriDish objects: `{name}_petri.json.gz`
- Individual objects: `individual_{i:02d}.json.gz`
- All paths use hierarchical parameter-based structure
- Use gzip compression for all JSON files

## Pipeline Decision Points

### When to Use Different Flags

**--plot-individuals**
- Use when: You need to understand individual growth trajectories
- Skip when: Running large-scale analyses (adds overhead)
- Output: Generates cell_JSD/methylation plots for each individual

**--uniform-mixing**
- Use when: Want perfect size consistency and eliminate sampling variation
- Skip when: Studying natural population variation or size effects
- Effect: **3-step normalization process** ensures identical final sizes
- **New Implementation (Fixed):**
  1. **Normalize**: All individuals reduced to minimum size via random sampling
  2. **Pool**: Create uniform snapshot pool sized for normalized individuals  
  3. **Mix**: Add identical pool to each individual → all same final size
- **All three groups** (mutant, control1, control2) end up identical size
- **Eliminates warnings**: No more "pool size mismatch" errors
- **Perfect comparison**: Size variation completely removed

**--normalize-size**
- Use when: Cell count variation affects your analysis
- Skip when: Natural size variation is important
- Effect: Trims/excludes to ensure uniform cell counts
- Note: Typically retains ~67% of individuals

**--individual-growth-phase**
- Small (3-5): Quick tests, small populations
- Medium (7-9): Standard analyses, ~100-500 cells
- Large (11-13): High-resolution, thousands of cells
- Trade-off: Larger = more biological realism but slower

**--n-quantiles and --cells-per-quantile**
- Testing: 4 quantiles × 1 cell = 4 individuals (fast)
- Standard: 10 quantiles × 3 cells = 30 individuals
- High-resolution: 20 quantiles × 5 cells = 100 individuals
- Trade-off: More individuals = better statistics but slower

## Performance Considerations

### Typical Runtimes
- Phase 1 simulation (100 years, 8192 cells): ~2-5 minutes
- Phase 2 pipeline (30 individuals, standard): ~5-10 minutes
- With --plot-individuals: Add ~2-3 minutes
- Large populations (32768 cells): Scale ~linearly

### Memory Usage
- Simulation history: ~100MB for 100 years, 8192 cells
- Pipeline intermediates: ~50MB per stage
- Cached snapshots: ~10MB per year
- Keep ~2GB RAM available for large runs

### Optimization Tips
- Use smaller --individual-growth-phase for testing
- Reduce --n-quantiles and --cells-per-quantile for prototyping
- Cached snapshots speed up repeated runs
- Parallelize multiple rate values externally

## Common Troubleshooting

### "File not found" for simulation
- Check path uses wildcards: `grow13-*/simulation.json.gz`
- Verify simulation completed successfully
- Ensure correct rate parameter matches directory

### Normalization excludes all individuals
- Reduce --individual-growth-phase (less variation)
- Increase --cells-per-quantile (more individuals)
- Check homeostasis isn't too aggressive

### ~~"Individual needs X but pool has Y" warnings (FIXED)~~
- ~~This was a bug in uniform mixing that has been fixed~~
- **Fixed in latest version** with 3-step normalization approach
- Old versions: Pool was created for median size but individuals varied
- **New behavior**: All individuals normalized first, then identical pool added

### Plots not generating
- Ensure plotly and kaleido are installed
- Check write permissions in output directory
- Verify --plot-individuals flag is set
- **Fixed**: "No history found" error resolved (attribute name fix)

### Gene rate compatibility
- Use `--gene-rate-groups` (not `--rate`) for gene-specific rate data
- Format: `"50:0.004,50:0.005,50:0.006"` (n_genes:rate pairs)
- Must match original Phase 1 simulation gene configuration

### Reproducibility issues
- Always set --seed explicitly
- Don't mix code versions
- Check numpy version consistency
- Avoid filesystem race conditions

### Memory errors
- Reduce population size (--growth-phase)
- Process fewer individuals
- Clear cached data periodically
- Use system with more RAM

## Dependencies
```bash
pip install -r requirements.txt
# Core: Python 3.7+ standard library (no external dependencies for basic simulation)
# Plotting: plotly>=5.0.0,<6.0.0, kaleido==0.2.1
# Analysis: scipy, numpy (imported dynamically if available)
# Note: scipy and numpy are optional but recommended for full functionality
```