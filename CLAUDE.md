# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Biologically realistic DNA methylation simulation modeling epigenetic drift through cell growth, division, and homeostasis. Uses object-oriented design with Cell and PetriDish classes to simulate how cells accumulate methylation patterns over time.

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

## Architecture

### Core Classes (phase1/cell.py)
```python
Cell(n=1000, rate=0.005)        # Individual cell with CpG sites
PetriDish(rate=0.005, growth_phase=13)  # Population manager
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
- **phase2/**: Analysis pipeline (quantile sampling → growth → mixing → statistics)
  - 8-stage pipeline from snapshot extraction to statistical analysis

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

### Run Simulation
```bash
cd phase1
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# With gene JSD tracking:
python run_simulation.py --rate 0.005 --years 100 --track-gene-jsd --seed 42

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
```bash
cd phase2
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 --seed 42
# Quick test: --n-quantiles 4 --cells-per-quantile 1 --individual-growth-phase 2

# With individual growth plots (NEW):
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/.../simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --plot-individuals --seed 42

# Generate plots for existing run:
python plot_individuals.py data/.../pipeline_output_dir/
```

Expected pipeline output:
```
=== PHASE 2 PIPELINE STARTING ===
Rate: 0.005
Seed: 42

Stage 1: Extracting year 50 snapshot...
  Found 8192 cells

Stage 2: Plotting cell_JSD distribution...
  Mean cell_JSD: 0.0456
  Median cell_JSD: 0.0448
  SD: 0.0089

Stage 3: Creating individuals from quantiles...
  Mutant: 30 individuals (3 from each decile)
  Control1: 30 individuals (uniform sampling)

Stage 4: Growing individuals for 10 years...
  Growing from year 50 to year 60
  Individual growth phase: 7 years (target: 128 cells)
  
Stage 5: Extracting year 60 snapshot...

Stage 6: Mixing populations (80% snapshot, 20% grown)...
  Target size per individual: ~400 cells
  
Stage 7: Creating Control2 (pure snapshot)...

Stage 8: Statistical analysis...
  T-test Mutant vs Control1: p=0.0234
  T-test Mutant vs Control2: p=0.0012
  T-test Control1 vs Control2: p=0.4567

Pipeline complete. Results in:
data/rate_0.00500-grow13-sites1000-years100/snap50to60-growth7-quant10x3-mix80-seed42-b8c9/
```

### Run Tests
```bash
# Phase 1 tests (all are standalone scripts)
cd phase1/tests
python test_small.py           # Quick validation
python test_comprehensive.py   # Full features
python test_edge_cases.py      # Edge cases

# Phase 2 tests  
cd phase2/tests
python test_reproducibility_robust.py   # Reproducibility
python test_dynamic_mix_year.py        # Dynamic year calculations
python test_history_tracking.py        # History tracking and plotting
python test_normalization_comprehensive.py  # Size normalization
```

### Compare Pipeline Runs
```bash
cd phase2/tools
python compare_two_runs.py --dir1 path1 --dir2 path2
```

## Key Implementation Details

### Biological Model
- **Methylation**: Random with rate per site per year (default 0.005 = 0.5%)
- **Inheritance**: Daughter cells inherit parent's methylation pattern exactly
- **Population dynamics**: Growth phase → steady state via division and culling
- **Gene-level distribution**: Tracks methylation per gene (5 sites/gene by default)

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
When using `--uniform-mixing` flag:
- **Stage 6**: Creates a shared pool of cells from the second snapshot (80% of target size)
- **Mutant/Control1**: Receive this shared pool + their own grown cells
- **Control2**: Receives the same shared pool + additional unique snapshot cells per individual
- **Key functions**:
  - `create_uniform_mixing_pool()`: Returns both cells and indices used
  - `create_control2_with_uniform_base()`: Combines shared pool with additional cells
- **Benefits**: Eliminates inter-individual sampling variation in the base 80%, allowing cleaner comparison of the growth effect in the remaining 20%

### History Tracking and Plotting
- **Cell history**: Renamed from `history` to `cell_history` for clarity
- **Gene JSD history**: New `gene_jsd_history` tracks population-level gene JSDs
- **--plot-individuals flag**: Enables growth trajectory tracking for phase2 individuals
- **--track-gene-jsd flag**: Enables gene JSD tracking in phase1
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
- Use when: Want to eliminate sampling variation between individuals
- Skip when: Studying natural population variation
- Effect: All individuals share the same base snapshot cells
- Details: With this flag, all three groups (mutant, control1, control2) share the same 80% base cells from the second snapshot. The remaining 20% differs:
  - Mutant/Control1: Their own grown cells (different per individual)
  - Control2: Additional random snapshot cells (different per individual)
- This isolates the effect of growth vs. pure snapshot in the 20% portion

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

### Plots not generating
- Ensure plotly and kaleido are installed
- Check write permissions in output directory
- Verify --plot-individuals flag is set

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
# Core: Python 3.7+ standard library
# Plotting: plotly>=5.0.0,<6.0.0, kaleido==0.2.1
# Analysis: scipy, numpy (usually pre-installed)
# Note: scipy and numpy typically come with scientific Python distributions
```