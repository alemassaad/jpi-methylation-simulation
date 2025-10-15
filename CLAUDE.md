# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands - Quick Start

### Complete Pipeline (Most Common Workflow)
```bash
# 1. Run Phase 1 simulation
cd phase1
python run_simulation.py --config config_default.yaml

# 2. Generate Phase 2 datasets (outputs to root data directory)
cd ../phase2
python phase2_pipeline.py --simulation ../data/gene_rates_*/simulation.json.gz

# 3. Extract and plot with Phase 3 analysis
cd ../phase3
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/
```

## Commands - Detailed

### Phase 1: Run Simulation
```bash
cd phase1

# Quick test (100 sites, 50 years, 512 cells)
python run_simulation.py --config config_default.yaml

# Full simulation with all options
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000 --seed 42

# Gene-specific methylation rates (default configuration)
python run_simulation.py --gene-rate-groups "5:0.004,5:0.005,5:0.006,5:0.007" --gene-size 5

# Output saved to data/gene_rates_*/size*-seed*-{timestamp}/simulation.json[.gz]
```

### Phase 2: Data Generation Pipeline
```bash
cd phase2

# Complete pipeline - outputs directly to root data directory
python phase2_pipeline.py --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json

# With custom parameters (override config defaults)
python phase2_pipeline.py --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json \
    --n-quantiles 10 --cells-per-quantile 3 --mix-ratio 80

# Quick test (smaller parameters)
python phase2_pipeline.py --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json \
    --n-quantiles 4 --cells-per-quantile 2 --individual-growth-phase 6
```

### Phase 3: Analysis Pipeline
```bash
cd phase3

# Run complete analysis pipeline (6 stages)
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/

# Useful options
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/ --skip-plots  # Extract CSVs only
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/ --skip-gene-comparison  # Skip 40 per-gene plots

# Run individual analysis tools
python extract_simulation_timeline.py --simulation ../data/*/simulation.json.gz --output-dir tables/
python extract_batch_comparison.py --individuals-dir ../data/{phase1_dir}/{phase2_subdir}/individuals
python plot_comparison_generic.py --csv tables/batch_comparison_cell.csv --column cell_jsd_mean --output plots/cell_jsd.png
python plot_comparison_by_gene.py --csv tables/batch_comparison_gene.csv --output-dir plots/per_gene/
python plot_simulation_timeline.py results/tables/
```

## High-Level Architecture

### Three-Phase Pipeline
1. **Phase 1**: Core simulation engine - generates cell populations over time with methylation dynamics
2. **Phase 2**: Single-file data generation pipeline - creates structured datasets from phase1 simulations (no subprocesses, in-memory processing)
3. **Phase 3**: Analysis pipeline - extracts data tables and generates comprehensive plots

### Key Classes
- `Cell`: Individual cell with methylation state and JSD calculations (phase1/cell.py)
- `PetriDish`: Population of cells with growth/homeostasis dynamics (phase1/cell.py)

### Directory Structure
```
data/gene_rates_*/size*-seed*-{timestamp}/                     # Phase 1 output (at repo root)
├── simulation.json[.gz]                                       # Complete simulation
└── snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}-seed{S}-{timestamp}/  # Phase 2 output
    ├── snapshots/
    │   ├── year{N}_snapshot.json[.gz]
    │   └── metadata.json
    └── individuals/
        ├── test2/           # Test 2: Quantile-sampled populations
        ├── test1/           # Test 1: Random-sampled populations
        ├── control/         # Control: Pure snapshot populations
        ├── common_pool.json[.gz]
        └── mixing_metadata.json
```

### Phase 3 Pipeline Stages
```
Stage 1: Extract timeline from simulation → 4 CSV files
Stage 2: Generate 4 histogram plots (2 years × 2 metrics)
Stage 3: Generate 4 timeline plots from CSVs
Stage 4: Extract batch comparisons → 2 CSV files
Stage 5: Calculate statistical tests → 2 CSV files (pairwise t-tests + ANOVA)
Stage 6: Generate 6 comparison plots (4 mean + 2 std)
Stage 7: Generate 40 per-gene plots (20 genes × 2 metrics)
```

## Key Implementation Patterns

### Data Loading Pattern
```python
# Phase 3 must use 'individual_final' for mixed population
def load_petri_dish(filepath):
    with smart_open(filepath, 'r') as f:
        data = json.load(f)

    # CRITICAL: Use 'individual_final' (~210 cells after mixing with default config)
    # NOT 'history' (only ~68 cells pre-mixing)
    if 'individual_final' in data:
        cells = data['individual_final']['cells']  # CORRECT
```

### Smart File I/O
```python
def smart_open(filepath, mode='r'):
    """Transparently handle .json and .json.gz files"""
    if filepath.endswith('.gz'):
        import gzip
        return gzip.open(filepath, mode + 't')
    return open(filepath, mode)
```

### Import Structure
```python
# Phase 2/3 must add phase1 to path
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish, Cell
```

### Year Key Pattern
```python
# ALWAYS use string keys for years in history dictionaries
year_str = str(year)
history[year_str] = data  # CORRECT
# history[year] = data  # WRONG - will fail

# When iterating over years
for year_str, cells in history.items():
    year = int(year_str)  # Convert back to int if needed for calculations
```

### Cell Creation Pattern
```python
# Creating cells with gene-specific rates
cell = Cell(
    n=1000,
    gene_rate_groups=[(5, 0.004), (5, 0.005), (5, 0.006), (5, 0.007)],
    gene_size=5
)

# Creating PetriDish from cells (phase2 pattern)
petri = PetriDish.from_cells(cells, growth_phase=7)
```

### CSV Data Formats
- **Cell CSV** (N rows, where N = total individuals across batches): `individual_id,batch,cell_jsd_mean,cell_methylation_mean,cell_count`
- **Gene CSV** (N × 20 rows): `individual_id,batch,gene_index,gene_jsd,gene_methylation`
- **Timeline CSV**: `year,value_0,value_1,...` (cells) or `year,gene_0,gene_1,...` (genes)
- **P-values CSV** (18 rows): `metric,comparison,students_t_pvalue,welchs_t_pvalue`
- **ANOVA CSV** (12 rows): `metric,test_type,f_statistic,df1,df2,p_value`

## Phase 3 Comparison Plots with P-Values

### Overview
Phase 3 generates 6 comparison plots in Stage 5 (phase3/run_pipeline.py:425-598):
1. Cell JSD mean comparison
2. Cell methylation mean comparison
3. Gene JSD mean comparison
4. Gene methylation mean comparison
5. Gene JSD standard deviation comparison
6. Gene methylation standard deviation comparison

### P-Value Implementation Location
P-values should be added in `phase3/plot_comparison_generic.py` where batch comparisons are visualized. The function `plot_comparison_generic()` handles all 6 plots.

### Key Implementation Points

#### Data Structure
- Each batch (control, test1, test2) contains a variable number of individuals
  - Initial: `n_quantiles × cells_per_quantile` per batch (from config)
  - Final: After normalization filtering (median - 0.5σ threshold)
  - Typical: ~40 individuals per batch with default config (10 quantiles × 5 cells)
- Cell metrics: `cell_jsd_mean`, `cell_methylation_mean` per individual
- Gene metrics: `gene_jsd`, `gene_methylation` (20 values per individual, aggregated as mean or std)

#### Where to Add P-Values
In `plot_comparison_generic.py`:
1. **After line 201** where statistics are calculated - add statistical tests
2. **Around lines 256-280** where annotations are displayed - add p-value display
3. Consider adding p-value calculation helper function

#### Statistical Test Recommendations
```python
from scipy import stats

# For comparing two batches (e.g., control vs test1)
def calculate_pvalue_ttest(batch1_data, batch2_data):
    """
    Calculate p-value using Welch's t-test.
    Appropriate for samples with potentially unequal variances.
    """
    statistic, pvalue = stats.ttest_ind(
        batch1_data,
        batch2_data,
        equal_var=False  # Welch's t-test
    )
    return pvalue

# Alternative: Mann-Whitney U test (non-parametric)
def calculate_pvalue_nonparametric(batch1_data, batch2_data):
    """
    Calculate p-value using Mann-Whitney U test.
    Non-parametric alternative when normality is questionable.
    """
    statistic, pvalue = stats.mannwhitneyu(
        batch1_data,
        batch2_data,
        alternative='two-sided'
    )
    return pvalue
```

#### Display Format
P-values should be shown:
1. As annotations on the plot (e.g., horizontal lines with asterisks)
2. In the statistics text below each batch
3. Format: `p < 0.001***`, `p = 0.023*`, `p = 0.456 (ns)`

#### Comparison Pairs
For each metric, calculate p-values for:
- Control vs Test 1
- Control vs Test 2
- Test 1 vs Test 2

### Implementation Checklist
- [ ] Import scipy.stats in plot_comparison_generic.py
- [ ] Add p-value calculation function
- [ ] Calculate p-values for all three batch pairs
- [ ] Add visual indicators (lines/brackets) between compared groups
- [ ] Include p-values in text annotations
- [ ] Handle edge cases (identical values, insufficient data)
- [ ] Add p-value significance thresholds (*, **, ***)

### P-Value Significance Conventions
- `***` : p < 0.001
- `**`  : p < 0.01
- `*`   : p < 0.05
- `ns`  : p ≥ 0.05 (not significant)

## Phase 3 Statistical Tests (Automatic)

### Overview
Phase 3 pipeline **automatically** generates comprehensive statistical tests in Stage 5:
- **Pairwise t-tests**: Compare each pair of batches (control vs test1, control vs test2, test1 vs test2)
- **ANOVA**: Overall test across all 3 batches

### Tests Performed
For each of the 6 metrics (Cell JSD Mean, Cell Methylation Mean, Gene JSD Mean, Gene Methylation Mean, Gene JSD Std Dev, Gene Methylation Std Dev):

#### Pairwise Tests (3 comparisons × 6 metrics = 18 rows)
1. **Student's t-test** - Assumes equal variance between groups
2. **Welch's t-test** - Does not assume equal variance (more robust)

#### Omnibus Tests (2 tests × 6 metrics = 12 rows)
1. **Fisher's ANOVA** - Standard one-way ANOVA, assumes equal variance
2. **Welch's ANOVA** - Heteroscedastic one-way ANOVA, does not assume equal variance

### Output Files
Generated automatically in `results/tables/`:
- **pvalues.csv** - All pairwise t-test results
- **anova_results.csv** - All ANOVA results with F-statistics and degrees of freedom

### Standalone Usage
The statistical tests can also be run independently:
```bash
cd phase3
python calculate_pvalues.py \
  --cell-csv path/to/batch_comparison_cell.csv \
  --gene-csv path/to/batch_comparison_gene.csv \
  --output-csv path/to/pvalues.csv \
  --anova-csv path/to/anova_results.csv
```

### Important Notes
- **Bonferroni correction**: Should be applied when interpreting p-values (multiply by number of comparisons)
- **Sample sizes**: Number of individuals per batch determined by config (`n_quantiles × cells_per_quantile`) and normalization
- **ANOVA interpretation**: Significant ANOVA indicates at least one batch differs, use pairwise tests to identify which ones

## Configuration System

### Default Configurations
- **Phase 1**: `config_default.yaml`
  - `gene_rate_groups="5:0.004,5:0.005,5:0.006,5:0.007"` (20 genes total)
  - `growth_phase=9` (512 cells)
  - `years=50`, `sites=100`
- **Phase 2**: Default parameters
  - Snapshots: years 23 and 38
  - Quantiles: 5, cells per quantile: 2
  - Growth: 6 years, mix ratio: 82%

### Loading Order
1. Default config loaded first
2. Custom config file overrides defaults
3. Command-line arguments override everything

## Critical Implementation Details

### Phase 1 Specifics
- History uses string year keys: `history["23"]`, not `history[23]`
- Growth phase: exponential growth to 2^N cells
- Homeostasis: divide → cull → methylate cycle
- JSD calculations track divergence from baseline
- Cell ages track years since creation

### Phase 2 Specifics
- **Single-file pipeline** (`phase2_pipeline.py`) - no subprocesses, direct function calls
- **In-memory processing** - data passed between stages without temporary files
- **Outputs directly to Phase 1 directory** (no separate phase2/data folder)
- Snapshot extraction preserves year wrapper: `{"23": {...}}`
- Test 2: quantile-based sampling (sorts by cell_jsd)
- Test 1: random uniform sampling
- Control: pure second snapshot populations
- Normalization: median - 0.5σ threshold applied
- All cells in a PetriDish must have identical gene_rate_groups
- ~10x faster than old multi-script version

### Phase 3 Specifics
- CSV-first workflow: extract once, plot many times
- Auto-detects Phase 1 simulation in parent directory
- Gene aggregation: mean or std across 20 genes per individual
- Batch comparison: Variable number of individuals per batch (determined by Phase 2 config and normalization)
- Timeline plots include percentile bands (5-95%, 25-75%)
- Must use `individual_final` for mixed populations, not `history`

## Development Workflow

### Debugging Data Issues
```python
# Check cell counts in phase3
with smart_open(filepath, 'r') as f:
    data = json.load(f)

# Check available keys
print("Keys:", data.keys())

# Check cell counts
if 'individual_final' in data:
    print(f"individual_final cells: {len(data['individual_final']['cells'])}")
if 'history' in data:
    for year_str in data['history']:
        print(f"history year {year_str}: {len(data['history'][year_str]['cells'])} cells")
```

### Running a Single Phase
```bash
# Test phase1 only
cd phase1
python run_simulation.py --years 5 --sites 20 --growth-phase 3  # Quick test

# Test phase2 with existing simulation (quick test)
cd phase2
python phase2_pipeline.py --simulation ../data/*/simulation.json.gz \
    --n-quantiles 2 --cells-per-quantile 1 --individual-growth-phase 3

# Test phase3 histogram generation
cd phase3
python plot_histogram_original.py --phase2-dir ../data/{dir}/{subdir}/ \
    --year 38 --output test_histogram.png
```

## Project-Specific Rules

### Critical Rules
- **ALWAYS use `python` instead of `python3`** for all commands
- **NO backward compatibility** unless explicitly requested
- Year keys are always strings: `year_str = str(year)`
- Deep copies when necessary to avoid reference issues
- **Dictation note**: "jean" or "gin" means "gene"
- Phase 2 outputs go directly into root data directory alongside Phase 1 simulation

## Common Issues & Solutions

### Phase 3 Issues
- **Wrong cell count in histograms**: Ensure using `individual_final` not `history`
  - `history` only has pre-mixing cells (~68)
  - `individual_final` has post-mixing cells (~210 with default config)
- **Missing timeline plots**: Check Phase 1 simulation exists in parent directory
- **Gene CSV rows**: Should be N × 20 where N = total individuals (varies by config)
- **KeyError in plot_histogram_original.py**: Missing or incorrect data structure

### Phase 2 Issues
- **Gene rate groups mismatch**: All cells must have identical gene_rate_groups
- **Insufficient cells**: Snapshot too small for sampling
- **Memory issues**: Reduce `n_quantiles` or `cells_per_quantile`
- **Output location confusion**: Phase 2 outputs to root data/ directory, not phase2/data
- **Pipeline faster but same results**: New single-file version ~10x faster, output format unchanged

### General Issues
- **ImportError**: Check sys.path additions for phase1
- **No plots**: Install `plotly kaleido` with correct versions
- **Slow simulation**: Reduce sites or years
- **File not found**: Check for .gz extension, use smart_open()

## Dependencies

Required packages (via `pip install -r requirements.txt`):
- `numpy` (>=1.19.0) - Optional but recommended for performance
- `scipy` (>=1.7.0) - Required for statistical analysis
- `pyyaml` (>=6.0) - Configuration file support
- `plotly` (>=5.0.0, <6.0.0) - Visualization
- `kaleido` (0.2.1) - PNG export from plotly

## Performance Tips

- Use compressed files (.json.gz) for large simulations to save ~90% disk space
- Numpy is optional but provides ~10x speedup for large populations
- Use `--skip-plots` in phase3 to quickly generate CSVs without plots
- Use `--skip-gene-comparison` to skip 40 per-gene plots when not needed