# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based methylation simulation project that models DNA methylation patterns over time using stochastic processes. The simulation tracks CpG site methylation across cell populations to study epigenetic drift and aging.

### Primary Workflow: Step 1 + Step23 Unified Pipeline
1. **Step 1**: Base methylation simulation over 100 years (10,000 cells, 1,000 CpG sites)
2. **Step23**: Unified 8-stage pipeline for cell sampling, growth, mixing, and statistical analysis

### Legacy Workflow (deprecated)
- **Step 2 & 3**: Original multi-step approach (now in `legacy/` directory)

## Architecture & Code Structure

### Core Simulation (step1/)
- **`cell.py`**: Core classes and mathematical functions
  - `Cell`: Models individual cell methylation states (N=1000 CpG sites)
  - `History`: Manages M=10,000 cells over T_MAX=100 years
  - Mathematical functions: `KL_div`, `JS_div` for divergence calculations
- **`run_large_simulation.py`**: Main CLI runner with progress tracking
- **`plot_history.py`**: Creates percentile-based plots with Plotly
- **`data/`**: Compressed `.json.gz` simulation outputs
- **`tests/`**: Reproducibility test suite

### Unified Pipeline (step23/)
- **`run_pipeline.py`**: 8-stage orchestrator with intelligent skip logic
- **`pipeline_utils.py`**: Sampling, growth, mixing utilities
- **`pipeline_analysis.py`**: Statistical analysis and visualization
- **`pipeline_checkpoint.py`**: JSON-based progress tracking
- **`data/rate_X.XXXXXX/`**: Rate-specific outputs
  - `pipeline_checkpoint.json`: Progress tracker
  - `snapshots/`: Cached year 50/60 extractions
  - `individuals/`: 30 mutant, 30 control1, 30 control2 files
  - `plots/` and `results/`: Analysis outputs

### Legacy (step2/, step3/)
Deprecated multi-step approach - use step23 instead. Files moved to `legacy/` directory.

## Key Commands

### Step 1: Run Base Simulation
```bash
cd step1
python run_large_simulation.py --rate 0.005  # Creates 10,000 cells, 100 years
python plot_history.py data/simulation_rate_0.005000_m10000_n1000_t100.json.gz --jsd-only
```

### Step23: Run Unified Pipeline V2
```bash
cd step23
# Standard run (30 individuals)
python run_pipeline_v2.py --rate 0.005 --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz

# Quick test (4 individuals)
python run_pipeline_v2.py --rate 0.005 --simulation ../step1/data/*.json.gz \
    --n-quantiles 4 --cells-per-quantile 1 --growth-years 2
```

**Pipeline Options (run_pipeline_v2.py):**
- `--n-quantiles 10`: Number of quantiles for sampling (10=deciles, 4=quartiles)
- `--cells-per-quantile 3`: Cells sampled per quantile (total = n-quantiles × cells-per-quantile)
- `--growth-years 10`: Years of growth (defaults to 10, matching 50→60 gap)
- `--mix-ratio 80`: Percentage of year 60 cells in mix
- `--bins 200`: Histogram bins for JSD plots
- `--seed 42`: Random seed for reproducibility

### 8-Stage Pipeline Process

1. **Extract year 50** → `snapshots/year50_snapshot.json.gz` (cached)
2. **Plot JSD distribution** → `plots/year50_jsd_distribution_200bins.png`
   - Step histogram with statistics: Mean, Median, SD, CV, MAD, 5%, 95%
3. **Create individuals**: 
   - Mutant: Sample from quantiles (default: 3 per decile = 30 total)
   - Control1: Uniform sampling (matches mutant count)
4. **Grow 10 years**: 1→2→4→...→1024 cells (divide then age each year)
5. **Extract year 60** → `snapshots/year60_snapshot.json.gz` (cached)
6. **Mix populations**: Add year 60 cells (80% ratio) → 5120 total
7. **Create control2**: Pure year 60 individuals (5120 cells each)
8. **Analysis**: 
   - Scatter plots with mean/quartile/percentile lines
   - Comprehensive statistics below each group
   - T-tests between all group pairs

### Skip Logic & Checkpoints

The pipeline intelligently skips completed stages by checking:
- File existence AND cell counts (not just first file)
- Handles partial completion (e.g., some grown, some mixed)
- Detects corrupted files gracefully

```bash
# Force recreation patterns
rm -rf data/rate_0.005000/                      # Complete restart
rm -rf data/rate_0.005000/individuals/          # Recreate all individuals
rm -rf data/rate_0.005000/individuals/mutant/   # Recreate mutants only
rm data/rate_0.005000/pipeline_checkpoint.json  # Reset checkpoint

# IMPORTANT: When changing parameters, clean first!
# Parameters requiring cleanup: n-quantiles, cells-per-quantile, mix-ratio, growth-years
```

### Quick Testing
```bash
# Minimal test run (3 individuals, 2 years growth → 4 cells)
python run_pipeline.py --rate 0.005 --simulation ../step1/data/*.json.gz \
    --n-individuals 3 --growth-years 2
```

### Testing & Validation

```bash
# Run reproducibility tests
cd step1/tests
python test_reproducibility.py              # Generate expected results
python test_reproducibility.py --check      # Verify against saved results

# Quick simulation test
cd step1 && python test_small.py  # rate=1%, n=15, m=5, t_max=5
```

### Dependencies
```bash
pip install -r requirements.txt  # plotly and kaleido for visualization only
```
Core simulation requires only Python 3.7+ standard library.

## Key Improvements in Pipeline V2

### Flexible Quantile Sampling
- `--n-quantiles`: Control number of quantiles (4 for quartiles, 10 for deciles, etc.)
- `--cells-per-quantile`: Control cells sampled per quantile
- Total individuals = n-quantiles × cells-per-quantile
- Example: `--n-quantiles 4 --cells-per-quantile 1` = 4 individuals for quick testing

### Improved Skip Logic
- Checks ALL files, not just the first one
- Detects partial states (some grown, some mixed, some corrupted)
- `check_individual_states()` function provides comprehensive file analysis
- Handles corrupted files gracefully without crashing

### Enhanced Visualizations
- Clean scatter plots with jittered points
- Horizontal lines for mean (solid), 25-75% (dashed), 5-95% (dotted)
- Comprehensive statistics: Mean, Median, SD, CV, MAD, 5%, 95%
- Statistics positioned below each group to avoid overlap
- Consistent colors: Mutant (blue), Control1 (orange), Control2 (green)

## Key Implementation Details

### Cell Division & Growth Model
```python
# Division creates identical copy, then both cells age
new_cells = [copy.deepcopy(cell) for cell in cells]  # Division
cells.extend(new_cells)
for cell in cells:
    cell.age()  # Stochastic methylation
```

### Sampling Without Replacement
```python
# Each individual gets unique random subset from year 60
sampled_indices = random.sample(range(len(year60_cells)), n_to_add)
```

### Skip Logic Checks
```python
# Check cell counts to determine stage completion
current_cells = len(data['cells'])
skip_growth = (current_cells == 2**args.growth_years)  # 1024 cells
skip_mixing = (current_cells > 2**args.growth_years)   # >1024 means mixed
```

### Important Constants (step1/cell.py)
- `N = 1000`: CpG sites per cell
- `M = 10_000`: Cells in simulation  
- `T_MAX = 100`: Years to simulate
- `GENE_SIZE = 5`: Sites per gene
- `BASELINE_METHYLATION_DISTRIBUTION`: Reference for JSD calculation

### Data Format
Compressed JSON (`.json.gz`) with structure:
```json
{
  "0": [  // Year
    {
      "cpg_sites": [0, 0, 1, ...],  // N methylation states
      "methylation_proportion": 0.0,
      "methylation_distribution": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      "jsd": 0.0,
      "age": 0,
      "gene_size": 5,
      "rate": 0.005
    }
    // ... M cells
  ]
  // ... T_MAX+1 years
}
```

## Common Development Tasks

### Debugging & Inspection
```bash
# Inspect compressed files
zcat file.json.gz | python -m json.tool | head -20

# Check cell counts in individuals  
python -c "import json, gzip; data=json.load(gzip.open('individual_00.json.gz', 'rt')); print(len(data['cells']))"

# Validate checkpoint status
python -c "import json; print(json.load(open('pipeline_checkpoint.json')))"
```

### Performance Optimization
- Standard simulation (M=10,000, T_MAX=100): ~10 minutes
- Complete pipeline: ~20-30 minutes
- Pipeline intelligently caches snapshots to avoid re-extraction
- Use `--n-individuals 3 --growth-years 2` for quick testing

### Critical Implementation Notes
1. **Always use `.copy()` on lists** to avoid reference sharing bugs
2. **Type consistency**: Use float literals (1.0, 0.0) in distributions
3. **Random seeds**: Set before History creation for reproducibility
4. **File lists**: Always populate file lists before checking (bug fix)
5. **Progress tracking**: Adaptive (frequent early, then every 5 years)

### Running Complete Analysis for a New Rate

```bash
# Step 1: Generate base simulation
cd step1
python run_large_simulation.py --rate 0.007

# Step23: Run unified pipeline
cd ../step23
python run_pipeline.py --rate 0.007 --simulation ../step1/data/simulation_rate_0.007000_m10000_n1000_t100.json.gz
```

### Modifying Simulation Parameters

Edit constants in `step1/cell.py`:
```python
M = 5_000        # Number of cells (default: 10_000)
N = 1000         # CpG sites per cell
GENE_SIZE = 10   # Must evenly divide N
T_MAX = 100      # Years to simulate
```