# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based methylation simulation project that models DNA methylation patterns over time using stochastic processes. The simulation tracks CpG site methylation across cell populations to study epigenetic drift and aging.

### Simulation Options

#### Step 1-Prime: Biologically Realistic (NEW - Recommended)
- Starts with 1 cell, grows exponentially to target population
- Growth phase: Population doubles yearly (1→2→4→8...2^growth_phase)
- Steady state: Maintains population via division + random culling
- Time-based state transitions (no population-based flip-flopping)
- Configurable target via `--growth-phase` parameter (default: 13 → 8192 cells)

#### Step 1: Traditional Approach
- Starts with 10,000 unmethylated cells
- All cells age in parallel without division
- Fixed population throughout simulation

### Primary Workflow: Step 1/1-Prime + Step23 Unified Pipeline
1. **Step 1 or 1-Prime**: Base methylation simulation
2. **Step23**: Unified 8-stage pipeline for cell sampling, growth, mixing, and statistical analysis

### Legacy Workflow (deprecated)
- **Step 2 & 3**: Original multi-step approach (now in `legacy/` directory)

## Architecture & Code Structure

### Biologically Realistic Simulation (step1-prime/)
- **`cell.py`**: Core classes with biological realism
  - `Cell`: Individual cell with methylation tracking
    - `methylate()`: Apply stochastic methylation
    - `create_daughter_cell()`: Mitosis (creates identical copy)
    - `to_dict()`: Serialization for output
  - `PetriDish`: Population manager with growth/homeostasis
    - `divide_cells()`: All cells undergo mitosis
    - `methylate_cells()`: Apply methylation to all cells
    - `random_cull_cells()`: ~50% survival for homeostasis
    - `age_cells_growth_phase()`: Growth phase logic (divide → methylate)
    - `age_cells_steady_state()`: Steady state logic (divide → cull → methylate)
    - `simulate_year()`: Time-based state management
  - Mathematical functions: `KL_div`, `JS_div` for divergence
- **`run_simulation.py`**: CLI with comprehensive parameters
  - `--growth-phase`: Control target population (1-20)
  - `--rate`: Methylation rate per site per year
  - `--years`: Total simulation time
  - `--seed`: Reproducibility (default: 42)
- **`tests/`**: Comprehensive test coverage
  - `test_comprehensive.py`: 26 tests covering all functionality
  - `test_edge_cases.py`: 10 edge case tests
- **`data/`**: Output format: `simulation_rate_X_gY_mZ_nN_tT_seedS.json.gz`

### Traditional Simulation (step1/)
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

### Step 1-Prime: Biologically Realistic Simulation (Recommended)
```bash
cd step1-prime

# Standard run (8192 cells, 100 years)
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Quick test (16 cells, 20 years)
python run_simulation.py --rate 0.01 --years 20 --growth-phase 4

# Different population sizes
python run_simulation.py --growth-phase 10  # 1024 cells
python run_simulation.py --growth-phase 13  # 8192 cells (default)
python run_simulation.py --growth-phase 15  # 32768 cells

# Run tests
python tests/test_comprehensive.py  # All 26 tests
python tests/test_edge_cases.py     # Edge cases
```

### Step 1: Traditional Simulation
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

#### Step 1-Prime Tests
```bash
cd step1-prime

# Comprehensive test suite (36 total tests)
python tests/test_comprehensive.py  # 26 tests:
  # - Growth phase tests (population doubling, division, methylation)
  # - Steady state tests (population maintenance, culling, transitions)
  # - Methylation mechanics (initial state, accumulation, distribution)
  # - Output format tests (JSON structure, cell format, file I/O)
  # - Reproducibility tests (same/different seeds)
  # - Edge cases (small populations, gene validation)
  # - Compatibility tests (step23 format, statistics)
  # - Growth phase parameter tests

python tests/test_edge_cases.py  # 10 tests:
  # - Single year simulation
  # - Ending during growth phase
  # - Immediate steady state (growth_phase=1)
  # - Extreme culling survival
  # - 100% and 0% methylation rates
  # - Minimum gene configuration
  # - Maximum growth_phase (20)
  # - Filename generation
  # - Population variability
```

#### Step 1 Tests
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

### Important Constants

#### Step 1-Prime (step1-prime/cell.py)
- `N = 1000`: CpG sites per cell
- `RATE = 0.005`: Default methylation rate (0.5% per site per year)
- `GENE_SIZE = 5`: Sites per gene
- `T_MAX = 100`: Default years to simulate
- `DEFAULT_GROWTH_PHASE = 13`: Default growth phase (2^13 = 8192 cells)
- `BASELINE_METHYLATION_DISTRIBUTION`: Reference for JSD calculation

#### Step 1 (step1/cell.py)
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

#### Step 1-Prime Specific
1. **Time-based state transitions**: Use `year <= growth_phase` not population checks
2. **Growth phase is deterministic**: Exactly 2^year cells each year
3. **Steady state is stochastic**: Population varies around target
4. **Cell division**: Use `create_daughter_cell()` for deep copies
5. **Filename format**: Includes growth phase (`_g13_`) and actual final population
6. **Seed handling**: Default 42, `_noseed` suffix when seed=None
7. **Parameter validation**: growth_phase must be 1-20 (max 2^20 = 1M cells)

#### General Notes
1. **Always use `.copy()` on lists** to avoid reference sharing bugs
2. **Type consistency**: Use float literals (1.0, 0.0) in distributions
3. **Random seeds**: Set before History/PetriDish creation for reproducibility
4. **File lists**: Always populate file lists before checking (bug fix)
5. **Progress tracking**: Adaptive (frequent early, then every 5 years)

### Running Complete Analysis for a New Rate

#### Using Step 1-Prime (Biologically Realistic)
```bash
# Generate simulation with custom rate
cd step1-prime
python run_simulation.py --rate 0.007 --years 100 --growth-phase 13 --seed 42

# Run unified pipeline
cd ../step23
python run_pipeline_v2.py --rate 0.007 --simulation ../step1-prime/data/simulation_rate_0.007000_g13_m*_n1000_t100_seed42.json.gz
```

#### Using Step 1 (Traditional)
```bash
# Generate base simulation
cd step1
python run_large_simulation.py --rate 0.007

# Run unified pipeline
cd ../step23
python run_pipeline_v2.py --rate 0.007 --simulation ../step1/data/simulation_rate_0.007000_m10000_n1000_t100.json.gz
```

### Modifying Simulation Parameters

#### Step 1-Prime Parameters
Edit constants in `step1-prime/cell.py` or use CLI:
```python
# In cell.py:
N = 1000                    # CpG sites per cell
RATE = 0.005               # Default methylation rate
GENE_SIZE = 5              # Must evenly divide N
T_MAX = 100                # Default years to simulate
DEFAULT_GROWTH_PHASE = 13  # Default growth phase (2^13 = 8192)

# Or use CLI arguments:
python run_simulation.py \
  --sites 2000 \
  --rate 0.01 \
  --gene-size 10 \
  --years 50 \
  --growth-phase 10  # 1024 cells
```

#### Step 1 Parameters
Edit constants in `step1/cell.py`:
```python
M = 10_000       # Number of cells
N = 1000         # CpG sites per cell
GENE_SIZE = 5    # Must evenly divide N
T_MAX = 100      # Years to simulate
RATE = 0.01      # Default methylation rate
```

## Step 1 vs Step 1-Prime: When to Use Each

### Use Step 1-Prime (Biologically Realistic) When:
- **Modeling real tissue dynamics**: Need to simulate actual cell division and population growth
- **Studying lineage effects**: Want to track how methylation patterns propagate through cell divisions
- **Variable population sizes**: Need to test different population sizes easily (via growth_phase)
- **Homeostasis modeling**: Want to study steady-state dynamics with cell turnover
- **More realistic variance**: Population bottlenecks and expansions affect methylation diversity

### Use Step 1 (Traditional) When:
- **Backward compatibility**: Need to reproduce previous results or compare with existing data
- **Simple aging model**: Only interested in methylation accumulation without division
- **Fixed population required**: Study design requires exactly 10,000 cells throughout
- **Computational simplicity**: Faster runtime without division/culling overhead
- **Pure methylation kinetics**: Want to isolate methylation dynamics from population dynamics

### Key Differences Summary

| Aspect | Step 1-Prime | Step 1 |
|--------|-------------|---------|
| Starting cells | 1 | 10,000 |
| Population dynamics | Growth → Homeostasis | Fixed |
| Cell division | Yes (mitosis) | No |
| Population control | Random culling | N/A |
| State transitions | Time-based | N/A |
| Methylation timing | After division | Each year |
| Population size | Configurable (2^growth_phase) | Fixed (M) |
| Biological realism | High | Low |
| Computational cost | Higher | Lower |
| Output filename | Includes growth_phase | Standard |

### Migration Guide

To convert Step 1 results to Step 1-Prime equivalent:
1. Choose growth_phase where 2^growth_phase ≈ original M value
   - M=10,000 → growth_phase=13 (8,192 cells)
   - M=5,000 → growth_phase=12 (4,096 cells)
2. Account for lineage effects in analysis (cells are related in Step 1-Prime)
3. Consider that Step 1-Prime steady state has ongoing turnover

### Performance Considerations

Step 1-Prime computational requirements:
- **Memory**: O(2^growth_phase × N) for cell storage
- **Time per year**: 
  - Growth phase: O(2^year) for division + methylation
  - Steady state: O(2^growth_phase) for division + culling + methylation
- **Typical runtime**: 
  - Growth_phase=13, 100 years: ~30-60 seconds
  - Growth_phase=10, 100 years: ~10-20 seconds

Step 1 computational requirements:
- **Memory**: O(M × N) constant
- **Time per year**: O(M) for methylation only
- **Typical runtime**: M=10,000, 100 years: ~10 minutes

### Troubleshooting Step 1-Prime

Common issues and solutions:

1. **Population dies out**: Shouldn't happen (always ≥1 survivor), but if it does, check seed
2. **Memory issues**: Reduce growth_phase (each unit doubles memory requirement)
3. **Unexpected population size**: Remember it's stochastic in steady state
4. **Test failures**: Run with small parameters first (growth_phase=3, years=10)
5. **Filename parsing**: Use glob patterns for varying population (*_m*_*)

### Development Guidelines for Step 1-Prime

When modifying Step 1-Prime:

1. **Maintain time-based logic**: Never use population size for state decisions
2. **Preserve deep copying**: Always use create_daughter_cell() for mitosis
3. **Test edge cases**: Single year, immediate steady state, maximum growth_phase
4. **Document stochastic elements**: Clearly mark what's deterministic vs random
5. **Keep backward compatibility**: Don't change default parameters without good reason
EOF < /dev/null