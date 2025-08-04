# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based methylation simulation project that models DNA methylation patterns over time. The project is organized into three sequential steps:

1. **Step 1**: Base methylation simulation over 100 years
2. **Step 2**: Cell division experiments with lineage tracking (year 50)
3. **Step 3**: Mixed population analysis (year 60)

### Directory Structure
- `step1/`: Base simulation with core classes and runners
  - `cell.py`: Core simulation classes (Cell and History) and mathematical functions
  - `run_large_simulation.py`: Main simulation runner with CLI arguments
  - `main.py`: Legacy multi-rate runner
  - `test_small.py`: Quick test script
  - `plot_history.py`: Visualization script
  - `data/`: Simulation outputs
- `step2/`: Cell division experiments
  - `scripts/`: All executable scripts
  - `data/`: Snapshots, lineages (mutant/control), plots
- `step3/`: Mixed population analysis
  - `scripts/`: All executable scripts
  - `data/`: Snapshots, individuals (mutant/control), plots, results
- `tests/`: Reproducibility test suite
- `requirements.txt`: Python dependencies (plotly and kaleido for visualization only)

## Key Commands

### Step 1: Running Base Simulations

```bash
cd step1
```

**Main simulation with CLI:**
```bash
python run_large_simulation.py --rate 0.005
```
Runs simulation with specified rate, n=1000, m=10,000, t_max=100:
- Progress tracking with time estimates
- Compressed JSON output to `data/`
- Supports fractional rates (e.g., 0.0025 for 0.25%)

**Legacy multi-rate simulation:**
```bash
python main.py
```
Runs simulations with rates 0.002, 0.005, 0.01.

**Quick test:**
```bash
python test_small.py
```
Runs with rate=1%, n=15, m=5, t_max=5 for testing.

**Data Visualization:**
```bash
# Create both JSD and methylation plots
python plot_history.py data/simulation_rate_0.005000_m10000_n1000_t100.json.gz

# Create only JSD plot
python plot_history.py data/simulation_rate_0.005000_m10000_n1000_t100.json.gz --jsd-only

# Custom output name
python plot_history.py data/simulation_rate_0.005000_m10000_n1000_t100.json.gz -o custom_name
```
Creates PNG files with mean line and 5-95/25-75 percentile bands.

### Step 2: Cell Division Experiments

```bash
cd step2/scripts
```

**Extract year 50 snapshot:**
```bash
python extract_snapshot.py --simulation ../../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz --year 50
```
Extracts all cells at specified year.

**Create lineages (both mutant and control):**
```bash
python create_lineages.py --type both
```
- **Mutant lineages**: Samples 3 cells from each JSD decile (30 total)
- **Control lineages**: Samples 30 cells uniformly from population
- Ages each cell for 10 years with division
- Creates 1,024 cells per lineage
- Outputs to `data/lineages/mutant/` and `data/lineages/control/`

**Plot JSD distribution:**
```bash
python plot_jsd_distribution.py ../data/snapshots/year50_snapshot.json.gz 200
```

### Step 3: Mixed Population Analysis

```bash
cd step3/scripts
```

**Extract year 60 from original simulation:**
```bash
python extract_year60_original.py --simulation ../../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz
```

**Create mixed individuals (both types):**
```bash
python create_individuals.py --type both
```
- Uses lineages from step2 (both mutant and control)
- Each individual: 80% original year 60 cells + 20% lineage cells
- Creates 30 individuals per type
- Outputs to `data/individuals/mutant/` and `data/individuals/control/`

**Plot comparison:**
```bash
python plot_distributions.py
```
- Compares mean JSD between mutant and control populations
- Saves plots and statistics to `data/plots/` and `data/results/`

### Testing

**Run reproducibility tests:**
```bash
cd tests
python test_reproducibility.py
```
Verifies that optimizations haven't changed simulation outputs using fixed random seeds.

**Test step pipelines:**
```bash
# Test step 2
cd step2/scripts
python test_pipeline.py

# Test step 3
cd step3/scripts
python test_pipeline.py
```
Runs complete pipelines with small test data.

### Dependencies
```bash
pip install -r requirements.txt
```
Installs:
- plotly>=5.0.0,<6.0.0 (for visualization)
- kaleido==0.2.1 (for PNG export)

Note: Core simulation requires only Python 3.7+ standard library.

## Architecture & Code Structure

### Simulation Architecture (step1/cell.py)

1. **Cell Class**: Models individual cell methylation states
   - Tracks methylation of N CpG sites (stored in `cpg_sites` array) grouped into genes of GENE_SIZE
   - Implements aging mechanism with stochastic methylation
   - Calculates methylation distribution and Jensen-Shannon Divergence (JSD)
   - Provides `to_dict()` method for serialization (MUST use .copy() on lists to avoid reference issues)
   - All distribution lists use float types for consistency

2. **History Class**: Manages cell population and records simulation over time
   - Creates and manages M cells directly
   - Coordinates aging across all cells
   - Stores complete cell state (including cpg_sites, distributions, etc.) at each time point
   - Handles data persistence to JSON in the `history/` directory
   - Custom JSON formatting keeps cpg_sites arrays on single lines for compactness
   - Provides methods: `run()`, `save()`, `load()` (static method)

3. **Key Mathematical Functions**:
   - `KL_div`: Kullback-Leibler divergence calculation
   - `JS_div`: Jensen-Shannon divergence (symmetrized KL divergence)

### Step 1 Scripts

**run_large_simulation.py**: Main simulation runner
- Command-line interface with argparse
- Supports custom methylation rates with 6 decimal precision
- Progress tracking with time estimates
- Compressed output to `data/` directory

**main.py**: Legacy multi-rate runner
- Runs simulations with rates 0.002, 0.005, 0.01
- Saves to `data/` directory

**plot_history.py**: Visualization tool
- Creates percentile-based plots (5-95 and 25-75 bands)
- Exports to PNG using kaleido backend
- Supports --jsd-only and --methylation-only flags

## Important Constants

Defined in `step1/cell.py`:
- `N = 1000`: Number of methylation sites per cell
- `M = 10_000`: Number of cells in simulation
- `T_MAX = 100`: Maximum age in years
- `GENE_SIZE = 5`: Size of methylation genes
- `BASELINE_METHYLATION_DISTRIBUTION`: Reference distribution for JSD calculations

## Data Format

Simulation output structure:
```json
{
  "0": [  // Time point (year)
    {
      "cpg_sites": [0, 0, 0, ...],  // Array of N methylation states
      "methylation_proportion": 0.0,
      "methylation_distribution": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      "jsd": 0.0,
      "age": 0,
      "gene_size": 5,
      "rate": 0.01
    },
    // ... M cells total
  ],
  // ... T_MAX+1 time points
}
```

## Code Style & Best Practices

1. **Type Hints**: All functions and methods have type annotations
2. **Docstrings**: Classes and key methods have descriptive docstrings
3. **Naming**: Use descriptive names (e.g., `cpg_sites` instead of `strain`)
4. **File Organization**: Each step has `scripts/` and `data/` subdirectories
5. **Output Organization**: All outputs go to appropriate `data/` subdirectories

## Common Tasks

### Adding a New Methylation Rate
Use command-line argument:
```bash
python step1/run_large_simulation.py --rate 0.015
```

### Changing Population Size
Edit constants in `step1/cell.py`:
```python
M = 5_000  # Reduced from 10_000 for faster testing
```

### Modifying Gene Size
Change `GENE_SIZE` in `step1/cell.py` and ensure `N` is divisible by it:
```python
GENE_SIZE = 10  # Changed from 5
N = 1000  # Still divisible
```

### Creating Custom Simulation Scripts
Copy structure from `step1/test_small.py`:
1. Import from `step1.cell`
2. Set custom parameters
3. Create History instance
4. Run simulation
5. Save to `step1/data/` with descriptive filename

## Known Issues & Solutions

1. **Reference Sharing Bug**: Always use `.copy()` when storing lists in dictionaries to avoid mutations affecting historical data
2. **Type Consistency**: Use float literals (1.0, 0.0) for distribution lists to match type hints
3. **Performance**: For large simulations (M=10,000, T_MAX=100), expect runtime of several minutes
4. **Memory Usage**: Large simulations can use significant memory; compressed output helps with disk storage

## Performance Notes

- Standard simulation (M=10,000, T_MAX=100): ~10 minutes
- Progress tracking: Every year for first 10 years, then every 5 years
- Compressed output reduces file size by ~90%
- Plotting large files may take 20-30 seconds
- Complete 3-step pipeline: ~20-30 minutes

## Development Workflow

### Before Making Changes
1. Run reproducibility tests to establish baseline:
   ```bash
   cd tests && python test_reproducibility.py --check
   ```
2. Create a test script based on `test_small.py` for quick iteration

### After Making Changes
1. Run your test script to verify basic functionality
2. Run reproducibility tests to check for unintended changes:
   ```bash
   cd tests && python test_reproducibility.py --check
   ```
3. If changes are intentional, update expected results:
   ```bash
   cd tests && python test_reproducibility.py  # Regenerates expected.json
   ```
4. Run a full simulation to verify performance hasn't degraded

### Debugging Tips
- Use `step1/test_small.py` as a template for minimal reproducible examples
- Check `cell.to_dict()` output for serialization issues (remember `.copy()` for lists)
- For performance issues, profile `step1/cell.py` methods
- Test pipelines available in `step2/scripts/test_pipeline.py` and `step3/scripts/test_pipeline.py`
- Compressed files can be inspected with: `zcat file.json.gz | python -m json.tool | head`

## Recent Changes

- **Major restructuring**: Organized project into 3 sequential steps
- **Step 1**: Moved all simulation files to dedicated directory with `data/` subdirectory
- **Step 2**: Unified lineage creation for both mutant (decile-based) and control (uniform)
- **Step 3**: Updated to use step2's control lineages instead of creating its own
- **CLI improvements**: Added command-line arguments to all major scripts
- **Precision update**: Increased rate precision from 3 to 6 decimal places
- **Directory structure**: Consistent `scripts/` and `data/` organization across steps
- **Test pipelines**: Added test scripts for steps 2 and 3
- **Gitignore**: Updated for new directory structure