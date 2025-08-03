# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based methylation simulation project that models DNA methylation patterns over time. The project consists of:
- `cell.py`: Core simulation classes (Cell and History) and mathematical functions
- `main.py`: Standard simulation runner for multiple methylation rates
- `run_large_fast_save.py`: Optimized script for large-scale simulations with compressed output
- `test_small.py`: Quick test script with small parameters
- `plot_history.py`: Visualization script creating PNG plots with percentile bands
- `tests/`: Reproducibility test suite with fixed random seeds
- `step2/`: Cell division experiments and analysis scripts
- `history/`: Directory containing simulation output files (compressed JSON `.json.gz`)
- `plots/`: Directory containing PNG visualizations
- `requirements.txt`: Python dependencies (plotly and kaleido for visualization only)

## Key Commands

### Running Simulations

**Standard simulation (multiple rates):**
```bash
python main.py
```
Runs simulations with rates 0.002, 0.005, 0.01 and saves uncompressed JSON files.

**Large-scale optimized simulation:**
```bash
python run_large_fast_save.py
```
Runs simulation with rate=0.5%, n=1000, m=10,000, t_max=100 with:
- Progress tracking every year/5 years with time estimates
- Fast compressed JSON output (.json.gz)
- ~10 minutes total runtime

**Quick test:**
```bash
python test_small.py
```
Runs with rate=1%, n=15, m=5, t_max=5 for testing.

### Data Visualization
```bash
# Create both JSD and methylation plots (default)
python plot_history.py simulation_rate_0.005_m10000_n1000_t100.json.gz

# Create only JSD plot
python plot_history.py simulation_rate_0.005_m10000_n1000_t100.json.gz --jsd-only

# Create only methylation plot  
python plot_history.py simulation_rate_0.005_m10000_n1000_t100.json.gz --methylation-only

# Custom output name (creates custom_name_jsd.png and custom_name_methylation.png)
python plot_history.py simulation_rate_0.005_m10000_n1000_t100.json.gz -o custom_name
```
Creates separate PNG files with mean line and 5-95/25-75 percentile bands.
All plots are automatically saved to the `plots/` directory.

### Step 2: Cell Division Analysis

**Extract a snapshot from existing simulation:**
```bash
cd step2
python extract_snapshot.py ../history/simulation_rate_0.005_m10000_n1000_t100.json.gz 50 year50_snapshot.json.gz
```
Extracts all cells at year 50 and saves them as a snapshot for further analysis.

**Plot JSD distribution:**
```bash
python plot_jsd_distribution.py year50_snapshot.json.gz 200  # 200 bins
```
Creates a step histogram showing JSD distribution across all cells.

**NEW: Run cell division with separate lineage tracking:**
```bash
python sample_divide_age_lineages.py
```
- Loads year 50 snapshot
- Samples 3 cells from each JSD decile (30 total)
- Ages each cell SEPARATELY for 10 years with division
- Creates 30 lineage files in `lineages/` directory
- Each file contains 1,024 cells from one lineage

**Plot individual lineage distributions:**
```bash
# Plot all lineages
for file in lineages/lineage_*.json.gz; do
    python plot_jsd_distribution.py "$file" 100
done
```

### Step 3: Mixed Population Analysis

**Extract year 60 from original simulation:**
```bash
cd step3
python extract_year60_original.py
```

**Create 30 mixed individuals:**
```bash
python create_individuals.py
```
- Loads 30 lineage files from step2
- For each lineage: mixes 1,024 lineage cells with 4,096 original cells (80-20 ratio)
- Creates 30 "individuals" in `individuals/` directory
- Each individual has 5,120 cells total

**Create 30 control individuals:**
```bash
python create_control_individuals.py
```
- Samples 5,120 cells from original year 60 (without replacement per individual)
- Creates 30 control individuals in `control_individuals/` directory

**Plot comparison:**
```bash
python plot_distributions.py
```
- Compares mean JSD distributions between mixed and control individuals
- Creates visualization showing both distributions with statistics

### Testing

**Run reproducibility tests:**
```bash
cd tests
python test_reproducibility.py
```
Verifies that optimizations haven't changed simulation outputs using fixed random seeds.

**Check against expected results:**
```bash
cd tests
python test_reproducibility.py --check
```
Compares current simulation outputs against saved expected values to detect unintended changes.

### Dependencies
```bash
pip install -r requirements.txt
```
Installs:
- plotly>=5.0.0,<6.0.0 (for visualization)
- kaleido==0.2.1 (for PNG export)

Note: Core simulation requires only Python 3.7+ standard library.

## Architecture & Code Structure

### Simulation Architecture (cell.py)

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

### Simulation Execution Scripts

**main.py**: Standard multi-rate runner
- Imports necessary classes and constants from cell.py
- Runs simulations with different methylation rates (0.002, 0.005, 0.01)
- Extracts and displays final statistics (average methylation proportion and JSD)
- Saves complete history data for each rate with decimal filenames (e.g., `simulation_rate_0.002.json`)

**run_large_fast_save.py**: Optimized for large simulations
- Hardcoded for rate=0.005 (0.5%)
- Uses gzip compression for output
- Includes progress tracking with time estimates
- Outputs to `simulation_rate_0.005_m10000_n1000_t100.json.gz`

**plot_history.py**: Visualization tool
- Supports both compressed (.json.gz) and uncompressed (.json) files
- Creates percentile-based plots (5-95 and 25-75 bands)
- Exports to PNG using kaleido backend
- Annotates plots with final statistics

## Important Constants

- `N = 1000`: Number of methylation sites per cell
- `M = 10_000`: Number of cells in simulation
- `T_MAX = 100`: Maximum age in years
- `GENE_SIZE = 5`: Size of methylation genes
- `BASELINE_METHYLATION_DISTRIBUTION`: Reference distribution for JSD calculations (list of floats)

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
4. **File Organization**: Separate concerns - `cell.py` for library code, execution scripts for specific use cases
5. **Output Organization**: All simulation outputs go to `history/` directory, plots to `plots/`

## Common Tasks

### Adding a New Methylation Rate
Edit `main.py` and add to the rates list:
```python
rates = [0.002, 0.005, 0.01, 0.015]  # Added 0.015
```

### Changing Population Size
Edit constants in `cell.py`:
```python
M = 5_000  # Reduced from 10_000 for faster testing
```

### Modifying Gene Size
Change `GENE_SIZE` in `cell.py` and ensure `N` is divisible by it:
```python
GENE_SIZE = 10  # Changed from 5
N = 1000  # Still divisible
```

### Creating Custom Simulation Scripts
Copy structure from `test_small.py` or `run_large_fast_save.py`:
1. Import from cell.py
2. Set custom parameters
3. Create History instance
4. Run simulation
5. Save with descriptive filename

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
- Use `test_small.py` as a template for minimal reproducible examples
- Check `cell.to_dict()` output for serialization issues (remember `.copy()` for lists)
- For performance issues, profile `cell.py` methods (methylation and distribution calculations)
- Compressed files can be inspected with: `zcat file.json.gz | python -m json.tool | head`

## Recent Changes

- Renamed `main.py` to `cell.py` and created new `main.py` for separation of concerns
- Changed `n_years` parameter to `t_max` throughout codebase
- Fixed typo: `distrbution` → `distribution`
- Added comprehensive type hints and docstrings
- Moved all output files to `history/` directory
- Fixed methylation distribution to use floats for type consistency
- Added `plots/` directory for visualization outputs
- Optimized `cell.py` methods for ~2x speedup:
  - Replaced `statistics.mean()` with counter-based calculation
  - Added early exit for fully methylated cells
  - Optimized `compute_methylation_distribution()` to avoid sublists
- Changed filenames from percentage format (0.5%) to decimal format (0.005)
- Added `tests/` directory with reproducibility test suite
- Created `step2/` for cell division experiments and analysis