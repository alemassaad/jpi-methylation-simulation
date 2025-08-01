# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based methylation simulation project that models DNA methylation patterns over time. The project consists of:
- `cell.py`: Core simulation classes (Cell and History) and mathematical functions
- `main.py`: Standard simulation runner for multiple methylation rates
- `run_large_fast_save.py`: Optimized script for large-scale simulations with compressed output
- `test_small.py`: Quick test script with small parameters
- `plot_history.py`: Visualization script creating PNG plots with percentile bands
- Simulation output files: Compressed JSON files (`.json.gz`) saved in `history/` directory

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
- Progress tracking every year/5 years
- Fast compressed JSON output (.json.gz)
- ~10 minutes total runtime

**Quick test:**
```bash
python test_small.py
```
Runs with rate=1%, n=15, m=5, t_max=5 for testing.

### Data Visualization
```bash
# Create both JSD and methylation plots
python plot_history.py baseline_rate_0.5%_m10000_n1000_t100.json.gz

# Create only JSD plot
python plot_history.py baseline_rate_0.5%_m10000_n1000_t100.json.gz --jsd-only

# Create only methylation plot  
python plot_history.py baseline_rate_0.5%_m10000_n1000_t100.json.gz --methylation-only
```
Creates separate PNG files with mean line and 5-95/25-75 percentile bands.

### Dependencies
```bash
pip install -r requirements.txt
```
Installs:
- plotly>=5.0.0,<6.0.0 (for visualization)
- kaleido==0.2.1 (for PNG export)

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

3. **Key Mathematical Functions**:
   - `KL_div`: Kullback-Leibler divergence calculation
   - `JS_div`: Jensen-Shannon divergence (symmetrized KL divergence)

### Simulation Execution (main.py)

- Imports necessary classes and constants from cell.py
- Runs simulations with different methylation rates (0.002, 0.005, 0.01)
- Extracts and displays final statistics (average methylation proportion and JSD)
- Saves complete history data for each rate with percentage-based filenames (e.g., `simulation_history_0.200%.json`)

## Important Constants

- `N = 1000`: Number of methylation sites per cell
- `M = 10_000`: Number of cells in simulation
- `T_MAX = 100`: Maximum age in years
- `GENE_SIZE = 5`: Size of methylation genes
- `BASELINE_METHYLATION_DISTRIBUTION`: Reference distribution for JSD calculations (list of floats)

## Code Style & Best Practices

1. **Type Hints**: All functions and methods have type annotations
2. **Docstrings**: Classes and key methods have descriptive docstrings
3. **Naming**: Use descriptive names (e.g., `cpg_sites` instead of `strain`)
4. **File Organization**: Separate concerns - `cell.py` for library code, `main.py` for execution
5. **Output Organization**: All simulation outputs go to `history/` directory

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

## Known Issues & Solutions

1. **Reference Sharing Bug**: Always use `.copy()` when storing lists in dictionaries to avoid mutations affecting historical data
2. **Type Consistency**: Use float literals (1.0, 0.0) for distribution lists to match type hints
3. **Performance**: For large simulations (M=10,000, T_MAX=100), expect runtime of several minutes

## Recent Changes

- Renamed `main.py` to `cell.py` and created new `main.py` for separation of concerns
- Changed `n_years` parameter to `t_max` throughout codebase
- Fixed typo: `distrbution` → `distribution`
- Added comprehensive type hints and docstrings
- Moved all output files to `history/` directory
- Fixed methylation distribution to use floats for type consistency