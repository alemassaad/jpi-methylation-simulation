# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a Python-based methylation simulation project that models DNA methylation patterns over time. The project supports two approaches:

### Recommended: Unified Pipeline
1. **Step 1**: Base methylation simulation over 100 years
2. **Step23**: Unified cell sampling, growth, mixing, and analysis pipeline

### Legacy: Three-Step Pipeline
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
  - `run_simulation.py`: Alternative simulation runner
  - `test_simulation.py`: Testing script
  - `data/`: Simulation outputs (compressed .json.gz files)
  - `tests/`: Reproducibility test suite
    - `test_reproducibility.py`: Verifies simulation determinism
    - `test_reproducibility_expected.json`: Expected results
    - `test_reproducibility_README.md`: Test documentation
- `step23/`: **RECOMMENDED** Unified pipeline
  - `run_pipeline.py`: Main pipeline orchestrator with 8-stage process
  - `pipeline_utils.py`: Core utilities for sampling, growth, and mixing
  - `pipeline_analysis.py`: Visualization and statistical analysis
  - `pipeline_checkpoint.py`: Checkpoint tracking system for skip logic
  - `__init__.py`: Package initialization
  - `PIPELINE_DESIGN.md`: Original design document
  - `data/`: All outputs organized by rate
    - `rate_X.XXXXXX/`: Rate-specific directory
      - `pipeline_checkpoint.json`: Progress tracking and skip logic
      - `snapshots/`: Year 50 and 60 snapshots (cached)
      - `individuals/`: All individual files
        - `mutant/`: 30 mutant individuals (decile-sampled)
        - `control1/`: 30 control individuals (uniform-sampled, grown)
        - `control2/`: 30 control individuals (pure year 60)
      - `plots/`: JSD distributions and comparisons
      - `results/`: Statistical analysis and metadata
- `step2/`: **LEGACY** Cell division experiments
  - `scripts/`: All executable scripts
    - `extract_snapshot.py`: Extract year from simulation
    - `create_lineages.py`: Create mutant/control lineages
    - `plot_jsd_distribution.py`: Visualize JSD distributions
    - `batch_processor.py`: Batch process multiple rates
    - `config/`: Configuration files for batch processing
    - `utils/`: Utility modules for file handling, logging, process management
  - `data/`: Snapshots, lineages (mutant/control), plots
  - `test_sample_divide_age_lineages.py`: Test lineage creation
- `step3/`: **LEGACY** Mixed population analysis
  - `scripts/`: All executable scripts
    - `extract_year60_original.py`: Extract year 60 from simulation
    - `create_individuals.py`: Create mixed populations
    - `plot_distributions.py`: Compare JSD distributions
    - `batch_processor_step3.py`: Batch process multiple rates
    - `test_pipeline.py`: Test complete pipeline
    - `config/`: Configuration files for batch processing
    - `utils/`: Extended utilities for step 3 processing
  - `data/`: Snapshots, individuals (mutant/control), plots, results
  - Legacy scripts for individual creation and testing
- `requirements.txt`: Python dependencies (plotly and kaleido for visualization only)

## Step23: Unified Pipeline (RECOMMENDED - Use This!)

### Quick Start
```bash
cd step23
python run_pipeline.py --rate 0.005 --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz
```

### Detailed Pipeline Documentation

#### Architecture
The pipeline consists of 8 sequential stages, each with intelligent skip logic:
1. Extract year 50 snapshot
2. Plot year 50 JSD distribution  
3. Create initial individuals (mutant & control1)
4. Grow individuals for 10 years
5. Extract year 60 snapshot
6. Mix populations with year 60 cells
7. Create control2 individuals
8. Analyze and compare all groups

#### Key Features
- **Checkpoint System**: Tracks progress in `pipeline_checkpoint.json`
- **Smart Skip Logic**: Automatically detects completed work
- **In-place Operations**: No intermediate lineage files
- **Efficient Caching**: Snapshots cached for reuse
- **Clean Visualizations**: Simple scatter plots, enhanced statistics

#### Command-Line Options
```bash
python run_pipeline.py \
    --rate 0.005 \                      # Methylation rate (MUST match simulation)
    --simulation path/to/sim.json.gz \  # Path to step1 simulation file
    --n-individuals 30 \                 # Number of individuals per group (default: 30)
    --growth-years 10 \                  # Years of growth 50→60 (default: 10)
    --mix-ratio 80 \                     # % of year 60 cells in mix (default: 80)
    --bins 200 \                         # Histogram bins for plots (default: 200)
    --seed 42 \                          # Random seed (default: 42)
    --output-dir data                    # Output directory (default: data)
```

#### Stage Details

**Stage 1: Extract Year 50 Snapshot**
- Extracts all cells at year 50
- Saves to `snapshots/year50_snapshot.json.gz`
- Skip: If file exists AND checkpoint confirms

**Stage 2: Plot JSD Distribution**
- Step histogram with filled area
- Statistics: Mean, Median, SD, CV, MAD, 5%, 95%
- X-axis: "JSD Score"
- Rate shown as percentage
- Bins in filename only

**Stage 3: Create Individuals**
- Mutant: 3 cells from each JSD decile (30 total)
- Control1: 30 cells uniformly sampled
- Skip: If all 30 files exist per group

**Stage 4: Grow Individuals**
- Model: Divide then age each year
- Growth: 1→2→4→8...→1024 cells over 10 years
- In-place file updates
- Skip: If already at 1024 cells

**Stage 5: Extract Year 60**
- Same as Stage 1 for year 60

**Stage 6: Mix Populations**
- Samples year 60 cells WITHOUT replacement
- Default: 80% year 60, 20% grown = 5120 cells
- Skip: If cells > 1024 (already mixed)

**Stage 7: Create Control2**
- 30 individuals of pure year 60 cells
- 5120 cells each (matches mixed size)
- Skip: If all 30 files exist

**Stage 8: Analysis**
- Calculates mean JSD per individual
- Simple scatter plots with jitter
- T-tests between groups
- Saves to `results/` and `plots/`

#### Skip Logic Implementation

**File-based checks:**
```python
# Check if individuals exist
mutant_files = sorted(glob.glob(os.path.join(mutant_dir, "*.json.gz")))
skip_mutant = len(mutant_files) == args.n_individuals

# Check if grown
with gzip.open(mutant_files[0], 'rt') as f:
    data = json.load(f)
    current_cells = len(data['cells'])
    skip_growth = (current_cells == 2**args.growth_years)

# Check if mixed
skip_mixing = (current_cells > 2**args.growth_years)
```

**Checkpoint tracking:**
```python
checkpoint = PipelineCheckpoint("pipeline_checkpoint.json")
checkpoint.mark_stage_complete("extract_year50", {"cells": 10000})
if checkpoint.is_stage_complete("extract_year50"):
    print("Skipping year 50 extraction")
```

#### Forcing Recreation
```bash
# Complete restart
rm -rf data/rate_0.005000/

# Recreate individuals only
rm -rf data/rate_0.005000/individuals/

# Recreate specific group
rm -rf data/rate_0.005000/individuals/mutant/

# Clear checkpoint to reprocess with existing files
rm data/rate_0.005000/pipeline_checkpoint.json
```

#### Common Issues & Solutions

**Pipeline hangs during growth:**
- Large computation (10,000 cells × 1024)
- Solution: Be patient or test with fewer individuals

**KeyError: 'cpg_sites' during growth:**
- Invalid cell structure
- Solution: Use real cells from snapshots

**Mixing fails "Need X cells but only Y available":**
- Not enough year 60 cells
- Solution: Adjust mix-ratio or check snapshot

**Files exist but pipeline recreates:**
- Checkpoint mismatch
- Solution: Check pipeline_checkpoint.json

#### Testing the Pipeline
```bash
# Quick test with minimal data
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz \
    --n-individuals 3 \
    --growth-years 2

# Creates 3 individuals, grows to 4 cells each
```

#### Key Implementation Notes

1. **Sampling without replacement:**
   ```python
   sampled_indices = random.sample(range(len(year60_cells)), n_to_add)
   ```

2. **Growth model:**
   ```python
   # Each year: divide then age
   new_cells = [cell.copy() for cell in cells]  # Division
   cells.extend(new_cells)
   for cell in cells:
       cell.age()  # Methylation
   ```

3. **File list bug fix:**
   ```python
   # Always get file lists (not conditional)
   mutant_files = sorted(glob.glob(os.path.join(mutant_dir, "*.json.gz")))
   control1_files = sorted(glob.glob(os.path.join(control1_dir, "*.json.gz")))
   ```

4. **Plot improvements:**
   - Step histogram (not bar)
   - CV and MAD statistics
   - "JSD Score" label
   - Rate as percentage
   - Simple scatter (no violin)

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

### Step23: Unified Pipeline (RECOMMENDED)

```bash
cd step23
```

**Run complete pipeline:**
```bash
python run_pipeline.py --rate 0.005 --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz
```

**With custom parameters:**
```bash
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz \
    --bins 300 \
    --mix-ratio 80 \
    --n-individuals 30 \
    --growth-years 10 \
    --seed 42 \
    --output-dir data
```

**Pipeline stages:**
1. **Extract year 50 snapshot** → `data/rate_X/snapshots/year50_snapshot.json.gz`
   - Skips extraction if file already exists (caching for faster reruns)
2. **Plot JSD distribution** → `data/rate_X/plots/year50_jsd_distribution_Nbins.png`
   - Step histogram with filled area (not bars)
   - Statistics box showing: Mean, Median, SD, CV, MAD, 5%, 95%
   - Subtitle displays cell count and methylation rate as percentage
   - Filename includes number of bins for reference
3. **Create individuals**:
   - Sample 30 mutant cells (3 per JSD decile)
   - Sample 30 control1 cells (uniform distribution)
   - Save as single-cell individuals
4. **Grow individuals** (10 years, in-place):
   - Each year: duplicate all cells, then methylate
   - Year 1: 1→2 cells, Year 2: 2→4 cells, ..., Year 10: 512→1024 cells
   - Updates individual files directly (no separate lineage files)
5. **Extract year 60 snapshot** → `data/rate_X/snapshots/year60_snapshot.json.gz`
   - Also cached if already exists
6. **Mix populations** (in-place):
   - Randomly samples year 60 cells WITHOUT replacement (4,096 from 10,000)
   - Appends to grown cells to reach mix ratio (default 80% year 60, 20% grown)
   - Total: 5,120 cells per individual at 80-20 ratio
   - Each individual gets different random sample from year 60
7. **Create control2 individuals**: 30 pure year 60 samples (5,120 cells each)
8. **Analysis**:
   - Calculate mean JSD per individual
   - Compare distributions: mutant vs control1 vs control2
   - Generate plots and statistics (PNG only, no HTML files)

### Step 2: Cell Division Experiments (LEGACY)

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
- **Control lineages**: Samples 30 cells uniformly from entire population (not decile-based)
- Cell division process:
  - Each selected cell undergoes 10 rounds of division
  - Division creates two daughter cells with identical methylation states
  - Cells age between each division round
  - Results in 2^10 = 1,024 cells per lineage
- Outputs to `data/lineages/mutant/` and `data/lineages/control/`

**Plot JSD distribution:**
```bash
python plot_jsd_distribution.py ../data/snapshots/year50_snapshot.json.gz 200
```

**Batch processing for multiple rates:**
```bash
python batch_processor.py
```
Processes multiple simulation files according to `config/batch_config.yaml`.
- Automates extraction, lineage creation, and plotting
- Maintains state across runs (can resume if interrupted)
- Logs progress to `logs/` directory

### Step 3: Mixed Population Analysis (LEGACY)

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

**Batch processing for multiple rates:**
```bash
python batch_processor_step3.py
```
Processes multiple rates according to `config/step3_config.yaml`.
- Automates year 60 extraction, individual creation, and analysis
- Aggregates results across all rates for comparison
- Maintains state for resumable processing

### Testing

**Run reproducibility tests:**
```bash
cd step1/tests
python test_reproducibility.py              # Generate expected results
python test_reproducibility.py --check      # Check against saved results
```
Verifies that optimizations haven't changed simulation outputs using fixed random seeds.

**Test step pipelines:**
```bash
# Test step 3 with small test data
cd step3/scripts
python test_pipeline.py
```
Runs complete pipeline with small test data.
- Creates test individuals from test lineages
- Plots distributions with test data
- Verifies pipeline functionality

**Test lineage creation:**
```bash
cd step2
python test_sample_divide_age_lineages.py
```

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
   - **Division process**: Creates exact copy, both copies then age independently

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

### Step23 Pipeline Architecture

The unified pipeline orchestrates the complete analysis workflow:

1. **Individual Growth Model**:
   - Starts with single cell sampled from year 50
   - Each year: all cells divide (duplicate), then age (methylate)
   - Growth follows 2^n pattern: 1→2→4→8...→1024 cells over 10 years
   - All operations happen in-place within individual files

2. **Sampling Strategies**:
   - **Mutant**: Decile-based sampling (3 cells per JSD decile)
   - **Control1**: Uniform random sampling from population
   - **Control2**: Pure year 60 samples (no growth phase)

3. **Mixing Process**:
   - Configurable ratio (default 80% year 60, 20% grown)
   - Ensures all individuals have same total cell count
   - Preserves individual identity throughout pipeline

4. **File Management**:
   - Rate-specific directories prevent conflicts
   - Individual files updated in-place during growth
   - Compressed JSON format for efficiency

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

**run_simulation.py**: Alternative simulation runner
- Similar to run_large_simulation.py

**test_simulation.py**: Testing script for simulation functionality

## Important Constants

Defined in `step1/cell.py`:
- `N = 1000`: Number of methylation sites per cell
- `M = 10_000`: Number of cells in simulation
- `T_MAX = 100`: Maximum age in years
- `GENE_SIZE = 5`: Size of methylation genes
- `BASELINE_METHYLATION_DISTRIBUTION`: Reference distribution for JSD calculations

## Data Format

All simulation outputs use gzip compression (`.json.gz` extension) to reduce file size by ~90%. The scripts can read these compressed files directly without manual decompression.

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

To inspect compressed files:
```bash
zcat file.json.gz | python -m json.tool | head
```

## Code Style & Best Practices

1. **Type Hints**: All functions and methods have type annotations
2. **Docstrings**: Classes and key methods have descriptive docstrings
3. **Naming**: Use descriptive names (e.g., `cpg_sites` instead of `strain`)
4. **File Organization**: Each step has `scripts/` and `data/` subdirectories
5. **Output Organization**: All outputs go to appropriate `data/` subdirectories
6. **Random Seeds**: For reproducible simulations, set seed before creating History:
   ```python
   import random
   random.seed(42)  # Or any fixed value
   history = History(m=10000, n=1000, t_max=100, rate=0.005)
   ```
7. **Linting/Type Checking**: No specific tools configured. User should provide commands if needed.

## Common Tasks

### Running Complete Analysis for a New Rate

```bash
# Step 1: Generate base simulation
cd step1
python run_large_simulation.py --rate 0.007

# Step23: Run unified pipeline
cd ../step23
python run_pipeline.py --rate 0.007 --simulation ../step1/data/simulation_rate_0.007000_m10000_n1000_t100.json.gz
```

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

### Adjusting Step23 Pipeline Parameters

Create custom configuration:
```yaml
# step23/config/custom_config.yaml
rate: 0.005
bins: 300
mix_ratio: 70  # 70% year 60, 30% grown
n_individuals: 50  # More individuals
growth_years: 10
seed: 12345
```

Run with custom config:
```bash
python run_pipeline.py --config config/custom_config.yaml --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz
```

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
- Test pipeline available in `step3/scripts/test_pipeline.py`
- Compressed files can be inspected with: `zcat file.json.gz | python -m json.tool | head`
- Batch processors maintain state in `.batch_state.json` files for resumable processing

## Batch Processing

Both Step 2 and Step 3 include batch processors for processing multiple methylation rates:

### Step 2 Batch Processing
```bash
cd step2/scripts
python batch_processor.py
```
- Processes simulation files based on `config/batch_config.yaml`
- Automatically extracts snapshots, creates lineages, and generates plots
- Maintains state in `.batch_state.json` for resumable processing
- Logs detailed progress to `logs/` directory

### Step 3 Batch Processing
```bash
cd step3/scripts
python batch_processor_step3.py
```
- Processes multiple rates based on `config/step3_config.yaml`
- Creates mixed populations and runs analysis for each rate
- Aggregates results across all rates for comparison
- State management for fault-tolerant processing

## Step23 Implementation Details

### Key Features & Improvements

1. **Efficient Caching**: Snapshots are cached - if year 50/60 files exist, they're loaded instead of re-extracted
2. **Beautiful Plots**: 
   - Step histograms with filled areas (matching original step2 style)
   - Statistics box with: Mean, Median, SD, CV (coefficient of variation), MAD (median absolute deviation), 5%, 95%
   - Methylation rate shown as percentage (e.g., "0.5% methylation rate")
   - Clean subtitles (bins only in filename, not plot)
   - High resolution (2400x1200 at scale=2)
   - PNG only (no HTML files)
3. **Smart Sampling**: 
   - Year 60 cells sampled WITHOUT replacement for each individual
   - Each individual gets different random sample
   - Uses `random.sample()` for efficiency
4. **In-Place Growth**: Individual files updated directly during growth (no separate lineage files)

### Pipeline Flow
```python
# Simplified pipeline flow with improvements
def run_pipeline(rate, simulation_file, config):
    # 1. Extract year 50 snapshot (with caching)
    if exists(year50_snapshot):
        year50_cells = load_snapshot(year50_snapshot)
    else:
        year50_cells = extract_snapshot(simulation_file, year=50)
        save_snapshot(year50_cells, year50_snapshot)
    
    # 2. Plot JSD distribution with enhanced stats
    plot_jsd_distribution(year50_cells, bins=config['bins'], 
                         filename=f"year50_jsd_distribution_{bins}bins.png",
                         rate=rate)  # Shows as percentage in plot
    
    # 3. Create individuals (single cells)
    mutant_cells = sample_by_deciles(year50_cells, n=30)  # 3 per decile
    control1_cells = sample_uniform(year50_cells, n=30)
    
    # 4. Grow individuals (in-place modification)
    for year in range(1, 11):  # 10 years of growth
        for individual_file in all_individual_files:
            cells = load_individual(individual_file)
            # Divide: each cell creates a copy
            new_cells = [copy.deepcopy(cell) for cell in cells] * 2
            # Age: all cells methylate
            for cell in new_cells:
                cell.age_1_year()
            # Save back to SAME file
            save_individual(new_cells, individual_file)
    
    # 5. Extract year 60 snapshot (with caching)
    if exists(year60_snapshot):
        year60_cells = load_snapshot(year60_snapshot)
    else:
        year60_cells = extract_snapshot(simulation_file, year=60)
        save_snapshot(year60_cells, year60_snapshot)
    
    # 6. Mix populations (without replacement sampling)
    for individual_file in mutant_and_control1_files:
        grown_cells = load_individual(individual_file)  # 1,024 cells
        n_to_add = calculate_cells_needed(grown_cells, mix_ratio)
        
        # Sample WITHOUT replacement from year 60
        added_cells = random.sample(year60_cells, n_to_add)  # 4,096 cells
        all_cells = grown_cells + added_cells  # Simple list concatenation
        random.shuffle(all_cells)  # Mix them
        
        save_individual(all_cells, individual_file)  # Update in-place
    
    # 7. Create control2 (pure year 60, also without replacement)
    for i in range(30):
        control2_cells = random.sample(year60_cells, 5120)
        save_individual(control2_cells, control2_file)
    
    # 8. Analysis with enhanced statistics
    results = analyze_populations()  # Calculates mean JSD per individual
    plot_comparison(results)  # Creates violin and box plots
```

### Key Implementation Considerations

1. **Memory Management**: Load/save individuals one at a time during growth
2. **Random Seeds**: Set before each sampling operation for reproducibility
3. **File Format**: Maintain compatibility with existing Cell.to_dict() format
4. **Progress Tracking**: Display progress for each pipeline stage
5. **Error Handling**: Graceful failure with informative messages
6. **Parallelization**: Consider parallel processing for individual growth

## Recent Changes (August 7, 2024 Session)

### 🚀 Major Addition: Step23 Unified Pipeline
- **Complete reimplementation** of steps 2 and 3 as single efficient pipeline
- **Files created**:
  - `run_pipeline.py`: 8-stage orchestrator with intelligent skip logic
  - `pipeline_utils.py`: Core utilities for sampling, growth, and mixing
  - `pipeline_analysis.py`: Visualization with enhanced statistics
  - `pipeline_checkpoint.py`: JSON-based progress tracking system

### 🎯 Skip Logic & Checkpoint System
- **Smart skip detection**: Checks file existence AND cell counts
- **Checkpoint tracking**: Records completed stages in JSON
- **Critical bug fix**: File lists now always populated (was causing skip failures)
- **Force recreation**: Can delete directories or checkpoint to restart

### 📊 Visualization Improvements
- **Step histograms**: Replaced bar plots with filled step plots
- **Statistics box**: Added CV and MAD, removed IQR and Skew
- **X-axis label**: Now "JSD Score" instead of verbose name
- **Rate display**: Shows as percentage (e.g., "0.5% methylation rate")
- **Clean subtitles**: Bins only in filename, not plot
- **Simple scatter**: Replaced complex violin plots

### 🔧 Implementation Details
- **Sampling**: WITHOUT replacement using `random.sample()`
- **In-place updates**: No intermediate lineage files
- **Snapshot caching**: Automatic reuse when files exist
- **Growth model**: Divide (copy) then age (methylate) each year
- **File structure**: All cells have required fields (cpg_sites, rate, etc.)

### 📝 Documentation
- **README.md**: Comprehensive step23 section with examples
- **CLAUDE.md**: This exhaustive implementation guide
- **Inline docs**: All functions properly documented

### 📁 Directory Reorganization
- ✅ Moved step2 and step3 to `legacy/` directory  
- ✅ Positioned step23 as the recommended approach

## Previous Changes (Earlier Sessions)
- **Major restructuring**: Organized into 3 sequential steps
- **Step 1**: Dedicated directory with `data/` subdirectory
- **Step 2**: Unified lineage creation
- **Step 3**: Uses step2's control lineages
- **CLI improvements**: Command-line arguments added
- **Batch processing**: Automated multi-rate processing
- **Precision**: Rate precision to 6 decimal places