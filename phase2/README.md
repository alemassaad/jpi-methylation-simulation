# Phase 2: Analysis Pipeline

A refactored analysis pipeline using object-oriented design with `PetriDish` and `Cell` classes from phase1.

## 🚨 Breaking Changes (Latest)

- **New Lean JSON Format Required**: Phase 1 simulations must use the new lean format (~90% smaller)
- **No Backward Compatibility**: Old JSON format is not supported
- **Config File Support**: YAML configuration files now preferred over complex CLI commands
- **Updated Defaults**: Optimized for testing (100 sites, 50 years, growth phase 6)

## Quick Start

```bash
# Using config file (recommended)
python run_pipeline.py --config configs/standard.yaml

# Standard analysis with CLI
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --individual-growth-phase 7 --seed 42

# Gene-specific rates
python run_pipeline.py \
    --gene-rate-groups "50:0.004,50:0.005,50:0.006,50:0.007" \
    --simulation ../phase1/data/gene_rates_*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 --seed 42

# With advanced features
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --individual-growth-phase 7 \
    --uniform-mixing \
    --normalize-size \
    --seed 42
```

## Configuration

### Config File Structure

Create a YAML file (e.g., `my_config.yaml`):

```yaml
simulation:
  rate: 0.005              # Or use gene_rate_groups
  first_snapshot: 50
  second_snapshot: 60
  individual_growth_phase: 7

analysis:
  n_quantiles: 10
  cells_per_quantile: 3
  mix_ratio: 80

options:
  uniform_mixing: false
  normalize_size: false

seed: 42
```

### Command-Line Parameters

#### Required
- `--config`: Path to YAML config file (optional, but recommended)
- `--rate` OR `--gene-rate-groups`: Methylation rate configuration (must match simulation)
- `--simulation`: Path to phase1 simulation file (supports wildcards)

#### Snapshot Configuration
- `--first-snapshot`: Year for initial cell extraction (default: 50)
- `--second-snapshot`: Year for mixing cells (default: 60)
- `--individual-growth-phase`: Years of exponential growth before homeostasis (default: 7)
  - Determines target population: 2^growth_phase cells
  - Example: 7 → 128 cells, 8 → 256 cells, 10 → 1024 cells

#### Sampling Configuration
- `--n-quantiles`: Number of quantiles for sampling (default: 10)
- `--cells-per-quantile`: Cells per quantile (default: 3)
- `--mix-ratio`: Percentage of snapshot cells in final mix (default: 80)

#### Advanced Features
- `--uniform-mixing`: All individuals receive same snapshot cells (3-step normalization)
- `--normalize-size`: Normalize to same size before mixing (median - 0.5σ)
- `--plot-individuals`: Generate growth trajectory plots for each individual

#### Other Options
- `--seed`: Random seed for reproducibility (default: 42)
- `--bins`: Histogram bins for JSD plot (default: 200)
- `--force-reload`: Force reload of cached snapshots
- `--force-recreate`: Force recreation of individuals

## Pipeline Stages

1. **Extract First Snapshot** (e.g., year 50)
   - Load cells from simulation at specified year
   - Cache in `snapshots/year{N}_snapshot.json.gz`

2. **Plot JSD Distribution**
   - Histogram with statistical overlays
   - Shows population methylation state

3. **Create Initial Individuals**
   - **Mutant**: Quantile-based stratified sampling
   - **Control1**: Uniform random sampling
   - Each starts with single cell

4. **Grow Individuals**
   - Exponential growth for `individual-growth-phase` years
   - Then homeostasis: divide → cull → methylate
   - Maintains ~2^growth_phase cells

5. **Extract Second Snapshot**
   - Load cells from year = first_snapshot + (second - first)

6. **Mix Populations** (with optional normalization)
   - If `--normalize-size`: Apply median - 0.5σ normalization
   - If `--uniform-mixing`: **NEW 3-step process** (fixed):
     1. Normalize all individuals to same size
     2. Create uniform pool sized for normalized individuals
     3. Add identical pool to each → perfect size consistency
   - Mix grown individuals with snapshot cells

7. **Create Control2**
   - Pure second snapshot cells
   - Count adjusted if normalization applied

8. **Analysis & Visualization**
   - Statistical comparisons (t-tests)
   - JSD distribution plots
   - Cell-level and individual-level analyses

## Output Structure

```
data/
└── rate_0.00500-grow13-sites1000-years100/  # Or gene_rates_50x0.004_50x0.005.../
    └── snap50to60-growth7-quant10x3-mix80[u][n]-seed42-YYYYMMDDHHMMSS/
        ├── snapshots/
        │   ├── year50_snapshot.json.gz
        │   └── year60_snapshot.json.gz
        ├── individuals/
        │   ├── mutant/     # Quantile-sampled individuals
        │   ├── control1/   # Uniformly-sampled individuals
        │   └── control2/   # Pure snapshot individuals
        ├── results/
        │   ├── year50_jsd_distribution_200bins.png
        │   ├── cell_jsd_comparison.png
        │   ├── cell_jsd_analysis.json         # Cell-level JSD analysis
        │   ├── gene_jsd_analysis.json          # Gene-level JSD distributions
        │   ├── mixing_statistics.json          # Only when --uniform-mixing
        │   └── pipeline_metadata.json          # Pipeline parameters only
        └── individual_plots/  # If --plot-individuals used
            ├── mutant_01_jsd.png
            ├── mutant_01_methylation.png
            └── mutant_01_combined.png
```

### JSON Output Files

#### cell_jsd_analysis.json (NEW - Consolidated Format)
Contains all cell JSD statistics and individual-level data:
```json
{
  "summary_statistics": {
    "mutant": {
      "mean": 0.234,
      "std": 0.045,
      "median": 0.229,
      "min": 0.156,
      "max": 0.341,
      "n_individuals": 7
    },
    "control1": { ... },
    "control2": { ... }
  },
  "statistical_tests": {
    "mutant_vs_control1": {
      "t_statistic": 2.34,
      "p_value": 0.023
    },
    "mutant_vs_control2": { ... },
    "control1_vs_control2": { ... }
  },
  "individual_means": {
    "mutant": [0.234, 0.189, ...],    // Mean cell JSD per individual
    "control1": [0.198, 0.211, ...],
    "control2": [0.165, 0.171, ...]
  }
}
```

#### gene_jsd_analysis.json (NEW)
Gene-level JSD distributions for each gene across all batches:
```json
{
  "gene_jsd_distributions": {
    "gene_0": {
      "mutant": [0.234, 0.189, ...],     // JSD values from mutant individuals
      "control1": [0.198, 0.211, ...],   // JSD values from control1 individuals
      "control2": [0.165, 0.171, ...]    // JSD values from control2 individuals
    },
    "gene_1": { ... },
    // ... for all genes
  },
  "gene_metadata": {
    "n_genes": 200,
    "n_individuals_per_batch": {
      "mutant": 8,
      "control1": 8,
      "control2": 8
    },
    "gene_rate_groups": [...]  // Only if using gene-specific rates
  }
}
```

#### pipeline_metadata.json
Complete pipeline configuration and runtime information (parameters only, no results).

#### mixing_statistics.json
Created only when `--uniform-mixing` is used. Contains pool size and normalization details.

### Directory Naming Convention
- Base: `snap{first}to{second}-growth{phase}-quant{n}x{m}-mix{ratio}`
- Suffixes:
  - No suffix: Standard (independent sampling)
  - `u`: Uniform mixing
  - `n`: Size normalization  
  - `un`: Both features
- Timestamp: `YYYYMMDDHHMMSS` format for chronological sorting and uniqueness

## Features

### New Lean JSON Format
- **~90% file size reduction** compared to old format
- Parameters stored once at top level (not per cell)
- Cells contain only essential data: methylated array and cell_JSD
- **No backward compatibility** - old format files not supported
- Automatic compression with .json.gz

### Object-Oriented Design
- Uses `Cell` and `PetriDish` objects throughout
- No dictionary conversions during processing
- Clean separation of concerns

### Homeostasis Modeling
- Biologically realistic growth dynamics
- Exponential phase followed by steady-state
- Stochastic population maintenance

### Uniform Mixing (FIXED)
- **NEW**: 3-step normalization process eliminates size mismatches
- **Step 1**: Normalize all individuals to minimum size
- **Step 2**: Create uniform pool sized exactly for normalized individuals  
- **Step 3**: Add identical pool → perfect size consistency
- **Benefits**: No warnings, fair comparisons, preserved diversity

### Size Normalization
- Addresses homeostasis variation
- Median - 0.5σ threshold
- Typically retains ~67% of individuals
- Ensures fair comparison

### Full Reproducibility
- Comprehensive random seeding
- Deterministic file operations
- Timestamp-based unique directories (YYYYMMDDHHMMSS)

## Analysis Tools

```bash
# Compare two runs for reproducibility
python tools/compare_two_runs.py --dir1 path1 --dir2 path2

# Analyze individual size distributions
python analyze_individual_sizes.py [output_directory]
```

## Testing

```bash
cd tests
# Run comprehensive tests
python test_normalization_comprehensive.py
python test_gene_rate_support.py
python test_final_integration.py

# Run all uniform mixing tests
python run_all_uniform_tests.py
```

## Troubleshooting

### "Old format detected" Error
- Phase 1 simulation must use the new lean JSON format
- Re-run phase 1 simulation with latest code
- Old format files are not supported (no backward compatibility)

### High CV Warning
If you see "High variation in individual sizes (CV > 20%)", consider:
- Using `--normalize-size` to equalize sizes
- Adjusting `--individual-growth-phase`
- Increasing initial sample size

### All Individuals Excluded
If normalization excludes all individuals:
- Check growth parameters
- Consider larger `--individual-growth-phase`
- Increase `--cells-per-quantile`

### Gene Rate Groups
- Use `--gene-rate-groups` instead of `--rate` for gene-specific rates
- Format: `"50:0.004,50:0.005,50:0.006,50:0.007"`
- Must match the original phase 1 simulation configuration

### Memory Issues
For large simulations:
- Process in smaller batches
- Use `--force-reload` sparingly
- Consider reducing `--n-quantiles`

### Config File Not Loading
- Install PyYAML: `pip install pyyaml`
- Check YAML syntax (proper indentation)
- CLI arguments override config values

## Project Structure

```
phase2/
├── README.md               # Main documentation
├── run_pipeline.py         # Main pipeline entry point
│
├── core/                   # Core pipeline modules
│   ├── __init__.py         # Package exports
│   ├── pipeline_utils.py  # Helper functions for pipeline
│   ├── pipeline_analysis.py # Analysis and visualization functions
│   └── path_utils.py       # Path generation utilities
│
├── visualization/          # Visualization modules
│   ├── __init__.py         # Package exports
│   ├── plot_individuals.py # Individual growth trajectory plots
│   └── analyze_individual_sizes.py # Size distribution analysis
│
├── configs/                # Configuration files
│   ├── config_default.yaml # Default configuration
│   ├── quick_test.yaml    # Quick testing configuration
│   ├── full_analysis.yaml # Full analysis configuration
│   └── uniform_mixing.yaml # Uniform mixing configuration
│
├── scripts/                # Standalone scripts
│   ├── complete_analysis.py # Run complete analysis
│   ├── create_control2.py  # Create control2 individuals
│   └── clean_individuals.sh # Clean individual files
│
├── tools/                  # Analysis tools
│   └── compare_two_runs.py # Compare two pipeline runs
│
├── tests/                  # Test suite (organized)
│   ├── config/            # Configuration tests
│   ├── integration/       # Integration tests
│   └── unit/              # Unit tests
│
├── docs/                   # Additional documentation
│   ├── NORMALIZATION_USAGE.md
│   ├── JSON_CONSOLIDATION_PLAN.md
│   └── PHASE2_REFACTORING_PLAN.md
│
├── data/                   # Output directory (gitignored)
└── tmp/                    # Temporary files (gitignored)