# Phase 2: Analysis Pipeline

A refactored analysis pipeline using object-oriented design with `PetriDish` and `Cell` classes from phase1.

## Quick Start

```bash
# Standard analysis
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --individual-growth-phase 7 --seed 42

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

## Parameters

### Required
- `--rate`: Methylation rate (must match simulation)
- `--simulation`: Path to phase1 simulation file

### Snapshot Configuration
- `--first-snapshot`: Year for initial cell extraction (default: 50)
- `--second-snapshot`: Year for mixing cells (default: 60)
- `--individual-growth-phase`: Years of exponential growth before homeostasis (default: 7)
  - Determines target population: 2^growth_phase cells
  - Example: 7 → 128 cells, 8 → 256 cells, 10 → 1024 cells

### Sampling Configuration
- `--n-quantiles`: Number of quantiles for sampling (default: 10)
- `--cells-per-quantile`: Cells per quantile (default: 3)
- `--mix-ratio`: Percentage of snapshot cells in final mix (default: 80)

### Advanced Features
- `--uniform-mixing`: All individuals receive same snapshot cells
- `--normalize-size`: Normalize to same size before mixing (median - 0.5σ)

### Other Options
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
└── rate_0.00500-grow13-sites1000-years100/
    └── snap50to60-growth7-quant10x3-mix80[u][n]-seed42-XXXX/
        ├── snapshots/
        │   ├── year50_snapshot.json.gz
        │   └── year60_snapshot.json.gz
        ├── individuals/
        │   ├── mutant/     # Quantile-sampled individuals
        │   ├── control1/   # Uniformly-sampled individuals
        │   └── control2/   # Pure snapshot individuals
        └── results/
            ├── year50_jsd_distribution_200bins.png
            ├── jsd_comparison.png
            ├── statistics.json
            └── pipeline_metadata.json
```

### Directory Naming Convention
- Base: `snap{first}to{second}-growth{phase}-quant{n}x{m}-mix{ratio}`
- Suffixes:
  - No suffix: Standard (independent sampling)
  - `u`: Uniform mixing
  - `n`: Size normalization  
  - `un`: Both features

## Features

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
- MD5 hashing for unique directories

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
```

## Troubleshooting

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

### Memory Issues
For large simulations:
- Process in smaller batches
- Use `--force-reload` sparingly
- Consider reducing `--n-quantiles`