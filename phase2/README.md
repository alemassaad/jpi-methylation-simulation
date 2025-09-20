# Phase 2: Simulation and Growth Pipeline

A modular pipeline for creating and simulating cell populations, split into independent scripts for flexibility and maintainability.

## 🆕 Major Architecture Change (2025-01-20)

Phase 2 has been completely reorganized from a monolithic script into 4 modular components:
- **`extract_snapshots.py`**: Extract cell snapshots from phase1 simulations
- **`simulate_individuals.py`**: Create, grow, and mix cell populations
- **`create_control2.py`**: Create control populations from snapshots
- **`analyze_and_plot.py`**: Generate all analysis and visualizations
- **`run_pipeline.py`**: Main driver that orchestrates all scripts

## Quick Start

```bash
# Run complete pipeline with defaults
cd phase2
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz

# Standard analysis
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 \
    --second-snapshot 50 \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --individual-growth-phase 7 \
    --mix-ratio 80 \
    --seed 42

# Quick test (smaller parameters)
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 \
    --second-snapshot 50 \
    --n-quantiles 4 \
    --cells-per-quantile 2 \
    --individual-growth-phase 6 \
    --mix-ratio 70

# With config file
python run_pipeline.py \
    --config configs/quick_test.yaml \
    --simulation ../phase1/data/*/simulation.json.gz
```

## Pipeline Stages

### Stage 1-2: Extract Snapshots (`extract_snapshots.py`)
- Extracts cells from two time points in the simulation
- Saves snapshots for reuse across multiple analyses
- Validates gene rate groups consistency

### Stage 3-5: Simulate Individuals (`simulate_individuals.py`)
- **Stage 3**: Create initial individuals
  - Mutant: Quantile-based sampling
  - Control1: Uniform random sampling
- **Stage 4**: Grow populations
  - Exponential growth phase
  - Homeostasis with random culling
- **Stage 5**: Mix with snapshot cells
  - Independent or uniform mixing modes
  - Optional size normalization

### Stage 6: Create Control2 (`create_control2.py`)
- Pure snapshot populations (no growth)
- Matches mutant/control1 population sizes
- Uses uniform pool if uniform mixing was applied

### Stage 7: Analysis and Plotting (`analyze_and_plot.py`)
- Cell-level and gene-level metrics
- Population comparisons
- Individual growth trajectories
- Timeline visualizations

## Command-Line Options

### Required Arguments
- `--simulation`: Path to phase1 simulation file (supports wildcards)

### Snapshot Parameters
- `--first-snapshot`: First year to extract (default: 30)
- `--second-snapshot`: Second year to extract (default: 50)

### Sampling Parameters
- `--n-quantiles`: Number of quantiles for mutant sampling (default: 10)
- `--cells-per-quantile`: Cells to sample per quantile (default: 3)

### Growth Parameters
- `--individual-growth-phase`: Years of exponential growth (default: 7)
  - 6 = 64 cells, 7 = 128 cells, 8 = 256 cells
- `--mix-ratio`: Percentage from second snapshot (default: 80)

### Advanced Options
- `--uniform-mixing`: Use same snapshot cells for all individuals
- `--normalize-size`: Normalize populations to same size before mixing
- `--bins`: Number of bins for histograms (default: 200)
- `--max-gene-plots`: Limit number of per-gene plots
- `--no-compress`: Save uncompressed JSON files
- `--force-reload`: Force re-extraction of snapshots
- `--force-recreate`: Force recreation of individuals
- `--seed`: Random seed for reproducibility (default: 42)

## Configuration Files

Create YAML configuration files in `configs/`:

```yaml
# configs/my_config.yaml
simulation: ../phase1/data/*/simulation.json.gz
first_snapshot: 30
second_snapshot: 50
n_quantiles: 10
cells_per_quantile: 3
individual_growth_phase: 7
mix_ratio: 80
uniform_mixing: true
normalize_size: true
seed: 42
```

Use with: `python run_pipeline.py --config configs/my_config.yaml`

## Running Individual Stages

Each script can be run independently for debugging or custom workflows:

```bash
# Extract snapshots only
python extract_snapshots.py \
    --simulation ../phase1/data/.../simulation.json.gz \
    --first-snapshot 30 \
    --second-snapshot 50 \
    --output-dir data/my_analysis

# Simulate individuals only
python simulate_individuals.py \
    --base-dir data/my_analysis \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --growth-phase 7 \
    --mix-ratio 80

# Create control2 only
python create_control2.py \
    --base-dir data/my_analysis \
    --seed 342

# Analysis only (can be run multiple times)
python analyze_and_plot.py \
    --base-dir data/my_analysis \
    --simulation ../phase1/data/.../simulation.json.gz \
    --bins 200
```

## Output Structure

```
data/{rate_info}/snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}-seed{S}-{timestamp}/
├── snapshots/
│   ├── year30_snapshot.json.gz
│   ├── year50_snapshot.json.gz
│   └── metadata.json
├── individuals/
│   ├── mutant/
│   │   └── individual_*.json.gz
│   ├── control1/
│   │   └── individual_*.json.gz
│   ├── control2/
│   │   └── individual_*.json.gz
│   └── mixing_metadata.json
└── results/
    ├── cell_metrics/
    │   ├── distributions/
    │   ├── comparisons/
    │   ├── individual_trajectories/
    │   └── timeline/
    └── gene_metrics/
        ├── distributions/
        ├── comparisons/
        ├── per_gene/
        └── timeline/
```

## Mixing Modes

### Independent Mixing (Default)
- Each individual randomly samples from snapshot
- Natural variation between individuals
- Standard statistical analysis

### Uniform Mixing (`--uniform-mixing`)
- All individuals use same snapshot cells
- Reduces sampling noise
- Better for comparing growth effects

## Size Normalization

### Without Normalization (Default)
- Natural size variation from growth stochasticity
- Some individuals may go extinct
- Realistic population dynamics

### With Normalization (`--normalize-size`)
- All populations normalized to median - 0.5σ
- Removes extinct individuals
- Better for comparing methylation patterns

## Gene-Specific Rates

For simulations with gene-specific methylation rates:

```bash
python run_pipeline.py \
    --simulation ../phase1/data/gene_rates_*/simulation.json.gz \
    --first-snapshot 30 \
    --second-snapshot 50
# Rate configuration auto-detected from simulation
```

## Testing

```bash
# Run unit tests
cd phase2/tests
python test_validation.py
python test_gene_jsd_extraction.py

# Run integration test
python test_pipeline_with_config.py

# Quick pipeline test
cd phase2
python run_pipeline.py \
    --config configs/quick_test.yaml \
    --simulation ../phase1/data/*/simulation.json.gz
```

## Troubleshooting

### Common Issues

1. **"No module named 'yaml'"**: Install with `pip install pyyaml`
2. **"Simulation not found"**: Check path and wildcards
3. **"Validation failed"**: Ensure consistent gene rate groups
4. **Memory issues**: Reduce `--n-quantiles` or `--cells-per-quantile`

### Performance Tips

- Use compressed files (`.json.gz`) for large simulations
- Run stages individually for debugging
- Reuse snapshots with `--force-reload` only when needed
- Limit gene plots with `--max-gene-plots` for faster analysis

## Migration from Old Pipeline

The new modular architecture is NOT backward compatible with old pipeline runs. To migrate:

1. Re-run phase1 simulations if using old JSON format
2. Use new `run_pipeline.py` instead of old monolithic script
3. Update any automation scripts to use new command structure
4. Config files remain compatible

## Future Plans

- Phase 3 will contain only analysis/plotting (current Stage 7)
- Phase 2 will be pure simulation (Stages 1-6)
- Potential parallelization of individual simulations
- Web interface for interactive analysis