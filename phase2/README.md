# Phase 2: Data Generation Pipeline

A modular pipeline for generating structured datasets from phase1 simulations. Phase2 creates cell populations with uniform mixing (all individuals receive the same snapshot cells) for reproducible analysis. All plotting and visualization has been moved to phase3.

## Current Structure (After Cleanup)

Phase 2 now contains only essential data generation components (10 files total):
- **Core Scripts**: `run_pipeline.py`, `extract_snapshots.py`, `simulate_individuals.py`, `create_control2.py`
- **Core Modules**: `pipeline_utils.py`, `individual_helpers.py`, `path_utils.py`, `validation.py`
- **Configuration**: `config_default.yaml` (loaded by default)

All testing and plotting functionality has been removed for a clean, focused pipeline.

## Quick Start

**Important Update (Sept 2025)**: Phase 2 now outputs directly to the Phase 1 simulation directory. No need to specify an output directory - it's automatic!

```bash
# Run complete pipeline with defaults
cd phase2
python run_pipeline.py --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json

# Standard analysis
python run_pipeline.py \
    --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json \
    --first-snapshot 30 \
    --second-snapshot 50 \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --individual-growth-phase 7 \
    --mix-ratio 80 \
    --seed 42

# Quick test (smaller parameters)
python run_pipeline.py \
    --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json \
    --first-snapshot 30 \
    --second-snapshot 50 \
    --n-quantiles 4 \
    --cells-per-quantile 2 \
    --individual-growth-phase 6 \
    --mix-ratio 70

# Config is loaded automatically, CLI args override
python run_pipeline.py --simulation ../data/gene_rates_*/size*-seed*-*/simulation.json
```

## Pipeline Stages

### Stage 1-2: Extract Snapshots (`extract_snapshots.py`)
- Extracts year data from two time points as direct copies
- Saves snapshots with year key wrapper (e.g., {"30": {...}})
- Preserves all fields including gene_jsd if present
- Creates metadata.json from simulation config

### Stage 3-5: Simulate Individuals (`simulate_individuals.py`)
- **Stage 3**: Create initial individuals
  - Mutant: Quantile-based sampling
  - Control1: Uniform random sampling
- **Stage 4**: Grow populations
  - Exponential growth phase
  - Homeostasis with random culling
- **Stage 5**: Mix with snapshot cells
  - Uniform mixing (all individuals get same cells)
  - Optional size normalization

### Stage 6: Create Control2 (`create_control2.py`)
- Pure snapshot populations (no growth)
- Matches mutant/control1 population sizes
- Uses uniform pool (always applied)

### Next Step: Analysis (Phase 3)
- Use phase3 to analyze the generated data
- All plotting and visualization moved to phase3
- Multiple analyses possible on same phase2 data

## Output Directory Structure

**New in Sept 2025**: Phase 2 outputs are now created as subdirectories within the Phase 1 simulation directory:

```
data/gene_rates_*/size*-seed*-TIMESTAMP/            # Phase 1 directory
├── simulation.json                                 # Phase 1 simulation
└── snap30to50-growth7-quant10x3-mix80u-seed42-TIMESTAMP/  # Phase 2 output
    ├── snapshots/
    │   ├── year30_snapshot.json
    │   ├── year50_snapshot.json
    │   └── metadata.json
    └── individuals/
        ├── mutant/
        ├── control1/
        ├── control2/
        ├── uniform_pool.json
        └── mixing_metadata.json
```

Multiple Phase 2 runs with the same or different parameters will each get their own timestamped subdirectory, preventing collisions.

## Command-Line Options

**Note**: Phase 2 no longer includes plotting options or output directory specification.

### Required Arguments
- `--simulation`: Path to phase1 simulation file (supports wildcards)

### Snapshot Parameters
- `--first-snapshot`: First year to extract (default: 20)
- `--second-snapshot`: Second year to extract (default: 30)

### Sampling Parameters
- `--n-quantiles`: Number of quantiles for mutant sampling (default: 3)
- `--cells-per-quantile`: Cells to sample per quantile (default: 2)

### Growth Parameters
- `--individual-growth-phase`: Years of exponential growth (default: 5)
  - 6 = 64 cells, 7 = 128 cells, 8 = 256 cells
- `--mix-ratio`: Percentage from second snapshot (default: 60)

### Advanced Options
- `--compress`: Compress output files (.json.gz)
- `--no-compress`: Don't compress output files (.json) - default
- `--seed`: Random seed for reproducibility (default: 42)

## Configuration Files

The default configuration (`config_default.yaml`) is always loaded automatically from the phase2 directory. You can override defaults with command-line arguments or a custom config file:

```yaml
# Custom config example
first_snapshot: 25
second_snapshot: 45
n_quantiles: 5
cells_per_quantile: 4
individual_growth_phase: 8
mix_ratio: 85
seed: 123
```

Use with: `python run_pipeline.py --config my_custom.yaml --simulation ../phase1/data/*/simulation.json.gz`

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

# For analysis, use phase3 instead:
cd ../phase3
python run_analysis.py \
    --phase2-dir ../phase2/data/my_analysis \
    --simulation ../phase1/data/.../simulation.json.gz
```

## Output Structure

```
data/{rate_info}/snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}-seed{S}-{timestamp}/
├── snapshots/
│   ├── year30_snapshot.json[.gz]  # Format: {"30": {"cells": [...], "gene_jsd": [...]}}
│   ├── year50_snapshot.json[.gz]  # Format: {"50": {"cells": [...], "gene_jsd": [...]}}
│   └── metadata.json               # Config from simulation
├── individuals/
│   ├── mutant/
│   │   └── individual_*.json.gz
│   ├── control1/
│   │   └── individual_*.json.gz
│   ├── control2/
│   │   └── individual_*.json.gz
│   └── mixing_metadata.json
# No results/ directory in phase2 - data only
# For analysis results, use phase3
```

## Mixing Mode

### Uniform Mixing (Always Enabled)
- All individuals receive the exact same snapshot cells
- Eliminates sampling variation between individuals
- Focuses analysis on biological differences rather than sampling noise
- Directory suffix always includes 'u' (e.g., mix80u)
- Provides reproducible, comparable results across batches

## Size Normalization

### Population Normalization (Always Applied)
- All populations normalized to median - 0.5σ threshold
- Removes outlier individuals below threshold
- Ensures fair comparison between individuals
- Better for analyzing methylation patterns

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
- Snapshots are cached within each run directory
- Generate data once, then run multiple analyses in phase3

## Migration from Old Pipeline

The new three-phase architecture is NOT backward compatible with old pipeline runs. To migrate:

1. Re-run phase1 simulations if using old JSON format
2. Use new phase2 `run_pipeline.py` for data generation only
3. Use new phase3 `run_analysis.py` for all analysis and plotting
4. Update automation scripts to use phase2 + phase3 workflow
5. Config files remain compatible within each phase

## Integration with Phase 3

Phase 2 generates data that phase3 analyzes:

```bash
# Generate data with phase2
cd phase2
python run_pipeline.py --simulation ../phase1/data/*/simulation.json.gz

# Analyze data with phase3
cd ../phase3
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/*/simulation.json.gz
```

### Benefits of Separation
- **Efficiency**: Generate data once, analyze multiple ways
- **Flexibility**: Different analysis parameters on same data
- **Scalability**: Batch analysis of multiple phase2 runs
- **Development**: Independent analysis development

## Future Plans

- Potential parallelization of individual population simulations
- Distributed processing for large-scale studies
- Enhanced validation and quality control
- Integration with external analysis tools