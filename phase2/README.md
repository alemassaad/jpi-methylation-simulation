# Phase 2: Data Generation Pipeline

A modular pipeline for generating structured datasets from phase1 simulations. Phase2 creates cell populations and snapshot data for downstream analysis, with no plotting or visualization (moved to phase3).

## 🆕 Major Architecture Changes

### Phase Split (2025-09-20)
Phase 2 is now **data generation only**. All analysis and plotting has been moved to **phase3** for clean separation of concerns.

### Modular Organization (2025-01-20)
Phase 2 consists of 3 modular data generation components:
- **`extract_snapshots.py`**: Extract cell snapshots from phase1 simulations
- **`simulate_individuals.py`**: Create, grow, and mix cell populations
- **`create_control2.py`**: Create control populations from snapshots
- **`run_pipeline.py`**: Main driver that orchestrates all scripts

For analysis and visualization, use **phase3** instead.

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

### Next Step: Analysis (Phase 3)
- Use phase3 to analyze the generated data
- All plotting and visualization moved to phase3
- Multiple analyses possible on same phase2 data

## Command-Line Options

**Note**: Phase 2 no longer includes plotting options. For visualization, use phase3 after running phase2.

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
- `--force-reload`: Force re-extraction of snapshots
- `--force-recreate`: Force recreation of individuals
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
# No results/ directory in phase2 - data only
# For analysis results, use phase3
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