# Phase 2: Data Generation Pipeline

A streamlined single-file pipeline for generating structured datasets from phase1 simulations. Phase2 creates cell populations with uniform mixing for reproducible analysis. All plotting and visualization is done in phase3.

## Quick Start

**Important**: Phase 2 outputs directly to the Phase 1 simulation directory automatically!

```bash
# Run complete pipeline with defaults
python phase2_pipeline.py --simulation ../phase1/data/gene_rates_*/simulation.json

# Standard parameters
python phase2_pipeline.py \
    --simulation ../phase1/data/gene_rates_*/simulation.json \
    --first-snapshot 30 \
    --second-snapshot 50 \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --individual-growth-phase 7 \
    --mix-ratio 80 \
    --seed 42

# Quick test (smaller parameters)
python phase2_pipeline.py \
    --simulation ../phase1/data/gene_rates_*/simulation.json \
    --n-quantiles 4 \
    --cells-per-quantile 2 \
    --individual-growth-phase 6

# Using config file
python phase2_pipeline.py \
    --config my_config.yaml \
    --simulation ../phase1/data/gene_rates_*/simulation.json
```

## Architecture

**Single-File Pipeline**: All functionality consolidated into `phase2_pipeline.py` (~800 lines).

### Pipeline Stages (executed sequentially in memory)

1. **Extract Snapshots**: Load simulation, extract two time points
2. **Create Individuals**: Sample cells (quantile-based for Test 2, uniform for Test 1)
3. **Grow Populations**: Exponential growth + homeostasis
4. **Normalize & Mix**: Size normalization + mix with common snapshot pool
5. **Create Control**: Pure snapshot populations
6. **Save Outputs**: Write all individuals to disk

### Key Features

- **No subprocesses**: Direct function calls for better performance
- **In-memory data**: Keep cells in memory between stages
- **Faster execution**: No process spawning or temporary file I/O
- **Easier debugging**: Step through entire pipeline
- **Same interface**: Command-line arguments unchanged

## Output Directory Structure

Phase 2 outputs are created as subdirectories within the Phase 1 simulation directory:

```
phase1/data/gene_rates_*/size*-seed*-TIMESTAMP/            # Phase 1 directory
├── simulation.json.gz                                     # Phase 1 simulation
└── snap30to50-growth7-quant10x3-mix80-seed42-TIMESTAMP/  # Phase 2 output
    ├── snapshots/
    │   ├── year30_snapshot.json
    │   ├── year50_snapshot.json
    │   └── metadata.json
    └── individuals/
        ├── test2/           # Test 2: Quantile-sampled populations
        ├── test1/           # Test 1: Random-sampled populations
        ├── control/         # Control: Pure snapshot populations
        ├── common_pool.json # Shared snapshot cells
        └── mixing_metadata.json
```

Multiple Phase 2 runs get their own timestamped subdirectories, preventing collisions.

## Command-Line Options

### Required
- `--simulation`: Path to phase1 simulation file (supports wildcards)

### Snapshot Parameters
- `--first-snapshot`: First year to extract (default: 30)
- `--second-snapshot`: Second year to extract (default: 50)

### Sampling Parameters
- `--n-quantiles`: Number of quantiles for Test 2 sampling (default: 4)
- `--cells-per-quantile`: Cells per quantile (default: 2)

### Growth Parameters
- `--individual-growth-phase`: Years of exponential growth (default: 6)
  - 6 = 64 cells, 7 = 128 cells, 8 = 256 cells
- `--mix-ratio`: Percentage from second snapshot (default: 70)

### Other Options
- `--config`: YAML configuration file
- `--seed`: Random seed for reproducibility (default: 42)
- `--compress` / `--no-compress`: Compress output files (default: false)

## Configuration Files

The default configuration (`config_default.yaml`) is loaded automatically. Override with command-line arguments or custom config:

```yaml
# Example custom config
first_snapshot: 25
second_snapshot: 45
n_quantiles: 5
cells_per_quantile: 4
individual_growth_phase: 8
mix_ratio: 85
seed: 123
```

Use with: `python phase2_pipeline.py --config my_config.yaml --simulation ../phase1/data/*/simulation.json`

## Key Concepts

### Uniform Mixing (Always Enabled)
- All individuals receive the exact same snapshot cells
- Eliminates sampling variation between individuals
- Focuses analysis on biological differences
- Provides reproducible, comparable results

### Size Normalization (Always Applied)
- Populations normalized to median - 0.5σ threshold
- Removes outlier individuals below threshold
- Ensures fair comparison between individuals
- Typically retains ~67% of individuals

### Three Population Types
1. **Test 2**: Quantile-based sampling (captures full JSD spectrum)
2. **Test 1**: Uniform random sampling (baseline comparison)
3. **Control**: Pure second snapshot (no growth, mixing reference)

## Integration with Phase 3

Phase 2 generates data that phase3 analyzes:

```bash
# Generate data with phase2
python phase2_pipeline.py --simulation ../phase1/data/*/simulation.json

# Analyze data with phase3
cd ../phase3
python run_pipeline.py \
    --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/ \
    --simulation ../phase1/data/{phase1_dir}/simulation.json
```

### Benefits of Separation
- **Efficiency**: Generate data once, analyze multiple ways
- **Flexibility**: Different analysis parameters on same data
- **Scalability**: Batch analysis of multiple phase2 runs
- **Development**: Independent analysis development

## Troubleshooting

### Common Issues

1. **"No module named 'yaml'"**: Install with `pip install pyyaml`
2. **"Simulation not found"**: Check path and wildcards
3. **Memory issues**: Reduce `--n-quantiles` or `--cells-per-quantile`
4. **All individuals excluded**: Increase growth phase or adjust normalization

### Performance Tips

- Use compressed files (`.json.gz`) for large simulations
- Default config provides good balance of speed and comprehensiveness
- Generate data once, then run multiple analyses in phase3
- Pipeline runs ~10x faster than old multi-script version

## Migration from Old Pipeline

The old multi-script architecture (`run_pipeline.py`, `extract_snapshots.py`, etc.) has been replaced with a single-file pipeline.

**Old backup available**: `run_pipeline_OLD.py`

### What Changed
- ✅ Command-line interface: **Same**
- ✅ Output format: **Same**
- ✅ Config files: **Compatible**
- ✅ Phase 3 integration: **No changes needed**
- 🚀 Performance: **Faster** (no subprocess overhead)
- 🧹 Code: **Simpler** (800 lines vs 3000 lines)

### Update Your Scripts
```bash
# Old
python run_pipeline.py --simulation ...

# New
python phase2_pipeline.py --simulation ...
```

## Dependencies

Required packages (via `pip install -r requirements.txt`):
- `numpy` (>=1.19.0)
- `scipy` (>=1.7.0) - used by phase1
- `pyyaml` (>=6.0) - for config files
- `plotly`, `kaleido` - for phase3 only

## Development

The single-file architecture makes phase2 easy to understand and modify:

- All code in one place (`phase2_pipeline.py`)
- Clear class structure with 6 main methods
- Essential utilities inlined (~200 lines)
- No complex abstractions or indirection

To add new features, simply extend the `Phase2Pipeline` class.
