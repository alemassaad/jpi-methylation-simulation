# Phase 2 Migration Guide

## Overview
On January 20, 2025, Phase 2 was completely reorganized from a monolithic script into a modular architecture with 4 independent scripts. This guide helps you migrate from the old pipeline to the new one.

## ⚠️ Breaking Changes

### 1. **No Backward Compatibility**
The new architecture is NOT compatible with old pipeline runs. You cannot:
- Resume old runs with the new scripts
- Mix old and new output directories
- Use old analysis results with new code

### 2. **Main Script Changed**
- **Old**: Single monolithic `run_pipeline.py` did everything
- **New**: `run_pipeline.py` is now a driver that calls 4 separate scripts

### 3. **Command-Line Arguments**
Some arguments have been renamed or removed:
- `--plot-individuals` → Removed (always generates plots now)
- `--rate` → Still works but auto-detected from simulation
- `--gene-rate-groups` → Still works but auto-detected

## Migration Steps

### For Basic Users

If you were running:
```bash
# OLD
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60
```

Simply change to:
```bash
# NEW
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60
# Rate auto-detected from simulation!
```

### For Advanced Users

If you were using advanced features:
```bash
# OLD
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --uniform-mixing --normalize-size \
    --plot-individuals --no-compress
```

Change to:
```bash
# NEW
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --uniform-mixing --normalize-size \
    --no-compress
# Note: --plot-individuals removed (always plots now)
```

### For Config File Users

Config files remain largely compatible. Just update:

```yaml
# OLD config
simulation:
  rate: 0.005
  # ... other simulation params

# NEW config (simpler!)
# Rate is auto-detected from simulation file
first_snapshot: 30
second_snapshot: 50
n_quantiles: 10
cells_per_quantile: 3
```

### For Script Automation

If you have automation scripts, update them:

```bash
# OLD automation
for rate in 0.004 0.005 0.006; do
    python run_pipeline.py --rate $rate ...
done

# NEW automation
for sim in ../phase1/data/*/simulation.json.gz; do
    python run_pipeline.py --simulation "$sim" ...
done
```

## New Features

### 1. **Modular Execution**
You can now run stages independently:

```bash
# Extract snapshots once, reuse many times
python extract_snapshots.py --simulation ../phase1/data/*/simulation.json.gz \
    --output-dir data/my_analysis

# Try different parameters without re-extracting
python simulate_individuals.py --base-dir data/my_analysis \
    --n-quantiles 10 --cells-per-quantile 3

# Re-run analysis with different plot settings
python analyze_and_plot.py --base-dir data/my_analysis \
    --bins 100 --max-gene-plots 5
```

### 2. **Snapshot Caching**
Snapshots are now cached and reused:
- First run extracts and saves snapshots
- Subsequent runs reuse cached snapshots
- Use `--force-reload` to re-extract

### 3. **Better Error Recovery**
If pipeline fails at stage 5, you can:
1. Fix the issue
2. Re-run just that stage
3. Continue from there

### 4. **Metadata Files**
New metadata files for debugging:
- `snapshots/metadata.json`: Extraction parameters
- `individuals/mixing_metadata.json`: Mixing configuration

## Troubleshooting

### Q: "Module not found" errors?
The new scripts need proper imports:
```bash
# Make sure you're in phase2 directory
cd phase2
python run_pipeline.py ...
```

### Q: Old output directories?
Old output directories are incompatible. Archive them:
```bash
mv data/old_runs data/archive_pre_20250120
```

### Q: Config file not working?
Check for removed parameters:
- Remove `plot_individuals` (always on now)
- Remove explicit `rate` if using auto-detection
- Update paths if needed

### Q: Performance different?
The new architecture may be slightly slower due to subprocess calls, but provides better reliability and debugging capabilities.

## Benefits of Migration

1. **Modularity**: Debug individual stages easily
2. **Caching**: Reuse expensive computations
3. **Flexibility**: Custom workflows possible
4. **Maintainability**: Cleaner, more organized code
5. **Future-proof**: Ready for Phase 3 split

## Need Help?

1. Check the updated README.md for examples
2. Run `python run_pipeline.py --help` for all options
3. Look at example configs in `configs/`
4. File issues on GitHub if you find bugs

## Rollback Option

If you need the old pipeline temporarily:
```bash
# Old pipeline backed up as:
python run_pipeline_old.py ...
# But this is deprecated and will be removed!
```

**Note**: The old pipeline is deprecated and will be removed in future versions. Please migrate as soon as possible.