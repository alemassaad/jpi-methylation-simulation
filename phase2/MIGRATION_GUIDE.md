# Pipeline Migration Guide

## Overview
This guide covers two major migrations:
1. **Modular Phase 2** (January 20, 2025): Reorganized phase2 into modular scripts
2. **Three-Phase Architecture** (September 20, 2025): Split analysis into separate phase3

This guide helps you migrate from old workflows to the new three-phase architecture.

## ⚠️ Breaking Changes

### 1. **Three-Phase Architecture (NEW - Sep 2025)**
- **Phase 2** now only generates data (NO plotting/analysis)
- **Phase 3** handles ALL analysis and visualization
- **Workflow**: Must run phase1 → phase2 → phase3 for complete analysis
- **No backward compatibility** with old two-phase workflows

### 2. **Modular Phase 2 (Jan 2025)**
- Phase 2 reorganized into modular scripts
- `run_pipeline.py` is now a driver for data generation only
- Analysis moved to separate phase3

### 3. **Command-Line Arguments**
Changes from old to new architecture:
- `--plot-individuals` → Removed (moved to phase3)
- `--bins` → Moved to phase3
- `--max-gene-plots` → Moved to phase3
- `--rate` → Auto-detected from simulation (phase2)
- `--gene-rate-groups` → Auto-detected from simulation (phase2)

## Migration Steps

### From Old Two-Phase to New Three-Phase

If you were running the old complete pipeline:
```bash
# OLD (single command for everything)
cd phase2
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 --second-snapshot 50
```

Now you need TWO commands:
```bash
# NEW (data generation)
cd phase2
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 --second-snapshot 50
# Rate auto-detected from simulation!

# NEW (analysis and visualization)
cd ../phase3
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz
```

### For Advanced Users

If you were using advanced features:
```bash
# OLD (everything in one command)
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 --second-snapshot 50 \
    --uniform-mixing --normalize-size \
    --plot-individuals --bins 200 --no-compress
```

Now split into two phases:
```bash
# NEW Phase 2 (data generation)
cd phase2
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 --second-snapshot 50 \
    --uniform-mixing --normalize-size \
    --no-compress

# NEW Phase 3 (analysis with custom settings)
cd ../phase3
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --bins 200 --max-gene-plots 10
```

### For Config File Users

Now you need separate config files for phase2 and phase3:

```yaml
# Phase 2 config (data generation)
# configs/my_phase2_config.yaml
simulation: ../phase1/data/*/simulation.json.gz
first_snapshot: 30
second_snapshot: 50
n_quantiles: 10
cells_per_quantile: 3
uniform_mixing: true
normalize_size: true
```

```yaml
# Phase 3 config (analysis)
# configs/my_phase3_config.yaml
bins: 200
max_gene_plots: 15
# output_dir: my_analysis
```

Usage:
```bash
# Run phase2 with config
cd phase2
python run_pipeline.py --config configs/my_phase2_config.yaml

# Run phase3 with config
cd ../phase3
python run_analysis.py \
    --phase2-dir ../phase2/data/{run}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --config configs/my_phase3_config.yaml
```

### For Script Automation

Update automation scripts for three-phase workflow:

```bash
# OLD automation (single phase)
for sim in ../phase1/data/*/simulation.json.gz; do
    python run_pipeline.py --simulation "$sim" ...
done

# NEW automation (two-step process)
# Step 1: Generate data with phase2
for sim in ../phase1/data/*/simulation.json.gz; do
    echo "Processing simulation: $sim"
    cd phase2
    python run_pipeline.py --simulation "$sim" --first-snapshot 30 --second-snapshot 50
    cd ..
done

# Step 2: Analyze all generated datasets with phase3
for phase2_dir in phase2/data/*/; do
    echo "Analyzing: $phase2_dir"
    cd phase3
    python run_analysis.py \
        --phase2-dir "../$phase2_dir" \
        --simulation ../phase1/data/*/simulation.json.gz \
        --config configs/default.yaml
    cd ..
done
```

## New Features

### 1. **Complete Phase Separation**
Clean separation of data generation and analysis:
- **Phase 2**: Only generates datasets (no plotting)
- **Phase 3**: Only analyzes and visualizes (no data generation)
- **Re-analysis**: Run phase3 multiple times on same phase2 data

### 2. **Modular Execution (Phase 2)**
You can run phase2 stages independently:

```bash
# Extract snapshots once, reuse many times
python extract_snapshots.py --simulation ../phase1/data/*/simulation.json.gz \
    --output-dir data/my_analysis

# Try different parameters without re-extracting
python simulate_individuals.py --base-dir data/my_analysis \
    --n-quantiles 10 --cells-per-quantile 3

# Create control populations
python create_control2.py --base-dir data/my_analysis
```

### 3. **Flexible Analysis (Phase 3)**
Analyze same data multiple ways:

```bash
# Quick analysis
python run_analysis.py ... --config configs/quick_analysis.yaml

# High-resolution analysis
python run_analysis.py ... --bins 400 --max-gene-plots 20

# Custom output location
python run_analysis.py ... --output-dir publication_plots
```

### 4. **Snapshot Caching**
Snapshots are cached and reused:
- First run extracts and saves snapshots
- Subsequent runs reuse cached snapshots
- Use `--force-reload` to re-extract

### 5. **Batch Processing**
Built-in support for analyzing multiple phase2 runs:

```bash
# Analyze all phase2 runs
for dir in ../phase2/data/*/; do
    python run_analysis.py --phase2-dir "$dir" --simulation ../phase1/data/*/simulation.json.gz
done
```

### 6. **Metadata Files**
Enhanced metadata for debugging:
- `snapshots/metadata.json`: Extraction parameters
- `individuals/mixing_metadata.json`: Mixing configuration
- `results/metadata/analysis_metadata.json`: Analysis settings

## Troubleshooting

### Q: "Phase2 directory not found" errors?
Make sure phase2 completed successfully:
```bash
# Check phase2 output
ls ../phase2/data/
# Should see directories like: snap30to50-growth7-quant10x3-mix80-seed42-YYYYMMDDHHMMSS/
```

### Q: "Module not found" errors?
Ensure you're in the correct directory:
```bash
# For phase2
cd phase2
python run_pipeline.py ...

# For phase3
cd phase3
python run_analysis.py ...
```

### Q: Old output directories?
Old output directories are incompatible. Archive them:
```bash
# Archive old two-phase outputs
mv phase2/data/old_runs phase2/data/archive_pre_20250920

# Archive old modular outputs (if any)
mv phase2/data/old_results phase2/data/archive_pre_modular
```

### Q: Config file not working?
Check for moved parameters:
- **Phase2 configs**: Remove `bins`, `max_gene_plots`, `plot_individuals`
- **Phase3 configs**: Remove `first_snapshot`, `n_quantiles`, etc.
- Update paths to be relative to correct phase directory

### Q: Performance different?
- **Phase2**: May be slightly slower due to modular architecture
- **Phase3**: Much faster since it only does analysis
- **Overall**: More efficient for iterative analysis (run phase3 multiple times)

### Q: How to find phase2 output directory?
```bash
# List recent phase2 runs
ls -lt ../phase2/data/

# Or use wildcard pattern
python run_analysis.py --phase2-dir ../phase2/data/snap*/ --simulation ...
```

## Benefits of Migration

1. **Separation of Concerns**: Data generation and analysis are completely separate
2. **Efficiency**: Re-analyze without expensive re-simulation
3. **Flexibility**: Multiple analysis approaches on same dataset
4. **Modularity**: Debug individual stages easily
5. **Caching**: Reuse expensive computations
6. **Scalability**: Batch process multiple datasets
7. **Maintainability**: Cleaner, more organized code
8. **Development**: Independent development of analysis features

## Complete Migration Example

```bash
# Complete workflow: phase1 → phase2 → phase3

# Step 1: Generate simulation (phase1)
cd phase1
python run_simulation.py --config configs/production.yaml

# Step 2: Generate datasets (phase2)
cd ../phase2
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --config configs/quick_test.yaml

# Step 3: Analyze and visualize (phase3)
cd ../phase3
python run_analysis.py \
    --phase2-dir ../phase2/data/{latest_run}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --config configs/default.yaml
```

## Need Help?

1. **Phase2**: Check `phase2/README.md` for data generation examples
2. **Phase3**: Check `phase3/README.md` for analysis examples
3. **Commands**: Run `--help` for all options:
   - `python phase2/run_pipeline.py --help`
   - `python phase3/run_analysis.py --help`
4. **Configs**: Look at example configs in each phase's `configs/` directory
5. **Issues**: File bugs on GitHub with specific phase information

## Rollback Option

If you need the old pipeline temporarily:
```bash
# Old monolithic pipeline backed up as:
cd phase2
python run_pipeline_old.py ...
# But this is deprecated and will be removed!
```

**Note**: The old pipeline is deprecated and will be removed in future versions. Please migrate to the three-phase architecture as soon as possible.

## Summary

- **Old**: Single command did everything
- **New**: phase1 → phase2 → phase3 workflow
- **Benefits**: Re-analysis, flexibility, separation of concerns
- **Cost**: Slightly more complex workflow, but much more powerful