# Phase 2: Analysis Pipeline

A refactored analysis pipeline using object-oriented design with `PetriDish` and `Cell` classes from phase1.

## Features

- **OOP Design**: Uses Cell and PetriDish objects instead of dictionaries
- **Dynamic Snapshots**: Configurable snapshot years (not hardcoded to 50/60)
- **Hierarchical Output**: Organized directory structure with parameter tracking
- **Full Reproducibility**: Proper random seeding throughout
- **Backward Compatible**: Works with both old and new simulation formats

## Usage

```bash
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-sites1000-years100-seed42-*/simulation.json.gz \
    --snapshot-year 50 \
    --growth-years 10 \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --mix-ratio 80 \
    --seed 42
```

## Parameters

- `--rate`: Methylation rate (must match simulation)
- `--simulation`: Path to phase1 simulation file
- `--snapshot-year`: Year for first snapshot (default: 50)
- `--growth-years`: Years to grow individuals (default: 10)
- `--n-quantiles`: Number of quantiles for sampling (default: 10)
- `--cells-per-quantile`: Cells per quantile (default: 3)
- `--mix-ratio`: Percentage of snapshot cells in final mix (default: 80)
- `--seed`: Random seed for reproducibility (default: 42)

## Output Structure

```
data/
└── rate_0.00500-grow13-sites1000-years100/
    └── snap50-quant10x3-grow10-mix80-seed42-a3f2/
        ├── snapshots/
        │   ├── year50_snapshot.json.gz
        │   └── year60_snapshot.json.gz
        ├── individuals/
        │   ├── mutant/
        │   ├── control1/
        │   └── control2/
        ├── plots/
        └── results/
```

## Directory Organization

- **Core Files** (root): Pipeline implementation
  - `run_pipeline.py`: Main entry point
  - `pipeline_utils.py`: Core utilities
  - `pipeline_analysis.py`: Analysis functions
  - `path_utils.py`: Path parsing/generation

- **tests/**: Test files for validation
- **tools/**: Comparison and analysis tools
- **scripts/**: Utility scripts
- **archive/**: Development artifacts and documentation

## Analysis Tools

Various comparison and analysis tools are available in the `tools/` directory:

```bash
python tools/compare_two_runs.py  # Compare reproducibility between runs
python tools/generate_comparison_report.py  # Generate comprehensive comparison report
```