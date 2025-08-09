# Step23-Prime Pipeline

A refactored version of the Step23 pipeline using object-oriented design with `PetriDish` and `Cell` classes from step1-prime.

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
    --simulation ../step1-prime/data/rate_0.00500/grow13-sites1000-years100-seed42-*/simulation.json.gz \
    --snapshot-year 50 \
    --growth-years 10 \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --mix-ratio 80 \
    --seed 42
```

## Parameters

- `--rate`: Methylation rate (must match simulation)
- `--simulation`: Path to step1/step1-prime simulation file
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

## Comparison with Original Step23

To validate results against the original step23:

```bash
python tools/generate_comparison_report.py
```

This compares snapshots, individuals, and statistical results between step23 and step23-prime.