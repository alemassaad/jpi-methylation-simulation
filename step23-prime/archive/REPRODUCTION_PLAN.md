# Step23 vs Step23-Prime Reproduction Plan

## Summary
This document details exactly how to reproduce step23 results using step23-prime to validate the refactoring.

## Key Parameters from Original Step23 Run

From `/step23/data/rate_0.005000/pipeline_checkpoint.json`:
- **Rate**: 0.005
- **Quantiles**: 10 (deciles)
- **Cells per quantile**: 3
- **Total individuals**: 30 per group (10 quantiles × 3 cells)
- **Growth years**: 10 (1→1024 cells)
- **Mix ratio**: 80% (80% year 60, 20% grown)
- **Seed**: 42
- **Bins**: 200 (for JSD histogram)

## Input Simulation
- **File**: `/step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz`
- **Cells**: 10,000
- **Sites**: 1,000 CpG sites per cell
- **Years**: 100 years total

## Expected Output Structure

```
step23-prime/data/rate_0.005000/
├── snapshots/
│   ├── year50_snapshot.json.gz  # 10,000 cells
│   └── year60_snapshot.json.gz  # 10,000 cells
├── individuals/
│   ├── mutant/
│   │   └── individual_00-29.json.gz  # 30 files, 5120 cells each
│   ├── control1/
│   │   └── individual_00-29.json.gz  # 30 files, 5120 cells each
│   └── control2/
│       └── individual_00-29.json.gz  # 30 files, 5120 cells each
├── plots/
│   ├── year50_jsd_distribution_200bins.png
│   └── cell_level_jsd_distributions.png
└── results/
    ├── analysis_results.json
    ├── statistics.json
    └── jsd_comparison.png
```

## Exact Commands to Run

### Step 1: Run Step23-Prime with Original Step1 Data

```bash
cd step23-prime

# Run with exact same parameters as original step23
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --growth-years 10 \
    --mix-ratio 80 \
    --seed 42 \
    --bins 200 \
    --output-dir data
```

### Step 2: Compare Results

The results should be identical because:

1. **Same random seed (42)**: Ensures same cells sampled
2. **Same quantile logic**: 10 deciles, 3 cells each
3. **Same growth**: 10 years (1→2→4→...→1024 cells)
4. **Same mixing**: 80% year 60 cells (1024→5120 total)
5. **Same input data**: Original step1 simulation

## Key Implementation Differences (Internal Only)

While results should be identical, the implementation differs:

### Step23 (Original)
- Works with dictionaries throughout
- Converts dict→Cell for operations, then Cell→dict for storage
- Manual division loops
- Complex file I/O with repeated conversions

### Step23-Prime (Refactored)
- Works with Cell/PetriDish objects natively
- Single dict→Cell conversion on load, Cell→dict on save
- Uses PetriDish.divide_cells() and methylate_cells()
- Cleaner separation of concerns

## Validation Checklist

After running both pipelines, verify:

- [ ] Same number of files in each directory (30 per group)
- [ ] Same cell counts (5120 per individual after mixing)
- [ ] Same mean JSD values for each group
- [ ] Same statistical test results (p-values)
- [ ] Visual comparison of plots

## Expected Results (from original step23)

Based on the typical output:
- **Mutant mean JSD**: ~0.76 ± 0.12
- **Control1 mean JSD**: ~0.75 ± 0.12  
- **Control2 mean JSD**: ~0.76 ± 0.12
- **P-values**: Generally >0.05 (no significant differences)

## Troubleshooting

If results differ:

1. **Check seed usage**: Ensure seed is set before each sampling operation
2. **Check quantile boundaries**: Verify same cells selected from each quantile
3. **Check growth logic**: Confirm division→age order matches
4. **Check mixing**: Verify same year 60 cells selected for mixing

## Running the Comparison

```bash
# Run original step23 (if not already done)
cd step23
python run_pipeline_v2.py --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz

# Run new step23-prime
cd ../step23-prime  
python run_pipeline.py --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz

# Compare results
diff -u ../step23/data/rate_0.005000/results/statistics.json \
        data/rate_0.005000/results/statistics.json
```

## Notes

- The `dict_to_cell()` function in step23-prime handles the conversion from step1's JSON format
- Both pipelines cache year 50/60 snapshots for efficiency
- Random seeds must be handled carefully for reproducibility
- Cell order in files may differ but statistics should match