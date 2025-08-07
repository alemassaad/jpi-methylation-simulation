# Step 3 Batch Processor

Automates the complete Step 3 pipeline for multiple methylation rates.

## What it does

1. **Auto-discovers** completed rates from Step 2
2. **Extracts year 60** from original simulations
3. **Creates mixed populations** (80% original + 20% lineage cells)
4. **Analyzes differences** between mutant and control groups
5. **Performs cross-rate analysis** to see trends

## Quick Start

```bash
cd step3/scripts

# Process all available rates
python batch_processor_step3.py

# Process specific rates only
python batch_processor_step3.py --rates 0.003000,0.005000

# Dry run to see what would happen
python batch_processor_step3.py --dry-run

# Resume if interrupted
python batch_processor_step3.py --resume
```

## Output Structure

```
step3/data/
├── rate_0.003000/
│   ├── snapshots/
│   │   └── year60_snapshot.json.gz
│   ├── individuals/
│   │   ├── mutant/     (30 mixed individuals)
│   │   └── control/    (30 mixed individuals)
│   ├── plots/
│   │   └── jsd_distribution_comparison.png
│   └── results/
│       └── jsd_distributions.json
├── rate_0.005000/
│   └── ... (same structure)
├── combined_analysis/
│   ├── plots/
│   │   ├── cross_rate_analysis.png
│   │   └── cross_rate_analysis.html
│   └── results/
│       └── summary_analysis.json
└── .batch_state_step3.json
```

## Features

- **Automatic discovery** of Step 2 outputs
- **Smart path resolution** - finds all required files automatically
- **Cross-rate analysis** - compares trends across different methylation rates
- **Interactive plots** - HTML outputs for detailed exploration
- **Statistical testing** - Mann-Whitney U tests for significance
- **Resume capability** - pick up where you left off

## Configuration

Edit `config/step3_config.yaml` to customize:
- Number of individuals (default: 30)
- Mixture ratio (default: 20% lineage cells)
- Statistical tests
- Plot formats

## Time Estimate

- Per rate: ~3-5 minutes
- Total for 2 rates: ~10 minutes
- Cross-rate analysis: ~1 minute

## Interpreting Results

1. **Individual plots**: Show JSD distribution differences between mutant/control
2. **Cross-rate plot**: Shows how the effect changes with methylation rate
3. **Statistical results**: P-values indicate if differences are significant

## Next Steps

After processing:
1. Check individual rate results in `data/rate_*/`
2. Review cross-rate analysis in `data/combined_analysis/`
3. Open HTML plots for interactive exploration