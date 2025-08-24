# Gene-Specific Methylation Simulation Usage Guide

## Overview
This guide explains how to use the new gene-specific methylation rates feature alongside the existing uniform rate functionality.

## Key Concepts

### Uniform vs Gene-Specific Rates
- **Uniform Rate**: All genes methylate at the same rate (original behavior)
- **Gene-Specific Rates**: Different gene groups can have different methylation rates (new feature)

### How Gene-Specific Rates Work
1. You define groups of genes with their methylation rates
2. Format: `"n1:rate1,n2:rate2,..."` where n = number of genes
3. Total genes must equal total_sites / gene_size
4. Each gene group methylates independently at its specified rate

## Phase 1: Running Simulations

### Basic Examples

#### 1. Uniform Rate (Traditional)
```bash
cd phase1

# Standard simulation with uniform 0.5% methylation rate
python run_simulation.py --rate 0.005 --years 100 --seed 42

# Quick test with fewer years and smaller population
python run_simulation.py --rate 0.005 --years 20 --growth-phase 4 --seed 42
```

#### 2. Gene-Specific Rates (New)
```bash
cd phase1

# Example 1: Four groups with increasing rates
# 50 genes at 0.4%, 50 at 0.45%, 50 at 0.5%, 50 at 0.55%
python run_simulation.py \
    --gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055" \
    --years 100 --seed 42

# Example 2: Two groups with very different rates
# 100 genes at 0.1%, 100 genes at 1.0% (10x difference)
python run_simulation.py \
    --gene-rate-groups "100:0.001,100:0.01" \
    --years 100 --seed 42

# Example 3: Many small groups with gradient
# 10 groups of 20 genes each, rates from 0.1% to 1.0%
python run_simulation.py \
    --gene-rate-groups "20:0.001,20:0.002,20:0.003,20:0.004,20:0.005,20:0.006,20:0.007,20:0.008,20:0.009,20:0.01" \
    --years 100 --seed 42
```

### Advanced Options

```bash
# Track gene-level JSD evolution
python run_simulation.py \
    --gene-rate-groups "50:0.002,50:0.004,50:0.006,50:0.008" \
    --years 100 --track-gene-jsd --seed 42

# Disable JSD calculations for maximum performance
python run_simulation.py \
    --gene-rate-groups "50:0.002,50:0.004,50:0.006,50:0.008" \
    --years 100 --no-jsds --seed 42

# Custom parameters
python run_simulation.py \
    --gene-rate-groups "40:0.003,40:0.006,40:0.009,40:0.012,40:0.015" \
    --sites 1000 \           # Total CpG sites (default: 1000)
    --gene-size 5 \          # Sites per gene (default: 5)
    --growth-phase 13 \      # Years of exponential growth (default: 13)
    --years 100 \            # Total simulation years
    --seed 42
```

### Output Structure

Data is saved in hierarchical directories:

**Uniform rate:**
```
phase1/data/
└── rate_0.00500/
    └── grow13-sites1000-years100-seed42-XXXX/
        ├── simulation.json.gz      # Complete cell history
        ├── jsd_history.png         # JSD over time plot
        └── methylation_history.png # Methylation over time plot
```

**Gene-specific rates:**
```
phase1/data/
└── gene_rates_50x0.00400_50x0.00450_50x0.00500_50x0.00550/
    └── grow13-sites1000-years100-seed42-XXXX/
        ├── simulation.json.gz
        ├── jsd_history.png
        └── methylation_history.png
```

## Phase 2: Analysis Pipeline

The phase2 pipeline works with both uniform and gene-specific rate simulations:

```bash
cd phase2

# Analyze uniform rate simulation
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --seed 42

# Analyze gene-specific rate simulation
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/gene_rates_*/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --seed 42

# With individual growth tracking
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/*/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --plot-individuals \
    --seed 42
```

## Understanding Gene-Specific Rates

### Biological Motivation
Different genomic regions methylate at different rates due to:
- Chromatin accessibility
- Transcriptional activity  
- Sequence context
- Proximity to regulatory elements

### Parameter Selection

#### Example 1: Modeling CpG Islands vs Non-Islands
```bash
# 30% of genes are CpG islands (low methylation)
# 70% of genes are non-islands (higher methylation)
--gene-rate-groups "60:0.001,140:0.007"
```

#### Example 2: Gradual Rate Variation
```bash
# Smooth gradient across genome
--gene-rate-groups "25:0.002,25:0.003,25:0.004,25:0.005,25:0.006,25:0.007,25:0.008,25:0.009"
```

#### Example 3: Hotspots and Coldspots
```bash
# Most genes normal, some hotspots
--gene-rate-groups "150:0.005,30:0.001,20:0.02"
```

### Validation Rules
1. **Cannot mix**: Can't specify both `--rate` and `--gene-rate-groups`
2. **Must specify one**: Either uniform or gene-specific required
3. **Gene count must match**: Total genes = n_sites / gene_size
   - Default: 1000 sites / 5 sites per gene = 200 genes total
   - Your groups must sum to exactly 200 genes

### Error Examples
```bash
# ERROR: Both rate types specified
python run_simulation.py --rate 0.005 --gene-rate-groups "100:0.004,100:0.006"

# ERROR: Gene count mismatch (150 ≠ 200)
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006"

# ERROR: Invalid format
python run_simulation.py --gene-rate-groups "50-0.004,50-0.005"  # Wrong separator
```

## Workflow Examples

### Workflow 1: Compare Uniform vs Variable Rates

```bash
# Step 1: Run uniform rate simulation
cd phase1
python run_simulation.py --rate 0.005 --years 100 --seed 42

# Step 2: Run variable rate simulation with same average
python run_simulation.py \
    --gene-rate-groups "50:0.003,50:0.004,50:0.006,50:0.007" \
    --years 100 --seed 42

# Step 3: Analyze both
cd ../phase2
# Analyze uniform
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 --seed 42

# Analyze variable  
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/gene_rates_*/*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 --seed 42

# Step 4: Compare results
cd tools
python compare_two_runs.py \
    --dir1 ../data/rate_0.00500*/snap50to60*/ \
    --dir2 ../data/gene_rates*/snap50to60*/
```

### Workflow 2: Explore Rate Heterogeneity Effects

```bash
cd phase1

# Low heterogeneity (rates: 0.4% to 0.6%)
python run_simulation.py \
    --gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055" \
    --years 100 --seed 42 --track-gene-jsd

# High heterogeneity (rates: 0.1% to 1.0%)  
python run_simulation.py \
    --gene-rate-groups "50:0.001,50:0.003,50:0.007,50:0.010" \
    --years 100 --seed 42 --track-gene-jsd

# Extreme heterogeneity (rates: 0.01% to 2.0%)
python run_simulation.py \
    --gene-rate-groups "50:0.0001,50:0.001,50:0.01,50:0.02" \
    --years 100 --seed 42 --track-gene-jsd
```

### Workflow 3: Quick Testing

```bash
cd phase1

# Quick test with small population and few years
python run_simulation.py \
    --gene-rate-groups "50:0.01,50:0.02,50:0.03,50:0.04" \
    --years 10 --growth-phase 3 --seed 42

# Verify it worked
ls -la data/gene_rates_*/grow3-sites1000-years10-seed42-*/
```

## Testing the Implementation

```bash
# Run comprehensive tests
cd phase1/tests
python test_comprehensive.py  # Includes gene rate tests
python test_edge_cases.py     # Tests edge cases
python test_gene_jsd.py       # Tests gene JSD with variable rates

# Run a specific test
cd phase1
python -c "
from cell import Cell, PetriDish
groups = [(50, 0.01), (50, 0.02), (50, 0.03), (50, 0.04)]
cell = Cell(n=1000, gene_rate_groups=groups, gene_size=5)
print(f'Created cell with {len(cell.site_rates)} site rates')
print(f'First 5 rates: {cell.site_rates[:5]}')
print(f'Rates at positions 250, 500, 750: {[cell.site_rates[i] for i in [250, 500, 750]]}')
"
```

## Tips and Best Practices

1. **Start small**: Test with `--years 10 --growth-phase 3` before long runs
2. **Use consistent seeds**: Same seed allows direct comparison
3. **Track gene JSD**: Use `--track-gene-jsd` to monitor gene-level heterogeneity
4. **Balance groups**: Similar group sizes often work better than extreme differences
5. **Document parameters**: Save your command in a script or README

## Troubleshooting

### Common Issues

**Issue**: "Cannot specify both rate and gene-rate-groups"
- **Solution**: Use either `--rate` OR `--gene-rate-groups`, not both

**Issue**: "gene-rate-groups specifies X genes, but simulation has Y genes"
- **Solution**: Ensure gene counts sum to n_sites/gene_size (default: 200)

**Issue**: Output directory name too long
- **Solution**: Gene rate directory names are truncated at 50 characters

**Issue**: Memory usage with large simulations
- **Solution**: Use `--no-jsds` flag to disable JSD calculations

## Next Steps

1. Start with the quick test examples
2. Run a comparison between uniform and variable rates
3. Explore how different rate distributions affect outcomes
4. Use phase2 pipeline to analyze population heterogeneity

## Questions?

Check the main documentation:
- `CLAUDE.md`: Technical architecture details
- `README.md`: Project overview
- `gene_rate_implementation_plan.md`: Implementation details