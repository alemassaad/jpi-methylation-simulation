# Phase 1: DNA Methylation Simulation Engine

A biologically realistic methylation simulation that models epigenetic drift through cell growth, division, and homeostasis.

## 🚨 Breaking Changes (Latest)

- **New Lean JSON Format**: ~90% file size reduction, no backward compatibility
- **Config File Support**: YAML configuration files now preferred over complex CLI commands
- **Renamed Flag**: `calculate_jsds` → `calculate_cell_jsds` for clarity
- **New Defaults**: Optimized for faster testing (100 sites, 50 years, 64 cells)

## Key Features

- **Single-cell origin**: Simulation starts with one unmethylated cell
- **Exponential growth phase**: Population doubles each year until reaching target size
- **Homeostatic maintenance**: After growth, population maintained through division and culling
- **Lean JSON format**: Efficient storage with ~90% size reduction
- **Configuration files**: YAML-based configuration system
- **Flexible methylation**: Uniform or gene-specific methylation rates
- **JSD tracking**: Both cell-level and gene-level divergence metrics
- **Compression control**: Choose between .json or .json.gz output

## Current Structure (After Cleanup)

Phase 1 now contains only core simulation components (5 files total):
- **`cell.py`**: Core Cell and PetriDish classes
- **`run_simulation.py`**: Main simulation script
- **`plot_history.py`**: Visualization tool for simulation results
- **`config_default.yaml`**: Default configuration template
- **`README.md`**: This documentation

All test files and example configs have been removed for simplicity.

## Quick Start

```bash
# Quick test with defaults (100 sites, 50 years, 64 cells)
python run_simulation.py --rate 0.005

# Standard production run
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000 --seed 42

# Gene-specific methylation rates
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006,50:0.007" --gene-size 5

# Create your own config from config_default.yaml template
python run_simulation.py --config my_config.yaml

# Visualize results
python plot_history.py data/rate_0.00500/*/simulation.json.gz
```

## Configuration

### Config File Structure

Create a YAML file (e.g., `my_config.yaml`):

```yaml
simulation:
  rate: 0.005           # Methylation rate (or use gene_rate_groups)
  sites: 1000           # CpG sites per cell
  years: 100            # Simulation duration
  gene_size: 5          # Sites per gene
  growth_phase: 13      # 2^13 = 8192 cells

output:
  directory: "data"
  compress: true        # Use .json.gz

performance:
  track_gene_jsd: true
  calculate_cell_jsds: true

seed: 42
```

### Command-Line Parameters

- `--config`: Path to YAML config file
- `--rate`: Uniform methylation rate per site per year (0.005 = 0.5%)
- `--gene-rate-groups`: Gene-specific rates (e.g., "50:0.004,50:0.006")
- `--sites`: Number of CpG sites per cell (default: 100)
- `--years`: Total simulation time (default: 50)
- `--gene-size`: Sites per gene (default: 5)
- `--growth-phase`: Years of exponential growth (default: 6 → 64 cells)
- `--seed`: Random seed (-1 for random)
- `--no-compress`: Save as .json instead of .json.gz
- `--no-cell-jsds`: Skip cell JSD calculations
- `--no-gene-jsd`: Skip gene JSD tracking
- `--no-jsds`: Skip ALL JSD calculations (maximum speed)
- `--no-save`: Run without saving results

### Performance Flags

```bash
# Default: Calculate everything
python run_simulation.py --rate 0.005

# Skip cell JSDs only (still track gene JSDs)
python run_simulation.py --rate 0.005 --no-cell-jsds

# Skip gene JSDs only (still calculate cell JSDs)
python run_simulation.py --rate 0.005 --no-gene-jsd

# Skip ALL JSDs for maximum performance
python run_simulation.py --rate 0.005 --no-jsds
```

## Output Format

### New Lean JSON Structure

```json
{
  "parameters": {
    "rate": 0.005,
    "gene_rate_groups": null,
    "n": 1000,
    "gene_size": 5,
    "growth_phase": 13,
    "years": 100,
    "seed": 42,
    "track_cell_history": true,
    "track_gene_jsd": true
  },
  "history": {
    "0": {
      "cells": [
        {
          "methylated": [0, 0, 1, 0, ...],  // Only essential data
          "cell_JSD": 0.0                   // Cell divergence
        }
      ],
      "gene_jsd": [0.0, 0.0, ...]          // Per-gene JSDs
    },
    "1": { ... }
  }
}
```

**Key improvements:**
- No redundant data per cell (rate, gene_size, site_rates removed)
- Parameters stored once at top level
- ~90% file size reduction
- Cleaner structure

### Output Directory Structure

```
data/
└── rate_0.00500/                           # Or gene_rates_50x0.004_50x0.006.../
    └── grow13-sites1000-years100-seed42-YYYYMMDDHHMMSS/  # Timestamp for uniqueness
        ├── simulation.json.gz               # Compressed (default for production)
        ├── simulation.json                  # Uncompressed (if --no-compress)
        ├── jsd_history.png                  # Cell JSD trajectory
        └── methylation_history.png          # Methylation accumulation
```

The directory name includes a timestamp (YYYYMMDDHHMMSS format) for chronological sorting and uniqueness.

## Biological Model

### Growth Phase (Years 0 to growth_phase)
- Population doubles synchronously each year
- Deterministic growth: exactly 2^year cells
- Models embryonic/developmental expansion
- All cells divide, then methylate

### Steady State (Years growth_phase+1 onwards)
- Population maintained via homeostasis
- Each year: divide (2x) → cull (~50% survival) → methylate
- Stochastic population varies around target
- Models adult tissue maintenance

### Methylation Process
- Each unmethylated site has probability `rate` of becoming methylated per year
- Once methylated, sites remain methylated (irreversible)
- Methylation patterns inherited perfectly during division
- No de novo methylation during mitosis

## Example Configurations

### Quick Test
```bash
python run_simulation.py --config configs/quick_test.yaml
# Or manually: 4 cells, 10 years, fast
python run_simulation.py --rate 0.01 --years 10 --growth-phase 2 --sites 100
```

### Production Run
```bash
python run_simulation.py --config configs/production.yaml
# Or manually: 8192 cells, 100 years, standard
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000
```

### Debug Mode
```bash
python run_simulation.py --config configs/debug.yaml
# Verbose output, uncompressed files
```

### Gene-Specific Rates
```bash
python run_simulation.py --config configs/gene_rates.yaml
# Or manually specify different rates for gene groups
python run_simulation.py --gene-rate-groups "100:0.003,100:0.007" --gene-size 5
```

## Metrics and Analysis

### Cell JSD (Jensen-Shannon Divergence)
- Measures how different each cell's methylation pattern is from baseline
- Range: 0 (identical to start) to 1 (maximally different)
- Tracked for every cell at every time point
- Useful for measuring "epigenetic age"

### Gene JSD
- Measures methylation heterogeneity across population for each gene
- Calculated per gene across all cells
- Tracks which genes show most variation
- Useful for identifying variable genomic regions

### Methylation Statistics
- Mean/median methylation fraction
- Standard deviation across cells
- Min/max values in population
- Distribution histograms in plots

## Performance

### Typical Runtimes
- **Quick test** (100 sites, 10 years, 4 cells): ~1 second
- **Standard test** (100 sites, 50 years, 64 cells): ~5 seconds
- **Production** (1000 sites, 100 years, 8192 cells): ~2-5 minutes
- **Large scale** (1000 sites, 200 years, 32768 cells): ~10-20 minutes

### File Sizes (New Lean Format)
- **Small** (100 sites, 50 years, 64 cells): 
  - Uncompressed: ~0.5 MB
  - Compressed: ~10 KB
- **Standard** (1000 sites, 100 years, 8192 cells):
  - Uncompressed: ~7 MB
  - Compressed: ~100 KB
- **Large** (1000 sites, 200 years, 32768 cells):
  - Uncompressed: ~30 MB
  - Compressed: ~400 KB

### Optimization Tips
- Use `--no-jsds` for maximum speed if JSDs not needed
- Use `--no-gene-jsd` if only cell-level metrics needed
- Use compressed output for production runs
- Smaller `sites` parameter for quick tests
- Config files load faster than parsing many CLI args

## Testing

```bash
# Run individual tests
cd tests
python test_small.py           # Quick validation
python test_comprehensive.py   # Full feature test
python test_edge_cases.py      # Edge case handling
python test_gene_jsd.py        # Gene JSD functionality

# Additional tests
python test_key_consistency.py # Dictionary key consistency
python test_rate_consistency.py # Rate consistency validation
```

## Troubleshooting

### Common Issues

**Large file sizes:**
- Ensure using new code (lean format)
- Use `--compress` or set `compress: true` in config
- Check not using legacy format

**JSD values all zero:**
- Check if `--no-cell-jsds` or `--no-jsds` was used
- Verify `calculate_cell_jsds: true` in config

**Config not loading:**
- Install PyYAML: `pip install pyyaml`
- Check YAML syntax (proper indentation)
- Verify file path is correct

**Memory issues:**
- Reduce `growth_phase` for smaller populations
- Use `--no-jsds` to save memory
- Process output in chunks for analysis

## API Usage

```python
from cell import Cell, PetriDish

# Create a petri dish
petri = PetriDish(rate=0.005, growth_phase=10, seed=42)

# Run simulation
petri.run_simulation(years=50)

# Access results
for cell in petri.cells:
    print(f"Cell methylation: {cell.methylation_proportion:.2%}")
    print(f"Cell JSD: {cell.cell_JSD:.4f}")

# Save results
petri.save_history(compress=True)
```

## Requirements

- Python 3.7+
- NumPy (optional, for performance)
- SciPy (optional, for JSD calculations)
- Plotly (for visualization)
- PyYAML (for config files)
- Kaleido (for static plot export)

Install all dependencies:
```bash
pip install numpy scipy plotly pyyaml kaleido
```