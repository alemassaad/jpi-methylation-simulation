# DNA Methylation Simulation

A biologically realistic simulation framework for modeling DNA methylation patterns and epigenetic drift over time. This project uses object-oriented design to simulate how cells accumulate methylation through growth, division, and homeostasis.

## 🚨 Breaking Changes (Latest)

- **New Lean JSON Format**: ~90% file size reduction, no backward compatibility with old format
- **Config File Support**: YAML configuration files now preferred over complex CLI commands
- **Renamed Flag**: `calculate_jsds` → `calculate_cell_jsds` for clarity
- **New Defaults**: Optimized for faster testing (100 sites, 50 years, 64 cells)

## Overview

The simulation models:
- Single-cell origin with exponential growth and homeostasis
- Stochastic methylation accumulation (epigenetic drift)
- Cell division with methylation pattern inheritance
- Population dynamics and steady-state maintenance
- Quantile-based stratified sampling
- Mixed population experiments
- Jensen-Shannon divergence (cell_JSD) from baseline patterns

## Pipeline: Phase 1 → Phase 2

The simulation and analysis pipeline features:
- **Biological Realism**: Single cell → exponential growth → homeostasis
- **Clean Architecture**: Object-oriented Cell and PetriDish classes
- **Full Reproducibility**: Comprehensive random seeding
- **Hierarchical Organization**: Parameter-based directory structure with timestamps
- **Dynamic Configuration**: Flexible snapshot years and growth parameters

## Installation

1. Clone the repository:
```bash
git clone <repository-url>
cd jpi-methylation-simulation
```

2. Install dependencies:
```bash
pip install -r requirements.txt
```

This installs:
- `plotly` (>=5.0.0, <6.0.0) - For data visualization
- `kaleido` (0.2.1) - For PNG export
- `scipy` - For statistical analysis
- `numpy` - For numerical operations
- `pyyaml` - For configuration file support

## Quick Start

### Phase 1: Biologically Realistic Simulation
```bash
cd phase1

# Using config file (recommended)
python run_simulation.py --config configs/production.yaml

# Or with CLI arguments
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Quick test with new defaults (100 sites, 50 years, 64 cells)
python run_simulation.py --rate 0.005

# Gene-specific methylation rates
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006,50:0.007" --gene-size 5

# Save uncompressed for inspection
python run_simulation.py --rate 0.005 --no-compress

# Output will be in hierarchical structure:
# data/gene_rates_200x0.00500/size8192-sites1000-genesize5-years100-seed42-YYYYMMDDHHMMSS/simulation.json.gz
# Or: data/gene_rates_50x0.004_50x0.005.../size8192-sites1000-genesize5-years100-seed42-YYYYMMDDHHMMSS/simulation.json.gz
```

### Phase 2: Analysis Pipeline
```bash
cd phase2

# Using config file (recommended)
python run_pipeline.py --config configs/standard.yaml

# Standard analysis (30 individuals from 10 deciles)
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/gene_rates_*/size8192-*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 \
    --individual-growth-phase 7 --seed 42

# With gene-specific rates
python run_pipeline.py --gene-rate-groups "50:0.004,50:0.005,50:0.006,50:0.007" \
    --simulation ../phase1/data/gene_rates_*/simulation.json.gz \
    --first-snapshot 50 --second-snapshot 60 --seed 42

# Quick test (4 individuals from quartiles, short growth)
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/gene_rates_*/size8192-*/simulation.json.gz \
    --first-snapshot 10 --second-snapshot 15 \
    --individual-growth-phase 3 \
    --n-quantiles 4 --cells-per-quantile 1

# Output structure:
# data/gene_rates_200x0.00500-size8192-sites1000-years100/
#   └── snap50to60-growth7-quant10x3-mix80-seed42-XXXX/
#       ├── snapshots/     # Cached year 50 and 60 extracts
#       ├── individuals/   # Mutant, Control1, Control2 populations
#       ├── plots/         # cell_JSD distributions and comparisons
#       └── results/       # Statistical analyses
#
# With advanced features:
#   mix80u   - Uniform mixing
#   mix80n   - Size normalization
#   mix80un  - Both features
```

## Detailed Usage

### Phase 1: Biologically Realistic Simulation

Simulates cellular population dynamics starting from a single cell:

```bash
cd phase1

# Standard simulation using config
python run_simulation.py --config configs/production.yaml

# Or with CLI parameters
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Parameters:
# --config: Path to YAML configuration file
# --rate: Uniform methylation rate (0.005 = 0.5%)
# --gene-rate-groups: Gene-specific rates "n1:rate1,n2:rate2,..."
# --years: Total simulation time (default 50)
# --growth-phase: Years of exponential growth (default 6 → 64 cells)
# --seed: Random seed (default 42, use -1 for no seed)
# --sites: Number of CpG sites per cell (default 100)
# --gene-size: Sites per gene (default 5)
# --no-compress: Save as .json instead of .json.gz
# --no-cell-jsds: Disable cell JSD calculations
# --no-gene-jsd: Disable gene JSD tracking
# --no-jsds: Disable ALL JSD calculations

# Examples with different population sizes:
python run_simulation.py --growth-phase 10  # 1024 cells (2^10)
python run_simulation.py --growth-phase 13  # 8192 cells (2^13) - default
python run_simulation.py --growth-phase 15  # 32768 cells (2^15)
```

#### Growth Phase vs Steady State

**Growth Phase (Years 1 to growth-phase):**
- Population doubles each year through cell division
- All cells divide synchronously
- Deterministic population size: exactly 2^year cells
- Models embryonic/developmental growth

**Steady State (Years growth-phase+1 onwards):**
- Population maintained around target size
- Each year: cells divide (double), then ~50% randomly culled
- Stochastic population size varies around 2^growth-phase
- Models adult tissue homeostasis

### Phase 2: Analysis Pipeline

Complete pipeline from cell sampling to analysis:

```bash
cd phase2

# Full run with all parameters
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../phase1/data/gene_rates_*/size8192-*/simulation.json.gz \
    --first-snapshot 50 \           # First snapshot year
    --second-snapshot 60 \          # Second snapshot year  
    --individual-growth-phase 7 \   # Growth before homeostasis (2^7=128 cells)
    --n-quantiles 10 \              # Number of quantiles for sampling
    --cells-per-quantile 3 \        # Cells sampled per quantile
    --mix-ratio 80 \                # % of snapshot cells in final mix
    --uniform-mixing \              # Use same snapshot cells for all (optional)
    --normalize-size \              # Normalize to same size before mixing (optional)
    --bins 200 \                    # Histogram bins for cell_JSD plot
    --seed 42                       # Random seed for reproducibility
```

**Pipeline Stages:**

1. **Extract First Snapshot** (e.g., year 50)
   - Cached in `snapshots/year{N}_snapshot.json.gz`
   
2. **Plot cell_JSD Distribution**
   - Step histogram with statistics overlay
   - Mean, Median, SD, CV, MAD, 5%, 95%
   
3. **Create Initial Individuals**
   - Mutant: Quantile-based sampling (e.g., 3 from each decile)
   - Control1: Uniform random sampling
   
4. **Grow Individuals**
   - Exponential growth for individual-growth-phase years
   - Then homeostasis (divide → cull → methylate) maintaining ~2^growth cells
   - More biologically realistic than pure exponential growth
   
5. **Extract Second Snapshot** (first_snapshot + growth_years)
   - Age-aligned with grown individuals
   
6. **Mix Populations**
   - Add snapshot cells to grown individuals
   - Default: 80% snapshot, 20% grown
   
7. **Create Control2**
   - Pure second snapshot cells
   - Same total count as mixed populations
   
8. **Analysis & Visualization**
   - Scatter plots with statistical overlays
   - T-tests between all group pairs
   - Cell-level and individual-level distributions

## Key Features

### Hierarchical Directory Structure
```
phase1 output:
data/gene_rates_200x0.00500/size8192-sites1000-genesize5-years100-seed42-YYYYMMDDHHMMSS/simulation.json.gz

phase2 output:
data/gene_rates_200x0.00500-size8192-sites1000-years100/
  └── snap50to60-growth7-quant10x3-mix80[u][n]-seed42-XXXX/
      ├── snapshots/
      ├── individuals/
      ├── plots/
      └── results/
      
Suffix indicators:
  mix80   - Standard mixing (independent sampling)
  mix80u  - Uniform mixing (shared snapshot cells)
  mix80n  - Size normalization (median - 0.5σ threshold)
  mix80un - Both uniform mixing and normalization
```

### Full Reproducibility
- Global random seeding at pipeline start
- NumPy seeding for all array operations
- Deterministic file I/O and sampling
- Timestamp-based unique directory names

### Object-Oriented Design
```python
# Cell class - individual cell methylation
cell = Cell(n=1000, rate=0.005)
cell.methylate()  # Apply stochastic methylation
daughter = cell.create_daughter_cell()  # Mitosis

# PetriDish class - population management
petri = PetriDish(rate=0.005, growth_phase=13)
petri.divide_cells()  # Population doubles
petri.methylate_cells()  # Apply methylation
petri.random_cull_cells()  # Homeostatic culling
```

## Advanced Features (Phase 2)

### Uniform Mixing (`--uniform-mixing`)
- All individuals receive the exact same set of snapshot cells
- Reduces inter-individual variation from random sampling
- Useful for focusing on biological variation rather than sampling noise
- Directory suffix: 'u' (e.g., mix80u)

### Size Normalization (`--normalize-size`)
- Normalizes all individuals to the same cell count before mixing
- Uses median - 0.5σ threshold to determine target size
- Excludes individuals below threshold, trims those above
- Addresses variation from homeostasis stochasticity
- Typically retains ~67% of individuals
- Directory suffix: 'n' (e.g., mix80n)
- Works best combined with uniform mixing (mix80un)

### Homeostasis in Individual Growth
- `--individual-growth-phase`: Years of exponential growth before homeostasis
- After growth phase: divide → cull → methylate each year
- Maintains population around 2^growth_phase cells
- More biologically realistic than pure exponential growth
- Example: `--individual-growth-phase 7` → ~128 cells

## Project Structure

```
jpi-methylation-simulation/
├── phase1/               # Simulation engine
│   ├── cell.py               # Core classes (Cell, PetriDish)
│   ├── run_simulation.py     # Main CLI runner
│   ├── plot_history.py       # Time-series visualization
│   ├── tests/                # Unit tests
│   └── data/                 # Hierarchical output structure
├── phase2/              # Analysis pipeline
│   ├── run_pipeline.py       # Pipeline orchestrator
│   ├── pipeline_utils.py     # Cell/PetriDish utilities
│   ├── pipeline_analysis.py  # Visualization functions
│   ├── path_utils.py         # Path parsing/generation
│   ├── tests/                # Pipeline tests
│   ├── tools/                # Comparison and analysis tools
│   └── data/                 # Hierarchical output structure
├── requirements.txt           # Python dependencies
├── README.md                  # This file
└── CLAUDE.md                  # AI assistant guidance
```

## Key Concepts

### CpG Sites
Cytosine-guanine dinucleotides where methylation occurs in DNA. Initially unmethylated (0), they can become methylated (1) with a certain probability each year.

### Gene-Level Distribution
CpG sites are grouped into genes. The distribution tracks how many genes have 0, 1, 2, ... up to GENE_SIZE methylated sites.

### Jensen-Shannon Divergence (cell_JSD)
A symmetric measure of the difference between the cell's methylation distribution and a baseline distribution. Used to quantify how far a cell has diverged from the expected pattern.

### Epigenetic Drift
The gradual accumulation of methylation changes over time, modeled through stochastic methylation events.

## Mathematical Background

The simulation uses:
- **Kullback-Leibler (KL) divergence**: Measures difference between probability distributions
- **Jensen-Shannon (JS) divergence**: Symmetric version of KL divergence, calculated as the average of KL divergences to the midpoint distribution

## Testing

### Phase 1 Tests
```bash
cd phase1/tests
python test_small.py           # Quick validation
python test_comprehensive.py   # Full feature tests
python test_edge_cases.py      # Edge case handling
python test_gene_jsd.py        # Gene JSD functionality
python test_key_consistency.py # Dictionary key consistency
python test_rate_consistency.py # Rate consistency validation
```

### Phase 2 Tests
```bash
cd phase2/tests
python test_validation.py              # Data validation functions
python test_cleanup_simple.py          # Test cleanup procedures
python test_static_petridish.py        # PetriDish functionality
python test_rate_consistency_phase2.py # Rate consistency in phase2
python test_gene_jsd_extraction.py     # Gene JSD extraction tests

# Test configs
cd phase2/tests/config
python test_config_simple.py           # Simple config tests
python test_config_phase2.py           # Phase2 config tests
python test_pipeline_with_config.py    # Pipeline with config
```

## License

[Add license information here]