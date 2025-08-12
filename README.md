# DNA Methylation Simulation

A biologically realistic simulation framework for modeling DNA methylation patterns and epigenetic drift over time. This project uses object-oriented design to simulate how cells accumulate methylation through growth, division, and homeostasis.

## Overview

The simulation models:
- Single-cell origin with exponential growth and homeostasis
- Stochastic methylation accumulation (epigenetic drift)
- Cell division with methylation pattern inheritance
- Population dynamics and steady-state maintenance
- Quantile-based stratified sampling
- Mixed population experiments
- Jensen-Shannon divergence (JSD) from baseline patterns

## Pipeline: Phase 1 → Phase 2

The simulation and analysis pipeline features:
- **Biological Realism**: Single cell → exponential growth → homeostasis
- **Clean Architecture**: Object-oriented Cell and PetriDish classes
- **Full Reproducibility**: Comprehensive random seeding
- **Hierarchical Organization**: Parameter-based directory structure with MD5 hashing
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

## Quick Start

### Phase 1: Biologically Realistic Simulation
```bash
# Run simulation starting from 1 cell, growing to 8192 cells (2^13)
cd phase1
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Output will be in hierarchical structure:
# data/rate_0.00500/grow13-sites1000-years100-seed42-XXXX/simulation.json.gz
```

### Phase 2: Analysis Pipeline
```bash
cd phase2

# Standard analysis (30 individuals from 10 deciles)
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --snapshot-year 50 --growth-years 10 --seed 42

# Quick test (4 individuals from quartiles, 2 years growth)
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --n-quantiles 4 --cells-per-quantile 1 --growth-years 2

# Output structure:
# data/rate_0.00500-grow13-sites1000-years100/
#   └── snap50-quant10x3-grow10-mix80-seed42-XXXX/
#       ├── snapshots/     # Cached year 50 and 60 extracts
#       ├── individuals/   # Mutant, Control1, Control2 populations
#       ├── plots/         # JSD distributions and comparisons
#       └── results/       # Statistical analyses
```

## Detailed Usage

### Phase 1: Biologically Realistic Simulation

Simulates cellular population dynamics starting from a single cell:

```bash
cd phase1

# Standard simulation (8192 cells, 100 years)
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Parameters:
# --rate: Methylation rate per site per year (0.005 = 0.5%)
# --years: Total simulation time
# --growth-phase: Years of exponential growth (final population = 2^growth-phase)
# --seed: Random seed (default 42, use -1 for no seed)
# --sites: Number of CpG sites per cell (default 1000)
# --gene-size: Sites per gene (default 5)

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
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --snapshot-year 50 \        # First snapshot year (default: 50)
    --growth-years 10 \         # Years to grow individuals (default: 10)
    --n-quantiles 10 \          # Number of quantiles for sampling
    --cells-per-quantile 3 \    # Cells sampled per quantile
    --mix-ratio 80 \            # % of snapshot cells in final mix
    --bins 200 \                # Histogram bins for JSD plot
    --seed 42                   # Random seed for reproducibility
```

**Pipeline Stages:**

1. **Extract First Snapshot** (e.g., year 50)
   - Cached in `snapshots/year{N}_snapshot.json.gz`
   
2. **Plot JSD Distribution**
   - Step histogram with statistics overlay
   - Mean, Median, SD, CV, MAD, 5%, 95%
   
3. **Create Initial Individuals**
   - Mutant: Quantile-based sampling (e.g., 3 from each decile)
   - Control1: Uniform random sampling
   
4. **Grow Individuals**
   - Cell division + methylation for N years
   - Growth: 1→2→4→...→2^N cells
   
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
step1-prime output:
data/rate_0.00500/grow13-sites1000-years100-seed42-XXXX/simulation.json.gz

step23-prime output:
data/rate_0.00500-grow13-sites1000-years100/
  └── snap50-quant10x3-grow10-mix80-seed42-XXXX/
      ├── snapshots/
      ├── individuals/
      ├── plots/
      └── results/
```

### Full Reproducibility
- Global random seeding at pipeline start
- NumPy seeding for all array operations
- Deterministic file I/O and sampling
- MD5 hashing for unique directory names

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

### Jensen-Shannon Divergence (JSD)
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
```

### Phase 2 Tests
```bash
cd phase2/tests
python test_reproducibility_robust.py   # Test full reproducibility
python test_dynamic_mix_year.py        # Test dynamic year calculations
```

## License

[Add license information here]