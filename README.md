# DNA Methylation Simulation

A biologically realistic simulation framework for modeling DNA methylation patterns and epigenetic drift over time. Simulates how cells accumulate methylation through growth, division, and homeostasis.

## Overview

The simulation models:
- Single-cell origin with exponential growth and homeostasis
- Stochastic methylation accumulation (epigenetic drift)
- Cell division with methylation pattern inheritance
- Population dynamics and steady-state maintenance
- Jensen-Shannon divergence (JSD) tracking from baseline patterns

## Quick Start

```bash
# 1. Run Phase 1 simulation
cd phase1
python run_simulation.py --config config_default.yaml

# 2. Generate Phase 2 datasets
cd ../phase2
python phase2_pipeline.py --simulation ../data/gene_rates_*/simulation.json.gz

# 3. Run Phase 3 analysis
cd ../phase3
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/
```

## Installation

```bash
git clone <repository-url>
cd jpi-methylation-simulation
pip install -r requirements.txt
```

**Dependencies:**
- `numpy` (>=1.19.0) - Numerical operations
- `scipy` (>=1.7.0) - Statistical analysis
- `plotly` (>=5.0.0, <6.0.0) - Visualization
- `kaleido` (0.2.1) - PNG export
- `pyyaml` (>=6.0) - Configuration files

## Architecture

### Three-Phase Pipeline

```
Phase 1: Simulation    →    Phase 2: Data Generation    →    Phase 3: Analysis
   (cell.py)                  (phase2_pipeline.py)           (run_pipeline.py)
```

**Phase 1: Biologically Realistic Simulation**
- Single cell origin with exponential growth and homeostasis
- Stochastic methylation accumulation and inheritance
- Object-oriented Cell and PetriDish classes
- Outputs to `data/gene_rates_*/size*-seed*-{timestamp}/simulation.json.gz`

**Phase 2: Data Generation Pipeline**
- Extracts snapshots from Phase 1 simulations
- Creates test populations via quantile-based and random sampling
- Grows individual cell populations with homeostasis
- Generates control populations from pure snapshots
- Outputs directly into Phase 1 simulation directory

**Phase 3: Analysis and Visualization**
- Extracts metrics to CSV files
- Generates histograms, timelines, and comparison plots
- Calculates statistical tests (t-tests, ANOVA)
- Re-analyzable without re-simulation

### Directory Structure

```
data/gene_rates_*/size*-seed*-{timestamp}/           # Phase 1 output
├── simulation.json.gz                               # Complete simulation
└── snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}-seed{S}-{timestamp}/
    ├── snapshots/                                   # Phase 2: Extracted snapshots
    │   ├── year{N}_snapshot.json
    │   └── metadata.json
    └── individuals/                                 # Phase 2: Generated populations
        ├── test2/                                   # Quantile-sampled
        ├── test1/                                   # Random-sampled
        ├── control/                                 # Pure snapshot
        └── mixing_metadata.json
```

## Usage

### Phase 1: Simulation

```bash
cd phase1

# Using config file (recommended)
python run_simulation.py --config config_default.yaml

# With CLI arguments
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Gene-specific methylation rates
python run_simulation.py --gene-rate-groups "5:0.004,5:0.005,5:0.006,5:0.007" --gene-size 5
```

**Key Parameters:**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--rate` | Uniform methylation rate per site per year | 0.005 |
| `--gene-rate-groups` | Gene-specific rates as `n1:rate1,n2:rate2,...` | - |
| `--years` | Total simulation time | 50 |
| `--growth-phase` | Years of exponential growth (final size = 2^N) | 9 |
| `--sites` | Number of CpG sites per cell | 100 |
| `--gene-size` | Sites per gene | 5 |
| `--seed` | Random seed for reproducibility | 42 |

**Growth Phases:**
- **Growth Phase** (Years 1 to N): Population doubles each year, deterministic 2^year cells
- **Steady State** (Years N+1 onwards): Divide → 50% cull → methylate, maintains ~2^N cells

### Phase 2: Data Generation

```bash
cd phase2

# Standard pipeline
python phase2_pipeline.py --simulation ../data/gene_rates_*/simulation.json.gz

# With custom parameters
python phase2_pipeline.py --simulation ../data/gene_rates_*/simulation.json.gz \
    --first-snapshot 30 --second-snapshot 50 \
    --n-quantiles 10 --cells-per-quantile 3 \
    --individual-growth-phase 7 --mix-ratio 80
```

**Key Parameters:**
| Parameter | Description | Default |
|-----------|-------------|---------|
| `--first-snapshot` | First snapshot year | 30 |
| `--second-snapshot` | Second snapshot year | 50 |
| `--n-quantiles` | Number of quantiles for stratified sampling | 10 |
| `--cells-per-quantile` | Cells sampled per quantile | 5 |
| `--individual-growth-phase` | Growth years per individual | 6 |
| `--mix-ratio` | % of snapshot cells in final mix | 80 |

**Pipeline Stages:**
1. Extract snapshots from two time points
2. Create Test 2 individuals (quantile-based sampling by JSD)
3. Create Test 1 individuals (uniform random sampling)
4. Grow populations with exponential growth then homeostasis
5. Normalize population sizes (median - 0.5σ threshold)
6. Mix grown cells with snapshot cells
7. Create Control populations (pure second snapshot)

### Phase 3: Analysis

```bash
cd phase3

# Complete analysis pipeline
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/

# Skip plots (CSV extraction only)
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/ --skip-plots

# Skip per-gene plots (faster)
python run_pipeline.py --phase2-dir ../data/{phase1_dir}/{phase2_subdir}/ --skip-gene-comparison
```

**Pipeline Stages:**
1. Extract timeline from simulation → 4 CSV files
2. Generate histogram plots (2 years × 2 metrics)
3. Generate timeline plots
4. Extract batch comparisons → 2 CSV files
5. Calculate statistical tests → p-values and ANOVA CSVs
6. Generate comparison plots (mean + std)
7. Generate per-gene plots (20 genes × 2 metrics)

## Key Concepts

### CpG Sites
Cytosine-guanine dinucleotides where methylation occurs. Initially unmethylated (0), they become methylated (1) with a probability each year.

### Genes
CpG sites are grouped into genes. Each gene can have a different methylation rate, modeling biological variation.

### Jensen-Shannon Divergence (JSD)
A symmetric measure of divergence between a cell's methylation distribution and the baseline (unmethylated). Quantifies how far a cell has drifted epigenetically.

### Epigenetic Drift
The gradual accumulation of methylation changes over time through stochastic methylation events.

### Quantile Sampling (Test 2)
Stratified sampling that selects cells from each JSD quantile, ensuring representation across the full distribution.

### Random Sampling (Test 1)
Uniform random sampling that may over/under-represent certain JSD ranges.

## Project Structure

```
jpi-methylation-simulation/
├── phase1/                        # Simulation Engine
│   ├── cell.py                    # Core Cell and PetriDish classes
│   ├── run_simulation.py          # Main simulation runner
│   └── config_default.yaml        # Default configuration
├── phase2/                        # Data Generation
│   ├── phase2_pipeline.py         # Single-file pipeline
│   └── config_default.yaml        # Default configuration
├── phase3/                        # Analysis
│   ├── run_pipeline.py            # Main analysis runner
│   ├── extract_*.py               # Data extraction scripts
│   ├── plot_*.py                  # Plotting scripts
│   └── calculate_pvalues.py       # Statistical tests
├── requirements.txt
├── README.md
└── CLAUDE.md                      # AI assistant guidance
```

## Testing

```bash
# Phase 1 tests
cd phase1/tests
python test_small.py               # Quick validation
python test_comprehensive.py       # Full feature tests

# Phase 2 tests
cd phase2/tests
python test_validation.py          # Data validation
python test_rate_consistency_phase2.py
```

## Reproducibility

Simulations are fully reproducible when using `--seed`:
- Both Python `random` and NumPy RNGs are seeded
- Phase 2 uses offset seeds for different operations
- Deterministic sampling and file I/O order

## License

[Add license information here]
