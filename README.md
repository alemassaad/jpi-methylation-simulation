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

## Three-Phase Architecture: Simulation → Data Generation → Analysis

The complete pipeline is now organized into three distinct phases:

### Phase 1: Biologically Realistic Simulation
- Single cell origin with exponential growth and homeostasis
- Stochastic methylation accumulation and inheritance
- Object-oriented Cell and PetriDish classes
- Lean JSON format with ~90% size reduction
- Hierarchical output directory structure

### Phase 2: Data Generation Pipeline
- Extract snapshots from phase1 simulations
- Create and grow individual cell populations
- Implement quantile-based sampling strategies
- Generate mutant and control populations
- Full reproducibility with comprehensive seeding

### Phase 3: Analysis and Visualization
- Pure analysis pipeline (no simulation)
- Generate all plots and statistical analysis
- Cell-level and gene-level metric visualization
- Population comparisons and timeline analysis
- Re-analyzable without re-simulation

This separation provides:
- **Clean Architecture**: Clear separation of concerns
- **Efficiency**: Re-analyze without re-simulation
- **Flexibility**: Multiple analysis approaches on same data
- **Scalability**: Batch processing and parallel analysis

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

# Using your own config (copy from config_default.yaml)
python run_simulation.py --config my_config.yaml

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

### Phase 2: Data Generation Pipeline
```bash
cd phase2

# Using config file (recommended)
python run_pipeline.py --config configs/quick_test.yaml \
    --simulation ../phase1/data/*/simulation.json.gz

# Standard data generation (30 individuals from 10 deciles)
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 --second-snapshot 50 \
    --individual-growth-phase 7 --seed 42

# Quick test (8 individuals from quartiles, short growth)
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 --second-snapshot 50 \
    --individual-growth-phase 4 \
    --n-quantiles 4 --cells-per-quantile 2

# Output structure:
# data/gene_rates_*/snap30to50-growth7-quant10x3-mix80-seed42-XXXX/
#   ├── snapshots/     # Extracted cell snapshots
#   └── individuals/   # Mutant, Control1, Control2 populations
#       ├── mutant/
#       ├── control1/
#       ├── control2/
#       └── mixing_metadata.json
```

### Phase 3: Analysis and Visualization
```bash
cd phase3

# Analyze phase2 data
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz

# Using config file
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --config configs/quick_analysis.yaml

# Custom analysis parameters
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --bins 150 --max-gene-plots 10

# Output structure:
# data/analysis_bins200_TIMESTAMP/results/
#   ├── cell_metrics/      # Cell-level analysis
#   │   ├── distributions/
#   │   ├── comparisons/
#   │   ├── individual_trajectories/
#   │   └── timeline/
#   ├── gene_metrics/      # Gene-level analysis
#   │   ├── distributions/
#   │   ├── comparisons/
#   │   ├── per_gene/
#   │   └── timeline/
#   └── metadata/
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

### Phase 2: Data Generation Pipeline

Generates structured datasets from phase1 simulations for downstream analysis:

```bash
cd phase2

# Complete data generation pipeline
python run_pipeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --first-snapshot 30 \           # First snapshot year
    --second-snapshot 50 \          # Second snapshot year  
    --individual-growth-phase 7 \   # Growth before homeostasis (2^7=128 cells)
    --n-quantiles 10 \              # Number of quantiles for sampling
    --cells-per-quantile 3 \        # Cells sampled per quantile
    --mix-ratio 80 \                # % of snapshot cells in final mix
    --seed 42                       # Random seed for reproducibility
```

**Data Generation Stages:**

1. **Extract Snapshots** (`extract_snapshots.py`)
   - Extract cells from two time points
   - Cache in `snapshots/year{N}_snapshot.json.gz`
   - Validate gene rate group consistency
   
2. **Create Initial Individuals** (`simulate_individuals.py`)
   - Mutant: Quantile-based sampling (e.g., 3 from each decile)
   - Control1: Uniform random sampling
   
3. **Grow Populations**
   - Exponential growth for individual-growth-phase years
   - Then homeostasis (divide → cull → methylate)
   - Maintains population around 2^growth_phase cells
   
4. **Mix with Snapshots**
   - Add second snapshot cells to grown individuals
   - Default: 80% snapshot, 20% grown
   - Optional uniform mixing and size normalization
   
5. **Create Control2** (`create_control2.py`)
   - Pure second snapshot populations
   - Same total count as mixed populations
   - Uses uniform pool if uniform mixing was applied

### Phase 3: Analysis and Visualization

Pure analysis pipeline that reads phase2 data and generates comprehensive visualizations:

```bash
cd phase3

# Complete analysis pipeline
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --bins 200 \                    # Histogram bins
    --max-gene-plots 20 \           # Limit per-gene plots
    --output-dir custom_analysis     # Custom output location
```

**Analysis Stages:**

1. **Snapshot Distributions**
   - Cell and gene JSD distributions
   - Methylation proportion histograms
   - Statistical overlays (mean, median, SD, CV, percentiles)
   
2. **Population Comparisons**
   - Compare mutant vs control1 vs control2
   - Statistical significance testing
   - Batch-level comparison plots
   
3. **Individual Trajectories**
   - Growth trajectories for each individual
   - Cell JSD and methylation progression over time
   - Gene-level evolution tracking
   
4. **Timeline Analysis**
   - Original simulation timeline plots
   - Full history visualization from phase1
   - Gene JSD heatmaps and progression

## Key Features

### Hierarchical Directory Structure
```
Phase 1 output:
data/gene_rates_200x0.00500/size8192-sites1000-genesize5-years100-seed42-YYYYMMDDHHMMSS/
├── simulation.json.gz           # Complete simulation history
├── jsd_history.png             # Cell JSD timeline
└── methylation_history.png     # Methylation timeline

Phase 2 output:
data/gene_rates_200x0.00500-grow9-sites1000-years100/
  └── snap30to50-growth7-quant10x3-mix80[u][n]-seed42-YYYYMMDDHHMMSS/
      ├── snapshots/            # Extracted cell snapshots
      │   ├── year30_snapshot.json.gz
      │   ├── year50_snapshot.json.gz
      │   └── metadata.json
      └── individuals/          # Generated populations
          ├── mutant/
          ├── control1/
          ├── control2/
          └── mixing_metadata.json

Phase 3 output:
data/analysis_bins200_YYYYMMDDHHMMSS/
  └── results/
      ├── cell_metrics/         # Cell-level analysis
      │   ├── distributions/
      │   ├── comparisons/
      │   ├── individual_trajectories/
      │   └── timeline/
      ├── gene_metrics/         # Gene-level analysis
      │   ├── distributions/
      │   ├── comparisons/
      │   ├── per_gene/
      │   └── timeline/
      └── metadata/
      
Phase 2 suffix indicators:
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

## Advanced Features

### Phase 2: Data Generation Options

#### Uniform Mixing (`--uniform-mixing`)
- All individuals receive the exact same set of snapshot cells
- Reduces inter-individual variation from random sampling
- Useful for focusing on biological variation rather than sampling noise
- Directory suffix: 'u' (e.g., mix80u)

#### Size Normalization (Always Applied)
- Normalizes all individuals to the same cell count before mixing
- Uses median - 0.5σ threshold to determine target size
- Excludes individuals below threshold, trims those above
- Addresses variation from homeostasis stochasticity
- Typically retains ~67% of individuals
- Ensures fair comparison between individuals

#### Homeostasis in Individual Growth
- `--individual-growth-phase`: Years of exponential growth before homeostasis
- After growth phase: divide → cull → methylate each year
- Maintains population around 2^growth_phase cells
- More biologically realistic than pure exponential growth
- Example: `--individual-growth-phase 7` → ~128 cells

### Phase 3: Analysis Options

#### Flexible Visualization
- Configurable histogram bins (`--bins`)
- Limit gene plots for faster analysis (`--max-gene-plots`)
- Custom output directories (`--output-dir`)
- YAML configuration files for consistent analysis

#### Re-analysis Capability
- Analyze same phase2 data with different parameters
- No need to re-run expensive simulations
- Batch processing of multiple phase2 runs
- Independent analysis timeline from data generation

## Project Structure

```
jpi-methylation-simulation/
├── phase1/                    # Simulation Engine
│   ├── cell.py                    # Core classes (Cell, PetriDish)
│   ├── run_simulation.py          # Main simulation runner
│   ├── plot_history.py            # Time-series visualization
│   ├── configs/                   # YAML configuration files
│   ├── tests/                     # Unit tests
│   └── data/                      # Simulation output
├── phase2/                    # Data Generation Pipeline
│   ├── run_pipeline.py            # Pipeline orchestrator
│   ├── extract_snapshots.py       # Extract cell snapshots
│   ├── simulate_individuals.py    # Create and grow populations
│   ├── create_control2.py         # Generate control populations
│   ├── core/                      # Pipeline utilities
│   │   ├── pipeline_utils.py      # Cell/PetriDish utilities
│   │   ├── path_utils.py          # Path parsing/generation
│   │   └── validation.py          # Data validation
│   ├── configs/                   # Configuration files
│   ├── tests/                     # Pipeline tests
│   └── data/                      # Generated datasets
├── phase3/                    # Analysis and Visualization
│   ├── run_analysis.py            # Main analysis runner
│   ├── core/                      # Analysis functions
│   │   ├── data_loader.py         # Data loading utilities
│   │   ├── analysis_functions.py  # Statistical analysis
│   │   └── plot_paths.py          # Plot organization
│   ├── configs/                   # Analysis configurations
│   └── data/                      # Analysis results
├── requirements.txt               # Python dependencies
├── README.md                      # This file
├── CLAUDE.md                      # AI assistant guidance
└── PLOT_DOCUMENTATION.md          # Plot descriptions
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

### Phase 2 Tests (Data Generation)
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

### Phase 3 Tests (Analysis)
```bash
cd phase3/tests  # (if test directory exists)
# Test analysis functions, plot generation, data loading
# Configuration validation, output format verification
```

## License

[Add license information here]