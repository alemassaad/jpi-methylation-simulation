# DNA Methylation Simulation

A Python-based simulation framework for modeling DNA methylation patterns in cells over time. This project simulates how CpG sites become methylated as cells age, tracking the distribution of methylation across genes and calculating divergence from baseline patterns.

## Overview

The simulation models:
- Individual cells with configurable numbers of CpG sites
- Stochastic methylation events over time
- Methylation distribution across genes
- Jensen-Shannon divergence (JSD) from baseline methylation patterns
- Population-level statistics across thousands of cells
- Cell division with lineage tracking
- Mixed population experiments

## Pipeline Overview

The project is organized into three sequential steps:

1. **Step 1**: Run base methylation simulation over 100 years
2. **Step 2**: Extract cells at year 50, perform cell division experiments with lineage tracking
3. **Step 3**: Mix aged lineages with original year 60 cells to analyze population dynamics

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

## Quick Start

```bash
# Step 1: Run simulation
cd step1
python run_large_simulation.py --rate 0.005

# Step 2: Create lineages from year 50
cd ../step2/scripts
python extract_snapshot.py
python create_lineages.py --type both

# Step 3: Create mixed populations and analyze
cd ../../step3/scripts
python extract_year60_original.py
python create_individuals.py --type both
python plot_distributions.py
```

## Detailed Usage

### Step 1: Base Simulation

Run methylation simulation over time:

```bash
cd step1

# Run with default parameters (rate=0.005, n=1000, m=10000, t=100)
python run_large_simulation.py

# Run with custom methylation rate
python run_large_simulation.py --rate 0.01

# Quick test with small parameters
python test_small.py

# Visualize results
python plot_history.py data/simulation_rate_0.005000_m10000_n1000_t100.json.gz
```

### Step 2: Cell Division Experiments

Extract cells and perform division experiments:

```bash
cd step2/scripts

# Extract year 50 snapshot from simulation
python extract_snapshot.py --year 50

# Create both mutant and control lineages
python create_lineages.py --type both

# Plot JSD distribution of a snapshot
python plot_jsd_distribution.py ../data/snapshots/year50_snapshot.json.gz 200
```

Lineage creation:
- **Mutant lineages**: Samples 3 cells from each JSD decile (30 total)
- **Control lineages**: Samples 30 cells uniformly from population
- Each lineage undergoes 10 years of division, resulting in 1,024 cells

### Step 3: Mixed Population Analysis

Analyze mixed populations:

```bash
cd step3/scripts

# Extract year 60 from original simulation
python extract_year60_original.py

# Create mixed individuals (both mutant and control)
python create_individuals.py --type both

# Plot comparison
python plot_distributions.py
```

Individual creation:
- Each individual: 80% original year 60 cells + 20% lineage cells
- Creates 30 individuals for both mutant and control groups

## Customizing Parameters

Edit constants in `step1/cell.py`:
- `N`: Number of CpG sites per cell (default: 1000)
- `M`: Number of cells in population (default: 10,000)  
- `T_MAX`: Maximum simulation time in years (default: 100)
- `GENE_SIZE`: Number of CpG sites per gene (default: 5)
- `RATE`: Default methylation rate (default: 0.01)

## Output Format

Simulations produce compressed JSON files (`.json.gz`) with structure:

Example structure:
```json
{
  "0": [
    {
      "cpg_sites": [0, 0, 0, ...],
      "methylation_proportion": 0.0,
      "methylation_distribution": [1.0, 0.0, 0.0, 0.0, 0.0, 0.0],
      "jsd": 0.0,
      "age": 0,
      "gene_size": 5,
      "rate": 0.01
    },
    ...
  ],
  ...
}
```

## Project Structure

```
jpi-methylation-simulation/
├── step1/                     # Base methylation simulation
│   ├── cell.py               # Core simulation classes (Cell, History)
│   ├── run_large_simulation.py # Main simulation runner with CLI
│   ├── main.py               # Legacy multi-rate runner
│   ├── test_small.py         # Quick test script
│   ├── plot_history.py       # Visualization script
│   └── data/                 # Simulation outputs
│       └── *.json.gz         # Compressed simulation data
├── step2/                     # Cell division experiments
│   ├── scripts/
│   │   ├── extract_snapshot.py    # Extract year from simulation
│   │   ├── create_lineages.py     # Create mutant & control lineages
│   │   ├── plot_jsd_distribution.py # Plot JSD distributions
│   │   └── test_pipeline.py       # Test the pipeline
│   └── data/
│       ├── snapshots/        # Year snapshots
│       ├── lineages/
│       │   ├── mutant/       # 30 decile-based lineages
│       │   └── control/      # 30 uniform lineages
│       └── plots/            # Visualization outputs
├── step3/                     # Mixed population analysis
│   ├── scripts/
│   │   ├── extract_year60_original.py # Extract year 60
│   │   ├── create_individuals.py      # Create mixed populations
│   │   ├── plot_distributions.py      # Compare distributions
│   │   └── test_pipeline.py          # Test the pipeline
│   └── data/
│       ├── snapshots/        # Year 60 snapshots
│       ├── individuals/
│       │   ├── mutant/       # 30 mixed mutant individuals
│       │   └── control/      # 30 mixed control individuals
│       ├── plots/            # Comparison visualizations
│       └── results/          # Analysis results
├── tests/                     # Reproducibility tests
│   ├── test_reproducibility.py
│   └── test_reproducibility_expected.json
├── requirements.txt           # Python dependencies
├── README.md                  # This file
├── CLAUDE.md                  # Guidance for AI assistants
└── .gitignore                # Git ignore patterns
```

## Key Concepts

### CpG Sites
Cytosine-guanine dinucleotides where methylation occurs in DNA. Initially unmethylated (0), they can become methylated (1) with a certain probability each year.

### Gene-Level Distribution
CpG sites are grouped into genes. The distribution tracks how many genes have 0, 1, 2, ... up to GENE_SIZE methylated sites.

### Jensen-Shannon Divergence (JSD)
A symmetric measure of the difference between the cell's methylation distribution and a baseline distribution. Used to quantify how far a cell has diverged from the expected pattern.

## Mathematical Background

The simulation uses:
- **Kullback-Leibler (KL) divergence**: Measures difference between probability distributions
- **Jensen-Shannon (JS) divergence**: Symmetric version of KL divergence, calculated as the average of KL divergences to the midpoint distribution

## Contributing

Feel free to submit issues or pull requests. Please ensure any changes maintain compatibility with the existing data format and include appropriate documentation.

## License

[Add license information here]