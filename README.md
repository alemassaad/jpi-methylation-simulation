# DNA Methylation Simulation

A Python-based simulation framework for modeling DNA methylation patterns in cells over time. This project simulates how CpG sites become methylated as cells age, tracking the distribution of methylation across genes and calculating divergence from baseline patterns.

## Overview

The simulation models:
- Individual cells with configurable numbers of CpG sites
- Stochastic methylation events over time
- Methylation distribution across genes
- Jensen-Shannon divergence (JSD) from baseline methylation patterns
- Population-level statistics across thousands of cells

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

## Usage

### Quick Start

```bash
# Run standard simulation with multiple rates
python main.py

# Run large-scale optimized simulation
python run_large_fast_save.py

# Plot results
python plot_history.py baseline_rate_0.5%_m10000_n1000_t100.json.gz
```

### Running Simulations

**Standard multi-rate simulation:**
```bash
python main.py
```
- Runs with rates: 0.002, 0.005, 0.01
- Saves uncompressed JSON files to `history/`

**Large-scale optimized simulation:**
```bash
python run_large_fast_save.py
```
- Parameters: rate=0.5%, n=1000, m=10,000, t_max=100
- Progress tracking with time estimates
- Compressed output (.json.gz)
- Runtime: ~10 minutes

**Quick test:**
```bash
python test_small.py
```
- Small parameters for testing: rate=1%, n=15, m=5, t_max=5

### Visualizing Results

Create plots from simulation output (automatically saved to `plots/` directory):

```bash
# Create both JSD and methylation plots (separate PNGs in plots/)
python plot_history.py history/simulation_file.json.gz

# Create only JSD plot
python plot_history.py history/simulation_file.json.gz --jsd-only

# Create only methylation plot
python plot_history.py history/simulation_file.json.gz --methylation-only

# Specify custom output name
python plot_history.py history/simulation_file.json.gz -o custom_name
```

Plots are saved to `plots/` directory and show:
- Mean line (solid)
- 25-75 percentile band (darker shading)
- 5-95 percentile band (lighter shading)
- Final statistics annotation

### Step 2: Cell Division Experiments

The `step2/` directory contains scripts for analyzing cell division effects:

**Extract a snapshot from simulation:**
```bash
cd step2
python extract_snapshot.py ../history/simulation_rate_0.005_m10000_n1000_t100.json.gz 50 year50_snapshot.json.gz
```

**Plot JSD distribution at a specific year:**
```bash
python plot_jsd_distribution.py year50_snapshot.json.gz 200  # 200 bins
```

**Run cell division simulation with lineage tracking:**
```bash
python sample_divide_age_lineages.py
```
This improved script:
- Samples 3 cells from each JSD decile (30 cells total)
- Ages each cell SEPARATELY for 10 years with division
- Creates 30 individual lineage files in `lineages/` directory
- Each lineage contains 1,024 cells from one original cell
- Enables tracking of individual cell lineages

### Step 3: Mixed Population Analysis

The `step3/` directory analyzes how lineage cells affect populations when mixed:

**Create mixed individuals:**
```bash
cd step3
python extract_year60_original.py  # Extract year 60 baseline
python create_individuals.py       # Create 30 mixed individuals (80% original, 20% lineage)
```

**Create control individuals:**
```bash
python create_control_individuals.py  # Create 30 control individuals (100% original)
```

**Compare distributions:**
```bash
python plot_distributions.py
```
Generates a visualization comparing mean JSD between mixed and control populations.

### Customizing Parameters

Edit constants in `cell.py`:
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
├── cell.py                    # Core simulation classes (Cell, History) and math functions
├── main.py                    # Standard simulation runner
├── run_large_fast_save.py     # Optimized large-scale simulation
├── test_small.py              # Quick test script
├── plot_history.py            # Visualization script
├── requirements.txt           # Python dependencies
├── history/                   # Simulation output directory
│   └── *.json, *.json.gz     # Simulation data files
├── plots/                     # Visualization output directory
│   └── *_jsd.png, *_methylation.png  # Generated plots
├── tests/                     # Reproducibility tests
│   ├── test_reproducibility.py
│   ├── test_reproducibility_expected.json
│   └── test_reproducibility_README.md
├── step2/                     # Cell division experiments
│   ├── extract_snapshot.py    # Extract specific year from simulation
│   ├── sample_divide_age_lineages.py # Age cells separately with lineage tracking
│   ├── plot_jsd_distribution.py # Plot JSD distributions
│   ├── lineages/              # 30 separate lineage files
│   └── plots/                 # Visualization outputs
├── step3/                     # Mixed population analysis
│   ├── extract_year60_original.py # Extract year 60 baseline
│   ├── create_individuals.py  # Create mixed individuals
│   ├── create_control_individuals.py # Create control individuals
│   ├── plot_distributions.py  # Compare distributions
│   ├── individuals/           # 30 mixed individuals
│   ├── control_individuals/   # 30 control individuals
│   └── plots/                 # Comparison visualizations
├── README.md                  # This file
└── CLAUDE.md                  # Guidance for AI assistants
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