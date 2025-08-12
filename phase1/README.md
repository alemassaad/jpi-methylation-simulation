# Phase 1: Biologically Realistic Methylation Simulation

A production-ready methylation simulation that models realistic cellular population dynamics starting from a single cell.

## Key Features

- **Single-cell origin**: Simulation starts with one unmethylated cell
- **Exponential growth phase**: Population doubles each year until reaching target size
- **Homeostatic maintenance**: After growth, population maintained through division and culling
- **Configurable population**: Target population = 2^growth_phase (default: 2^13 = 8192 cells)
- **Object-oriented design**: Clean Cell and PetriDish classes
- **Hierarchical output**: Organized directory structure with parameter tracking

## Quick Start

```bash
# Standard simulation (1 cell → 8192 cells over 100 years)
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Quick test (1 cell → 16 cells over 20 years)
python run_simulation.py --rate 0.01 --years 20 --growth-phase 4 --seed 42

# Visualize results
python plot_history.py data/rate_0.00500/grow13-*/simulation.json.gz
```

## Parameters

- `--rate`: Methylation rate per site per year (default: 0.005 = 0.5%)
- `--sites`: Number of CpG sites per cell (default: 1000)
- `--years`: Total simulation time (default: 100)
- `--gene-size`: Sites per gene (default: 5, must divide sites evenly)
- `--growth-phase`: Years of exponential growth (default: 13 → 8192 cells)
- `--seed`: Random seed for reproducibility (default: 42, use -1 for no seed)
- `--output`: Custom output filename (optional)
- `--no-save`: Run without saving results

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

## Output Structure

```
data/
└── rate_0.00500/
    └── grow13-sites1000-years100-seed42-XXXX/
        ├── simulation.json.gz    # Complete simulation history
        └── plots/                 # Visualization outputs
            ├── *_jsd.png         # JSD over time
            ├── *_methylation.png # Methylation over time
            └── *_combined.png    # Both plots
```

## Visualization

The `plot_history.py` script creates time-series plots showing:
- Jensen-Shannon divergence over time
- Methylation percentage over time
- Population size tracking
- Growth phase transition markers
- Mean with 5-95% and 25-75% percentile bands

```bash
# Create JSD plot only
python plot_history.py data/*/grow13-*/simulation.json.gz --jsd-only

# Create combined plot
python plot_history.py data/*/grow13-*/simulation.json.gz --combined

# Specify output directory
python plot_history.py data/*/grow13-*/simulation.json.gz -o my_plots/
```

## Core Classes

### Cell
Represents an individual cell with methylation sites:
- `methylate()`: Apply stochastic methylation
- `create_daughter_cell()`: Create identical copy (mitosis)
- `to_dict()`: Serialize for JSON output

### PetriDish
Manages cell population with biological dynamics:
- `divide_cells()`: All cells undergo mitosis
- `methylate_cells()`: Apply methylation to all cells
- `random_cull_cells()`: Stochastic removal (~50% survival)
- `simulate_year()`: Run one year of simulation
- `run_simulation()`: Complete multi-year simulation
- `save_history()`: Save to hierarchical directory structure

## Examples

### Different Population Sizes
```bash
# Small population (256 cells)
python run_simulation.py --growth-phase 8

# Medium population (1024 cells)
python run_simulation.py --growth-phase 10

# Large population (32768 cells)
python run_simulation.py --growth-phase 15
```

### Different Methylation Rates
```bash
# Low rate (0.1% per year)
python run_simulation.py --rate 0.001

# Standard rate (0.5% per year)
python run_simulation.py --rate 0.005

# High rate (1% per year)
python run_simulation.py --rate 0.01
```

## Integration with Phase 2

Output from phase1 serves as input for the phase2 analysis pipeline:

```bash
# Run simulation
python run_simulation.py --rate 0.005 --growth-phase 13 --seed 42

# Analyze results with phase2
cd ../phase2
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz
```

## Scientific Background

This simulation models epigenetic drift through:
- **Stochastic methylation**: Random CpG site methylation over time
- **Mitotic inheritance**: Methylation patterns preserved through division
- **Population dynamics**: Realistic growth and homeostasis
- **Jensen-Shannon divergence**: Quantifies deviation from baseline