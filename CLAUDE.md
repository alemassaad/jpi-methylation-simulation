# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

Biologically realistic DNA methylation simulation modeling epigenetic drift through cell growth, division, and homeostasis. Uses object-oriented design with Cell and PetriDish classes.

## Architecture

### Main Pipeline (USE THESE)
- **phase1/**: Simulation engine (single cell → population growth → homeostasis)
- **phase2/**: Analysis pipeline (quantile sampling → growth → mixing → statistics)

### Core Classes (phase1/cell.py)
```python
Cell(n=1000, rate=0.005)        # Individual cell with CpG sites
PetriDish(rate=0.005, growth_phase=13)  # Population manager
```

Key methods:
- `Cell.methylate()`: Apply stochastic methylation
- `Cell.create_daughter_cell()`: Mitosis (identical copy)
- `PetriDish.divide_cells()`: Population doubling
- `PetriDish.random_cull_cells()`: Homeostatic ~50% survival

### Directory Structure
```
# Hierarchical output with MD5 hash suffixes
phase1/data/rate_0.00500/grow13-sites1000-years100-seed42-XXXX/
phase2/data/rate_0.00500-grow13-sites1000-years100/snap50-quant10x3-grow10-mix80-seed42-XXXX/
```

## Commands

### Run Simulation
```bash
cd phase1
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42
# Quick test: --years 20 --growth-phase 4
```

### Run Analysis Pipeline
```bash
cd phase2
python run_pipeline.py --rate 0.005 \
    --simulation ../phase1/data/rate_0.00500/grow13-*/simulation.json.gz \
    --snapshot-year 50 --growth-years 10 --seed 42
# Quick test: --n-quantiles 4 --cells-per-quantile 1 --growth-years 2
```

### Run Tests
```bash
# Phase 1
cd phase1/tests
python test_small.py           # Quick validation
python test_comprehensive.py   # Full features
python test_edge_cases.py      # Edge cases

# Phase 2
cd phase2/tests
python test_reproducibility_robust.py   # Reproducibility
python test_dynamic_mix_year.py        # Dynamic year calculations
```

## Key Implementation Details

### Biological Model
- **Growth phase** (years 0-growth_phase): Deterministic 2^year cells
- **Steady state** (years growth_phase+1-T): Stochastic ~2^growth_phase cells
- **Methylation**: Random with rate per site per year
- **Inheritance**: Daughter cells inherit parent's methylation pattern

### Pipeline Stages (8 total)
1. Extract first snapshot (e.g., year 50)
2. Plot JSD distribution with statistics
3. Create individuals (mutant: quantile-based, control1: uniform)
4. Grow N years (PetriDish simulation)
5. Extract second snapshot (first + growth_years)
6. Mix populations (default 80% snapshot)
7. Create control2 (pure second snapshot)
8. Analysis (t-tests, scatter plots)

### Reproducibility
- Set seeds globally AND in each function
- Use deep copy for Cell objects
- MD5 hashing for unique directories
- Hierarchical caching of intermediate results

## Important Constants
```python
N = 1000                    # CpG sites per cell
RATE = 0.005               # Methylation rate (0.5%)
GENE_SIZE = 5              # Sites per gene
DEFAULT_GROWTH_PHASE = 13  # → 8192 cells
```

## Development Guidelines

### Always Use Current Versions
- phase1 and phase2 are production pipelines
- legacy/ directory for historical reference only

### Object-Oriented Design
```python
# Direct manipulation - no dict conversions
cells = load_snapshot_as_cells(simulation_file, year=50)
petri = PetriDish(rate=0.005)
petri.cells = [cell]
grow_petri_for_years(petri, years=10)
```

### Type Consistency
- Use float literals in distributions: `[1.0, 0.0, 0.0]`
- Deep copy cells to avoid reference bugs
- Dynamic years, not hardcoded 50/60

### Testing & Validation
```bash
# Quick validation of changes
cd phase1/tests && python test_small.py
cd phase2/tests && python test_reproducibility_robust.py

# Check reproducibility between runs
cd phase2/tools
python compare_two_runs.py --dir1 path1 --dir2 path2
```

## Dependencies
```bash
pip install -r requirements.txt
# Installs: plotly>=5.0.0,<6.0.0, kaleido==0.2.1
# Also needs: scipy, numpy (usually pre-installed)
```