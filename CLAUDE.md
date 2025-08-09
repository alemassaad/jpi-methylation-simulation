# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a biologically realistic methylation simulation project that models DNA methylation patterns over time using object-oriented design. The simulation tracks epigenetic drift in cell populations through growth, division, and homeostasis.

## Main Pipeline: Step1-Prime → Step23-Prime

**IMPORTANT: The prime versions (step1-prime and step23-prime) are the main production pipelines. Use these for all new work.**

### Primary Workflow
1. **Step1-Prime**: Biologically realistic simulation starting from single cell
2. **Step23-Prime**: Object-oriented analysis pipeline using Cell/PetriDish classes

## Architecture & Code Structure

### Core Simulation (step1-prime/) - MAIN
- **`cell.py`**: Core OOP classes
  - `Cell`: Individual cell with methylation tracking and mitosis
  - `PetriDish`: Population manager with growth/homeostasis logic
  - Mathematical functions: `KL_div`, `JS_div` for divergence calculations
- **`run_simulation.py`**: Main CLI runner with progress tracking
- **Hierarchical output**: `data/rate_X.XXXXX/growG-sitesN-yearsT-seedS-HASH/simulation.json.gz`

### Analysis Pipeline (step23-prime/) - MAIN
- **`run_pipeline.py`**: 8-stage orchestrator with dynamic snapshot years
- **`pipeline_utils.py`**: Cell/PetriDish manipulation utilities
- **`pipeline_analysis.py`**: Statistical analysis and visualization
- **`path_utils.py`**: Hierarchical path parsing/generation
- **Hierarchical output**: `data/rate_X-growG-sitesN-yearsT/snapS-quantQxC-growG-mixM-seedS-HASH/`

## Key Commands

### Step1-Prime: Run Biologically Realistic Simulation
```bash
cd step1-prime
# Standard run (1 cell → 8192 cells, 100 years)
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --seed 42

# Quick test (1 cell → 16 cells, 20 years)
python run_simulation.py --rate 0.01 --years 20 --growth-phase 4 --seed 42
```

### Step23-Prime: Run Analysis Pipeline
```bash
cd step23-prime
# Standard analysis (30 individuals from 10 deciles)
python run_pipeline.py --rate 0.005 \
    --simulation ../step1-prime/data/rate_0.00500/grow13-*/simulation.json.gz \
    --snapshot-year 50 --growth-years 10 --seed 42

# Quick test (4 individuals from quartiles, 2 years growth)
python run_pipeline.py --rate 0.005 \
    --simulation ../step1-prime/data/rate_0.00500/grow13-*/simulation.json.gz \
    --n-quantiles 4 --cells-per-quantile 1 --growth-years 2
```

**Pipeline Parameters:**
- `--snapshot-year 50`: First extraction year (dynamic, not hardcoded)
- `--growth-years 10`: Years to grow individuals (second snapshot = first + growth)
- `--n-quantiles 10`: Number of quantiles for stratified sampling
- `--cells-per-quantile 3`: Cells sampled per quantile
- `--mix-ratio 80`: Percentage of snapshot cells in final mix
- `--seed 42`: Random seed for full reproducibility

## Key Implementation Details

### Biological Realism (step1-prime)
```python
# Growth phase: Years 0-13 (configurable)
if self.year <= self.growth_phase:
    self.divide_cells()     # Population doubles deterministically
    self.methylate_cells()  # Apply stochastic methylation
    # Result: exactly 2^year cells

# Steady state: Years 14+ 
else:
    self.divide_cells()      # Double population
    self.random_cull_cells() # ~50% survival (stochastic)
    self.methylate_cells()   # Apply methylation
    # Result: ~8192 cells (varies randomly)
```

### Object-Oriented Design (step23-prime)
```python
# Direct Cell/PetriDish manipulation - no dict conversions!
cells = load_snapshot_as_cells(simulation_file, year=50)  # Returns List[Cell]
petri = PetriDish(rate=0.005, n=1000)
petri.cells = [cell]  # Start with single cell
grow_petri_for_years(petri, years=10)  # Uses PetriDish methods
```

### Hierarchical Directory Structure
```python
# step1-prime: data/rate_0.00500/grow13-sites1000-years100-seed42-a3f2/
level1 = f"rate_{rate:.5f}"
level2 = f"grow{growth}-sites{n}-years{t}-seed{seed}-{hash[:4]}"

# step23-prime: data/rate_0.00500-grow13-sites1000-years100/snap50-quant10x3-grow10-mix80-seed42-b5e1/
level1 = f"rate_{rate:.5f}-grow{growth}-sites{n}-years{t}"
level2 = f"snap{snap}-quant{q}x{c}-grow{g}-mix{m}-seed{s}-{hash[:4]}"
```

### Full Reproducibility
```python
# Set at pipeline start for deterministic results
import random
import numpy as np
random.seed(args.seed)
np.random.seed(args.seed)

# Also seed in every sampling function
def sample_by_quantiles(cells, n_quantiles, cells_per_quantile, seed):
    random.seed(seed)
    np.random.seed(seed)  # Both for full reproducibility
```

## Pipeline Process (8 Stages)

1. **Extract first snapshot** → Cell objects (cached)
2. **Plot JSD distribution** → Histogram with stats overlay
3. **Create individuals**: 
   - Mutant: Quantile-based stratified sampling
   - Control1: Uniform random sampling
4. **Grow N years**: PetriDish.divide_cells() + methylate_cells()
5. **Extract second snapshot** → Age-aligned cells
6. **Mix populations**: Add snapshot cells maintaining ratio
7. **Create control2**: Pure second snapshot individuals
8. **Analysis**: Scatter plots, t-tests, comprehensive statistics

## Important Constants (step1-prime/cell.py)
- `N = 1000`: CpG sites per cell
- `RATE = 0.005`: Default methylation rate (0.5%)
- `GENE_SIZE = 5`: Sites per gene
- `T_MAX = 100`: Default simulation years
- `DEFAULT_GROWTH_PHASE = 13`: Default growth years (→ 8192 cells)
- `BASELINE_METHYLATION_DISTRIBUTION`: Reference for JSD calculation

## Critical Implementation Notes

1. **Always use prime versions** for new work
2. **Cell.methylate()** is the core aging mechanism (formerly age_1_year)
3. **PetriDish** manages populations with biologically realistic dynamics
4. **Dynamic snapshot years** ensure age alignment (not hardcoded 50/60)
5. **Use deep copy** for cells to avoid reference bugs
6. **Type consistency**: Use float literals (1.0, 0.0) in distributions
7. **Random seeds**: Set globally AND in each function for reproducibility
8. **MD5 hashing**: 4-char hashes ensure unique directory names

## Testing & Validation

```bash
# Test reproducibility (step23-prime)
cd step23-prime/tests
python test_reproducibility_robust.py

# Compare prime vs legacy pipelines
cd step23-prime/tools
python generate_comparison_report.py
```

## Common Development Tasks

### Debugging
```bash
# Inspect compressed files
zcat file.json.gz | python -m json.tool | head -20

# Check cell counts
python -c "import json, gzip; data=json.load(gzip.open('file.json.gz', 'rt')); print(len(data['cells']))"

# Validate hierarchical paths
ls -la step1-prime/data/rate_*/grow*/
ls -la step23-prime/data/rate_*/snap*/
```

### Performance
- Step1-prime simulation: ~5-10 minutes (8192 cells, 100 years)
- Step23-prime pipeline: ~15-20 minutes (complete analysis)
- Use `--n-quantiles 4 --cells-per-quantile 1 --growth-years 2` for quick tests

---

## Legacy Pipelines (Historical Reference Only)

### Legacy Directories
- `step1/`: Traditional simulation (10,000 cells, no growth)
- `step2/`, `step3/`: Original multi-step approach
- `step23/`: Pre-OOP unified pipeline
- `legacy/`: Archived old implementations

These are maintained only for:
- Backward compatibility
- Validation of prime pipeline results
- Historical documentation

**DO NOT use legacy pipelines for new work. Always use step1-prime and step23-prime.**