# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Commands - Quick Start

### Complete Pipeline (Most Common Workflow)
```bash
# 1. Run Phase 1 simulation
cd phase1
python run_simulation.py --config config_default.yaml

# 2. Generate Phase 2 datasets (outputs to Phase 1 directory)
cd ../phase2
python run_pipeline.py --simulation ../phase1/data/gene_rates_*/simulation.json.gz

# 3. Extract and plot with Phase 3 analysis
cd ../phase3
python run_pipeline.py --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/
```

## Commands - Detailed

### Phase 1: Run Simulation
```bash
cd phase1

# Quick test (100 sites, 50 years, 512 cells)
python run_simulation.py --config config_default.yaml

# Full simulation with all options
python run_simulation.py --rate 0.005 --years 100 --growth-phase 13 --sites 1000 --seed 42

# Gene-specific methylation rates (default configuration)
python run_simulation.py --gene-rate-groups "5:0.004,5:0.005,5:0.006,5:0.007" --gene-size 5

# Output saved to data/gene_rates_*/size*-seed*-{timestamp}/simulation.json[.gz]
```

### Phase 2: Data Generation Pipeline
```bash
cd phase2

# Complete pipeline - outputs directly to Phase 1 directory
python run_pipeline.py --simulation ../phase1/data/gene_rates_*/size*-seed*-*/simulation.json

# With custom parameters (override config defaults)
python run_pipeline.py --simulation ../phase1/data/gene_rates_*/size*-seed*-*/simulation.json \
    --n-quantiles 10 --cells-per-quantile 3 --mix-ratio 80

# Run individual stages (for debugging/custom workflows)
python extract_snapshots.py --simulation ../phase1/data/*/simulation.json.gz --output-dir data/my_run
python simulate_individuals.py --base-dir data/my_run --n-quantiles 10 --cells-per-quantile 3
python create_control2.py --base-dir data/my_run
```

### Phase 3: Analysis Pipeline

Phase 3 is the streamlined CSV-first analysis pipeline that prioritizes modularity and reusability:

```bash
cd phase3

# Run complete analysis pipeline (6 stages)
python run_pipeline.py --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/

# Pipeline stages:
# Stage 1: Extract timeline from Phase 1 simulation → 4 CSV files
# Stage 2: Generate 4 histogram plots (2 years × 2 metrics)
# Stage 3: Generate 4 timeline plots from extracted CSVs
# Stage 4: Extract batch comparisons → batch_comparison_{cell,gene}.csv
# Stage 5: Generate 4 comparison plots
# Stage 6: Generate 40 per-gene plots (optional)

# Useful options
python run_pipeline.py --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/ --skip-plots  # Extract CSVs only

# Extract full timeline from Phase 1 simulation (standalone)
python extract_simulation_timeline.py \
    --simulation ../phase1/data/*/simulation.json.gz \
    --output-dir tables/

# Extract comparison data separately (split into cell and gene files)
python extract_batch_comparison.py \
    --individuals-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/individuals

# Generate single comparison plot from CSV
python plot_comparison_generic.py \
    --csv tables/batch_comparison_cell.csv \
    --column cell_jsd_mean \
    --title "Cell JSD Comparison" \
    --ylabel "Mean Cell JSD" \
    --output plots/cell_jsd_comparison.png

# Generate per-gene comparison plots (40 plots: 20 genes × 2 metrics)
python plot_comparison_by_gene.py \
    --csv tables/batch_comparison_gene.csv \
    --output-dir plots/per_gene/

# Plot simulation timeline from extracted CSV files
python plot_simulation_timeline.py results/tables/

# Custom results directory with plots
python run_pipeline.py --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/ \
    --results-dir /path/to/results
```

## High-Level Architecture

### Three-Phase Pipeline
1. **Phase 1**: Core simulation engine - generates cell populations over time
2. **Phase 2**: Data generation pipeline - creates structured datasets from phase1
3. **Phase 3**: Analysis pipeline - extracts data tables and generates plots

### Phase 3 Architecture

Phase 3 is a streamlined CSV-first analysis pipeline that prioritizes modularity and reusability:

```
phase3/
├── run_pipeline.py                # Main orchestrator for all stages
├── extract_simulation_timeline.py # Extract timeline from Phase 1 simulation
├── data_extractors.py             # Core data extraction utilities
├── plot_histogram_original.py     # Generate histogram plots
├── plot_simulation_timeline.py    # Plot timeline from CSVs
├── extract_batch_comparison.py    # Extract batch comparison data
├── plot_comparison_generic.py     # Generic comparison plotting
├── plot_comparison_by_gene.py     # Per-gene comparison plots
└── test_gene_extraction.py        # Test file
```

**Pipeline Execution Flow**:
```
Phase 2 data directory
    ↓
Stage 1: Extract timeline from ../simulation.json.gz
    → simulation_cell_jsd.csv
    → simulation_cell_methylation.csv
    → simulation_gene_jsd.csv
    → simulation_gene_methylation.csv
    ↓
Stage 2: Generate histogram plots (4 plots)
    → year30_cell_jsd_histogram.png
    → year30_cell_methylation_histogram.png
    → year50_cell_jsd_histogram.png
    → year50_cell_methylation_histogram.png
    ↓
Stage 3: Generate timeline plots (4 plots)
    → simulation_cell_jsd_timeline.png
    → simulation_cell_methylation_timeline.png
    → simulation_gene_jsd_timeline.png
    → simulation_gene_methylation_timeline.png
    ↓
Stage 4: Extract batch comparisons
    → batch_comparison_cell.csv (9 rows)
    → batch_comparison_gene.csv (180 rows)
    ↓
Stage 5: Generate comparison plots (4 plots)
    → cell_jsd_comparison.png
    → cell_methylation_comparison.png
    → gene_jsd_comparison.png
    → gene_methylation_comparison.png
    ↓
Stage 6: Generate per-gene plots (40 plots optional)
    → per_gene/gene_{0-19}_{jsd,methylation}_comparison.png
```

**Key Design Principles**:
- **CSV-first workflow**: Extract once, plot many times - avoids reprocessing large JSON files
- **Generic functions**: Single plotting function handles all metrics via column parameters
- **Dual CSV architecture**: Separate files for cell-level (9 rows) and gene-level (180 rows) data
- **Modular stages**: Each stage can run independently for debugging/customization
- **Auto-detection**: Automatically finds Phase 1 simulation in parent directory
- **Column consistency**:
  - Cell metrics: `value_0`, `value_1`, ... (cells can die/divide)
  - Gene metrics: `gene_0`, `gene_1`, ... (fixed gene identities)

**Output Structure**:
```
<phase2-dir>/results/
├── tables/                    # All extracted CSV data
│   ├── simulation_cell_jsd.csv
│   ├── simulation_cell_methylation.csv
│   ├── simulation_gene_jsd.csv
│   ├── simulation_gene_methylation.csv
│   ├── year30_cells.csv
│   ├── year50_cells.csv
│   ├── batch_comparison_cell.csv
│   └── batch_comparison_gene.csv
└── plots/                     # All generated plots
    ├── year30_cell_jsd_histogram.png
    ├── year30_cell_methylation_histogram.png
    ├── year50_cell_jsd_histogram.png
    ├── year50_cell_methylation_histogram.png
    ├── simulation_cell_jsd_timeline.png
    ├── simulation_cell_methylation_timeline.png
    ├── simulation_gene_jsd_timeline.png
    ├── simulation_gene_methylation_timeline.png
    ├── cell_jsd_comparison.png
    ├── cell_methylation_comparison.png
    ├── gene_jsd_comparison.png
    ├── gene_methylation_comparison.png
    └── per_gene/               # Optional Stage 6
        ├── gene_0_jsd_comparison.png
        ├── gene_0_methylation_comparison.png
        └── ... (38 more plots)
```

### Key Classes
- `Cell`: Individual cell with methylation state and JSD calculations
- `PetriDish`: Population of cells with growth/homeostasis dynamics

### Directory Structure
```
phase1/data/gene_rates_*/size*-seed*-{phase1_timestamp}/      # Phase 1 directory
├── simulation.json                                           # Phase 1 simulation
└── snap{Y1}to{Y2}-growth{G}-quant{Q}x{C}-mix{M}u-seed{S}-{phase2_timestamp}/
    ├── snapshots/
    │   ├── year{N}_snapshot.json[.gz]     # Snapshots with year key wrapper
    │   └── metadata.json                   # Extraction metadata
    └── individuals/
        ├── mutant/                         # Quantile-sampled populations
        ├── control1/                       # Random-sampled populations
        ├── control2/                       # Pure snapshot populations
        ├── common_pool.json[.gz]          # Shared snapshot cells (array)
        └── mixing_metadata.json           # Mixing parameters
```

## Key Implementation Patterns

### Phase 3 Data Extraction Pattern
```python
# All phase3 scripts use this pattern to access individual_final
def load_petri_dish(filepath):
    with smart_open(filepath, 'r') as f:
        data = json.load(f)

    # CRITICAL: Use 'individual_final' for mixed population (~288 cells)
    # NOT 'history' which only contains pre-mixing growth (~68 cells)
    if 'individual_final' in data:
        cells = data['individual_final']['cells']  # CORRECT
    elif 'history' in data:
        # WARNING: This is pre-mixing data!
        last_year = max(history.keys())
        cells = history[last_year]['cells']  # INCOMPLETE
```

### Gene Metric Aggregation Pattern (phase3)
```python
# In extract_batch_comparison.py
# Extract per-gene metrics for all 20 genes
for gene_idx in range(n_genes):
    gene_jsds = petri.calculate_gene_jsd()
    gene_jsd = float(gene_jsds[gene_idx])
    gene_methylation = calculate_gene_methylation(cells, gene_idx)
    # Output one row per gene per individual

# In plot_comparison_generic.py
# Auto-aggregate when plotting gene data
if 'gene_index' in df.columns and df['gene_index'].nunique() > 1:
    # Multiple genes: aggregate by individual
    df = df.groupby(['batch', 'individual_id'])[value_column].mean().reset_index()
    # Result: 9 rows (one mean value per individual)
```

### CSV Column Naming Convention
```python
# Cell-level CSV (9 rows - 1 per individual)
"individual_id,batch,cell_jsd_mean,cell_methylation_mean"

# Gene-level CSV (180 rows - 20 per individual)
"individual_id,batch,gene_index,gene_jsd,gene_methylation"

# Timeline CSVs (variable columns based on cell count)
"year,value_0,value_1,..."  # Cells (count varies)
"year,gene_0,gene_1,..."    # Genes (always 20)
```

### Import Structure
```python
# Phase 2/3 must add phase1 to path
import sys
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'phase1'))
from cell import PetriDish, Cell
```

### Smart File I/O Pattern
```python
# Universal smart_open() function
def smart_open(filepath, mode='r'):
    """Open .json or .json.gz files transparently"""
    if filepath.endswith('.gz'):
        import gzip
        return gzip.open(filepath, mode + 't')  # Text mode for JSON
    return open(filepath, mode)

# Usage:
with smart_open('file.json.gz', 'r') as f:
    data = json.load(f)
```

## Phase 2 Data Flow - Key Implementation Details

### Snapshot Extraction (Stage 1-2)
- Loads Phase 1 simulation using `smart_open()` (handles .json/.json.gz)
- Accesses history with string keys: `history["23"]`, `history["38"]`
- Preserves year key wrapper in snapshots: `{"23": {...}}`
- Direct copy - NO transformation of cell data

### Individual Creation (Stage 3-5)
- **Mutant**: Quantile-based sampling (sorts by cell_jsd)
- **Control1**: Random uniform sampling
- Growth phase: exponential for N years, then homeostasis
- Common mixing: all individuals get identical snapshot cells
- Normalization: median - 0.5σ threshold, excludes/trims populations

### Control2 Creation (Stage 6)
- Pure second snapshot populations
- Uses common pool from Stage 5
- No growth phase (no cell_history)

### Critical JSON Structures

**Phase 1 History (string year keys):**
```json
{
    "history": {
        "0": {"cells": [...], "gene_jsd": [...]},
        "23": {...},
        "38": {...}
    }
}
```

**Snapshot Format (preserves year wrapper):**
```json
{
    "23": {
        "cells": [...],
        "gene_jsd": [...]
    }
}
```

**Common Pool (direct array):**
```json
[
    {"cpg_sites": [...], "gene_rate_groups": [...], "cell_jsd": 0.456, ...}
]
```

## Phase 3 Comparison Plot Customization

### plot_comparison_generic.py Parameters

The generic comparison plotting function supports extensive customization:

```python
plot_comparison_generic(
    csv_path='data.csv',           # Input CSV or DataFrame
    value_column='metric',         # Column to plot
    title='Plot Title',           # Plot title
    ylabel='Y-axis Label',        # Y-axis label
    output_path='plot.png',       # Output file path

    # Visual customization
    batch_colors={'mutant': '#1f77b4', 'control1': '#ff7f0e', 'control2': '#2ca02c'},
    figsize=(1200, 600),          # Width x Height in pixels
    scale=2,                      # DPI scale factor
    seed=42,                      # Random seed for consistent jitter

    # Y-axis range customization
    y_range_padding=0.1,          # Padding factor (0.1 = 10% padding on each side)
    y_range_min=None,             # Manual override for y-axis minimum
    y_range_max=None,             # Manual override for y-axis maximum

    # Y-axis tick customization
    y_tick_spacing=0.005,         # Manual tick spacing override
    auto_tick_density=True,       # Enable adaptive tick density

    # Gene data aggregation
    aggregation='mean',           # 'mean' or 'std' for gene data

    # Output control
    verbose=True                  # Print progress messages
)
```

### Y-Axis Range Padding

The y-axis range is automatically padded to prevent data points from being clipped at the edges:

- **Default padding**: 10% (`y_range_padding=0.1`)
- **How it works**: Adds padding equal to `data_range * padding_factor` on each side
- **Manual override**: Set `y_range_min` and/or `y_range_max` for exact control
- **Disable padding**: Set `y_range_padding=0`

### Adaptive Y-Axis Tick Spacing

The function automatically adjusts tick density based on data range:

| Data Range | Auto Tick Spacing | Example Use Case |
|------------|-------------------|------------------|
| < 0.01     | 0.001            | Very tight JSD clustering |
| 0.01-0.02  | 0.002            | Very tight methylation |
| 0.02-0.05  | 0.005            | Tight methylation range |
| 0.05-0.1   | 0.01             | Moderate spread |
| 0.1-0.2    | 0.02             | Common methylation (0.1-0.25) |
| 0.2-0.3    | 0.025            | Medium-wide range |
| 0.3-0.5    | 0.05             | Wide distribution |
| 0.5-1.0    | 0.1              | Full 0-1 range |
| > 1.0      | Auto (Plotly)    | Large values |

Override with `y_tick_spacing` for manual control or set `auto_tick_density=False` to disable.

## Plot Generation Workflow

### Quick Plot Generation Commands

```bash
# Generate ALL plots from an existing phase2 dataset
cd phase3
python run_pipeline.py --phase2-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/

# Generate only histogram plots
python plot_histogram_original.py \
    --snapshot-dir ../phase1/data/{phase1_dir}/{phase2_subdir}/snapshots/ \
    --output-dir plots/

# Plot simulation timeline from extracted CSV files
python plot_simulation_timeline.py results/tables/
```

### Plot Output Locations

All plots are saved as high-resolution PNG files (1200x600, scale=2):

```
# Phase 3 outputs
results/
├── plots/
│   ├── year30_cell_jsd_histogram.png
│   ├── year30_cell_methylation_histogram.png
│   ├── year50_cell_jsd_histogram.png
│   ├── year50_cell_methylation_histogram.png
│   ├── simulation_cell_jsd_timeline.png
│   ├── simulation_cell_methylation_timeline.png
│   ├── simulation_gene_jsd_timeline.png
│   ├── simulation_gene_methylation_timeline.png
│   ├── cell_jsd_comparison.png
│   ├── cell_methylation_comparison.png
│   ├── gene_jsd_comparison.png
│   ├── gene_methylation_comparison.png
│   └── per_gene/
│       ├── gene_0_jsd_comparison.png
│       ├── gene_0_methylation_comparison.png
│       └── ... (38 more plots)
└── tables/
    ├── simulation_cell_jsd.csv
    ├── simulation_cell_methylation.csv
    ├── simulation_gene_jsd.csv
    ├── simulation_gene_methylation.csv
    ├── year30_cells.csv
    ├── year50_cells.csv
    ├── batch_comparison_cell.csv
    └── batch_comparison_gene.csv
```

## Plot Types Generated (Phase 3 Pipeline)

The pipeline generates comprehensive analysis plots across two metrics (JSD and methylation) at two levels (cell and gene):

### Stage 2: Snapshot Histograms (4 plots)
- **Year 30 Cell JSD**: Blue histogram with mean line
- **Year 30 Cell Methylation**: Red histogram with mean line
- **Year 50 Cell JSD**: Blue histogram with mean line
- **Year 50 Cell Methylation**: Red histogram with mean line

### Stage 3: Timeline Plots (4 plots)
- **Cell JSD Timeline**: Evolution with percentile bands
- **Cell Methylation Timeline**: Evolution with percentile bands
- **Gene JSD Timeline**: 20 gene trajectories with percentile bands
- **Gene Methylation Timeline**: 20 gene trajectories with percentile bands

### Stage 5: Comparison Plots (4 plots)
- **Cell JSD Comparison**: Scatter plot across batches
- **Cell Methylation Comparison**: Scatter plot across batches
- **Gene JSD Comparison**: Mean across 20 genes per individual
- **Gene Methylation Comparison**: Mean across 20 genes per individual

### Stage 6: Per-Gene Plots (40 plots, optional)
- 20 genes × 2 metrics (JSD, methylation)
- Individual gene behavior across batches

### Plot Styling Conventions
- **Color by metric type**:
  - JSD plots: Blue (#1f77b4) with red mean lines
  - Methylation plots: Red (#d62728) with dark blue mean lines
- **Batch colors**: Mutant (blue), Control1 (orange), Control2 (green)
- **Figure size**: 1200x600 pixels, scale=2 for high DPI
- **Template**: plotly_white with gray gridlines
- **Percentile bands**: 5-95% (light), 25-75% (dark) for trajectories

## Gene-Level Comparison Plot Architecture

The phase3 pipeline generates gene-level comparison plots with flexible aggregation options:

### Data Flow for Gene Metrics
1. **Extraction** (`extract_batch_comparison.py`):
   - Loads each individual's PetriDish from `individual_final` (mixed population ~288 cells)
   - Calculates gene JSD using `petri.calculate_gene_jsd()` → 20 values per individual
   - Calculates gene methylation by averaging across cells for each gene region
   - Outputs 180 rows to `batch_comparison_gene.csv` (9 individuals × 20 genes)

2. **Plotting** (`plot_comparison_generic.py`):
   - Detects gene format by presence of `gene_index` column
   - Aggregates by taking mean OR std across all genes per individual
   - Results in 9 data points (3 per batch) for scatter plot
   - Displays mean line, 25-75% bands, 5-95% bands

### Generating Standard Deviation Plots

```bash
# Generated automatically by run_pipeline.py:
# Mean plots (default)
gene_jsd_mean_comparison.png       # Mean of gene JSDs across 20 genes
gene_methylation_mean_comparison.png  # Mean of gene methylations

# Standard deviation plots
gene_jsd_std_comparison.png        # Std dev of gene JSDs across 20 genes
gene_methylation_std_comparison.png   # Std dev of gene methylations

# Or manually with plot_comparison_generic.py:
python plot_comparison_generic.py \
    --csv batch_comparison_gene.csv \
    --column gene_jsd \
    --title "Gene JSD Standard Deviation Across Batches" \
    --ylabel "Std Dev of Gene JSD" \
    --output gene_jsd_std.png \
    --aggregation std  # Key parameter for std aggregation
```

### Understanding the Aggregation Logic

For gene-level data (180 rows: 9 individuals × 20 genes):
- **Mean aggregation**: For each individual, computes mean(gene_0, gene_1, ..., gene_19)
- **Std aggregation**: For each individual, computes std(gene_0, gene_1, ..., gene_19)

This measures within-individual variability across genes:
- High std → genes have divergent methylation/JSD patterns within that individual
- Low std → genes have uniform methylation/JSD patterns within that individual

## Project-Specific Instructions

### Critical Rules
- **ALWAYS use `python` instead of `python3`** for all commands
- **NO backward compatibility** unless explicitly requested - prefer clean breaks
- **Dictation note**: "jean" or "gin" means "gene"

### Configuration System
1. Default config loaded first (e.g., `config_default.yaml`)
2. Custom config file overrides defaults
3. Command-line arguments override everything

**Default Configurations:**
- Phase 1: `gene_rate_groups="5:0.004,5:0.005,5:0.006,5:0.007"`, growth_phase=9, years=50, sites=100
- Phase 2: snapshots (23, 38), quantiles=5, cells_per_quantile=2, growth=6, mix_ratio=82

### File I/O Patterns
- Always use `smart_open()` for .json/.json.gz handling
- Text mode for JSON operations (`'rt'`/`'wt'` for gzip)
- Deep copies when necessary to avoid reference issues
- Year keys always strings: `year_str = str(year)`

### Validation Points
- Cell compatibility: check `gene_rate_groups` consistency
- Population size: normalization handles homeostasis variation
- Mixing validation: ensure pool compatibility before mixing

## Installation & Dependencies

```bash
pip install -r requirements.txt
```

Required packages:
- `numpy` (>=1.19.0): Performance optimization
- `scipy` (>=1.7.0): Statistical analysis
- `pyyaml` (>=6.0): Configuration files
- `plotly` (>=5.0.0, <6.0.0): Interactive visualization
- `kaleido` (0.2.1): PNG export from plotly

## Common Issues & Solutions

### Phase 3 Specific Issues
- **Wrong cell count in comparisons**: Ensure using `individual_final` not `history`
  - `individual_final`: ~288 cells (mixed population) ✓
  - `history[last_year]`: ~68 cells (pre-mixing) ✗
- **Missing timeline plots**: Check Phase 1 simulation exists in parent directory
- **Gene CSV has wrong rows**: Should be 180 (9 individuals × 20 genes)
- **Comparison plots look sparse**: Normal - only 9 data points (3 per batch)

### Phase 2 Issues
- **Memory**: Reduce `n_quantiles` or `cells_per_quantile`
- **"Gene rate groups mismatch"**: Inconsistent simulation
- **"Insufficient cells"**: Snapshot too small for sampling
- **High variation**: Normalization automatically applied

### General Issues
- **ImportError**: Check sys.path additions
- **No plots**: Install `plotly kaleido`
- **Slow simulation**: Reduce sites or years
- **Out of memory**: Use compression or smaller parameters

## Performance Tips
- Extract tables once, plot multiple times
- Use CSV for large tables (faster than JSON)
- Cache extracted data in tables/ directory
- Parallelize extraction for multiple individuals

## Testing

### Limited Test Infrastructure
```bash
# Phase 3 test (only test file in codebase)
cd phase3
python test_gene_extraction.py
```

**Note**: Most test files were removed in the 2025-01-20 cleanup. Only `test_gene_extraction.py` remains, which tests:
- Mock data generation (9 individuals: 3 per batch)
- CSV extraction and validation
- Gene index completeness
- Plot generation verification

## Recent Updates

### January 2025
- Phase 3 is the primary analysis pipeline (renamed from phase3-prime)
- Major refactoring: ~100 files reduced to 29 Python files
- Clear phase separation: simulation → data generation → analysis
- Lean JSON format with ~90% size reduction
- Test infrastructure largely removed (only 1 test file remains)

### September 2025
- Phase 2 outputs directly to Phase 1 simulation directory
- Common mixing always enabled (uniform mixing mode)
- Gene JSD trajectory plots improved with percentile bands
- Phase 3 streamlined as CSV-first analysis pipeline