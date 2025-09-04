# Implementation Plan for New Plots

## Overview
This document details the implementation plan for adding 5 new plot types to the phase2 pipeline.

## 1. Gene-level JSD Snapshot Histogram

### Implementation Approach
**Standalone function** in `phase2/core/pipeline_analysis.py`

### Function to Create
```python
def plot_gene_jsd_snapshot_histogram(snapshot_cells, output_path, year=None, bins=20)
```

### Changes Required

#### In `phase2/core/pipeline_analysis.py`:
- **Add new function** after existing `plot_cell_jsd_distribution()`
- Function will:
  1. Create temporary PetriDish from snapshot cells
  2. Call `petri.calculate_gene_jsd()` to get 20 gene-level JSD values
  3. Create histogram with 20 bins (not 200 like cell-level)
  4. Add statistics box with mean, median, std of the 20 values
  5. Use similar styling to existing histograms

#### In `phase2/run_pipeline.py`:
- **After Stage 2** (line ~645): Add call to plot year 30 gene-level JSD histogram
- **After Stage 5.5** (line ~818): Add call to plot year 50 gene-level JSD histogram
- Import the new function at top

### Output Files
- `results/year30_gene_jsd_histogram.png`
- `results/year50_gene_jsd_histogram.png`

---

## 2. Gene-level JSD Individual Trajectories

### Implementation Approach
**Extend PetriDishPlotter class** in `phase1/cell.py`

### Method to Add
```python
class PetriDishPlotter:
    def plot_gene_jsd_trajectory(self, title, output_path)
```

### Changes Required

#### In `phase1/cell.py`:
- **Add new method** to PetriDishPlotter class (~line 800)
- Method will:
  1. Loop through `self.petri.cell_history`
  2. For each year, create temporary PetriDish with those cells
  3. Call `calculate_gene_jsd()` to get 20 values
  4. Calculate mean and median of the 20 values
  5. Store both for plotting
  6. Create plot with two lines (mean in blue, median in red/dashed)
  7. Use same styling as existing trajectory plots

#### In `phase2/run_pipeline.py`:
- **In trajectory plotting section** (~line 1185): Add gene-level JSD trajectory generation
- Loop through mutant and control1 individuals
- For each, call `plotter.plot_gene_jsd_trajectory()`
- Skip control2 (no trajectory)

### Output Files
- `results/individual_trajectories/mutant_XX_gene_jsd.png`
- `results/individual_trajectories/control1_XX_gene_jsd.png`

---

## 3. Gene-level JSD Original Timeline

### Implementation Approach
**Extend PetriDishPlotter class** in `phase1/cell.py`

### Method to Add
```python
class PetriDishPlotter:
    def plot_gene_jsd_timeline(self, title, output_path)
```

### Changes Required

#### In `phase1/cell.py`:
- **Add new method** to PetriDishPlotter class
- Method will:
  1. Check if `self.petri.gene_jsd_history` exists
  2. If yes: Use pre-calculated values directly
  3. If no: Loop through `self.petri.cell_history`, calculate gene JSDs
  4. For each year, store mean and median of 20 gene JSDs
  5. Create timeline plot with two lines
  6. Include percentile bands if desired (like existing timeline plots)

#### In `phase2/run_pipeline.py`:
- **In original simulation timeline section** (~line 1465): Add gene-level timeline
- After existing timeline plots, add:
  ```python
  gene_jsd_timeline_path = os.path.join(results_dir, "original_simulation_gene_jsd_timeline.png")
  plotter.plot_gene_jsd_timeline("Original Simulation Gene JSD Timeline", gene_jsd_timeline_path)
  ```

### Output Files
- `results/original_simulation_gene_jsd_timeline.png`

---

## 4. Cell Methylation Snapshot Histogram

### Implementation Approach
**Standalone function** in `phase2/core/pipeline_analysis.py`

### Function to Create
```python
def plot_cell_methylation_histogram(snapshot_cells, bins, output_path, year=None)
```

### Changes Required

#### In `phase2/core/pipeline_analysis.py`:
- **Add new function** after gene JSD functions
- Function will:
  1. For each cell, calculate `sum(cell.methylated) / len(cell.methylated)`
  2. Create list of methylation proportions
  3. Create histogram with 200 bins
  4. Add statistics box (mean, median, std, percentiles)
  5. Use same styling as cell JSD histogram
  6. X-axis will likely be narrower range (0.2-0.4 for year 30, 0.4-0.6 for year 50)

#### In `phase2/run_pipeline.py`:
- **After Stage 2**: Add call to plot year 30 methylation histogram
- **After Stage 5.5**: Add call to plot year 50 methylation histogram
- Place right after the gene JSD histogram calls

### Output Files
- `results/year30_methylation_histogram.png`
- `results/year50_methylation_histogram.png`

---

## 5. Cell Methylation Batch Comparison

### Implementation Approach
**Standalone function** in `phase2/core/pipeline_analysis.py`

### Function to Create
```python
def analyze_cell_methylation_comparison(mutant_dishes, control1_dishes, control2_dishes, output_dir)
```

### Changes Required

#### In `phase2/core/pipeline_analysis.py`:
- **Add new function** similar to `analyze_populations_from_dishes()`
- Function will:
  1. For each PetriDish, calculate mean methylation proportion across all cells
  2. Create scatter plot using existing `create_comparison_plot_from_jsds()` as template
  3. Rename to `create_comparison_plot_from_methylations()` or make generic
  4. Calculate t-tests between batches
  5. Save statistics to JSON
  6. Save plot

#### In `phase2/run_pipeline.py`:
- **In Stage 8 Analysis section** (~line 1270): Add methylation comparison
- Right after existing cell JSD analysis:
  ```python
  methylation_results = analyze_cell_methylation_comparison(
      mutant_dishes, control1_dishes, control2_dishes, results_dir
  )
  ```

### Output Files
- `results/cell_methylation_comparison.png`
- `results/cell_methylation_analysis.json`

---

## Implementation Order

1. **Cell Methylation Snapshot Histogram** - Easiest, similar to existing histogram
2. **Cell Methylation Batch Comparison** - Similar to existing batch comparison
3. **Gene-level JSD Snapshot Histogram** - Similar to #1 but with gene calculation
4. **Gene-level JSD Individual Trajectories** - Requires PetriDishPlotter extension
5. **Gene-level JSD Original Timeline** - Most complex, requires history handling

## Testing Strategy

After implementing each plot:
1. Run with small test data first
2. Verify output file is created
3. Check plot visually for correctness
4. Ensure no existing plots are broken
5. Update PLOT_DOCUMENTATION.md if needed

## Files Modified Summary

- `phase1/cell.py` - Add 2 methods to PetriDishPlotter
- `phase2/core/pipeline_analysis.py` - Add 3 new functions
- `phase2/run_pipeline.py` - Add 5 plot generation calls
- `PLOT_DOCUMENTATION.md` - Update status from 🆕 to ✅ after implementation