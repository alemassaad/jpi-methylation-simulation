# Gene JSD Implementation Plan

This document outlines the comprehensive plan for implementing gene-level JSD plotting functionality in the DNA methylation simulation pipeline.

## Overview

We need to implement gene JSD plots to complement the existing cell JSD plots, ensuring clear distinction between:
- **Cell JSD**: Per-cell measure of methylation pattern divergence from baseline
- **Gene JSD**: Population-level measure of methylation heterogeneity for each gene

## PART 1: Fix Missing Cell JSD Labels

### 1.1 Fix `jsd_comparison.png` (in `pipeline_analysis.py`)
**Location:** `phase2/pipeline_analysis.py` - `create_comparison_plot_from_jsds()` function
**Changes Needed:**
- Update title from "Mean JSD per Individual" to "Mean Cell JSD per Individual"
- Update y-axis label to "Mean Cell JSD per Individual"
- Update hover text to specify "Cell JSD"

### 1.2 Fix Individual History Plots (in `cell.py`)
**Location:** `phase1/cell.py` - `PetriDishPlotter` class methods
**Changes Needed:**
- Already partially done, but need to verify all labels
- Check `plot_jsd()`, `plot_combined()` methods
- Ensure legend entries say "Cell JSD" not just "JSD"

## PART 2: Gene JSD Plot Implementations

All plots will follow the existing phase 2 design pattern:
- **Color scheme:** Same blues (#1f77b4), oranges (#ff7f0e), greens (#2ca02c)
- **Template:** `plotly_white`
- **Dimensions:** 1200x600 px (or 1200x800 for heatmaps)
- **Font:** Arial for annotations, size 16 for titles
- **Grid:** Light gray (rgba(0,0,0,0.1))
- **Annotations:** Top-right corner with statistics box

### Plot 1: Gene JSD Heatmap Over Time

**Function Name:** `plot_gene_jsd_heatmap()`
**Location:** Add to `phase1/cell.py` in `PetriDishPlotter` class
**Design Specs:**
```python
def plot_gene_jsd_heatmap(self, title: str = None, output_path: str = None,
                          width: int = 1200, height: int = 800):
    """
    Create heatmap showing each gene's JSD evolution over time.
    
    Layout:
    - X-axis: "Age (years)" 
    - Y-axis: "Gene Index"
    - Colorbar: "Gene JSD Score"
    - Title: "Gene JSD Evolution Heatmap"
    - Colorscale: 'Blues' (0=white, max=dark blue) to match existing plots
    """
```

**Data Structure:**
- 2D array: rows = genes (0-199), columns = years
- Values: JSD for each gene at each time point

**Special Features:**
- Add horizontal lines every 50 genes if using gene-rate-groups
- Annotation box showing: "Genes: 200, Years: X, Max JSD: X.XX"

### Plot 2: Gene JSD Distribution Comparison

**Function Name:** `plot_gene_jsd_distribution_comparison()`
**Location:** Add to `phase2/pipeline_analysis.py`
**Design Specs:**
```python
def plot_gene_jsd_distribution_comparison(snapshot1_cells, snapshot2_cells, 
                                          year1, year2, output_path):
    """
    Overlapping histograms comparing gene JSD distributions at two time points.
    
    Layout:
    - X-axis: "Gene JSD Score"
    - Y-axis: "Number of Genes"
    - Title: "Gene JSD Distribution: Year X vs Year Y"
    - Two overlapping histograms with transparency
    """
```

**Style Matching:**
- Use same histogram style as `plot_jsd_distribution_from_cells()`
- Step histogram with fill
- Year 1: Blue (#1f77b4) with 0.5 opacity
- Year 2: Orange (#ff7f0e) with 0.5 opacity
- 50 bins (fewer than cell JSD since only 200 genes)

**Statistics Box:**
```
Year X Stats:
Mean: X.XXX
Median: X.XXX
Max: X.XXX

Year Y Stats:
Mean: X.XXX
Median: X.XXX
Max: X.XXX
```

### Plot 3: Gene Rate Group Comparison

**Function Name:** `plot_gene_jsd_by_rate_group()`
**Location:** Add to `phase1/cell.py` in `PetriDishPlotter` class
**Design Specs:**
```python
def plot_gene_jsd_by_rate_group(self, title: str = None, output_path: str = None,
                                width: int = 1200, height: int = 600):
    """
    Line plot showing mean JSD for each gene rate group over time.
    
    Layout:
    - X-axis: "Age (years)"
    - Y-axis: "Mean Gene JSD Score"
    - Title: "Gene JSD by Methylation Rate Group"
    - One line per rate group
    """
```

**Style Details:**
- Line colors: Use color gradient (blue→red) for low→high rates
- Line width: 2px
- Include confidence bands (25-75 percentile) with 0.2 opacity
- Legend: "Group 1: X.X% rate" format

**Annotation Box:**
```
Gene Groups:
Group 1: 50 genes @ 0.4%
Group 2: 50 genes @ 0.45%
Group 3: 50 genes @ 0.5%
Group 4: 50 genes @ 0.55%
```

### Plot 4: Gene vs Cell JSD Scatter Plot

**Function Name:** `plot_gene_vs_cell_jsd_comparison()`
**Location:** Add to `phase2/pipeline_analysis.py`
**Design Specs:**
```python
def plot_gene_vs_cell_jsd_comparison(mutant_dishes, control1_dishes, control2_dishes,
                                     output_path):
    """
    Scatter plot comparing mean gene JSD vs mean cell JSD.
    
    Layout:
    - X-axis: "Mean Cell JSD per Individual"
    - Y-axis: "Mean Gene JSD per Individual"
    - Title: "Gene JSD vs Cell JSD Comparison"
    - Three groups with different colors/shapes
    """
```

**Style Details:**
- Mutant: Blue circles (#1f77b4)
- Control1: Orange triangles (#ff7f0e)
- Control2: Green squares (#2ca02c)
- Add diagonal reference line (y=x) in light gray
- Point size: 8px

**Statistics Box:**
```
Correlations:
Mutant: r=X.XX
Control1: r=X.XX
Control2: r=X.XX
Overall: r=X.XX
```

### Plot 5: Top Variable Genes Bar Chart

**Function Name:** `plot_top_variable_genes()`
**Location:** Add to `phase2/pipeline_analysis.py`
**Design Specs:**
```python
def plot_top_variable_genes(petri_dishes, n_top=20, output_path=None):
    """
    Bar chart showing the most heterogeneous genes.
    
    Layout:
    - X-axis: "Gene Index"
    - Y-axis: "Gene JSD Score"
    - Title: "Top 20 Most Heterogeneous Genes"
    - Horizontal bar chart for readability
    """
```

**Style Details:**
- Bars colored by gradient (light→dark blue based on JSD value)
- Include gene rate group info if applicable
- Sort by JSD value (highest at top)

## PART 3: Integration Points

### 3.1 Phase 1 Integration
**Files to modify:**
- `phase1/cell.py`: Add heatmap and rate group methods to `PetriDishPlotter`
- `phase1/run_simulation.py`: Ensure `--track-gene-jsd` flag works
- `phase1/plot_history.py`: Add option to generate gene JSD plots

### 3.2 Phase 2 Integration  
**Files to modify:**
- `phase2/pipeline_analysis.py`: Add comparison functions
- `phase2/run_pipeline.py`: Add calls to generate gene JSD plots
- Add new flag: `--plot-gene-jsd` to enable gene JSD plotting

### 3.3 File Naming Convention
Follow existing pattern:
- Phase 1: `{base_name}_gene_jsd_heatmap.png`, `{base_name}_gene_rate_comparison.png`
- Phase 2: `gene_jsd_distribution.png`, `gene_vs_cell_jsd.png`, `top_variable_genes.png`

## PART 4: Implementation Order

1. **First:** Fix existing Cell JSD labels (quick fixes)
2. **Second:** Implement Gene JSD Heatmap (most informative single plot)
3. **Third:** Implement Gene JSD Distribution Comparison (key for phase 2)
4. **Fourth:** Implement Gene Rate Group Comparison (if using gene-rate-groups)
5. **Fifth:** Implement Gene vs Cell JSD Scatter (correlation analysis)
6. **Sixth:** Implement Top Variable Genes (optional but insightful)

## PART 5: Testing Strategy

For each plot:
1. Test with small simulation (100 sites, 5 years)
2. Test with gene-rate-groups
3. Test with uniform rate
4. Verify file outputs match naming convention
5. Check plot readability and statistics accuracy

## Implementation Notes

### Key Design Principles
- Maintain visual consistency with existing phase 2 plots
- Clearly distinguish between Cell JSD and Gene JSD in all titles/labels
- Use consistent color schemes and layout patterns
- Include informative statistics boxes
- Follow existing file naming conventions

### Data Requirements
- Gene JSD tracking must be enabled with `--track-gene-jsd` flag
- Requires `gene_jsd_history` attribute in PetriDish objects
- For phase 2, need gene JSD data from grown individuals

### Color Palette
- Primary blue: #1f77b4
- Secondary orange: #ff7f0e  
- Tertiary green: #2ca02c
- Grid/background: rgba(0,0,0,0.1)
- Statistics box background: rgba(255,255,255,0.9)

## Expected Output Files

### Phase 1 (Single Simulation)
- `simulation_gene_jsd.png` (already implemented - evolution over time)
- `simulation_gene_jsd_heatmap.png` (new - detailed heatmap)
- `simulation_gene_rate_comparison.png` (new - if using gene-rate-groups)

### Phase 2 (Pipeline Analysis) 
- `gene_jsd_distribution.png` (new - year comparison histograms)
- `gene_vs_cell_jsd.png` (new - correlation scatter plot)
- `top_variable_genes.png` (new - bar chart of most heterogeneous genes)

This comprehensive plan ensures all gene JSD visualizations maintain consistency with the existing pipeline design while providing clear scientific insights into population-level methylation heterogeneity.