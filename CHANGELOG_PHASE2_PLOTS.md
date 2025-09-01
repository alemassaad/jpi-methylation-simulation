# Phase 2 Plot Generation Fixes and Enhancements

## Date: 2025-09-02

### Overview
Fixed critical issues preventing timeline plot generation and added new gene-level JSD individual comparison visualization.

### New Features
1. **Gene JSD Individual Comparison Plot** (`gene_jsd_individual_comparison.png`)
   - Shows individual-averaged gene JSD values across batches
   - Matches exact design of cell-level JSD comparison plots
   - Includes comprehensive statistics and quantile lines

### Fixed Issues

#### 1. Individual Growth Trajectory Plots
- **Problem**: Plots failed with "local variable 'PetriDishPlotter' referenced before assignment"
- **Cause**: Missing import statement within try blocks
- **Fix**: Added `from cell import PetriDishPlotter` where needed

#### 2. Original Simulation Timeline Plots
- **Problem**: Failed with cryptic '0' error for all simulations
- **Cause**: Multiple issues:
  - PetriDishPlotter expected string keys in cell_history but got integers
  - Lean JSON format incompatibility (uses 'methylated' not 'cpg_sites')
  - Missing n_sites and gene_size parameters for gene-specific rate simulations
- **Fixes**:
  - Keep year keys as strings in cell_history
  - Convert lean format cells to expected dictionary format
  - Pass n_sites and gene_size when creating PetriDish objects

#### 3. PetriDishPlotter Compatibility
- **Problem**: 'cell_jsd' KeyError when plotting
- **Cause**: Naming inconsistency between 'cell_jsd' and 'cell_JSD'
- **Fix**: Updated PetriDishPlotter to handle both naming conventions

#### 4. Critical Indentation Errors
- **Problem**: Entire run_pipeline function body incorrectly indented
- **Cause**: Previous editing mistake
- **Fix**: Corrected indentation for lines 384-1434

### Generated Plots
All plots now generate automatically without special flags:

**Timeline Plots:**
- `original_simulation_jsd_timeline.png` - Full simulation JSD evolution
- `original_simulation_methylation_timeline.png` - Methylation over time
- `original_simulation_combined_timeline.png` - Combined view
- `individual_trajectories/` - Growth trajectories for each individual

**Comparison Plots:**
- `cell_jsd_comparison.png` - Cell-level JSD comparison
- `gene_jsd_individual_comparison.png` - NEW: Gene-level individual averages
- `year{X}_jsd_distribution_{bins}bins.png` - Snapshot distributions

**Gene Plots:**
- `gene_jsd_plots/` - Per-gene distributions
- `top_variable_genes.png` - Most variable genes
- `gene_jsd_distribution.png` - Snapshot comparisons

### Files Modified
- `phase2/run_pipeline.py` - Fixed indentation, imports, and timeline plot generation
- `phase2/pipeline_analysis.py` - Added gene JSD individual comparison functions
- `phase1/cell.py` - Updated PetriDishPlotter for compatibility
- `CLAUDE.md` - Updated documentation

### Testing
All plots verified to generate correctly for:
- ✅ Uniform methylation rate simulations
- ✅ Gene-specific methylation rate simulations
- ✅ With and without normalization
- ✅ Various quantile and cell count configurations