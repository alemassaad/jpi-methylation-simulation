# Summary of Fixes Applied - September 7, 2025

## Overview
Fixed 5 additional issues discovered during final review of the pipeline output.

## Issues Fixed

### 1. ✅ File Overwrite Bug (CRITICAL)
**Problem**: Two different plots were saved to the same filename `gene_jsd_comparison.png`, causing the second to overwrite the first.

**Solution**: 
- Added new method `get_gene_jsd_snapshot_comparison_path()` to PlotPaths class
- Now saves snapshot comparison as `gene_jsd_snapshot_comparison.png`
- Individual-averaged comparison remains as `gene_jsd_comparison.png`

**Files Modified**:
- `core/plot_paths.py`: Added new path method
- `run_pipeline.py`: Updated to use new path for snapshot comparison

### 2. ✅ Confusing Timeline Description
**Problem**: Timeline mixed simulation years with individual growth years, stating "Years 30→36: Exponential growth" which was misleading.

**Solution**: Clarified the timeline to show:
- Growth happens over 20 years in individual's reference frame (years 0-20)
- Ages progress from 30→50 (simulation time)
- Clear separation between growth phases

**Files Modified**:
- `run_pipeline.py`: Rewrote timeline description section

**Before**:
```
📊 Individual Simulation Timeline:
   Year 30: Sample individual cells
   Years 30→36: Exponential growth (1 → 64 cells)
   Years 36→50: Homeostasis (~64 cells)
```

**After**:
```
📊 Individual Growth Timeline:
   Sampling: Individual cells from year 30 snapshot
   Growth period: 20 years (ages 30→50)
     • Years 0-6: Exponential growth (1 → 64 cells)
     • Years 7-20: Homeostasis (~64 cells)
```

### 3. ✅ Directory Naming Confusion
**Problem**: Output directory used "grow9" which referred to the original simulation's growth phase, not the pipeline's growth phase.

**Solution**: Changed to "simgrow9" to clarify it refers to the simulation parameter.

**Files Modified**:
- `core/path_utils.py`: Changed `grow{sim_params['growth_phase']}` to `simgrow{sim_params['growth_phase']}`

**Example**: `data/gene_rates_...-simgrow9-sites100-years50/...`

### 4. ⏭️ Individual Renumbering After Exclusion (SKIPPED)
**Note**: User requested to leave this as-is, so no changes were made to the renumbering behavior.

### 5. ✅ Control2 Count Explanation
**Problem**: Message "Adjusted control2 count after normalization: 6 (Based on 7 mutant + 6 control1)" didn't explain the logic clearly.

**Solution**: Clarified to show it's the average of the two groups.

**Files Modified**:
- `run_pipeline.py`: Improved message clarity

**Before**:
```
Adjusted control2 count after normalization: 6
  (Based on 7 mutant + 6 control1)
```

**After**:
```
Control2 count set to: 6
  (Average of 7 mutant + 6 control1 = 6.5)
```

## Testing
All fixes tested successfully:
- ✓ Two separate gene JSD comparison files created
- ✓ Timeline description is clearer
- ✓ Directory name shows "simgrow9" 
- ✓ Control2 count message is clearer

## No Breaking Changes
All changes are backward compatible and only affect:
- Output messages and plot filenames
- Directory naming for new runs
- No changes to data processing or analysis logic