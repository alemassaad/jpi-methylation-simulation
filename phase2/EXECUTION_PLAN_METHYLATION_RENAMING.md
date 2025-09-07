# Execution Plan: Cell Methylation Proportion Renaming

## Overview
Rename all references to "methylation" that represent cell-level methylation proportions to "cell_methylation_proportion" to distinguish from gene-level metrics.

## Phase 1: Core Cell Class Changes

### File: `phase1/cell.py`

#### 1. Cell Class Attributes (Line ~133)
**Change:**
```python
self.methylation_proportion = 0.0
```
**To:**
```python
self.cell_methylation_proportion = 0.0
```
**Why:** Primary attribute storing cell's methylation fraction

#### 2. Methylate Method Check (Line ~159)
**Change:**
```python
if self.methylation_proportion >= 1.0:
```
**To:**
```python
if self.cell_methylation_proportion >= 1.0:
```
**Why:** Check if cell is fully methylated

#### 3. Methylate Method Update (Line ~184)
**Change:**
```python
self.methylation_proportion = methylated_count / self.n
```
**To:**
```python
self.cell_methylation_proportion = methylated_count / self.n
```
**Why:** Update the proportion after methylation

#### 4. To Dict Method (Line ~241)
**Change:**
```python
'methylation_proportion': n_methylated / self.n,  # Proportion methylated
```
**To:**
```python
'cell_methylation_proportion': n_methylated / self.n,  # Proportion methylated
```
**Why:** Serialization key for saving to JSON

#### 5. From Dict Method (Line ~282)
**Change:**
```python
if 'methylation_proportion' in data:
    stored_prop = data['methylation_proportion']
```
**To:**
```python
# Handle both old and new keys for backward compatibility during transition
if 'cell_methylation_proportion' in data:
    stored_prop = data['cell_methylation_proportion']
elif 'methylation_proportion' in data:
    # Legacy support - can remove after transition
    stored_prop = data['methylation_proportion']
```
**Why:** Deserialization with temporary backward compatibility

#### 6. PetriDish Method (Line ~1119)
**Change:**
```python
def _calculate_mean_methylation_proportion(self) -> float:
```
**To:**
```python
def _calculate_mean_cell_methylation_proportion(self) -> float:
```
**Why:** Helper method for calculating mean across all cells

#### 7. PetriDish Method Implementation (Line ~1120-1124)
**Change:**
```python
return sum(cell.methylation_proportion for cell in self.cells) / len(self.cells)
```
**To:**
```python
return sum(cell.cell_methylation_proportion for cell in self.cells) / len(self.cells)
```
**Why:** Access renamed attribute

#### 8. PetriDish Statistics (Line ~1144)
**Change:**
```python
'mean_methylation': self._calculate_mean_methylation_proportion(),
```
**To:**
```python
'mean_cell_methylation_proportion': self._calculate_mean_cell_methylation_proportion(),
```
**Why:** Statistics key for JSON output

#### 9. History Stats Processing (Lines ~1290-1303)
**Change:**
```python
if 'methylation_proportion' in cell:
    meth_values.append(cell['methylation_proportion'] * 100)
...
meth_values.append(getattr(cell, 'methylation_proportion', 0.0) * 100)
```
**To:**
```python
# Handle both old and new keys
if 'cell_methylation_proportion' in cell:
    meth_values.append(cell['cell_methylation_proportion'] * 100)
elif 'methylation_proportion' in cell:
    # Legacy support
    meth_values.append(cell['methylation_proportion'] * 100)
...
meth_values.append(getattr(cell, 'cell_methylation_proportion', 
                           getattr(cell, 'methylation_proportion', 0.0)) * 100)
```
**Why:** Handle both dictionary and object access with backward compatibility

#### 10. Plot Method Names (Line ~1528)
**Change:**
```python
def plot_methylation(self, title: str = None, output_path: str = None,
```
**To:**
```python
def plot_cell_methylation_proportion(self, title: str = None, output_path: str = None,
```
**Why:** Rename plotting method for clarity

#### 11. Plot Method Title (Line ~1545)
**Change:**
```python
title = title or "Methylation Timeline"
```
**To:**
```python
title = title or "Cell Methylation Proportion Timeline"
```
**Why:** Default title clarity

#### 12. Plot Method Y-Axis (Line ~1560)
**Change:**
```python
yaxis_title='Mean Methylation (%)',
```
**To:**
```python
yaxis_title='Mean Cell Methylation Proportion (%)',
```
**Why:** Y-axis label clarity

#### 13. Combined Plot Method (Line ~1914)
**Change:**
```python
self.plot_methylation(title, f"{base_path}_methylation.png")
```
**To:**
```python
self.plot_cell_methylation_proportion(title, f"{base_path}_cell_methylation_proportion.png")
```
**Why:** Call renamed method with updated filename

## Phase 2: Phase1 Scripts and Tests

### File: `phase1/run_simulation.py`

#### 1. Import Comment (Line ~393)
**Change:**
```python
meth_props = [cell.methylation_proportion for cell in petri_dish.cells]
```
**To:**
```python
meth_props = [cell.cell_methylation_proportion for cell in petri_dish.cells]
```

#### 2. Print Statements (Lines ~395-396)
**Change:**
```python
print(f"  Mean methylation: {statistics.mean(meth_props):.2%}")
print(f"  Median methylation: {statistics.median(meth_props):.2%}")
```
**To:**
```python
print(f"  Mean cell methylation proportion: {statistics.mean(meth_props):.2%}")
print(f"  Median cell methylation proportion: {statistics.median(meth_props):.2%}")
```

### File: `phase1/plot_history.py`

#### 1. Extract Values (Line ~84)
**Change:**
```python
meth_values = np.array([cell['methylation_proportion'] for cell in year_data]) * 100
```
**To:**
```python
# Handle both old and new keys
meth_values = []
for cell in year_data:
    if 'cell_methylation_proportion' in cell:
        meth_values.append(cell['cell_methylation_proportion'])
    elif 'methylation_proportion' in cell:
        meth_values.append(cell['methylation_proportion'])
    else:
        meth_values.append(0.0)
meth_values = np.array(meth_values) * 100
```

#### 2. Function Name (Line ~260)
**Change:**
```python
def create_methylation_plot(stats, filename):
```
**To:**
```python
def create_cell_methylation_proportion_plot(stats, filename):
```

#### 3. Plot Title (Line ~275)
**Change:**
```python
title='Methylation History',
```
**To:**
```python
title='Cell Methylation Proportion History',
```

#### 4. Y-Axis Title (Line ~288)
**Change:**
```python
yaxis_title='Methylation (%)',
```
**To:**
```python
yaxis_title='Cell Methylation Proportion (%)',
```

#### 5. Method Call (Lines ~707-708)
**Change:**
```python
meth_output = os.path.join(plots_dir, f"{base_name}_methylation.png")
plotter.plot_methylation(filename, meth_output)
```
**To:**
```python
meth_output = os.path.join(plots_dir, f"{base_name}_cell_methylation_proportion.png")
plotter.plot_cell_methylation_proportion(filename, meth_output)
```

### File: `phase1/test_new_format.py`

#### 1. Test Assertions (Lines ~27, 44, 208)
**Change:** Update all references to check for new key name
```python
_ = cell.cell_methylation_proportion  # Trigger calculation
assert 'cell_methylation_proportion' not in cell_dict, "Redundant field present"
# Comment updates about cell_methylation_proportion
```

### File: `phase1/tests/test_comprehensive.py`

#### All occurrences of `methylation_proportion` → `cell_methylation_proportion`
- Lines 283-284, 417, 430-431, 615-620, 708

### File: `phase1/tests/test_edge_cases.py`

#### All occurrences (Lines 78, 94)
**Change:**
```python
assert cell.methylation_proportion == 1.0
assert cell.methylation_proportion == 0.0
```
**To:**
```python
assert cell.cell_methylation_proportion == 1.0
assert cell.cell_methylation_proportion == 0.0
```

### File: `phase1/tests/test_small.py`

#### Print statement (Line 61)
**Change:**
```python
print("\nMethylation verification:")
```
**To:**
```python
print("\nCell methylation proportion verification:")
```

## Phase 3: Phase2 Pipeline Changes

### File: `phase2/run_pipeline.py`

#### 1. Import (Line ~52)
**Change:**
```python
plot_cell_methylation_histogram,
```
**To:**
```python
plot_cell_methylation_proportion_histogram,
```

#### 2. Import (Line ~54)
**Change:**
```python
analyze_cell_methylation_comparison,
```
**To:**
```python
analyze_cell_methylation_proportion_comparison,
```

#### 3. First Snapshot Plot (Lines 648-649)
**Change:**
```python
methylation_plot_path = os.path.join(results_dir, f"year{args.first_snapshot}_methylation_histogram.png")
plot_cell_methylation_histogram(first_snapshot_cells, args.bins, methylation_plot_path,
```
**To:**
```python
methylation_plot_path = os.path.join(results_dir, f"year{args.first_snapshot}_cell_methylation_proportion_histogram.png")
plot_cell_methylation_proportion_histogram(first_snapshot_cells, args.bins, methylation_plot_path,
```

#### 4. Second Snapshot Plot (Lines 831-832)
**Change:**
```python
methylation_plot_path_50 = os.path.join(results_dir, f"year{args.second_snapshot}_methylation_histogram.png")
plot_cell_methylation_histogram(second_snapshot_cells, args.bins, methylation_plot_path_50,
```
**To:**
```python
methylation_plot_path_50 = os.path.join(results_dir, f"year{args.second_snapshot}_cell_methylation_proportion_histogram.png")
plot_cell_methylation_proportion_histogram(second_snapshot_cells, args.bins, methylation_plot_path_50,
```

#### 5. Individual Trajectories - Mutant (Lines 1234, 1238)
**Change:**
```python
meth_path = os.path.join(individual_trajectories_dir, f"mutant_{i:02d}_methylation.png")
plotter.plot_methylation(f"Mutant Individual {i:02d}", meth_path)
```
**To:**
```python
meth_path = os.path.join(individual_trajectories_dir, f"mutant_{i:02d}_cell_methylation_proportion.png")
plotter.plot_cell_methylation_proportion(f"Mutant Individual {i:02d} Cell Methylation Proportion", meth_path)
```

#### 6. Individual Trajectories - Control1 (Lines 1263, 1267)
**Change:**
```python
meth_path = os.path.join(individual_trajectories_dir, f"control1_{i:02d}_methylation.png")
plotter.plot_methylation(f"Control1 Individual {i:02d}", meth_path)
```
**To:**
```python
meth_path = os.path.join(individual_trajectories_dir, f"control1_{i:02d}_cell_methylation_proportion.png")
plotter.plot_cell_methylation_proportion(f"Control1 Individual {i:02d} Cell Methylation Proportion", meth_path)
```

#### 7. Analysis Call (Line 1301)
**Change:**
```python
methylation_results = analyze_cell_methylation_comparison(
```
**To:**
```python
methylation_results = analyze_cell_methylation_proportion_comparison(
```

#### 8. Dict Creation for Compatibility (Line 1469)
**Change:**
```python
'methylation_proportion': sum(cell_dict['methylated']) / len(cell_dict['methylated']),
```
**To:**
```python
'cell_methylation_proportion': sum(cell_dict['methylated']) / len(cell_dict['methylated']),
```

#### 9. Original Timeline (Lines 1493-1495)
**Change:**
```python
meth_timeline_path = os.path.join(results_dir, "original_simulation_methylation_timeline.png")
plotter.plot_methylation("Original Simulation Methylation Timeline", meth_timeline_path)
print(f"  ✓ Generated original_simulation_methylation_timeline.png")
```
**To:**
```python
meth_timeline_path = os.path.join(results_dir, "original_simulation_cell_methylation_proportion_timeline.png")
plotter.plot_cell_methylation_proportion("Original Simulation Cell Methylation Proportion Timeline", meth_timeline_path)
print(f"  ✓ Generated original_simulation_cell_methylation_proportion_timeline.png")
```

### File: `phase2/core/pipeline_analysis.py`

#### 1. Function Definition (Line 171)
**Change:**
```python
def plot_cell_methylation_histogram(cells: List[Cell], bins: int, output_path: str,
```
**To:**
```python
def plot_cell_methylation_proportion_histogram(cells: List[Cell], bins: int, output_path: str,
```

#### 2. Extract Proportions (Lines 187-196)
**Change:**
```python
methylation_props = []
for cell in cells:
    if isinstance(cell, dict):
        prop = sum(cell.get('methylated', [])) / len(cell.get('methylated', [1]))
    else:
        prop = cell.methylation_proportion
    methylation_props.append(prop)
```
**To:**
```python
cell_methylation_props = []
for cell in cells:
    if isinstance(cell, dict):
        prop = sum(cell.get('methylated', [])) / len(cell.get('methylated', [1]))
    else:
        # Try new attribute first, fall back to old
        prop = getattr(cell, 'cell_methylation_proportion', 
                      getattr(cell, 'methylation_proportion', 0.0))
    cell_methylation_props.append(prop)
```

#### 3. Statistics (Lines 199-202)
**Change:**
```python
mean_meth = np.mean(methylation_props)
std_meth = np.std(methylation_props)
median_meth = np.median(methylation_props)
min_meth = np.min(methylation_props)
max_meth = np.max(methylation_props)
```
**To:**
```python
mean_meth = np.mean(cell_methylation_props)
std_meth = np.std(cell_methylation_props)
median_meth = np.median(cell_methylation_props)
min_meth = np.min(cell_methylation_props)
max_meth = np.max(cell_methylation_props)
```

#### 4. Plot Labels (Lines 212-215)
**Change:**
```python
title = title or f"Cell Methylation Distribution (Year {year})" if year else "Cell Methylation Distribution"
...
xaxis_title='Methylation Proportion',
```
**To:**
```python
title = title or f"Cell Methylation Proportion Distribution (Year {year})" if year else "Cell Methylation Proportion Distribution"
...
xaxis_title='Cell Methylation Proportion',
```

#### 5. Function Definition (Line 440)
**Change:**
```python
def analyze_cell_methylation_comparison(mutant_dishes: List[PetriDish],
```
**To:**
```python
def analyze_cell_methylation_proportion_comparison(mutant_dishes: List[PetriDish],
```

#### 6. Helper Function (Line 459)
**Change:**
```python
def get_mean_methylation_from_dishes(dishes: List[PetriDish]) -> np.ndarray:
```
**To:**
```python
def get_mean_cell_methylation_proportion_from_dishes(dishes: List[PetriDish]) -> np.ndarray:
```

#### 7. Helper Implementation (Lines 467-472)
**Change:**
```python
mean_vals.append(dish._calculate_mean_methylation_proportion())
```
**To:**
```python
mean_vals.append(dish._calculate_mean_cell_methylation_proportion())
```

#### 8. Function Calls (Lines 475-477)
**Change:**
```python
mutant_means = get_mean_methylation_from_dishes(mutant_dishes)
control1_means = get_mean_methylation_from_dishes(control1_dishes)
control2_means = get_mean_methylation_from_dishes(control2_dishes)
```
**To:**
```python
mutant_means = get_mean_cell_methylation_proportion_from_dishes(mutant_dishes)
control1_means = get_mean_cell_methylation_proportion_from_dishes(control1_dishes)
control2_means = get_mean_cell_methylation_proportion_from_dishes(control2_dishes)
```

#### 9. Results Dictionary Keys (Lines 487-502)
**Change all keys:**
```python
'mean_methylation' → 'mean_cell_methylation_proportion'
'std_methylation' → 'std_cell_methylation_proportion'
'median_methylation' → 'median_cell_methylation_proportion'
'min_methylation' → 'min_cell_methylation_proportion'
'max_methylation' → 'max_cell_methylation_proportion'
```

#### 10. Plot Title and Labels (Lines 525, 541)
**Change:**
```python
title='Cell Methylation Comparison'
yaxis_title='Methylation Proportion'
```
**To:**
```python
title='Cell Methylation Proportion Comparison'
yaxis_title='Cell Methylation Proportion'
```

#### 11. Output Path (Line 554)
**Change:**
```python
plot_path = os.path.join(output_dir, "cell_methylation_comparison.png")
```
**To:**
```python
plot_path = os.path.join(output_dir, "cell_methylation_proportion_comparison.png")
```

#### 12. JSON Output Path (Line 560)
**Change:**
```python
json_path = os.path.join(output_dir, "cell_methylation_analysis.json")
```
**To:**
```python
json_path = os.path.join(output_dir, "cell_methylation_proportion_analysis.json")
```

### File: `phase2/visualization/plot_individuals.py`

#### Method Calls (Lines 58, 708)
**Change:**
```python
plotter.plot_methylation(title, meth_path)
```
**To:**
```python
plotter.plot_cell_methylation_proportion(title, meth_path)
```

## Phase 4: Test File Updates

### All Test Files in `phase2/tests/`

Update all test files that reference methylation:
- `test_all_new_plots.py`: Update expected filenames
- `test_new_plots_implementation.py`: Update function names and expected files
- `test_validation.py`: Update parameter name in create_test_cells
- Unit tests in `phase2/tests/unit/`: Update any methylation references

## Phase 5: Documentation Updates

### File: `CLAUDE.md`

#### Update all sections mentioning methylation:
1. Plot documentation section
2. Generated files section
3. JSON structure examples
4. Method documentation

**Key changes:**
- "methylation_history.png" → "cell_methylation_proportion_history.png"
- "Mean methylation" → "Mean cell methylation proportion"
- Clarify distinction between cell and gene metrics

## Phase 6: Testing Plan

### 1. Unit Tests to Run
```bash
# Phase 1 tests
cd phase1/tests
python test_small.py
python test_comprehensive.py
python test_edge_cases.py
python test_new_format.py

# Phase 2 tests
cd phase2/tests
python test_all_new_plots.py
python test_new_plots_implementation.py
python test_validation.py

cd phase2/tests/unit
python test_gene_proportions.py
python test_gene_metrics.py
```

### 2. Integration Test
Run a small simulation and pipeline to verify:
```bash
# Phase 1
cd phase1
python run_simulation.py --rate 0.005 --years 10 --growth-phase 2 --sites 100

# Phase 2
cd phase2
python run_pipeline.py --simulation ../phase1/data/.../simulation.json --first-snapshot 5 --second-snapshot 10
```

### 3. Backward Compatibility Test
Test that old JSON files can still be read with the temporary compatibility code.

## Phase 7: Cleanup (After Verification)

Once all tests pass, remove backward compatibility code:
1. Remove legacy key handling in `from_dict` methods
2. Remove fallback attribute access
3. Update any remaining documentation

## Risk Assessment

**Low Risk Changes:**
- Plot titles and labels
- Print statements
- Documentation

**Medium Risk Changes:**
- Function/method names (need to update all calls)
- JSON keys (need backward compatibility)

**High Risk Changes:**
- Core Cell class attribute (affects entire codebase)
- Serialization/deserialization logic

## Rollback Plan

If issues arise:
1. Git revert the commit
2. Re-add compatibility shims if needed
3. Gradually migrate one component at a time

## Success Criteria

1. All tests pass
2. New simulations generate correctly named files
3. Old simulations can be read (compatibility mode)
4. Plots display correct labels
5. JSON files contain new keys
6. No ambiguity between cell and gene metrics