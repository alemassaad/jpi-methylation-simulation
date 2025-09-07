# Phase 2 Pipeline Investigation Findings

## Investigation Date: 2025-09-07

## Summary
Investigated issues found in pipeline output to determine root causes and fix approaches.

## Issue 1: PlotPaths Integration - ✅ WORKING CORRECTLY

### Finding:
The new PlotPaths batch subdirectory structure IS working correctly. Files are being saved in the proper locations:
- `results/cell_metrics/individual_trajectories/jsd/mutant/individual_01.png`
- `results/cell_metrics/individual_trajectories/jsd/control1/individual_01.png`

### Evidence:
- run_pipeline.py correctly calls PlotPaths methods (lines 1246-1248, 1275-1277)
- Actual file output confirms correct structure
- The verbose output just doesn't show full paths

### Status: **Non-issue** - Working as intended

---

## Issue 2: Gene JSD Reshape Error - 🐛 CRITICAL BUG FOUND

### Root Cause:
In `pipeline_analysis.py` line 1333-1339, PetriDish objects are created with empty constructor:
```python
petri1 = PetriDish()  # ❌ No parameters!
petri1.cells = snapshot1_cells
gene_jsds1 = petri1.calculate_gene_jsds()
```

### Problem:
- PetriDish() with no parameters defaults to some configuration
- It doesn't know the actual gene_rate_groups or gene configuration
- When calculate_gene_jsds() tries to reshape, it expects wrong dimensions
- Error: "cannot reshape array of size 46100 into shape (461,200,5)"
  - 46100 = 461 cells × 100 sites/cell ✓
  - But expects (461, 200, 5) = 461,000 ❌
  - The 200 is wrong (should be 20 genes)

### Fix Required:
Must pass proper parameters when creating PetriDish:
```python
# Need to get gene_rate_groups from cells or global config
petri1 = PetriDish(gene_rate_groups=gene_rate_groups)
petri1.cells = snapshot1_cells
```

---

## Issue 3: ID/Indexing Inconsistency - 📊 DESIGN INCONSISTENCY

### Finding:
Multiple indexing schemes are used throughout the pipeline:

1. **Growth Phase** (individual_helpers.py line 140):
   - Uses `individual_id` from metadata, defaults to `i + 1` (1-based)
   - Output: "Individual 01", "Individual 02", etc.

2. **Normalization** (pipeline_utils.py line 1234):
   - Uses `enumerate(mutant_dishes)` directly (0-based)
   - Output: "Mutant 00", "Mutant 01", etc.

3. **File Saving**:
   - Uses 1-based indexing for filenames
   - individual_01.json, individual_02.json, etc.

### Problem:
When individuals are excluded during normalization:
- Original IDs: 01, 02, 03, 04, 05, 06, 07, 08
- After excluding 04: Still saves as 01, 02, 03, 05, 06, 07, 08
- This creates gaps in numbering and confusion

### Design Question:
Should we:
1. Preserve original IDs (with gaps) for traceability?
2. Renumber consecutively after exclusion?

---

## Issue 4: Misleading Growth Messages - 📝 MINOR

### Issue:
Output shows "Individual 01: 1 → ~64 cells" but actual final counts vary wildly (37, 94, 111, etc.)

### Root Cause:
The "~64" is the expected population, but homeostasis creates high variance.

### Fix:
Change message to be less specific:
- "Individual 01: Growing for 20 years..."
- Or: "Individual 01: 1 cell → homeostasis"

---

## Issue 5: Verbose Message Terminology - 📝 MINOR

### Issues Found:
1. "Plotting Gene-level JSD histogram..." should be "Plotting gene JSD histogram..."
2. Several messages still use old terminology

### Locations:
- Various print statements throughout the code
- Easy search and replace fix

---

## Issue 6: Validation Messages Lack Detail - 📝 MINOR

### Issue:
Messages like "✓ Validation passed" don't show what was validated

### Fix:
Add counts and details:
- "✓ Validation passed: 7 mutant, 6 control1 individuals"
- "✓ Validation passed: 461 cells, all with valid gene_rate_groups"

---

## Priority for Fixes

### Critical (Breaking):
1. **Gene JSD reshape error** - Causes analysis to fail

### Important (Confusing):
2. **ID/indexing inconsistency** - Creates user confusion

### Minor (Cosmetic):
3. Misleading growth messages
4. Verbose message terminology
5. Validation message details

---

## Recommendations

1. **Immediate Fix**: The Gene JSD reshape bug needs fixing ASAP as it breaks functionality

2. **Design Decision Needed**: The ID/indexing issue needs a decision:
   - Option A: Keep original IDs with gaps (current behavior)
   - Option B: Renumber consecutively after exclusion
   - Option C: Use a mapping system to track both

3. **Minor Fixes**: Can be done as a batch of small improvements

---

## Additional Notes

- The new PlotPaths directory structure is working perfectly
- The pipeline completes successfully despite the Gene JSD error (it's caught and skipped)
- Control2 numbering appears correct (adjusts count based on excluded individuals)