# Detailed Fixing Plan for Phase 2 Pipeline Issues

## Date: 2025-09-07

## Executive Summary
This document provides detailed, step-by-step fixing plans for all identified issues in the Phase 2 pipeline.

**Key Decision**: Use 1-based indexing consistently throughout the pipeline.

---

## Fix 1: Gene JSD Reshape Error [CRITICAL]

### Problem
PetriDish objects created without parameters don't know the gene configuration, causing reshape errors.

### Root Cause Location
- File: `phase2/core/pipeline_analysis.py`
- Lines: 1333-1339
- Function: `plot_gene_jsd_distribution_comparison()`

### Detailed Fix Plan

#### Step 1: Extract gene configuration from cells
The cells already contain the gene configuration. We need to extract it.

**Location**: `phase2/core/pipeline_analysis.py`, before line 1333

**Add this code**:
```python
# Extract gene configuration from first cell
if snapshot1_cells and hasattr(snapshot1_cells[0], 'gene_rate_groups'):
    gene_rate_groups = snapshot1_cells[0].gene_rate_groups
    n_sites = len(snapshot1_cells[0].cpg_sites)
    gene_size = snapshot1_cells[0].gene_size if hasattr(snapshot1_cells[0], 'gene_size') else 5
else:
    # Fallback: calculate from cell data
    n_sites = len(snapshot1_cells[0].cpg_sites) if snapshot1_cells else 100
    gene_size = 5  # Default gene size
    gene_rate_groups = None
```

#### Step 2: Create PetriDish with proper parameters
**Replace lines 1333-1339** with:
```python
# Create PetriDish with proper configuration
if gene_rate_groups:
    petri1 = PetriDish(gene_rate_groups=gene_rate_groups, n=n_sites)
else:
    # Uniform rate - need to determine rate from somewhere
    # This might need to be passed as a parameter to the function
    petri1 = PetriDish(n=n_sites, gene_size=gene_size)

petri1.cells = snapshot1_cells
gene_jsds1 = petri1.calculate_gene_jsds()

# Same for petri2
if gene_rate_groups:
    petri2 = PetriDish(gene_rate_groups=gene_rate_groups, n=n_sites)
else:
    petri2 = PetriDish(n=n_sites, gene_size=gene_size)
    
petri2.cells = snapshot2_cells
gene_jsds2 = petri2.calculate_gene_jsds()
```

#### Step 3: Alternative approach - Pass configuration to function
**Better solution**: Modify function signature to accept gene configuration

**Change function signature** (line 1309):
```python
def plot_gene_jsd_distribution_comparison(snapshot1_cells, snapshot2_cells, 
                                          year1, year2, output_path,
                                          gene_rate_groups=None, gene_size=5):
```

**Update function call** in `run_pipeline.py` (around line 1351):
```python
plot_gene_jsd_distribution_comparison(
    snapshot1_cells, snapshot2_cells, 
    args.first_snapshot, args.second_snapshot, gene_dist_path,
    gene_rate_groups=gene_rate_groups, gene_size=gene_size
)
```

### Testing
1. Run pipeline with gene-specific rates
2. Verify no reshape error occurs
3. Check that gene JSD distribution plot is generated

---

## Fix 2: Consistent 1-Based Indexing [IMPORTANT]

### Problem
Mixed 0-based and 1-based indexing throughout pipeline causes confusion.

### Decision
**Use 1-based indexing everywhere** for user-facing output and file names.

### Detailed Fix Plan

#### Location 1: Normalization output
**File**: `phase2/core/pipeline_utils.py`
**Lines**: 1240, 1245, 1254, 1268, 1273, 1282

**Change all occurrences from**:
```python
print(f"    Mutant {i:02d}: ...")  # 0-based
```

**To**:
```python
print(f"    Mutant {i+1:02d}: ...")  # 1-based
```

**Specific changes**:
- Line 1240: `f"    Mutant {i+1:02d}: {current_size} cells - EXCLUDED (below threshold)"`
- Line 1245: `f"    Mutant {i+1:02d}: {current_size} cells - kept as is"`
- Line 1254: `f"    Mutant {i+1:02d}: {current_size} → {threshold_size} cells (trimmed)"`
- Line 1268: `f"    Control1 {i+1:02d}: {current_size} cells - EXCLUDED (below threshold)"`
- Line 1273: `f"    Control1 {i+1:02d}: {current_size} cells - kept as is"`
- Line 1282: `f"    Control1 {i+1:02d}: {current_size} → {threshold_size} cells (trimmed)"`

#### Location 2: File saving after exclusion
**File**: `phase2/run_pipeline.py`
**Around lines 960-990** (where mixed individuals are saved)

**Current issue**: When individuals are excluded, file numbers have gaps.

**Solution**: Track mapping between original and new IDs

**Add before the mixing loop**:
```python
# Create ID mapping for excluded individuals
mutant_id_map = {}  # Maps new sequential ID to original ID
new_id = 1
for i, petri in enumerate(normalized_mutant_dishes):
    original_id = petri.metadata.get('individual_id', i + 1)
    mutant_id_map[new_id] = original_id
    petri.metadata['original_id'] = original_id
    petri.metadata['individual_id'] = new_id
    new_id += 1

# Same for control1
control1_id_map = {}
new_id = 1
for i, petri in enumerate(normalized_control1_dishes):
    original_id = petri.metadata.get('individual_id', i + 1)
    control1_id_map[new_id] = original_id
    petri.metadata['original_id'] = original_id
    petri.metadata['individual_id'] = new_id
    new_id += 1
```

**Update file saving** to use new sequential IDs:
```python
# Use petri.metadata['individual_id'] for file naming (sequential)
# Store petri.metadata['original_id'] for traceability
```

#### Location 3: Growth phase output
**Already uses 1-based** - No change needed

### Testing
1. Run pipeline with normalization that excludes some individuals
2. Verify output shows consistent 1-based indexing
3. Check file names are sequential without gaps
4. Verify metadata contains both original and new IDs

---

## Fix 3: Improve Growth Messages [MINOR]

### Problem
Message "Individual 01: 1 → ~64 cells" is misleading as final count varies.

### Fix Plan

**File**: `phase2/core/individual_helpers.py`
**Line**: 144

**Change from**:
```python
print(f"    Individual {individual_id:02d}: {current_cells} → ~{expected_population} cells")
```

**To**:
```python
print(f"    Individual {individual_id:02d}: Starting growth from {current_cells} cell(s)")
```

**Alternative option** (more informative):
```python
print(f"    Individual {individual_id:02d}: Growing for {years} years (target ~{expected_population} cells)")
```

---

## Fix 4: Update Verbose Message Terminology [MINOR]

### Problem
Some messages use outdated or inconsistent terminology.

### Fix Plan

#### Location 1: Gene JSD histogram message
**File**: Search for "Plotting Gene-level JSD histogram" in all files

**Find and replace**:
- From: "Plotting Gene-level JSD histogram..."
- To: "Plotting gene JSD histogram..."

#### Location 2: Other terminology
**Search for and update**:
1. "Gene-level" → "gene" (when referring to JSD)
2. "Cell-level" → "cell" (when referring to JSD)
3. Any remaining "methylation" that should be "cell_methylation_proportion"

**Use grep to find all occurrences**:
```bash
grep -r "Gene-level" phase2/
grep -r "Cell-level" phase2/
```

---

## Fix 5: Enhance Validation Messages [MINOR]

### Problem
Validation messages don't show what was validated.

### Fix Plan

#### Location 1: Snapshot validation
**File**: `phase2/core/validation.py`
**Function**: `validate_snapshot()`

**Current**:
```python
print(f"  [INFO] ✓ Snapshot validation passed")
```

**Change to**:
```python
print(f"  [INFO] ✓ Snapshot validation passed: {len(cells)} valid cells with gene_rate_groups: {gene_rate_groups}")
```

#### Location 2: Individual validation
**Various locations in validation.py**

**Add details to all validation success messages**:
- Show counts
- Show key parameters verified
- Show any important metrics

**Example changes**:
```python
# Instead of:
print(f"  [INFO] ✓ Initial individuals validation passed")

# Use:
print(f"  [INFO] ✓ Initial individuals validation passed: {mutant_count} mutant, {control1_count} control1")
```

---

## Implementation Order

### Phase 1: Critical Fix
1. **Fix Gene JSD reshape error** - Prevents functionality

### Phase 2: Important Fixes  
2. **Implement consistent 1-based indexing** - Reduces confusion
3. **Update validation messages** - Helps debugging

### Phase 3: Minor Improvements
4. **Improve growth messages** - Clarity
5. **Update verbose terminology** - Consistency

---

## Testing Plan

### Test Scenario 1: Gene JSD Fix
```bash
python run_pipeline.py \
    --gene-rate-groups "5:0.004,5:0.005,5:0.006,5:0.007" \
    --simulation [path] \
    --first-snapshot 30 --second-snapshot 50
```
- Verify no reshape error
- Check gene JSD distribution plot is created

### Test Scenario 2: Indexing Consistency
```bash
# Run with settings that cause some individuals to be excluded
python run_pipeline.py \
    --rate 0.005 \
    --simulation [path] \
    --normalize-size  # This often excludes some individuals
```
- Check all output shows 1-based indexing
- Verify file names are sequential (no gaps)
- Check metadata for ID mapping

### Test Scenario 3: Message Improvements
- Run standard pipeline
- Verify all messages are clear and informative
- Check terminology is consistent

---

## Rollback Plan

If any fix causes issues:
1. All changes are isolated and can be reverted independently
2. Gene JSD fix is most critical - test thoroughly
3. Indexing changes affect multiple files but are straightforward
4. Message changes are cosmetic and low risk

---

## Questions Resolved

Q: How to handle individual IDs after exclusion?
A: Use sequential 1-based IDs for files, preserve original IDs in metadata for traceability.

Q: What parameters to use for PetriDish initialization?
A: Extract from cells or pass as function parameters.

Q: How detailed should validation messages be?
A: Show counts and key parameters without being overly verbose.