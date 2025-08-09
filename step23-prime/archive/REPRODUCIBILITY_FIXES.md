# Reproducibility Fixes Applied to step23-prime

## Issues Found and Fixed

### 1. ✅ No Global Random Seed (FIXED)
**Problem**: Pipeline started with undefined random state
**Fix**: Added to `run_pipeline.py` line 52-56:
```python
# Set global random seed for reproducibility
import random
import numpy as np
random.seed(args.seed)
np.random.seed(args.seed)
print(f"Random seed set to {args.seed} for reproducibility")
```

### 2. ✅ NumPy Random Not Seeded (FIXED)
**Problem**: Functions only seeded `random`, not `np.random`
**Fixes Applied**:

#### `pipeline_utils.py` - sample_by_quantiles (line 243)
```python
random.seed(seed)
np.random.seed(seed)  # Also seed numpy for full reproducibility
```

#### `pipeline_utils.py` - sample_uniform (line 298)
```python
random.seed(seed)
np.random.seed(seed)  # Also seed numpy for full reproducibility
```

#### `pipeline_utils.py` - mix_petri_with_snapshot (line 327-328)
```python
if seed is not None:
    random.seed(seed)
    np.random.seed(seed)  # Also seed numpy for full reproducibility
```

#### `pipeline_utils.py` - create_pure_snapshot_petri (line 370-371)
```python
if seed is not None:
    random.seed(seed)
    np.random.seed(seed)  # Also seed numpy for full reproducibility
```

## Places Where Randomness Occurs

### 1. Initial Sampling (Seeded ✅)
- `sample_by_quantiles()` - Selects cells from quantiles
- `sample_uniform()` - Selects cells uniformly
- Both now properly seed random and np.random

### 2. Growth Phase (Uses Global Seed ✅)
- `grow_petri_for_years()` → `PetriDish.divide_cells()` → `PetriDish.methylate_cells()`
- These use `Cell.methylate()` which calls `random.random()`
- Now controlled by global seed set at pipeline start

### 3. Mixing Phase (Seeded ✅)
- `mix_petri_with_snapshot()` - Samples year 60 cells
- `random.shuffle()` to mix cells
- Now properly seeds both random generators

### 4. Control2 Creation (Seeded ✅)
- `create_pure_snapshot_petri()` - Samples pure year 60
- Now properly seeds both random generators

### 5. Visualization (Already Seeded ✅)
- `pipeline_analysis.py` line 308: `np.random.seed(42)`
- Used for jittering points in plots (cosmetic only)

## Testing Reproducibility

### Quick Test
```bash
# Run twice with same seed
python run_pipeline.py --rate 0.005 --simulation ... --seed 42
mv data/rate_0.005000 data/rate_0.005000_run1

python run_pipeline.py --rate 0.005 --simulation ... --seed 42
mv data/rate_0.005000 data/rate_0.005000_run2

# Compare
diff -r data/rate_0.005000_run1 data/rate_0.005000_run2
```

### Full Test
```bash
python test_full_reproducibility.py
```

## What This Achieves

### Before Fixes:
- ❌ Each run produced different results even with same seed
- ❌ Could not reproduce bugs or verify results
- ❌ Tests would be non-deterministic

### After Fixes:
- ✅ Same seed → Same results (reproducible)
- ✅ Can write deterministic tests
- ✅ Can reproduce and debug issues
- ✅ Scientific validity through reproducibility

## Remaining Considerations

### 1. File System Order
Python 3.7+ preserves dict order, but file listing might vary by OS.
We use `sorted()` when listing files to ensure consistency.

### 2. Floating Point
Minor differences possible across platforms, but should be negligible.

### 3. Parallel Processing
Current implementation is serial, so no race conditions.

## Summary

**step23-prime is now fully reproducible!**

The pipeline will produce identical results when run with:
- Same input simulation
- Same parameters
- Same random seed

This ensures scientific validity and enables reliable testing.