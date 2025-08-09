# Random Seed Usage Analysis: step23 vs step23-prime

## Critical Differences Found

### 1. Main Pipeline Seeding

#### step23 (run_pipeline_v2.py)
```python
# Line 183-184: Seeds BOTH random and np.random at start
random.seed(args.seed)
np.random.seed(args.seed)
```

#### step23-prime (run_pipeline.py)
```python
# NO GLOBAL SEEDING in main pipeline!
# Seeds are only passed to individual functions
```

**Impact**: step23-prime doesn't set a global seed, relying on functions to seed themselves

### 2. Quantile Sampling

#### step23 (run_pipeline_v2.py)
```python
def sample_by_quantiles(..., seed=42):
    random.seed(seed)  # Line 48
    # ...
    sampled = random.sample(quantile_cells, cells_per_quantile)  # Line 73
```

#### step23-prime (pipeline_utils.py)
```python
def sample_by_quantiles(..., seed=42):
    random.seed(seed)  # Line 242
    # ...
    sample_indices = random.sample(range(len(quantile_cells)), cells_per_quantile)  # Line 269
```

**Difference**: 
- step23: `random.sample(quantile_cells, cells_per_quantile)` - samples cells directly
- step23-prime: `random.sample(range(len(quantile_cells)), cells_per_quantile)` - samples indices

This is subtle but could change the random sequence!

### 3. Uniform Sampling

#### step23 (pipeline_utils.py)
```python
def sample_uniform(cells, n_cells=30, seed=None):
    if seed is not None:
        random.seed(seed)      # Line 165
        np.random.seed(seed)   # Line 166 - ALSO seeds numpy!
    
    sampled_indices = random.sample(range(len(cells)), n_cells)
```

#### step23-prime (pipeline_utils.py)
```python
def sample_uniform(cells, n_samples=30, seed=42):
    random.seed(seed)  # Line 296 - ONLY seeds random, not numpy!
    
    sample_indices = random.sample(range(len(cells)), n_samples)
```

**Critical Difference**: step23 seeds BOTH `random` and `np.random`, step23-prime only seeds `random`

### 4. Growth Phase

#### step23 (pipeline_utils.py)
```python
def grow_individual_inplace(filepath, years, rate=None):
    # NO SEEDING during growth!
    # Uses whatever random state is current
    for year in range(years):
        # ... division and aging using Cell.age_1_year()
```

#### step23-prime (pipeline_utils.py)
```python
def grow_petri_for_years(petri: PetriDish, years: int, verbose: bool = True):
    # NO SEEDING during growth!
    # Also uses whatever random state is current
    for year in range(years):
        petri.divide_cells()
        petri.methylate_cells()
```

### 5. Mixing with Year 60

#### step23 (pipeline_utils.py)
```python
def mix_with_year60_inplace(filepath, year60_cells, mix_ratio=80, seed=None):
    if seed is not None:
        random.seed(seed)      # Line 320
        np.random.seed(seed)   # Line 321 - Seeds both!
    
    sampled_indices = random.sample(range(len(year60_cells)), n_to_add)
```

#### step23-prime (pipeline_utils.py)
```python
def mix_petri_with_snapshot(petri, snapshot_cells, mix_ratio=0.8, seed=None):
    if seed is not None:
        random.seed(seed)  # Line 325 - Only seeds random!
    
    sample_indices = random.sample(range(len(snapshot_cells)), n_to_add)
```

## Key Issues That Affect Reproducibility

### Issue 1: Global vs Local Seeding
- **step23**: Sets global seeds at start (`random.seed(42)` and `np.random.seed(42)`)
- **step23-prime**: No global seeding, only per-function

### Issue 2: NumPy Random Not Seeded
- **step23**: Always seeds both `random` and `np.random` together
- **step23-prime**: Only seeds `random`, never `np.random`

### Issue 3: Different Sampling Methods
- **step23**: `random.sample(cells, n)` in quantile sampling
- **step23-prime**: `random.sample(range(len(cells)), n)` then index

### Issue 4: Seed Propagation
In step23-prime, seeds are passed like:
```python
# mutant sampling
sampled = sample_by_quantiles(year50_cells, seed=args.seed)

# control1 sampling  
sampled_cells = sample_uniform(year50_cells, seed=args.seed + 1000)

# mixing
mix_petri_with_snapshot(petri, year60_cells, seed=args.seed + 100 + i)
```

But without global seeding, the random state between calls is undefined!

## Why This Matters

The lack of global seeding in step23-prime means:
1. Random state is undefined at start
2. Each function call might start with different state
3. Any numpy random calls use unseeded generator

## Fix to Make Them Identical

To make step23-prime produce identical results, add at the start of `run_pipeline()`:

```python
# In run_pipeline() function, right after setup:
import random
import numpy as np

random.seed(args.seed)
np.random.seed(args.seed)
```

And ensure all utility functions seed both random generators:
```python
if seed is not None:
    random.seed(seed)
    np.random.seed(seed)  # Add this line!
```