# Comparison Plan: step23 vs step23-prime

## What Are Individuals?

**Individuals = PetriDish objects** that represent experimental samples.

### Object Structure:
```python
# Each individual is a PetriDish containing:
individual = PetriDish(
    cells=[Cell, Cell, ...],     # List of Cell objects (5120 after mixing)
    year=60,                      # Current year
    rate=0.005,                   # Methylation rate
    metadata={                    # Extra tracking info
        'individual_id': 0,
        'individual_type': 'mutant',
        'source_quantile': 3,
        'mixed': True
    }
)
```

### Individual Life Cycle:
1. **Birth**: 1 cell sampled from year 50
2. **Growth**: 10 years of division (1→2→4→8→...→1024)
3. **Mixing**: Add year 60 cells (1024 + 4096 = 5120 total)
4. **Analysis**: Calculate mean JSD across all cells

## Comparison Strategy

### 1. Snapshots (Should Be IDENTICAL)
Snapshots are just extracted cells from the simulation. They should be exactly the same.

```bash
# Run snapshot comparison
python compare_snapshots.py
```

Expected result: ✅ Identical cells, same JSD values

### 2. Individuals (Should Be STATISTICALLY SIMILAR)
Individuals involve random sampling, so exact matches are unlikely, but statistics should be close.

```bash
# Run individual comparison
python compare_individuals.py
```

Expected results:
- ✅ Same number of individuals (30 per group)
- ✅ Same cell counts (5120 per individual)
- ✅ Similar mean JSD (within ~0.01)
- ⚠️ Individual cells will differ (random sampling)

## Quick Verification Commands

### Check Snapshots Are Identical:
```python
# Quick check that year 50 snapshots match
python -c "
import json, gzip, numpy as np

# Load both
with gzip.open('../step23/data/rate_0.005000/snapshots/year50_snapshot.json.gz', 'rt') as f:
    step23 = json.load(f)
with gzip.open('data/rate_0.005000/snapshots/year50_snapshot.json.gz', 'rt') as f:
    prime = json.load(f)

# Extract cells (handle format differences)
cells23 = step23 if isinstance(step23, list) else step23.get('cells', step23)
cells_prime = prime['cells'] if 'cells' in prime else prime

# Compare
jsds23 = sorted([c['jsd'] for c in cells23])
jsds_prime = sorted([c['jsd'] for c in cells_prime])

print(f'Cells match: {len(cells23) == len(cells_prime)}')
print(f'JSD match: {jsds23 == jsds_prime}')
"
```

### Check Individual Statistics:
```python
# Compare mean JSD for first mutant individual
python -c "
import json, gzip, numpy as np

# Load first mutant from each
with gzip.open('../step23/data/rate_0.005000/individuals/mutant/individual_00.json.gz', 'rt') as f:
    step23 = json.load(f)
with gzip.open('data/rate_0.005000/individuals/mutant/individual_00.json.gz', 'rt') as f:
    prime = json.load(f)

# Get cells
cells23 = step23['cells'] if 'cells' in step23 else step23
cells_prime = prime['cells'] if 'cells' in prime else prime

# Compare
print(f'Step23:       {len(cells23)} cells, mean JSD = {np.mean([c[\"jsd\"] for c in cells23]):.6f}')
print(f'Step23-prime: {len(cells_prime)} cells, mean JSD = {np.mean([c[\"jsd\"] for c in cells_prime]):.6f}')
"
```

## What to Expect

### Snapshots
- **Year 50 & 60**: Should be IDENTICAL (same cells, same order, same JSD values)
- If different: Check if same simulation file was used

### Individuals  
- **Cell counts**: Should all be 5120
- **Mean JSD per group**: Should be within ~0.01
- **Individual cells**: Will differ (random sampling with seed)

### Why Individuals Might Differ Slightly

Even with same seed (42), individuals can differ because:

1. **Random state management**: Python's random vs numpy.random
2. **Sampling order**: Which cells get sampled first
3. **Growth randomness**: Methylation during growth is stochastic
4. **Mixing randomness**: Which year 60 cells get selected

But the **statistical properties should be preserved**:
- Same mean JSD per group (±0.01)
- Same variance patterns
- Same p-values in comparisons

## Full Comparison Workflow

```bash
# 1. Wait for pipeline to complete
python run_pipeline.py --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz

# 2. Compare snapshots (should be identical)
python compare_snapshots.py

# 3. Compare individuals (should be statistically similar)
python compare_individuals.py

# 4. Compare final statistics
diff -u ../step23/data/rate_0.005000/results/statistics.json \
        data/rate_0.005000/results/statistics.json
```

## Success Criteria

✅ **Implementation is correct if:**
1. Snapshots are identical (same cells extracted)
2. Individual counts match (30 per group, 5120 cells each)
3. Group mean JSD within 0.01 of original
4. Overall statistical patterns preserved