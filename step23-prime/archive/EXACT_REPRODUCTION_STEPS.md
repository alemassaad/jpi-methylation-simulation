# Exact Steps to Reproduce Step23 with Step23-Prime (Including Timeout Workarounds)

## What Actually Happened (Including All Issues)

### Step 1: Initial Run (TIMED OUT)
```bash
cd step23-prime

python run_pipeline.py \
    --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --growth-years 10 \
    --mix-ratio 80 \
    --seed 42 \
    --bins 200
```
**Result**: Timed out after 10 minutes

### Step 2: Check What Was Created
```bash
# Check if directories were created
ls -la data/rate_0.005000/

# Found these were created:
# - snapshots/ (with year50 and year60 files)
# - individuals/mutant/ (30 files)
# - individuals/control1/ (30 files)
# - individuals/control2/ (EMPTY!)
```

### Step 3: Verify Partial Success
```bash
# Check snapshots
ls -la data/rate_0.005000/snapshots/
# ✓ year50_snapshot.json.gz (2.2MB)
# ✓ year60_snapshot.json.gz (2.4MB)

# Count individuals
ls data/rate_0.005000/individuals/mutant/ | wc -l     # 30 ✓
ls data/rate_0.005000/individuals/control1/ | wc -l   # 30 ✓
ls data/rate_0.005000/individuals/control2/ | wc -l   # 0 ✗

# Check if individuals were fully processed (grown and mixed)
python -c "import json, gzip; \
data=json.load(gzip.open('data/rate_0.005000/individuals/mutant/individual_00.json.gz', 'rt')); \
print(f'Cells: {len(data[\"cells\"])}')"
# Output: Cells: 5120 ✓ (correctly mixed)
```

### Step 4: Pipeline Got Stuck Creating Control2
The pipeline completed:
- ✓ Stage 1: Extract year 50 (cached)
- ✓ Stage 2: Plot JSD distribution 
- ✓ Stage 3: Create mutant and control1 individuals
- ✓ Stage 4: Grow individuals (1→1024 cells)
- ✓ Stage 5: Extract year 60 (cached)
- ✓ Stage 6: Mix populations (1024→5120 cells)
- ✗ Stage 7: Create control2 (TIMED OUT HERE)
- ✗ Stage 8: Analysis

### Step 5: Manually Create Control2 Individuals
Created a helper script `create_control2.py` but it also timed out. So did it in batches:

```python
# First attempt - timed out after creating 5
python create_control2.py  # Created 0-4, then timeout

# Second attempt - create more in batches
python -c "
import os, sys
sys.path.insert(0, '..')
sys.path.append('../step1-prime')
from pipeline_utils import load_snapshot_cells, create_pure_snapshot_petri, save_petri_dish

year60_cells = load_snapshot_cells('data/rate_0.005000/snapshots/year60_snapshot.json.gz')

# Create individuals 5-27 (timed out at 28)
for i in range(5, 30):
    print(f'Creating individual {i}...')
    petri = create_pure_snapshot_petri(year60_cells, n_cells=5120, rate=0.005, seed=42+300+i)
    petri.metadata = {'individual_id': i, 'individual_type': 'control2'}
    save_petri_dish(petri, f'data/rate_0.005000/individuals/control2/individual_{i:02d}.json.gz')
"
# This created up to individual 27, then timed out

# Final attempt - create last 2
python -c "
import os, sys
sys.path.insert(0, '..')
sys.path.append('../step1-prime')
from pipeline_utils import load_snapshot_cells, create_pure_snapshot_petri, save_petri_dish

year60_cells = load_snapshot_cells('data/rate_0.005000/snapshots/year60_snapshot.json.gz')

# Create last 2 control2 individuals
for i in range(28, 30):
    print(f'Creating individual {i}...')
    petri = create_pure_snapshot_petri(year60_cells, n_cells=5120, rate=0.005, seed=42+300+i)
    petri.metadata = {'individual_id': i, 'individual_type': 'control2'}
    save_petri_dish(petri, f'data/rate_0.005000/individuals/control2/individual_{i:02d}.json.gz')
print('Done!')
"
# Success! Created individuals 28-29
```

### Step 6: Verify All Individuals Created
```bash
ls data/rate_0.005000/individuals/mutant/ | wc -l     # 30 ✓
ls data/rate_0.005000/individuals/control1/ | wc -l   # 30 ✓  
ls data/rate_0.005000/individuals/control2/ | wc -l   # 30 ✓
```

### Step 7: Run Analysis (Also Timed Out)
```python
# Created complete_analysis.py script but it timed out too
python complete_analysis.py  # TIMEOUT

# So just did a quick verification instead:
python -c "
import json, gzip, numpy as np

with gzip.open('data/rate_0.005000/individuals/mutant/individual_00.json.gz', 'rt') as f:
    data = json.load(f)
    
print(f'Cells in file: {len(data[\"cells\"])}')
jsds = [c['jsd'] for c in data['cells']]
print(f'Mean JSD: {np.mean(jsds):.6f}')
"
# Output:
# Cells in file: 5120
# Mean JSD: 0.575080
```

## Why The Timeouts?

1. **10,000 cells is a LOT**: Each snapshot has 10,000 Cell objects to process
2. **90 individuals total**: 30 mutant + 30 control1 + 30 control2
3. **5120 cells per individual**: After mixing, each file has 5120 cells
4. **Total cells processed**: ~460,800 cells (90 individuals × 5120 cells)

## Clean Reproduction Steps (If Starting Fresh)

### Option A: Full Run (May Take Hours)
```bash
cd step23-prime

# Clean any previous attempts
rm -rf data/rate_0.005000/

# Run with extended timeout or in background
nohup python run_pipeline.py \
    --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --growth-years 10 \
    --mix-ratio 80 \
    --seed 42 \
    --bins 200 > pipeline.log 2>&1 &

# Monitor progress
tail -f pipeline.log
```

### Option B: Run in Stages (Recommended)

1. **First run to create snapshots and individuals**:
```bash
# This will create year 50/60 snapshots and start processing
python run_pipeline.py --rate 0.005 --simulation ../step1/data/*.json.gz
# Let it run until timeout, it will create snapshots and some individuals
```

2. **Check progress**:
```bash
ls -la data/rate_0.005000/snapshots/  # Should have 2 files
ls data/rate_0.005000/individuals/*/ | wc -l  # Check each directory
```

3. **If control2 missing, create manually**:
```bash
# Use the create_control2.py script in batches
# Or run in background with nohup
```

4. **Run analysis separately**:
```bash
# Use complete_analysis.py or do selective analysis on subset
```

### Option C: Test with Smaller Dataset First
```bash
# Test with fewer individuals to verify it works
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz \
    --n-quantiles 2 \
    --cells-per-quantile 1 \
    --growth-years 2 \
    --mix-ratio 80 \
    --seed 42
# This creates only 2 individuals per group with 20 cells each (much faster)
```

## Key Lessons Learned

1. **The pipeline is stateful**: It checks what exists and skips completed stages
2. **Snapshots are cached**: Year 50/60 extracts are reused (saves time)
3. **Control2 is the bottleneck**: Creating 30 × 5120 pure year-60 individuals is slow
4. **Analysis is memory-intensive**: Loading 90 individuals with 5120 cells each uses lots of RAM

## What Was NOT Cleaned/Restarted

**Important**: I did NOT clean and restart. The pipeline's skip logic allowed it to:
- Use cached snapshots
- Keep already-created mutant/control1 individuals  
- Only create missing control2 individuals
- This is actually a FEATURE - the pipeline can resume from interruptions!

## Verification Commands
```bash
# Quick check that reproduction worked
python -c "
import json, gzip
# Check one file from each group
for group in ['mutant', 'control1', 'control2']:
    with gzip.open(f'data/rate_0.005000/individuals/{group}/individual_00.json.gz', 'rt') as f:
        data = json.load(f)
        print(f'{group}: {len(data[\"cells\"])} cells')
"
# Should show 5120 cells for each group
```