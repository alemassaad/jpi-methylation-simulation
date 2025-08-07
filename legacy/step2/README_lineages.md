# Step 2: Cell Division with Lineage Tracking

This directory now contains updated scripts that track cell lineages separately for more biologically meaningful analysis.

## New Approach: Separate Lineage Tracking

### Overview

Instead of mixing all 30,720 cells together, we now:
1. Take 30 cells from year 50 (3 from each JSD decile)
2. Age each cell SEPARATELY with division for 10 years
3. Create 30 separate files, each containing 1,024 cells from one lineage
4. Each lineage represents cells descended from a single original cell

### Scripts

#### `sample_divide_age_lineages.py`
The main script that creates 30 separate lineage files:
```bash
cd step2
python sample_divide_age_lineages.py
```

This will:
- Load the year 50 snapshot
- Sample 3 cells from each JSD decile (30 total)
- Age each cell separately for 10 years with division
- Save 30 files to `lineages/` directory
- Each file contains 1,024 cells (2^10 from division)

Output files are named: `lineage_XX_decile_YY.json.gz` where:
- XX = lineage ID (00-29)
- YY = source decile (01-10)

#### `test_sample_divide_age_lineages.py`
Test version with smaller parameters for quick verification:
```bash
python test_sample_divide_age_lineages.py
```

Test parameters:
- 1 cell per decile (10 total instead of 30)
- 3 years of aging (8 cells per lineage instead of 1,024)
- Output to `test_lineages/` directory
- Takes only seconds to run

### File Format

Each lineage file contains:
```json
{
  "metadata": {
    "lineage_id": 0,
    "source_decile": 1,
    "within_decile_index": 0,
    "source_year": 50,
    "final_year": 60,
    "generations": 10,
    "final_cell_count": 1024,
    "expected_cell_count": 1024,
    "source_cell": {
      // Complete data of the original cell at year 50
    }
  },
  "cells": [
    // 1,024 cells at year 60, all descendants of source_cell
  ]
}
```

### Why This Approach?

This allows us to:
1. Track which cells came from which original cell
2. Model each lineage as representing a different "individual"
3. Mix each lineage separately with control cells in step 3
4. Calculate individual-level statistics

### Original Scripts (Deprecated)

The original `sample_divide_age.py` that creates one combined file is still available but no longer recommended for use.

## Next Steps

In step 3, we will:
1. Load each of the 30 lineage files
2. Mix each lineage with original cells (e.g., 80% original, 20% lineage)
3. Create 30 "individuals" with mixed populations
4. Calculate statistics for each individual
5. Analyze the distribution of these 30 data points