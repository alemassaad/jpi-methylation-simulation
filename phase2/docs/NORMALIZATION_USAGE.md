# Size Normalization Feature - Usage Guide

## Overview
The size normalization feature ensures all individuals have exactly the same number of cells before mixing, addressing variation from homeostasis. It uses a **median - 0.5σ** threshold to determine the target size.

## Command Line Usage

### Basic Normalization
```bash
python run_pipeline.py \
  --rate 0.005 \
  --simulation ../phase1/data/.../simulation.json.gz \
  --first-snapshot 50 \
  --second-snapshot 60 \
  --individual-growth-phase 7 \
  --normalize-size \
  --seed 42
```

### Normalization + Uniform Mixing (Recommended)
```bash
python run_pipeline.py \
  --rate 0.005 \
  --simulation ../phase1/data/.../simulation.json.gz \
  --first-snapshot 50 \
  --second-snapshot 60 \
  --individual-growth-phase 7 \
  --normalize-size \
  --uniform-mixing \
  --seed 42
```

## How It Works

1. **After Growth Phase**: Individuals have varying sizes due to homeostasis stochasticity
2. **Calculate Threshold**: `threshold = median(all_sizes) - 0.5 * stdev(all_sizes)`
3. **Normalize**:
   - Individuals **below** threshold → EXCLUDED
   - Individuals **at** threshold → KEPT AS-IS
   - Individuals **above** threshold → RANDOMLY TRIMMED to threshold
4. **Result**: All remaining individuals have exactly the same cell count

## Directory Structure

Output directories include suffix indicators:
- No flags: `mix80-seed42-xxxx`
- Uniform only: `mix80u-seed42-xxxx`
- Normalized only: `mix80n-seed42-xxxx`  
- Both: `mix80un-seed42-xxxx`

## Example Output

```
=== APPLYING SIZE NORMALIZATION ===
Using median - 0.5σ threshold

=== NORMALIZATION STATISTICS ===
All individual sizes: min=180, max=342
Median: 308.0 cells
Std Dev: 32.8 cells
Threshold (median - 0.5σ): 292 cells

  Mutant 00: 312 → 292 cells (trimmed)
  Mutant 01: 285 cells - EXCLUDED (below threshold)
  Mutant 02: 292 cells - kept as is
  ...

Mutant summary:
  Kept: 24 (8 as-is, 16 trimmed)
  Excluded: 6

Control1 summary:
  Kept: 22 (7 as-is, 15 trimmed)
  Excluded: 8

Overall retention: 46/60 (76.7%)
```

## Testing

Run comprehensive tests:
```bash
cd phase2/tests
python test_normalization_comprehensive.py

# With cleanup of test data:
python test_normalization_comprehensive.py --cleanup
```

## Edge Cases Handled

1. **All individuals excluded**: Error message with suggestions
2. **Single individual**: Uses that individual's size as threshold
3. **100% mix ratio**: Replaces all cells with snapshot cells
4. **Threshold < 1**: Uses minimum observed size instead
5. **No individuals**: Returns empty lists

## Recommendations

1. **Use with Uniform Mixing**: **NEW - Uniform mixing now includes automatic normalization!** In latest version, `--uniform-mixing` automatically normalizes individuals first, then applies uniform pool. You can still combine with `--normalize-size` for the traditional median-0.5σ approach, but it's no longer required to avoid pool size warnings.

2. **Check Retention Rate**: If < 50% of individuals are retained, consider:
   - Adjusting growth phase parameters
   - Increasing initial sample size
   - Using default pipeline without normalization

3. **Reproducibility**: Always set `--seed` for reproducible results

## Statistical Notes

- **Median - 0.5σ** typically retains 60-70% of individuals
- More conservative than percentile methods
- Adapts to population variance
- Robust to outliers

## Quick Test Example

```bash
# Small test with existing simulation
python run_pipeline.py \
  --rate 0.005 \
  --simulation ../phase1/data/rate_0.00500/grow5-sites100-years20-seed47-93ce/simulation.json.gz \
  --first-snapshot 10 \
  --second-snapshot 15 \
  --individual-growth-phase 3 \
  --n-quantiles 4 \
  --cells-per-quantile 2 \
  --normalize-size \
  --uniform-mixing \
  --seed 42
```