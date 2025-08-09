# Reproducibility Test

This test suite ensures that the methylation simulation produces deterministic, reproducible results when using fixed random seeds.

## Purpose

- Verify that code changes don't alter simulation behavior
- Ensure results are reproducible across different runs
- Detect unintended changes to the simulation logic

## Files

- `test_reproducibility.py` - The test script
- `test_reproducibility_expected.json` - Expected results (generated with current code)

## Usage

### Generate/Update Expected Results
```bash
cd tests
python test_reproducibility.py
```
This runs all tests and saves the results to `test_reproducibility_expected.json`.

### Verify Against Expected Results
```bash
cd tests
python test_reproducibility.py --check
```
This runs the tests and compares results against the saved expected values.

Or from the project root:
```bash
python tests/test_reproducibility.py --check
```

## Tests Included

1. **Single Cell Aging** (seed=42)
   - Ages one cell for 10 years
   - Tracks methylation proportion, distribution, and JSD at each year

2. **Full Simulation** (seed=12345)
   - 100 cells for 20 years
   - Verifies aggregate statistics and specific cell states

3. **Edge Cases** (seed=999)
   - High methylation rate (20%) to test fully methylated cells
   - Verifies early exit optimization

4. **Different Parameters** (various seeds)
   - Tests with different n, rate, and years combinations
   - Ensures robustness across parameter space

## How It Works

The tests use fixed random seeds to make the stochastic simulation deterministic. Any change to the simulation logic (intentional or accidental) will cause the tests to fail, alerting you to investigate.

## When to Update Expected Results

Update the expected results file **only** when:
1. You intentionally change the simulation logic
2. You've verified the changes are correct
3. You want the new behavior to be the reference

To update:
```bash
cd tests
python test_reproducibility.py  # Regenerate expected results
git add test_reproducibility_expected.json
git commit -m "Update expected results after [describe change]"
```