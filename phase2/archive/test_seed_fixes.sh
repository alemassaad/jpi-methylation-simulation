#!/bin/bash
# Test that seed fixes produce reproducible results

echo "=========================================="
echo "TESTING SEED FIXES IN STEP23-PRIME"
echo "=========================================="
echo ""

# Step 1: Save current results (post-fix)
echo "Step 1: Saving current results (with seed fixes)..."
if [ -d "data/rate_0.005000" ]; then
    mv data/rate_0.005000 data/rate_0.005000_postfix
    echo "  ✓ Saved to data/rate_0.005000_postfix"
else
    echo "  ❌ No current results found"
    exit 1
fi

# Step 2: Clean and re-run pipeline with same parameters
echo ""
echo "Step 2: Re-running pipeline with seed=42..."
echo "  (This tests if results are reproducible with seed fixes)"
echo ""

# Run with exact same parameters as original
python run_pipeline.py \
    --rate 0.005 \
    --simulation ../step1/data/simulation_rate_0.005000_m10000_n1000_t100.json.gz \
    --n-quantiles 10 \
    --cells-per-quantile 3 \
    --growth-years 10 \
    --mix-ratio 80 \
    --seed 42 \
    --bins 200

echo ""
echo "Step 3: Comparing post-fix results..."
echo ""

# Compare the two runs
python compare_two_runs.py data/rate_0.005000_postfix data/rate_0.005000

echo ""
echo "=========================================="
echo "SEED FIX TEST COMPLETE"
echo "=========================================="
echo ""
echo "If the runs are identical, the seed fixes work!"
echo "If they differ, we still have reproducibility issues."