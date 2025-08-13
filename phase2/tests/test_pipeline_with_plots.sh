#!/bin/bash
# Test script to demonstrate running phase2 pipeline with individual plots

echo "=========================================="
echo "Testing Phase2 Pipeline with History Plots"
echo "=========================================="

# Change to phase2 directory
cd ../

# Check if a phase1 simulation exists for testing
TEST_SIM="../phase1/data/rate_0.00500/grow4-sites100-years20-seed42*/simulation.json.gz"

if ! ls $TEST_SIM 1> /dev/null 2>&1; then
    echo ""
    echo "No test simulation found. Creating a small one..."
    echo ""
    
    # Create a small phase1 simulation for testing
    cd ../phase1
    python run_simulation.py --rate 0.005 --years 20 --growth-phase 4 --sites 100 --seed 42
    cd ../phase2
fi

# Find the simulation file
TEST_SIM=$(ls ../phase1/data/rate_0.00500/grow4-sites100-years20-seed42*/simulation.json.gz 2>/dev/null | head -1)

if [ -z "$TEST_SIM" ]; then
    echo "Error: Could not find or create test simulation"
    exit 1
fi

echo ""
echo "Using test simulation: $TEST_SIM"
echo ""
echo "Running phase2 pipeline with --plot-individuals flag..."
echo ""

# Run the pipeline with history tracking and plotting
python run_pipeline.py \
    --rate 0.005 \
    --simulation "$TEST_SIM" \
    --first-snapshot 10 \
    --second-snapshot 15 \
    --n-quantiles 3 \
    --cells-per-quantile 2 \
    --individual-growth-phase 3 \
    --mix-ratio 70 \
    --normalize-size \
    --uniform-mixing \
    --plot-individuals \
    --seed 42

# Check if plots were generated
OUTPUT_DIR=$(ls -td data/rate_0.00500-grow4-sites100-years20/snap10-quant3x2-grow3-mix70-seed42-* 2>/dev/null | head -1)

if [ -z "$OUTPUT_DIR" ]; then
    echo ""
    echo "Error: Pipeline output directory not found"
    exit 1
fi

echo ""
echo "Pipeline completed. Checking for generated plots..."
echo ""

# Check for individual plots
PLOT_DIR="$OUTPUT_DIR/individual_plots"

if [ -d "$PLOT_DIR" ]; then
    echo "✓ Individual plots directory created: $PLOT_DIR"
    
    # Count plots
    JSD_COUNT=$(ls $PLOT_DIR/*_jsd.png 2>/dev/null | wc -l)
    METH_COUNT=$(ls $PLOT_DIR/*_methylation.png 2>/dev/null | wc -l)
    COMBINED_COUNT=$(ls $PLOT_DIR/*_combined.png 2>/dev/null | wc -l)
    
    echo "  - JSD plots: $JSD_COUNT"
    echo "  - Methylation plots: $METH_COUNT"
    echo "  - Combined plots: $COMBINED_COUNT"
    
    if [ $JSD_COUNT -gt 0 ]; then
        echo ""
        echo "✓ SUCCESS: Individual growth trajectory plots generated!"
    else
        echo ""
        echo "⚠ WARNING: No plots found, but directory exists"
    fi
else
    echo "✗ Individual plots directory not found"
    echo "  Expected: $PLOT_DIR"
fi

echo ""
echo "Test complete."
echo ""
echo "To view the plots, navigate to:"
echo "  $PLOT_DIR"
echo ""