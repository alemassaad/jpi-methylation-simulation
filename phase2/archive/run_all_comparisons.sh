#!/bin/bash
# Master script to run all comparisons between step23 and step23-prime

echo "=========================================="
echo "COMPLETE COMPARISON: step23 vs step23-prime"
echo "=========================================="
echo ""

# Check if both directories exist
if [ ! -d "../step23/data/rate_0.005000" ]; then
    echo "❌ Error: step23 results not found at ../step23/data/rate_0.005000"
    echo "Please ensure step23 pipeline has completed"
    exit 1
fi

if [ ! -d "data/rate_0.005000" ]; then
    echo "❌ Error: step23-prime results not found at data/rate_0.005000"
    echo "Please ensure step23-prime pipeline has completed"
    exit 1
fi

echo "✓ Both directories found"
echo ""

# 1. Compare snapshots
echo "=========================================="
echo "1. SNAPSHOT COMPARISON (Cell-by-cell)"
echo "=========================================="
echo "Comparing year 50 and 60 snapshots..."
python compare_snapshots_detailed.py
echo ""

# 2. Quick individual statistics
echo "=========================================="
echo "2. QUICK INDIVIDUAL STATISTICS"
echo "=========================================="
echo "Getting basic statistics for all individuals..."
python quick_individual_stats.py
echo ""

# 3. Detailed individual comparison
echo "=========================================="
echo "3. DETAILED INDIVIDUAL COMPARISON"
echo "=========================================="
echo "This may take a minute..."
python compare_individuals_detailed.py
echo ""

# 4. Statistical test comparison
echo "=========================================="
echo "4. STATISTICAL TEST RESULTS"
echo "=========================================="
echo "Comparing p-values and test results..."
python compare_statistical_tests.py
echo ""

# 5. Generate final report
echo "=========================================="
echo "5. GENERATING FINAL REPORT"
echo "=========================================="
python generate_comparison_report.py

echo ""
echo "=========================================="
echo "COMPARISON COMPLETE!"
echo "Check COMPARISON_REPORT.txt for full results"
echo "=========================================="