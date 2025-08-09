# Cleanup Plan for jpi-methylation-simulation

## Items to Clean

### 1. Temporary Test Data Directories
```bash
# In step23/ - test data from development
rm -rf step23/data_test      # 3.3M of test data
rm -rf step23/data_test10    # 4K empty directory
```

### 2. Python Cache Files
```bash
# Remove all __pycache__ and .pyc files
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete
```

### 3. macOS System Files
```bash
# Remove .DS_Store files
find . -name ".DS_Store" -delete
```

### 4. step23-prime Development/Debug Scripts
These scripts were created during development and debugging. Consider archiving or removing:

**Debug/Fix Scripts** (used to identify and fix issues, may not be needed anymore):
- `step23-prime/debug_snapshots.py` - Used to debug snapshot differences
- `step23-prime/check_snapshot_issue.py` - Identified file hash vs content issue
- `step23-prime/fix_random_seeding.py` - Documents seeding fixes (already applied)

**Redundant Test Scripts** (functionality covered by test_reproducibility_robust.py):
- `step23-prime/test_reproducibility.py` - Basic version
- `step23-prime/test_full_reproducibility.py` - Another variant
- `step23-prime/test_seed_fixes.sh` - Shell script for testing

**Detailed Comparison Scripts** (keep generate_comparison_report.py, remove verbose variants):
- `step23-prime/compare_snapshots_detailed.py` - Verbose version
- `step23-prime/compare_individuals_detailed.py` - Verbose version

### 5. Development Documentation
Consider archiving these planning/analysis documents:
- `step23-prime/RANDOM_SEED_ANALYSIS.md` - Analysis done, fixes applied
- `step23-prime/REPRODUCIBILITY_FIXES.md` - Fixes already implemented
- `step23-prime/REPRODUCTION_PLAN.md` - Plan executed
- `step23-prime/EXACT_REPRODUCTION_STEPS.md` - Steps completed
- `step23-prime/GIT_COMMIT_PLAN.md` - Commits done
- `step23-prime/COMPARISON_REPORT.txt` - Old report output

## Items to Keep

### Essential Scripts
- **Core Pipeline**: `run_pipeline.py`, `pipeline_*.py`
- **Main Test**: `test_reproducibility_robust.py`
- **Main Comparison**: `generate_comparison_report.py`
- **Basic Comparisons**: `compare_snapshots.py`, `compare_individuals.py`, `compare_statistical_tests.py`
- **Utilities**: `clean_individuals.sh`, `complete_analysis.py`, `create_control2.py`

### Important Documentation
- `step23-prime/COMPARISON_PLAN.md` - Useful reference
- `step23-prime/REPRODUCTION_RESULTS.md` - Results documentation
- `step23-prime/CLEANING_GUIDE.md` - User guide

## Cleanup Commands

### Quick Cleanup (safe)
```bash
# Remove test data and cache files
rm -rf step23/data_test step23/data_test10
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete
find . -name ".DS_Store" -delete
```

### Full Cleanup (archives debug/development files)
```bash
# Create archive directory
mkdir -p step23-prime/archive_dev

# Move debug and development files
mv step23-prime/debug_snapshots.py step23-prime/archive_dev/
mv step23-prime/check_snapshot_issue.py step23-prime/archive_dev/
mv step23-prime/fix_random_seeding.py step23-prime/archive_dev/
mv step23-prime/test_reproducibility.py step23-prime/archive_dev/
mv step23-prime/test_full_reproducibility.py step23-prime/archive_dev/
mv step23-prime/test_seed_fixes.sh step23-prime/archive_dev/
mv step23-prime/compare_*_detailed.py step23-prime/archive_dev/
mv step23-prime/RANDOM_SEED_ANALYSIS.md step23-prime/archive_dev/
mv step23-prime/REPRODUCIBILITY_FIXES.md step23-prime/archive_dev/
mv step23-prime/REPRODUCTION_PLAN.md step23-prime/archive_dev/
mv step23-prime/EXACT_REPRODUCTION_STEPS.md step23-prime/archive_dev/
mv step23-prime/GIT_COMMIT_PLAN.md step23-prime/archive_dev/
mv step23-prime/COMPARISON_REPORT.txt step23-prime/archive_dev/ 2>/dev/null || true

# Remove test data and cache
rm -rf step23/data_test step23/data_test10
find . -type d -name "__pycache__" -exec rm -rf {} + 2>/dev/null
find . -name "*.pyc" -delete
find . -name ".DS_Store" -delete
```

## Summary

- **Total space to free**: ~3.3MB (mostly test data)
- **Files to archive**: 14 development/debug scripts and docs
- **Files to keep**: 12 essential scripts and 3 user guides
- **Recommendation**: Run the "Full Cleanup" to archive development files while keeping essential functionality