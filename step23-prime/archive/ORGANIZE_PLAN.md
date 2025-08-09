# Step23-Prime Organization Plan

## Current State: 23 Python files + multiple docs

## Proposed Structure:

```
step23-prime/
├── run_pipeline.py              # Main entry point
├── pipeline_utils.py            # Core utilities
├── pipeline_analysis.py         # Analysis functions
├── path_utils.py                # Path parsing/generation
├── __init__.py                  # Package init
├── README.md                    # Main documentation
│
├── tests/                       # All test files
│   ├── test_reproducibility.py
│   ├── test_reproducibility_robust.py
│   ├── test_full_reproducibility.py
│   ├── test_dynamic_mix_year.py
│   └── verify_determinism.py
│
├── tools/                       # Comparison and debug tools
│   ├── compare_individuals.py
│   ├── compare_snapshots.py
│   ├── compare_statistical_tests.py
│   ├── compare_two_runs.py
│   ├── generate_comparison_report.py
│   └── quick_individual_stats.py
│
├── scripts/                     # Utility scripts
│   ├── complete_analysis.py
│   ├── create_control2.py
│   └── clean_individuals.sh
│
└── archive/                     # Development artifacts (can be deleted)
    ├── debug_snapshots.py
    ├── check_snapshot_issue.py
    ├── fix_random_seeding.py
    ├── compare_individuals_detailed.py
    ├── compare_snapshots_detailed.py
    ├── test_seed_fixes.sh
    ├── run_all_comparisons.sh
    ├── COMPARISON_REPORT.txt
    ├── add_global_seeding.patch
    └── [various .md planning docs]
```

## Files to Keep in Root:
- Core pipeline files (run_pipeline.py, pipeline_*.py, path_utils.py)
- Main documentation (README.md)
- Package files (__init__.py)

## Files to Move:
- **tests/**: All test_*.py and verify_*.py files
- **tools/**: Comparison and analysis tools
- **scripts/**: Utility scripts for specific tasks
- **archive/**: Debug files and development artifacts

## Files to Consider Deleting:
- Detailed comparison variants (keep only main versions)
- Debug files that were for specific issues (now fixed)
- Old planning documents (already in git history)
- Patch files (fixes already applied)

## Commands to Execute:

```bash
# Create directories
mkdir -p step23-prime/{tests,tools,scripts,archive}

# Move test files
mv step23-prime/test*.py step23-prime/tests/
mv step23-prime/verify*.py step23-prime/tests/

# Move comparison tools
mv step23-prime/compare*.py step23-prime/tools/
mv step23-prime/generate_comparison_report.py step23-prime/tools/
mv step23-prime/quick_individual_stats.py step23-prime/tools/

# Move utility scripts
mv step23-prime/complete_analysis.py step23-prime/scripts/
mv step23-prime/create_control2.py step23-prime/scripts/
mv step23-prime/clean_individuals.sh step23-prime/scripts/

# Archive development files
mv step23-prime/{debug,check,fix}*.py step23-prime/archive/
mv step23-prime/*.patch step23-prime/archive/
mv step23-prime/COMPARISON_REPORT.txt step23-prime/archive/
mv step23-prime/*PLAN*.md step23-prime/archive/
mv step23-prime/*STEPS*.md step23-prime/archive/
mv step23-prime/*RESULTS*.md step23-prime/archive/
mv step23-prime/*ANALYSIS*.md step23-prime/archive/
mv step23-prime/*FIXES*.md step23-prime/archive/
```

## Alternative: Aggressive Cleanup

If you want to be more aggressive and only keep essentials:

```bash
# Keep only core files and main comparison tool
mkdir -p step23-prime/tools

# Move the main comparison tool
mv step23-prime/generate_comparison_report.py step23-prime/tools/

# Delete everything else that's not core
rm -f step23-prime/test*.py
rm -f step23-prime/compare*.py
rm -f step23-prime/debug*.py
rm -f step23-prime/check*.py
rm -f step23-prime/fix*.py
rm -f step23-prime/verify*.py
rm -f step23-prime/quick*.py
rm -f step23-prime/*.md
rm -f step23-prime/*.txt
rm -f step23-prime/*.patch
rm -f step23-prime/*.sh
```

## Benefits of Organization:
1. Clear separation of concerns
2. Easy to find files by purpose
3. Can easily exclude tests/tools from production
4. Development artifacts archived but available
5. Clean root directory with only essentials