# Git Commit Plan for step23-prime

## Commit 1: Core Pipeline Implementation
**Message**: `feat(step23-prime): Implement refactored pipeline with PetriDish/Cell objects`

**Files**:
```bash
git add step23-prime/__init__.py
git add step23-prime/run_pipeline.py
git add step23-prime/pipeline_utils.py
git add step23-prime/pipeline_analysis.py
```

**Description**: Core implementation of step23-prime using cleaner OOP design with PetriDish and Cell objects instead of dictionaries.

---

## Commit 2: Reproducibility Fixes
**Message**: `fix(step23-prime): Add proper random seeding for reproducibility`

**Files**:
```bash
# The fixes are already in run_pipeline.py and pipeline_utils.py
# Just document the analysis and fixes
git add step23-prime/RANDOM_SEED_ANALYSIS.md
git add step23-prime/REPRODUCIBILITY_FIXES.md
```

**Description**: Fixed missing global seeding and numpy random seeding to ensure reproducible results.

---

## Commit 3: Testing and Validation Suite
**Message**: `test(step23-prime): Add comprehensive testing and comparison tools`

**Files**:
```bash
# Comparison scripts
git add step23-prime/compare_snapshots.py
git add step23-prime/compare_snapshots_detailed.py
git add step23-prime/compare_individuals.py
git add step23-prime/compare_individuals_detailed.py
git add step23-prime/compare_statistical_tests.py
git add step23-prime/compare_two_runs.py
git add step23-prime/generate_comparison_report.py
git add step23-prime/quick_individual_stats.py
git add step23-prime/run_all_comparisons.sh

# Reproducibility tests
git add step23-prime/test_reproducibility.py
git add step23-prime/test_reproducibility_robust.py
git add step23-prime/verify_determinism.py
git add step23-prime/test_full_reproducibility.py
git add step23-prime/test_seed_fixes.sh
```

**Description**: Comprehensive test suite to validate statistical equivalence with step23 and ensure reproducibility.

---

## Commit 4: Utility Scripts
**Message**: `feat(step23-prime): Add utility scripts for pipeline management`

**Files**:
```bash
git add step23-prime/clean_individuals.sh
git add step23-prime/complete_analysis.py
git add step23-prime/create_control2.py
```

**Description**: Helper scripts for cleaning, completing analysis, and creating control groups.

---

## Commit 5: Documentation
**Message**: `docs(step23-prime): Add comprehensive documentation and guides`

**Files**:
```bash
git add step23-prime/REPRODUCTION_PLAN.md
git add step23-prime/REPRODUCTION_RESULTS.md
git add step23-prime/EXACT_REPRODUCTION_STEPS.md
git add step23-prime/COMPARISON_PLAN.md
git add step23-prime/CLEANING_GUIDE.md
git add step23-prime/GIT_COMMIT_PLAN.md
```

**Description**: Documentation covering reproduction plans, results, cleaning guides, and comparison strategies.

---

## Commit 6: Debug and Fix Scripts
**Message**: `debug(step23-prime): Add debugging tools for troubleshooting`

**Files**:
```bash
git add step23-prime/debug_snapshots.py
git add step23-prime/check_snapshot_issue.py
git add step23-prime/fix_random_seeding.py
```

**Description**: Debugging tools created during development to identify and fix issues.

---

## Commands to Execute (in order):

```bash
# Commit 1: Core implementation
git add step23-prime/__init__.py step23-prime/run_pipeline.py step23-prime/pipeline_utils.py step23-prime/pipeline_analysis.py
git commit -m "feat(step23-prime): Implement refactored pipeline with PetriDish/Cell objects"

# Commit 2: Reproducibility fixes (documentation)
git add step23-prime/RANDOM_SEED_ANALYSIS.md step23-prime/REPRODUCIBILITY_FIXES.md
git commit -m "fix(step23-prime): Add proper random seeding for reproducibility"

# Commit 3: Testing suite
git add step23-prime/compare_*.py step23-prime/test_*.py step23-prime/*.sh step23-prime/quick_individual_stats.py step23-prime/generate_comparison_report.py step23-prime/verify_determinism.py
git commit -m "test(step23-prime): Add comprehensive testing and comparison tools"

# Commit 4: Utilities
git add step23-prime/complete_analysis.py step23-prime/create_control2.py
git commit -m "feat(step23-prime): Add utility scripts for pipeline management"

# Commit 5: Documentation
git add step23-prime/*.md
git commit -m "docs(step23-prime): Add comprehensive documentation and guides"

# Commit 6: Debug tools
git add step23-prime/debug_snapshots.py step23-prime/check_snapshot_issue.py step23-prime/fix_random_seeding.py
git commit -m "debug(step23-prime): Add debugging tools for troubleshooting"
```

---

## Alternative: Single Comprehensive Commit

If you prefer one commit:

```bash
git add step23-prime/
git commit -m "feat(step23-prime): Refactor step23 with PetriDish/Cell objects and reproducibility fixes

- Implement cleaner pipeline using OOP design with PetriDish and Cell classes
- Fix reproducibility issues with proper random seeding
- Add comprehensive testing and comparison suite
- Include utility scripts for pipeline management
- Document reproduction process and results
- Validate statistical equivalence with original step23"
```

---

## Notes:
- All data files (*.json.gz) are excluded
- __pycache__ directories are excluded
- The changes maintain scientific validity while improving code quality
- Reproducibility is now guaranteed with proper seeding