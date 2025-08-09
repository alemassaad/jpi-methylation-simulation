# Git Commit Plan - Final Session Changes

## Summary of All Changes

### Major Features Implemented:
1. **Dynamic snapshot years**: Replaced hardcoded year 50/60 with configurable snapshot years
2. **Hierarchical directory structure**: Both step1-prime and step23-prime now use organized hierarchical paths
3. **Improved filename formats**: Removed variable `m` parameter, fixed rate precision, added hashing
4. **Pipeline enhancements**: Better parameter handling and path generation

### Files Changed:
- **step1-prime/cell.py**: New hierarchical save structure, removed `m` from filename
- **step23-prime/run_pipeline.py**: Dynamic snapshot years, new directory structure
- **step23-prime/path_utils.py**: NEW - Utilities for parsing and generating paths
- **step23-prime/test_dynamic_mix_year.py**: NEW - Tests for dynamic snapshot feature
- **CLAUDE.md**: Documentation updates for step23-prime
- **README.md**: Documentation updates for step23-prime
- **CLEANUP_PLAN.md**: NEW - Cleanup recommendations (can be excluded from commits)

---

## Commit Structure (5 commits)

### Commit 1: Dynamic Snapshot Years Feature
**Message**: 
```
feat(step23-prime): Add configurable snapshot years for flexible pipeline execution

- Add --snapshot-year parameter to control first snapshot (default: 50)
- Make second snapshot always dynamic (snapshot-year + growth-years)
- Remove hardcoded year 60 references throughout pipeline
- Update all stage titles and messages to use calculated years
- Add validation for snapshot year parameters
- Ensure age-aligned mixing (grown cells mix with same-age snapshot)

This change makes the pipeline more flexible and scientifically sound by ensuring
all mixed populations are age-matched.
```

**Files**:
- `step23-prime/run_pipeline.py` (partial - only snapshot-related changes)
- `step23-prime/test_dynamic_mix_year.py`

---

### Commit 2: Hierarchical Directory Structure for step1-prime
**Message**:
```
refactor(step1-prime): Implement hierarchical directory structure with improved naming

- Change directory structure to: data/rate_X.XXXXX/params-hash/simulation.json.gz
- Remove variable 'm' parameter from filenames (population varies with culling)
- Change rate precision from 6 to 5 decimal places (0.005000 → 0.00500)
- Add 4-character MD5 hash for uniqueness
- Use fixed filename 'simulation.json.gz' instead of long names
- Maintain backward compatibility for existing code

Example: data/rate_0.00500/grow13-sites1000-years100-seed42-a3f2/simulation.json.gz
```

**Files**:
- `step1-prime/cell.py`

---

### Commit 3: Hierarchical Directory Structure for step23-prime
**Message**:
```
refactor(step23-prime): Implement hierarchical directory structure with parameter tracking

- Add path_utils.py for parsing and generating hierarchical paths
- Parse source simulation parameters from both old and new formats
- Generate descriptive directory names with all parameters
- Structure: data/rate_X-sourceParams/pipelineParams-hash/
- Use hyphens as separators for readability
- Support both old flat and new hierarchical simulation formats

Example: data/rate_0.00500-grow13-sites1000-years100/snap50-quant10x3-grow10-mix80-seed42-7b3f/
```

**Files**:
- `step23-prime/run_pipeline.py` (partial - only structure-related changes)
- `step23-prime/path_utils.py`

---

### Commit 4: Documentation Updates
**Message**:
```
docs: Update documentation for step23-prime and new directory structures

- Document step23-prime in CLAUDE.md with features and usage
- Update README.md with step23-prime pipeline information
- Add examples of new hierarchical directory structures
- Document the dynamic snapshot year feature
- Update command examples with new paths
```

**Files**:
- `CLAUDE.md`
- `README.md`

---

### Commit 5: Development Artifacts (Optional)
**Message**:
```
chore: Add development and testing artifacts from step23-prime implementation

- Add comparison report from validation testing
- Add cleanup plan for organizing directories
- Add patch file documenting seeding fixes
- These files document the development process and validation

Note: These files can be excluded if only production code is desired
```

**Files**:
- `CLEANUP_PLAN.md`
- `step23-prime/COMPARISON_REPORT.txt`
- `step23-prime/add_global_seeding.patch`

---

## Commands to Execute (DO NOT RUN YET)

```bash
# Commit 1: Dynamic Snapshot Years
git add step23-prime/test_dynamic_mix_year.py
git add -p step23-prime/run_pipeline.py  # Select only snapshot-related changes
git commit -m "feat(step23-prime): Add configurable snapshot years for flexible pipeline execution"

# Commit 2: step1-prime Structure
git add step1-prime/cell.py
git commit -m "refactor(step1-prime): Implement hierarchical directory structure with improved naming"

# Commit 3: step23-prime Structure  
git add step23-prime/path_utils.py
git add step23-prime/run_pipeline.py  # Add remaining changes
git commit -m "refactor(step23-prime): Implement hierarchical directory structure with parameter tracking"

# Commit 4: Documentation
git add CLAUDE.md README.md
git commit -m "docs: Update documentation for step23-prime and new directory structures"

# Commit 5: Development Artifacts (Optional)
git add CLEANUP_PLAN.md step23-prime/COMPARISON_REPORT.txt step23-prime/add_global_seeding.patch
git commit -m "chore: Add development and testing artifacts from step23-prime implementation"
```

---

## Notes

- The commits are ordered logically: features → refactoring → documentation → artifacts
- Each commit is atomic and focuses on a specific aspect
- Commit messages follow conventional commit format
- The development artifacts commit is optional and can be skipped
- Consider squashing commits 2 & 3 if you prefer fewer commits

## Alternative: 3 Commits Instead

If you prefer fewer commits:

1. **feat: Add dynamic snapshot years and hierarchical directory structures**
   - Combine commits 1, 2, and 3
   - All code changes in one commit

2. **docs: Update documentation for new features**
   - Commit 4 as-is

3. **chore: Add development artifacts**
   - Commit 5 as-is (or skip)