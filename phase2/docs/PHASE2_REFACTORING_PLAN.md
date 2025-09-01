# Phase 2 Refactoring Plan

## Current Issues
1. **Root clutter**: 10 Python files mixed in root directory
2. **Test files scattered**: Test files in root (`test_*.py`) AND in `tests/` directory
3. **No clear module structure**: Core pipeline code mixed with utilities and analysis
4. **Inconsistent naming**: Some files use underscores, some don't
5. **Documentation scattered**: Multiple .md files in root
6. **Temporary files**: test_run.log, test_output/, test_timeline_plots/ in root

## Proposed New Structure

```
phase2/
├── README.md                     # Main documentation (keep in root)
├── __init__.py                   # Package init
├── run_pipeline.py              # Main entry point (keep in root for easy access)
│
├── core/                        # Core pipeline modules
│   ├── __init__.py
│   ├── pipeline_utils.py       # Pipeline utilities
│   ├── pipeline_analysis.py    # Analysis functions
│   └── path_utils.py           # Path generation utilities
│
├── visualization/               # All plotting/visualization code
│   ├── __init__.py
│   ├── plot_individuals.py     # Individual trajectory plots
│   └── analyze_individual_sizes.py  # Size distribution analysis
│
├── scripts/                     # Standalone scripts (already exists)
│   ├── clean_individuals.sh
│   ├── complete_analysis.py
│   └── create_control2.py
│
├── tools/                       # Analysis tools (already exists)
│   └── compare_two_runs.py
│
├── configs/                     # Configuration files (already exists)
│   ├── config_default.yaml     # Move from root
│   ├── debug.yaml
│   ├── full_analysis.yaml
│   ├── quick_test.yaml
│   └── uniform_mixing.yaml
│
├── tests/                       # ALL test files
│   ├── __init__.py
│   ├── config/                 # Config-related tests
│   │   ├── test_config_phase2.py
│   │   ├── test_config_simple.py
│   │   └── test_pipeline_with_config.py
│   ├── integration/            # Integration tests
│   │   └── test_final_integration.py
│   └── unit/                   # Unit tests
│       └── [other test files]
│
├── docs/                        # Documentation
│   ├── NORMALIZATION_USAGE.md
│   ├── JSON_CONSOLIDATION_PLAN.md
│   └── PHASE2_REFACTORING_PLAN.md
│
├── data/                        # Output data (already exists, gitignored)
└── tmp/                         # Temporary files (gitignored)
    ├── test_output/
    ├── test_timeline_plots/
    └── *.log
```

## Migration Steps

### Step 1: Create New Directory Structure
```bash
mkdir -p core visualization docs tmp
mkdir -p tests/config tests/integration tests/unit
```

### Step 2: Move Files
```bash
# Core modules
mv pipeline_utils.py pipeline_analysis.py path_utils.py core/

# Visualization modules  
mv plot_individuals.py analyze_individual_sizes.py visualization/

# Config files
mv config_default.yaml configs/

# Test files
mv test_config*.py test_pipeline_with_config.py tests/config/
mv tests/test_final_integration.py tests/integration/
mv tests/test_*.py tests/unit/

# Documentation
mv NORMALIZATION_USAGE.md JSON_CONSOLIDATION_PLAN.md docs/

# Temporary files
mv test_output test_timeline_plots tmp/
mv *.log tmp/
```

### Step 3: Update Imports
All files will need import updates:

**In run_pipeline.py:**
```python
from core.pipeline_utils import *
from core.pipeline_analysis import *
from core.path_utils import *
from visualization.plot_individuals import plot_all_individuals
```

**In core modules:**
```python
# Add relative imports between core modules
from .pipeline_utils import ...
from .path_utils import ...
```

### Step 4: Update sys.path additions
Remove hardcoded sys.path manipulations and use proper package imports

### Step 5: Create __init__.py files
Add proper __init__.py files to make subdirectories proper Python packages

### Step 6: Update .gitignore
```
# Add to .gitignore
phase2/tmp/
phase2/data/
*.log
__pycache__/
.DS_Store
```

## Benefits
1. **Clear separation of concerns**: Core logic, visualization, testing, tools
2. **Easier navigation**: Know exactly where to find specific functionality
3. **Better imports**: Use proper package imports instead of sys.path hacks
4. **Cleaner root**: Only README and main entry point in root
5. **Professional structure**: Follows Python package best practices

## Testing Plan
After refactoring:
1. Run main pipeline with various configurations
2. Run all tests from new locations
3. Verify all imports work correctly
4. Check that relative paths still work
5. Ensure backward compatibility for command-line usage

## Rollback Plan
If issues arise:
1. Git stash changes
2. Revert to previous commit
3. Apply refactoring incrementally