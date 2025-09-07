# Priority 3 Issues - Detailed Fixing Plan

## Overview
This document contains the complete fixing plan for all Priority 3 issues and remaining minor issues identified in the Phase 2 pipeline.

## Issues to Fix

### 1. PlotPaths Class Documentation Enhancement
**Location**: `phase2/core/plot_paths.py`
**Changes**:

#### A. Add comprehensive class docstring (after line 6, before class definition)
```python
"""Centralized path management for pipeline output organization.

This class provides a structured approach to organizing pipeline outputs
into metric-based subdirectories, ensuring consistent file placement and
easy navigation of results.

Directory Structure:
    results/
    ├── cell_metrics/          # Cell-level measurements
    │   ├── timeline/          # Time series plots
    │   ├── distributions/     # Histogram distributions
    │   ├── comparisons/       # Batch comparisons
    │   ├── individual_trajectories/  # Per-individual growth
    │   │   ├── jsd/          # Cell JSD trajectories
    │   │   │   ├── mutant/
    │   │   │   └── control1/
    │   │   └── methylation_proportion/  # Methylation trajectories
    │   │       ├── mutant/
    │   │       └── control1/
    │   └── analysis/         # JSON analysis results
    ├── gene_metrics/         # Gene-level measurements
    │   ├── timeline/         # Gene JSD over time
    │   ├── distributions/    # Gene JSD distributions
    │   ├── comparisons/      # Batch comparisons
    │   ├── individual_trajectories/  # Per-individual gene JSD
    │   │   ├── mutant/
    │   │   └── control1/
    │   ├── per_gene/         # Individual gene distributions
    │   └── analysis/         # JSON analysis results
    └── metadata/             # Pipeline configuration and metadata

Attributes:
    results_dir: Base results directory path
    cell_metrics_dir: Root for cell-level metrics
    gene_metrics_dir: Root for gene-level metrics
    metadata_dir: Root for metadata files
    All subdirectory attributes for organized structure

Example:
    >>> plot_paths = PlotPaths('/path/to/results')
    >>> plot_paths.create_all_directories()
    >>> jsd_plot = plot_paths.get_cell_jsd_timeline_path()
"""
```

#### B. Update create_all_directories method docstring (line 35)
```python
def create_all_directories(self):
    """Create all necessary subdirectories for organized output.
    
    Creates the complete directory structure for pipeline outputs,
    ensuring all paths are available before file generation begins.
    
    Raises:
        OSError: If directory creation fails due to permissions.
    """
```

#### C. Add validation method (after create_all_directories method)
```python
def validate_structure(self) -> bool:
    """Validate that all expected directories exist.
    
    Returns:
        bool: True if all directories exist, False otherwise.
    
    Example:
        >>> if not plot_paths.validate_structure():
        ...     plot_paths.create_all_directories()
    """
    dirs = [
        self.cell_metrics_dir, self.gene_metrics_dir, self.metadata_dir,
        self.cell_timeline_dir, self.cell_distributions_dir, self.cell_comparisons_dir,
        self.cell_trajectories_dir, self.cell_jsd_trajectories_dir, self.cell_meth_trajectories_dir,
        self.cell_analysis_dir,
        self.gene_timeline_dir, self.gene_distributions_dir, self.gene_comparisons_dir,
        self.gene_trajectories_dir, self.gene_per_gene_dir, self.gene_analysis_dir
    ]
    return all(os.path.exists(d) for d in dirs)
```

#### D. Add get_all_paths method (after validate_structure)
```python
def get_all_paths(self) -> dict:
    """Return dictionary of all configured paths.
    
    Returns:
        Dict mapping path names to their full paths.
        
    Example:
        >>> paths = plot_paths.get_all_paths()
        >>> print(f"Cell metrics at: {paths['cell_metrics_dir']}")
    """
    return {
        'results_dir': self.results_dir,
        'cell_metrics_dir': self.cell_metrics_dir,
        'gene_metrics_dir': self.gene_metrics_dir,
        'metadata_dir': self.metadata_dir,
        'cell_timeline_dir': self.cell_timeline_dir,
        'cell_distributions_dir': self.cell_distributions_dir,
        'cell_comparisons_dir': self.cell_comparisons_dir,
        'cell_trajectories_dir': self.cell_trajectories_dir,
        'cell_jsd_trajectories_dir': self.cell_jsd_trajectories_dir,
        'cell_meth_trajectories_dir': self.cell_meth_trajectories_dir,
        'cell_analysis_dir': self.cell_analysis_dir,
        'gene_timeline_dir': self.gene_timeline_dir,
        'gene_distributions_dir': self.gene_distributions_dir,
        'gene_comparisons_dir': self.gene_comparisons_dir,
        'gene_trajectories_dir': self.gene_trajectories_dir,
        'gene_per_gene_dir': self.gene_per_gene_dir,
        'gene_analysis_dir': self.gene_analysis_dir
    }
```

### 2. Individual Trajectory Subdirectory Organization
**Location**: `phase2/core/plot_paths.py`
**Changes**:

#### A. Update __init__ method to add batch-specific subdirectories (lines 23-32)
Replace current trajectory directory definitions with:
```python
# Cell metrics subdirectories
self.cell_timeline_dir = os.path.join(self.cell_metrics_dir, 'timeline')
self.cell_distributions_dir = os.path.join(self.cell_metrics_dir, 'distributions')
self.cell_comparisons_dir = os.path.join(self.cell_metrics_dir, 'comparisons')
self.cell_trajectories_dir = os.path.join(self.cell_metrics_dir, 'individual_trajectories')

# Cell JSD trajectory subdirectories by batch
self.cell_jsd_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'jsd')
self.cell_jsd_mutant_dir = os.path.join(self.cell_jsd_trajectories_dir, 'mutant')
self.cell_jsd_control1_dir = os.path.join(self.cell_jsd_trajectories_dir, 'control1')

# Cell methylation trajectory subdirectories by batch
self.cell_meth_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'methylation_proportion')
self.cell_meth_mutant_dir = os.path.join(self.cell_meth_trajectories_dir, 'mutant')
self.cell_meth_control1_dir = os.path.join(self.cell_meth_trajectories_dir, 'control1')

self.cell_analysis_dir = os.path.join(self.cell_metrics_dir, 'analysis')

# Gene metrics subdirectories
self.gene_timeline_dir = os.path.join(self.gene_metrics_dir, 'timeline')
self.gene_distributions_dir = os.path.join(self.gene_metrics_dir, 'distributions')
self.gene_comparisons_dir = os.path.join(self.gene_metrics_dir, 'comparisons')
self.gene_trajectories_dir = os.path.join(self.gene_metrics_dir, 'individual_trajectories')

# Gene trajectory subdirectories by batch
self.gene_trajectories_mutant_dir = os.path.join(self.gene_trajectories_dir, 'mutant')
self.gene_trajectories_control1_dir = os.path.join(self.gene_trajectories_dir, 'control1')

self.gene_per_gene_dir = os.path.join(self.gene_metrics_dir, 'per_gene')
self.gene_analysis_dir = os.path.join(self.gene_metrics_dir, 'analysis')
```

#### B. Update create_all_directories method (line 37-46)
```python
def create_all_directories(self):
    """Create all necessary subdirectories for organized output."""
    dirs = [
        # Main directories
        self.cell_metrics_dir, self.gene_metrics_dir, self.metadata_dir,
        # Cell metric subdirectories
        self.cell_timeline_dir, self.cell_distributions_dir, self.cell_comparisons_dir,
        self.cell_trajectories_dir, self.cell_jsd_trajectories_dir, 
        self.cell_jsd_mutant_dir, self.cell_jsd_control1_dir,
        self.cell_meth_trajectories_dir,
        self.cell_meth_mutant_dir, self.cell_meth_control1_dir,
        self.cell_analysis_dir,
        # Gene metric subdirectories
        self.gene_timeline_dir, self.gene_distributions_dir, self.gene_comparisons_dir,
        self.gene_trajectories_dir, 
        self.gene_trajectories_mutant_dir, self.gene_trajectories_control1_dir,
        self.gene_per_gene_dir, self.gene_analysis_dir
    ]
    for dir_path in dirs:
        os.makedirs(dir_path, exist_ok=True)
```

#### C. Update path getter methods (lines 67-71, 89-91)
```python
def get_individual_cell_jsd_path(self, batch: str, individual_id: int) -> str:
    """Get path for individual cell JSD trajectory plot."""
    if batch == 'mutant':
        return os.path.join(self.cell_jsd_mutant_dir, f'individual_{individual_id:02d}.png')
    elif batch == 'control1':
        return os.path.join(self.cell_jsd_control1_dir, f'individual_{individual_id:02d}.png')
    else:
        raise ValueError(f"Unexpected batch type: {batch}")

def get_individual_cell_methylation_path(self, batch: str, individual_id: int) -> str:
    """Get path for individual cell methylation proportion trajectory plot."""
    if batch == 'mutant':
        return os.path.join(self.cell_meth_mutant_dir, f'individual_{individual_id:02d}.png')
    elif batch == 'control1':
        return os.path.join(self.cell_meth_control1_dir, f'individual_{individual_id:02d}.png')
    else:
        raise ValueError(f"Unexpected batch type: {batch}")

def get_individual_gene_jsd_path(self, batch: str, individual_id: int) -> str:
    """Get path for individual gene JSD trajectory plot."""
    if batch == 'mutant':
        return os.path.join(self.gene_trajectories_mutant_dir, f'individual_{individual_id:02d}.png')
    elif batch == 'control1':
        return os.path.join(self.gene_trajectories_control1_dir, f'individual_{individual_id:02d}.png')
    else:
        raise ValueError(f"Unexpected batch type: {batch}")
```

### 3. Create Test File for PlotPaths
**Location**: New file `phase2/tests/unit/test_plot_paths.py`
**Content**: Complete test file following project patterns

### 4. Update PLOT_DOCUMENTATION.md
**Location**: `/Users/alessandromassaad/jpi-methylation-simulation/PLOT_DOCUMENTATION.md`
**Changes**:
- Update all file paths to reflect new subdirectory structure
- Change "methylation" to "cell_methylation_proportion" where applicable
- Add note about PlotPaths class usage
- Update directory structure examples

### 5. Fix Verbose Output Messages
**Location**: `phase2/run_pipeline.py`
**Changes**:

#### Line 654:
Change: `print("STAGE 2: Plot JSD Distribution")`
To: `print("STAGE 2: Plot cell JSD Distribution")`

#### Line 837:
Change: `print(f"STAGE 5.5: Plot Year {args.second_snapshot} JSD Distribution")`
To: `print(f"STAGE 5.5: Plot Year {args.second_snapshot} cell JSD Distribution")`

#### Lines 1556-1558:
Change:
```python
print(f"  Mutant mean JSD: {summary_stats['mutant']['mean']:.6f} ± {summary_stats['mutant']['std']:.6f}")
print(f"  Control1 mean JSD: {summary_stats['control1']['mean']:.6f} ± {summary_stats['control1']['std']:.6f}")
print(f"  Control2 mean JSD: {summary_stats['control2']['mean']:.6f} ± {summary_stats['control2']['std']:.6f}")
```
To:
```python
print(f"  Mutant mean cell JSD: {summary_stats['mutant']['mean']:.6f} ± {summary_stats['mutant']['std']:.6f}")
print(f"  Control1 mean cell JSD: {summary_stats['control1']['mean']:.6f} ± {summary_stats['control1']['std']:.6f}")
print(f"  Control2 mean cell JSD: {summary_stats['control2']['mean']:.6f} ± {summary_stats['control2']['std']:.6f}")
```

## Testing Plan

1. After making all changes, run:
   - `python test_plot_paths.py` to verify PlotPaths functionality
   - Run a quick test pipeline to verify directory creation
   - Check that all plots are generated in correct subdirectories

2. Verify:
   - All directories are created with correct structure
   - Individual trajectory plots are organized by batch
   - Path generation methods work correctly
   - Documentation is clear and accurate

## Execution Order

1. Fix PlotPaths class documentation
2. Implement subdirectory organization in PlotPaths
3. Create test file for PlotPaths
4. Update PLOT_DOCUMENTATION.md
5. Fix verbose output messages in run_pipeline.py
6. Run all tests to verify changes