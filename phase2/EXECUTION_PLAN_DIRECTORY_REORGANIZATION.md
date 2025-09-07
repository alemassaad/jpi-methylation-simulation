# Execution Plan: Results Directory Reorganization

## Overview
Reorganize the results directory to group plots by metric type (cell vs gene) with clear subdirectory structure for better organization and discoverability.

## Current Structure
```
results/
├── cell_jsd_analysis.json
├── cell_jsd_comparison.png
├── cell_methylation_analysis.json
├── cell_methylation_comparison.png
├── gene_jsd_analysis.json
├── gene_jsd_individual_comparison.png
├── gene_jsd_plots/
│   └── gene_XXX_jsd_comparison.png (20 files)
├── individual_trajectories/
│   ├── mutant_XX_jsd.png
│   ├── mutant_XX_methylation.png
│   ├── mutant_XX_gene_jsd.png
│   ├── control1_XX_jsd.png
│   ├── control1_XX_methylation.png
│   └── control1_XX_gene_jsd.png
├── mixing_statistics.json
├── original_simulation_combined_timeline.png
├── original_simulation_gene_jsd_timeline.png
├── original_simulation_jsd_timeline.png
├── original_simulation_methylation_timeline.png
├── pipeline_metadata.json
├── year30_gene_jsd_histogram.png
├── year30_jsd_distribution_200bins.png
├── year30_methylation_histogram.png
├── year50_gene_jsd_histogram.png
├── year50_jsd_distribution_200bins.png
└── year50_methylation_histogram.png
```

## Target Structure
```
results/
├── cell_metrics/
│   ├── timeline/
│   │   ├── cell_jsd_timeline.png
│   │   └── cell_methylation_proportion_timeline.png
│   ├── distributions/
│   │   ├── year30_cell_jsd.png
│   │   ├── year50_cell_jsd.png
│   │   ├── year30_cell_methylation_proportion.png
│   │   └── year50_cell_methylation_proportion.png
│   ├── comparisons/
│   │   ├── cell_jsd_comparison.png
│   │   └── cell_methylation_proportion_comparison.png
│   ├── individual_trajectories/
│   │   ├── jsd/
│   │   │   ├── mutant_01_cell_jsd.png
│   │   │   ├── mutant_02_cell_jsd.png
│   │   │   ├── control1_01_cell_jsd.png
│   │   │   └── ...
│   │   └── methylation_proportion/
│   │       ├── mutant_01_cell_methylation_proportion.png
│   │       ├── mutant_02_cell_methylation_proportion.png
│   │       ├── control1_01_cell_methylation_proportion.png
│   │       └── ...
│   └── analysis/
│       ├── cell_jsd_analysis.json
│       └── cell_methylation_proportion_analysis.json
│
├── gene_metrics/
│   ├── timeline/
│   │   └── gene_jsd_timeline.png
│   ├── distributions/
│   │   ├── year30_gene_jsd.png
│   │   └── year50_gene_jsd.png
│   ├── comparisons/
│   │   └── gene_jsd_comparison.png
│   ├── individual_trajectories/
│   │   ├── mutant_01_gene_jsd.png
│   │   ├── mutant_02_gene_jsd.png
│   │   ├── control1_01_gene_jsd.png
│   │   └── ...
│   ├── per_gene/
│   │   ├── gene_000_jsd.png
│   │   ├── gene_001_jsd.png
│   │   └── ... (20 total)
│   └── analysis/
│       └── gene_jsd_analysis.json
│
└── metadata/
    ├── pipeline_metadata.json
    └── mixing_statistics.json (if --uniform-mixing)
```

## Implementation Steps

### Step 1: Update Path Generation Functions

Create new file: `phase2/core/plot_paths.py`
```python
"""
Centralized path generation for organized plot output structure.
"""
import os
from typing import Optional

class PlotPaths:
    """Generate organized paths for plot outputs."""
    
    def __init__(self, results_dir: str):
        self.results_dir = results_dir
        
        # Create main metric directories
        self.cell_metrics_dir = os.path.join(results_dir, 'cell_metrics')
        self.gene_metrics_dir = os.path.join(results_dir, 'gene_metrics')
        self.metadata_dir = os.path.join(results_dir, 'metadata')
        
        # Cell metrics subdirectories
        self.cell_timeline_dir = os.path.join(self.cell_metrics_dir, 'timeline')
        self.cell_distributions_dir = os.path.join(self.cell_metrics_dir, 'distributions')
        self.cell_comparisons_dir = os.path.join(self.cell_metrics_dir, 'comparisons')
        self.cell_trajectories_dir = os.path.join(self.cell_metrics_dir, 'individual_trajectories')
        self.cell_jsd_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'jsd')
        self.cell_meth_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'methylation_proportion')
        self.cell_analysis_dir = os.path.join(self.cell_metrics_dir, 'analysis')
        
        # Gene metrics subdirectories
        self.gene_timeline_dir = os.path.join(self.gene_metrics_dir, 'timeline')
        self.gene_distributions_dir = os.path.join(self.gene_metrics_dir, 'distributions')
        self.gene_comparisons_dir = os.path.join(self.gene_metrics_dir, 'comparisons')
        self.gene_trajectories_dir = os.path.join(self.gene_metrics_dir, 'individual_trajectories')
        self.gene_per_gene_dir = os.path.join(self.gene_metrics_dir, 'per_gene')
        self.gene_analysis_dir = os.path.join(self.gene_metrics_dir, 'analysis')
    
    def create_all_directories(self):
        """Create all necessary directories."""
        dirs = [
            self.cell_metrics_dir, self.gene_metrics_dir, self.metadata_dir,
            self.cell_timeline_dir, self.cell_distributions_dir, self.cell_comparisons_dir,
            self.cell_trajectories_dir, self.cell_jsd_trajectories_dir, self.cell_meth_trajectories_dir,
            self.cell_analysis_dir,
            self.gene_timeline_dir, self.gene_distributions_dir, self.gene_comparisons_dir,
            self.gene_trajectories_dir, self.gene_per_gene_dir, self.gene_analysis_dir
        ]
        for dir_path in dirs:
            os.makedirs(dir_path, exist_ok=True)
    
    # Cell metric paths
    def get_cell_jsd_timeline_path(self) -> str:
        return os.path.join(self.cell_timeline_dir, 'cell_jsd_timeline.png')
    
    def get_cell_methylation_timeline_path(self) -> str:
        return os.path.join(self.cell_timeline_dir, 'cell_methylation_proportion_timeline.png')
    
    def get_cell_jsd_distribution_path(self, year: int) -> str:
        return os.path.join(self.cell_distributions_dir, f'year{year}_cell_jsd.png')
    
    def get_cell_methylation_distribution_path(self, year: int) -> str:
        return os.path.join(self.cell_distributions_dir, f'year{year}_cell_methylation_proportion.png')
    
    def get_cell_jsd_comparison_path(self) -> str:
        return os.path.join(self.cell_comparisons_dir, 'cell_jsd_comparison.png')
    
    def get_cell_methylation_comparison_path(self) -> str:
        return os.path.join(self.cell_comparisons_dir, 'cell_methylation_proportion_comparison.png')
    
    def get_individual_cell_jsd_path(self, batch: str, individual_id: int) -> str:
        return os.path.join(self.cell_jsd_trajectories_dir, f'{batch}_{individual_id:02d}_cell_jsd.png')
    
    def get_individual_cell_methylation_path(self, batch: str, individual_id: int) -> str:
        return os.path.join(self.cell_meth_trajectories_dir, f'{batch}_{individual_id:02d}_cell_methylation_proportion.png')
    
    def get_cell_jsd_analysis_path(self) -> str:
        return os.path.join(self.cell_analysis_dir, 'cell_jsd_analysis.json')
    
    def get_cell_methylation_analysis_path(self) -> str:
        return os.path.join(self.cell_analysis_dir, 'cell_methylation_proportion_analysis.json')
    
    # Gene metric paths
    def get_gene_jsd_timeline_path(self) -> str:
        return os.path.join(self.gene_timeline_dir, 'gene_jsd_timeline.png')
    
    def get_gene_jsd_distribution_path(self, year: int) -> str:
        return os.path.join(self.gene_distributions_dir, f'year{year}_gene_jsd.png')
    
    def get_gene_jsd_comparison_path(self) -> str:
        return os.path.join(self.gene_comparisons_dir, 'gene_jsd_comparison.png')
    
    def get_individual_gene_jsd_path(self, batch: str, individual_id: int) -> str:
        return os.path.join(self.gene_trajectories_dir, f'{batch}_{individual_id:02d}_gene_jsd.png')
    
    def get_per_gene_jsd_path(self, gene_idx: int) -> str:
        return os.path.join(self.gene_per_gene_dir, f'gene_{gene_idx:03d}_jsd.png')
    
    def get_gene_jsd_analysis_path(self) -> str:
        return os.path.join(self.gene_analysis_dir, 'gene_jsd_analysis.json')
    
    # Metadata paths
    def get_pipeline_metadata_path(self) -> str:
        return os.path.join(self.metadata_dir, 'pipeline_metadata.json')
    
    def get_mixing_statistics_path(self) -> str:
        return os.path.join(self.metadata_dir, 'mixing_statistics.json')
```

### Step 2: Update run_pipeline.py

#### Import PlotPaths (Add after other imports)
```python
from core.plot_paths import PlotPaths
```

#### Initialize PlotPaths (After creating results_dir, ~Line 600)
```python
# Initialize organized path structure
plot_paths = PlotPaths(results_dir)
plot_paths.create_all_directories()
```

#### Update All Path Generations

**Stage 5: First Snapshot Plots (~Lines 642-653)**
```python
# OLD
plot_path = os.path.join(results_dir, f"year{args.first_snapshot}_jsd_distribution_{args.bins}bins.png")
# NEW
plot_path = plot_paths.get_cell_jsd_distribution_path(args.first_snapshot)

# OLD
methylation_plot_path = os.path.join(results_dir, f"year{args.first_snapshot}_methylation_histogram.png")
# NEW
methylation_plot_path = plot_paths.get_cell_methylation_distribution_path(args.first_snapshot)

# OLD
gene_jsd_plot_path = os.path.join(results_dir, f"year{args.first_snapshot}_gene_jsd_histogram.png")
# NEW
gene_jsd_plot_path = plot_paths.get_gene_jsd_distribution_path(args.first_snapshot)
```

**Stage 6: Second Snapshot Plots (~Lines 825-836)**
```python
# OLD
second_plot_path = os.path.join(results_dir, f"year{args.second_snapshot}_jsd_distribution_{args.bins}bins.png")
# NEW
second_plot_path = plot_paths.get_cell_jsd_distribution_path(args.second_snapshot)

# Similar updates for methylation and gene JSD paths
```

**Stage 7: Individual Trajectories (~Lines 1233-1269)**
```python
# Mutant trajectories
# OLD
jsd_path = os.path.join(individual_trajectories_dir, f"mutant_{i:02d}_jsd.png")
meth_path = os.path.join(individual_trajectories_dir, f"mutant_{i:02d}_methylation.png")
gene_jsd_path = os.path.join(individual_trajectories_dir, f"mutant_{i:02d}_gene_jsd.png")

# NEW
jsd_path = plot_paths.get_individual_cell_jsd_path('mutant', i)
meth_path = plot_paths.get_individual_cell_methylation_path('mutant', i)
gene_jsd_path = plot_paths.get_individual_gene_jsd_path('mutant', i)

# Similar for control1
```

**Stage 8: Comparison Plots**
```python
# Cell JSD comparison (update analyze_cell_jsd_comparison to use new path)
# Cell methylation comparison (update analyze_cell_methylation_proportion_comparison to use new path)
# Gene JSD comparison (~Line 1337)
# OLD
gene_dist_path = os.path.join(results_dir, "gene_jsd_distribution.png")
# NEW
gene_dist_path = plot_paths.get_gene_jsd_comparison_path()
```

**Pipeline Metadata (~Line 1343)**
```python
# OLD
metadata_path = os.path.join(results_dir, 'pipeline_metadata.json')
# NEW
metadata_path = plot_paths.get_pipeline_metadata_path()
```

**Mixing Statistics (if uniform mixing)**
```python
# OLD
stats_path = os.path.join(results_dir, 'mixing_statistics.json')
# NEW
stats_path = plot_paths.get_mixing_statistics_path()
```

**Original Simulation Timelines (~Lines 1488-1505)**
```python
# OLD
jsd_timeline_path = os.path.join(results_dir, "original_simulation_jsd_timeline.png")
# NEW
jsd_timeline_path = plot_paths.get_cell_jsd_timeline_path()

# OLD
meth_timeline_path = os.path.join(results_dir, "original_simulation_methylation_timeline.png")
# NEW
meth_timeline_path = plot_paths.get_cell_methylation_timeline_path()

# OLD (DELETE THIS)
combined_timeline_path = os.path.join(results_dir, "original_simulation_combined_timeline.png")
# Remove the combined timeline generation code

# OLD
gene_jsd_timeline_path = os.path.join(results_dir, "original_simulation_gene_jsd_timeline.png")
# NEW
gene_jsd_timeline_path = plot_paths.get_gene_jsd_timeline_path()
```

### Step 3: Update pipeline_analysis.py

#### Update plot_individual_gene_jsds function (~Line 1630)
```python
# OLD
plot_path = os.path.join(gene_plots_dir, f'gene_{gene_idx:03d}_jsd_comparison.png')
# NEW (pass plot_paths object)
plot_path = plot_paths.get_per_gene_jsd_path(gene_idx)
```

#### Update analyze_cell_jsd_comparison (~Line 837)
```python
# OLD
plot_path = os.path.join(output_dir, "cell_jsd_comparison.png")
# NEW (pass plot_paths object)
plot_path = plot_paths.get_cell_jsd_comparison_path()

# OLD
json_path = os.path.join(output_dir, "cell_jsd_analysis.json")
# NEW
json_path = plot_paths.get_cell_jsd_analysis_path()
```

#### Update analyze_cell_methylation_proportion_comparison (~Line 554)
```python
# OLD
plot_path = os.path.join(output_dir, "cell_methylation_comparison.png")
# NEW (pass plot_paths object)
plot_path = plot_paths.get_cell_methylation_comparison_path()

# OLD
json_path = os.path.join(output_dir, "cell_methylation_analysis.json")
# NEW
json_path = plot_paths.get_cell_methylation_analysis_path()
```

#### Update plot_gene_jsd_comparison (~Line 1926)
```python
# OLD
output_path = os.path.join(output_dir, 'gene_jsd_individual_comparison.png')
# NEW (pass plot_paths object)
output_path = plot_paths.get_gene_jsd_comparison_path()
```

### Step 4: Update Function Signatures

Since we're passing plot_paths object, update function signatures:

```python
def analyze_cell_jsd_comparison(mutant_dishes, control1_dishes, control2_dishes, 
                                plot_paths: PlotPaths, verbose=False):
    # Use plot_paths.get_cell_jsd_comparison_path() internally
    
def analyze_cell_methylation_proportion_comparison(mutant_dishes, control1_dishes, control2_dishes,
                                                   plot_paths: PlotPaths, verbose=False):
    # Use plot_paths.get_cell_methylation_comparison_path() internally

def plot_individual_gene_jsds(individual_gene_jsds, plot_paths: PlotPaths):
    # Use plot_paths.get_per_gene_jsd_path(gene_idx) internally
    
def plot_gene_jsd_comparison(batches_data, plot_paths: PlotPaths):
    # Use plot_paths.get_gene_jsd_comparison_path() internally
```

### Step 5: Remove Combined Timeline

In `phase1/cell.py`, remove or comment out:
1. The `plot_combined()` method (~Line 1850-1914)
2. Any calls to this method

In `run_pipeline.py`:
1. Remove combined timeline generation (~Lines 1497-1500)
2. Remove the print statement for combined timeline

### Step 6: Update Test Files

#### Update Expected Paths in Tests

File: `phase2/tests/test_all_new_plots.py`
```python
# Update expected plot locations to match new structure
expected_plots = {
    "cell_metrics/timeline/cell_jsd_timeline.png",
    "cell_metrics/distributions/year30_cell_jsd.png",
    "cell_metrics/comparisons/cell_jsd_comparison.png",
    # ... etc
}
```

File: `phase2/tests/test_new_plots_implementation.py`
```python
# Update all path assertions to match new structure
```

### Step 7: Update plot_individuals.py Standalone Script

```python
# Import PlotPaths
from core.plot_paths import PlotPaths

# Initialize when processing
plot_paths = PlotPaths(output_dir)
plot_paths.create_all_directories()

# Use plot_paths for generating output paths
```

### Step 8: Update CLAUDE.md Documentation

Update the documentation to reflect new directory structure:

```markdown
### 📊 Generated Plots (Phase 2)
All plots are organized by metric type:

**Cell Metrics (`results/cell_metrics/`):**
- `timeline/`: Evolution over time
  - `cell_jsd_timeline.png`: Full JSD evolution from original simulation
  - `cell_methylation_proportion_timeline.png`: Methylation proportion evolution
- `distributions/`: Snapshot distributions
  - `year{X}_cell_jsd.png`: Cell JSD distribution at year X
  - `year{X}_cell_methylation_proportion.png`: Methylation distribution at year X
- `comparisons/`: Batch comparisons
  - `cell_jsd_comparison.png`: Box plots comparing batches
  - `cell_methylation_proportion_comparison.png`: Methylation comparison
- `individual_trajectories/`: Individual growth trajectories
  - `jsd/`: Cell JSD trajectories for each individual
  - `methylation_proportion/`: Methylation trajectories
- `analysis/`: JSON analysis files
  - `cell_jsd_analysis.json`: Statistical analysis
  - `cell_methylation_proportion_analysis.json`: Methylation statistics

**Gene Metrics (`results/gene_metrics/`):**
- `timeline/`: Gene JSD evolution
  - `gene_jsd_timeline.png`: Per-gene JSD over time
- `distributions/`: Gene JSD distributions
  - `year{X}_gene_jsd.png`: Gene JSD at snapshot years
- `comparisons/`: Gene-level comparisons
  - `gene_jsd_comparison.png`: Individual averages comparison
- `individual_trajectories/`: Gene JSD trajectories
  - `{batch}_{id:02d}_gene_jsd.png`: Per-individual gene trajectories
- `per_gene/`: Individual gene analysis
  - `gene_{idx:03d}_jsd.png`: JSD for each specific gene
- `analysis/`: Gene-level analysis
  - `gene_jsd_analysis.json`: Gene JSD statistics

**Metadata (`results/metadata/`):**
- `pipeline_metadata.json`: Pipeline configuration
- `mixing_statistics.json`: Mixing stage statistics (if --uniform-mixing)
```

## Testing Plan

### 1. Create Test Script
Create `phase2/tests/test_directory_reorganization.py`:
```python
"""Test that new directory structure is created correctly."""
import os
import tempfile
import shutil
from core.plot_paths import PlotPaths

def test_directory_creation():
    """Test all directories are created."""
    with tempfile.TemporaryDirectory() as tmpdir:
        results_dir = os.path.join(tmpdir, 'results')
        plot_paths = PlotPaths(results_dir)
        plot_paths.create_all_directories()
        
        # Check all directories exist
        assert os.path.exists(plot_paths.cell_metrics_dir)
        assert os.path.exists(plot_paths.gene_metrics_dir)
        assert os.path.exists(plot_paths.metadata_dir)
        # ... check all subdirectories
        
def test_path_generation():
    """Test path generation methods."""
    plot_paths = PlotPaths('/tmp/results')
    
    # Test cell paths
    assert plot_paths.get_cell_jsd_timeline_path() == '/tmp/results/cell_metrics/timeline/cell_jsd_timeline.png'
    assert plot_paths.get_individual_cell_jsd_path('mutant', 1) == '/tmp/results/cell_metrics/individual_trajectories/jsd/mutant_01_cell_jsd.png'
    
    # Test gene paths
    assert plot_paths.get_gene_jsd_timeline_path() == '/tmp/results/gene_metrics/timeline/gene_jsd_timeline.png'
    assert plot_paths.get_per_gene_jsd_path(5) == '/tmp/results/gene_metrics/per_gene/gene_005_jsd.png'
    
    print("✓ All path generation tests passed")

if __name__ == "__main__":
    test_directory_creation()
    test_path_generation()
```

### 2. Integration Test
Run small pipeline to verify new structure:
```bash
cd phase2
python run_pipeline.py \
    --simulation ../phase1/data/.../simulation.json \
    --first-snapshot 5 \
    --second-snapshot 10 \
    --n-quantiles 2 \
    --cells-per-quantile 1

# Verify directory structure
find data/.../results -type f -name "*.png" | sort
```

### 3. Verification Checklist
- [ ] All directories created
- [ ] Plots saved to correct locations
- [ ] JSON files in analysis subdirectories
- [ ] No files in root results directory (except metadata/)
- [ ] Individual trajectories properly separated
- [ ] Combined timeline removed
- [ ] Tests pass with new structure

## Migration Support

For existing results, create migration script:
```python
"""Migrate old results structure to new organized structure."""
import os
import shutil
from core.plot_paths import PlotPaths

def migrate_results(old_results_dir: str):
    """Migrate old flat structure to new organized structure."""
    # Create new structure
    plot_paths = PlotPaths(old_results_dir)
    plot_paths.create_all_directories()
    
    # Define migration mappings
    migrations = {
        'cell_jsd_comparison.png': plot_paths.get_cell_jsd_comparison_path(),
        'original_simulation_jsd_timeline.png': plot_paths.get_cell_jsd_timeline_path(),
        # ... complete mapping
    }
    
    # Move files
    for old_name, new_path in migrations.items():
        old_path = os.path.join(old_results_dir, old_name)
        if os.path.exists(old_path):
            shutil.move(old_path, new_path)
            print(f"Moved {old_name} → {new_path}")
    
    # Delete combined timeline if exists
    combined_path = os.path.join(old_results_dir, 'original_simulation_combined_timeline.png')
    if os.path.exists(combined_path):
        os.remove(combined_path)
        print("Removed combined timeline")
```

## Rollback Plan

If issues arise:
1. Git revert the commit
2. Restore flat structure
3. Remove PlotPaths class
4. Revert path generation changes

## Success Criteria

1. All plots generate in correct subdirectories
2. No plots in root results directory
3. Clear separation between cell and gene metrics
4. Individual trajectories organized by metric type
5. Combined timeline successfully removed
6. All tests pass with new structure
7. Documentation updated to reflect new structure

## Benefits of New Structure

1. **Clear Organization**: Instantly know what type of metric each plot shows
2. **Easier Navigation**: Related plots grouped together
3. **Scalability**: Easy to add new plot types in appropriate locations
4. **Reduced Clutter**: Root directory only has metadata
5. **Logical Grouping**: Individual trajectories with their metric type
6. **Better Discoverability**: Users can find plots by browsing metric categories