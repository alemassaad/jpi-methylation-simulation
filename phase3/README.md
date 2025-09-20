# Phase 3: Analysis and Visualization Pipeline

Pure analysis pipeline that reads data from phase2 and generates all plots and statistical analysis.

## Overview

Phase 3 is the **pure analysis pipeline** that reads structured datasets from phase2 and generates comprehensive visualizations and statistical analysis. This complete separation from simulation provides:

### Key Benefits
- **Efficiency**: Re-analyze without expensive re-simulation
- **Flexibility**: Multiple analysis approaches on same dataset
- **Scalability**: Batch processing of multiple phase2 runs
- **Development**: Independent analysis workflow iteration
- **Reproducibility**: Consistent analysis across different datasets

### What Phase 3 Does
- Load phase2 snapshots and individual populations
- Generate cell-level and gene-level distribution plots
- Create population comparison visualizations
- Perform statistical analysis and significance testing
- Generate individual growth trajectory plots
- Create timeline visualizations from original phase1 data
- Organize all outputs in structured directory hierarchy

### What Phase 3 Does NOT Do
- Does not run simulations (phase1)
- Does not generate or modify cell populations (phase2)
- Does not create new datasets
- Only reads and analyzes existing phase2 data

## Quick Start

### Prerequisites
Ensure you have:
1. **Phase1 simulation**: Complete simulation file (`.json.gz`)
2. **Phase2 data**: Complete phase2 pipeline output with `snapshots/` and `individuals/` directories
3. **Dependencies**: Install with `pip install -r requirements.txt`

### Basic Usage
```bash
cd phase3

# Standard analysis
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{simulation_directory}/simulation.json.gz

# With wildcard pattern matching
python run_analysis.py \
    --phase2-dir ../phase2/data/gene_rates_*/snap30to50*/ \
    --simulation ../phase1/data/gene_rates_*/simulation.json.gz
```

### Custom Analysis Parameters
```bash
# Fast analysis (fewer bins, limited gene plots)
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --bins 100 \
    --max-gene-plots 5

# Comprehensive analysis (high resolution)
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --bins 300 \
    --output-dir detailed_analysis_$(date +%Y%m%d)
```

### Using Configuration Files
```bash
# Quick analysis preset
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --config configs/quick_analysis.yaml

# Custom config
python run_analysis.py \
    --phase2-dir ../phase2/data/{run_directory}/ \
    --simulation ../phase1/data/{sim}/simulation.json.gz \
    --config my_custom_config.yaml
```

## Command-Line Options

### Required Arguments
- `--phase2-dir`: Path to phase2 output directory
  - Must contain `snapshots/` and `individuals/` subdirectories
  - Example: `../phase2/data/gene_rates_*/snap30to50-growth7-*/`
  - Supports glob patterns for batch processing

- `--simulation`: Path to original phase1 simulation file
  - Used for timeline plots and parameter validation
  - Must be the same simulation used for phase2 data generation
  - Example: `../phase1/data/gene_rates_*/simulation.json.gz`

### Analysis Control Options
- `--bins`: Number of histogram bins (default: 200)
  - Higher values: more detailed distributions, slower analysis
  - Lower values: faster analysis, less detail
  - Typical range: 50-500

- `--max-gene-plots`: Limit per-gene plots (default: all genes)
  - Set to small number (5-10) for quick analysis
  - `null` or omit for all genes (typically 20 genes)
  - Useful for large datasets or quick testing

### Output Control Options
- `--output-dir`: Custom output directory name
  - Default: auto-generated `analysis_bins{N}_{timestamp}`
  - Relative to `phase3/data/` directory
  - Use for organized batch processing

- `--config`: YAML configuration file path
  - Overrides command-line defaults
  - Useful for standardized analysis protocols
  - See `configs/` directory for examples

## Configuration Files

YAML configuration files provide reproducible analysis settings and are recommended for production workflows.

### Available Configurations

#### `configs/default.yaml` (Standard Analysis)
```yaml
# Default configuration for comprehensive analysis
bins: 200                    # High-resolution histograms
max_gene_plots: null         # Plot all genes
# output_dir: null           # Auto-generated directory
```

#### `configs/quick_analysis.yaml` (Fast Analysis)
```yaml
# Quick analysis for testing and iteration
bins: 100                    # Lower resolution for speed
max_gene_plots: 5            # Only first 5 genes
# output_dir: quick_results  # Optional custom directory
```

### Custom Configuration Example
```yaml
# configs/publication_analysis.yaml
bins: 300                    # High resolution for publication
max_gene_plots: 20           # All genes
output_dir: publication_2024 # Organized output
```

### Configuration Priority
1. **Command-line arguments** (highest priority)
2. **User config file** (`--config`)
3. **Default config** (`configs/default.yaml`)
4. **Hardcoded defaults** (lowest priority)

### Using Configurations
```bash
# Use default config (recommended)
python run_analysis.py --phase2-dir ... --simulation ... --config configs/default.yaml

# Override specific parameters
python run_analysis.py --phase2-dir ... --simulation ... --config configs/quick_analysis.yaml --bins 150
```

## Output Structure

By default, results are saved to:
```
phase3/data/analysis_bins{N}_{timestamp}/results/
```

Or with custom `--output-dir`:
```
{output_dir}/results/
```

Directory structure:
```
results/
├── cell_metrics/
│   ├── distributions/          # Snapshot histograms
│   │   ├── year30_cell_jsd.png
│   │   ├── year30_cell_methylation_proportion.png
│   │   ├── year50_cell_jsd.png
│   │   └── year50_cell_methylation_proportion.png
│   ├── comparisons/            # Batch comparisons
│   │   ├── cell_jsd_comparison.png
│   │   └── cell_methylation_proportion_comparison.png
│   ├── individual_trajectories/ # Growth trajectories
│   │   ├── jsd/
│   │   │   ├── mutant/
│   │   │   └── control1/
│   │   └── methylation_proportion/
│   │       ├── mutant/
│   │       └── control1/
│   └── timeline/               # Original simulation timelines
│       ├── cell_jsd_timeline.png
│       └── cell_methylation_proportion_timeline.png
├── gene_metrics/
│   ├── distributions/          # Gene-level distributions
│   ├── comparisons/            # Gene-level comparisons
│   ├── per_gene/              # Individual gene plots
│   │   ├── gene_000_jsd.png
│   │   └── ... (up to gene_019_jsd.png)
│   ├── individual_trajectories/
│   └── timeline/
│       └── gene_jsd_timeline.png
├── metadata/
│   └── analysis_metadata.json
└── simulation_gene_jsd_heatmap.png
```

## Analysis Stages

Phase 3 performs comprehensive analysis in four main stages:

### Stage 1: Snapshot Distribution Analysis
**Purpose**: Analyze cell populations at specific time points

**Cell-Level Metrics**:
- Cell JSD distributions with statistical overlays
- Cell methylation proportion histograms
- Population statistics (mean, median, SD, CV, MAD, percentiles)

**Gene-Level Metrics**:
- Gene JSD distributions across population
- Gene methylation proportion analysis
- Per-gene heterogeneity measurements

**Outputs**:
- `year{N}_cell_jsd.png`: Cell JSD distribution for snapshot
- `year{N}_cell_methylation_proportion.png`: Cell methylation histogram
- `year{N}_gene_jsd.png`: Gene JSD distribution
- `year{N}_gene_methylation.png`: Gene methylation distribution

### Stage 2: Population Comparison Analysis
**Purpose**: Statistical comparison between experimental groups

**Comparisons Performed**:
- **Mutant vs Control1**: Quantile sampling vs random sampling effects
- **Mutant vs Control2**: Growth effects vs pure snapshot effects
- **Control1 vs Control2**: Random sampling vs snapshot composition

**Statistical Tests**:
- T-tests for significance testing
- Effect size calculations
- Distribution overlap analysis

**Outputs**:
- `cell_jsd_comparison.png`: Three-group cell JSD comparison
- `cell_methylation_proportion_comparison.png`: Three-group methylation comparison
- `gene_jsd_comparison.png`: Gene-level population comparisons
- `gene_methylation_comparison.png`: Gene methylation comparisons

### Stage 3: Individual Trajectory Analysis
**Purpose**: Track individual population dynamics over time

**Trajectory Types**:
- **Cell JSD trajectories**: Individual epigenetic divergence over time
- **Methylation progression**: Accumulation patterns in each individual
- **Gene JSD evolution**: Per-gene heterogeneity development
- **Population size dynamics**: Growth and homeostasis patterns

**Grouping**:
- Separate trajectories for mutant and control1 groups
- Control2 has no trajectories (snapshot-only)

**Outputs**:
- `individual_trajectories/jsd/mutant/`: Individual cell JSD time series
- `individual_trajectories/jsd/control1/`: Control1 trajectories
- `individual_trajectories/methylation_proportion/`: Methylation progression

### Stage 4: Timeline and Historical Analysis
**Purpose**: Visualize complete simulation history from phase1

**Timeline Plots**:
- **Cell JSD timeline**: Population-wide JSD evolution
- **Methylation timeline**: Aggregate methylation accumulation
- **Gene JSD timeline**: Per-gene heterogeneity development

**Historical Analysis**:
- **Gene JSD heatmaps**: Spatial-temporal patterns
- **Rate group comparisons**: Effects of different methylation rates
- **Population statistics over time**: Mean, variance, distribution changes

**Outputs**:
- `timeline/cell_jsd_timeline.png`: Complete cell JSD history
- `timeline/cell_methylation_proportion_timeline.png`: Methylation history
- `timeline/gene_jsd_timeline.png`: Gene JSD evolution
- `simulation_gene_jsd_heatmap.png`: Comprehensive gene heatmap

### Per-Gene Detailed Analysis
**Purpose**: Individual gene characterization

**For Each Gene (0-19)**:
- Distribution analysis across populations
- Statistical comparisons between groups
- Heterogeneity measurements
- Rate group effects (if applicable)

**Outputs**:
- `per_gene/gene_000_jsd.png` through `gene_019_jsd.png`
- Limited by `--max-gene-plots` parameter for performance

## Batch Processing and Automation

### Batch Analysis Script
Process multiple phase2 runs efficiently:

```python
#!/usr/bin/env python3
"""
Batch analysis script for processing multiple phase2 runs.
Run from phase3 directory.
"""
import glob
import subprocess
import os
import sys
from pathlib import Path

def find_phase2_runs(pattern="../phase2/data/*/"):
    """Find all phase2 output directories."""
    return glob.glob(pattern)

def find_matching_simulation(phase2_dir):
    """Find corresponding phase1 simulation for a phase2 run."""
    # Extract parameters from phase2 directory name
    phase2_name = os.path.basename(phase2_dir.rstrip('/'))
    
    # Look for matching simulation
    sim_pattern = "../phase1/data/*/simulation.json.gz"
    simulations = glob.glob(sim_pattern)
    
    # Return first match (could be enhanced with smarter matching)
    return simulations[0] if simulations else None

def run_analysis(phase2_dir, simulation_path, config="configs/quick_analysis.yaml"):
    """Run analysis for a single phase2 run."""
    cmd = [
        "python", "run_analysis.py",
        "--phase2-dir", phase2_dir,
        "--simulation", simulation_path,
        "--config", config
    ]
    
    print(f"\nAnalyzing: {phase2_dir}")
    print(f"Simulation: {simulation_path}")
    print(f"Command: {' '.join(cmd)}")
    
    try:
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        print(f"✓ Success: {phase2_dir}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"❌ Failed: {phase2_dir}")
        print(f"Error: {e.stderr}")
        return False

if __name__ == "__main__":
    # Find all phase2 runs
    phase2_dirs = find_phase2_runs()
    
    if not phase2_dirs:
        print("No phase2 runs found. Check path pattern.")
        sys.exit(1)
    
    print(f"Found {len(phase2_dirs)} phase2 runs to analyze")
    
    # Process each run
    successes = 0
    for phase2_dir in phase2_dirs:
        simulation = find_matching_simulation(phase2_dir)
        if simulation:
            if run_analysis(phase2_dir, simulation):
                successes += 1
        else:
            print(f"⚠️  No matching simulation found for {phase2_dir}")
    
    print(f"\nCompleted: {successes}/{len(phase2_dirs)} analyses")
```

### Shell Script Alternative
```bash
#!/bin/bash
# batch_analysis.sh - Simple batch processing

for phase2_dir in ../phase2/data/*/; do
    echo "Analyzing: $phase2_dir"
    
    python run_analysis.py \
        --phase2-dir "$phase2_dir" \
        --simulation "../phase1/data/*/simulation.json.gz" \
        --config configs/quick_analysis.yaml
    
    if [ $? -eq 0 ]; then
        echo "✓ Success: $phase2_dir"
    else
        echo "❌ Failed: $phase2_dir"
    fi
done
```

### Parallel Processing
For large-scale analysis, use GNU parallel:

```bash
# Create list of phase2 directories
find ../phase2/data -maxdepth 1 -type d -name "snap*" > phase2_dirs.txt

# Run analyses in parallel (adjust -j for CPU cores)
parallel -j 4 python run_analysis.py --phase2-dir {} --simulation ../phase1/data/*/simulation.json.gz --config configs/quick_analysis.yaml :::: phase2_dirs.txt
```

### Automated Quality Control
```python
#!/usr/bin/env python3
"""
Quality control script to validate batch analysis results.
"""
import glob
import json
from pathlib import Path

def validate_analysis_output(analysis_dir):
    """Check if analysis completed successfully."""
    results_dir = Path(analysis_dir) / "results"
    
    # Required directories
    required_dirs = [
        "cell_metrics", "gene_metrics", "metadata"
    ]
    
    # Check directory structure
    for req_dir in required_dirs:
        if not (results_dir / req_dir).exists():
            return False, f"Missing directory: {req_dir}"
    
    # Check for key plots
    cell_dist = results_dir / "cell_metrics" / "distributions"
    if not any(cell_dist.glob("*.png")):
        return False, "No cell distribution plots found"
    
    # Check metadata
    metadata_file = results_dir / "metadata" / "analysis_metadata.json"
    if not metadata_file.exists():
        return False, "Missing analysis metadata"
    
    return True, "Analysis complete"

if __name__ == "__main__":
    # Find all analysis outputs
    analysis_dirs = glob.glob("data/analysis_*/")
    
    print(f"Validating {len(analysis_dirs)} analysis outputs...")
    
    for analysis_dir in analysis_dirs:
        valid, message = validate_analysis_output(analysis_dir)
        status = "✓" if valid else "❌"
        print(f"{status} {analysis_dir}: {message}")
```

## Installation and Dependencies

### System Requirements
- **Python**: 3.7 or higher (tested on 3.8-3.11)
- **Operating System**: Linux, macOS, or Windows
- **Memory**: Minimum 4GB RAM (8GB+ recommended for large datasets)
- **Storage**: 2-5x the size of phase2 data for analysis outputs
- **CPU**: Single-threaded, benefits from high clock speed

### Installation

#### Quick Installation
```bash
# From the phase3 directory
pip install -r requirements.txt
```

#### Manual Installation
```bash
# Core numerical libraries
pip install "numpy>=1.19.0" "scipy>=1.7.0"

# Visualization libraries
pip install "plotly>=5.0.0,<6.0.0" "kaleido==0.2.1"

# Configuration support
pip install "pyyaml>=6.0"
```

#### Development Installation
```bash
# Include development and testing dependencies
pip install -r requirements.txt
pip install pytest pytest-cov jupyter  # Optional for development
```

### Required Packages

| Package | Version | Purpose |
|---------|---------|----------|
| **numpy** | ≥ 1.19.0 | Numerical computations, array operations |
| **scipy** | ≥ 1.7.0 | Statistical functions, JSD calculations |
| **plotly** | ≥ 5.0.0, < 6.0.0 | Interactive and static plot generation |
| **kaleido** | == 0.2.1 | PNG export from plotly (exact version required) |
| **pyyaml** | ≥ 6.0 | YAML configuration file support |

### Optional Packages

| Package | Purpose | Installation |
|---------|---------|-------------|
| **psutil** | System monitoring | `pip install psutil` |
| **jupyter** | Interactive analysis | `pip install jupyter` |
| **pytest** | Running tests | `pip install pytest` |
| **tqdm** | Progress bars | `pip install tqdm` |

### Version Compatibility

#### Tested Combinations
- **Python 3.8** + numpy 1.21 + scipy 1.7 + plotly 5.15
- **Python 3.9** + numpy 1.23 + scipy 1.9 + plotly 5.17
- **Python 3.10** + numpy 1.24 + scipy 1.10 + plotly 5.18
- **Python 3.11** + numpy 1.24 + scipy 1.10 + plotly 5.18

#### Known Issues
- **kaleido 0.2.1**: Exact version required (newer versions have compatibility issues)
- **plotly 6.x**: Not supported (breaking changes in API)
- **numpy < 1.19**: Missing some required array functions
- **scipy < 1.7**: Incompatible statistical function signatures

### Environment Setup

#### Virtual Environment (Recommended)
```bash
# Create virtual environment
python -m venv phase3_env

# Activate environment
source phase3_env/bin/activate  # Linux/macOS
# or
phase3_env\Scripts\activate  # Windows

# Install dependencies
pip install -r requirements.txt
```

#### Conda Environment
```bash
# Create conda environment
conda create -n phase3 python=3.9
conda activate phase3

# Install packages
conda install numpy scipy pyyaml
pip install plotly==5.18.0 kaleido==0.2.1
```

### Installation Verification

#### Quick Test
```bash
# Test imports
python -c "import numpy, scipy, plotly, kaleido, yaml; print('All imports successful')"

# Test plotly export capability
python -c "import plotly.graph_objects as go; fig = go.Figure(); fig.write_image('test.png'); print('Plot export working')"
```

#### Comprehensive Test
```python
#!/usr/bin/env python3
"""
Installation verification script for phase3.
"""
import sys
import importlib
from packaging import version

# Required packages with versions
REQUIREMENTS = {
    'numpy': '1.19.0',
    'scipy': '1.7.0', 
    'plotly': '5.0.0',
    'kaleido': '0.2.1',
    'yaml': '6.0'  # pyyaml package imports as 'yaml'
}

def check_package(name, min_version):
    """Check if package is installed with correct version."""
    try:
        module = importlib.import_module(name)
        if hasattr(module, '__version__'):
            installed_version = module.__version__
            if version.parse(installed_version) >= version.parse(min_version):
                print(f"✓ {name}: {installed_version} (>= {min_version})")
                return True
            else:
                print(f"❌ {name}: {installed_version} (need >= {min_version})")
                return False
        else:
            print(f"⚠️  {name}: version unknown")
            return True
    except ImportError:
        print(f"❌ {name}: not installed")
        return False

def test_functionality():
    """Test basic functionality."""
    try:
        import numpy as np
        import plotly.graph_objects as go
        
        # Test numpy
        arr = np.array([1, 2, 3])
        assert arr.sum() == 6
        
        # Test plotly
        fig = go.Figure(data=go.Scatter(x=[1, 2], y=[1, 2]))
        assert fig.data[0].x[0] == 1
        
        print("✓ Basic functionality test passed")
        return True
    except Exception as e:
        print(f"❌ Functionality test failed: {e}")
        return False

if __name__ == "__main__":
    print(f"Python version: {sys.version}")
    print("\nChecking package versions...")
    
    all_good = True
    for package, min_ver in REQUIREMENTS.items():
        if not check_package(package, min_ver):
            all_good = False
    
    print("\nTesting functionality...")
    if not test_functionality():
        all_good = False
    
    if all_good:
        print("\n🎉 Installation verification successful!")
        sys.exit(0)
    else:
        print("\n❌ Installation verification failed. Please fix issues above.")
        sys.exit(1)
```

#### Run Verification
```bash
# Save script as verify_installation.py and run
python verify_installation.py
```

### Troubleshooting Installation

#### Common Issues

**"No module named 'kaleido'"**
```bash
# Kaleido sometimes needs special installation
pip uninstall kaleido
pip install kaleido==0.2.1 --force-reinstall
```

**"plotly export not working"**
```bash
# Test kaleido separately
python -c "import kaleido; print('Kaleido working')"

# Alternative: install via conda-forge
conda install -c conda-forge kaleido
```

**"Version conflicts"**
```bash
# Check for conflicts
pip check

# Start fresh if needed
pip freeze > requirements_backup.txt
pip uninstall -y -r requirements_backup.txt
pip install -r requirements.txt
```

#### Platform-Specific Notes

**macOS with Apple Silicon**
```bash
# May need specific numpy installation
conda install numpy  # Often more reliable than pip
```

**Windows**
```bash
# May need Visual C++ build tools for some packages
# Download from Microsoft if compilation errors occur
```

**Linux (older distributions)**
```bash
# May need to update pip first
pip install --upgrade pip setuptools wheel
```

## Data Integration and Dependencies

### Phase2 Data Requirements
Phase3 requires specific files from phase2 output:

#### Required Directory Structure
```
phase2-output-directory/
├── snapshots/
│   ├── year{N}_snapshot.json[.gz]     # Cell snapshots (e.g., year30, year50)
│   └── metadata.json                  # Snapshot extraction metadata
└── individuals/
    ├── mutant/
    │   └── individual_*.json[.gz]      # Mutant population files
    ├── control1/
    │   └── individual_*.json[.gz]      # Control1 population files
    ├── control2/
    │   └── individual_*.json[.gz]      # Control2 population files
    └── mixing_metadata.json            # Population mixing configuration
```

#### File Content Requirements
- **Snapshots**: Must contain extracted cells with JSD and methylation data
- **Individuals**: Must be PetriDish objects with complete history tracking
- **Metadata**: Must include gene rate groups, snapshot years, and mixing parameters

### Phase1 Data Requirements
Phase3 also requires the original phase1 simulation for timeline analysis:

- **File**: `simulation.json[.gz]` from phase1 output
- **Content**: Complete simulation history with parameters and cell evolution
- **Purpose**: Generate timeline plots and validate consistency with phase2 data

### Data Validation
Phase3 performs extensive validation:

1. **Directory Structure**: Verifies all required directories and files exist
2. **File Format**: Validates JSON structure and compression formats
3. **Parameter Consistency**: Ensures phase1 and phase2 parameters match
4. **Gene Rate Groups**: Validates gene rate group consistency across phases
5. **Data Integrity**: Checks for complete individual histories and valid JSD values

### Read-Only Data Access
**Important**: Phase3 is completely read-only:
- Never modifies phase1 or phase2 data
- Only reads and analyzes existing datasets
- Creates new analysis outputs in separate directory structure
- Safe to run multiple times on same data
- Multiple analyses can run concurrently on same phase2 data

## Performance Optimization

### Analysis Speed Optimization

#### Quick Analysis Settings
```bash
# Fastest analysis (development/testing)
python run_analysis.py ... --bins 50 --max-gene-plots 3

# Balanced analysis (standard workflow)
python run_analysis.py ... --bins 150 --max-gene-plots 10

# High-quality analysis (publication)
python run_analysis.py ... --bins 300  # All genes
```

#### Memory Management
- **Large datasets**: Use fewer bins to reduce memory usage
- **Many individuals**: Limit gene plots to essential genes
- **Batch processing**: Process runs sequentially to avoid memory conflicts
- **System monitoring**: Monitor RAM usage during analysis

### I/O Optimization

#### File System Performance
- **SSD storage**: Use SSD for phase2 data and analysis output
- **Local processing**: Copy data locally for network-attached storage
- **Compression**: Ensure phase2 data is compressed (`.json.gz`)
- **Cleanup**: Remove old analysis outputs to free disk space

#### Data Loading Strategies
```bash
# Pre-validate data before expensive analysis
python -c "from core.data_loader import validate_phase2_data; validate_phase2_data('path')"

# Check available memory before large analysis
free -h  # Linux
top -l 1 | head -n 10 | grep PhysMem  # macOS
```

### Parallel Processing

#### System-Level Parallelization
```bash
# Use GNU parallel for multiple phase2 runs
ls -d ../phase2/data/*/ | parallel python run_analysis.py --phase2-dir {} --simulation ../phase1/data/*/simulation.json.gz

# Limit parallel jobs based on CPU cores
parallel -j 4 python run_analysis.py --phase2-dir {} --simulation ../phase1/data/*/simulation.json.gz :::: phase2_dirs.txt
```

#### Process Monitoring
```bash
# Monitor analysis progress
watch -n 5 'ps aux | grep run_analysis'

# Monitor I/O usage
iostat -x 1  # Linux
```

### Resource Requirements

#### Typical Resource Usage
| Dataset Size | RAM Usage | Analysis Time | Storage Output |
|--------------|-----------|---------------|----------------|
| Small (64 cells, 50 years) | 0.5 GB | 30 seconds | 10 MB |
| Medium (512 cells, 100 years) | 2 GB | 2 minutes | 50 MB |
| Large (8192 cells, 200 years) | 8 GB | 10 minutes | 200 MB |
| Very Large (32768 cells, 500 years) | 32 GB | 45 minutes | 1 GB |

#### Scaling Guidelines
- **RAM**: ~4x the compressed phase2 data size
- **CPU**: Single-threaded, benefits from high clock speed
- **Storage**: ~5-10x phase2 data size for output
- **Time**: Roughly linear with number of cells × time points

### Performance Monitoring

#### Built-in Timing
Phase3 includes automatic timing for major stages:
```
Stage 1: Snapshot distributions... 15.3s
Stage 2: Population comparisons... 8.7s  
Stage 3: Individual trajectories... 45.2s
Stage 4: Timeline analysis... 12.1s
Total analysis time: 81.3s
```

#### Custom Benchmarking
```bash
# Time complete analysis
time python run_analysis.py --phase2-dir ... --simulation ...

# Profile memory usage (requires psutil)
python -c "import psutil; print(f'Available RAM: {psutil.virtual_memory().available / 1e9:.1f} GB')"
```

### Troubleshooting Performance Issues

#### Common Bottlenecks
1. **Memory exhaustion**: Reduce bins or limit gene plots
2. **Slow I/O**: Move data to faster storage
3. **CPU saturation**: Avoid concurrent analyses
4. **Large outputs**: Use custom output directory on fast storage

#### Performance Diagnostics
```bash
# Check disk space
df -h

# Monitor real-time resource usage
htop  # or top

# Check for memory leaks in long-running batch jobs
while true; do ps aux | grep python | awk '{print $6}' | sort -n | tail -1; sleep 10; done
```

## Troubleshooting and Diagnostics

### Common Issues and Solutions

#### File and Path Issues

**"Phase2 directory not found"**
```bash
# Verify phase2 directory exists and has correct structure
ls -la ../phase2/data/
ls -la {phase2_dir}/snapshots/
ls -la {phase2_dir}/individuals/

# Check for typos in path
find ../phase2/data -name "*snap*" -type d
```

**"Simulation file not found"**
```bash
# Find available simulation files
find ../phase1/data -name "simulation.json*" -type f

# Verify file is readable
file {simulation_path}
zcat {simulation_path} | head -n 5  # For .gz files
```

**"Permission denied"**
```bash
# Check file permissions
ls -l {phase2_dir}/snapshots/
ls -l {simulation_path}

# Fix permissions if needed
chmod 644 {files}
chmod 755 {directories}
```

#### Memory and Performance Issues

**"Memory error" or "Out of memory"**
```bash
# Check available RAM
free -h  # Linux
top -l 1 | head -n 10 | grep PhysMem  # macOS

# Reduce memory usage
python run_analysis.py ... --bins 50 --max-gene-plots 5
```

**"Analysis taking too long"**
```bash
# Monitor progress
tail -f analysis.log  # If logging enabled
ps aux | grep run_analysis

# Use faster settings
python run_analysis.py ... --config configs/quick_analysis.yaml
```

#### Data Integrity Issues

**"Gene rate groups inconsistent"**
```bash
# Check phase1 simulation parameters
zcat {simulation_path} | jq '.parameters.gene_rate_groups'

# Check phase2 metadata
cat {phase2_dir}/snapshots/metadata.json | jq '.gene_rate_groups'
```

**"Missing individuals or snapshots"**
```bash
# Count individuals
find {phase2_dir}/individuals -name "individual_*.json*" | wc -l

# Check snapshot years
ls {phase2_dir}/snapshots/year*_snapshot.json*
```

### Comprehensive Validation

Phase3 performs extensive validation before analysis:

#### Pre-Analysis Validation
```python
# Manual validation (useful for debugging)
from core.data_loader import validate_phase2_data, validate_simulation_file

# Validate phase2 data structure
try:
    validate_phase2_data("/path/to/phase2/dir")
    print("✓ Phase2 data valid")
except Exception as e:
    print(f"❌ Phase2 validation failed: {e}")

# Validate simulation file
try:
    validate_simulation_file("/path/to/simulation.json.gz")
    print("✓ Simulation file valid")
except Exception as e:
    print(f"❌ Simulation validation failed: {e}")
```

#### Validation Checklist
- [ ] **Directory structure**: All required subdirectories exist
- [ ] **File formats**: JSON files are valid and readable
- [ ] **Compression**: Compressed files decompress correctly
- [ ] **Metadata consistency**: Parameters match between phase1 and phase2
- [ ] **Data completeness**: All expected individuals and snapshots present
- [ ] **Gene rate groups**: Consistent across all data sources
- [ ] **Time points**: Snapshot years are valid and consistent
- [ ] **File permissions**: All files are readable

### Diagnostic Tools

#### Built-in Diagnostics
```bash
# Dry-run mode (validate without analysis)
python run_analysis.py ... --validate-only  # If implemented

# Verbose output for debugging
python run_analysis.py ... --verbose
```

#### Manual Inspection Tools
```bash
# Check JSON file structure
jq '.' {json_file} | head -20

# Verify compression
file {compressed_file}
gzip -t {gz_file}  # Test compression integrity

# Count data elements
jq '.history | keys | length' {simulation_file}  # Time points
jq '.cells | length' {snapshot_file}  # Cells in snapshot
```

#### System Diagnostics
```bash
# Check system resources
df -h  # Disk space
free -h  # Memory usage
nproc  # CPU cores

# Monitor during analysis
watch -n 2 'ps aux | grep python; free -h'
```

### Error Recovery

#### Partial Analysis Recovery
If analysis fails partway through:

1. **Check partial outputs**: Some stages may have completed
2. **Identify failure point**: Look at error messages and partial results
3. **Resume if possible**: Some operations can be repeated safely
4. **Clean restart**: Remove partial outputs and restart

#### Data Corruption Recovery
```bash
# Verify data integrity
find {phase2_dir} -name "*.json.gz" -exec gzip -t {} \;

# Re-extract from phase2 if needed
cd ../phase2
python extract_snapshots.py --force-reload ...
```

### Getting Help

#### Debug Information Collection
When reporting issues, collect:

```bash
# System information
python --version
pip list | grep -E "(numpy|scipy|plotly|yaml)"
uname -a

# Data structure information
ls -la {phase2_dir}
du -sh {phase2_dir}/*
jq '.parameters' {simulation_file}

# Error logs
# Run with verbose output and save logs
python run_analysis.py ... 2>&1 | tee analysis_debug.log
```

#### Support Resources
1. **Documentation**: Check CLAUDE.md and README files
2. **Test cases**: Run on smaller datasets first
3. **Issue tracking**: Document reproducible error cases
4. **Version information**: Ensure all components are compatible

## Future Enhancements and Roadmap

### Near-Term Improvements
- **Interactive plots**: HTML output with plotly for web viewing
- **Statistical exports**: CSV/Excel reports for external analysis
- **Enhanced batch processing**: Built-in parallel processing
- **Quality metrics**: Automated analysis quality assessment
- **Plot customization**: User-defined color schemes and layouts

### Medium-Term Features
- **Automated comparisons**: Cross-run comparative analysis
- **Machine learning integration**: Pattern recognition in methylation data
- **Custom analysis pipelines**: User-defined analysis workflows
- **Database integration**: Store and query analysis results
- **Real-time monitoring**: Progress tracking for long-running analyses

### Long-Term Vision
- **Web dashboard**: Interactive results browser
- **Cloud integration**: Distributed analysis across clusters
- **GPU acceleration**: CUDA-based analysis for large datasets
- **API interface**: RESTful API for programmatic access
- **Collaborative features**: Shared analysis environments

### Performance Optimization
- **Memory efficiency**: Streaming analysis for large datasets
- **Parallel processing**: Multi-core utilization
- **Caching strategies**: Intermediate result caching
- **Incremental analysis**: Update analysis with new data
- **Format optimization**: HDF5 support for large numerical data

### Integration Improvements
- **External tool connectors**: R, MATLAB, ImageJ integration
- **Workflow engines**: Nextflow, Snakemake compatibility
- **Container support**: Docker/Singularity deployment
- **Cloud storage**: S3, GCS data access
- **Version control**: Analysis provenance tracking

### Research Features
- **Advanced statistics**: Bayesian analysis, mixed-effects models
- **Pathway analysis**: Gene set enrichment analysis
- **Spatial analysis**: Tissue-level methylation patterns
- **Temporal modeling**: Time-series analysis methods
- **Comparative genomics**: Cross-species analysis support

## Examples and Use Cases

### Quick Testing and Development
```bash
# Fast analysis for development and testing
python run_analysis.py \
    --phase2-dir ../phase2/data/latest_run/ \
    --simulation ../phase1/data/latest/simulation.json.gz \
    --config configs/quick_analysis.yaml

# Custom quick analysis
python run_analysis.py \
    --phase2-dir ../phase2/data/test_run/ \
    --simulation ../phase1/data/test/simulation.json.gz \
    --bins 50 --max-gene-plots 3
```

### Production Analysis
```bash
# High-resolution analysis for publication
python run_analysis.py \
    --phase2-dir ../phase2/data/production_run/ \
    --simulation ../phase1/data/production/simulation.json.gz \
    --bins 300 \
    --output-dir publication_analysis_$(date +%Y%m%d_%H%M%S)

# Standard production analysis
python run_analysis.py \
    --phase2-dir ../phase2/data/experiment_001/ \
    --simulation ../phase1/data/exp001/simulation.json.gz \
    --config configs/default.yaml
```

### Comparative Studies
```bash
# Analyze multiple experimental conditions
for condition in low_rate medium_rate high_rate; do
    python run_analysis.py \
        --phase2-dir ../phase2/data/gene_rates_*${condition}*/ \
        --simulation ../phase1/data/gene_rates_*${condition}*/simulation.json.gz \
        --output-dir comparative_${condition}_$(date +%Y%m%d)
done
```

### Time Series Analysis
```bash
# Analyze different time points
for timepoint in snap20to40 snap30to50 snap40to60; do
    python run_analysis.py \
        --phase2-dir ../phase2/data/*${timepoint}*/ \
        --simulation ../phase1/data/*/simulation.json.gz \
        --output-dir timeseries_${timepoint}
done
```

### Parameter Sensitivity Analysis
```bash
# Test different analysis parameters on same data
for bins in 100 200 400; do
    python run_analysis.py \
        --phase2-dir ../phase2/data/reference_run/ \
        --simulation ../phase1/data/reference/simulation.json.gz \
        --bins $bins \
        --output-dir sensitivity_bins${bins}
done
```

### Large-Scale Batch Processing
```bash
# Process all available phase2 runs
#!/bin/bash
log_file="batch_analysis_$(date +%Y%m%d_%H%M%S).log"
echo "Starting batch analysis at $(date)" > $log_file

count=0
success=0

for phase2_dir in ../phase2/data/*/; do
    count=$((count + 1))
    echo "[$count] Processing: $phase2_dir" | tee -a $log_file
    
    if python run_analysis.py \
        --phase2-dir "$phase2_dir" \
        --simulation "../phase1/data/*/simulation.json.gz" \
        --config configs/default.yaml >> $log_file 2>&1; then
        success=$((success + 1))
        echo "[$count] ✓ Success" | tee -a $log_file
    else
        echo "[$count] ❌ Failed" | tee -a $log_file
    fi
done

echo "Completed: $success/$count analyses" | tee -a $log_file
echo "Log saved to: $log_file"
```

### Integration with External Tools
```bash
# Export results for external analysis
python run_analysis.py \
    --phase2-dir ../phase2/data/export_study/ \
    --simulation ../phase1/data/export/simulation.json.gz \
    --output-dir external_analysis

# Convert plots to different formats if needed
cd data/external_analysis/results
find . -name "*.png" -exec convert {} {}.pdf \;
```