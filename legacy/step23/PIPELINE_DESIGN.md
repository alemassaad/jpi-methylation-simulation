# Step23 Unified Pipeline Design Document

## Overview

Step23 is a unified pipeline that combines the functionality of steps 2 and 3 into a single, efficient process. This document provides detailed implementation specifications for the pipeline.

## Pipeline Architecture

### Input Requirements
- **Primary Input**: Step 1 simulation file (e.g., `simulation_rate_0.005000_m10000_n1000_t100.json.gz`)
- **Parameters**:
  - `rate`: Methylation rate (must match simulation file)
  - `bins`: Number of bins for JSD distribution plots (default: 200)
  - `mix_ratio`: Percentage of year 60 cells in final mixture (default: 80)
  - `n_individuals`: Number of individuals per group (default: 30)
  - `growth_years`: Years to grow individuals (default: 10)
  - `seed`: Random seed for reproducibility (default: 42)

### Output Structure
```
step23/data/rate_0.XXXXXX/
├── snapshots/
│   ├── year50_snapshot.json.gz       # All 10,000 cells at year 50
│   └── year60_snapshot.json.gz       # All 10,000 cells at year 60
├── individuals/
│   ├── mutant/
│   │   └── individual_XX.json.gz     # 30 files, each with 5,120 cells
│   ├── control1/
│   │   └── individual_XX.json.gz     # 30 files, each with 5,120 cells
│   └── control2/
│       └── individual_XX.json.gz     # 30 files, each with 5,120 cells
├── plots/
│   ├── year50_jsd_distribution.png   # JSD distribution at year 50
│   ├── year60_jsd_distribution.png   # JSD distribution at year 60
│   ├── jsd_comparison.png            # Final comparison plot
│   └── jsd_distributions_boxplot.png # Box plot comparison
└── results/
    ├── statistics.json                # Statistical analysis results
    ├── jsd_distributions.json         # Raw JSD values for all individuals
    └── pipeline_metadata.json         # Pipeline parameters and timing
```

## Detailed Pipeline Stages

### Stage 1: Extract Year 50 Snapshot
**Purpose**: Extract all cells at year 50 from the original simulation.

**Process**:
1. Load simulation file
2. Extract year 50 data (index 50 in the time series)
3. Save as compressed JSON to `snapshots/year50_snapshot.json.gz`

**Output Format**:
```json
[
    {
        "cpg_sites": [0, 1, 0, ...],
        "methylation_proportion": 0.234,
        "methylation_distribution": [0.4, 0.3, 0.2, 0.1, 0.0, 0.0],
        "jsd": 0.1234,
        "age": 50,
        "gene_size": 5,
        "rate": 0.005
    },
    // ... 9,999 more cells
]
```

### Stage 2: Plot JSD Distribution
**Purpose**: Visualize the JSD distribution of year 50 cells.

**Process**:
1. Extract JSD values from all cells
2. Create histogram with specified number of bins
3. Add statistics (mean, median, std)
4. Save plot to `plots/year50_jsd_distribution.png`

### Stage 3: Create Initial Individuals
**Purpose**: Sample cells to create initial single-cell individuals.

**Sampling Strategies**:
- **Mutant**: 
  1. Sort cells by JSD value
  2. Divide into 10 deciles
  3. Randomly sample 3 cells from each decile
  4. Total: 30 cells → 30 individual files
  
- **Control1**:
  1. Uniformly sample 30 cells from entire population
  2. No JSD-based selection
  3. Total: 30 cells → 30 individual files

**File Creation**:
Each sampled cell becomes a single-cell individual file:
```json
[
    {
        "cpg_sites": [0, 1, 0, ...],
        "methylation_proportion": 0.234,
        "methylation_distribution": [0.4, 0.3, 0.2, 0.1, 0.0, 0.0],
        "jsd": 0.1234,
        "age": 50,
        "gene_size": 5,
        "rate": 0.005
    }
]
```

### Stage 4: Grow Individuals
**Purpose**: Simulate 10 years of cell division and aging.

**Growth Process** (for each individual file):
```python
for year in range(1, 11):  # Years 51-60
    # Load current cells
    cells = load_individual(file_path)
    
    # Division phase: each cell divides
    new_cells = []
    for cell in cells:
        # Create two daughter cells (exact copies)
        daughter1 = cell.copy()
        daughter2 = cell.copy()
        new_cells.extend([daughter1, daughter2])
    
    # Aging phase: all cells age by 1 year
    for cell in new_cells:
        cell.age()  # Stochastic methylation
    
    # Save updated population back to same file
    save_individual(new_cells, file_path)
```

**Growth Progression**:
- Year 50: 1 cell (initial)
- Year 51: 2 cells (after 1 division + aging)
- Year 52: 4 cells
- Year 53: 8 cells
- Year 54: 16 cells
- Year 55: 32 cells
- Year 56: 64 cells
- Year 57: 128 cells
- Year 58: 256 cells
- Year 59: 512 cells
- Year 60: 1,024 cells (final)

### Stage 5: Extract Year 60 Snapshot
**Purpose**: Extract all cells at year 60 from the original simulation.

**Process**:
1. Load simulation file
2. Extract year 60 data (index 60 in the time series)
3. Save as compressed JSON to `snapshots/year60_snapshot.json.gz`

### Stage 6: Mix Populations
**Purpose**: Add year 60 cells to grown individuals to create mixed populations.

**Mixing Formula**:
```python
# For default 80-20 ratio:
n_grown = 1024  # Cells after growth
percentage_grown = 20  # 20% of final population
percentage_year60 = 80  # 80% of final population

n_total = n_grown * 100 / percentage_grown  # 5,120 cells
n_to_add = n_total - n_grown  # 4,096 cells from year 60
```

**Process** (for each mutant and control1 individual):
1. Load grown individual (1,024 cells)
2. Randomly sample 4,096 cells from year 60 snapshot
3. Combine both populations
4. Save mixed population (5,120 cells) back to same file

### Stage 7: Create Control2 Individuals
**Purpose**: Create pure year 60 control group.

**Process**:
1. For each of 30 individuals:
   - Randomly sample 5,120 cells from year 60 snapshot
   - Save as new individual file in `control2/` directory

### Stage 8: Analysis and Visualization
**Purpose**: Compare JSD distributions across all three groups.

**Analysis Steps**:
1. For each individual in each group:
   - Load all cells
   - Calculate mean JSD across all cells
   
2. Statistical comparison:
   - Calculate mean, median, std for each group
   - Perform statistical tests (t-test, Mann-Whitney U)
   - Generate summary statistics

3. Visualization:
   - Histogram overlay of three distributions
   - Box plot comparison
   - Statistical significance annotations

**Output Files**:
- `statistics.json`: Summary statistics and test results
- `jsd_distributions.json`: Raw JSD values for all individuals
- `jsd_comparison.png`: Main comparison plot
- `pipeline_metadata.json`: Pipeline configuration and timing

## Implementation Guidelines

### Cell Division Implementation
```python
def divide_cell(cell):
    """Create a copy of the cell for division."""
    # Must create deep copy to avoid reference issues
    new_cell = Cell(n=cell.n, gene_size=cell.gene_size, rate=cell.rate)
    new_cell.cpg_sites = cell.cpg_sites.copy()
    new_cell.age_val = cell.age_val
    # Recalculate derived properties
    new_cell.methylation_proportion = sum(new_cell.cpg_sites) / new_cell.n
    new_cell.methylation_distribution = new_cell.calculate_methylation_distribution()
    new_cell.jsd = new_cell.calculate_jsd()
    return new_cell
```

### Memory Management
- Process individuals one at a time during growth
- Use generators for large cell populations
- Clear unused variables after each stage

### Random Seed Management
```python
def set_stage_seed(base_seed, stage_num):
    """Set reproducible seed for each pipeline stage."""
    random.seed(base_seed + stage_num)
    np.random.seed(base_seed + stage_num)
```

### Progress Tracking
```python
def progress_bar(current, total, stage_name):
    """Display progress for current stage."""
    percentage = (current / total) * 100
    bar_length = 50
    filled = int(bar_length * current / total)
    bar = '█' * filled + '-' * (bar_length - filled)
    print(f'\r{stage_name}: |{bar}| {percentage:.1f}% ({current}/{total})', end='')
```

### Error Handling
- Validate input simulation file exists and contains required years
- Check rate parameter matches simulation file
- Ensure output directories are writable
- Graceful handling of interrupted growth process

## Performance Considerations

### Expected Runtime
- Stage 1 (Extract year 50): ~5 seconds
- Stage 2 (Plot JSD): ~2 seconds
- Stage 3 (Create individuals): ~1 second
- Stage 4 (Grow individuals): ~5-10 minutes (main bottleneck)
- Stage 5 (Extract year 60): ~5 seconds
- Stage 6 (Mix populations): ~30 seconds
- Stage 7 (Create control2): ~10 seconds
- Stage 8 (Analysis): ~5 seconds

**Total**: ~6-12 minutes per rate

### Optimization Opportunities
1. **Parallel growth**: Process multiple individuals simultaneously
2. **Batch operations**: Load/save multiple cells at once
3. **Numpy vectorization**: Use numpy for methylation calculations
4. **Caching**: Cache year 50/60 snapshots if running multiple analyses

## Testing Strategy

### Unit Tests
- Cell division produces correct number of cells
- Decile sampling selects from correct ranges
- Mixing produces correct proportions
- JSD calculations match expected values

### Integration Tests
- Full pipeline runs without errors
- Output files have correct structure
- Reproducible with same seed
- Results match legacy pipeline for same inputs

### Validation Tests
- Growth produces 2^n cells after n years
- All cells age correctly during growth
- Mixed populations have correct total size
- Statistical tests produce meaningful results

## Migration from Legacy Pipeline

### Equivalence Mapping
- Step2 lineages → Step23 grown individuals (before mixing)
- Step3 mixed individuals → Step23 mixed individuals (after stage 6)
- Step3 control → Step23 control1
- New in Step23 → control2 (pure year 60)

### Key Differences
1. **No separate lineage files**: Individuals grow in place
2. **Direct pipeline**: No intermediate manual steps
3. **Additional control group**: Control2 for baseline comparison
4. **Unified configuration**: Single config file for all parameters

## Future Enhancements

### Planned Features
1. **Batch processing**: Run multiple rates in parallel
2. **Resume capability**: Checkpoint after each stage
3. **Custom growth models**: Support different division patterns
4. **Extended analysis**: Additional statistical metrics
5. **Interactive plots**: HTML outputs with plotly

### Extensibility Points
- Custom sampling strategies
- Alternative mixing ratios
- Different growth periods
- Additional visualization types
- Export to different formats

## Appendix: Mathematical Details

### Cell Division and Aging
1. **Division**: Creates exact copy with identical methylation state
2. **Aging**: Each CpG site has probability `rate` of becoming methylated
3. **Independence**: Each cell ages independently after division

### JSD Calculation
```python
def calculate_jsd(distribution, baseline):
    """Calculate Jensen-Shannon Divergence."""
    # Ensure distributions sum to 1
    dist = np.array(distribution) / np.sum(distribution)
    base = np.array(baseline) / np.sum(baseline)
    
    # Calculate midpoint distribution
    midpoint = (dist + base) / 2
    
    # Calculate KL divergences
    kl_dist = np.sum(dist * np.log2(dist / midpoint + 1e-10))
    kl_base = np.sum(base * np.log2(base / midpoint + 1e-10))
    
    # JSD is average of KL divergences
    jsd = (kl_dist + kl_base) / 2
    return jsd
```

### Statistical Tests
- **T-test**: Compare means between groups (assumes normality)
- **Mann-Whitney U**: Non-parametric comparison (no normality assumption)
- **Kolmogorov-Smirnov**: Compare entire distributions
- **Effect size**: Cohen's d for practical significance