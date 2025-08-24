# Complete Simulation Guide - All Flags and Options

## Phase 1: Running Simulations

### Basic Command Structure
```bash
cd phase1
python run_simulation.py [OPTIONS]
```

### ALL Available Flags

#### 1. METHYLATION RATE FLAGS (Choose ONE)

##### `--rate` or `-r` (Uniform methylation)
- **Purpose**: Sets uniform methylation rate for ALL genes
- **Format**: Decimal number (e.g., 0.005 = 0.5% per year)
- **Default**: 0.005 if flag used without value
- **Examples**:
```bash
# Explicit rate
python run_simulation.py --rate 0.005

# Using short form
python run_simulation.py -r 0.005

# Using default (0.005)
python run_simulation.py -r
```

##### `--gene-rate-groups` (Variable methylation) 
- **Purpose**: Different gene groups methylate at different rates
- **Format**: `"n1:rate1,n2:rate2,n3:rate3,..."`
  - n = number of genes in group
  - rate = methylation rate for that group
  - Total genes MUST equal sites/gene_size (default: 200)
- **Examples**:
```bash
# 4 groups of 50 genes each
python run_simulation.py --gene-rate-groups "50:0.002,50:0.004,50:0.006,50:0.008"

# 2 groups: slow and fast
python run_simulation.py --gene-rate-groups "100:0.001,100:0.01"

# 10 groups with gradient
python run_simulation.py --gene-rate-groups "20:0.001,20:0.002,20:0.003,20:0.004,20:0.005,20:0.006,20:0.007,20:0.008,20:0.009,20:0.010"
```

#### 2. SIMULATION PARAMETERS

##### `--years` or `-t`
- **Purpose**: Total number of years to simulate
- **Default**: 100
- **Range**: Any positive integer
- **Examples**:
```bash
# Quick test
python run_simulation.py --rate 0.005 --years 10

# Standard run
python run_simulation.py --rate 0.005 --years 100

# Long simulation
python run_simulation.py --rate 0.005 --years 200
```

##### `--sites` or `-n`
- **Purpose**: Number of CpG sites per cell
- **Default**: 1000
- **Constraint**: Must be divisible by gene_size
- **Examples**:
```bash
# Small genome
python run_simulation.py --rate 0.005 --sites 500

# Large genome
python run_simulation.py --rate 0.005 --sites 5000

# With gene groups (must adjust counts)
python run_simulation.py --sites 2000 --gene-rate-groups "200:0.004,200:0.006"  # 2000/5 = 400 genes
```

##### `--gene-size` or `-g`
- **Purpose**: Number of CpG sites per gene
- **Default**: 5
- **Constraint**: Must divide sites evenly
- **Examples**:
```bash
# Smaller genes
python run_simulation.py --rate 0.005 --gene-size 2

# Larger genes  
python run_simulation.py --rate 0.005 --gene-size 10

# With 1000 sites and gene-size 10 = 100 genes total
python run_simulation.py --gene-size 10 --gene-rate-groups "50:0.004,50:0.006"
```

##### `--growth-phase`
- **Purpose**: Duration of exponential growth phase (years)
- **Default**: 13 (gives 8192 cells)
- **Range**: 1-20
- **Formula**: Final population = 2^growth_phase
- **Examples**:
```bash
# Small population (2^3 = 8 cells)
python run_simulation.py --rate 0.005 --growth-phase 3

# Medium population (2^7 = 128 cells)
python run_simulation.py --rate 0.005 --growth-phase 7

# Default large population (2^13 = 8192 cells)
python run_simulation.py --rate 0.005 --growth-phase 13

# Maximum (2^20 = 1,048,576 cells) - VERY SLOW!
python run_simulation.py --rate 0.005 --growth-phase 20
```

#### 3. TRACKING AND PERFORMANCE FLAGS

##### `--track-gene-jsd`
- **Purpose**: Track gene-level JSD evolution during simulation
- **Default**: False
- **Output**: Stores gene_jsd_history in simulation file
- **Performance**: Slightly slower
- **Example**:
```bash
python run_simulation.py --rate 0.005 --years 100 --track-gene-jsd
```

##### `--no-jsds`
- **Purpose**: Disable ALL JSD calculations for maximum speed
- **Default**: False (JSDs are calculated)
- **Use case**: Very large simulations where JSD not needed
- **Performance**: Significantly faster
- **Example**:
```bash
python run_simulation.py --rate 0.005 --years 100 --no-jsds
```

#### 4. OUTPUT CONTROL FLAGS

##### `--output` or `-o`
- **Purpose**: Custom output filename
- **Default**: Auto-generated based on parameters
- **Example**:
```bash
python run_simulation.py --rate 0.005 --output my_custom_simulation.json.gz
```

##### `--no-save`
- **Purpose**: Run simulation WITHOUT saving results
- **Default**: False (results are saved)
- **Use case**: Testing or when only need console output
- **Example**:
```bash
python run_simulation.py --rate 0.005 --years 10 --no-save
```

#### 5. REPRODUCIBILITY FLAG

##### `--seed`
- **Purpose**: Random seed for reproducibility
- **Default**: 42
- **Special value**: -1 means no seed (truly random)
- **Examples**:
```bash
# Reproducible run (default)
python run_simulation.py --rate 0.005 --seed 42

# Different seed
python run_simulation.py --rate 0.005 --seed 12345

# Truly random (different results each time)
python run_simulation.py --rate 0.005 --seed -1
```

### COMPLETE EXAMPLE COMMANDS

#### Example 1: Minimal Command
```bash
# Must specify rate method, everything else uses defaults
python run_simulation.py --rate 0.005
# OR
python run_simulation.py --gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055"
```

#### Example 2: Quick Test Run
```bash
python run_simulation.py \
    --rate 0.005 \
    --years 10 \
    --growth-phase 3 \
    --seed 42
```

#### Example 3: Full Custom Parameters
```bash
python run_simulation.py \
    --gene-rate-groups "40:0.002,40:0.004,40:0.006,40:0.008,40:0.010" \
    --sites 1000 \
    --gene-size 5 \
    --years 100 \
    --growth-phase 13 \
    --track-gene-jsd \
    --seed 12345 \
    --output my_simulation.json.gz
```

#### Example 4: Performance-Optimized Large Run
```bash
python run_simulation.py \
    --rate 0.005 \
    --years 200 \
    --growth-phase 15 \
    --no-jsds \
    --seed 42
```

#### Example 5: Biological Scenario - CpG Islands
```bash
# 30% CpG islands (low methylation), 70% non-islands (normal methylation)
python run_simulation.py \
    --gene-rate-groups "60:0.001,140:0.007" \
    --years 100 \
    --growth-phase 13 \
    --track-gene-jsd \
    --seed 42
```

### OUTPUT STRUCTURE

#### Console Output During Run
```
============================================================
STEP1-PRIME SIMULATION
============================================================
Parameters:
  Methylation rate: 0.500%           # OR Gene-specific rates: 4 groups
  CpG sites per cell: 1000
  Gene size: 5
  Growth phase: 13 years
  Target population: 8192 (2^13)
  Max years: 100
  Random seed: 42
============================================================

Year 0: Starting with 1 unmethylated cell

Year 1:
  Growth phase (year 1 of 13)
  Division: 1 → 2 cells
  Methylation applied to 2 cells
  Final count: 2 cells (predictable: 2^1)
  Mean JSD: 0.0045

[... continues for each year ...]

============================================================
Simulation complete!
Final population: 8145 cells
============================================================

Simulation runtime: 45.23 seconds

Results saved to: data/rate_0.00500/grow13-sites1000-years100-seed42-a3f7/simulation.json.gz

Final Statistics:
  Total cells: 8145
  Mean JSD: 0.0234
  Median JSD: 0.0221
  Std Dev JSD: 0.0089
  Min JSD: 0.0098
  Max JSD: 0.0567
  Mean methylation: 4.85%
  Median methylation: 4.70%
```

#### File Output Structure
```
phase1/data/
├── rate_0.00500/                    # For uniform rate
│   └── grow13-sites1000-years100-seed42-XXXX/
│       └── simulation.json.gz       # Complete history
│
└── gene_rates_50x0.00400_50x0.00600.../  # For gene-specific rates
    └── grow13-sites1000-years100-seed42-XXXX/
        └── simulation.json.gz
```

#### Simulation.json.gz Contents
```json
{
  "0": [                              # Year 0
    {
      "cpg_sites": [0,0,0,...],      # Methylation state per site
      "methylation_proportion": 0.0,  # Fraction methylated
      "methylation_distribution": [1.0, 0.0, ...],  # Per-gene distribution
      "cell_jsd": 0.0,                # Cell's JSD score
      "age": 0,                        # Cell age
      "gene_size": 5,                  # Sites per gene
      "rate": 0.005,                   # If uniform rate
      "gene_rate_groups": [[50, 0.004], ...],  # If gene-specific
      "site_rates": [0.004, 0.004, ...]  # Pre-computed rates per site
    }
  ],
  "1": [...],                         # Year 1
  "2": [...],                         # Year 2
  ...
  "100": [...]                        # Year 100
}
```

### ERROR MESSAGES AND SOLUTIONS

#### Error: "Cannot specify both 'rate' and 'gene_rate_groups'"
```bash
# WRONG
python run_simulation.py --rate 0.005 --gene-rate-groups "100:0.004,100:0.006"

# CORRECT - Choose one:
python run_simulation.py --rate 0.005
# OR
python run_simulation.py --gene-rate-groups "100:0.004,100:0.006"
```

#### Error: "Must specify either --rate or --gene-rate-groups"
```bash
# WRONG
python run_simulation.py --years 100

# CORRECT
python run_simulation.py --rate 0.005 --years 100
```

#### Error: "gene-rate-groups specifies X genes, but simulation has Y genes"
```bash
# WRONG (only 150 genes, need 200)
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006"

# CORRECT (200 genes total)
python run_simulation.py --gene-rate-groups "50:0.004,50:0.005,50:0.006,50:0.007"
```

#### Error: "Number of sites (X) must be divisible by gene size (Y)"
```bash
# WRONG (1000 not divisible by 7)
python run_simulation.py --rate 0.005 --sites 1000 --gene-size 7

# CORRECT
python run_simulation.py --rate 0.005 --sites 1001 --gene-size 7
```

#### Error: "growth-phase must be between 1 and 20"
```bash
# WRONG
python run_simulation.py --rate 0.005 --growth-phase 25

# CORRECT
python run_simulation.py --rate 0.005 --growth-phase 15
```

### PERFORMANCE TIPS

1. **Quick Tests**: Use `--years 10 --growth-phase 3` for testing
2. **Speed**: Use `--no-jsds` for large simulations if JSD not needed
3. **Memory**: Growth phase >15 creates >32,768 cells (uses lots of RAM)
4. **Disk Space**: Long simulations create large files (100 years ≈ 10-100 MB compressed)

### BIOLOGICAL INTERPRETATIONS

| Growth Phase | Final Cells | Biological Context |
|--------------|-------------|-------------------|
| 1 | 2 | Minimal tissue |
| 3 | 8 | Small organoid |
| 7 | 128 | Tissue sample |
| 10 | 1,024 | Small organ region |
| 13 | 8,192 | Standard tissue (default) |
| 15 | 32,768 | Large tissue sample |
| 20 | 1,048,576 | Full organ (VERY SLOW) |

### NEXT STEPS

After running simulations, analyze with Phase 2:
```bash
cd ../phase2
python run_pipeline.py \
    --simulation ../phase1/data/*/grow*/simulation.json.gz \
    --first-snapshot 50 \
    --second-snapshot 60 \
    --seed 42
```