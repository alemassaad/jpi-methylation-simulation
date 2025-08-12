# Cleaning Guide for Step23-Prime

## Selective Cleaning Options

### 1. Clean ONLY Individuals (Keep Snapshots)
```bash
# This is the most common case - keep expensive snapshots, redo individuals
rm -rf data/rate_0.005000/individuals/

# Or clean specific groups only:
rm -rf data/rate_0.005000/individuals/mutant/    # Clean only mutants
rm -rf data/rate_0.005000/individuals/control1/  # Clean only control1
rm -rf data/rate_0.005000/individuals/control2/  # Clean only control2
```

### 2. Clean Everything EXCEPT Snapshots
```bash
# Remove individuals, plots, and results but keep snapshots
rm -rf data/rate_0.005000/individuals/
rm -rf data/rate_0.005000/plots/
rm -rf data/rate_0.005000/results/
# Keep: data/rate_0.005000/snapshots/
```

### 3. Clean Specific Individuals
```bash
# Remove specific individuals if corrupted
rm data/rate_0.005000/individuals/mutant/individual_00.json.gz
rm data/rate_0.005000/individuals/mutant/individual_01.json.gz
# Pipeline will recreate only missing ones
```

### 4. Full Clean (Complete Restart)
```bash
# Remove everything including cached snapshots
rm -rf data/rate_0.005000/
```

### 5. Clean Multiple Rates
```bash
# If you have multiple rates
rm -rf data/rate_*/individuals/  # Clean all individuals
# Keeps all snapshots for all rates
```

## Using the Clean Script
```bash
# Use the provided script
./clean_individuals.sh 0.005  # Cleans individuals for rate 0.005
./clean_individuals.sh 0.007  # Cleans individuals for rate 0.007
```

## What Gets Cached and Why

### Expensive to Recreate (Keep These):
- **snapshots/year50_snapshot.json.gz** - Extracting 10,000 cells from year 50
- **snapshots/year60_snapshot.json.gz** - Extracting 10,000 cells from year 60

### Cheap to Recreate (OK to Clean):
- **plots/** - Quick visualization generation
- **results/** - Quick statistics calculation

### Moderate Cost (Main Processing):
- **individuals/** - The actual experimental pipeline

## Pipeline Skip Logic

The pipeline checks in this order:
1. **Snapshots exist?** → Skip extraction (saves ~1-2 minutes)
2. **Individuals exist?** → Check cell counts
   - If correct count → Skip that stage
   - If wrong count → Recreate
3. **Plots exist?** → Skip if file exists
4. **Results exist?** → Skip if file exists

## Example Workflow

```bash
# First run - creates everything
python run_pipeline.py --rate 0.005 --simulation ../step1/data/*.json.gz

# Something went wrong with individuals? Clean and retry:
rm -rf data/rate_0.005000/individuals/
python run_pipeline.py --rate 0.005 --simulation ../step1/data/*.json.gz
# This time it will skip snapshot extraction (fast!)

# Want to try different parameters? Clean individuals only:
./clean_individuals.sh 0.005
python run_pipeline.py --rate 0.005 --simulation ../step1/data/*.json.gz \
    --n-quantiles 4 --cells-per-quantile 2  # Different parameters

# Corrupted control2? Clean just that:
rm -rf data/rate_0.005000/individuals/control2/
python run_pipeline.py --rate 0.005 --simulation ../step1/data/*.json.gz
# Will only recreate control2
```

## Pro Tips

1. **Always keep snapshots** unless you suspect they're corrupted
2. **Check file sizes** before cleaning:
   ```bash
   du -sh data/rate_0.005000/*
   # Snapshots are usually 2-3MB each
   # Individuals are ~100MB per group (30 files)
   ```

3. **Partial runs are OK** - The pipeline can resume from any stage

4. **Test with small parameters first**:
   ```bash
   # Quick test run
   python run_pipeline.py --n-quantiles 2 --cells-per-quantile 1 --growth-years 2
   # Then clean and run full:
   ./clean_individuals.sh 0.005
   python run_pipeline.py --n-quantiles 10 --cells-per-quantile 3 --growth-years 10
   ```

5. **Monitor disk space** - Full run creates ~300MB of data per rate