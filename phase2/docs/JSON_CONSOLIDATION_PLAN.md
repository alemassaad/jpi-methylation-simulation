# JSON File Consolidation Plan

## Current State (4 separate files in results/)

1. **`statistics.json`**
   - Statistical comparison of cell JSDs between batches
   - Mean, std, median, min, max for each batch
   - T-test p-values between batch pairs

2. **`jsd_distributions.json`**
   - Raw arrays of mean cell JSD per individual
   - One array per batch (mutant, control1, control2)
   - Each value is the mean JSD across all cells in that individual

3. **`mixing_statistics.json`**
   - Statistics about the uniform mixing pool (only when --uniform-mixing used)
   - Pool size, individual sizes before/after mixing
   - Created in Stage 6 of pipeline

4. **`pipeline_metadata.json`**
   - Complete pipeline configuration and parameters
   - All args, paths, timestamps, versions
   - Full record of how the pipeline was run

## Issues with Current Naming

- `jsd_distributions.json` is ambiguous - doesn't specify cell vs gene JSD
- `statistics.json` is too generic - what statistics?
- Files are scattered without clear relationships
- Hard to add gene JSD data with parallel structure

## Proposed Solution: Option 4 (Recommended)

### New Structure

**1. `cell_jsd_analysis.json`** (combines statistics.json + jsd_distributions.json)
```json
{
  "summary_statistics": {
    "mutant": {
      "mean": 0.234,
      "std": 0.045,
      "median": 0.229,
      "min": 0.156,
      "max": 0.341,
      "n_individuals": 7
    },
    "control1": { ... },
    "control2": { ... }
  },
  "statistical_tests": {
    "mutant_vs_control1": {
      "t_statistic": 2.34,
      "p_value": 0.023
    },
    "mutant_vs_control2": { ... },
    "control1_vs_control2": { ... }
  },
  "individual_means": {  // Renamed from "distributions" for clarity
    "mutant": [0.234, 0.189, 0.267, ...],  // One value per individual
    "control1": [0.198, 0.211, ...],
    "control2": [0.165, 0.171, ...]
  }
}
```

**2. `pipeline_metadata.json`** (unchanged)
- Keep separate as it serves different purpose
- Contains configuration, not results

**3. `mixing_statistics.json`** (unchanged, conditional)
- Only created when --uniform-mixing is used
- Specific to mixing process, not analysis results

### Future Gene JSD Structure (parallel structure)

**`gene_jsd_analysis.json`** (when implemented)
```json
{
  "summary_statistics": {
    "by_gene": {
      "gene_0": {"mean": ..., "std": ...},
      "gene_1": { ... }
    },
    "overall": {
      "mean_across_genes": ...,
      "std_across_genes": ...
    }
  },
  "statistical_tests": { ... },
  "gene_distributions": {
    "gene_0": {
      "mutant": [...],
      "control1": [...],
      "control2": [...]
    }
  }
}
```

## Implementation Steps

1. **Modify `pipeline_analysis.py`**:
   - Change `analyze_populations_from_dishes()` to create single `cell_jsd_analysis.json`
   - Remove separate statistics.json and jsd_distributions.json
   - Use clear field names (individual_means vs distributions)

2. **Update plot filename**:
   - Already done: `jsd_comparison.png` → `cell_jsd_comparison.png`

3. **Update documentation**:
   - README.md file structure
   - Code comments
   - Any analysis scripts that load these files

4. **Backward compatibility**:
   - Consider adding a small script to convert old format to new
   - Or document the breaking change clearly

## Benefits

1. **Clarity**: File names explicitly state "cell_jsd"
2. **Organization**: Related data grouped together
3. **Extensibility**: Easy to add gene_jsd_analysis.json with parallel structure
4. **Efficiency**: One file read for all cell JSD analysis data
5. **Semantics**: "individual_means" clearer than "distributions"

## Alternative Options Considered

### Option 1: Single mega-file
- Too large, mixes different data types
- Hard to load just what you need

### Option 2: Just rename existing files
- `cell_jsd_statistics.json`, `cell_jsd_distributions.json`
- Better than current but still scattered

### Option 3: Group by data type
- `analysis_results.json`, `raw_distributions.json`
- Less clear what's cell vs gene JSD

## Decision

**Go with Option 4**: Create `cell_jsd_analysis.json` that consolidates the cell JSD statistics and distributions with clear, semantic field names.

## TODO After Compacting

- [ ] Implement the consolidation in `pipeline_analysis.py`
- [ ] Update any scripts that read the old files
- [ ] Test with existing pipelines
- [ ] Update documentation
- [ ] Consider migration script for old data