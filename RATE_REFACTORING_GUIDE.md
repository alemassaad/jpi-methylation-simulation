# Rate Refactoring Guide: Phase 1 → Phase 2

## Overview
This document details the comprehensive refactoring performed on Phase 1 to unify rate handling through `gene_rate_groups`, and provides an exhaustive guide for applying the same refactoring to Phase 2.

## Phase 1 Refactoring Summary

### Core Concept
We eliminated dual code paths by making `gene_rate_groups` the single internal representation for methylation rates. The `rate` parameter is now just syntactic sugar that gets immediately converted to `gene_rate_groups` format.

### Key Changes Made

#### 1. **Added Conversion Utility** (`phase1/cell.py`)
```python
def rate_to_gene_rate_groups(rate: float, n: int, gene_size: int) -> List[Tuple[int, float]]:
    """Convert uniform methylation rate to gene_rate_groups format."""
    n_genes = n // gene_size
    return [(n_genes, rate)]
```

#### 2. **Cell Class Refactoring**

##### Cell.__init__ Changes:
- **Before**: Stored both `self.rate` and `self.gene_rate_groups`, checked which was provided
- **After**: Only stores `self.gene_rate_groups`, converts `rate` to gene_rate_groups format if provided
- **Key code**:
```python
# Convert rate to gene_rate_groups if provided
if rate is not None:
    if gene_rate_groups is not None:
        raise ValueError("Cannot specify both 'rate' and 'gene_rate_groups'")
    gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
```

##### Cell._build_site_rates() Simplification:
- **Before**: Two branches - one for uniform rate, one for gene_rate_groups
- **After**: Single pathway using only gene_rate_groups
- **Removed**: All `if self.rate is not None` checks

##### Cell.from_dict() Updates:
- Now checks for `gene_rate_groups` in the data dictionary
- Falls back to parameters for backward compatibility
- Converts lists to tuples as needed

##### Added Backward Compatibility Property:
```python
@property
def rate(self) -> Optional[float]:
    """Get uniform rate if applicable (backward compatibility)."""
    rates = set(rate for _, rate in self.gene_rate_groups)
    if len(rates) == 1:
        return rates.pop()
    return None
```

#### 3. **PetriDish Class Refactoring**

##### PetriDish.__init__ Changes:
- **Before**: Stored both `self.rate` and `self.gene_rate_groups`
- **After**: Only stores `self.gene_rate_groups`
- Converts rate to gene_rate_groups at initialization
- Uses default rate if neither specified

##### Removed Dual Paths:
- **Before**: Two branches for creating initial cell (uniform vs gene-specific)
- **After**: Single path using gene_rate_groups

##### Updated Validation:
- Renamed `_validate_cell_rate_consistency()` → `validate_cell_consistency()` (made public)
- Added static method `validate_cells_compatible()` for pre-validation
- Simplified to only check gene_rate_groups equality

##### Display Changes:
- Always shows rates as gene groups format in simulation output
- Example: "Gene-specific rates: 1 group(s)" even for uniform rates

#### 4. **JSON Serialization Changes**

##### Parameters Section:
- **Before**: Stored both `'rate'` and `'gene_rate_groups'`
- **After**: Only stores `'gene_rate_groups'`

##### Cell Serialization Enhanced:
```python
# Now includes in each cell:
{
    'methylated': [...],
    'cell_JSD': 0.123,
    'age': 50,
    'gene_rate_groups': [(200, 0.005)],  # NEW
    'methylation_proportion': 0.5,        # NEW
    'n_methylated': 500                   # NEW
}
```

##### Directory Naming:
- **Before**: `data/rate_0.00500/grow13-sites1000-years100-seed42-XXXX/`
- **After**: `data/gene_rates_200x0.00500/size8192-sites1000-genesize5-years100-seed42-YYYYMMDD-HHMMSS/`
- Always uses gene_rates format, even for uniform rates

#### 5. **run_simulation.py Changes**

##### Early Conversion:
```python
# Parse gene rate groups if specified as string
if gene_rate_groups_str:
    gene_rate_groups = parse_gene_rate_groups(gene_rate_groups_str, n, gene_size)
# Convert rate to gene_rate_groups if specified
elif rate is not None:
    gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
else:
    gene_rate_groups = None  # Will use default in PetriDish
```

##### Single PetriDish Creation:
- **Before**: Two branches (rate vs gene_rate_groups)
- **After**: Always passes gene_rate_groups

#### 6. **Config Handling**
- CLI `--rate` overrides config gene_rate_groups (clears it)
- CLI `--gene-rate-groups` overrides config rate (clears it)
- Validation ensures both aren't set simultaneously

## Phase 2 Refactoring Requirements

### Files to Modify

#### 1. **phase2/core/pipeline_utils.py**

##### dict_to_cell() Function:
- **Current Issue**: Still checks for 'rate' in cell dictionaries (lines 58-73)
- **Required Change**: Remove backward compatibility code for old format
- **Update to**:
```python
def dict_to_cell(cell_dict: Dict[str, Any]) -> Cell:
    # Only use gene_rate_groups from cell data
    gene_rate_groups = cell_dict.get('gene_rate_groups')
    if gene_rate_groups and isinstance(gene_rate_groups[0], list):
        gene_rate_groups = [tuple(group) for group in gene_rate_groups]
    
    n = len(cell_dict.get('methylated', cell_dict.get('cpg_sites', [])))
    gene_size = cell_dict.get('gene_size', GENE_SIZE)
    
    # Create cell with gene_rate_groups
    cell = Cell(
        n=n,
        gene_rate_groups=gene_rate_groups,
        gene_size=gene_size
    )
    # ... rest of attribute setting
```

##### load_snapshot_as_cells() Function:
- **Current**: Extracts both rate and gene_rate_groups from parameters
- **Update to**: Only use gene_rate_groups
```python
# Extract rate configuration from parameters
gene_rate_groups = params.get('gene_rate_groups')
if not gene_rate_groups:
    # Backward compatibility: convert rate if found
    rate = params.get('rate')
    if rate:
        gene_rate_groups = rate_to_gene_rate_groups(rate, n, gene_size)
```

##### save_snapshot_cells() Function:
- Ensure it only saves gene_rate_groups in parameters
- Remove any rate-related code

##### save_individual() Function:
- Update metadata to only include gene_rate_groups
- Remove any 'rate' field from metadata

#### 2. **phase2/run_pipeline.py**

##### Import Changes:
```python
from phase1.cell import PetriDish, Cell, rate_to_gene_rate_groups
```

##### Rate Configuration Handling:
- Convert rate to gene_rate_groups early in the pipeline
- Similar to phase1/run_simulation.py approach

##### Remove Dual Paths:
- Search for any `if rate is not None` checks
- Replace with single gene_rate_groups pathway

##### Stage-Specific Updates:

**Stage 3 (Create Initial Individuals)**:
- Ensure cells are created with gene_rate_groups only
- Validate all cells have same gene_rate_groups

**Stage 4 (Grow Individuals)**:
- Pass gene_rate_groups when creating new PetriDish instances
- Remove any rate parameter passing

**Stage 6 (Mix Populations)**:
- Use `PetriDish.validate_cells_compatible()` before mixing
- Ensure mixed cells have consistent gene_rate_groups

**Stage 7 (Create Control2)**:
- Create with gene_rate_groups from snapshot

#### 3. **phase2/core/pipeline_analysis.py**

##### Visualization Functions:
- Update any rate-related labels to show gene_rate_groups format
- Ensure plots show actual gene group structure

#### 4. **phase2/path_utils.py**

##### parse_path_parameters() Function:
- Update to parse new directory format with gene_rates
- Handle size instead of grow
- Parse genesize parameter

##### generate_output_path() Function:
- Generate paths using new format
- Always use gene_rates prefix

#### 5. **phase2/visualization/plot_individuals.py**

##### Any Rate References:
- Update to use gene_rate_groups
- Show gene group structure in plots

### Validation Points to Add in Phase 2

1. **Before Mixing Cells** (Stage 6):
```python
# Validate all cells are compatible before mixing
all_cells = mutant_cells + snapshot_cells
PetriDish.validate_cells_compatible(all_cells)
```

2. **After Loading Snapshots**:
```python
# Verify all loaded cells have consistent rates
PetriDish.validate_cells_compatible(snapshot_cells)
```

3. **When Creating Individuals**:
```python
# Ensure individual's cells are consistent
individual_petri.validate_cell_consistency()
```

### Testing Strategy for Phase 2

1. **Create test_rate_consistency_phase2.py**:
   - Test loading phase1 simulations with new format
   - Test mixing cells from different sources
   - Test control group creation
   - Verify gene_rate_groups consistency throughout pipeline

2. **Update Existing Tests**:
   - Remove any tests expecting 'rate' in metadata
   - Update path expectations for new directory format
   - Verify gene_rate_groups in all saved individuals

3. **Integration Tests**:
   - Full pipeline with uniform rate (converted to gene_rate_groups)
   - Full pipeline with gene-specific rates
   - Verify output format consistency

### Config File Updates for Phase 2

Update phase2 config files to:
- Use gene_rate_groups in examples
- Remove rate from templates
- Update documentation

### Backward Compatibility Considerations

For Phase 2, add temporary compatibility layer:
```python
def load_legacy_simulation(filepath):
    """Load old format simulations and convert to new format."""
    data = load_json(filepath)
    
    # Convert old format
    if 'parameters' in data:
        params = data['parameters']
        if 'rate' in params and 'gene_rate_groups' not in params:
            # Convert rate to gene_rate_groups
            rate = params['rate']
            n = params.get('n', 1000)
            gene_size = params.get('gene_size', 5)
            params['gene_rate_groups'] = rate_to_gene_rate_groups(rate, n, gene_size)
            del params['rate']
    
    return data
```

### Critical Files to Search and Update

Run these searches in phase2 directory:
```bash
# Find all rate references
grep -r "\.rate" phase2/
grep -r "['rate']" phase2/
grep -r '["rate"]' phase2/
grep -r "args.rate" phase2/
grep -r "config.get('rate')" phase2/

# Find dual pathway code
grep -r "if.*rate.*is not None" phase2/
grep -r "if cell.rate" phase2/
grep -r "has_uniform" phase2/
```

### Expected Outcomes

After Phase 2 refactoring:
1. All internal code uses only gene_rate_groups
2. No dual pathways remain
3. Consistent JSON format between phases
4. Rate validation at all boundaries
5. Clear error messages for rate mismatches
6. Backward compatibility for loading old simulations

### Benefits

1. **Consistency**: Same architecture as Phase 1
2. **Simplicity**: Single code path throughout
3. **Validation**: Strong rate consistency checking
4. **Clarity**: Directory names show actual configuration
5. **Maintainability**: Easier to debug and extend

## Implementation Order

1. Start with pipeline_utils.py (core loading/saving)
2. Update run_pipeline.py (main pipeline)
3. Fix visualization and analysis
4. Update tests
5. Run full integration tests
6. Update documentation

## Key Testing Commands

```bash
# Test with uniform rate
python run_pipeline.py --rate 0.005 --simulation ../phase1/data/gene_rates*/simulation.json.gz ...

# Test with gene-specific rates  
python run_pipeline.py --gene-rate-groups "10:0.004,10:0.006" --simulation ...

# Test loading old format (if compatibility layer added)
python run_pipeline.py --simulation ../old_format/simulation.json.gz ...
```

## Success Criteria

- [ ] No references to `self.rate` in Cell objects
- [ ] No references to `petri.rate` in PetriDish objects  
- [ ] All JSON files only contain gene_rate_groups in parameters
- [ ] Directory names always use gene_rates format
- [ ] All cells in same PetriDish have identical gene_rate_groups
- [ ] Tests pass with both uniform and gene-specific rates
- [ ] Clear error messages for rate mismatches