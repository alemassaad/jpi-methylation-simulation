# In-Memory Validation Implementation Plan for Phase 2 Pipeline

## Architecture Overview

A **centralized validation module** with checkpoint-based validation throughout the pipeline to ensure data integrity and scientific validity for in-memory processing.

## 1. Create New Validation Module

**File**: `phase2/core/validation.py`

**Why here?**
- Centralized validation logic (DRY principle)
- Easy to test independently
- Can be imported anywhere needed
- Keeps run_pipeline.py cleaner

**Structure**:
```python
# phase2/core/validation.py
"""
Comprehensive validation system for in-memory pipeline operations.
Ensures data integrity and scientific validity throughout processing.
"""

import numpy as np
from typing import List, Dict, Optional, Tuple
from cell import PetriDish, Cell
import logging

# Configure logging for validation issues
logger = logging.getLogger(__name__)

class ValidationError(Exception):
    """Custom exception for validation failures."""
    pass

class ValidationWarning:
    """Container for non-fatal validation issues."""
    pass

# Main validation orchestrator
class PipelineValidator:
    def __init__(self, verbose: bool = True):
        self.verbose = verbose
        self.warnings = []
        self.metrics = {}
    
    # Individual validation methods...
```

## 2. Integration Points in Pipeline

### Where to add validation calls in `run_pipeline.py`:

#### STAGE 1 - After Snapshot Extraction
```python
# Line ~615 (after loading first_snapshot_cells)
from core.validation import PipelineValidator
validator = PipelineValidator(verbose=args.verbose)

try:
    validator.validate_snapshot(
        cells=first_snapshot_cells,
        expected_year=args.first_snapshot,
        expected_gene_rate_groups=gene_rate_groups,
        min_cells=args.n_quantiles * args.cells_per_quantile * 2  # Need enough for sampling
    )
except ValidationError as e:
    print(f"❌ Snapshot validation failed: {e}")
    sys.exit(1)
```

#### STAGE 3 - After Individual Creation
```python
# Line ~675 (after creating both mutant and control1)
validator.validate_initial_individuals(
    mutant_dishes=mutant_dishes,
    control1_dishes=control1_dishes,
    expected_count=expected_individuals,
    expected_gene_rate_groups=gene_rate_groups,
    snapshot_year=args.first_snapshot
)
```

#### STAGE 4 - After Growth
```python
# Line ~735 (after growing both batches)
extinction_report = validator.validate_grown_individuals(
    mutant_dishes=mutant_dishes,
    control1_dishes=control1_dishes,
    expected_population=expected_population,
    growth_years=timeline_duration,
    allow_extinction=True  # Flag to handle gracefully
)

if extinction_report['total_extinct'] > 0:
    print(f"⚠️ Warning: {extinction_report['total_extinct']} individuals went extinct during growth")
    if extinction_report['complete_batch_extinction']:
        print("❌ Fatal: Entire batch extinct. Consider different parameters.")
        sys.exit(1)
```

#### STAGE 6 - Multiple Checkpoints
```python
# After normalization (line ~815)
if args.normalize_size:
    validator.validate_normalized_populations(
        mutant_before=mutant_dishes_before_norm,  # Keep reference
        control1_before=control1_dishes_before_norm,
        mutant_after=mutant_dishes,
        control1_after=control1_dishes,
        threshold=normalization_threshold
    )

# After mixing but before save (line ~925)
validator.validate_mixed_populations(
    mutant_dishes=mutant_dishes,
    control1_dishes=control1_dishes,
    mix_ratio=args.mix_ratio,
    uniform_mixing=args.uniform_mixing,
    second_snapshot_year=args.second_snapshot
)

# Right before save (line ~940)
validator.validate_before_save(
    mutant_dishes=mutant_dishes,
    control1_dishes=control1_dishes
)
```

## 3. Detailed Validation Functions

### Critical Validations (Fail Pipeline)
```python
def validate_snapshot(self, cells: List[Cell], expected_year: int, 
                     expected_gene_rate_groups: List[Tuple[int, float]],
                     min_cells: int) -> None:
    """
    Validate snapshot cells for biological and technical correctness.
    Raises ValidationError on critical issues.
    """
    # Critical checks (fail pipeline)
    if not cells:
        raise ValidationError("Snapshot has no cells")
    
    if len(cells) < min_cells:
        raise ValidationError(
            f"Insufficient cells for sampling: {len(cells)} < {min_cells} required"
        )
    
    # Check all cells have correct age
    ages = set(cell.age for cell in cells)
    if len(ages) > 1:
        raise ValidationError(f"Inconsistent cell ages in snapshot: {ages}")
    if ages.pop() != expected_year:
        raise ValidationError(f"Cell age doesn't match snapshot year")
    
    # Verify gene_rate_groups
    for i, cell in enumerate(cells):
        if cell.gene_rate_groups != expected_gene_rate_groups:
            raise ValidationError(
                f"Cell {i} has wrong gene_rate_groups: "
                f"{cell.gene_rate_groups} != {expected_gene_rate_groups}"
            )
    
    # Check methylation validity
    for cell in cells:
        methylation = np.mean([site for site in cell.cpg_sites])
        max_possible = expected_year * 0.005 * 3  # 3x safety factor
        if methylation > max_possible:
            raise ValidationError(
                f"Cell has impossible methylation level: {methylation:.3f} > {max_possible:.3f}"
            )
    
    # Warning checks (don't fail)
    if len(cells) == len(set(id(cell) for cell in cells)):
        pass  # All unique objects, good
    else:
        self.warnings.append("Duplicate cell objects detected in snapshot")
```

### Growth Validation with Extinction Handling
```python
def validate_grown_individuals(self, mutant_dishes: List[PetriDish],
                              control1_dishes: List[PetriDish],
                              expected_population: int,
                              growth_years: int,
                              allow_extinction: bool = True) -> Dict:
    """
    Validate grown individuals, handling extinction gracefully.
    Returns extinction report.
    """
    report = {
        'mutant_extinct': [],
        'control1_extinct': [],
        'mutant_low_pop': [],
        'control1_low_pop': [],
        'total_extinct': 0,
        'complete_batch_extinction': False
    }
    
    # Check each individual
    for i, dish in enumerate(mutant_dishes):
        cell_count = len(dish.cells)
        if cell_count == 0:
            report['mutant_extinct'].append(dish.metadata.get('individual_id', i+1))
        elif cell_count < expected_population * 0.3:
            report['mutant_low_pop'].append({
                'id': dish.metadata.get('individual_id', i+1),
                'cells': cell_count
            })
        
        # Validate growth history
        if hasattr(dish, 'cell_history'):
            if len(dish.cell_history) != growth_years + 1:  # +1 for year 0
                raise ValidationError(
                    f"Individual {i+1} has incomplete history: "
                    f"{len(dish.cell_history)} entries for {growth_years} years"
                )
    
    # Similar for control1...
    
    # Check for complete batch extinction
    if all(len(d.cells) == 0 for d in mutant_dishes):
        report['complete_batch_extinction'] = True
        if not allow_extinction:
            raise ValidationError("Complete mutant batch extinction")
    
    report['total_extinct'] = len(report['mutant_extinct']) + len(report['control1_extinct'])
    
    return report
```

## 4. Metadata Preservation Validation
```python
def validate_metadata_integrity(self, dishes: List[PetriDish], 
                               stage: str) -> None:
    """
    Ensure metadata is preserved and accumulated correctly.
    """
    required_fields = {
        'initial': ['individual_id', 'individual_type'],
        'grown': ['individual_id', 'individual_type', 'initial_year'],
        'normalized': ['individual_id', 'individual_type', 'normalized'],
        'mixed': ['individual_id', 'individual_type', 'mixed', 'mix_ratio']
    }
    
    for dish in dishes:
        if not hasattr(dish, 'metadata'):
            raise ValidationError(f"No metadata at {stage}")
        
        # Check required fields for this stage
        for field in required_fields.get(stage, []):
            if field not in dish.metadata:
                raise ValidationError(
                    f"Missing {field} in metadata at {stage}"
                )
        
        # Check ID uniqueness
        ids = [d.metadata.get('individual_id') for d in dishes]
        if len(ids) != len(set(ids)):
            raise ValidationError(f"Duplicate individual_ids at {stage}")
```

## 5. Scientific Validity Checks
```python
def validate_biological_consistency(self, dish: PetriDish) -> None:
    """
    Validate biological/scientific consistency.
    """
    # Methylation should only increase
    if hasattr(dish, 'cell_history'):
        prev_methylation = 0
        for year, cells_data in sorted(dish.cell_history.items()):
            year_methylation = np.mean([
                np.mean(cell.get('methylated', [])) 
                for cell in cells_data
            ])
            if year_methylation < prev_methylation * 0.95:  # 5% tolerance
                self.warnings.append(
                    f"Methylation decreased at year {year}: "
                    f"{prev_methylation:.3f} -> {year_methylation:.3f}"
                )
            prev_methylation = year_methylation
    
    # Check gene JSD values are valid
    if hasattr(dish, 'metadata') and 'gene_jsds' in dish.metadata:
        jsds = dish.metadata['gene_jsds']
        if any(j < 0 or j > 1 for j in jsds):
            raise ValidationError("Invalid gene JSD values outside [0,1]")
```

## 6. Pre-Save Final Validation
```python
def validate_before_save(self, mutant_dishes: List[PetriDish],
                        control1_dishes: List[PetriDish]) -> None:
    """
    Final validation before saving to disk.
    Ensures everything is ready for serialization.
    """
    all_dishes = mutant_dishes + control1_dishes
    
    for dish in all_dishes:
        # Must have metadata
        if not hasattr(dish, 'metadata'):
            raise ValidationError("Missing metadata before save")
        
        # Must have individual_id
        if 'individual_id' not in dish.metadata:
            raise ValidationError("Missing individual_id before save")
        
        # Check for None values in critical fields
        critical = ['individual_id', 'individual_type']
        for field in critical:
            if dish.metadata.get(field) is None:
                raise ValidationError(f"None value in {field}")
        
        # Verify cells are valid
        if not dish.cells:
            self.warnings.append(
                f"Individual {dish.metadata.get('individual_id')} "
                f"has no cells (extinct?)"
            )
```

## 7. Configuration and Flexibility

Add validation configuration to allow different strictness levels:

```python
class ValidationConfig:
    """Configuration for validation strictness."""
    
    # Extinction handling
    ALLOW_EXTINCTION = True
    MAX_EXTINCTION_RATE = 0.5  # Fail if >50% extinct
    
    # Population thresholds  
    MIN_POPULATION_FRACTION = 0.3  # Warn if <30% of expected
    MAX_POPULATION_FRACTION = 3.0  # Warn if >300% of expected
    
    # Methylation validation
    MAX_METHYLATION_FACTOR = 3.0  # Max 3x expected methylation
    
    # Sampling requirements
    MIN_CELLS_PER_QUANTILE = 10  # Need enough for statistics
    
    @classmethod
    def strict(cls):
        """Strict validation for production."""
        cls.ALLOW_EXTINCTION = False
        cls.MAX_EXTINCTION_RATE = 0.1
        return cls
    
    @classmethod  
    def relaxed(cls):
        """Relaxed validation for experimentation."""
        cls.ALLOW_EXTINCTION = True
        cls.MAX_EXTINCTION_RATE = 0.9
        return cls
```

## 8. Integration with Logging

Add structured logging for debugging:

```python
import logging
import json

class ValidationLogger:
    """Structured logging for validation events."""
    
    def __init__(self, log_file: Optional[str] = None):
        self.logger = logging.getLogger('validation')
        if log_file:
            handler = logging.FileHandler(log_file)
            handler.setFormatter(
                logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
            )
            self.logger.addHandler(handler)
    
    def log_validation_point(self, stage: str, metrics: Dict):
        """Log validation metrics at each stage."""
        self.logger.info(f"Validation at {stage}: {json.dumps(metrics)}")
    
    def log_warning(self, warning: str):
        self.logger.warning(warning)
    
    def log_error(self, error: str):
        self.logger.error(error)
```

## 9. Error Recovery Suggestions

Include helpful error messages:

```python
class ValidationErrorWithSuggestion(ValidationError):
    """Validation error with recovery suggestions."""
    
    def __init__(self, message: str, suggestion: str = None):
        super().__init__(message)
        self.suggestion = suggestion

# Usage:
if extinction_rate > 0.8:
    raise ValidationErrorWithSuggestion(
        f"Extinction rate too high: {extinction_rate:.1%}",
        suggestion=(
            "Try: 1) Increasing growth_phase, "
            "2) Reducing growth years, "
            "3) Using different seed, "
            "4) Adjusting methylation rate"
        )
    )
```

## 10. Testing Strategy

Create test file: `phase2/tests/test_validation.py`
- Test each validation function independently
- Test edge cases (extinction, no cells, etc.)
- Test metadata preservation
- Test warning accumulation
- Test with mock data for speed

## Implementation Priority

1. **First**: Core validation module structure
2. **Second**: Critical validation functions (snapshot, initial, save)
3. **Third**: Integration points in run_pipeline.py
4. **Fourth**: Growth and mixing validation
5. **Fifth**: Logging and error recovery
6. **Last**: Configuration flexibility

## Benefits

This architecture provides:
- **Comprehensive validation** without cluttering main pipeline
- **Flexible configuration** for different use cases
- **Clear error messages** with recovery suggestions
- **Performance tracking** through metrics
- **Backward compatibility** (can disable if needed)
- **Testability** through modular design

## Scientific Goals Addressed

### Core Scientific Goals:
1. Compare epigenetic drift between different sampling strategies (quantile vs uniform)
2. Ensure statistical validity of comparisons between mutant/control1/control2
3. Track methylation accumulation accurately through growth and mixing
4. Preserve gene-specific rate effects through all transformations

### Critical Validation Points:

#### POST-STAGE 1 (Snapshot Extraction)
- All cells must have correct age (should match snapshot year)
- All cells must have valid methylation patterns (0s and 1s only)
- Gene_rate_groups must match the simulation exactly
- Cell JSD values must be in valid range [0, 1]
- No cells should have more methylation than mathematically possible for their age
- Methylation should follow expected distribution (not all 0s or all 1s)
- Minimum cell count (need enough for sampling)
- No duplicate cells (each should be unique object)
- All cells must have required attributes (cpg_sites, age, cell_JSD, gene_rate_groups)

#### POST-STAGE 3 (Individual Creation)
- Mutant: Exactly n_quantiles * cells_per_quantile individuals
- Control1: Same count as mutant
- Each individual must have exactly 1 cell (founding cell)
- All founding cells must be from the snapshot (age = snapshot year)
- Each individual must have unique individual_id (1, 2, 3, ...)
- Individual_ids must be sequential with no gaps
- Must have individual_type ('mutant' or 'control1')
- Mutant: must have source_quantile metadata
- Control1: must have source='uniform' metadata
- No two individuals should have the same founding cell (unique sampling)
- Quantile distribution must be preserved (correct cells in each quantile)
- All PetriDish objects must have growth_phase matching args

#### POST-STAGE 4 (Growth)
- Cell count must be in range [target * 0.3, target * 3.0] (homeostasis variation)
- BUT: Handle extinction case (0 cells is valid but needs flagging)
- All cells in each individual must have age = snapshot_year + growth_years
- Cell history must have entries for each year of growth
- Year 0 of history = single cell from snapshot
- Final year of history = current cells
- Mean methylation must increase monotonically in history
- No cell should lose methylation (it's irreversible)
- Methylation rate should match expected accumulation
- All cells must maintain original gene_rate_groups
- No gene_rate_groups should be modified during growth

#### PRE-STAGE 6 (Before Mixing)
- Warn if >50% of individuals have <50% of target population
- Flag if any batch has high extinction rate
- Ensure both mutant and control1 have at least 1 viable individual
- Handle if ALL individuals went extinct
- Handle if one batch extinct but other survived
- Validate enough cells for meaningful mixing

#### POST-NORMALIZATION (If Applied)
- All retained individuals must have exactly threshold_size cells
- Excluded individuals must not be in the lists
- Metadata must reflect normalization (normalized=True, threshold recorded)
- Individual_ids must be preserved (not renumbered!)
- May have gaps in IDs (e.g., 1,3,4,7 if 2,5,6 excluded)
- Original metadata must be intact
- Warn if one batch loses >80% of individuals
- Error if either batch completely excluded
- Track and report retention rates

#### POST-MIXING
- Expected final cell count per individual (based on mix_ratio)
- Verify mix actually happened (not all cells from one source)
- Check age distribution matches expected (mix of ages)
- All individuals must have identical final cell count (uniform mixing)
- Verify uniform pool used consistently
- mixed=True, mix_ratio recorded
- final_cells count recorded
- All original metadata preserved

#### PRE-SAVE VALIDATION
- Every individual has required metadata fields
- No null/None values in critical fields
- Cell histories intact and complete
- No duplicate individual_ids within batch
- File names will be unique
- Gene JSDs calculated if requested
- Mean methylation values reasonable
- No data corruption indicators

#### CONTROL2 VALIDATION
- growth_phase must be None
- All cells must have age = second_snapshot_year
- No growth history (or single timepoint only)
- Count matches expected based on other batches

## Edge Cases to Handle

**Extinction scenarios:**
- Individual with 0 cells after growth
- Entire batch extinct
- Mix with 0 cells (what to mix?)

**Extreme distributions:**
- All cells in one quantile
- Snapshot with very few cells
- Mix ratio 0% or 100%

**Data anomalies:**
- Cells with corrupted methylation data
- Missing or malformed metadata
- Memory pressure with large populations

**User errors:**
- Inconsistent parameters between runs
- Wrong snapshot years
- Incompatible gene_rate_groups

## Cross-Stage Consistency Tracking
- Total individual count per batch
- Individual_id mapping and uniqueness  
- Gene_rate_groups never change
- Metadata accumulation (only add, never remove)
- Memory objects not accidentally shared/aliased