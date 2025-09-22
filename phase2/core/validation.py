"""
Comprehensive validation system for in-memory pipeline operations.
Ensures data integrity and scientific validity throughout processing.
"""

import numpy as np
from typing import List, Dict, Optional, Tuple, Any
import sys
import os
import logging
import copy

# Add paths for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import PetriDish, Cell

# Configure logging for validation issues
logger = logging.getLogger(__name__)


class ValidationError(Exception):
    """Custom exception for validation failures."""
    pass


class ValidationWarning:
    """Container for non-fatal validation issues."""
    def __init__(self, message: str):
        self.message = message
    
    def __str__(self):
        return self.message


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
        config = copy.copy(cls)
        config.ALLOW_EXTINCTION = False
        config.MAX_EXTINCTION_RATE = 0.1
        return config
    
    @classmethod  
    def relaxed(cls):
        """Relaxed validation for experimentation."""
        config = copy.copy(cls)
        config.ALLOW_EXTINCTION = True
        config.MAX_EXTINCTION_RATE = 0.9
        return config


class PipelineValidator:
    """Main validation orchestrator for the pipeline."""
    
    def __init__(self, config: Optional[ValidationConfig] = None, verbose: bool = True):
        self.config = config or ValidationConfig()
        self.verbose = verbose
        self.warnings = []
        self.metrics = {}
        
    def _log(self, message: str, level: str = "info"):
        """Internal logging helper."""
        if self.verbose:
            print(f"  [{level.upper()}] {message}")
        
        if level == "warning":
            self.warnings.append(ValidationWarning(message))
            logger.warning(message)
        elif level == "error":
            logger.error(message)
        else:
            logger.info(message)
    
    def validate_snapshot(self, cells: List[Cell], expected_year: int, 
                         expected_gene_rate_groups: List[Tuple[int, float]],
                         min_cells: int) -> None:
        """
        Validate snapshot cells for biological and technical correctness.
        Raises ValidationError on critical issues.
        """
        self._log(f"Validating snapshot: {len(cells)} cells from year {expected_year}")
        
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
        
        actual_age = ages.pop() if ages else 0
        if actual_age != expected_year:
            raise ValidationError(
                f"Cell age doesn't match snapshot year: {actual_age} != {expected_year}"
            )
        
        # Verify gene_rate_groups (check first 10 cells for efficiency)
        sample_size = min(10, len(cells))
        for i in range(sample_size):
            cell = cells[i]
            if not hasattr(cell, 'gene_rate_groups'):
                raise ValidationError(f"Cell {i} missing gene_rate_groups attribute")
            if cell.gene_rate_groups != expected_gene_rate_groups:
                raise ValidationError(
                    f"Cell {i} has wrong gene_rate_groups: "
                    f"{cell.gene_rate_groups} != {expected_gene_rate_groups}"
                )
        
        # Check methylation validity
        methylation_values = []
        for i, cell in enumerate(cells[:100]):  # Check first 100 for efficiency
            methylation = np.mean(cell.cpg_sites)
            methylation_values.append(methylation)
            
            # Maximum possible methylation (with safety factor)
            max_possible = min(1.0, expected_year * 0.005 * self.config.MAX_METHYLATION_FACTOR)
            if methylation > max_possible:
                raise ValidationError(
                    f"Cell {i} has impossible methylation level: {methylation:.3f} > {max_possible:.3f}"
                )
        
        # Check for data anomalies
        mean_meth = np.mean(methylation_values)
        if mean_meth == 0:
            self._log("WARNING: All cells completely unmethylated", "warning")
        elif mean_meth == 1:
            self._log("WARNING: All cells completely methylated", "warning")
        
        # Check for duplicate objects (memory aliasing)
        unique_ids = len(set(id(cell) for cell in cells[:100]))
        if unique_ids < min(100, len(cells)):
            self._log("WARNING: Duplicate cell objects detected in snapshot", "warning")
        
        # Get gene rate groups info for display
        if cells and hasattr(cells[0], 'gene_rate_groups') and cells[0].gene_rate_groups:
            grg_info = f", gene_rate_groups: {cells[0].gene_rate_groups}"
        else:
            grg_info = ""
        self._log(f"✓ Snapshot validation passed: {len(cells)} valid cells{grg_info}")
        
    def validate_initial_individuals(self, 
                                   mutant_dishes: List[PetriDish],
                                   control1_dishes: List[PetriDish],
                                   expected_count: int,
                                   expected_gene_rate_groups: List[Tuple[int, float]],
                                   snapshot_year: int) -> None:
        """
        Validate newly created individuals before growth.
        """
        self._log(f"Validating initial individuals: {len(mutant_dishes)} mutant, {len(control1_dishes)} control1")
        
        # Check counts
        if len(mutant_dishes) != expected_count:
            raise ValidationError(
                f"Wrong mutant count: {len(mutant_dishes)} != {expected_count}"
            )
        if len(control1_dishes) != expected_count:
            raise ValidationError(
                f"Wrong control1 count: {len(control1_dishes)} != {expected_count}"
            )
        
        # Validate each batch
        for batch_name, dishes in [("mutant", mutant_dishes), ("control1", control1_dishes)]:
            seen_ids = set()
            
            for i, dish in enumerate(dishes):
                # Check single founding cell
                if len(dish.cells) != 1:
                    raise ValidationError(
                        f"{batch_name} individual {i} has {len(dish.cells)} cells, expected 1"
                    )
                
                # Check founding cell age
                cell = dish.cells[0]
                if cell.age != snapshot_year:
                    raise ValidationError(
                        f"{batch_name} individual {i} founding cell has age {cell.age}, "
                        f"expected {snapshot_year}"
                    )
                
                # Check metadata
                if not hasattr(dish, 'metadata'):
                    raise ValidationError(f"{batch_name} individual {i} missing metadata")
                
                metadata = dish.metadata
                
                # Check individual_id
                if 'individual_id' not in metadata:
                    raise ValidationError(f"{batch_name} individual {i} missing individual_id")
                
                ind_id = metadata['individual_id']
                if ind_id in seen_ids:
                    raise ValidationError(f"Duplicate individual_id: {ind_id}")
                seen_ids.add(ind_id)
                
                # Check individual_type
                if metadata.get('individual_type') != batch_name:
                    raise ValidationError(
                        f"Wrong individual_type: {metadata.get('individual_type')} != {batch_name}"
                    )
                
                # Check batch-specific metadata
                if batch_name == "mutant":
                    if 'source_quantile' not in metadata:
                        self._log(f"Mutant {ind_id} missing source_quantile", "warning")
                elif batch_name == "control1":
                    if metadata.get('source') != 'uniform':
                        self._log(f"Control1 {ind_id} missing source='uniform'", "warning")
                
                # Check gene_rate_groups
                if cell.gene_rate_groups != expected_gene_rate_groups:
                    raise ValidationError(
                        f"{batch_name} {ind_id} has wrong gene_rate_groups"
                    )
        
        self._log(f"✓ Initial individuals validation passed: {len(mutant_dishes)} mutant, {len(control1_dishes)} control1")
    
    def validate_grown_individuals(self, 
                                  mutant_dishes: List[PetriDish],
                                  control1_dishes: List[PetriDish],
                                  expected_population: int,
                                  growth_years: int,
                                  allow_extinction: Optional[bool] = None) -> Dict:
        """
        Validate grown individuals, handling extinction gracefully.
        Returns extinction report.
        """
        if allow_extinction is None:
            allow_extinction = self.config.ALLOW_EXTINCTION
            
        self._log(f"Validating grown individuals after {growth_years} years")
        
        report = {
            'mutant_extinct': [],
            'control1_extinct': [],
            'mutant_low_pop': [],
            'control1_low_pop': [],
            'total_extinct': 0,
            'complete_batch_extinction': False
        }
        
        # Process each batch
        for batch_name, dishes in [("mutant", mutant_dishes), ("control1", control1_dishes)]:
            extinct_ids = []
            low_pop_ids = []
            
            for dish in dishes:
                ind_id = dish.metadata.get('individual_id', 0)
                cell_count = len(dish.cells)
                
                # Check for extinction
                if cell_count == 0:
                    extinct_ids.append(ind_id)
                    self._log(f"{batch_name} {ind_id}: EXTINCT", "warning")
                    continue
                
                # Check population size
                min_expected = int(expected_population * self.config.MIN_POPULATION_FRACTION)
                max_expected = int(expected_population * self.config.MAX_POPULATION_FRACTION)
                
                if cell_count < min_expected:
                    low_pop_ids.append({'id': ind_id, 'cells': cell_count})
                    self._log(
                        f"{batch_name} {ind_id}: Low population {cell_count} < {min_expected}", 
                        "warning"
                    )
                elif cell_count > max_expected:
                    self._log(
                        f"{batch_name} {ind_id}: High population {cell_count} > {max_expected}", 
                        "warning"
                    )
                
                # Validate growth history
                if hasattr(dish, 'cell_history'):
                    history_years = len(dish.cell_history)
                    expected_history = growth_years + 1  # +1 for year 0
                    
                    if history_years != expected_history:
                        self._log(
                            f"{batch_name} {ind_id}: History mismatch {history_years} != {expected_history}",
                            "warning"
                        )
                    
                    # Check methylation progression
                    prev_meth = 0
                    for year_str in sorted(dish.cell_history.keys(), key=int):
                        cells_data = dish.cell_history[year_str]
                        if cells_data:  # Skip if no cells (extinction)
                            # Get mean methylation, handling both empty lists and missing keys
                            methylation_values = []
                            for cell in cells_data:
                                # Support both 'cpg_sites' (correct) and 'methylated' (legacy) keys
                                cpg_sites = cell.get('cpg_sites') or cell.get('methylated', [])
                                if cpg_sites:  # Only include if not empty
                                    methylation_values.append(np.mean(cpg_sites))
                            
                            if methylation_values:
                                year_meth = np.mean(methylation_values)
                            else:
                                year_meth = 0.0  # Default to 0 if no valid data
                            
                            if year_meth < prev_meth * 0.95:  # 5% tolerance
                                self._log(
                                    f"{batch_name} {ind_id}: Methylation decreased at year {year_str}",
                                    "warning"
                                )
                            prev_meth = year_meth
            
            # Store in report
            if batch_name == "mutant":
                report['mutant_extinct'] = extinct_ids
                report['mutant_low_pop'] = low_pop_ids
            else:
                report['control1_extinct'] = extinct_ids
                report['control1_low_pop'] = low_pop_ids
        
        # Calculate totals
        report['total_extinct'] = len(report['mutant_extinct']) + len(report['control1_extinct'])
        total_individuals = len(mutant_dishes) + len(control1_dishes)
        
        # Check extinction rates
        if total_individuals > 0:
            extinction_rate = report['total_extinct'] / total_individuals
            
            if extinction_rate > self.config.MAX_EXTINCTION_RATE:
                msg = f"Extinction rate too high: {extinction_rate:.1%}"
                if not allow_extinction:
                    raise ValidationError(msg)
                else:
                    self._log(msg, "warning")
        
        # Check for complete batch extinction
        if all(len(d.cells) == 0 for d in mutant_dishes):
            report['complete_batch_extinction'] = True
            if not allow_extinction:
                raise ValidationError("Complete mutant batch extinction")
        
        if all(len(d.cells) == 0 for d in control1_dishes):
            report['complete_batch_extinction'] = True
            if not allow_extinction:
                raise ValidationError("Complete control1 batch extinction")
        
        total_individuals = len(mutant_dishes) + len(control1_dishes)
        self._log(f"✓ Growth validation complete: {total_individuals} individuals, {report['total_extinct']} extinct")
        return report
    
    def validate_normalized_populations(self,
                                       mutant_after: List[PetriDish],
                                       control1_after: List[PetriDish],
                                       threshold: int) -> None:
        """
        Validate populations after normalization.
        """
        self._log(f"Validating normalized populations (threshold: {threshold} cells)")
        
        # Check all individuals have threshold size
        for batch_name, dishes in [("mutant", mutant_after), ("control1", control1_after)]:
            for dish in dishes:
                cell_count = len(dish.cells)
                if cell_count != threshold:
                    raise ValidationError(
                        f"{batch_name} {dish.metadata.get('individual_id')}: "
                        f"Size {cell_count} != threshold {threshold}"
                    )
                
                # Check metadata updated
                if not dish.metadata.get('normalized'):
                    self._log(
                        f"{batch_name} {dish.metadata.get('individual_id')}: "
                        f"Missing 'normalized' flag",
                        "warning"
                    )
        
        # Check retention
        if len(mutant_after) == 0 and len(control1_after) == 0:
            raise ValidationError("All individuals excluded by normalization!")
        
        total_kept = len(mutant_after) + len(control1_after)
        self._log(f"✓ Normalization validation passed: {total_kept} individuals kept (threshold: {threshold} cells)")
    
    def validate_mixed_populations(self,
                                  mutant_dishes: List[PetriDish],
                                  control1_dishes: List[PetriDish],
                                  mix_ratio: int,
                                  uniform_mixing: bool,  # Kept for compatibility, always True now
                                  second_snapshot_year: int) -> None:
        """
        Validate populations after mixing.
        """
        self._log(f"Validating mixed populations (ratio: {mix_ratio}%, uniform mixing)")
        
        # Calculate expected composition
        grown_fraction = (100 - mix_ratio) / 100.0
        
        for batch_name, dishes in [("mutant", mutant_dishes), ("control1", control1_dishes)]:
            if not dishes:
                self._log(f"No {batch_name} individuals to validate", "warning")
                continue
                
            cell_counts = []
            
            for dish in dishes:
                cell_count = len(dish.cells)
                cell_counts.append(cell_count)
                
                # Check metadata
                if not dish.metadata.get('mixed'):
                    raise ValidationError(
                        f"{batch_name} {dish.metadata.get('individual_id')}: Not marked as mixed"
                    )
                
                if dish.metadata.get('mix_ratio') != mix_ratio:
                    self._log(
                        f"{batch_name} {dish.metadata.get('individual_id')}: "
                        f"Wrong mix_ratio in metadata",
                        "warning"
                    )
                
                # Check age distribution 
                # NOTE: In the current experimental design, grown cells (aged from year 30 to 50)
                # and snapshot cells (from year 50) both end up at age 50, so uniform age is expected.
                # This validation is disabled as it produces false positives.
                #
                # TODO: Re-enable with smarter logic that accounts for experimental design:
                # - If first_snapshot + growth_years == second_snapshot, expect uniform age
                # - Otherwise, expect mixed ages
                #
                # if dish.cells:
                #     ages = [cell.age for cell in dish.cells[:10]]  # Sample first 10
                #     unique_ages = set(ages)
                #     
                #     # Should have cells from snapshot (older) and grown (younger)
                #     if len(unique_ages) < 2 and mix_ratio not in [0, 100]:
                #         self._log(
                #             f"{batch_name} {dish.metadata.get('individual_id')}: "
                #             f"Suspicious age distribution: {unique_ages}",
                #             "warning"
                #         )
            
            # Check uniform mixing consistency (always enabled now)
            if len(set(cell_counts)) > 1:
                self._log(
                    f"{batch_name}: Non-uniform sizes in uniform mixing: {set(cell_counts)}",
                    "warning"
                )
        
        total_mixed = len(mutant_dishes) + len(control1_dishes)
        self._log(f"✓ Mixing validation passed: {len(mutant_dishes)} mutant, {len(control1_dishes)} control1 (ratio: {mix_ratio}%)")
    
    def validate_before_save(self,
                           mutant_dishes: List[PetriDish],
                           control1_dishes: List[PetriDish]) -> None:
        """
        Final validation before saving to disk.
        """
        self._log("Performing final pre-save validation")
        
        # Check each batch separately (they're saved in different directories)
        for batch_name, dishes in [("mutant", mutant_dishes), ("control1", control1_dishes)]:
            batch_ids = []
            
            for dish in dishes:
                # Must have metadata
                if not hasattr(dish, 'metadata'):
                    raise ValidationError(f"{batch_name}: Missing metadata before save")
                
                # Must have individual_id
                if 'individual_id' not in dish.metadata:
                    raise ValidationError(f"{batch_name}: Missing individual_id before save")
                
                ind_id = dish.metadata['individual_id']
                batch_ids.append(ind_id)
                
                # Check for None values in critical fields
                critical = ['individual_id', 'individual_type']
                for field in critical:
                    if dish.metadata.get(field) is None:
                        raise ValidationError(f"{batch_name}: None value in {field} for individual {ind_id}")
                
                # Verify cells are valid (unless extinct)
                if not dish.cells:
                    self._log(
                        f"{batch_name} {ind_id} has no cells (extinct?)",
                        "warning"
                    )
                else:
                    # Check cells have required attributes
                    sample_cell = dish.cells[0]
                    if not hasattr(sample_cell, 'cpg_sites'):
                        raise ValidationError(f"{batch_name} {ind_id}: cells missing cpg_sites")
                    if not hasattr(sample_cell, 'age'):
                        raise ValidationError(f"{batch_name} {ind_id}: cells missing age")
            
            # Check for duplicate IDs within batch
            if len(batch_ids) != len(set(batch_ids)):
                raise ValidationError(f"{batch_name}: Duplicate individual_ids detected: {batch_ids}")
        
        total_dishes = len(mutant_dishes) + len(control1_dishes)
        self._log(f"✓ Pre-save validation passed: {len(mutant_dishes)} mutant, {len(control1_dishes)} control1 ready to save")
    
    def validate_control2(self,
                         control2_dishes: List[PetriDish],
                         second_snapshot_year: int,
                         expected_count: int) -> None:
        """
        Validate control2 individuals (pure snapshot).
        """
        self._log(f"Validating {len(control2_dishes)} control2 individuals")
        
        if len(control2_dishes) != expected_count:
            self._log(
                f"Control2 count mismatch: {len(control2_dishes)} != {expected_count}",
                "warning"
            )
        
        for dish in control2_dishes:
            # Check growth_phase is None
            if dish.growth_phase is not None:
                raise ValidationError(
                    f"Control2 {dish.metadata.get('individual_id')}: "
                    f"growth_phase should be None, got {dish.growth_phase}"
                )
            
            # Check all cells have correct age
            if dish.cells:
                ages = set(cell.age for cell in dish.cells[:10])  # Sample
                if len(ages) > 1:
                    raise ValidationError(
                        f"Control2 {dish.metadata.get('individual_id')}: Mixed ages {ages}"
                    )
                if ages and ages.pop() != second_snapshot_year:
                    raise ValidationError(
                        f"Control2 {dish.metadata.get('individual_id')}: "
                        f"Wrong cell age, expected {second_snapshot_year}"
                    )
            
            # Check metadata
            if dish.metadata.get('individual_type') != 'control2':
                raise ValidationError(
                    f"Wrong individual_type: {dish.metadata.get('individual_type')}"
                )
        
        self._log(f"✓ Control2 validation passed: {len(control2_dishes)} individuals created")
    
    def get_summary(self) -> Dict[str, Any]:
        """Get summary of validation results."""
        return {
            'warnings_count': len(self.warnings),
            'warnings': [str(w) for w in self.warnings],
            'metrics': self.metrics
        }