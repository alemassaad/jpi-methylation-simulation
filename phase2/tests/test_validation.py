#!/usr/bin/env python3
"""
Comprehensive tests for the validation system.
Tests all validation functions and edge cases.
"""

import sys
import os
import copy
import numpy as np

# Add paths for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from core.validation import PipelineValidator, ValidationError, ValidationConfig


def create_test_cells(n_cells: int, age: int, gene_rate_groups=None, methylation_fraction=0.1):
    """Create test cells with specified properties."""
    if gene_rate_groups is None:
        gene_rate_groups = [(20, 0.005)]
    
    cells = []
    for i in range(n_cells):
        cell = Cell(n=100, gene_rate_groups=gene_rate_groups, gene_size=5)
        cell.age = age
        # Add some methylation
        n_methylated = int(len(cell.cpg_sites) * methylation_fraction)
        for j in range(n_methylated):
            cell.cpg_sites[j] = 1
        cells.append(cell)
    return cells


def test_snapshot_validation():
    """Test snapshot validation with various scenarios."""
    print("\n" + "="*60)
    print("TEST: Snapshot Validation")
    print("="*60)
    
    validator = PipelineValidator(verbose=False)
    gene_rate_groups = [(20, 0.005)]
    
    # Test 1: Valid snapshot
    print("\n1. Valid snapshot...")
    cells = create_test_cells(100, age=30, gene_rate_groups=gene_rate_groups)
    try:
        validator.validate_snapshot(cells, expected_year=30, 
                                   expected_gene_rate_groups=gene_rate_groups,
                                   min_cells=50)
        print("   ✓ Valid snapshot passed")
    except ValidationError as e:
        print(f"   ✗ Unexpected error: {e}")
        return False
    
    # Test 2: Empty snapshot
    print("\n2. Empty snapshot...")
    try:
        validator.validate_snapshot([], expected_year=30,
                                   expected_gene_rate_groups=gene_rate_groups,
                                   min_cells=50)
        print("   ✗ Should have raised error for empty snapshot")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 3: Insufficient cells
    print("\n3. Insufficient cells...")
    cells = create_test_cells(10, age=30)
    try:
        validator.validate_snapshot(cells, expected_year=30,
                                   expected_gene_rate_groups=gene_rate_groups,
                                   min_cells=50)
        print("   ✗ Should have raised error for insufficient cells")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 4: Wrong age
    print("\n4. Wrong cell age...")
    cells = create_test_cells(100, age=40)  # Wrong age
    try:
        validator.validate_snapshot(cells, expected_year=30,  # Expected 30
                                   expected_gene_rate_groups=gene_rate_groups,
                                   min_cells=50)
        print("   ✗ Should have raised error for wrong age")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 5: Mixed ages
    print("\n5. Mixed cell ages...")
    cells = create_test_cells(50, age=30) + create_test_cells(50, age=40)
    try:
        validator.validate_snapshot(cells, expected_year=30,
                                   expected_gene_rate_groups=gene_rate_groups,
                                   min_cells=50)
        print("   ✗ Should have raised error for mixed ages")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 6: Wrong gene_rate_groups
    print("\n6. Wrong gene_rate_groups...")
    cells = create_test_cells(100, age=30, gene_rate_groups=[(10, 0.004), (10, 0.006)])
    try:
        validator.validate_snapshot(cells, expected_year=30,
                                   expected_gene_rate_groups=[(20, 0.005)],  # Different
                                   min_cells=50)
        print("   ✗ Should have raised error for wrong gene_rate_groups")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 7: Excessive methylation
    print("\n7. Excessive methylation...")
    cells = create_test_cells(100, age=30, methylation_fraction=0.8)  # Too much
    try:
        validator.validate_snapshot(cells, expected_year=30,
                                   expected_gene_rate_groups=gene_rate_groups,
                                   min_cells=50)
        print("   ✗ Should have raised error for excessive methylation")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    return True


def test_initial_individuals_validation():
    """Test validation of newly created individuals."""
    print("\n" + "="*60)
    print("TEST: Initial Individuals Validation")
    print("="*60)
    
    validator = PipelineValidator(verbose=False)
    gene_rate_groups = [(20, 0.005)]
    
    # Create valid test individuals
    def create_test_individual(individual_id, individual_type, snapshot_year=30):
        cell = create_test_cells(1, age=snapshot_year, gene_rate_groups=gene_rate_groups)[0]
        petri = PetriDish.from_cells(
            cell,
            growth_phase=7,
            calculate_cell_jsds=True,
            metadata={
                'individual_id': individual_id,
                'individual_type': individual_type,
                'source_quantile': 1 if individual_type == 'mutant' else None,
                'source': 'uniform' if individual_type == 'control1' else None
            }
        )
        return petri
    
    # Test 1: Valid individuals
    print("\n1. Valid individuals...")
    mutant_dishes = [create_test_individual(i+1, 'mutant') for i in range(4)]
    control1_dishes = [create_test_individual(i+1, 'control1') for i in range(4)]
    
    try:
        validator.validate_initial_individuals(
            mutant_dishes, control1_dishes,
            expected_count=4,
            expected_gene_rate_groups=gene_rate_groups,
            snapshot_year=30
        )
        print("   ✓ Valid individuals passed")
    except ValidationError as e:
        print(f"   ✗ Unexpected error: {e}")
        return False
    
    # Test 2: Wrong count
    print("\n2. Wrong individual count...")
    try:
        validator.validate_initial_individuals(
            mutant_dishes[:2], control1_dishes,  # Only 2 mutants
            expected_count=4,
            expected_gene_rate_groups=gene_rate_groups,
            snapshot_year=30
        )
        print("   ✗ Should have raised error for wrong count")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 3: Duplicate IDs
    print("\n3. Duplicate individual IDs...")
    mutant_dishes_dup = [create_test_individual(1, 'mutant') for _ in range(4)]  # All ID=1
    try:
        validator.validate_initial_individuals(
            mutant_dishes_dup, control1_dishes,
            expected_count=4,
            expected_gene_rate_groups=gene_rate_groups,
            snapshot_year=30
        )
        print("   ✗ Should have raised error for duplicate IDs")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 4: Multiple founding cells
    print("\n4. Multiple founding cells...")
    bad_petri = create_test_individual(1, 'mutant')
    bad_petri.cells = create_test_cells(2, age=30)  # 2 cells instead of 1
    mutant_dishes_bad = [bad_petri] + mutant_dishes[1:]
    
    try:
        validator.validate_initial_individuals(
            mutant_dishes_bad, control1_dishes,
            expected_count=4,
            expected_gene_rate_groups=gene_rate_groups,
            snapshot_year=30
        )
        print("   ✗ Should have raised error for multiple founding cells")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    return True


def test_grown_individuals_validation():
    """Test validation of grown individuals with extinction handling."""
    print("\n" + "="*60)
    print("TEST: Grown Individuals Validation")
    print("="*60)
    
    # Test with strict config (no extinction allowed)
    config_strict = ValidationConfig.strict()
    validator_strict = PipelineValidator(config=config_strict, verbose=False)
    
    # Test with relaxed config (extinction allowed)
    config_relaxed = ValidationConfig.relaxed()
    validator_relaxed = PipelineValidator(config=config_relaxed, verbose=False)
    
    gene_rate_groups = [(20, 0.005)]
    
    # Create grown individuals
    def create_grown_individual(individual_id, n_cells, with_history=True):
        cells = create_test_cells(n_cells, age=40, gene_rate_groups=gene_rate_groups)
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=gene_rate_groups,
            n=100,
            gene_size=5,
            growth_phase=7
        )
        petri.metadata = {'individual_id': individual_id}
        
        if with_history:
            # Add fake history
            petri.cell_history = {}
            for year in range(11):  # 0-10 years
                petri.cell_history[str(year)] = [
                    {'methylated': [0] * 100} for _ in range(min(2**year, n_cells))
                ]
        
        return petri
    
    # Test 1: Normal population
    print("\n1. Normal population sizes...")
    mutant_dishes = [create_grown_individual(i+1, 100) for i in range(4)]
    control1_dishes = [create_grown_individual(i+1, 110) for i in range(4)]
    
    report = validator_relaxed.validate_grown_individuals(
        mutant_dishes, control1_dishes,
        expected_population=128,
        growth_years=10,
        allow_extinction=True
    )
    
    if report['total_extinct'] == 0:
        print("   ✓ Normal populations validated correctly")
    else:
        print(f"   ✗ Unexpected extinctions: {report}")
        return False
    
    # Test 2: Some extinctions
    print("\n2. Some extinct individuals...")
    mutant_dishes[0] = create_grown_individual(1, 0)  # Extinct
    mutant_dishes[1] = create_grown_individual(2, 10)  # Low pop
    
    report = validator_relaxed.validate_grown_individuals(
        mutant_dishes, control1_dishes,
        expected_population=128,
        growth_years=10,
        allow_extinction=True
    )
    
    if report['total_extinct'] == 1 and len(report['mutant_low_pop']) == 1:
        print(f"   ✓ Correctly identified {report['total_extinct']} extinct, "
              f"{len(report['mutant_low_pop'])} low population")
    else:
        print(f"   ✗ Wrong extinction detection: {report}")
        return False
    
    # Test 3: Complete batch extinction (should fail with strict config)
    print("\n3. Complete batch extinction...")
    all_extinct = [create_grown_individual(i+1, 0) for i in range(4)]
    
    try:
        report = validator_strict.validate_grown_individuals(
            all_extinct, control1_dishes,
            expected_population=128,
            growth_years=10,
            allow_extinction=False  # Strict mode
        )
        print("   ✗ Should have raised error for complete extinction")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 4: High extinction rate
    print("\n4. High extinction rate...")
    mostly_extinct = [create_grown_individual(i+1, 0 if i < 3 else 100) for i in range(4)]
    
    report = validator_relaxed.validate_grown_individuals(
        mostly_extinct, control1_dishes,
        expected_population=128,
        growth_years=10,
        allow_extinction=True
    )
    
    if report['total_extinct'] == 3:
        print(f"   ✓ High extinction rate handled: {report['total_extinct']}/{8} extinct")
    else:
        print(f"   ✗ Wrong extinction count: {report}")
        return False
    
    return True


def test_normalization_validation():
    """Test validation after normalization."""
    print("\n" + "="*60)
    print("TEST: Normalization Validation")
    print("="*60)
    
    validator = PipelineValidator(verbose=False)
    gene_rate_groups = [(20, 0.005)]
    
    def create_normalized_individual(individual_id, n_cells):
        cells = create_test_cells(n_cells, age=40, gene_rate_groups=gene_rate_groups)
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=gene_rate_groups,
            n=100,
            gene_size=5,
            growth_phase=7
        )
        petri.metadata = {
            'individual_id': individual_id,
            'normalized': True,
            'normalization_threshold': n_cells
        }
        return petri
    
    # Test 1: Correctly normalized
    print("\n1. Correctly normalized populations...")
    threshold = 50
    mutant_dishes = [create_normalized_individual(i+1, threshold) for i in range(3)]
    control1_dishes = [create_normalized_individual(i+1, threshold) for i in range(2)]
    
    try:
        validator.validate_normalized_populations(
            mutant_after=mutant_dishes,
            control1_after=control1_dishes,
            threshold=threshold
        )
        print("   ✓ Correctly normalized populations passed")
    except ValidationError as e:
        print(f"   ✗ Unexpected error: {e}")
        return False
    
    # Test 2: Wrong size after normalization
    print("\n2. Wrong size after normalization...")
    mutant_dishes[0].cells = create_test_cells(45, age=40)  # Wrong size
    
    try:
        validator.validate_normalized_populations(
            mutant_after=mutant_dishes,
            control1_after=control1_dishes,
            threshold=threshold
        )
        print("   ✗ Should have raised error for wrong size")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 3: All excluded
    print("\n3. All individuals excluded...")
    try:
        validator.validate_normalized_populations(
            mutant_after=[],
            control1_after=[],
            threshold=threshold
        )
        print("   ✗ Should have raised error for all excluded")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    return True


def test_mixed_populations_validation():
    """Test validation after mixing."""
    print("\n" + "="*60)
    print("TEST: Mixed Populations Validation")
    print("="*60)
    
    validator = PipelineValidator(verbose=False)
    gene_rate_groups = [(20, 0.005)]
    
    def create_mixed_individual(individual_id, n_cells, mix_ratio, single_age=50):
        # After mixing, cells are at the same age (snapshot year)
        # But they came from different sources (some grown from year 40, some from year 50 snapshot)
        cells = create_test_cells(n_cells, age=single_age, gene_rate_groups=gene_rate_groups)
        
        # Mark some cells to simulate they came from different sources 
        # (this would be reflected in their methylation patterns, not ages)
        
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=gene_rate_groups,
            n=100,
            gene_size=5,
            growth_phase=7
        )
        petri.metadata = {
            'individual_id': individual_id,
            'mixed': True,
            'mix_ratio': mix_ratio,
            'final_cells': n_cells
        }
        return petri
    
    # Test 1: Valid mixed populations
    print("\n1. Valid mixed populations...")
    mutant_dishes = [create_mixed_individual(i+1, 640, 80) for i in range(4)]
    control1_dishes = [create_mixed_individual(i+1, 640, 80) for i in range(4)]
    
    try:
        validator.validate_mixed_populations(
            mutant_dishes, control1_dishes,
            mix_ratio=80,
            uniform_mixing=False,
            second_snapshot_year=50
        )
        print("   ✓ Valid mixed populations passed")
    except ValidationError as e:
        print(f"   ✗ Unexpected error: {e}")
        return False
    
    # Test 2: Not marked as mixed
    print("\n2. Individual not marked as mixed...")
    mutant_dishes[0].metadata['mixed'] = False
    
    try:
        validator.validate_mixed_populations(
            mutant_dishes, control1_dishes,
            mix_ratio=80,
            uniform_mixing=False,
            second_snapshot_year=50
        )
        print("   ✗ Should have raised error for not mixed")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Reset for next test
    mutant_dishes[0].metadata['mixed'] = True
    
    # Test 3: Uniform mixing with different sizes
    print("\n3. Uniform mixing with non-uniform sizes...")
    mutant_dishes[1].cells = create_test_cells(500, age=40)  # Different size
    
    validator.validate_mixed_populations(
        mutant_dishes, control1_dishes,
        mix_ratio=80,
        uniform_mixing=True,  # Should all be same size
        second_snapshot_year=50
    )
    
    # Check for warning
    if len(validator.warnings) > 0:
        print(f"   ✓ Warning generated for non-uniform sizes: {len(validator.warnings)} warnings")
    else:
        print("   ✗ Should have generated warning for non-uniform sizes")
        return False
    
    return True


def test_control2_validation():
    """Test control2 validation."""
    print("\n" + "="*60)
    print("TEST: Control2 Validation")
    print("="*60)
    
    validator = PipelineValidator(verbose=False)
    gene_rate_groups = [(20, 0.005)]
    
    # Test 1: Valid control2
    print("\n1. Valid control2 individuals...")
    control2_dishes = []
    for i in range(3):
        cells = create_test_cells(640, age=50, gene_rate_groups=gene_rate_groups)
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=gene_rate_groups,
            n=100,
            gene_size=5,
            growth_phase=None  # Static population
        )
        petri.metadata = {
            'individual_id': i+1,
            'individual_type': 'control2'
        }
        control2_dishes.append(petri)
    
    try:
        validator.validate_control2(
            control2_dishes,
            second_snapshot_year=50,
            expected_count=3
        )
        print("   ✓ Valid control2 passed")
    except ValidationError as e:
        print(f"   ✗ Unexpected error: {e}")
        return False
    
    # Test 2: Wrong growth_phase
    print("\n2. Control2 with growth_phase...")
    control2_dishes[0].growth_phase = 7  # Should be None
    
    try:
        validator.validate_control2(
            control2_dishes,
            second_snapshot_year=50,
            expected_count=3
        )
        print("   ✗ Should have raised error for growth_phase")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Reset
    control2_dishes[0].growth_phase = None
    
    # Test 3: Wrong cell age
    print("\n3. Control2 with wrong cell age...")
    control2_dishes[1].cells = create_test_cells(640, age=40)  # Wrong age
    
    try:
        validator.validate_control2(
            control2_dishes,
            second_snapshot_year=50,
            expected_count=3
        )
        print("   ✗ Should have raised error for wrong age")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    return True


def test_pre_save_validation():
    """Test pre-save validation."""
    print("\n" + "="*60)
    print("TEST: Pre-Save Validation")
    print("="*60)
    
    validator = PipelineValidator(verbose=False)
    gene_rate_groups = [(20, 0.005)]
    
    def create_final_individual(individual_id, individual_type='mutant'):
        cells = create_test_cells(640, age=50, gene_rate_groups=gene_rate_groups)
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=gene_rate_groups,
            n=100,
            gene_size=5,
            growth_phase=7
        )
        petri.metadata = {
            'individual_id': individual_id,
            'individual_type': individual_type,
            'mixed': True,
            'final_cells': 640
        }
        return petri
    
    # Test 1: Valid pre-save
    print("\n1. Valid individuals ready to save...")
    mutant_dishes = [create_final_individual(i+1, 'mutant') for i in range(4)]
    control1_dishes = [create_final_individual(i+5, 'control1') for i in range(4)]  # Different IDs (5-8)
    
    try:
        validator.validate_before_save(mutant_dishes, control1_dishes)
        print("   ✓ Valid pre-save passed")
    except ValidationError as e:
        print(f"   ✗ Unexpected error: {e}")
        return False
    
    # Test 2: Missing metadata
    print("\n2. Missing metadata...")
    bad_petri = create_final_individual(5, 'mutant')
    del bad_petri.metadata
    mutant_dishes_bad = mutant_dishes + [bad_petri]
    
    try:
        validator.validate_before_save(mutant_dishes_bad, control1_dishes)
        print("   ✗ Should have raised error for missing metadata")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 3: Duplicate IDs within batch
    print("\n3. Duplicate individual IDs within batch...")
    control1_dishes[0].metadata['individual_id'] = 6  # Duplicate with control1[1] (also ID 6)
    control1_dishes[1].metadata['individual_id'] = 6
    
    try:
        validator.validate_before_save(mutant_dishes, control1_dishes)
        print("   ✗ Should have raised error for duplicate IDs")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    # Test 4: None values in critical fields
    print("\n4. None values in critical fields...")
    control1_dishes[1].metadata['individual_id'] = None
    
    try:
        validator.validate_before_save(mutant_dishes, control1_dishes)
        print("   ✗ Should have raised error for None value")
        return False
    except ValidationError as e:
        print(f"   ✓ Correctly raised error: {str(e)[:50]}...")
    
    return True


def test_edge_cases():
    """Test various edge cases."""
    print("\n" + "="*60)
    print("TEST: Edge Cases")
    print("="*60)
    
    validator = PipelineValidator(verbose=False)
    gene_rate_groups = [(20, 0.005)]
    
    # Test 1: Empty population (all extinct)
    print("\n1. Completely extinct population...")
    extinct_dishes = []
    for i in range(4):
        petri = PetriDish(
            cells=[],  # No cells
            gene_rate_groups=gene_rate_groups,
            n=100,
            gene_size=5,
            growth_phase=7
        )
        petri.metadata = {'individual_id': i+1, 'individual_type': 'mutant'}
        extinct_dishes.append(petri)
    
    # Should generate warnings but not fail
    validator.warnings = []  # Clear previous warnings
    validator.validate_before_save(extinct_dishes, [])
    
    if len(validator.warnings) > 0:
        print(f"   ✓ Generated {len(validator.warnings)} warnings for extinct individuals")
    else:
        print("   ✗ Should have generated warnings")
        return False
    
    # Test 2: Very large population
    print("\n2. Very large population...")
    large_cells = create_test_cells(10000, age=40, gene_rate_groups=gene_rate_groups)
    large_petri = PetriDish(
        cells=large_cells,
        gene_rate_groups=gene_rate_groups,
        n=100,
        gene_size=5,
        growth_phase=7
    )
    large_petri.metadata = {'individual_id': 1}
    
    report = validator.validate_grown_individuals(
        [large_petri], [],
        expected_population=128,
        growth_years=10
    )
    
    # Should flag as high population
    if len(validator.warnings) > 0:
        print(f"   ✓ Flagged large population: {len(large_petri.cells)} cells")
    else:
        print("   ✗ Should have warned about large population")
        return False
    
    # Test 3: Methylation regression
    print("\n3. Methylation regression detection...")
    petri = PetriDish(
        cells=create_test_cells(100, age=40),
        gene_rate_groups=gene_rate_groups,
        n=100,
        gene_size=5,
        growth_phase=7
    )
    petri.metadata = {'individual_id': 1}
    
    # Create history with decreasing methylation
    petri.cell_history = {
        '0': [{'methylated': [0] * 50 + [1] * 50}],  # 50% methylated
        '1': [{'methylated': [0] * 80 + [1] * 20}],  # 20% methylated (regression!)
        '2': [{'methylated': [0] * 70 + [1] * 30}],  # 30% methylated
    }
    
    validator.warnings = []
    report = validator.validate_grown_individuals(
        [petri], [],
        expected_population=128,
        growth_years=2
    )
    
    if len(validator.warnings) > 0:
        print(f"   ✓ Detected methylation regression: {len(validator.warnings)} warnings")
    else:
        print("   ✗ Should have detected methylation regression")
        return False
    
    return True


def main():
    """Run all validation tests."""
    print("="*60)
    print("VALIDATION SYSTEM TESTS")
    print("="*60)
    
    tests = [
        ("Snapshot validation", test_snapshot_validation),
        ("Initial individuals validation", test_initial_individuals_validation),
        ("Grown individuals validation", test_grown_individuals_validation),
        ("Normalization validation", test_normalization_validation),
        ("Mixed populations validation", test_mixed_populations_validation),
        ("Control2 validation", test_control2_validation),
        ("Pre-save validation", test_pre_save_validation),
        ("Edge cases", test_edge_cases),
    ]
    
    results = []
    for test_name, test_func in tests:
        try:
            passed = test_func()
            results.append((test_name, passed))
        except Exception as e:
            print(f"\n✗ Test '{test_name}' crashed: {e}")
            import traceback
            traceback.print_exc()
            results.append((test_name, False))
    
    # Summary
    print("\n" + "="*60)
    print("TEST SUMMARY")
    print("="*60)
    
    for test_name, passed in results:
        status = "✓ PASSED" if passed else "✗ FAILED"
        print(f"{status}: {test_name}")
    
    passed_count = sum(1 for _, passed in results if passed)
    total_count = len(results)
    
    print(f"\nTotal: {passed_count}/{total_count} tests passed")
    
    if passed_count == total_count:
        print("\n🎉 All validation tests passed!")
        return 0
    else:
        print(f"\n❌ {total_count - passed_count} test(s) failed")
        return 1


if __name__ == "__main__":
    exit(main())