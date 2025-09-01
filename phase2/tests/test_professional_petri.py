#!/usr/bin/env python3
"""
Comprehensive test suite for the professional PetriDish implementation.
Tests the new from_snapshot_cell() method, metadata management, and save/load cycle.
"""

import os
import sys
import json
import tempfile
import traceback

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish, load_petri_dish


def test_from_snapshot_cell():
    """Test the new from_snapshot_cell() class method."""
    print("\nTest 1: from_snapshot_cell() method")
    print("=" * 50)
    
    # Create a snapshot cell with some methylation
    cell = Cell(n=100, rate=0.005)
    for _ in range(10):
        cell.methylate()
    
    # Test metadata
    metadata = {
        'individual_id': 42,
        'individual_type': 'test',
        'source': 'unit_test'
    }
    
    # Create PetriDish using professional method
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=7,
        calculate_cell_jsds=True,
        metadata=metadata
    )
    
    # Verify properties
    assert len(petri.cells) == 1, f"Should have 1 cell, got {len(petri.cells)}"
    assert petri.cells[0].n == 100, f"Cell should have 100 sites"
    assert petri.cells[0].rate == 0.005, f"Cell rate should be 0.005"
    assert petri.growth_phase == 7, f"Growth phase should be 7"
    
    # Verify metadata
    assert hasattr(petri, 'metadata'), "PetriDish should have metadata"
    assert petri.metadata['individual_id'] == 42, "Individual ID mismatch"
    assert petri.metadata['individual_type'] == 'test', "Individual type mismatch"
    assert petri.metadata['source'] == 'unit_test', "Source mismatch"
    assert petri.metadata['num_cells'] == 1, "num_cells should be 1"
    assert petri.metadata['creation_method'] == 'from_snapshot_cell', "Creation method not set"
    
    print("  ✅ from_snapshot_cell() works correctly")
    return True


def test_metadata_management():
    """Test metadata management methods."""
    print("\nTest 2: Metadata Management")
    print("=" * 50)
    
    # Create a simple PetriDish
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(cell, growth_phase=5)
    
    # Test update_metadata
    print("  Testing update_metadata()...")
    petri.update_metadata({
        'experiment': 'test_run',
        'batch': 3,
        'custom_field': 'value'
    })
    
    assert petri.metadata['experiment'] == 'test_run', "Experiment not updated"
    assert petri.metadata['batch'] == 3, "Batch not updated"
    assert petri.metadata['num_cells'] == 1, "num_cells should stay at 1"
    
    # Grow the PetriDish
    petri.divide_cells()
    petri.update_metadata({'after_division': True})
    
    assert petri.metadata['num_cells'] == 2, "num_cells should update to 2 after division"
    assert petri.metadata['after_division'] == True, "New field not added"
    
    print("  ✅ update_metadata() maintains num_cells correctly")
    
    # Test get_current_stats
    print("  Testing get_current_stats()...")
    stats = petri.get_current_stats()
    
    assert stats['num_cells'] == 2, "Stats should show 2 cells"
    assert stats['year'] == 0, "Year should be 0"
    assert 'mean_methylation' in stats, "Should have mean_methylation"
    assert 'mean_cell_jsd' in stats, "Should have mean_cell_jsd"
    
    # Should have gene metrics if calculate_cell_jsds is True
    assert 'gene_jsds' in stats, "Should have gene_jsds"
    assert 'gene_mean_methylation' in stats, "Should have gene_mean_methylation"
    assert stats['n_genes'] == 20, "Should have 20 genes (100 sites / 5 sites per gene)"
    
    print("  ✅ get_current_stats() returns accurate data")
    
    # Test merge_metadata
    print("  Testing merge_metadata()...")
    
    # Set some initial values
    petri.metadata['gene_jsds'] = [0.1, 0.2, 0.3]
    petri.metadata['to_override'] = 'old_value'
    
    # Merge with new metadata
    new_metadata = {
        'to_override': 'new_value',
        'gene_jsds': [0.9, 0.9, 0.9],  # Should be preserved
        'new_field': 'added'
    }
    
    petri.merge_metadata(new_metadata)
    
    assert petri.metadata['to_override'] == 'new_value', "Field should be overridden"
    assert petri.metadata['gene_jsds'] == [0.1, 0.2, 0.3], "gene_jsds should be preserved"
    assert petri.metadata['new_field'] == 'added', "New field should be added"
    assert petri.metadata['num_cells'] == 2, "num_cells should stay accurate"
    
    print("  ✅ merge_metadata() preserves critical fields")
    
    return True


def test_prepare_for_save():
    """Test the prepare_for_save() method."""
    print("\nTest 3: prepare_for_save() method")
    print("=" * 50)
    
    # Create PetriDish with some cells
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=7,
        metadata={'test_id': 123}
    )
    
    # Grow to get multiple cells
    for _ in range(3):
        petri.divide_cells()
        petri.methylate_cells()
    
    # Test prepare_for_save
    save_data = petri.prepare_for_save(include_gene_metrics=True)
    
    # Verify structure
    assert 'cells' in save_data, "Should have cells array"
    assert 'metadata' in save_data, "Should have metadata"
    
    # Verify cell count
    assert len(save_data['cells']) == 8, "Should have 8 cells (2^3)"
    assert save_data['metadata']['num_cells'] == 8, "Metadata should show 8 cells"
    
    # Verify metadata contents
    metadata = save_data['metadata']
    assert metadata['test_id'] == 123, "Custom metadata preserved"
    assert metadata['year'] == 0, "Year should be 0"
    assert metadata['growth_phase'] == 7, "Growth phase preserved"
    assert metadata['rate'] == 0.005, "Rate preserved"
    
    # Verify gene metrics included
    assert 'gene_jsds' in metadata, "Should have gene JSDs"
    assert 'gene_mean_methylation' in metadata, "Should have gene means"
    assert len(metadata['gene_jsds']) == 20, "Should have 20 gene JSDs"
    assert len(metadata['gene_mean_methylation']) == 20, "Should have 20 gene means"
    
    print("  ✅ prepare_for_save() creates correct structure")
    return True


def test_save_load_cycle():
    """Test saving and loading with new professional methods."""
    print("\nTest 4: Save/Load Cycle")
    print("=" * 50)
    
    # Create a PetriDish with custom metadata
    cell = Cell(n=100, rate=0.005)
    original_petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=6,
        metadata={
            'individual_id': 99,
            'individual_type': 'mutant',
            'source_quantile': 0.95
        }
    )
    
    # Grow and add some history
    original_petri.enable_history_tracking(start_year=50, track_gene_jsd=True)
    for i in range(3):
        original_petri.divide_cells()
        original_petri.methylate_cells()
        original_petri.increment_year()
    
    # Save to temporary file
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        # Save with all features
        save_petri_dish(
            original_petri, 
            temp_file,
            include_cell_history=True,
            include_gene_jsd=True,
            include_gene_metrics=True,
            compress=False
        )
        
        # Verify file contents
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        # Check metadata
        assert data['metadata']['num_cells'] == 8, "Should have 8 cells"
        assert data['metadata']['individual_id'] == 99, "Individual ID lost"
        assert data['metadata']['individual_type'] == 'mutant', "Type lost"
        assert data['metadata']['source_quantile'] == 0.95, "Quantile lost"
        
        # Check gene metrics
        assert 'gene_jsds' in data['metadata'], "Gene JSDs missing"
        assert 'gene_mean_methylation' in data['metadata'], "Gene means missing"
        
        # Check histories
        assert 'cell_history' in data, "Cell history missing"
        assert 'gene_jsd_history' in data, "Gene JSD history missing"
        
        print("  ✅ Save preserves all metadata and metrics")
        
        # Load and verify
        loaded_petri = load_petri_dish(temp_file, include_cell_history=True, include_gene_jsd=True)
        
        assert len(loaded_petri.cells) == 8, "Loaded wrong number of cells"
        assert loaded_petri.year == 3, "Year not preserved"
        assert hasattr(loaded_petri, 'metadata'), "Metadata not loaded"
        assert loaded_petri.metadata['individual_id'] == 99, "Metadata not preserved"
        
        print("  ✅ Load restores PetriDish correctly")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def test_gene_specific_rates():
    """Test from_snapshot_cell with gene-specific rates."""
    print("\nTest 5: Gene-Specific Rates")
    print("=" * 50)
    
    # Create cell with gene-specific rates
    gene_rate_groups = [(10, 0.001), (10, 0.01)]  # 2 groups, 10 genes each
    cell = Cell(n=100, gene_rate_groups=gene_rate_groups)
    
    # Create PetriDish
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=5,
        metadata={'rate_type': 'gene_specific'}
    )
    
    # Verify configuration
    assert petri.rate is None, "Should not have uniform rate"
    assert petri.gene_rate_groups == gene_rate_groups, "Gene rate groups not preserved"
    assert petri.cells[0].gene_rate_groups == gene_rate_groups, "Cell should have gene rates"
    
    # Test that it grows correctly
    petri.divide_cells()
    assert len(petri.cells) == 2, "Should have 2 cells after division"
    
    for cell in petri.cells:
        assert cell.gene_rate_groups == gene_rate_groups, "Daughter cells should inherit gene rates"
    
    print("  ✅ Gene-specific rates handled correctly")
    return True


def test_edge_cases():
    """Test edge cases and error conditions."""
    print("\nTest 6: Edge Cases")
    print("=" * 50)
    
    # Test 1: Empty metadata
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(cell, growth_phase=5, metadata=None)
    
    assert hasattr(petri, 'metadata'), "Should create metadata even if None provided"
    assert petri.metadata['num_cells'] == 1, "Should set num_cells"
    assert petri.metadata['creation_method'] == 'from_snapshot_cell', "Should set creation method"
    print("  ✅ Handles None metadata")
    
    # Test 2: Large cell population
    petri = PetriDish.from_snapshot_cell(cell, growth_phase=10)
    
    # Grow to target
    for _ in range(10):
        petri.divide_cells()
    
    assert len(petri.cells) == 1024, "Should have 2^10 cells"
    
    stats = petri.get_current_stats()
    assert stats['num_cells'] == 1024, "Stats should show correct count"
    
    # Verify metadata stays in sync
    petri.update_metadata({'large_pop': True})
    assert petri.metadata['num_cells'] == 1024, "num_cells should stay accurate"
    
    print("  ✅ Handles large populations")
    
    # Test 3: Zero cells edge case
    petri = PetriDish.from_snapshot_cell(cell, growth_phase=5)
    petri.cells = []  # Artificially empty
    
    stats = petri.get_current_stats()
    assert stats['num_cells'] == 0, "Should handle empty dish"
    assert stats['mean_cell_jsd'] == 0.0, "Should return 0 for empty dish"
    
    petri.update_metadata({'empty_test': True})
    assert petri.metadata['num_cells'] == 0, "Should update to 0"
    
    print("  ✅ Handles empty dish edge case")
    
    # Test 4: Metadata override protection
    petri = PetriDish.from_snapshot_cell(cell, growth_phase=5)
    petri.metadata['num_cells'] = 999  # Try to set wrong value
    
    petri.update_metadata({})  # This should fix it
    assert petri.metadata['num_cells'] == 1, "Should auto-correct num_cells"
    
    print("  ✅ Auto-corrects incorrect num_cells")
    
    return True


def test_enforcement():
    """Test that only professional methods work correctly."""
    print("\nTest 7: Professional Method Enforcement")
    print("=" * 50)
    
    # Create PetriDish professionally
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(cell, growth_phase=5)
    
    # Test that prepare_for_save always exists
    assert hasattr(petri, 'prepare_for_save'), "Should have prepare_for_save"
    assert hasattr(petri, 'update_metadata'), "Should have update_metadata"
    assert hasattr(petri, 'get_current_stats'), "Should have get_current_stats"
    
    # Test that direct manipulation still gets corrected
    petri.cells = [cell, cell, cell]  # Direct assignment (not recommended but possible)
    
    # Any metadata operation should fix the count
    petri.update_metadata({})
    assert petri.metadata['num_cells'] == 3, "Should auto-correct after direct assignment"
    
    # Test saving always works
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        assert data['metadata']['num_cells'] == 3, "Save should have correct count"
        assert len(data['cells']) == 3, "Should have 3 cells"
        
        print("  ✅ Professional methods enforced correctly")
        print("  ✅ Auto-correction works even with direct manipulation")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def test_incremental_updates():
    """Test that metadata updates correctly through pipeline stages."""
    print("\nTest 8: Incremental Pipeline Updates")
    print("=" * 50)
    
    # Simulate pipeline stages
    
    # Stage 1: Create individual
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=7,
        metadata={
            'individual_id': 5,
            'individual_type': 'mutant',
            'source_quantile': 0.75,
            'initial_year': 50
        }
    )
    
    assert petri.metadata['num_cells'] == 1, "Should start with 1 cell"
    print("  Stage 1: Created individual with 1 cell")
    
    # Stage 2: Growth
    for _ in range(7):
        petri.divide_cells()
        petri.methylate_cells()
    
    petri.update_metadata({'growth_complete': True})
    assert petri.metadata['num_cells'] == 128, "Should have 128 cells after growth"
    assert petri.metadata['growth_complete'] == True, "Growth flag set"
    print(f"  Stage 2: Grown to {petri.metadata['num_cells']} cells")
    
    # Stage 3: Mixing (simulate adding cells)
    mix_cells = [Cell(n=100, rate=0.005) for _ in range(72)]
    petri.cells.extend(mix_cells)
    
    petri.update_metadata({
        'mixed': True,
        'mix_ratio': 80,
        'mix_mode': 'test'
    })
    
    assert petri.metadata['num_cells'] == 200, "Should have 200 cells after mixing"
    assert petri.metadata['mixed'] == True, "Mix flag set"
    print(f"  Stage 3: Mixed to {petri.metadata['num_cells']} cells")
    
    # Verify all metadata preserved
    assert petri.metadata['individual_id'] == 5, "ID preserved"
    assert petri.metadata['individual_type'] == 'mutant', "Type preserved"
    assert petri.metadata['source_quantile'] == 0.75, "Quantile preserved"
    assert petri.metadata['initial_year'] == 50, "Year preserved"
    
    print("  ✅ Metadata preserved through all stages")
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("Professional PetriDish Implementation Tests")
    print("=" * 60)
    
    tests = [
        test_from_snapshot_cell,
        test_metadata_management,
        test_prepare_for_save,
        test_save_load_cycle,
        test_gene_specific_rates,
        test_edge_cases,
        test_enforcement,
        test_incremental_updates
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"  ❌ Test failed: {e}")
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Test Summary: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("🎉 All tests passed!")
        return 0
    else:
        print(f"⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())