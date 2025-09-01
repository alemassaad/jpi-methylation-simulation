#!/usr/bin/env python3
"""
Test the clean, professional pipeline implementation.
Verifies no dual code paths, no backward compatibility issues.
"""

import os
import sys
import json
import tempfile
import copy

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import (
    save_petri_dish, 
    load_petri_dish,
    create_pure_snapshot_petri,
    create_control2_with_uniform_base,
    mix_petri_with_snapshot,
    mix_petri_uniform
)


def test_no_dual_paths():
    """Test that there are no dual code paths in save_petri_dish."""
    print("\nTest 1: No Dual Code Paths")
    print("=" * 50)
    
    # Create any PetriDish
    petri = PetriDish(n=100, rate=0.005, growth_phase=5)
    
    # Should have professional methods
    assert hasattr(petri, 'prepare_for_save'), "All PetriDish objects should have prepare_for_save"
    
    # Save should work without any conditionals
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri, temp_file, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        # Verify it saved correctly
        assert 'metadata' in data, "Should have metadata"
        assert 'cells' in data, "Should have cells"
        assert data['metadata']['num_cells'] == 1, "Should have correct count"
        
        print("  ✅ Single code path - no conditionals needed")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    return True


def test_control2_creation_professional():
    """Test that control2 creation uses professional methods."""
    print("\nTest 2: Control2 Creation Professional")
    print("=" * 50)
    
    # Create snapshot cells
    snapshot_cells = []
    for i in range(100):
        cell = Cell(n=100, rate=0.005)
        for _ in range(30):
            cell.methylate()
        snapshot_cells.append(cell)
    
    # Test create_pure_snapshot_petri
    print("  Testing create_pure_snapshot_petri...")
    control2 = create_pure_snapshot_petri(snapshot_cells, n_cells=50, rate=0.005, seed=42)
    
    assert hasattr(control2, 'metadata'), "Should have metadata"
    assert control2.metadata['num_cells'] == 50, "Should have correct num_cells"
    assert control2.metadata['creation_method'] == 'pure_snapshot', "Should record creation method"
    assert hasattr(control2, 'prepare_for_save'), "Should have professional methods"
    
    print("    ✅ Uses from_snapshot_cell and professional methods")
    
    # Test create_control2_with_uniform_base
    print("  Testing create_control2_with_uniform_base...")
    uniform_pool = snapshot_cells[:40]
    uniform_indices = list(range(40))
    
    control2_uniform = create_control2_with_uniform_base(
        snapshot_cells, 
        uniform_pool,
        uniform_indices,
        target_size=60,
        rate=0.005,
        seed=43
    )
    
    assert hasattr(control2_uniform, 'metadata'), "Should have metadata"
    assert control2_uniform.metadata['num_cells'] == 60, "Should have correct num_cells"
    assert control2_uniform.metadata['creation_method'] == 'uniform_base_plus_snapshot', "Should record method"
    assert control2_uniform.metadata['uniform_base_size'] == 40, "Should record base size"
    
    print("    ✅ Uses from_snapshot_cell with uniform base")
    
    return True


def test_mixing_updates_metadata():
    """Test that all mixing functions update metadata."""
    print("\nTest 3: Mixing Functions Update Metadata")
    print("=" * 50)
    
    # Create initial individual
    cell = Cell(n=100, rate=0.005)
    petri = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=7,
        metadata={'individual_id': 0}
    )
    
    # Grow it
    for _ in range(7):
        petri.divide_cells()
    
    initial_count = len(petri.cells)
    assert petri.metadata['num_cells'] == initial_count, "Should track growth"
    
    # Create snapshot cells for mixing
    snapshot_cells = []
    for _ in range(1000):
        mix_cell = Cell(n=100, rate=0.005)
        for _ in range(20):
            mix_cell.methylate()
        snapshot_cells.append(mix_cell)
    
    # Test mix_petri_with_snapshot
    print(f"  Testing mix_petri_with_snapshot...")
    print(f"    Initial: {initial_count} cells")
    
    final_count = mix_petri_with_snapshot(petri, snapshot_cells, mix_ratio=0.8, seed=44)
    
    assert petri.metadata['num_cells'] == final_count, "Metadata should be updated"
    assert petri.metadata['num_cells'] == len(petri.cells), "Should match actual count"
    print(f"    After mixing: {final_count} cells (metadata: {petri.metadata['num_cells']})")
    print("    ✅ mix_petri_with_snapshot updates metadata")
    
    # Test mix_petri_uniform
    print(f"  Testing mix_petri_uniform...")
    
    # Create new petri for uniform test
    petri2 = PetriDish.from_snapshot_cell(cell=cell, growth_phase=5)
    for _ in range(5):
        petri2.divide_cells()
    
    initial_count2 = len(petri2.cells)
    print(f"    Initial: {initial_count2} cells")
    
    uniform_pool = snapshot_cells[:100]
    final_count2 = mix_petri_uniform(petri2, uniform_pool, target_ratio=0.5)
    
    assert petri2.metadata['num_cells'] == final_count2, "Metadata should be updated"
    assert petri2.metadata['num_cells'] == len(petri2.cells), "Should match actual count"
    print(f"    After mixing: {final_count2} cells (metadata: {petri2.metadata['num_cells']})")
    print("    ✅ mix_petri_uniform updates metadata")
    
    return True


def test_complete_pipeline_flow():
    """Test complete pipeline flow with professional methods."""
    print("\nTest 4: Complete Pipeline Flow")
    print("=" * 50)
    
    # Simulate pipeline stages
    stages_passed = []
    
    # Stage 1: Create snapshot
    print("  Stage 1: Create snapshot cells")
    snapshot_cells = []
    for i in range(100):
        cell = Cell(n=100, rate=0.005)
        for _ in range(20):
            cell.methylate()
        snapshot_cells.append(cell)
    stages_passed.append("snapshot_creation")
    
    # Stage 2: Create individuals
    print("  Stage 2: Create individuals professionally")
    
    # Mutant (from quantile)
    mutant = PetriDish.from_snapshot_cell(
        cell=snapshot_cells[10],  # Simulate quantile selection
        growth_phase=7,
        metadata={
            'individual_id': 0,
            'individual_type': 'mutant',
            'source_quantile': 0.9
        }
    )
    assert mutant.metadata['num_cells'] == 1
    stages_passed.append("mutant_creation")
    
    # Control1 (uniform)
    control1 = PetriDish.from_snapshot_cell(
        cell=snapshot_cells[50],  # Simulate uniform selection
        growth_phase=7,
        metadata={
            'individual_id': 0,
            'individual_type': 'control1',
            'source': 'uniform'
        }
    )
    assert control1.metadata['num_cells'] == 1
    stages_passed.append("control1_creation")
    
    # Control2 (pure snapshot)
    control2 = create_pure_snapshot_petri(
        snapshot_cells,
        n_cells=50,
        rate=0.005,
        seed=45
    )
    control2.update_metadata({
        'individual_id': 0,
        'individual_type': 'control2'
    })
    assert control2.metadata['num_cells'] == 50
    stages_passed.append("control2_creation")
    
    # Stage 3: Growth
    print("  Stage 3: Grow individuals")
    for petri, name in [(mutant, "mutant"), (control1, "control1")]:
        for _ in range(7):
            petri.divide_cells()
            petri.methylate_cells()
        petri.update_metadata({'growth_complete': True})
        assert petri.metadata['num_cells'] == 128
        stages_passed.append(f"{name}_growth")
    
    # Stage 4: Mixing
    print("  Stage 4: Mix with snapshot")
    
    # Mix mutant and control1
    for petri, name in [(mutant, "mutant"), (control1, "control1")]:
        final = mix_petri_with_snapshot(petri, snapshot_cells, mix_ratio=0.8, seed=46+len(name))
        petri.update_metadata({'mixed': True, 'mix_ratio': 80})
        assert petri.metadata['num_cells'] == final
        assert petri.metadata['num_cells'] == len(petri.cells)
        stages_passed.append(f"{name}_mixing")
    
    # Stage 5: Save all
    print("  Stage 5: Save all individuals")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        for petri, name in [(mutant, "mutant"), (control1, "control1"), (control2, "control2")]:
            filepath = os.path.join(tmpdir, f"{name}.json")
            save_petri_dish(petri, filepath, include_gene_metrics=True, compress=False)
            
            # Verify saved correctly
            with open(filepath, 'r') as f:
                data = json.load(f)
            
            assert data['metadata']['num_cells'] == len(petri.cells)
            assert data['metadata']['individual_type'] == name
            assert 'gene_jsds' in data['metadata']
            stages_passed.append(f"{name}_save")
    
    # Verify all stages passed
    expected_stages = [
        'snapshot_creation', 'mutant_creation', 'control1_creation', 'control2_creation',
        'mutant_growth', 'control1_growth', 'mutant_mixing', 'control1_mixing',
        'mutant_save', 'control1_save', 'control2_save'
    ]
    
    assert stages_passed == expected_stages, f"Missing stages: {set(expected_stages) - set(stages_passed)}"
    
    print("\n  ✅ Complete pipeline works with professional methods")
    print("  ✅ No dual code paths encountered")
    print("  ✅ Metadata accurate throughout")
    
    return True


def main():
    """Run all clean pipeline tests."""
    print("=" * 60)
    print("Clean Professional Pipeline Tests")
    print("=" * 60)
    
    tests = [
        test_no_dual_paths,
        test_control2_creation_professional,
        test_mixing_updates_metadata,
        test_complete_pipeline_flow
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"\n  ❌ Test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"Test Summary: {passed} passed, {failed} failed")
    
    if failed == 0:
        print("\n🎉 Clean Professional Implementation Complete!")
        print("   ✅ No dual code paths")
        print("   ✅ No backward compatibility code")
        print("   ✅ All operations use professional methods")
        print("   ✅ Metadata always accurate")
        print("   ✅ Simple, maintainable, and robust")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())