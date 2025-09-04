#!/usr/bin/env python3
"""
Test static PetriDish functionality (growth_phase=None).
Verifies that control2 individuals are created correctly without fake metadata.
"""

import sys
import os
import json
import tempfile

# Add paths for imports
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from core.pipeline_utils import create_pure_snapshot_petri, save_petri_dish, load_petri_dish


def test_static_petridish_creation():
    """Test creating a static PetriDish with growth_phase=None."""
    print("\n" + "="*60)
    print("TEST: Static PetriDish Creation")
    print("="*60)
    
    # Create cells with specific age
    cells = []
    for i in range(10):
        cell = Cell(n=100, gene_rate_groups=[(20, 0.005)], gene_size=5)
        cell.age = 50  # Simulate year 50 cells
        # Add some methylation
        for j in range(25):
            cell.cpg_sites[j] = 1
        cells.append(cell)
    
    # Create static PetriDish
    try:
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=[(20, 0.005)],
            n=100,
            gene_size=5,
            growth_phase=None,  # Static population
            calculate_cell_jsds=True
        )
        print("✓ Created static PetriDish with growth_phase=None")
        print(f"  PetriDish year: {petri.year} (should be 50)")
        print(f"  Growth phase: {petri.growth_phase}")
        print(f"  Target population: {petri.target_population}")
        
        # Verify year matches cell ages
        if petri.year == 50:
            print("✓ PetriDish year correctly set from cell ages")
        else:
            print(f"✗ Wrong year: expected 50, got {petri.year}")
            return False
            
        # Check cell_history key
        if hasattr(petri, 'cell_history') and petri.cell_history:
            history_keys = list(petri.cell_history.keys())
            print(f"  Cell history keys: {history_keys}")
            if '50' in history_keys:
                print("✓ Cell history uses correct year as key")
            else:
                print(f"✗ Cell history has wrong key: {history_keys}")
                return False
                
        return True
        
    except Exception as e:
        print(f"✗ Failed to create static PetriDish: {e}")
        return False


def test_static_petridish_cannot_age():
    """Test that static PetriDish cannot be aged."""
    print("\n" + "="*60)
    print("TEST: Static PetriDish Cannot Be Aged")
    print("="*60)
    
    # Create a static PetriDish
    cells = [Cell(n=100, gene_rate_groups=[(20, 0.005)], gene_size=5) for _ in range(5)]
    for cell in cells:
        cell.age = 50
    
    petri = PetriDish(
        cells=cells,
        gene_rate_groups=[(20, 0.005)],
        n=100,
        gene_size=5,
        growth_phase=None
    )
    
    # Try to age it
    try:
        petri.simulate_year()
        print("✗ Should have raised error when aging static PetriDish")
        return False
    except ValueError as e:
        if "static PetriDish" in str(e):
            print("✓ Correctly prevented aging of static PetriDish")
            print(f"  Error: {str(e)[:100]}...")
            return True
        else:
            print(f"✗ Wrong error: {e}")
            return False


def test_control2_creation():
    """Test control2 creation with correct metadata."""
    print("\n" + "="*60)
    print("TEST: Control2 Creation and Metadata")
    print("="*60)
    
    # Create snapshot cells at year 50
    snapshot_cells = []
    for i in range(100):
        cell = Cell(n=100, gene_rate_groups=[(20, 0.005)], gene_size=5)
        cell.age = 50
        for j in range(30):
            cell.cpg_sites[j] = 1
        cell.cell_JSD = 0.25
        snapshot_cells.append(cell)
    
    # Create control2 using pure_snapshot_petri
    try:
        control2 = create_pure_snapshot_petri(snapshot_cells, n_cells=50, seed=42)
        
        print(f"✓ Created control2 PetriDish")
        print(f"  Year: {control2.year} (should be 50)")
        print(f"  Growth phase: {control2.growth_phase} (should be None)")
        print(f"  Number of cells: {len(control2.cells)}")
        
        # Check metadata
        if hasattr(control2, 'metadata'):
            print(f"  Metadata keys: {list(control2.metadata.keys())}")
            
        # Save and reload to check JSON structure
        with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
            temp_path = f.name
            
        save_petri_dish(control2, temp_path, include_cell_history=False, 
                       include_gene_metrics=True, compress=False)
        
        # Read the JSON directly
        with open(temp_path, 'r') as f:
            data = json.load(f)
        
        print("\nJSON structure:")
        print(f"  Outer keys: {list(data.keys())}")
        if 'metadata' in data:
            print(f"  Metadata keys: {list(data['metadata'].keys())}")
            print(f"    year: {data['metadata'].get('year')} (should be 50)")
            print(f"    growth_phase: {data['metadata'].get('growth_phase')} (should be null)")
            
            # Check for absence of redundant fields
            if 'initial_year' in data['metadata']:
                print(f"  ✗ Found redundant 'initial_year' field")
                return False
            else:
                print(f"  ✓ No redundant 'initial_year' field")
                
        # Clean up
        os.unlink(temp_path)
        
        # Verify correct values
        if control2.year == 50 and control2.growth_phase is None:
            print("\n✓ Control2 has correct year and growth_phase")
            return True
        else:
            print(f"\n✗ Wrong values: year={control2.year}, growth_phase={control2.growth_phase}")
            return False
            
    except Exception as e:
        print(f"✗ Control2 creation failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def test_mixed_age_cells_rejected():
    """Test that cells with different ages are rejected."""
    print("\n" + "="*60)
    print("TEST: Mixed Age Cells Rejected")
    print("="*60)
    
    # Create cells with different ages
    cells = []
    for i in range(5):
        cell = Cell(n=100, gene_rate_groups=[(20, 0.005)], gene_size=5)
        cell.age = 30  # First 5 cells age 30
        cells.append(cell)
    
    for i in range(5):
        cell = Cell(n=100, gene_rate_groups=[(20, 0.005)], gene_size=5)
        cell.age = 50  # Next 5 cells age 50
        cells.append(cell)
    
    # Try to create PetriDish
    try:
        petri = PetriDish(
            cells=cells,
            gene_rate_groups=[(20, 0.005)],
            n=100,
            gene_size=5
        )
        print("✗ Should have rejected cells with different ages")
        return False
    except ValueError as e:
        if "same age" in str(e):
            print("✓ Correctly rejected cells with different ages")
            print(f"  Error: {str(e)}")
            return True
        else:
            print(f"✗ Wrong error: {e}")
            return False


def main():
    """Run all tests."""
    print("="*60)
    print("STATIC PETRIDISH TESTS")
    print("="*60)
    
    tests = [
        ("Static PetriDish creation", test_static_petridish_creation),
        ("Static PetriDish cannot age", test_static_petridish_cannot_age),
        ("Control2 creation and metadata", test_control2_creation),
        ("Mixed age cells rejected", test_mixed_age_cells_rejected),
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
        print("\n🎉 All tests passed!")
        return 0
    else:
        print(f"\n❌ {total_count - passed_count} test(s) failed")
        return 1


if __name__ == "__main__":
    exit(main())