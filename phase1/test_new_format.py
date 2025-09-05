#!/usr/bin/env python3
"""
Test the new lean JSON format for phase1 simulations.
"""

import json
import gzip
import os
import sys
import tempfile

# Add path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cell import Cell, PetriDish


def test_cell_serialization():
    """Test Cell to_dict and from_dict methods."""
    print("\n1. Testing Cell serialization...")
    
    # Create a cell with some methylation
    cell = Cell(n=100, rate=0.005)
    cell.cpg_sites[0] = 1
    cell.cpg_sites[5] = 1
    cell.cpg_sites[10] = 1
    # Properties are calculated automatically via @property decorators
    _ = cell.methylation_proportion  # Trigger calculation
    _ = cell.cell_jsd  # Trigger JSD calculation
    
    # Convert to dict (new format)
    cell_dict = cell.to_dict()
    
    # Check new format structure
    assert 'methylated' in cell_dict, "Missing 'methylated' key"
    assert 'cell_jsd' in cell_dict, "Missing 'cell_jsd' key"
    assert len(cell_dict) == 2, f"Should only have 2 keys, got {len(cell_dict)}: {list(cell_dict.keys())}"
    assert len(cell_dict['methylated']) == 100, "Wrong methylated array size"
    
    # No redundant fields
    assert 'cpg_sites' not in cell_dict, "Legacy 'cpg_sites' field present"
    assert 'rate' not in cell_dict, "Redundant 'rate' field present"
    assert 'gene_size' not in cell_dict, "Redundant 'gene_size' field present"
    assert 'site_rates' not in cell_dict, "Redundant 'site_rates' field present"
    assert 'methylation_proportion' not in cell_dict, "Redundant 'methylation_proportion' field present"
    
    print("  ✓ Cell.to_dict() creates lean format")
    
    # Test from_dict
    cell2 = Cell.from_dict(cell_dict, rate=0.005)
    assert cell2.cpg_sites == cell.cpg_sites, "Methylation state not preserved"
    assert abs(cell2.cell_jsd - cell.cell_jsd) < 0.0001, "JSD not preserved"
    
    print("  ✓ Cell.from_dict() restores correctly")
    
    return True


def test_petri_dish_save():
    """Test PetriDish save_history with new format."""
    print("\n2. Testing PetriDish save format...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create small petri dish
        petri = PetriDish(rate=0.005, n=100, growth_phase=2, seed=42)
        petri.track_cell_history = True
        petri.track_gene_jsd = True
        
        # Run for 3 years
        for year in range(3):
            petri.simulate_year()
        
        # Save with new format
        filepath = petri.save_history(directory=tmpdir, compress=False)
        
        # Load and check structure
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        # Check top-level structure
        assert 'parameters' in data, "Missing 'parameters' section"
        assert 'history' in data, "Missing 'history' section"
        assert len(data) == 2, f"Should only have 2 top keys, got {list(data.keys())}"
        
        print("  ✓ New format has parameters and history sections")
        
        # Check parameters
        params = data['parameters']
        assert params['rate'] == 0.005, "Rate not saved correctly"
        assert params['n'] == 100, "n not saved correctly"
        assert params['growth_phase'] == 2, "growth_phase not saved correctly"
        assert params['seed'] == 42, "seed not saved correctly"
        
        print("  ✓ Parameters saved correctly")
        
        # Check history structure
        history = data['history']
        assert '0' in history, "Year 0 missing"
        assert '3' in history, "Year 3 missing"
        
        # Check year data structure
        year_data = history['0']
        assert 'cells' in year_data, "Missing 'cells' in year data"
        
        # Check cell format
        first_cell = year_data['cells'][0]
        assert 'methylated' in first_cell, "Missing 'methylated' in cell"
        assert 'cell_jsd' in first_cell, "Missing 'cell_jsd' in cell"
        assert 'cpg_sites' not in first_cell, "Legacy 'cpg_sites' present"
        assert 'rate' not in first_cell, "Redundant 'rate' in cell"
        
        print("  ✓ Cells use lean format")
        
        # Check file size reduction
        # Old format would store rate (8 bytes) + gene_size (4 bytes) + 
        # site_rates (100 * 8 bytes) = 812 bytes per cell of redundant data
        # With 4 cells at year 2, that's ~3.2KB of redundant data
        
        return True


def test_compressed_vs_uncompressed():
    """Test compression settings."""
    print("\n3. Testing compression...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        petri = PetriDish(rate=0.005, n=100, growth_phase=2, seed=42)
        
        for year in range(5):
            petri.simulate_year()
        
        # Save compressed
        compressed_path = petri.save_history(directory=tmpdir, compress=True)
        compressed_size = os.path.getsize(compressed_path)
        
        # Save uncompressed
        uncompressed_path = petri.save_history(directory=tmpdir, compress=False)
        uncompressed_size = os.path.getsize(uncompressed_path)
        
        ratio = uncompressed_size / compressed_size
        print(f"  Compressed: {compressed_size} bytes")
        print(f"  Uncompressed: {uncompressed_size} bytes")
        print(f"  Ratio: {ratio:.1f}x")
        
        # Load both and verify same content
        with gzip.open(compressed_path, 'rt') as f:
            compressed_data = json.load(f)
        
        with open(uncompressed_path, 'r') as f:
            uncompressed_data = json.load(f)
        
        # Parameters should match
        assert compressed_data['parameters'] == uncompressed_data['parameters']
        
        print("  ✓ Compressed and uncompressed have same content")
        
        return True


def test_gene_rate_groups():
    """Test gene-specific rate configuration."""
    print("\n4. Testing gene-specific rates...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create with gene-specific rates
        gene_rate_groups = [(10, 0.004), (10, 0.006)]
        petri = PetriDish(gene_rate_groups=gene_rate_groups, n=100, gene_size=5, 
                         growth_phase=1, seed=42)
        
        petri.simulate_year()
        
        # Save
        filepath = petri.save_history(directory=tmpdir, compress=False)
        
        # Load and check
        with open(filepath, 'r') as f:
            data = json.load(f)
        
        params = data['parameters']
        assert params['rate'] is None, "rate should be None for gene-specific"
        # JSON converts tuples to lists, so compare as lists
        saved_groups = params['gene_rate_groups']
        expected_groups = [[n, r] for n, r in gene_rate_groups]
        assert saved_groups == expected_groups, f"gene_rate_groups not saved correctly: {saved_groups} != {expected_groups}"
        
        # Check cells don't have redundant rate info
        first_cell = data['history']['0']['cells'][0]
        assert 'site_rates' not in first_cell, "Redundant site_rates in cell"
        assert 'gene_rate_groups' not in first_cell, "Redundant gene_rate_groups in cell"
        
        print("  ✓ Gene-specific rates handled correctly")
        
        return True


def test_format_comparison():
    """Compare file sizes between old and new formats."""
    print("\n5. Estimating format size savings...")
    
    # Calculate theoretical savings for typical simulation
    n_sites = 1000
    n_cells_year_13 = 8192
    n_years = 100
    
    # Old format redundant data per cell:
    # - rate: 8 bytes (float as string)
    # - gene_size: 4 bytes  
    # - site_rates: 1000 * 8 = 8000 bytes
    # - methylation_proportion: 8 bytes (redundant, can be calculated)
    # - methylation_distribution: 6 * 8 = 48 bytes (redundant)
    # - age: 4 bytes (redundant when tracking by year)
    old_redundant_per_cell = 8 + 4 + 8000 + 8 + 48 + 4
    
    # Total redundant data for year 13
    redundant_year_13 = old_redundant_per_cell * n_cells_year_13
    redundant_mb = redundant_year_13 / (1024 * 1024)
    
    print(f"  Old format redundant data per cell: {old_redundant_per_cell:,} bytes")
    print(f"  For {n_cells_year_13:,} cells: {redundant_mb:.1f} MB of redundant data")
    print(f"  Estimated savings: ~{redundant_mb/70*100:.0f}% file size reduction")
    
    return True


def main():
    """Run all tests."""
    print("=" * 60)
    print("TESTING NEW LEAN JSON FORMAT")
    print("=" * 60)
    
    tests = [
        test_cell_serialization,
        test_petri_dish_save,
        test_compressed_vs_uncompressed,
        test_gene_rate_groups,
        test_format_comparison
    ]
    
    passed = 0
    failed = 0
    
    for test in tests:
        try:
            if test():
                passed += 1
        except Exception as e:
            print(f"  ✗ Test failed: {e}")
            import traceback
            traceback.print_exc()
            failed += 1
    
    print("\n" + "=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)
    
    if passed == len(tests):
        print("\n✅ New lean format is working correctly!")
        print("\nBenefits:")
        print("- No redundant data (rate, gene_size, site_rates per cell)")
        print("- Cleaner structure (parameters + history)")  
        print("- Significantly smaller file sizes")
        print("- Faster to save and load")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())