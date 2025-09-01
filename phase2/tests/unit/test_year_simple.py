#!/usr/bin/env python3
"""
Simple test for year tracking without absolute_year.
"""

import os
import sys
import json
import tempfile

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_utils import save_petri_dish


def main():
    print("=" * 60)
    print("Testing Simplified Year System")
    print("=" * 60)
    
    # Test 1: No absolute_year
    print("\n1. Checking absolute_year is removed...")
    petri = PetriDish(n=100, rate=0.005, growth_phase=5)
    
    assert not hasattr(petri, 'absolute_year'), "absolute_year should not exist"
    assert petri.year == 0, "year should start at 0"
    print("   ✅ No absolute_year field")
    print(f"   ✅ year = {petri.year}")
    
    # Test 2: Year increments
    print("\n2. Testing year increments...")
    petri.increment_year()
    assert petri.year == 1
    petri.increment_year()
    assert petri.year == 2
    print(f"   ✅ Year increments: 0 → 1 → 2")
    
    # Test 3: Metadata tracking
    print("\n3. Testing metadata...")
    cell = Cell(n=100, rate=0.005)
    petri2 = PetriDish.from_snapshot_cell(
        cell=cell,
        growth_phase=7,
        metadata={
            'initial_year': 50,
            'individual_type': 'test'
        }
    )
    
    assert petri2.year == 0, "New PetriDish starts at year 0"
    assert petri2.metadata['initial_year'] == 50
    print(f"   ✅ PetriDish age: {petri2.year}")
    print(f"   ✅ Snapshot origin: {petri2.metadata['initial_year']}")
    
    # Simulate growth
    for _ in range(20):
        petri2.year += 1
    
    print(f"   ✅ After 20 years: age = {petri2.year}")
    print(f"   ✅ Simulation year = {petri2.metadata['initial_year'] + petri2.year}")
    
    # Test 4: Save without absolute_year
    print("\n4. Testing save format...")
    with tempfile.NamedTemporaryFile(suffix='.json', delete=False) as f:
        temp_file = f.name
    
    try:
        save_petri_dish(petri2, temp_file, compress=False)
        
        with open(temp_file, 'r') as f:
            data = json.load(f)
        
        assert 'absolute_year' not in data['metadata']
        assert data['metadata']['year'] == 20
        assert data['metadata']['initial_year'] == 50
        
        print("   ✅ No absolute_year in saved data")
        print(f"   ✅ Saved year: {data['metadata']['year']}")
        print(f"   ✅ Saved initial_year: {data['metadata']['initial_year']}")
        
    finally:
        if os.path.exists(temp_file):
            os.remove(temp_file)
    
    print("\n" + "=" * 60)
    print("🎉 SUCCESS: Year system simplified!")
    print("   - No more confusing absolute_year")
    print("   - Simple year counter (0, 1, 2, ...)")
    print("   - Clear metadata (initial_year for origin)")
    print("=" * 60)
    
    return 0


if __name__ == "__main__":
    sys.exit(main())