#!/usr/bin/env python3
"""
Comprehensive test suite for Step1-Prime simulation.
Tests all aspects of the simulation with small, fast-running parameters.
"""

import json
import gzip
import os
import sys
import copy
import tempfile
import shutil
from typing import List, Dict, Any

# Add parent directory to path to import cell module
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cell import PetriDish, Cell, JS_div, BASELINE_METHYLATION_DISTRIBUTION


class TestResults:
    """Track test results and provide summary."""
    def __init__(self):
        self.passed = []
        self.failed = []
        self.total = 0
    
    def record(self, test_name: str, passed: bool, message: str = ""):
        self.total += 1
        if passed:
            self.passed.append(test_name)
            print(f"✓ {test_name}")
        else:
            self.failed.append((test_name, message))
            print(f"✗ {test_name}: {message}")
    
    def summary(self):
        print("\n" + "="*60)
        print("TEST SUMMARY")
        print("="*60)
        print(f"Total tests: {self.total}")
        print(f"Passed: {len(self.passed)} ({len(self.passed)/self.total*100:.1f}%)")
        print(f"Failed: {len(self.failed)} ({len(self.failed)/self.total*100:.1f}%)")
        
        if self.failed:
            print("\nFailed tests:")
            for test_name, message in self.failed:
                print(f"  - {test_name}: {message}")
        
        return len(self.failed) == 0


# Initialize test results tracker
results = TestResults()


# ==============================================================================
# SECTION 1: GROWTH PHASE TESTS
# ==============================================================================

def test_population_doubling():
    """Test 1.1: Verify population doubles correctly during growth phase."""
    test_name = "Population Doubling"
    try:
        petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, target_population=16)
        
        # Expected populations for years 0-4
        expected = {0: 1, 1: 2, 2: 4, 3: 8, 4: 16}
        
        for year in range(5):
            if year > 0:
                petri.simulate_year()
            
            actual = len(petri.cells)
            if actual != expected[year]:
                results.record(test_name, False, 
                             f"Year {year}: expected {expected[year]}, got {actual}")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_division_creates_identical_copies():
    """Test 1.2: Verify cell division creates identical copies."""
    test_name = "Division Creates Identical Copies"
    try:
        # Create a cell with some methylation
        cell = Cell(n=100, rate=0.01, gene_size=5)
        # Manually methylate some sites
        for i in range(10):
            cell.cpg_sites[i] = 1
        cell.compute_methylation_distribution()
        
        # Create daughter cell
        daughter = cell.create_daughter_cell()
        
        # Check they're identical but different objects
        if cell.cpg_sites != daughter.cpg_sites:
            results.record(test_name, False, "Daughter has different cpg_sites")
            return
        
        if id(cell.cpg_sites) == id(daughter.cpg_sites):
            results.record(test_name, False, "Daughter shares same cpg_sites object (not a copy)")
            return
        
        # Modify daughter and verify parent unchanged
        daughter.cpg_sites[50] = 1
        if cell.cpg_sites[50] == 1:
            results.record(test_name, False, "Modifying daughter affected parent")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_methylation_during_growth():
    """Test 1.3: Verify methylation accumulates during growth phase."""
    test_name = "Methylation During Growth"
    try:
        petri = PetriDish(rate=0.1, n=100, gene_size=5, seed=42, target_population=8)
        
        # Run for 3 years (1->2->4->8 cells)
        jsd_values = []
        for year in range(4):
            if year > 0:
                petri.simulate_year()
            # Calculate mean JSD
            mean_jsd = sum(c.JSD for c in petri.cells) / len(petri.cells)
            jsd_values.append(mean_jsd)
        
        # Check JSD increases monotonically
        for i in range(1, len(jsd_values)):
            if jsd_values[i] <= jsd_values[i-1]:
                results.record(test_name, False, 
                             f"JSD not increasing: year {i-1}={jsd_values[i-1]:.4f}, "
                             f"year {i}={jsd_values[i]:.4f}")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


# ==============================================================================
# SECTION 2: STEADY STATE TESTS  
# ==============================================================================

def test_population_maintenance():
    """Test 2.1: Verify population stays around target in steady state."""
    test_name = "Population Maintenance"
    try:
        petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, target_population=16)
        
        # Run to year 15 (steady state starts at year 4)
        for _ in range(15):
            petri.simulate_year()
        
        # Check years 5-15 are in steady state
        for year in range(5, 16):
            pop = len(petri.history[str(year)])
            # Allow more variation due to randomness (6-40 for target of 16)
            if pop < 6 or pop > 40:
                results.record(test_name, False, 
                             f"Year {year}: population {pop} outside range [6, 40]")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_culling_statistics():
    """Test 2.2: Verify culling gives approximately 50% survival."""
    test_name = "Culling Statistics"
    try:
        petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, target_population=16)
        
        # Get to steady state
        for _ in range(4):
            petri.simulate_year()
        
        # Track survival rates over multiple cullings
        survival_rates = []
        for _ in range(10):
            initial = len(petri.cells)
            petri.divide_cells()
            doubled = len(petri.cells)
            petri.random_cull_cells()
            final = len(petri.cells)
            
            survival_rate = final / doubled
            survival_rates.append(survival_rate)
            
            # Reset for next test
            petri.cells = petri.cells[:initial]
        
        # Check mean survival rate is close to 0.5
        mean_survival = sum(survival_rates) / len(survival_rates)
        if abs(mean_survival - 0.5) > 0.1:  # Allow 10% deviation
            results.record(test_name, False, 
                         f"Mean survival rate {mean_survival:.2f} not close to 0.5")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_phase_transition():
    """Test 2.3: Verify transition from growth to steady state."""
    test_name = "Phase Transition"
    try:
        petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, target_population=16)
        
        # Year 3: should be growth (starts with 4, ends with 8 cells)
        for _ in range(3):
            petri.simulate_year()
        
        # After year 3, we have 8 cells, still in growth
        if petri.reached_target:
            results.record(test_name, False, "reached_target True before reaching target")
            return
        
        # Year 4: starts with 8, divides to 16, should trigger reached_target
        petri.simulate_year()
        
        # Now we should have reached target (16 cells)
        if not petri.reached_target:
            # Check current population to debug
            current_pop = len(petri.cells)
            results.record(test_name, False, 
                         f"reached_target False after year 4 (pop={current_pop}, target=16)")
            return
        
        # Manually reduce population to test persistence
        petri.cells = petri.cells[:10]
        
        # Year 5: should stay in steady state despite low population
        petri.simulate_year()
        
        if not petri.reached_target:
            results.record(test_name, False, "reached_target reset after population drop")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


# ==============================================================================
# SECTION 3: METHYLATION MECHANICS TESTS
# ==============================================================================

def test_initial_state():
    """Test 3.1: Verify cells start unmethylated."""
    test_name = "Initial State"
    try:
        cell = Cell(n=100, rate=0.01, gene_size=5)
        
        # Check all sites unmethylated
        if any(site != 0 for site in cell.cpg_sites):
            results.record(test_name, False, "Initial cpg_sites not all zeros")
            return
        
        # Check methylation proportion is 0
        if cell.methylation_proportion != 0:
            results.record(test_name, False, f"Initial methylation_proportion = {cell.methylation_proportion}")
            return
        
        # Check initial JSD
        if cell.JSD != 0:
            results.record(test_name, False, f"Initial JSD = {cell.JSD}")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_methylation_accumulation():
    """Test 3.2: Verify methylation only increases (no demethylation)."""
    test_name = "Methylation Accumulation"
    try:
        cell = Cell(n=100, rate=0.5, gene_size=5)  # High rate for testing
        
        previous_methylated = 0
        for _ in range(5):
            cell.methylate()
            
            # Count methylated sites
            current_methylated = sum(cell.cpg_sites)
            
            # Check never decreases
            if current_methylated < previous_methylated:
                results.record(test_name, False, "Methylation decreased")
                return
            
            # Check all values are 0 or 1
            if any(site not in [0, 1] for site in cell.cpg_sites):
                results.record(test_name, False, "Invalid cpg_site value")
                return
            
            previous_methylated = current_methylated
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_methylation_distribution():
    """Test 3.3: Verify methylation distribution calculated correctly."""
    test_name = "Methylation Distribution"
    try:
        cell = Cell(n=100, rate=0.01, gene_size=5)
        
        # Manually set some methylation pattern
        # First gene (sites 0-4): 2 methylated
        cell.cpg_sites[0] = 1
        cell.cpg_sites[1] = 1
        # Second gene (sites 5-9): 0 methylated
        # Third gene (sites 10-14): 5 methylated
        for i in range(10, 15):
            cell.cpg_sites[i] = 1
        
        cell.compute_methylation_distribution()
        
        # Should have distribution[0]=18, distribution[2]=1, distribution[5]=1
        # Normalized: [0.9, 0.0, 0.05, 0.0, 0.0, 0.05]
        expected = [18/20, 0, 1/20, 0, 0, 1/20]
        
        for i, (actual, exp) in enumerate(zip(cell.methylation_distribution, expected)):
            if abs(actual - exp) > 1e-6:
                results.record(test_name, False, 
                             f"Distribution[{i}]: expected {exp}, got {actual}")
                return
        
        # Check sum equals 1
        total = sum(cell.methylation_distribution)
        if abs(total - 1.0) > 1e-6:
            results.record(test_name, False, f"Distribution sum = {total}, not 1.0")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


# ==============================================================================
# SECTION 4: OUTPUT FORMAT TESTS
# ==============================================================================

def test_json_structure():
    """Test 4.1: Verify JSON output structure is correct."""
    test_name = "JSON Structure"
    try:
        petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, target_population=8)
        
        # Run for 5 years
        for _ in range(5):
            petri.simulate_year()
        
        # Check history structure
        for year in range(6):  # 0-5
            year_str = str(year)
            if year_str not in petri.history:
                results.record(test_name, False, f"Year {year} missing from history")
                return
            
            if not isinstance(petri.history[year_str], list):
                results.record(test_name, False, f"Year {year} not a list")
                return
            
            # Check expected cell count
            if year == 0:
                expected = 1
            elif year <= 3:
                expected = 2 ** year
            else:
                expected = 8  # target reached
            
            if year < 4 and len(petri.history[year_str]) != expected:
                results.record(test_name, False, 
                             f"Year {year}: expected {expected} cells, got {len(petri.history[year_str])}")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_cell_format():
    """Test 4.2: Verify cell dictionary format is correct."""
    test_name = "Cell Dictionary Format"
    try:
        cell = Cell(n=100, rate=0.01, gene_size=5)
        cell.methylate()
        cell_dict = cell.to_dict()
        
        # Check required keys
        required_keys = ['cpg_sites', 'methylation_proportion', 'methylation_distribution',
                        'jsd', 'age', 'gene_size', 'rate']
        
        for key in required_keys:
            if key not in cell_dict:
                results.record(test_name, False, f"Missing key: {key}")
                return
        
        # Check data types
        if not isinstance(cell_dict['cpg_sites'], list):
            results.record(test_name, False, "cpg_sites not a list")
            return
        
        if not isinstance(cell_dict['methylation_proportion'], (int, float)):
            results.record(test_name, False, "methylation_proportion not numeric")
            return
        
        if not isinstance(cell_dict['methylation_distribution'], list):
            results.record(test_name, False, "methylation_distribution not a list")
            return
        
        if len(cell_dict['methylation_distribution']) != 6:  # gene_size + 1
            results.record(test_name, False, 
                         f"methylation_distribution wrong length: {len(cell_dict['methylation_distribution'])}")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_file_saving_loading():
    """Test 4.3: Verify file saving and loading works correctly."""
    test_name = "File Saving/Loading"
    try:
        # Create temporary directory
        with tempfile.TemporaryDirectory() as tmpdir:
            petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, target_population=4)
            
            # Run for 3 years
            for _ in range(3):
                petri.simulate_year()
            
            # Save to file
            filepath = petri.save_history(directory=tmpdir)
            
            # Check file exists
            if not os.path.exists(filepath):
                results.record(test_name, False, "File not created")
                return
            
            # Check it's compressed
            if not filepath.endswith('.json.gz'):
                results.record(test_name, False, "File not gzipped")
                return
            
            # Try to load it
            with gzip.open(filepath, 'rt') as f:
                loaded = json.load(f)
            
            # Verify content matches
            if str(3) not in loaded:
                results.record(test_name, False, "Year 3 not in loaded data")
                return
            
            # Year 3 should have 4 cells for target_population=4
            if len(loaded['3']) != 4:
                results.record(test_name, False, f"Year 3 has {len(loaded['3'])} cells, expected 4")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


# ==============================================================================
# SECTION 5: REPRODUCIBILITY TESTS
# ==============================================================================

def test_same_seed_same_results():
    """Test 5.1: Verify same seed produces identical results."""
    test_name = "Same Seed Same Results"
    try:
        # Important: Must use separate random instances or reset seed
        import random
        
        # Run simulation twice with same seed - reset seed each time
        random.seed(12345)
        petri1 = PetriDish(rate=0.1, n=100, gene_size=5, seed=12345, target_population=8)
        for _ in range(5):
            petri1.simulate_year()
            
        # Reset and run again
        random.seed(12345)
        petri2 = PetriDish(rate=0.1, n=100, gene_size=5, seed=12345, target_population=8)
        for _ in range(5):
            petri2.simulate_year()
        
        # Compare year 5 cells
        cells1 = petri1.history['5']
        cells2 = petri2.history['5']
        
        if len(cells1) != len(cells2):
            results.record(test_name, False, 
                         f"Different cell counts: {len(cells1)} vs {len(cells2)}")
            return
        
        # Check methylation patterns are identical
        for i, (c1, c2) in enumerate(zip(cells1, cells2)):
            if c1['cpg_sites'] != c2['cpg_sites']:
                results.record(test_name, False, f"Cell {i} has different methylation")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_different_seeds_different_results():
    """Test 5.2: Verify different seeds produce different results."""
    test_name = "Different Seeds Different Results"
    try:
        # Run simulation twice with different seeds
        petri1 = PetriDish(rate=0.1, n=100, gene_size=5, seed=111, target_population=16)
        petri2 = PetriDish(rate=0.1, n=100, gene_size=5, seed=999, target_population=16)
        
        # Run to steady state and beyond
        for _ in range(6):
            petri1.simulate_year()
            petri2.simulate_year()
        
        # Year 5 should be in steady state with culling
        # Cell counts should differ due to random culling
        cells1 = petri1.history['5']
        cells2 = petri2.history['5']
        
        # Check if populations are different (very likely with different seeds)
        if len(cells1) == len(cells2):
            # Even if same count, methylation should differ
            same_methylation = True
            for c1, c2 in zip(cells1[:min(len(cells1), len(cells2))], 
                            cells2[:min(len(cells1), len(cells2))]):
                if c1['cpg_sites'] != c2['cpg_sites']:
                    same_methylation = False
                    break
            
            if same_methylation:
                results.record(test_name, False, "Different seeds produced identical results")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


# ==============================================================================
# SECTION 6: EDGE CASES
# ==============================================================================

def test_very_small_population():
    """Test 6.1: Verify simulation works with very small target population."""
    test_name = "Very Small Population"
    try:
        petri = PetriDish(rate=0.01, n=50, gene_size=5, seed=42, target_population=2)
        
        # Run for 5 years
        for _ in range(5):
            petri.simulate_year()
        
        # Check year 1 has 2 cells (reached target)
        if len(petri.history['1']) != 2:
            results.record(test_name, False, f"Year 1 has {len(petri.history['1'])} cells, expected 2")
            return
        
        # Check steady state maintains around 2
        for year in range(2, 6):
            pop = len(petri.history[str(year)])
            if pop < 1 or pop > 4:  # Allow some variation
                results.record(test_name, False, f"Year {year} population {pop} out of range")
                return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_high_methylation_rate():
    """Test 6.2: Verify simulation handles high methylation rate."""
    test_name = "High Methylation Rate"
    try:
        petri = PetriDish(rate=0.9, n=100, gene_size=5, seed=42, target_population=4)
        
        # Run for 3 years
        for _ in range(3):
            petri.simulate_year()
        
        # Check cells are heavily methylated
        year3_cells = petri.history['3']
        methylation_props = [c['methylation_proportion'] for c in year3_cells]
        mean_methylation = sum(methylation_props) / len(methylation_props)
        
        if mean_methylation < 0.8:  # Should be very high with 90% rate
            results.record(test_name, False, 
                         f"Mean methylation {mean_methylation:.2%} too low for 90% rate")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_gene_size_validation():
    """Test 6.3: Verify gene size must divide n evenly."""
    test_name = "Gene Size Validation"
    try:
        # This should raise an error
        try:
            cell = Cell(n=100, rate=0.01, gene_size=7)  # 100 not divisible by 7
            results.record(test_name, False, "No error raised for invalid gene_size")
        except ValueError as e:
            if "divide" in str(e).lower():
                results.record(test_name, True)
            else:
                results.record(test_name, False, f"Wrong error message: {e}")
    except Exception as e:
        results.record(test_name, False, str(e))


# ==============================================================================
# SECTION 7: COMPATIBILITY TESTS
# ==============================================================================

def test_step23_compatible_format():
    """Test 7.1: Verify output is compatible with step23 pipeline."""
    test_name = "Step23 Compatible Format"
    try:
        with tempfile.TemporaryDirectory() as tmpdir:
            petri = PetriDish(rate=0.01, n=100, gene_size=5, seed=42, target_population=8)
            
            # Run simulation
            for _ in range(10):
                petri.simulate_year()
            
            # Save file
            filepath = petri.save_history(directory=tmpdir)
            
            # Load like step23 would
            with gzip.open(filepath, 'rt') as f:
                history = json.load(f)
            
            # Try to extract a year (like step23 does)
            year_5 = history.get('5')
            if not year_5:
                results.record(test_name, False, "Cannot extract year 5")
                return
            
            # Check cells have required fields for step23
            for cell in year_5:
                if 'jsd' not in cell:
                    results.record(test_name, False, "Cell missing 'jsd' field")
                    return
                if 'cpg_sites' not in cell:
                    results.record(test_name, False, "Cell missing 'cpg_sites' field")
                    return
                if 'rate' not in cell:
                    results.record(test_name, False, "Cell missing 'rate' field")
                    return
            
            results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


def test_statistics_extraction():
    """Test 7.2: Verify statistics can be calculated from output."""
    test_name = "Statistics Extraction"
    try:
        petri = PetriDish(rate=0.05, n=100, gene_size=5, seed=42, target_population=8)
        
        # Run simulation
        for _ in range(5):
            petri.simulate_year()
        
        # Extract statistics for year 5
        year5_cells = petri.history['5']
        
        # Calculate JSD statistics
        jsd_values = [c['jsd'] for c in year5_cells]
        mean_jsd = sum(jsd_values) / len(jsd_values)
        
        # Calculate methylation statistics
        meth_values = [c['methylation_proportion'] for c in year5_cells]
        mean_meth = sum(meth_values) / len(meth_values)
        
        # Basic sanity checks
        if mean_jsd <= 0:
            results.record(test_name, False, f"Mean JSD {mean_jsd} should be positive")
            return
        
        if mean_meth <= 0 or mean_meth >= 1:
            results.record(test_name, False, f"Mean methylation {mean_meth} out of range")
            return
        
        # Check can calculate median
        import statistics
        median_jsd = statistics.median(jsd_values)
        if median_jsd <= 0:
            results.record(test_name, False, "Median JSD calculation failed")
            return
        
        results.record(test_name, True)
    except Exception as e:
        results.record(test_name, False, str(e))


# ==============================================================================
# MAIN TEST RUNNER
# ==============================================================================

def run_all_tests():
    """Run all tests and report results."""
    print("="*60)
    print("STEP1-PRIME COMPREHENSIVE TEST SUITE")
    print("="*60)
    
    print("\n--- GROWTH PHASE TESTS ---")
    test_population_doubling()
    test_division_creates_identical_copies()
    test_methylation_during_growth()
    
    print("\n--- STEADY STATE TESTS ---")
    test_population_maintenance()
    test_culling_statistics()
    test_phase_transition()
    
    print("\n--- METHYLATION MECHANICS TESTS ---")
    test_initial_state()
    test_methylation_accumulation()
    test_methylation_distribution()
    
    print("\n--- OUTPUT FORMAT TESTS ---")
    test_json_structure()
    test_cell_format()
    test_file_saving_loading()
    
    print("\n--- REPRODUCIBILITY TESTS ---")
    test_same_seed_same_results()
    test_different_seeds_different_results()
    
    print("\n--- EDGE CASE TESTS ---")
    test_very_small_population()
    test_high_methylation_rate()
    test_gene_size_validation()
    
    print("\n--- COMPATIBILITY TESTS ---")
    test_step23_compatible_format()
    test_statistics_extraction()
    
    # Print summary
    all_passed = results.summary()
    
    if all_passed:
        print("\n🎉 ALL TESTS PASSED! 🎉")
        return 0
    else:
        print("\n❌ SOME TESTS FAILED")
        return 1


if __name__ == "__main__":
    exit(run_all_tests())