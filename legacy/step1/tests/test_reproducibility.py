#!/usr/bin/env python3
"""
Reproducibility test to verify simulation results are deterministic with fixed random seed.
This test can be used to verify that code changes don't alter the simulation behavior.

Usage:
    python test_reproducibility.py              # Run test and save results
    python test_reproducibility.py --check      # Check against saved results
"""

import random
import json
import sys
import os

# Add parent directory to path to import cell module
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from cell import Cell, History

def test_single_cell_aging():
    """Test that a single cell ages identically with fixed seed."""
    print("Test 1: Single cell aging with fixed seed")
    
    # Fix random seed for reproducibility
    random.seed(42)
    
    # Create and age a cell
    cell = Cell(n=100, rate=0.05, gene_size=5)
    
    # Store states at each year
    states = []
    for year in range(10):
        cell.age_1_year()
        state = {
            'year': year + 1,
            'methylation_proportion': cell.methylation_proportion,
            'distribution': cell.methylation_distribution[:],
            'jsd': cell.JSD,
            'cpg_sites_sum': sum(cell.cpg_sites)  # Verify actual methylation
        }
        states.append(state)
    
    # Print results
    print(f"  Final methylation: {cell.methylation_proportion:.4f}")
    print(f"  Final JSD: {cell.JSD:.4f}")
    print(f"  Final distribution: {[f'{x:.3f}' for x in cell.methylation_distribution]}")
    
    return states

def test_full_simulation():
    """Test a small but complete simulation."""
    print("\nTest 2: Full simulation with 100 cells for 20 years")
    
    # Fix random seed
    random.seed(12345)
    
    # Run simulation
    m = 100  # 100 cells
    n = 100  # 100 sites per cell
    t_max = 20  # 20 years
    rate = 0.01
    
    history = History(m=m, n=n, rate=rate, gene_size=5, t_max=t_max)
    history.age_multiple_years(t_max)
    
    # Extract final year statistics
    final_year = history.history[t_max]
    methylation_props = [cell['methylation_proportion'] for cell in final_year]
    jsd_scores = [cell['jsd'] for cell in final_year]
    
    # Calculate aggregate statistics
    avg_methylation = sum(methylation_props) / len(methylation_props)
    avg_jsd = sum(jsd_scores) / len(jsd_scores)
    min_methylation = min(methylation_props)
    max_methylation = max(methylation_props)
    
    print(f"  Average methylation: {avg_methylation:.6f}")
    print(f"  Min methylation: {min_methylation:.6f}")
    print(f"  Max methylation: {max_methylation:.6f}")
    print(f"  Average JSD: {avg_jsd:.6f}")
    
    # Check specific cells
    cells_to_check = [0, 50, 99]
    cell_states = {}
    for idx in cells_to_check:
        cell = final_year[idx]
        cell_states[f'cell_{idx}'] = {
            'methylation': cell['methylation_proportion'],
            'jsd': cell['jsd'],
            'distribution': cell['methylation_distribution'],
            'cpg_sum': sum(cell['cpg_sites'])
        }
    
    return {
        'avg_methylation': avg_methylation,
        'avg_jsd': avg_jsd,
        'min_methylation': min_methylation,
        'max_methylation': max_methylation,
        'cell_states': cell_states
    }

def test_edge_cases():
    """Test edge cases like fully methylated cells."""
    print("\nTest 3: Edge case - high methylation rate")
    
    random.seed(999)
    
    # High rate to ensure some cells become fully methylated
    cell = Cell(n=50, rate=0.2, gene_size=5)
    
    methylation_history = []
    fully_methylated_year = None
    
    for year in range(30):
        cell.age_1_year()
        methylation_history.append(cell.methylation_proportion)
        
        if cell.methylation_proportion >= 1.0 and fully_methylated_year is None:
            fully_methylated_year = year + 1
            print(f"  Cell fully methylated at year {fully_methylated_year}")
            # Continue aging to verify early exit works
            for extra_year in range(5):
                cell.age_1_year()
                methylation_history.append(cell.methylation_proportion)
            break
    
    return {
        'methylation_history': methylation_history,
        'fully_methylated_year': fully_methylated_year
    }

def test_different_parameters():
    """Test with various parameter combinations."""
    print("\nTest 4: Different parameter combinations")
    
    test_configs = [
        {'n': 50, 'rate': 0.001, 'gene_size': 5, 'years': 10, 'seed': 111},
        {'n': 200, 'rate': 0.05, 'gene_size': 5, 'years': 5, 'seed': 222},
        {'n': 75, 'rate': 0.1, 'gene_size': 5, 'years': 15, 'seed': 333},
    ]
    
    results = []
    for i, config in enumerate(test_configs):
        random.seed(config['seed'])
        cell = Cell(n=config['n'], rate=config['rate'], gene_size=config['gene_size'])
        
        for _ in range(config['years']):
            cell.age_1_year()
        
        result = {
            'config': config,
            'final_methylation': cell.methylation_proportion,
            'final_jsd': cell.JSD,
            'final_distribution': cell.methylation_distribution[:]
        }
        results.append(result)
        
        print(f"  Config {i+1}: n={config['n']}, rate={config['rate']}, "
              f"final_methylation={cell.methylation_proportion:.4f}")
    
    return results

def run_all_tests():
    """Run all tests and return results."""
    results = {
        'test_single_cell': test_single_cell_aging(),
        'test_full_simulation': test_full_simulation(),
        'test_edge_cases': test_edge_cases(),
        'test_different_parameters': test_different_parameters()
    }
    return results

def save_results(results):
    """Save test results to file."""
    # Save in the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(script_dir, 'test_reproducibility_expected.json')
    with open(filename, 'w') as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to {os.path.basename(filename)}")

def check_against_expected():
    """Compare current results against saved expected results."""
    # Load from the same directory as this script
    script_dir = os.path.dirname(os.path.abspath(__file__))
    filename = os.path.join(script_dir, 'test_reproducibility_expected.json')
    
    if not os.path.exists(filename):
        print(f"Error: {os.path.basename(filename)} not found. Run without --check first to generate expected results.")
        return False
    
    print("Running tests and comparing against expected results...\n")
    
    # Load expected results
    with open(filename, 'r') as f:
        expected = json.load(f)
    
    # Run tests
    actual = run_all_tests()
    
    # Compare results
    print("\n" + "="*60)
    print("COMPARISON RESULTS:")
    print("="*60)
    
    def compare_values(expected, actual, tolerance=1e-10):
        """Compare two values with tolerance for floats."""
        if type(expected) != type(actual):
            # Handle int/float type differences
            if isinstance(expected, (int, float)) and isinstance(actual, (int, float)):
                return abs(float(expected) - float(actual)) < tolerance
            return False
        
        if isinstance(expected, float):
            return abs(expected - actual) < tolerance
        elif isinstance(expected, list):
            if len(expected) != len(actual):
                return False
            return all(compare_values(e, a, tolerance) for e, a in zip(expected, actual))
        elif isinstance(expected, dict):
            if set(expected.keys()) != set(actual.keys()):
                return False
            return all(compare_values(expected[k], actual[k], tolerance) for k in expected)
        else:
            return expected == actual
    
    all_match = True
    for test_name in expected:
        if test_name in actual:
            if compare_values(expected[test_name], actual[test_name]):
                print(f"✓ {test_name}: PASS")
            else:
                print(f"✗ {test_name}: FAIL")
                all_match = False
        else:
            print(f"✗ {test_name}: MISSING")
            all_match = False
    
    print("="*60)
    if all_match:
        print("✓ ALL TESTS PASS - Results match expected values!")
    else:
        print("✗ SOME TESTS FAILED - Results differ from expected values!")
    print("="*60)
    
    return all_match

def main():
    """Main function."""
    if len(sys.argv) > 1 and sys.argv[1] == '--check':
        success = check_against_expected()
        sys.exit(0 if success else 1)
    else:
        print("Running reproducibility tests...\n")
        results = run_all_tests()
        save_results(results)
        print("\nTo verify results in the future, run:")
        print("  python test_reproducibility.py --check")

if __name__ == "__main__":
    main()