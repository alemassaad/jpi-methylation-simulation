#!/usr/bin/env python3
"""
Test the JSON consolidation functionality.
"""

import os
import sys
import json
import tempfile
import numpy as np

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from cell import Cell, PetriDish
from pipeline_analysis import analyze_populations_from_dishes


def test_json_consolidation():
    """Test that cell_jsd_analysis.json is created with correct structure."""
    print("\nTest: JSON Consolidation")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create mock PetriDish objects for testing
        print("Creating mock PetriDish objects...")
        
        # Create mutant dishes
        mutant_dishes = []
        for i in range(3):
            petri = PetriDish()
            # Add cells with varying JSD values
            for j in range(10):
                cell = Cell(n=100, rate=0.005)
                cell.cell_jsd = 0.2 + i * 0.01 + j * 0.001  # Varying values
                petri.cells.append(cell)
            mutant_dishes.append(petri)
        
        # Create control1 dishes
        control1_dishes = []
        for i in range(3):
            petri = PetriDish()
            for j in range(10):
                cell = Cell(n=100, rate=0.005)
                cell.cell_jsd = 0.15 + i * 0.01 + j * 0.001
                petri.cells.append(cell)
            control1_dishes.append(petri)
        
        # Create control2 dishes
        control2_dishes = []
        for i in range(3):
            petri = PetriDish()
            for j in range(10):
                cell = Cell(n=100, rate=0.005)
                cell.cell_jsd = 0.1 + i * 0.01 + j * 0.001
                petri.cells.append(cell)
            control2_dishes.append(petri)
        
        # Run analysis
        print("Running analysis...")
        results_dir = os.path.join(tmpdir, "results")
        results = analyze_populations_from_dishes(
            mutant_dishes, control1_dishes, control2_dishes, results_dir
        )
        
        # Check that new consolidated file exists
        cell_jsd_path = os.path.join(results_dir, "cell_jsd_analysis.json")
        assert os.path.exists(cell_jsd_path), "cell_jsd_analysis.json not created"
        print("  ✓ cell_jsd_analysis.json created")
        
        # Load and verify structure
        with open(cell_jsd_path, 'r') as f:
            cell_jsd_data = json.load(f)
        
        # Check top-level keys
        assert "summary_statistics" in cell_jsd_data, "Missing summary_statistics"
        assert "statistical_tests" in cell_jsd_data, "Missing statistical_tests"
        assert "individual_means" in cell_jsd_data, "Missing individual_means"
        print("  ✓ All top-level keys present")
        
        # Check summary_statistics structure
        for batch in ["mutant", "control1", "control2"]:
            assert batch in cell_jsd_data["summary_statistics"], f"Missing {batch} in summary_statistics"
            batch_stats = cell_jsd_data["summary_statistics"][batch]
            
            required_fields = ["mean", "std", "median", "min", "max", "n_individuals"]
            for field in required_fields:
                assert field in batch_stats, f"Missing {field} in {batch} statistics"
            
            # Verify n_individuals matches our test data
            assert batch_stats["n_individuals"] == 3, f"Wrong n_individuals for {batch}"
        print("  ✓ summary_statistics structure correct")
        
        # Check statistical_tests structure
        test_pairs = ["mutant_vs_control1", "mutant_vs_control2", "control1_vs_control2"]
        for pair in test_pairs:
            assert pair in cell_jsd_data["statistical_tests"], f"Missing {pair} in statistical_tests"
            test_data = cell_jsd_data["statistical_tests"][pair]
            assert "t_statistic" in test_data, f"Missing t_statistic in {pair}"
            assert "p_value" in test_data, f"Missing p_value in {pair}"
        print("  ✓ statistical_tests structure correct")
        
        # Check individual_means structure
        for batch in ["mutant", "control1", "control2"]:
            assert batch in cell_jsd_data["individual_means"], f"Missing {batch} in individual_means"
            means = cell_jsd_data["individual_means"][batch]
            assert isinstance(means, list), f"{batch} means not a list"
            assert len(means) == 3, f"Wrong number of individuals in {batch}"
            
            # Verify values are reasonable
            for mean_val in means:
                assert 0 <= mean_val <= 1, f"Invalid JSD value: {mean_val}"
        print("  ✓ individual_means structure correct")
        
        # Check that old files are NOT created (no backward compatibility)
        stats_path = os.path.join(results_dir, "statistics.json")
        distributions_path = os.path.join(results_dir, "jsd_distributions.json")
        
        assert not os.path.exists(stats_path), "statistics.json should not be created"
        assert not os.path.exists(distributions_path), "jsd_distributions.json should not be created"
        print("  ✓ Old files not created (no backward compatibility)")
        
        print("\n✅ JSON consolidation test passed!")
        return True


def test_semantic_improvements():
    """Test that semantic improvements are working."""
    print("\nTest: Semantic Improvements")
    print("=" * 50)
    
    with tempfile.TemporaryDirectory() as tmpdir:
        # Create simple test data
        mutant_dishes = [create_test_petri_dish(0.2) for _ in range(2)]
        control1_dishes = [create_test_petri_dish(0.15) for _ in range(2)]
        control2_dishes = [create_test_petri_dish(0.1) for _ in range(2)]
        
        results_dir = os.path.join(tmpdir, "results")
        results = analyze_populations_from_dishes(
            mutant_dishes, control1_dishes, control2_dishes, results_dir
        )
        
        # Load consolidated file
        with open(results["cell_jsd_path"], 'r') as f:
            data = json.load(f)
        
        # Check semantic improvements
        # 1. "n_individuals" instead of "n"
        assert "n_individuals" in data["summary_statistics"]["mutant"]
        assert "n" not in data["summary_statistics"]["mutant"]
        print("  ✓ Using 'n_individuals' instead of 'n'")
        
        # 2. "individual_means" instead of raw array
        assert "individual_means" in data
        assert isinstance(data["individual_means"], dict)
        print("  ✓ Using 'individual_means' for clarity")
        
        # 3. Clear separation of statistics vs tests
        assert "summary_statistics" in data
        assert "statistical_tests" in data
        assert "comparisons" not in data  # Old name not used
        print("  ✓ Clear separation of statistics and tests")
        
        print("\n✅ Semantic improvements test passed!")
        return True


def create_test_petri_dish(base_jsd):
    """Helper to create a test PetriDish with cells."""
    petri = PetriDish()
    for i in range(5):
        cell = Cell(n=100, rate=0.005)
        cell.cell_jsd = base_jsd + i * 0.01
        petri.cells.append(cell)
    return petri


def main():
    """Run all tests."""
    print("=" * 60)
    print("Testing JSON Consolidation")
    print("=" * 60)
    
    tests = [
        test_json_consolidation,
        test_semantic_improvements
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
        print("\n🎉 All JSON consolidation tests passed!")
        print("   ✅ cell_jsd_analysis.json created successfully")
        print("   ✅ Structure matches specification")
        print("   ✅ No backward compatibility files created")
        print("   ✅ Semantic improvements working")
        return 0
    else:
        print(f"\n⚠️  {failed} test(s) failed")
        return 1


if __name__ == "__main__":
    sys.exit(main())