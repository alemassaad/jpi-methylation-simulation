#!/usr/bin/env python3
"""
Comprehensive test script for gene JSD extraction from different data formats.
Tests the extract_gene_jsd_from_history function with all possible scenarios.
"""

import sys
import os
import json
import tempfile
import traceback

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1'))

from core.pipeline_utils import extract_gene_jsd_from_history

def test_correct_format():
    """Test extraction from correct format (nested in history)"""
    print("\n" + "="*60)
    print("TEST 1: Correct Format (nested in history)")
    print("="*60)
    
    # Create test data in the correct format
    sim_data = {
        "parameters": {"n": 100, "gene_size": 5},
        "history": {
            "0": {
                "cells": [],
                "gene_jsd": [0.0] * 20,
                "mean_gene_jsd": 0.0,
                "median_gene_jsd": 0.0
            },
            "10": {
                "cells": [],
                "gene_jsd": [0.1, 0.2, 0.15] + [0.1] * 17,
                "mean_gene_jsd": 0.11,
                "median_gene_jsd": 0.1
            },
            "20": {
                "cells": [],
                "gene_jsd": [0.2, 0.3, 0.25] + [0.2] * 17,
                "mean_gene_jsd": 0.21,
                "median_gene_jsd": 0.2
            }
        }
    }
    
    # Extract data
    gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(sim_data)
    
    # Verify extraction
    assert len(gene_jsd_hist) == 3, f"Expected 3 years, got {len(gene_jsd_hist)}"
    assert '0' in gene_jsd_hist, "Year 0 missing"
    assert '10' in gene_jsd_hist, "Year 10 missing"
    assert '20' in gene_jsd_hist, "Year 20 missing"
    
    # Check values
    assert gene_jsd_hist['0'] == [0.0] * 20, "Year 0 values incorrect"
    assert gene_jsd_hist['10'][0] == 0.1, "Year 10 first value incorrect"
    assert mean_hist['10'] == 0.11, "Year 10 mean incorrect"
    assert median_hist['20'] == 0.2, "Year 20 median incorrect"
    
    print("✓ Successfully extracted gene JSD from correct format")
    print(f"  - Found {len(gene_jsd_hist)} years of data")
    print(f"  - Years: {sorted(gene_jsd_hist.keys())}")
    print(f"  - Mean values: {mean_hist}")
    print(f"  - Median values: {median_hist}")
    
    return True


def test_legacy_format():
    """Test extraction from legacy format (top-level)"""
    print("\n" + "="*60)
    print("TEST 2: Legacy Format (top-level)")
    print("="*60)
    
    # Create test data in the legacy format
    sim_data = {
        "parameters": {"n": 100, "gene_size": 5},
        "gene_jsd_history": {
            "0": [0.0] * 20,
            "5": [0.05, 0.06, 0.04] + [0.05] * 17,
            "10": [0.1, 0.12, 0.08] + [0.1] * 17
        },
        "mean_gene_jsd_history": {
            "0": 0.0,
            "5": 0.05,
            "10": 0.1
        },
        "median_gene_jsd_history": {
            "0": 0.0,
            "5": 0.05,
            "10": 0.095
        },
        "history": {
            "0": {"cells": []},
            "5": {"cells": []},
            "10": {"cells": []}
        }
    }
    
    # Extract data
    gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(sim_data)
    
    # Verify extraction
    assert len(gene_jsd_hist) == 3, f"Expected 3 years, got {len(gene_jsd_hist)}"
    assert '0' in gene_jsd_hist, "Year 0 missing"
    assert '5' in gene_jsd_hist, "Year 5 missing"
    assert '10' in gene_jsd_hist, "Year 10 missing"
    
    # Check values
    assert gene_jsd_hist['5'][0] == 0.05, "Year 5 first value incorrect"
    assert mean_hist['5'] == 0.05, "Year 5 mean incorrect"
    assert median_hist['10'] == 0.095, "Year 10 median incorrect"
    
    print("✓ Successfully extracted gene JSD from legacy format")
    print(f"  - Found {len(gene_jsd_hist)} years of data")
    print(f"  - Years: {sorted(gene_jsd_hist.keys())}")
    
    return True


def test_mixed_format():
    """Test extraction when data exists in both locations (correct format takes precedence)"""
    print("\n" + "="*60)
    print("TEST 3: Mixed Format (both locations)")
    print("="*60)
    
    # Create test data with BOTH formats
    sim_data = {
        "parameters": {"n": 100, "gene_size": 5},
        # This is the OLD location (should be ignored if new exists)
        "gene_jsd_history": {
            "0": [0.999] * 20,  # Wrong values to test precedence
        },
        # This is the CORRECT location (should take precedence)
        "history": {
            "0": {
                "cells": [],
                "gene_jsd": [0.0] * 20,  # Correct values
                "mean_gene_jsd": 0.0,
                "median_gene_jsd": 0.0
            },
            "15": {
                "cells": [],
                "gene_jsd": [0.15] * 20,
                "mean_gene_jsd": 0.15,
                "median_gene_jsd": 0.15
            }
        }
    }
    
    # Extract data
    gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(sim_data)
    
    # Verify correct location takes precedence
    assert len(gene_jsd_hist) == 2, f"Expected 2 years from correct location, got {len(gene_jsd_hist)}"
    assert gene_jsd_hist['0'][0] == 0.0, "Should use correct location value (0.0), not legacy (0.999)"
    assert '15' in gene_jsd_hist, "Year 15 should be present from correct location"
    
    print("✓ Correct location takes precedence when both exist")
    print(f"  - Used data from history structure, not top-level")
    print(f"  - Years: {sorted(gene_jsd_hist.keys())}")
    
    return True


def test_no_gene_jsd():
    """Test extraction when no gene JSD data exists"""
    print("\n" + "="*60)
    print("TEST 4: No Gene JSD Data")
    print("="*60)
    
    # Create test data with NO gene JSD anywhere
    sim_data = {
        "parameters": {"n": 100, "gene_size": 5, "track_gene_jsd": False},
        "history": {
            "0": {"cells": []},
            "10": {"cells": []},
            "20": {"cells": []}
        }
    }
    
    # Extract data
    gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(sim_data)
    
    # Verify empty results
    assert len(gene_jsd_hist) == 0, f"Expected empty dict, got {len(gene_jsd_hist)} entries"
    assert len(mean_hist) == 0, "Mean history should be empty"
    assert len(median_hist) == 0, "Median history should be empty"
    
    print("✓ Correctly handles missing gene JSD data")
    print(f"  - Returned empty dictionaries")
    
    return True


def test_partial_data():
    """Test extraction when only some fields exist"""
    print("\n" + "="*60)
    print("TEST 5: Partial Data (gene_jsd without mean/median)")
    print("="*60)
    
    sim_data = {
        "parameters": {"n": 100, "gene_size": 5},
        "history": {
            "0": {
                "cells": [],
                "gene_jsd": [0.0] * 20
                # Note: no mean_gene_jsd or median_gene_jsd
            },
            "25": {
                "cells": [],
                "gene_jsd": [0.25] * 20,
                "mean_gene_jsd": 0.25  # Only mean, no median
            }
        }
    }
    
    # Extract data
    gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(sim_data)
    
    # Verify partial extraction
    assert len(gene_jsd_hist) == 2, "Should extract gene_jsd for both years"
    assert len(mean_hist) == 1, "Should only have mean for year 25"
    assert len(median_hist) == 0, "Should have no median data"
    assert '25' in mean_hist, "Year 25 should have mean"
    
    print("✓ Correctly handles partial data")
    print(f"  - Gene JSD years: {sorted(gene_jsd_hist.keys())}")
    print(f"  - Mean years: {sorted(mean_hist.keys())}")
    print(f"  - Median years: {sorted(median_hist.keys())}")
    
    return True


def test_real_simulation():
    """Test with a real simulation file if available"""
    print("\n" + "="*60)
    print("TEST 6: Real Simulation File")
    print("="*60)
    
    # Look for a real simulation file
    test_sim_path = "../phase1/data/gene_rates_5x0.00400_5x0.00500_5x0.00600_5x0.00700/size512-sites100-genesize5-years50-seed42-20250917-130044/simulation.json"
    
    if not os.path.exists(test_sim_path):
        print("  ⚠ No real simulation file found, skipping")
        return True
    
    try:
        # Load real simulation
        with open(test_sim_path, 'r') as f:
            sim_data = json.load(f)
        
        # Extract data
        gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(sim_data)
        
        # Check what we found
        if gene_jsd_hist:
            print(f"✓ Found gene JSD data in real simulation")
            print(f"  - Years with data: {sorted(gene_jsd_hist.keys())[:10]}...")
            print(f"  - Number of genes: {len(gene_jsd_hist['0']) if '0' in gene_jsd_hist else 'N/A'}")
            print(f"  - Has mean data: {'Yes' if mean_hist else 'No'}")
            print(f"  - Has median data: {'Yes' if median_hist else 'No'}")
        else:
            print(f"  ℹ Real simulation has no gene JSD data (was run with --no-gene-jsd)")
            print(f"    Parameters: {sim_data.get('parameters', {})}")
        
    except Exception as e:
        print(f"  ⚠ Error loading real simulation: {e}")
        return False
    
    return True


def test_integration_with_pipeline():
    """Test that the updated pipeline can process files correctly"""
    print("\n" + "="*60)
    print("TEST 7: Integration Test (mini pipeline run)")
    print("="*60)
    
    # Create a minimal simulation file
    with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False) as f:
        sim_data = {
            "parameters": {
                "n": 100,
                "gene_size": 5,
                "years": 20,
                "growth_phase": 3,
                "seed": 123,
                "track_gene_jsd": True
            },
            "history": {
                "0": {
                    "cells": [{"methylated": [0]*100, "cell_jsd": 0.0, "age": 0}] * 10,
                    "gene_jsd": [0.0] * 20,
                    "mean_gene_jsd": 0.0,
                    "median_gene_jsd": 0.0
                },
                "10": {
                    "cells": [{"methylated": [0]*50 + [1]*50, "cell_jsd": 0.5, "age": 10}] * 20,
                    "gene_jsd": [0.1] * 20,
                    "mean_gene_jsd": 0.1,
                    "median_gene_jsd": 0.1
                },
                "20": {
                    "cells": [{"methylated": [1]*100, "cell_jsd": 1.0, "age": 20}] * 30,
                    "gene_jsd": [0.2] * 20,
                    "mean_gene_jsd": 0.2,
                    "median_gene_jsd": 0.2
                }
            }
        }
        json.dump(sim_data, f)
        test_file = f.name
    
    try:
        # Test extraction directly
        with open(test_file, 'r') as f:
            loaded_data = json.load(f)
        
        gene_jsd_hist, mean_hist, median_hist = extract_gene_jsd_from_history(loaded_data)
        
        assert len(gene_jsd_hist) == 3, "Should extract all 3 years"
        assert gene_jsd_hist['10'][0] == 0.1, "Should extract correct values"
        
        print("✓ Integration test passed")
        print(f"  - Created test file: {test_file}")
        print(f"  - Successfully extracted gene JSD data")
        
    finally:
        # Clean up
        if os.path.exists(test_file):
            os.remove(test_file)
    
    return True


def run_all_tests():
    """Run all tests and report results"""
    print("\n" + "="*80)
    print("GENE JSD EXTRACTION TEST SUITE")
    print("="*80)
    
    tests = [
        ("Correct Format", test_correct_format),
        ("Legacy Format", test_legacy_format),
        ("Mixed Format", test_mixed_format),
        ("No Gene JSD", test_no_gene_jsd),
        ("Partial Data", test_partial_data),
        ("Real Simulation", test_real_simulation),
        ("Integration", test_integration_with_pipeline)
    ]
    
    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, success))
        except Exception as e:
            print(f"\n❌ {name} failed with error: {e}")
            print(traceback.format_exc())
            results.append((name, False))
    
    # Summary
    print("\n" + "="*80)
    print("TEST SUMMARY")
    print("="*80)
    
    passed = sum(1 for _, success in results if success)
    total = len(results)
    
    for name, success in results:
        status = "✓ PASSED" if success else "✗ FAILED"
        print(f"  {status}: {name}")
    
    print(f"\nTotal: {passed}/{total} tests passed")
    
    return passed == total


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1)