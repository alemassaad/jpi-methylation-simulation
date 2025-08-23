#!/usr/bin/env python3
"""
Master test runner for uniform mixing implementation.
Runs all test suites and reports comprehensive results.
"""

import sys
import os
import time
import subprocess
from typing import Dict, List, Tuple

def run_test_file(test_file: str) -> Tuple[bool, str, float]:
    """
    Run a single test file and capture results.
    Returns (success, output, duration).
    """
    start_time = time.time()
    
    try:
        result = subprocess.run(
            [sys.executable, test_file],
            capture_output=True,
            text=True,
            timeout=30  # 30 second timeout per test file
        )
        
        duration = time.time() - start_time
        output = result.stdout + result.stderr
        success = result.returncode == 0
        
        return success, output, duration
        
    except subprocess.TimeoutExpired:
        duration = time.time() - start_time
        return False, f"TIMEOUT after {duration:.1f} seconds", duration
        
    except Exception as e:
        duration = time.time() - start_time
        return False, f"ERROR: {str(e)}", duration


def extract_test_results(output: str) -> Dict[str, int]:
    """Extract pass/fail counts from test output."""
    results = {"passed": 0, "failed": 0}
    
    # Look for various result patterns
    patterns = [
        "RESULTS: ",
        "passed, ",
        "Test passed",
        "✓",
        "✗"
    ]
    
    # Count check marks
    results["passed"] = output.count("✓")
    results["failed"] = output.count("✗")
    
    # Also look for explicit counts
    for line in output.split('\n'):
        if "passed" in line and "failed" in line:
            try:
                parts = line.split()
                for i, part in enumerate(parts):
                    if "passed" in part and i > 0:
                        results["passed"] = max(results["passed"], int(parts[i-1]))
                    if "failed" in part and i > 0:
                        results["failed"] = max(results["failed"], int(parts[i-1]))
            except:
                pass
    
    return results


def main():
    """Run all test suites and report results."""
    print("="*70)
    print("COMPREHENSIVE TEST SUITE FOR UNIFORM MIXING IMPLEMENTATION")
    print("="*70)
    print()
    
    # Define test files
    test_files = [
        ("Unit Tests", "test_uniform_mixing_units.py"),
        ("Edge Cases", "test_uniform_mixing_edge_cases.py"),
        ("Comparison", "test_uniform_mixing_comparison.py"),
    ]
    
    # Overall statistics
    total_passed = 0
    total_failed = 0
    total_duration = 0
    all_success = True
    
    # Run each test suite
    for test_name, test_file in test_files:
        print(f"\n{'='*70}")
        print(f"Running {test_name}: {test_file}")
        print("-"*70)
        
        test_path = os.path.join(os.path.dirname(__file__), test_file)
        
        if not os.path.exists(test_path):
            print(f"  ⚠️  Test file not found: {test_file}")
            continue
        
        success, output, duration = run_test_file(test_path)
        total_duration += duration
        
        # Extract detailed results
        results = extract_test_results(output)
        total_passed += results["passed"]
        total_failed += results["failed"]
        
        if not success:
            all_success = False
        
        # Print summary for this test
        status_icon = "✅" if success else "❌"
        print(f"\n{status_icon} {test_name} Results:")
        print(f"   Duration: {duration:.2f} seconds")
        print(f"   Tests passed: {results['passed']}")
        print(f"   Tests failed: {results['failed']}")
        
        # Show key output lines
        print("\n   Key output:")
        for line in output.split('\n'):
            if any(keyword in line for keyword in ["✓", "✗", "PASS", "FAIL", "ERROR", "Test"]):
                print(f"     {line.strip()}")
                if len(line.strip()) > 100:  # Truncate long lines
                    break
    
    # Final summary
    print("\n" + "="*70)
    print("FINAL TEST SUMMARY")
    print("="*70)
    
    print(f"\nOverall Results:")
    print(f"  Total Tests Passed: {total_passed}")
    print(f"  Total Tests Failed: {total_failed}")
    print(f"  Total Duration: {total_duration:.2f} seconds")
    print(f"  Average per Test: {total_duration/len(test_files):.2f} seconds")
    
    if all_success:
        print("\n✅ ALL TEST SUITES PASSED! ✅")
        print("\nThe uniform mixing implementation is working correctly:")
        print("  • Core functions return correct data structures")
        print("  • Edge cases are handled gracefully")
        print("  • Statistical properties are as expected")
        print("  • All three groups share the same base cells")
        print("  • Control2 has unique additional cells per individual")
    else:
        print("\n❌ SOME TESTS FAILED ❌")
        print("\nIssues detected in the implementation.")
        print("Please review the failed tests above.")
    
    print("\n" + "="*70)
    
    return 0 if all_success else 1


if __name__ == "__main__":
    sys.exit(main())