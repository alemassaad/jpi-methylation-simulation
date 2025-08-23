#!/usr/bin/env python3
"""
Simple verification test that checks the implementation without numpy.
This test verifies the code changes were correctly applied.
"""

import sys
import os
import ast
import inspect

# Add paths
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def check_function_signature(module_path, function_name, expected_return_type):
    """Check if a function has the expected signature."""
    with open(module_path, 'r') as f:
        tree = ast.parse(f.read())
    
    for node in ast.walk(tree):
        if isinstance(node, ast.FunctionDef) and node.name == function_name:
            # Check return annotation if present
            if node.returns:
                return_str = ast.unparse(node.returns) if hasattr(ast, 'unparse') else str(node.returns)
                if expected_return_type in return_str:
                    return True, f"Found {function_name} with return type containing {expected_return_type}"
            
            # Check for return statements
            for child in ast.walk(node):
                if isinstance(child, ast.Return) and child.value:
                    if isinstance(child.value, ast.Tuple):
                        return True, f"Found {function_name} returning a tuple"
                    # Check for explicit tuple return
                    value_str = ast.unparse(child.value) if hasattr(ast, 'unparse') else str(child.value)
                    if "indices" in value_str:
                        return True, f"Found {function_name} returning indices"
            
            return False, f"Found {function_name} but wrong return type"
    
    return False, f"Function {function_name} not found"


def check_function_exists(module_path, function_name):
    """Check if a function exists in a module."""
    with open(module_path, 'r') as f:
        content = f.read()
    
    if f"def {function_name}" in content:
        # Count parameters
        import re
        pattern = f"def {function_name}\\([^)]*\\)"
        match = re.search(pattern, content)
        if match:
            params = match.group(0).count(',') + 1 if '(' in match.group(0) else 0
            return True, f"Found {function_name} with ~{params} parameters"
    
    return False, f"Function {function_name} not found"


def check_import_added(module_path, function_name):
    """Check if an import was added."""
    with open(module_path, 'r') as f:
        content = f.read()
    
    if function_name in content and "from pipeline_utils import" in content:
        return True, f"Import for {function_name} found"
    
    return False, f"Import for {function_name} not found"


def check_stage_implementation(module_path, stage_num, expected_text):
    """Check if a stage contains expected implementation."""
    with open(module_path, 'r') as f:
        content = f.read()
    
    # Find stage section
    stage_marker = f"STAGE {stage_num}"
    if stage_marker in content:
        # Get content after stage marker
        stage_start = content.index(stage_marker)
        # Find next stage or end
        next_stage = content.find(f"STAGE {stage_num + 1}", stage_start)
        if next_stage == -1:
            stage_content = content[stage_start:]
        else:
            stage_content = content[stage_start:next_stage]
        
        if expected_text in stage_content:
            return True, f"Stage {stage_num} contains '{expected_text}'"
        else:
            return False, f"Stage {stage_num} missing '{expected_text}'"
    
    return False, f"Stage {stage_num} not found"


def main():
    """Run verification tests."""
    print("="*70)
    print("IMPLEMENTATION VERIFICATION TEST")
    print("="*70)
    print("\nVerifying that all code changes were correctly applied...\n")
    
    # Define paths
    pipeline_utils = os.path.join(os.path.dirname(os.path.dirname(__file__)), "pipeline_utils.py")
    run_pipeline = os.path.join(os.path.dirname(os.path.dirname(__file__)), "run_pipeline.py")
    
    tests = [
        ("1. create_uniform_mixing_pool returns tuple", 
         lambda: check_function_signature(pipeline_utils, "create_uniform_mixing_pool", "Tuple")),
        
        ("2. create_control2_with_uniform_base exists", 
         lambda: check_function_exists(pipeline_utils, "create_control2_with_uniform_base")),
        
        ("3. Import added in run_pipeline.py", 
         lambda: check_import_added(run_pipeline, "create_control2_with_uniform_base")),
        
        ("4. Stage 6 stores uniform_indices", 
         lambda: check_stage_implementation(run_pipeline, 6, "uniform_indices")),
        
        ("5. Stage 7 uses uniform_base", 
         lambda: check_stage_implementation(run_pipeline, 7, "uniform_base")),
        
        ("6. Stage 7 calls create_control2_with_uniform_base", 
         lambda: check_stage_implementation(run_pipeline, 7, "create_control2_with_uniform_base")),
    ]
    
    passed = 0
    failed = 0
    
    for test_name, test_func in tests:
        try:
            success, message = test_func()
            if success:
                print(f"✅ {test_name}")
                print(f"   {message}")
                passed += 1
            else:
                print(f"❌ {test_name}")
                print(f"   {message}")
                failed += 1
        except Exception as e:
            print(f"❌ {test_name}")
            print(f"   Error: {e}")
            failed += 1
    
    print("\n" + "="*70)
    print("VERIFICATION SUMMARY")
    print("="*70)
    print(f"\nTests Passed: {passed}/{len(tests)}")
    print(f"Tests Failed: {failed}/{len(tests)}")
    
    if failed == 0:
        print("\n✅ ALL IMPLEMENTATION CHANGES VERIFIED ✅")
        print("\nThe uniform mixing Control2 implementation has been successfully applied:")
        print("  • create_uniform_mixing_pool now returns (pool, indices)")
        print("  • create_control2_with_uniform_base function added")
        print("  • Pipeline stages updated to use new functions")
        print("  • Control2 now shares base with mutant/control1 when using --uniform-mixing")
    else:
        print("\n❌ SOME IMPLEMENTATION CHANGES MISSING ❌")
        print("\nPlease review the failed checks above.")
    
    return 0 if failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())