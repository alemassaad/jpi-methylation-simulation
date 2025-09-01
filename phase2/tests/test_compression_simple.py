#!/usr/bin/env python3
"""
Simple test of compression logic without running full pipeline.
"""

import os
import sys
import argparse

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

def test_compression_logic():
    """Test the compression decision logic."""
    
    class Args:
        """Mock args object."""
        pass
    
    print("Testing Compression Logic")
    print("=" * 50)
    
    # Test 1: .json input, no flag
    print("\nTest 1: .json input, no --no-compress flag")
    args = Args()
    args.simulation = "test.json"
    
    # Simulate the new logic
    if hasattr(args, 'no_compress') and args.no_compress:
        args.use_compression = False
        source = "from --no-compress flag"
    elif args.simulation.endswith('.gz'):
        args.use_compression = True
        source = "matching compressed input"
    elif args.simulation.endswith('.json'):
        args.use_compression = False
        source = "matching uncompressed input"
    else:
        args.use_compression = True
        source = "default"
    
    print(f"  Input: {args.simulation}")
    print(f"  Result: use_compression = {args.use_compression} [{source}]")
    assert args.use_compression == False, "Should match uncompressed input"
    print("  ✅ Pass")
    
    # Test 2: .json.gz input, no flag
    print("\nTest 2: .json.gz input, no --no-compress flag")
    args = Args()
    args.simulation = "test.json.gz"
    
    if hasattr(args, 'no_compress') and args.no_compress:
        args.use_compression = False
        source = "from --no-compress flag"
    elif args.simulation.endswith('.gz'):
        args.use_compression = True
        source = "matching compressed input"
    elif args.simulation.endswith('.json'):
        args.use_compression = False
        source = "matching uncompressed input"
    else:
        args.use_compression = True
        source = "default"
    
    print(f"  Input: {args.simulation}")
    print(f"  Result: use_compression = {args.use_compression} [{source}]")
    assert args.use_compression == True, "Should match compressed input"
    print("  ✅ Pass")
    
    # Test 3: .json.gz input with --no-compress flag
    print("\nTest 3: .json.gz input WITH --no-compress flag")
    args = Args()
    args.simulation = "test.json.gz"
    args.no_compress = True
    
    if hasattr(args, 'no_compress') and args.no_compress:
        args.use_compression = False
        source = "from --no-compress flag"
    elif args.simulation.endswith('.gz'):
        args.use_compression = True
        source = "matching compressed input"
    elif args.simulation.endswith('.json'):
        args.use_compression = False
        source = "matching uncompressed input"
    else:
        args.use_compression = True
        source = "default"
    
    print(f"  Input: {args.simulation}")
    print(f"  Flag: --no-compress")
    print(f"  Result: use_compression = {args.use_compression} [{source}]")
    assert args.use_compression == False, "Flag should override input format"
    print("  ✅ Pass")
    
    # Test 4: Extension calculation
    print("\nTest 4: File extension calculation")
    for use_comp, expected_ext in [(True, ".json.gz"), (False, ".json")]:
        args.use_compression = use_comp
        snapshot_ext = ".json.gz" if args.use_compression else ".json"
        individual_ext = ".json.gz" if args.use_compression else ".json"
        print(f"  use_compression={use_comp} -> ext={snapshot_ext}")
        assert snapshot_ext == expected_ext
        assert individual_ext == expected_ext
    print("  ✅ Pass")
    
    print("\n" + "=" * 50)
    print("All compression logic tests passed! 🎉")

if __name__ == "__main__":
    test_compression_logic()