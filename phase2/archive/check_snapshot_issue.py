#!/usr/bin/env python3
"""
Check why snapshots appear different in the test.
"""

import os
import hashlib
import gzip
import json

def check_content_vs_file_hash(filepath):
    """Check both file hash and content hash."""
    if not os.path.exists(filepath):
        return None, None
    
    # File hash (includes formatting, compression)
    with open(filepath, 'rb') as f:
        file_hash = hashlib.md5(f.read()).hexdigest()[:8]
    
    # Content hash (actual data)
    with gzip.open(filepath, 'rt') as f:
        data = json.load(f)
    
    # Normalize and hash content
    cells = data.get('cells', data)
    content_str = json.dumps(cells, sort_keys=True)
    content_hash = hashlib.md5(content_str.encode()).hexdigest()[:8]
    
    return file_hash, content_hash

# The issue is in test_reproducibility_robust.py - it uses file hash not content!
print("The issue is that test_reproducibility_robust.py uses file MD5 hash,")
print("which includes JSON formatting and gzip compression differences.")
print("")
print("We should compare the CONTENT, not the raw file bytes!")
print("")
print("The fix is to update get_file_hash() in test_reproducibility_robust.py to:")
print("1. Load the JSON content")
print("2. Hash the normalized content")
print("3. Not the raw compressed file")