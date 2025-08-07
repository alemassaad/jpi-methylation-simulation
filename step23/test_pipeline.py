#!/usr/bin/env python3
"""Test the pipeline with existing snapshots."""

import json
import gzip
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from pipeline_utils import sample_uniform

print("Loading year 50 snapshot...")
with gzip.open('data/rate_0.005000/snapshots/year50_snapshot.json.gz', 'rt') as f:
    snapshot_data = json.load(f)
year50_cells = snapshot_data['cells']
print(f"Loaded {len(year50_cells)} cells")

print("\nSampling 4 cells uniformly...")
sampled = sample_uniform(year50_cells, n_cells=4, seed=42)
print(f"Sampled {len(sampled)} cells")

for i, cell in enumerate(sampled):
    print(f"  Cell {i}: JSD = {cell['jsd']:.4f}, age = {cell['age']}")

print("\nTest complete!")