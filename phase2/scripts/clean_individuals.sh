#!/bin/bash
# Clean only individuals, keeping snapshots cached

if [ $# -eq 0 ]; then
    echo "Usage: ./clean_individuals.sh <rate>"
    echo "Example: ./clean_individuals.sh 0.005"
    exit 1
fi

RATE=$1
RATE_DIR="data/rate_${RATE}000"

echo "Cleaning individuals for rate $RATE..."

# Remove only individuals directories
if [ -d "$RATE_DIR/individuals" ]; then
    echo "Removing $RATE_DIR/individuals/"
    rm -rf "$RATE_DIR/individuals"
    echo "✓ Individuals cleaned"
else
    echo "No individuals directory found at $RATE_DIR/individuals"
fi

# Keep these:
echo ""
echo "Preserved cached data:"
if [ -d "$RATE_DIR/snapshots" ]; then
    echo "✓ Snapshots preserved:"
    ls -lh "$RATE_DIR/snapshots/"
fi

if [ -d "$RATE_DIR/plots" ]; then
    echo "✓ Plots preserved"
fi

if [ -d "$RATE_DIR/results" ]; then
    echo "✓ Results preserved" 
fi

echo ""
echo "Ready to re-run pipeline. Snapshots will be reused!"