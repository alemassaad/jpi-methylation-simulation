#!/usr/bin/env python3
"""
Checkpoint system for tracking pipeline progress.
Maintains a JSON file that records which stages have been completed.
"""

import json
import os
from datetime import datetime


class PipelineCheckpoint:
    """Manages pipeline execution checkpoints."""
    
    def __init__(self, checkpoint_file):
        """
        Initialize checkpoint manager.
        
        Args:
            checkpoint_file: Path to checkpoint JSON file
        """
        self.checkpoint_file = checkpoint_file
        self.checkpoints = self.load()
    
    def load(self):
        """Load existing checkpoints or create new."""
        if os.path.exists(self.checkpoint_file):
            with open(self.checkpoint_file, 'r') as f:
                return json.load(f)
        return {
            "created": datetime.now().isoformat(),
            "stages": {},
            "parameters": {},
            "files": {}
        }
    
    def save(self):
        """Save checkpoints to file."""
        dir_path = os.path.dirname(self.checkpoint_file)
        if dir_path:
            os.makedirs(dir_path, exist_ok=True)
        with open(self.checkpoint_file, 'w') as f:
            json.dump(self.checkpoints, f, indent=2)
    
    def mark_stage_complete(self, stage_name, details=None):
        """
        Mark a stage as complete.
        
        Args:
            stage_name: Name of the stage (e.g., "create_mutant", "grow_mutant")
            details: Optional details about the stage completion
        """
        self.checkpoints["stages"][stage_name] = {
            "completed": True,
            "timestamp": datetime.now().isoformat(),
            "details": details
        }
        self.save()
    
    def is_stage_complete(self, stage_name):
        """Check if a stage has been completed."""
        return self.checkpoints["stages"].get(stage_name, {}).get("completed", False)
    
    def set_parameter(self, key, value):
        """Store a parameter value."""
        self.checkpoints["parameters"][key] = value
        self.save()
    
    def get_parameter(self, key, default=None):
        """Get a stored parameter value."""
        return self.checkpoints["parameters"].get(key, default)
    
    def record_file(self, file_type, file_path, cell_count=None):
        """
        Record information about a file.
        
        Args:
            file_type: Type of file (e.g., "mutant_individual_00")
            file_path: Path to the file
            cell_count: Optional number of cells in the file
        """
        self.checkpoints["files"][file_type] = {
            "path": file_path,
            "cell_count": cell_count,
            "timestamp": datetime.now().isoformat()
        }
        self.save()
    
    def get_file_info(self, file_type):
        """Get information about a recorded file."""
        return self.checkpoints["files"].get(file_type, {})
    
    def check_growth_status(self, group, expected_cells):
        """
        Check if individuals in a group are grown.
        
        Args:
            group: "mutant" or "control1"
            expected_cells: Expected number of cells after growth
            
        Returns:
            True if all individuals are grown, False otherwise
        """
        stage_key = f"grow_{group}"
        if not self.is_stage_complete(stage_key):
            return False
        
        # Additional check: verify expected cell count
        stored_count = self.checkpoints["stages"][stage_key].get("details", {}).get("cell_count")
        return stored_count == expected_cells
    
    def check_mixing_status(self, group):
        """
        Check if individuals in a group are mixed.
        
        Args:
            group: "mutant" or "control1"
            
        Returns:
            True if mixing is complete, False otherwise
        """
        return self.is_stage_complete(f"mix_{group}")
    
    def reset(self):
        """Reset all checkpoints."""
        self.checkpoints = {
            "created": datetime.now().isoformat(),
            "stages": {},
            "parameters": {},
            "files": {}
        }
        self.save()
    
    def summary(self):
        """Get a summary of completed stages."""
        completed = [stage for stage, info in self.checkpoints["stages"].items() 
                    if info.get("completed")]
        return {
            "completed_stages": completed,
            "total_stages": len(self.checkpoints["stages"]),
            "parameters": self.checkpoints["parameters"],
            "file_count": len(self.checkpoints["files"])
        }


# Stage names for the pipeline
PIPELINE_STAGES = [
    "extract_year50",
    "plot_year50_jsd",
    "create_mutant",
    "create_control1",
    "grow_mutant",
    "grow_control1",
    "extract_year60",
    "plot_year60_jsd",
    "mix_mutant",
    "mix_control1",
    "create_control2",
    "analysis"
]


def demonstrate_usage():
    """Demonstrate how to use the checkpoint system."""
    
    # Create checkpoint manager
    checkpoint = PipelineCheckpoint("test_checkpoint.json")
    
    # Store parameters
    checkpoint.set_parameter("rate", 0.005)
    checkpoint.set_parameter("n_individuals", 30)
    checkpoint.set_parameter("growth_years", 10)
    
    # Simulate pipeline stages
    print("Simulating pipeline with checkpoints...\n")
    
    # Stage 1: Extract year 50
    if checkpoint.is_stage_complete("extract_year50"):
        print("⏭ Skipping year 50 extraction (already complete)")
    else:
        print("✓ Extracting year 50 snapshot...")
        checkpoint.mark_stage_complete("extract_year50", {"cells": 10000})
    
    # Stage 2: Create mutant individuals
    if checkpoint.is_stage_complete("create_mutant"):
        print("⏭ Skipping mutant creation (already complete)")
    else:
        print("✓ Creating mutant individuals...")
        checkpoint.mark_stage_complete("create_mutant", {"count": 30})
        for i in range(3):  # Record a few files
            checkpoint.record_file(f"mutant_{i:02d}", f"individuals/mutant/individual_{i:02d}.json.gz", 1)
    
    # Stage 3: Grow mutant individuals
    if checkpoint.check_growth_status("mutant", 1024):
        print("⏭ Skipping mutant growth (already at 1024 cells)")
    else:
        print("✓ Growing mutant individuals...")
        checkpoint.mark_stage_complete("grow_mutant", {"cell_count": 1024})
    
    # Print summary
    print("\n" + "="*60)
    print("CHECKPOINT SUMMARY")
    print("="*60)
    summary = checkpoint.summary()
    print(f"Completed stages: {', '.join(summary['completed_stages'])}")
    print(f"Total stages: {summary['total_stages']}")
    print(f"Parameters: {summary['parameters']}")
    print(f"Files recorded: {summary['file_count']}")
    
    # Clean up
    os.remove("test_checkpoint.json")
    print("\nTest checkpoint file removed.")


if __name__ == "__main__":
    demonstrate_usage()