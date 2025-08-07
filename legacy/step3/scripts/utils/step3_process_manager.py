"""Process management for Step 3 scripts."""

import os
import time
from typing import Tuple, Optional, List, Dict
from .process_manager import ProcessManager


class Step3ProcessManager(ProcessManager):
    """Manages subprocess execution of Step 3 scripts."""
    
    def extract_year60(
        self,
        sim_file: str,
        output_path: str,
        timeout: int = 300
    ) -> Tuple[bool, str, float]:
        """
        Run extract_year60_original.py script.
        
        Args:
            sim_file: Path to original simulation file
            output_path: Output path for year 60 snapshot
            timeout: Timeout in seconds
            
        Returns:
            Tuple of (success, output/error, duration)
        """
        # extract_year60_original.py expects positional arguments
        cmd = [
            'python', 'extract_year60_original.py',
            sim_file,      # First positional argument
            output_path    # Second positional argument
        ]
        
        return self.run_command(
            cmd,
            f"Extracting year 60 snapshot",
            timeout=timeout
        )
    
    def create_individuals(
        self,
        year60_snapshot: str,
        lineage_dir: str,
        output_dir: str,
        individual_type: str = "both",
        num_individuals: int = 30,
        mixture_ratio: float = 0.2,
        timeout: int = 600
    ) -> Tuple[bool, str, float]:
        """
        Run create_individuals.py script.
        
        Args:
            year60_snapshot: Path to year 60 snapshot
            lineage_dir: Directory containing lineages from Step 2
            output_dir: Output directory for individuals
            individual_type: Type to create ("mutant", "control", "both")
            num_individuals: Number of individuals to create
            mixture_ratio: Proportion of lineage cells
            timeout: Timeout in seconds
            
        Returns:
            Tuple of (success, output/error, duration)
        """
        cmd = [
            'python', 'create_individuals.py',
            '--type', individual_type,
            '--year60-snapshot', year60_snapshot,
            '--lineage-base-dir', lineage_dir,
            '--output-base-dir', output_dir,
            '--ratio', str(1 - mixture_ratio)  # Script expects fraction of original cells
        ]
        
        return self.run_command(
            cmd,
            f"Creating {individual_type} individuals",
            timeout=timeout
        )
    
    def plot_distributions(
        self,
        mutant_dir: str,
        control_dir: str,
        output_dir: str,
        timeout: int = 300
    ) -> Tuple[bool, str, float]:
        """
        Run plot_distributions.py script.
        
        Args:
            mutant_dir: Directory containing mutant individuals
            control_dir: Directory containing control individuals
            output_dir: Output directory for plots and results
            timeout: Timeout in seconds
            
        Returns:
            Tuple of (success, output/error, duration)
        """
        cmd = [
            'python', 'plot_distributions.py',
            '--mutant-dir', mutant_dir,
            '--control-dir', control_dir,
            '--output-dir', output_dir
        ]
        
        return self.run_command(
            cmd,
            "Plotting distributions",
            timeout=timeout
        )
    
    def check_script_dependencies(self) -> Dict[str, bool]:
        """Check if all required Step 3 scripts exist."""
        scripts = {
            'extract_year60_original.py': 'Extract year 60 snapshots',
            'create_individuals.py': 'Create mixed individuals',
            'plot_distributions.py': 'Plot comparison distributions'
        }
        
        results = {}
        for script, description in scripts.items():
            exists = self.check_script_exists(script)
            results[script] = exists
            if exists:
                self.logger.info(f"  ├── Found: {script} ({description})")
            else:
                self.logger.error(f"  ├── Missing: {script} ({description})")
        
        return results