"""File handling utilities for batch processing."""

import os
import re
import glob
import json
import shutil
from typing import List, Dict, Tuple, Optional
from pathlib import Path


class FileHandler:
    """Handles file operations for batch processing."""
    
    @staticmethod
    def find_simulation_files(input_dir: str, pattern: str = "simulation_rate_*.json.gz") -> List[str]:
        """
        Find all simulation files matching the pattern.
        
        Args:
            input_dir: Directory to search in
            pattern: Glob pattern for simulation files
            
        Returns:
            List of absolute paths to simulation files
        """
        search_path = os.path.join(input_dir, pattern)
        files = glob.glob(search_path)
        return sorted([os.path.abspath(f) for f in files])
    
    @staticmethod
    def extract_rate_from_filename(filename: str) -> Optional[str]:
        """
        Extract rate value from simulation filename.
        
        Args:
            filename: Simulation filename
            
        Returns:
            Rate as string (e.g., "0.005000") or None if not found
        """
        basename = os.path.basename(filename)
        match = re.search(r'simulation_rate_(\d+\.\d+)_', basename)
        if match:
            return match.group(1)
        return None
    
    @staticmethod
    def create_rate_directories(base_dir: str, rate: str) -> Dict[str, str]:
        """
        Create directory structure for a specific rate.
        
        Args:
            base_dir: Base output directory
            rate: Rate value as string
            
        Returns:
            Dictionary of created directory paths
        """
        rate_dir = os.path.join(base_dir, f"rate_{rate}")
        
        dirs = {
            'base': rate_dir,
            'snapshots': os.path.join(rate_dir, 'snapshots'),
            'lineages_mutant': os.path.join(rate_dir, 'lineages', 'mutant'),
            'lineages_control': os.path.join(rate_dir, 'lineages', 'control'),
            'plots': os.path.join(rate_dir, 'plots')
        }
        
        for path in dirs.values():
            os.makedirs(path, exist_ok=True)
        
        return dirs
    
    @staticmethod
    def check_disk_space(path: str, required_gb: float = 20.0) -> Tuple[bool, float]:
        """
        Check if there's enough disk space.
        
        Args:
            path: Path to check
            required_gb: Required space in GB
            
        Returns:
            Tuple of (has_enough_space, available_gb)
        """
        stat = shutil.disk_usage(path)
        available_gb = stat.free / (1024 ** 3)
        return available_gb >= required_gb, available_gb
    
    @staticmethod
    def save_state(state_file: str, state: Dict):
        """Save processing state to JSON file."""
        temp_file = f"{state_file}.tmp"
        with open(temp_file, 'w') as f:
            json.dump(state, f, indent=2)
        # Atomic rename
        os.replace(temp_file, state_file)
    
    @staticmethod
    def load_state(state_file: str) -> Dict:
        """Load processing state from JSON file."""
        if os.path.exists(state_file):
            with open(state_file, 'r') as f:
                return json.load(f)
        return {
            'completed': {},
            'failed': {},
            'in_progress': {},
            'last_run': None
        }
    
    @staticmethod
    def get_output_paths(dirs: Dict[str, str], year: int) -> Dict[str, str]:
        """
        Generate output file paths for a specific rate.
        
        Args:
            dirs: Directory structure from create_rate_directories
            year: Year being extracted
            
        Returns:
            Dictionary of output file paths
        """
        return {
            'snapshot': os.path.join(dirs['snapshots'], f'year{year}_snapshot.json.gz'),
            'plot': os.path.join(dirs['plots'], f'year{year}_jsd_distribution.png'),
            'lineages_mutant': dirs['lineages_mutant'],
            'lineages_control': dirs['lineages_control']
        }
    
    @staticmethod
    def check_outputs_exist(paths: Dict[str, str]) -> Dict[str, bool]:
        """
        Check which outputs already exist.
        
        Args:
            paths: Dictionary of output paths
            
        Returns:
            Dictionary indicating which outputs exist
        """
        results = {}
        
        # Check snapshot
        results['snapshot'] = os.path.exists(paths['snapshot'])
        
        # Check plot
        results['plot'] = os.path.exists(paths['plot'])
        
        # Check lineages (at least one file in each directory)
        results['lineages_mutant'] = bool(
            glob.glob(os.path.join(paths['lineages_mutant'], '*.json.gz'))
        )
        results['lineages_control'] = bool(
            glob.glob(os.path.join(paths['lineages_control'], '*.json.gz'))
        )
        
        return results