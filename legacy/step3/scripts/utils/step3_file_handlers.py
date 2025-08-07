"""Step 3 specific file handling utilities."""

import os
import glob
from typing import List, Dict, Tuple, Optional
from pathlib import Path


class Step3FileHandler:
    """Handles file operations specific to Step 3."""
    
    @staticmethod
    def find_completed_rates(step2_dir: str) -> List[str]:
        """
        Find all rates that have completed Step 2 processing.
        
        Args:
            step2_dir: Path to Step 2 data directory
            
        Returns:
            List of rate strings that have complete lineages
        """
        completed_rates = []
        
        # Look for rate directories
        rate_dirs = glob.glob(os.path.join(step2_dir, "rate_*"))
        
        for rate_dir in sorted(rate_dirs):
            # Extract rate from directory name
            rate = os.path.basename(rate_dir).replace("rate_", "")
            
            # Check if lineages exist
            mutant_files = glob.glob(os.path.join(rate_dir, "lineages", "mutant", "*.json.gz"))
            control_files = glob.glob(os.path.join(rate_dir, "lineages", "control", "*.json.gz"))
            
            if mutant_files and control_files:
                completed_rates.append(rate)
        
        return completed_rates
    
    @staticmethod
    def get_step2_paths(step2_dir: str, rate: str) -> Dict[str, str]:
        """
        Get all relevant Step 2 output paths for a given rate.
        
        Args:
            step2_dir: Path to Step 2 data directory
            rate: Rate string
            
        Returns:
            Dictionary of paths
        """
        rate_dir = os.path.join(step2_dir, f"rate_{rate}")
        
        return {
            'snapshot': os.path.join(rate_dir, 'snapshots', 'year50_snapshot.json.gz'),
            'lineages_mutant': os.path.join(rate_dir, 'lineages', 'mutant'),
            'lineages_control': os.path.join(rate_dir, 'lineages', 'control'),
            'plot': os.path.join(rate_dir, 'plots', 'year50_jsd_distribution.png')
        }
    
    @staticmethod
    def get_original_simulation_path(step1_dir: str, rate: str) -> Optional[str]:
        """
        Find the original simulation file for a given rate.
        
        Args:
            step1_dir: Path to Step 1 data directory
            rate: Rate string
            
        Returns:
            Path to simulation file or None if not found
        """
        pattern = os.path.join(step1_dir, f"simulation_rate_{rate}_*.json.gz")
        files = glob.glob(pattern)
        
        if files:
            return files[0]  # Should only be one
        return None
    
    @staticmethod
    def create_step3_directories(base_dir: str, rate: str) -> Dict[str, str]:
        """
        Create directory structure for Step 3 outputs.
        
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
            'individuals_mutant': os.path.join(rate_dir, 'individuals', 'mutant'),
            'individuals_control': os.path.join(rate_dir, 'individuals', 'control'),
            'plots': os.path.join(rate_dir, 'plots'),
            'results': os.path.join(rate_dir, 'results')
        }
        
        for path in dirs.values():
            os.makedirs(path, exist_ok=True)
        
        return dirs
    
    @staticmethod
    def create_combined_analysis_dirs(base_dir: str) -> Dict[str, str]:
        """
        Create directories for combined analysis across rates.
        
        Args:
            base_dir: Base output directory
            
        Returns:
            Dictionary of created directory paths
        """
        combined_dir = os.path.join(base_dir, 'combined_analysis')
        
        dirs = {
            'base': combined_dir,
            'plots': os.path.join(combined_dir, 'plots'),
            'results': os.path.join(combined_dir, 'results')
        }
        
        for path in dirs.values():
            os.makedirs(path, exist_ok=True)
        
        return dirs
    
    @staticmethod
    def check_step3_outputs_exist(dirs: Dict[str, str], year: int = 60) -> Dict[str, bool]:
        """
        Check which Step 3 outputs already exist.
        
        Args:
            dirs: Directory structure from create_step3_directories
            year: Year being extracted
            
        Returns:
            Dictionary indicating which outputs exist
        """
        results = {}
        
        # Check snapshot
        snapshot_path = os.path.join(dirs['snapshots'], f'year{year}_snapshot.json.gz')
        results['snapshot'] = os.path.exists(snapshot_path)
        
        # Check individuals
        results['individuals_mutant'] = bool(
            glob.glob(os.path.join(dirs['individuals_mutant'], '*.json.gz'))
        )
        results['individuals_control'] = bool(
            glob.glob(os.path.join(dirs['individuals_control'], '*.json.gz'))
        )
        
        # Check plots
        plot_path = os.path.join(dirs['plots'], 'jsd_distribution_comparison.png')
        results['plot'] = os.path.exists(plot_path)
        
        # Check results
        results_path = os.path.join(dirs['results'], 'jsd_distributions.json')
        results['statistics'] = os.path.exists(results_path)
        
        return results