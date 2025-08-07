"""Process management utilities for running external scripts."""

import subprocess
import os
import time
from typing import Tuple, Optional, List, Dict
from .logger import get_logger


class ProcessManager:
    """Manages subprocess execution of existing scripts."""
    
    def __init__(self, script_dir: str = "."):
        """
        Initialize ProcessManager.
        
        Args:
            script_dir: Directory containing the scripts to run
        """
        self.script_dir = script_dir
        self.logger = get_logger()
    
    def run_command(
        self, 
        cmd: List[str], 
        description: str,
        timeout: Optional[int] = None,
        cwd: Optional[str] = None
    ) -> Tuple[bool, str, float]:
        """
        Execute a command with proper error handling.
        
        Args:
            cmd: Command as list of strings
            description: Description for logging
            timeout: Timeout in seconds (None for no timeout)
            cwd: Working directory for the command
            
        Returns:
            Tuple of (success, output/error, duration)
        """
        start_time = time.time()
        
        try:
            self.logger.debug(f"Running: {' '.join(cmd)}")
            
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=timeout,
                cwd=cwd or self.script_dir
            )
            
            duration = time.time() - start_time
            
            if result.returncode == 0:
                return True, result.stdout, duration
            else:
                error_msg = result.stderr or result.stdout or "Unknown error"
                self.logger.error(f"{description} failed: {error_msg}")
                return False, error_msg, duration
                
        except subprocess.TimeoutExpired:
            duration = time.time() - start_time
            error_msg = f"Timeout after {timeout} seconds"
            self.logger.error(f"{description} timeout: {error_msg}")
            return False, error_msg, duration
            
        except Exception as e:
            duration = time.time() - start_time
            error_msg = str(e)
            self.logger.error(f"{description} exception: {error_msg}")
            return False, error_msg, duration
    
    def extract_snapshot(
        self,
        sim_file: str,
        year: int,
        output_path: str,
        timeout: int = 300
    ) -> Tuple[bool, str, float]:
        """
        Run extract_snapshot.py script.
        
        Args:
            sim_file: Path to simulation file
            year: Year to extract
            output_path: Output path for snapshot
            timeout: Timeout in seconds
            
        Returns:
            Tuple of (success, output/error, duration)
        """
        cmd = [
            'python', 'extract_snapshot.py',
            sim_file,
            str(year),
            output_path
        ]
        
        return self.run_command(
            cmd,
            f"Extracting year {year} snapshot",
            timeout=timeout
        )
    
    def create_lineages(
        self,
        snapshot_path: str,
        output_dir: str,
        lineage_type: str = "both",
        timeout: int = 3600
    ) -> Tuple[bool, str, float]:
        """
        Run create_lineages.py script.
        
        Args:
            snapshot_path: Path to snapshot file
            output_dir: Output directory for lineages
            lineage_type: Type of lineages to create ("mutant", "control", "both")
            timeout: Timeout in seconds
            
        Returns:
            Tuple of (success, output/error, duration)
        """
        # create_lineages.py expects snapshot file as positional argument
        cmd = [
            'python', 'create_lineages.py',
            snapshot_path,  # Positional argument
            '--type', lineage_type,
            '--output-dir', output_dir
        ]
        
        return self.run_command(
            cmd,
            f"Creating {lineage_type} lineages",
            timeout=timeout
        )
    
    def plot_distribution(
        self,
        snapshot_path: str,
        output_path: str,
        bins: int = 200,
        timeout: int = 300
    ) -> Tuple[bool, str, float]:
        """
        Run plot_jsd_distribution.py script.
        
        Args:
            snapshot_path: Path to snapshot file
            output_path: Output path for plot
            bins: Number of bins for histogram
            timeout: Timeout in seconds
            
        Returns:
            Tuple of (success, output/error, duration)
        """
        # Extract directory and filename
        output_dir = os.path.dirname(output_path)
        output_name = os.path.splitext(os.path.basename(output_path))[0]
        
        cmd = [
            'python', 'plot_jsd_distribution.py',
            snapshot_path,
            str(bins)
        ]
        
        # The script saves to a default location, we'll need to move it
        success, output, duration = self.run_command(
            cmd,
            f"Plotting JSD distribution",
            timeout=timeout
        )
        
        if success:
            # Move the generated plot to the desired location
            # The script generates files like: year50_snapshot_jsd_distribution_200bins.png
            snapshot_base = os.path.splitext(os.path.basename(snapshot_path))[0]
            generated_plot = os.path.join(
                self.script_dir, '..', 'plots',
                f'{snapshot_base}_jsd_distribution_{bins}bins.png'
            )
            
            if os.path.exists(generated_plot):
                os.makedirs(output_dir, exist_ok=True)
                import shutil
                shutil.move(generated_plot, output_path)
        
        return success, output, duration
    
    def check_script_exists(self, script_name: str) -> bool:
        """Check if a script exists in the script directory."""
        script_path = os.path.join(self.script_dir, script_name)
        return os.path.exists(script_path)
    
    def get_python_version(self) -> str:
        """Get the Python version being used."""
        result = subprocess.run(
            ['python', '--version'],
            capture_output=True,
            text=True
        )
        return result.stdout.strip()