#!/usr/bin/env python3
"""
Batch processor for Step 2 of the methylation simulation pipeline.

This script processes multiple simulation files from Step 1, extracting snapshots,
creating lineages, and generating plots for each methylation rate.
"""

import os
import sys
import yaml
import time
import argparse
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple

from utils import setup_logger, get_logger, FileHandler, ProcessManager
from utils.logger import ProgressLogger


class BatchProcessor:
    """Orchestrates batch processing of simulation files."""
    
    def __init__(self, config_path: str = "config/batch_config.yaml"):
        """
        Initialize BatchProcessor with configuration.
        
        Args:
            config_path: Path to configuration file
        """
        self.config = self._load_config(config_path)
        self.logger = setup_logger(
            log_dir=self.config['logging']['log_dir'],
            level=self.config['logging']['level'],
            console=self.config['logging']['console']
        )
        self.progress_logger = ProgressLogger(self.logger)
        self.file_handler = FileHandler()
        self.process_manager = ProcessManager(script_dir=os.getcwd())
        
        # Load or initialize state
        self.state_file = self.config['paths']['state_file']
        self.state = self.file_handler.load_state(self.state_file)
        
        # Track current run
        self.start_time = None
        self.processed_count = 0
        self.failed_count = 0
    
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration from YAML file."""
        if not os.path.exists(config_path):
            raise FileNotFoundError(f"Configuration file not found: {config_path}")
        
        with open(config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def run(self, 
            rates: Optional[List[str]] = None,
            dry_run: bool = False,
            resume: bool = False) -> bool:
        """
        Run batch processing.
        
        Args:
            rates: Specific rates to process (None for all)
            dry_run: If True, only show what would be done
            resume: If True, skip completed simulations
            
        Returns:
            True if all simulations processed successfully
        """
        self.start_time = datetime.now()
        self.logger.info("="*60)
        self.logger.info("BATCH PROCESSING STARTED")
        self.logger.info(f"Start time: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("="*60)
        
        try:
            # Phase 1: Discovery
            simulation_files = self._discover_simulations(rates)
            if not simulation_files:
                self.logger.warning("No simulation files found!")
                return False
            
            # Phase 2: Validation
            if not self._validate_environment(dry_run):
                return False
            
            # Phase 3: Processing
            success = self._process_simulations(simulation_files, dry_run, resume)
            
            # Phase 4: Reporting
            self._generate_report()
            
            return success
            
        except KeyboardInterrupt:
            self.logger.warning("\nBatch processing interrupted by user!")
            self._save_state()
            return False
            
        except Exception as e:
            self.logger.error(f"Unexpected error: {e}", exc_info=True)
            self._save_state()
            return False
    
    def _discover_simulations(self, rates: Optional[List[str]] = None) -> List[Tuple[str, str]]:
        """
        Discover simulation files to process.
        
        Returns:
            List of (file_path, rate) tuples
        """
        self.logger.info("Phase 1: Discovery")
        self.progress_logger.indent()
        
        # Find all simulation files
        input_dir = self.config['paths']['input_dir']
        all_files = self.file_handler.find_simulation_files(input_dir)
        self.progress_logger.info(f"Found {len(all_files)} simulation files in {input_dir}")
        
        # Extract rates and filter
        simulations = []
        for file_path in all_files:
            rate = self.file_handler.extract_rate_from_filename(file_path)
            if rate:
                # Apply filters
                if rates and rate not in rates:
                    continue
                if self.config['filters']['rates'] != "all":
                    if rate not in self.config['filters']['rates']:
                        continue
                
                simulations.append((file_path, rate))
                self.progress_logger.info(f"Rate {rate}: {os.path.basename(file_path)}")
        
        self.progress_logger.dedent()
        self.logger.info(f"Selected {len(simulations)} simulations to process")
        
        return simulations
    
    def _validate_environment(self, dry_run: bool) -> bool:
        """Validate environment before processing."""
        self.logger.info("\nPhase 2: Validation")
        self.progress_logger.indent()
        
        # Check disk space
        output_base = self.config['paths']['output_base']
        has_space, available_gb = self.file_handler.check_disk_space(
            output_base,
            self.config['advanced']['min_disk_space']
        )
        
        if not has_space:
            self.progress_logger.error(
                f"Insufficient disk space! Available: {available_gb:.1f} GB, "
                f"Required: {self.config['advanced']['min_disk_space']} GB"
            )
            self.progress_logger.dedent()
            return False
        
        self.progress_logger.info(f"Disk space: {available_gb:.1f} GB available")
        
        # Check required scripts
        required_scripts = [
            'extract_snapshot.py',
            'create_lineages.py',
            'plot_jsd_distribution.py'
        ]
        
        all_found = True
        for script in required_scripts:
            if self.process_manager.check_script_exists(script):
                self.progress_logger.info(f"Script found: {script}")
            else:
                self.progress_logger.error(f"Script missing: {script}")
                all_found = False
        
        if not all_found and not dry_run:
            self.progress_logger.dedent()
            return False
        
        # Check Python version
        python_version = self.process_manager.get_python_version()
        self.progress_logger.info(f"Python version: {python_version}")
        
        self.progress_logger.dedent()
        return True
    
    def _process_simulations(self, 
                           simulations: List[Tuple[str, str]], 
                           dry_run: bool,
                           resume: bool) -> bool:
        """Process all simulations."""
        self.logger.info("\nPhase 3: Processing")
        
        total = len(simulations)
        all_success = True
        
        for idx, (sim_file, rate) in enumerate(simulations, 1):
            # Check if already completed
            if resume and rate in self.state['completed']:
                self.logger.info(f"\nSkipping rate {rate} ({idx}/{total}) - already completed")
                continue
            
            # Check if previously failed and should skip
            if self.config['filters']['skip_failed'] and rate in self.state['failed']:
                self.logger.info(f"\nSkipping rate {rate} ({idx}/{total}) - previously failed")
                continue
            
            self.logger.info(f"\nProcessing rate {rate} ({idx}/{total})")
            
            if dry_run:
                self.logger.info("  [DRY RUN] Would process this simulation")
                continue
            
            # Process this simulation
            success = self._process_single_simulation(sim_file, rate)
            
            if success:
                self.processed_count += 1
            else:
                self.failed_count += 1
                all_success = False
            
            # Save state periodically
            if idx % self.config['advanced']['checkpoint_interval'] == 0:
                self._save_state()
            
            # Show progress
            elapsed = (datetime.now() - self.start_time).total_seconds()
            if self.processed_count > 0:
                avg_time = elapsed / self.processed_count
                remaining = (total - idx) * avg_time
                eta = datetime.now() + timedelta(seconds=remaining)
                self.logger.info(
                    f"Progress: {idx}/{total} | "
                    f"Processed: {self.processed_count} | "
                    f"Failed: {self.failed_count} | "
                    f"ETA: {eta.strftime('%H:%M:%S')}"
                )
        
        # Final state save
        self._save_state()
        
        return all_success
    
    def _process_single_simulation(self, sim_file: str, rate: str) -> bool:
        """
        Process a single simulation file.
        
        Returns:
            True if successful, False otherwise
        """
        self.progress_logger.indent()
        
        # Mark as in progress
        self.state['in_progress'][rate] = {
            'start_time': datetime.now().isoformat(),
            'sim_file': sim_file
        }
        
        try:
            # Create directories
            dirs = self.file_handler.create_rate_directories(
                self.config['paths']['output_base'], 
                rate
            )
            
            # Get output paths
            year = self.config['processing']['year_to_extract']
            paths = self.file_handler.get_output_paths(dirs, year)
            
            # Check what already exists
            existing = self.file_handler.check_outputs_exist(paths)
            
            # Step 1: Extract snapshot
            if self.config['processing']['create_snapshot']:
                if existing['snapshot'] and self.config['filters']['skip_existing']:
                    self.progress_logger.info("Snapshot already exists, skipping")
                else:
                    success, output, duration = self.process_manager.extract_snapshot(
                        sim_file, year, paths['snapshot'],
                        timeout=self.config['processing']['timeouts']['snapshot']
                    )
                    if success:
                        self.progress_logger.success("Extracted snapshot", duration)
                    else:
                        self.progress_logger.error(f"Failed to extract snapshot: {output}")
                        raise Exception("Snapshot extraction failed")
            
            # Step 2: Create lineages
            if self.config['processing']['create_lineages']:
                if (existing['lineages_mutant'] and existing['lineages_control'] 
                    and self.config['filters']['skip_existing']):
                    self.progress_logger.info("Lineages already exist, skipping")
                else:
                    # Create lineages with correct output directory
                    # The script expects a 'lineages' subdirectory
                    lineage_output_dir = os.path.join(dirs['base'], 'lineages')
                    success, output, duration = self.process_manager.create_lineages(
                        paths['snapshot'], 
                        lineage_output_dir,
                        lineage_type=self.config['lineages']['type'],
                        timeout=self.config['processing']['timeouts']['lineages']
                    )
                    if success:
                        self.progress_logger.success("Created lineages", duration)
                    else:
                        self.progress_logger.error(f"Failed to create lineages: {output}")
                        raise Exception("Lineage creation failed")
            
            # Step 3: Create plots
            if self.config['processing']['create_plots']:
                if existing['plot'] and self.config['filters']['skip_existing']:
                    self.progress_logger.info("Plot already exists, skipping")
                else:
                    success, output, duration = self.process_manager.plot_distribution(
                        paths['snapshot'], paths['plot'],
                        bins=self.config['plots']['jsd_bins'],
                        timeout=self.config['processing']['timeouts']['plots']
                    )
                    if success:
                        self.progress_logger.success("Generated plot", duration)
                    else:
                        self.progress_logger.error(f"Failed to create plot: {output}")
                        # Don't fail the whole simulation for plot errors
            
            # Mark as completed
            self.state['completed'][rate] = {
                'timestamp': datetime.now().isoformat(),
                'sim_file': sim_file,
                'outputs': paths
            }
            
            # Remove from in_progress
            if rate in self.state['in_progress']:
                del self.state['in_progress'][rate]
            
            self.progress_logger.dedent()
            self.logger.info(f"Rate {rate} completed successfully")
            return True
            
        except Exception as e:
            # Mark as failed
            self.state['failed'][rate] = {
                'timestamp': datetime.now().isoformat(),
                'error': str(e),
                'attempts': self.state['failed'].get(rate, {}).get('attempts', 0) + 1
            }
            
            # Remove from in_progress
            if rate in self.state['in_progress']:
                del self.state['in_progress'][rate]
            
            self.progress_logger.dedent()
            self.logger.error(f"Rate {rate} failed: {e}")
            return False
    
    
    def _save_state(self):
        """Save current processing state."""
        self.state['last_run'] = datetime.now().isoformat()
        self.file_handler.save_state(self.state_file, self.state)
    
    def _generate_report(self):
        """Generate final summary report."""
        duration = datetime.now() - self.start_time
        
        self.logger.info("\n" + "="*60)
        self.logger.info("BATCH PROCESSING SUMMARY")
        self.logger.info("="*60)
        self.logger.info(f"Start time: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info(f"Duration: {duration}")
        self.logger.info(f"\nProcessed: {self.processed_count}")
        self.logger.info(f"Failed: {self.failed_count}")
        self.logger.info(f"Completed rates: {sorted(self.state['completed'].keys())}")
        
        if self.state['failed']:
            self.logger.info(f"\nFailed rates: {sorted(self.state['failed'].keys())}")
            for rate, info in self.state['failed'].items():
                self.logger.info(f"  - {rate}: {info['error']} (attempts: {info['attempts']})")
        
        self.logger.info(f"\nOutput directory: {self.config['paths']['output_base']}")
        
        # Suggest next steps
        if self.processed_count > 0:
            self.logger.info("\nNext steps:")
            self.logger.info("1. Review the generated outputs in each rate_* directory")
            self.logger.info("2. Check plots to compare JSD distributions across rates")
            self.logger.info("3. Run Step 3 processing for mixed population analysis")
        
        if self.failed_count > 0:
            self.logger.info("\nTo retry failed simulations:")
            self.logger.info("  python batch_processor.py --resume")


def main():
    """Main entry point with command-line interface."""
    parser = argparse.ArgumentParser(
        description="Batch process Step 2 for multiple simulation files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all simulations
  python batch_processor.py
  
  # Process specific rates
  python batch_processor.py --rates 0.001,0.005,0.010
  
  # Dry run to see what would be done
  python batch_processor.py --dry-run
  
  # Resume from previous run
  python batch_processor.py --resume
  
  # Use custom configuration
  python batch_processor.py --config my_config.yaml
        """
    )
    
    parser.add_argument(
        '--config',
        default='config/batch_config.yaml',
        help='Path to configuration file (default: config/batch_config.yaml)'
    )
    
    parser.add_argument(
        '--rates',
        help='Comma-separated list of rates to process (e.g., 0.001,0.005,0.010)'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without actually processing'
    )
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from previous run, skipping completed simulations'
    )
    
    args = parser.parse_args()
    
    # Parse rates if provided
    rates = None
    if args.rates:
        rates = [r.strip() for r in args.rates.split(',')]
    
    # Create and run processor
    try:
        processor = BatchProcessor(args.config)
        success = processor.run(
            rates=rates,
            dry_run=args.dry_run,
            resume=args.resume
        )
        sys.exit(0 if success else 1)
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()