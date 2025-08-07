#!/usr/bin/env python3
"""
Batch processor for Step 3 of the methylation simulation pipeline.

This script processes outputs from Step 2, creating mixed populations
and analyzing the effects of lineage mixing on JSD distributions.
"""

import os
import sys
import yaml
import json
import argparse
import numpy as np
from datetime import datetime, timedelta
from typing import Dict, List, Optional, Tuple
from collections import defaultdict

from utils import setup_logger, get_logger, FileHandler
from utils.logger import ProgressLogger
from utils.step3_file_handlers import Step3FileHandler
from utils.step3_process_manager import Step3ProcessManager


class Step3BatchProcessor:
    """Orchestrates batch processing for Step 3."""
    
    def __init__(self, config_path: str = "config/step3_config.yaml"):
        """Initialize processor with configuration."""
        self.config = self._load_config(config_path)
        self.logger = setup_logger(
            name="step3_batch",
            log_dir=self.config['logging']['log_dir'],
            level=self.config['logging']['level'],
            console=self.config['logging']['console']
        )
        self.progress_logger = ProgressLogger(self.logger)
        self.file_handler = FileHandler()
        self.step3_handler = Step3FileHandler()
        self.process_manager = Step3ProcessManager(script_dir=os.getcwd())
        
        # Load or initialize state
        self.state_file = self.config['paths']['state_file']
        self.state = self.file_handler.load_state(self.state_file)
        
        # Track current run
        self.start_time = None
        self.processed_count = 0
        self.failed_count = 0
        self.all_results = {}  # Store results for cross-rate analysis
    
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
        Run batch processing for Step 3.
        
        Args:
            rates: Specific rates to process (None for all)
            dry_run: If True, only show what would be done
            resume: If True, skip completed rates
            
        Returns:
            True if all rates processed successfully
        """
        self.start_time = datetime.now()
        self.logger.info("="*60)
        self.logger.info("STEP 3 BATCH PROCESSING STARTED")
        self.logger.info(f"Start time: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("="*60)
        
        try:
            # Phase 1: Discovery
            available_rates = self._discover_rates(rates)
            if not available_rates:
                self.logger.warning("No completed Step 2 rates found!")
                return False
            
            # Phase 2: Validation
            if not self._validate_environment(dry_run):
                return False
            
            # Phase 3: Processing
            success = self._process_rates(available_rates, dry_run, resume)
            
            # Phase 4: Cross-rate analysis
            if success and self.config['processing']['cross_rate_analysis'] and not dry_run:
                self._perform_cross_rate_analysis()
            
            # Phase 5: Reporting
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
    
    def _discover_rates(self, rates: Optional[List[str]] = None) -> List[str]:
        """Discover available rates from Step 2."""
        self.logger.info("Phase 1: Discovery")
        self.progress_logger.indent()
        
        # Find completed Step 2 rates
        step2_dir = self.config['paths']['step2_data']
        completed_rates = self.step3_handler.find_completed_rates(step2_dir)
        self.progress_logger.info(f"Found {len(completed_rates)} completed Step 2 rates")
        
        # Filter based on user input
        if rates:
            available_rates = [r for r in rates if r in completed_rates]
            missing = [r for r in rates if r not in completed_rates]
            if missing:
                self.progress_logger.warning(f"Rates not found in Step 2: {missing}")
        else:
            available_rates = completed_rates
        
        # Apply config filters
        if self.config['filters']['rates'] != "all":
            available_rates = [r for r in available_rates if r in self.config['filters']['rates']]
        
        for rate in available_rates:
            # Check if original simulation exists
            sim_path = self.step3_handler.get_original_simulation_path(
                self.config['paths']['step1_data'], rate
            )
            if sim_path:
                self.progress_logger.info(f"Rate {rate}: Found simulation and lineages")
            else:
                self.progress_logger.warning(f"Rate {rate}: Original simulation not found!")
                available_rates.remove(rate)
        
        self.progress_logger.dedent()
        self.logger.info(f"Selected {len(available_rates)} rates to process: {available_rates}")
        
        return available_rates
    
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
        all_found = all(self.process_manager.check_script_dependencies().values())
        
        if not all_found and not dry_run:
            self.progress_logger.dedent()
            return False
        
        self.progress_logger.dedent()
        return True
    
    def _process_rates(self, rates: List[str], dry_run: bool, resume: bool) -> bool:
        """Process all rates through Step 3."""
        self.logger.info("\nPhase 3: Processing")
        
        total = len(rates)
        all_success = True
        
        for idx, rate in enumerate(rates, 1):
            # Check if already completed
            if resume and rate in self.state['completed']:
                self.logger.info(f"\nSkipping rate {rate} ({idx}/{total}) - already completed")
                # Load previous results for cross-rate analysis
                results_file = os.path.join(
                    self.config['paths']['output_base'],
                    f"rate_{rate}",
                    "results",
                    "jsd_distributions.json"
                )
                if os.path.exists(results_file):
                    with open(results_file, 'r') as f:
                        self.all_results[rate] = json.load(f)
                continue
            
            self.logger.info(f"\nProcessing rate {rate} ({idx}/{total})")
            
            if dry_run:
                self.logger.info("  [DRY RUN] Would process this rate")
                continue
            
            # Process this rate
            success, results = self._process_single_rate(rate)
            
            if success:
                self.processed_count += 1
                self.all_results[rate] = results
            else:
                self.failed_count += 1
                all_success = False
            
            # Save state periodically
            if idx % self.config['advanced']['checkpoint_interval'] == 0:
                self._save_state()
            
            # Show progress
            self._show_progress(idx, total)
        
        # Final state save
        self._save_state()
        
        return all_success
    
    def _process_single_rate(self, rate: str) -> Tuple[bool, Optional[Dict]]:
        """Process a single rate through Step 3."""
        self.progress_logger.indent()
        
        # Mark as in progress
        self.state['in_progress'][rate] = {
            'start_time': datetime.now().isoformat()
        }
        
        try:
            # Get paths
            step1_dir = self.config['paths']['step1_data']
            step2_dir = self.config['paths']['step2_data']
            
            # Find original simulation
            sim_path = self.step3_handler.get_original_simulation_path(step1_dir, rate)
            if not sim_path:
                raise Exception(f"Original simulation not found for rate {rate}")
            
            # Get Step 2 paths
            step2_paths = self.step3_handler.get_step2_paths(step2_dir, rate)
            
            # Create Step 3 directories
            dirs = self.step3_handler.create_step3_directories(
                self.config['paths']['output_base'], rate
            )
            
            # Check existing outputs
            existing = self.step3_handler.check_step3_outputs_exist(dirs)
            
            # Step 1: Extract year 60 snapshot
            if self.config['processing']['extract_snapshots']:
                if existing['snapshot'] and self.config['filters']['skip_existing']:
                    self.progress_logger.info("Year 60 snapshot already exists, skipping")
                else:
                    snapshot_path = os.path.join(
                        dirs['snapshots'], 
                        f"year{self.config['processing']['year_to_extract']}_snapshot.json.gz"
                    )
                    success, output, duration = self.process_manager.extract_year60(
                        sim_path, snapshot_path,
                        timeout=self.config['processing']['timeouts']['snapshot']
                    )
                    if success:
                        self.progress_logger.success("Extracted year 60 snapshot", duration)
                    else:
                        raise Exception(f"Failed to extract snapshot: {output}")
            
            # Step 2: Create individuals
            if self.config['processing']['create_individuals']:
                if (existing['individuals_mutant'] and existing['individuals_control'] 
                    and self.config['filters']['skip_existing']):
                    self.progress_logger.info("Individuals already exist, skipping")
                else:
                    snapshot_path = os.path.join(
                        dirs['snapshots'],
                        f"year{self.config['processing']['year_to_extract']}_snapshot.json.gz"
                    )
                    
                    # The create_individuals script expects specific directory structure
                    lineage_base = os.path.dirname(step2_paths['lineages_mutant'])
                    
                    success, output, duration = self.process_manager.create_individuals(
                        snapshot_path,
                        lineage_base,
                        dirs['base'],
                        individual_type="both",
                        num_individuals=self.config['processing']['num_individuals'],
                        mixture_ratio=self.config['processing']['mixture_ratio'],
                        timeout=self.config['processing']['timeouts']['individuals']
                    )
                    if success:
                        self.progress_logger.success("Created individuals", duration)
                    else:
                        raise Exception(f"Failed to create individuals: {output}")
            
            # Step 3: Create plots and analysis
            results = None
            if self.config['processing']['create_plots']:
                if existing['plot'] and self.config['filters']['skip_existing']:
                    self.progress_logger.info("Plots already exist, skipping")
                    # Load existing results
                    results_path = os.path.join(dirs['results'], 'jsd_distributions.json')
                    if os.path.exists(results_path):
                        with open(results_path, 'r') as f:
                            results = json.load(f)
                else:
                    success, output, duration = self.process_manager.plot_distributions(
                        dirs['individuals_mutant'],
                        dirs['individuals_control'],
                        dirs['base'],
                        timeout=self.config['processing']['timeouts']['plots']
                    )
                    if success:
                        self.progress_logger.success("Generated plots and analysis", duration)
                        # Load the results
                        results_path = os.path.join(dirs['results'], 'jsd_distributions.json')
                        if os.path.exists(results_path):
                            with open(results_path, 'r') as f:
                                results = json.load(f)
                    else:
                        self.progress_logger.warning(f"Failed to create plots: {output}")
            
            # Mark as completed
            self.state['completed'][rate] = {
                'timestamp': datetime.now().isoformat(),
                'outputs': dirs
            }
            
            # Remove from in_progress
            if rate in self.state['in_progress']:
                del self.state['in_progress'][rate]
            
            self.progress_logger.dedent()
            self.logger.info(f"Rate {rate} completed successfully")
            return True, results
            
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
            return False, None
    
    def _perform_cross_rate_analysis(self):
        """Perform analysis across all processed rates."""
        self.logger.info("\nPhase 4: Cross-Rate Analysis")
        self.progress_logger.indent()
        
        if len(self.all_results) < 2:
            self.progress_logger.warning("Need at least 2 rates for cross-rate analysis")
            self.progress_logger.dedent()
            return
        
        try:
            # Create combined analysis directories
            combined_dirs = self.step3_handler.create_combined_analysis_dirs(
                self.config['paths']['output_base']
            )
            
            # Prepare data for analysis
            rates_sorted = sorted(self.all_results.keys())
            mutant_means = []
            control_means = []
            
            for rate in rates_sorted:
                if rate in self.all_results and self.all_results[rate]:
                    mutant_means.append(self.all_results[rate].get('mutant_mean', []))
                    control_means.append(self.all_results[rate].get('control_mean', []))
            
            # Generate combined plot
            self._create_combined_plot(rates_sorted, mutant_means, control_means, combined_dirs)
            
            # Perform statistical analysis
            summary = self._create_summary_statistics(rates_sorted, self.all_results)
            
            # Save summary
            summary_path = os.path.join(combined_dirs['results'], 'summary_analysis.json')
            with open(summary_path, 'w') as f:
                json.dump(summary, f, indent=2)
            
            self.progress_logger.success("Completed cross-rate analysis")
            
        except Exception as e:
            self.progress_logger.error(f"Cross-rate analysis failed: {e}")
        
        self.progress_logger.dedent()
    
    def _create_combined_plot(self, rates: List[str], mutant_data: List, control_data: List, dirs: Dict):
        """Create combined visualization across rates."""
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
            
            # Create subplots
            fig = make_subplots(
                rows=1, cols=2,
                subplot_titles=("Mean JSD by Rate", "Mutant vs Control Difference"),
                horizontal_spacing=0.15
            )
            
            # Convert rates to percentages for display
            rate_labels = [f"{float(r)*100:.1f}%" for r in rates]
            
            # Calculate means and differences
            mutant_means_avg = [np.mean(m) if m else 0 for m in mutant_data]
            control_means_avg = [np.mean(c) if c else 0 for c in control_data]
            differences = [m - c for m, c in zip(mutant_means_avg, control_means_avg)]
            
            # Plot 1: Mean JSD by rate
            fig.add_trace(
                go.Scatter(
                    x=rate_labels, y=mutant_means_avg,
                    mode='lines+markers',
                    name='Mutant',
                    line=dict(color='red', width=2),
                    marker=dict(size=8)
                ),
                row=1, col=1
            )
            
            fig.add_trace(
                go.Scatter(
                    x=rate_labels, y=control_means_avg,
                    mode='lines+markers',
                    name='Control',
                    line=dict(color='blue', width=2),
                    marker=dict(size=8)
                ),
                row=1, col=1
            )
            
            # Plot 2: Difference
            fig.add_trace(
                go.Bar(
                    x=rate_labels, y=differences,
                    name='Difference',
                    marker_color=['red' if d > 0 else 'blue' for d in differences]
                ),
                row=1, col=2
            )
            
            # Update layout
            fig.update_xaxes(title_text="Methylation Rate", row=1, col=1)
            fig.update_xaxes(title_text="Methylation Rate", row=1, col=2)
            fig.update_yaxes(title_text="Mean JSD", row=1, col=1)
            fig.update_yaxes(title_text="Difference (Mutant - Control)", row=1, col=2)
            
            fig.update_layout(
                title="Cross-Rate Analysis: Effect of Lineage Mixing",
                height=500,
                width=1200,
                showlegend=True
            )
            
            # Save plots
            output_path = os.path.join(dirs['plots'], 'cross_rate_analysis')
            fig.write_html(f"{output_path}.html")
            fig.write_image(f"{output_path}.png")
            
            self.progress_logger.info("Created combined analysis plots")
            
        except Exception as e:
            self.progress_logger.error(f"Failed to create combined plot: {e}")
    
    def _create_summary_statistics(self, rates: List[str], results: Dict) -> Dict:
        """Create summary statistics across all rates."""
        summary = {
            'analysis_date': datetime.now().isoformat(),
            'rates_analyzed': rates,
            'rate_summaries': {}
        }
        
        for rate in rates:
            if rate in results and results[rate]:
                rate_data = results[rate]
                summary['rate_summaries'][rate] = {
                    'mutant_mean': rate_data.get('mutant_mean_of_means', 0),
                    'control_mean': rate_data.get('control_mean_of_means', 0),
                    'difference': rate_data.get('mutant_mean_of_means', 0) - rate_data.get('control_mean_of_means', 0),
                    'p_value': rate_data.get('p_value', None),
                    'statistic': rate_data.get('statistic', None)
                }
        
        # Overall trend analysis
        differences = [s['difference'] for s in summary['rate_summaries'].values()]
        summary['overall_trend'] = {
            'mean_difference': np.mean(differences),
            'std_difference': np.std(differences),
            'trend': 'increasing' if differences[-1] > differences[0] else 'decreasing'
        }
        
        return summary
    
    def _show_progress(self, current: int, total: int):
        """Show progress with ETA."""
        elapsed = (datetime.now() - self.start_time).total_seconds()
        if self.processed_count > 0:
            avg_time = elapsed / self.processed_count
            remaining = (total - current) * avg_time
            eta = datetime.now() + timedelta(seconds=remaining)
            self.logger.info(
                f"Progress: {current}/{total} | "
                f"Processed: {self.processed_count} | "
                f"Failed: {self.failed_count} | "
                f"ETA: {eta.strftime('%H:%M:%S')}"
            )
    
    def _save_state(self):
        """Save current processing state."""
        self.state['last_run'] = datetime.now().isoformat()
        self.file_handler.save_state(self.state_file, self.state)
    
    def _generate_report(self):
        """Generate final summary report."""
        duration = datetime.now() - self.start_time
        
        self.logger.info("\n" + "="*60)
        self.logger.info("STEP 3 BATCH PROCESSING SUMMARY")
        self.logger.info("="*60)
        self.logger.info(f"Start time: {self.start_time.strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info(f"End time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info(f"Duration: {duration}")
        self.logger.info(f"\nProcessed: {self.processed_count}")
        self.logger.info(f"Failed: {self.failed_count}")
        self.logger.info(f"Completed rates: {sorted(self.state['completed'].keys())}")
        
        if self.state['failed']:
            self.logger.info(f"\nFailed rates: {sorted(self.state['failed'].keys())}")
        
        self.logger.info(f"\nOutput directory: {self.config['paths']['output_base']}")
        
        if self.config['processing']['cross_rate_analysis'] and len(self.all_results) > 1:
            self.logger.info("\nCross-rate analysis: COMPLETED")
            self.logger.info("  See combined_analysis/ for results")
        
        # Next steps
        if self.processed_count > 0:
            self.logger.info("\nNext steps:")
            self.logger.info("1. Review individual rate results in each rate_* directory")
            self.logger.info("2. Check combined_analysis/ for cross-rate comparisons")
            self.logger.info("3. Examine statistical results in results/ subdirectories")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Batch process Step 3 for multiple rates",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Process all rates from Step 2
  python batch_processor_step3.py
  
  # Process specific rates
  python batch_processor_step3.py --rates 0.003000,0.005000
  
  # Dry run
  python batch_processor_step3.py --dry-run
  
  # Resume from previous run
  python batch_processor_step3.py --resume
        """
    )
    
    parser.add_argument(
        '--config',
        default='config/step3_config.yaml',
        help='Configuration file path'
    )
    
    parser.add_argument(
        '--rates',
        help='Comma-separated list of rates to process'
    )
    
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='Show what would be done without processing'
    )
    
    parser.add_argument(
        '--resume',
        action='store_true',
        help='Resume from previous run'
    )
    
    args = parser.parse_args()
    
    # Parse rates
    rates = None
    if args.rates:
        rates = [r.strip() for r in args.rates.split(',')]
    
    # Run processor
    try:
        processor = Step3BatchProcessor(args.config)
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