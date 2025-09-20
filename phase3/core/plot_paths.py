"""
Centralized path generation for organized plot output structure.
"""
import os
from typing import Optional

class PlotPaths:
    """Centralized path management for pipeline output organization.

    This class provides a structured approach to organizing pipeline outputs
    into metric-based subdirectories, ensuring consistent file placement and
    easy navigation of results.

    Directory Structure:
        results/
        ├── cell_metrics/          # Cell-level measurements
        │   ├── timeline/          # Time series plots
        │   ├── distributions/     # Histogram distributions
        │   ├── comparisons/       # Batch comparisons
        │   ├── individual_trajectories/  # Per-individual growth
        │   │   ├── jsd/          # Cell JSD trajectories
        │   │   │   ├── mutant/
        │   │   │   └── control1/
        │   │   └── methylation_proportion/  # Methylation trajectories
        │   │       ├── mutant/
        │   │       └── control1/
        │   └── analysis/         # JSON analysis results
        ├── gene_metrics/         # Gene-level measurements
        │   ├── timeline/         # Gene JSD over time
        │   ├── distributions/    # Gene JSD distributions
        │   ├── comparisons/      # Batch comparisons
        │   ├── individual_trajectories/  # Per-individual gene JSD
        │   │   ├── mutant/
        │   │   └── control1/
        │   ├── per_gene/         # Individual gene distributions
        │   └── analysis/         # JSON analysis results
        └── metadata/             # Pipeline configuration and metadata

    Attributes:
        results_dir: Base results directory path
        cell_metrics_dir: Root for cell-level metrics
        gene_metrics_dir: Root for gene-level metrics
        metadata_dir: Root for metadata files
        All subdirectory attributes for organized structure

    Example:
        >>> plot_paths = PlotPaths('/path/to/results')
        >>> plot_paths.create_all_directories()
        >>> jsd_plot = plot_paths.get_cell_jsd_timeline_path()
    """
    
    def __init__(self, results_dir: str):
        self.results_dir = results_dir
        
        # Create main metric directories
        self.cell_metrics_dir = os.path.join(results_dir, 'cell_metrics')
        self.gene_metrics_dir = os.path.join(results_dir, 'gene_metrics')
        self.metadata_dir = os.path.join(results_dir, 'metadata')
        
        # Cell metrics subdirectories
        self.cell_timeline_dir = os.path.join(self.cell_metrics_dir, 'timeline')
        self.cell_distributions_dir = os.path.join(self.cell_metrics_dir, 'distributions')
        self.cell_comparisons_dir = os.path.join(self.cell_metrics_dir, 'comparisons')
        self.cell_trajectories_dir = os.path.join(self.cell_metrics_dir, 'individual_trajectories')
        
        # Cell JSD trajectory subdirectories by batch
        self.cell_jsd_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'jsd')
        self.cell_jsd_mutant_dir = os.path.join(self.cell_jsd_trajectories_dir, 'mutant')
        self.cell_jsd_control1_dir = os.path.join(self.cell_jsd_trajectories_dir, 'control1')
        
        # Cell methylation trajectory subdirectories by batch
        self.cell_meth_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'methylation_proportion')
        self.cell_meth_mutant_dir = os.path.join(self.cell_meth_trajectories_dir, 'mutant')
        self.cell_meth_control1_dir = os.path.join(self.cell_meth_trajectories_dir, 'control1')
        
        self.cell_analysis_dir = os.path.join(self.cell_metrics_dir, 'analysis')
        
        # Gene metrics subdirectories
        self.gene_timeline_dir = os.path.join(self.gene_metrics_dir, 'timeline')
        self.gene_distributions_dir = os.path.join(self.gene_metrics_dir, 'distributions')
        self.gene_comparisons_dir = os.path.join(self.gene_metrics_dir, 'comparisons')
        self.gene_trajectories_dir = os.path.join(self.gene_metrics_dir, 'individual_trajectories')
        
        # Gene trajectory subdirectories by batch
        self.gene_trajectories_mutant_dir = os.path.join(self.gene_trajectories_dir, 'mutant')
        self.gene_trajectories_control1_dir = os.path.join(self.gene_trajectories_dir, 'control1')
        
        self.gene_per_gene_dir = os.path.join(self.gene_metrics_dir, 'per_gene')
        self.gene_analysis_dir = os.path.join(self.gene_metrics_dir, 'analysis')
    
    def create_all_directories(self):
        """Create all necessary subdirectories for organized output.
        
        Creates the complete directory structure for pipeline outputs,
        ensuring all paths are available before file generation begins.
        
        Raises:
            OSError: If directory creation fails due to permissions.
        """
        dirs = [
            # Main directories
            self.cell_metrics_dir, self.gene_metrics_dir, self.metadata_dir,
            # Cell metric subdirectories
            self.cell_timeline_dir, self.cell_distributions_dir, self.cell_comparisons_dir,
            self.cell_trajectories_dir, self.cell_jsd_trajectories_dir, 
            self.cell_jsd_mutant_dir, self.cell_jsd_control1_dir,
            self.cell_meth_trajectories_dir,
            self.cell_meth_mutant_dir, self.cell_meth_control1_dir,
            self.cell_analysis_dir,
            # Gene metric subdirectories
            self.gene_timeline_dir, self.gene_distributions_dir, self.gene_comparisons_dir,
            self.gene_trajectories_dir, 
            self.gene_trajectories_mutant_dir, self.gene_trajectories_control1_dir,
            self.gene_per_gene_dir, self.gene_analysis_dir
        ]
        for dir_path in dirs:
            os.makedirs(dir_path, exist_ok=True)
    
    # Cell metric paths
    def get_cell_jsd_timeline_path(self) -> str:
        return os.path.join(self.cell_timeline_dir, 'cell_jsd_timeline.png')
    
    def get_cell_methylation_timeline_path(self) -> str:
        return os.path.join(self.cell_timeline_dir, 'cell_methylation_proportion_timeline.png')
    
    def get_cell_jsd_distribution_path(self, year: int) -> str:
        return os.path.join(self.cell_distributions_dir, f'year{year}_cell_jsd.png')
    
    def get_cell_methylation_distribution_path(self, year: int) -> str:
        return os.path.join(self.cell_distributions_dir, f'year{year}_cell_methylation_proportion.png')
    
    def get_cell_jsd_comparison_path(self) -> str:
        return os.path.join(self.cell_comparisons_dir, 'cell_jsd_comparison.png')
    
    def get_cell_methylation_comparison_path(self) -> str:
        return os.path.join(self.cell_comparisons_dir, 'cell_methylation_proportion_comparison.png')
    
    def get_individual_cell_jsd_path(self, batch: str, individual_id: int) -> str:
        """Get path for individual cell JSD trajectory plot."""
        if batch == 'mutant':
            return os.path.join(self.cell_jsd_mutant_dir, f'individual_{individual_id:02d}.png')
        elif batch == 'control1':
            return os.path.join(self.cell_jsd_control1_dir, f'individual_{individual_id:02d}.png')
        else:
            raise ValueError(f"Unexpected batch type: {batch}")
    
    def get_individual_cell_methylation_path(self, batch: str, individual_id: int) -> str:
        """Get path for individual cell methylation proportion trajectory plot."""
        if batch == 'mutant':
            return os.path.join(self.cell_meth_mutant_dir, f'individual_{individual_id:02d}.png')
        elif batch == 'control1':
            return os.path.join(self.cell_meth_control1_dir, f'individual_{individual_id:02d}.png')
        else:
            raise ValueError(f"Unexpected batch type: {batch}")
    
    def get_cell_jsd_analysis_path(self) -> str:
        return os.path.join(self.cell_analysis_dir, 'cell_jsd_analysis.json')
    
    def get_cell_methylation_analysis_path(self) -> str:
        return os.path.join(self.cell_analysis_dir, 'cell_methylation_proportion_analysis.json')
    
    # Gene metric paths
    def get_gene_jsd_timeline_path(self) -> str:
        return os.path.join(self.gene_timeline_dir, 'gene_jsd_timeline.png')
    
    def get_gene_jsd_distribution_path(self, year: int) -> str:
        return os.path.join(self.gene_distributions_dir, f'year{year}_gene_jsd.png')
    
    def get_gene_jsd_comparison_path(self) -> str:
        return os.path.join(self.gene_comparisons_dir, 'gene_jsd_comparison.png')
    
    def get_individual_gene_jsd_path(self, batch: str, individual_id: int) -> str:
        """Get path for individual gene JSD trajectory plot."""
        if batch == 'mutant':
            return os.path.join(self.gene_trajectories_mutant_dir, f'individual_{individual_id:02d}.png')
        elif batch == 'control1':
            return os.path.join(self.gene_trajectories_control1_dir, f'individual_{individual_id:02d}.png')
        else:
            raise ValueError(f"Unexpected batch type: {batch}")
    
    def get_per_gene_jsd_path(self, gene_idx: int) -> str:
        return os.path.join(self.gene_per_gene_dir, f'gene_{gene_idx:03d}_jsd.png')
    
    def get_gene_jsd_analysis_path(self) -> str:
        return os.path.join(self.gene_analysis_dir, 'gene_jsd_analysis.json')
    
    def get_gene_methylation_analysis_path(self) -> str:
        return os.path.join(self.gene_analysis_dir, 'gene_methylation_analysis.json')
    
    def get_gene_methylation_comparison_path(self) -> str:
        return os.path.join(self.gene_comparisons_dir, 'gene_methylation_comparison.png')
    
    def get_gene_methylation_distribution_path(self, year: int) -> str:
        return os.path.join(self.gene_distributions_dir, f'year{year}_gene_methylation.png')
    
    def get_gene_methylation_timeline_path(self) -> str:
        return os.path.join(self.gene_timeline_dir, 'gene_methylation_timeline.png')
    
    def get_individual_gene_methylation_path(self, batch: str, individual_id: int) -> str:
        """Get path for individual gene methylation trajectory plot."""
        if batch == 'mutant':
            return os.path.join(self.gene_trajectories_mutant_dir, f'individual_{individual_id:02d}_methylation.png')
        elif batch == 'control1':
            return os.path.join(self.gene_trajectories_control1_dir, f'individual_{individual_id:02d}_methylation.png')
        else:
            return None  # Control2 doesn't have trajectories
    
    # Metadata paths
    def get_pipeline_metadata_path(self) -> str:
        return os.path.join(self.metadata_dir, 'pipeline_metadata.json')
    
    def get_mixing_statistics_path(self) -> str:
        return os.path.join(self.metadata_dir, 'mixing_statistics.json')
    
    def validate_structure(self) -> bool:
        """Validate that all expected directories exist.
        
        Returns:
            bool: True if all directories exist, False otherwise.
        
        Example:
            >>> if not plot_paths.validate_structure():
            ...     plot_paths.create_all_directories()
        """
        dirs = [
            self.cell_metrics_dir, self.gene_metrics_dir, self.metadata_dir,
            self.cell_timeline_dir, self.cell_distributions_dir, self.cell_comparisons_dir,
            self.cell_trajectories_dir, self.cell_jsd_trajectories_dir, 
            self.cell_jsd_mutant_dir, self.cell_jsd_control1_dir,
            self.cell_meth_trajectories_dir,
            self.cell_meth_mutant_dir, self.cell_meth_control1_dir,
            self.cell_analysis_dir,
            self.gene_timeline_dir, self.gene_distributions_dir, self.gene_comparisons_dir,
            self.gene_trajectories_dir,
            self.gene_trajectories_mutant_dir, self.gene_trajectories_control1_dir,
            self.gene_per_gene_dir, self.gene_analysis_dir
        ]
        return all(os.path.exists(d) for d in dirs)
    
    def get_all_paths(self) -> dict:
        """Return dictionary of all configured paths.
        
        Returns:
            Dict mapping path names to their full paths.
            
        Example:
            >>> paths = plot_paths.get_all_paths()
            >>> print(f"Cell metrics at: {paths['cell_metrics_dir']}")
        """
        return {
            'results_dir': self.results_dir,
            'cell_metrics_dir': self.cell_metrics_dir,
            'gene_metrics_dir': self.gene_metrics_dir,
            'metadata_dir': self.metadata_dir,
            'cell_timeline_dir': self.cell_timeline_dir,
            'cell_distributions_dir': self.cell_distributions_dir,
            'cell_comparisons_dir': self.cell_comparisons_dir,
            'cell_trajectories_dir': self.cell_trajectories_dir,
            'cell_jsd_trajectories_dir': self.cell_jsd_trajectories_dir,
            'cell_jsd_mutant_dir': self.cell_jsd_mutant_dir,
            'cell_jsd_control1_dir': self.cell_jsd_control1_dir,
            'cell_meth_trajectories_dir': self.cell_meth_trajectories_dir,
            'cell_meth_mutant_dir': self.cell_meth_mutant_dir,
            'cell_meth_control1_dir': self.cell_meth_control1_dir,
            'cell_analysis_dir': self.cell_analysis_dir,
            'gene_timeline_dir': self.gene_timeline_dir,
            'gene_distributions_dir': self.gene_distributions_dir,
            'gene_comparisons_dir': self.gene_comparisons_dir,
            'gene_trajectories_dir': self.gene_trajectories_dir,
            'gene_trajectories_mutant_dir': self.gene_trajectories_mutant_dir,
            'gene_trajectories_control1_dir': self.gene_trajectories_control1_dir,
            'gene_per_gene_dir': self.gene_per_gene_dir,
            'gene_analysis_dir': self.gene_analysis_dir
        }