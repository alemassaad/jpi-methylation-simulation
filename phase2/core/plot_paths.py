"""
Centralized path generation for organized plot output structure.
"""
import os
from typing import Optional

class PlotPaths:
    """Generate organized paths for plot outputs."""
    
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
        self.cell_jsd_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'jsd')
        self.cell_meth_trajectories_dir = os.path.join(self.cell_trajectories_dir, 'methylation_proportion')
        self.cell_analysis_dir = os.path.join(self.cell_metrics_dir, 'analysis')
        
        # Gene metrics subdirectories
        self.gene_timeline_dir = os.path.join(self.gene_metrics_dir, 'timeline')
        self.gene_distributions_dir = os.path.join(self.gene_metrics_dir, 'distributions')
        self.gene_comparisons_dir = os.path.join(self.gene_metrics_dir, 'comparisons')
        self.gene_trajectories_dir = os.path.join(self.gene_metrics_dir, 'individual_trajectories')
        self.gene_per_gene_dir = os.path.join(self.gene_metrics_dir, 'per_gene')
        self.gene_analysis_dir = os.path.join(self.gene_metrics_dir, 'analysis')
    
    def create_all_directories(self):
        """Create all necessary directories."""
        dirs = [
            self.cell_metrics_dir, self.gene_metrics_dir, self.metadata_dir,
            self.cell_timeline_dir, self.cell_distributions_dir, self.cell_comparisons_dir,
            self.cell_trajectories_dir, self.cell_jsd_trajectories_dir, self.cell_meth_trajectories_dir,
            self.cell_analysis_dir,
            self.gene_timeline_dir, self.gene_distributions_dir, self.gene_comparisons_dir,
            self.gene_trajectories_dir, self.gene_per_gene_dir, self.gene_analysis_dir
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
        return os.path.join(self.cell_jsd_trajectories_dir, f'{batch}_{individual_id:02d}_cell_jsd.png')
    
    def get_individual_cell_methylation_path(self, batch: str, individual_id: int) -> str:
        return os.path.join(self.cell_meth_trajectories_dir, f'{batch}_{individual_id:02d}_cell_methylation_proportion.png')
    
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
        return os.path.join(self.gene_trajectories_dir, f'{batch}_{individual_id:02d}_gene_jsd.png')
    
    def get_per_gene_jsd_path(self, gene_idx: int) -> str:
        return os.path.join(self.gene_per_gene_dir, f'gene_{gene_idx:03d}_jsd.png')
    
    def get_gene_jsd_analysis_path(self) -> str:
        return os.path.join(self.gene_analysis_dir, 'gene_jsd_analysis.json')
    
    # Metadata paths
    def get_pipeline_metadata_path(self) -> str:
        return os.path.join(self.metadata_dir, 'pipeline_metadata.json')
    
    def get_mixing_statistics_path(self) -> str:
        return os.path.join(self.metadata_dir, 'mixing_statistics.json')