"""
Core modules for phase3 analysis.
"""

from .data_loader import (
    load_snapshot_cells,
    load_petri_dish,
    load_all_petri_dishes,
    extract_gene_jsd_from_history,
    load_phase2_metadata,
    smart_open
)

from .analysis_functions import (
    plot_cell_jsd_distribution,
    plot_cell_methylation_proportion_histogram,
    plot_gene_jsd_distribution,
    plot_gene_methylation_proportion_histogram,
    analyze_cell_methylation_proportion_comparison,
    analyze_populations_from_dishes,
    plot_gene_jsd_distributions,
    plot_gene_jsd_individual_comparison,
    generate_gene_methylation_analysis,
    plot_gene_methylation_individual_comparison,
    generate_gene_jsd_analysis,
    get_cell_jsd_array
)

from .plot_paths import PlotPaths

__all__ = [
    # Data loading
    'load_snapshot_cells',
    'load_petri_dish',
    'load_all_petri_dishes',
    'extract_gene_jsd_from_history',
    'load_phase2_metadata',
    'smart_open',
    # Analysis functions
    'plot_cell_jsd_distribution',
    'plot_cell_methylation_proportion_histogram',
    'plot_gene_jsd_distribution',
    'plot_gene_methylation_proportion_histogram',
    'analyze_cell_methylation_proportion_comparison',
    'analyze_populations_from_dishes',
    'plot_gene_jsd_distributions',
    'plot_gene_jsd_individual_comparison',
    'generate_gene_methylation_analysis',
    'plot_gene_methylation_individual_comparison',
    'generate_gene_jsd_analysis',
    'get_cell_jsd_array',
    # Plot paths
    'PlotPaths'
]