"""
Core pipeline modules for Phase 2.
"""

from .pipeline_utils import (
    load_snapshot_as_cells,
    save_snapshot_cells,
    load_snapshot_cells,
    sample_by_quantiles,
    sample_uniform,
    mix_petri_with_snapshot,
    mix_petri_uniform,
    create_uniform_mixing_pool,
    normalize_populations,
    normalize_individuals_for_uniform_mixing,
    create_pure_snapshot_petri,
    create_control2_with_uniform_base,
    load_petri_dish,
    save_petri_dish,
    load_all_petri_dishes,
    check_petri_files_state,
    grow_petri_for_years,
    calculate_population_statistics,
    print_mixing_statistics,
    dict_to_cell,
    get_petri_statistics,
    get_cell_jsd_array
)

from .pipeline_analysis import (
    plot_cell_jsd_distribution,
    analyze_populations_from_dishes,
    generate_gene_jsd_analysis,
    plot_gene_jsd_distributions,
    plot_gene_jsd_individual_comparison
)

from .path_utils import (
    parse_step1_simulation_path,
    generate_step23_output_dir
)

__all__ = [
    # pipeline_utils
    'load_snapshot_as_cells',
    'save_snapshot_cells',
    'load_snapshot_cells',
    'sample_by_quantiles',
    'sample_uniform',
    'mix_petri_with_snapshot',
    'mix_petri_uniform',
    'create_uniform_mixing_pool',
    'normalize_populations',
    'normalize_individuals_for_uniform_mixing',
    'create_pure_snapshot_petri',
    'create_control2_with_uniform_base',
    'load_petri_dish',
    'save_petri_dish',
    'load_all_petri_dishes',
    'check_petri_files_state',
    'grow_petri_for_years',
    'calculate_population_statistics',
    'print_mixing_statistics',
    'dict_to_cell',
    'get_petri_statistics',
    'get_cell_jsd_array',
    # pipeline_analysis
    'plot_cell_jsd_distribution',
    'analyze_populations_from_dishes',
    'generate_gene_jsd_analysis',
    'plot_gene_jsd_distributions',
    'plot_gene_jsd_individual_comparison',
    # path_utils
    'parse_step1_simulation_path',
    'generate_step23_output_dir'
]