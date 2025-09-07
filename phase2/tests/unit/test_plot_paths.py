#!/usr/bin/env python3
"""Test PlotPaths class functionality."""

import os
import sys
import tempfile
import shutil

# Add parent directories to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

# Add phase1 to path for cell module
phase1_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))), 'phase1')
sys.path.append(phase1_path)

from core.plot_paths import PlotPaths

def test_directory_creation():
    """Test that all directories are created correctly."""
    print("Testing directory creation...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        results_dir = os.path.join(tmpdir, "results")
        plot_paths = PlotPaths(results_dir)
        plot_paths.create_all_directories()
        
        # Check main directories
        assert os.path.exists(plot_paths.cell_metrics_dir), "cell_metrics_dir not created"
        assert os.path.exists(plot_paths.gene_metrics_dir), "gene_metrics_dir not created"
        assert os.path.exists(plot_paths.metadata_dir), "metadata_dir not created"
        print("  ✓ Main directories created")
        
        # Check cell metric subdirectories
        assert os.path.exists(plot_paths.cell_timeline_dir), "cell_timeline_dir not created"
        assert os.path.exists(plot_paths.cell_distributions_dir), "cell_distributions_dir not created"
        assert os.path.exists(plot_paths.cell_comparisons_dir), "cell_comparisons_dir not created"
        assert os.path.exists(plot_paths.cell_trajectories_dir), "cell_trajectories_dir not created"
        assert os.path.exists(plot_paths.cell_analysis_dir), "cell_analysis_dir not created"
        print("  ✓ Cell metric subdirectories created")
        
        # Check cell trajectory batch subdirectories
        assert os.path.exists(plot_paths.cell_jsd_mutant_dir), "cell_jsd_mutant_dir not created"
        assert os.path.exists(plot_paths.cell_jsd_control1_dir), "cell_jsd_control1_dir not created"
        assert os.path.exists(plot_paths.cell_meth_mutant_dir), "cell_meth_mutant_dir not created"
        assert os.path.exists(plot_paths.cell_meth_control1_dir), "cell_meth_control1_dir not created"
        print("  ✓ Cell trajectory batch subdirectories created")
        
        # Check gene metric subdirectories
        assert os.path.exists(plot_paths.gene_timeline_dir), "gene_timeline_dir not created"
        assert os.path.exists(plot_paths.gene_distributions_dir), "gene_distributions_dir not created"
        assert os.path.exists(plot_paths.gene_comparisons_dir), "gene_comparisons_dir not created"
        assert os.path.exists(plot_paths.gene_trajectories_dir), "gene_trajectories_dir not created"
        assert os.path.exists(plot_paths.gene_per_gene_dir), "gene_per_gene_dir not created"
        assert os.path.exists(plot_paths.gene_analysis_dir), "gene_analysis_dir not created"
        print("  ✓ Gene metric subdirectories created")
        
        # Check gene trajectory batch subdirectories
        assert os.path.exists(plot_paths.gene_trajectories_mutant_dir), "gene_trajectories_mutant_dir not created"
        assert os.path.exists(plot_paths.gene_trajectories_control1_dir), "gene_trajectories_control1_dir not created"
        print("  ✓ Gene trajectory batch subdirectories created")
        
        print("  ✓ All directories created successfully")

def test_path_generation():
    """Test that paths are generated correctly."""
    print("Testing path generation...")
    
    plot_paths = PlotPaths("/test/results")
    
    # Test cell metric paths
    jsd_path = plot_paths.get_cell_jsd_timeline_path()
    assert jsd_path == "/test/results/cell_metrics/timeline/cell_jsd_timeline.png"
    print("  ✓ Cell JSD timeline path correct")
    
    meth_path = plot_paths.get_cell_methylation_timeline_path()
    assert meth_path == "/test/results/cell_metrics/timeline/cell_methylation_proportion_timeline.png"
    print("  ✓ Cell methylation timeline path correct")
    
    # Test distribution paths
    dist_path = plot_paths.get_cell_jsd_distribution_path(30)
    assert dist_path == "/test/results/cell_metrics/distributions/year30_cell_jsd.png"
    print("  ✓ Cell JSD distribution path correct")
    
    # Test individual paths with batch subdirectories
    mutant_jsd_path = plot_paths.get_individual_cell_jsd_path("mutant", 5)
    assert mutant_jsd_path == "/test/results/cell_metrics/individual_trajectories/jsd/mutant/individual_05.png"
    print("  ✓ Mutant cell JSD trajectory path correct")
    
    control1_jsd_path = plot_paths.get_individual_cell_jsd_path("control1", 3)
    assert control1_jsd_path == "/test/results/cell_metrics/individual_trajectories/jsd/control1/individual_03.png"
    print("  ✓ Control1 cell JSD trajectory path correct")
    
    mutant_meth_path = plot_paths.get_individual_cell_methylation_path("mutant", 2)
    assert mutant_meth_path == "/test/results/cell_metrics/individual_trajectories/methylation_proportion/mutant/individual_02.png"
    print("  ✓ Mutant cell methylation trajectory path correct")
    
    # Test gene paths
    gene_jsd_path = plot_paths.get_gene_jsd_timeline_path()
    assert gene_jsd_path == "/test/results/gene_metrics/timeline/gene_jsd_timeline.png"
    print("  ✓ Gene JSD timeline path correct")
    
    mutant_gene_path = plot_paths.get_individual_gene_jsd_path("mutant", 7)
    assert mutant_gene_path == "/test/results/gene_metrics/individual_trajectories/mutant/individual_07.png"
    print("  ✓ Mutant gene JSD trajectory path correct")
    
    control1_gene_path = plot_paths.get_individual_gene_jsd_path("control1", 4)
    assert control1_gene_path == "/test/results/gene_metrics/individual_trajectories/control1/individual_04.png"
    print("  ✓ Control1 gene JSD trajectory path correct")
    
    # Test per-gene path
    per_gene_path = plot_paths.get_per_gene_jsd_path(12)
    assert per_gene_path == "/test/results/gene_metrics/per_gene/gene_012_jsd.png"
    print("  ✓ Per-gene JSD path correct")
    
    # Test analysis paths
    cell_analysis_path = plot_paths.get_cell_jsd_analysis_path()
    assert cell_analysis_path == "/test/results/cell_metrics/analysis/cell_jsd_analysis.json"
    print("  ✓ Cell JSD analysis path correct")
    
    # Test metadata paths
    metadata_path = plot_paths.get_pipeline_metadata_path()
    assert metadata_path == "/test/results/metadata/pipeline_metadata.json"
    print("  ✓ Pipeline metadata path correct")
    
    print("  ✓ Path generation works correctly")

def test_validation():
    """Test directory structure validation."""
    print("Testing structure validation...")
    
    with tempfile.TemporaryDirectory() as tmpdir:
        plot_paths = PlotPaths(tmpdir)
        
        # Should be invalid before creation
        assert not plot_paths.validate_structure(), "Structure incorrectly validated before creation"
        print("  ✓ Validation correctly returns False before directory creation")
        
        # Create and validate
        plot_paths.create_all_directories()
        assert plot_paths.validate_structure(), "Structure not validated after creation"
        print("  ✓ Validation correctly returns True after directory creation")
        
        # Remove a directory and check validation fails
        os.rmdir(plot_paths.cell_jsd_mutant_dir)
        assert not plot_paths.validate_structure(), "Structure incorrectly validated with missing directory"
        print("  ✓ Validation correctly detects missing directory")
        
        print("  ✓ Validation works correctly")

def test_get_all_paths():
    """Test get_all_paths method."""
    print("Testing get_all_paths method...")
    
    plot_paths = PlotPaths("/test/results")
    all_paths = plot_paths.get_all_paths()
    
    # Check that dictionary contains expected keys
    expected_keys = [
        'results_dir', 'cell_metrics_dir', 'gene_metrics_dir', 'metadata_dir',
        'cell_timeline_dir', 'cell_distributions_dir', 'cell_comparisons_dir',
        'cell_trajectories_dir', 'cell_jsd_trajectories_dir', 'cell_jsd_mutant_dir',
        'cell_jsd_control1_dir', 'cell_meth_trajectories_dir', 'cell_meth_mutant_dir',
        'cell_meth_control1_dir', 'cell_analysis_dir', 'gene_timeline_dir',
        'gene_distributions_dir', 'gene_comparisons_dir', 'gene_trajectories_dir',
        'gene_trajectories_mutant_dir', 'gene_trajectories_control1_dir',
        'gene_per_gene_dir', 'gene_analysis_dir'
    ]
    
    for key in expected_keys:
        assert key in all_paths, f"Missing key: {key}"
    print(f"  ✓ All {len(expected_keys)} expected keys present")
    
    # Check that paths are correctly formed
    assert all_paths['results_dir'] == "/test/results"
    assert all_paths['cell_metrics_dir'] == "/test/results/cell_metrics"
    assert all_paths['cell_jsd_mutant_dir'] == "/test/results/cell_metrics/individual_trajectories/jsd/mutant"
    print("  ✓ Path values correctly formed")
    
    print("  ✓ get_all_paths works correctly")

def test_error_handling():
    """Test error handling for invalid batch types."""
    print("Testing error handling...")
    
    plot_paths = PlotPaths("/test/results")
    
    # Test invalid batch type for cell JSD path
    try:
        plot_paths.get_individual_cell_jsd_path("invalid_batch", 1)
        assert False, "Should have raised ValueError for invalid batch"
    except ValueError as e:
        assert "Unexpected batch type: invalid_batch" in str(e)
        print("  ✓ Cell JSD path raises error for invalid batch")
    
    # Test invalid batch type for cell methylation path
    try:
        plot_paths.get_individual_cell_methylation_path("control2", 1)
        assert False, "Should have raised ValueError for control2 batch"
    except ValueError as e:
        assert "Unexpected batch type: control2" in str(e)
        print("  ✓ Cell methylation path raises error for control2 batch")
    
    # Test invalid batch type for gene JSD path
    try:
        plot_paths.get_individual_gene_jsd_path("unknown", 1)
        assert False, "Should have raised ValueError for unknown batch"
    except ValueError as e:
        assert "Unexpected batch type: unknown" in str(e)
        print("  ✓ Gene JSD path raises error for unknown batch")
    
    print("  ✓ Error handling works correctly")

if __name__ == "__main__":
    print("=" * 60)
    print("Testing PlotPaths class")
    print("=" * 60)
    
    test_directory_creation()
    test_path_generation()
    test_validation()
    test_get_all_paths()
    test_error_handling()
    
    print("\n" + "=" * 60)
    print("✅ All PlotPaths tests passed!")
    print("=" * 60)