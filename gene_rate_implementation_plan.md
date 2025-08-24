# Gene-Specific Methylation Rates Implementation Plan

## Overview
Implementation of gene-specific methylation rates allowing different groups of genes to have different methylation rates.

## Key Design Decisions
- Cells can have either uniform rate (backward compatible) OR gene-specific rates
- Cannot specify both rate types simultaneously  
- Gene rate groups format: `[(n_genes, rate), (n_genes, rate), ...]`
- Site rates are pre-computed and stored for efficiency
- All cells in a PetriDish must have same rate configuration

## Implementation Status

### ✅ Phase 1: Core Cell Class Updates (phase1/cell.py)
- [x] Updated Cell.__init__ to accept gene_rate_groups parameter
- [x] Added validation for rate specifications
- [x] Created _build_site_rates() method
- [x] Updated methylate() to use per-site rates
- [x] Updated to_dict() to save rate configuration and site_rates

### ✅ Phase 2: PetriDish Class Updates (phase1/cell.py)
- [x] Updated PetriDish.__init__ to accept gene_rate_groups
- [x] Added _validate_cell_rate_consistency() method
- [x] Updated cell creation logic for both rate types

### ✅ Phase 3: Pipeline Utils Updates (phase2/pipeline_utils.py)
- [x] Updated dict_to_cell() to handle both rate configurations
- [x] Updated load_petri_dish() to detect and propagate rate configuration
- [x] Added validation after loading cells

### ⏳ Phase 4: Command Line Updates (phase1/run_simulation.py)
Need to implement:
```python
# Add to argument parser
parser.add_argument('-r', '--rate', type=float, default=None, nargs='?', const=RATE,
                    help='Uniform methylation rate (default: 0.005 if flag used without value)')
parser.add_argument('--gene-rate-groups', type=str, default=None,
                    help='Gene-specific rates as "n1:rate1,n2:rate2,..."')

# Add parsing function
def parse_gene_rate_groups(groups_str: str, n_sites: int, gene_size: int) -> List[Tuple[int, float]]:
    """Parse gene rate groups string into list of tuples."""
    groups = []
    for i, group in enumerate(groups_str.split(',')):
        if ':' not in group:
            raise ValueError(f"Group {i} missing ':' separator: '{group}'")
        parts = group.split(':')
        if len(parts) != 2:
            raise ValueError(f"Group {i} should have format 'n:rate', got: '{group}'")
        n_genes = int(parts[0])
        rate = float(parts[1])
        if n_genes <= 0:
            raise ValueError(f"Group {i}: number of genes must be positive, got {n_genes}")
        if rate <= 0:
            raise ValueError(f"Group {i}: rate must be positive, got {rate}")
        groups.append((n_genes, rate))
    
    # Validate total
    total_genes = sum(n for n, _ in groups)
    expected_genes = n_sites // gene_size
    if total_genes != expected_genes:
        raise ValueError(f"gene-rate-groups specifies {total_genes} genes, but simulation has {expected_genes} genes")
    return groups

# Update main() to use new arguments
```

### ⏳ Phase 5: Test Updates
Need to add tests in:

#### test_comprehensive.py
```python
def test_gene_rate_groups():
    """Test gene-specific methylation rates."""
    # Test basic functionality
    groups = [(50, 0.01), (50, 0.02), (50, 0.03), (50, 0.04)]
    cell = Cell(n=1000, gene_rate_groups=groups, gene_size=5)
    
    # Verify site_rates built correctly
    assert len(cell.site_rates) == 1000
    for i in range(250):
        assert cell.site_rates[i] == 0.01
    for i in range(250, 500):
        assert cell.site_rates[i] == 0.02
    
    # Test daughter cells inherit rate structure
    daughter = cell.create_daughter_cell()
    assert daughter.gene_rate_groups == groups

def test_rate_validation():
    """Test validation of rate specifications."""
    # Can't specify both
    with pytest.raises(ValueError, match="Cannot specify both"):
        cell = Cell(n=1000, rate=0.01, gene_rate_groups=[(200, 0.01)])
    
    # Must specify one
    with pytest.raises(ValueError, match="Must specify either"):
        cell = Cell(n=1000)
    
    # Gene counts must match
    with pytest.raises(ValueError, match="150 genes.*200 genes"):
        cell = Cell(n=1000, gene_rate_groups=[(100, 0.01), (50, 0.02)])
```

#### test_edge_cases.py
```python
def test_gene_rate_edge_cases():
    """Test edge cases for gene rate groups."""
    # Single gene group (effectively uniform)
    cell = Cell(n=100, gene_rate_groups=[(20, 0.01)], gene_size=5)
    assert all(r == 0.01 for r in cell.site_rates)
    
    # Many small groups
    groups = [(1, 0.001 * i) for i in range(1, 21)]
    cell = Cell(n=100, gene_rate_groups=groups, gene_size=5)
    
    # Very high rates (>1 is allowed)
    cell = Cell(n=100, gene_rate_groups=[(20, 2.0)], gene_size=5)
    cell.methylate()
    
    # Mixed very low and very high
    groups = [(10, 0.0001), (10, 10.0)]
    cell = Cell(n=100, gene_rate_groups=groups, gene_size=5)
```

#### test_gene_jsd.py
```python
def test_gene_jsd_with_rate_groups():
    """Test gene JSD calculation with variable gene rates."""
    # Create PetriDish with different rates per gene group
    groups = [(5, 0.001), (5, 0.01), (5, 0.1), (5, 0.5)]
    petri = PetriDish(gene_rate_groups=groups, n=100, gene_size=5)
    petri.enable_history_tracking(start_year=0, track_gene_jsd=True)
    petri.grow_with_homeostasis(years=5, growth_phase=2, verbose=False)
    
    # Calculate gene JSDs
    gene_jsds = petri.calculate_gene_jsd()
    
    # High-rate genes should have higher JSDs
    low_rate_jsds = gene_jsds[0:5]   # First 5 genes (rate 0.001)
    high_rate_jsds = gene_jsds[15:20] # Last 5 genes (rate 0.5)
    assert max(low_rate_jsds) < min(high_rate_jsds)
```

### ⏳ Phase 6: Documentation Updates

#### CLAUDE.md
Add section:
```markdown
### Variable Methylation Rates
- Cells can have uniform rate (all genes methylate at same rate)
- OR gene-specific rates (different gene groups methylate at different rates)
- Specified via `--rate` (uniform) or `--gene-rate-groups` (variable)
- Example: `--gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055"`
- Cannot specify both flags simultaneously
- Total genes specified must equal n_sites/gene_size
```

## Command Examples

```bash
# Uniform rate (backward compatible)
python run_simulation.py --rate 0.005 --years 100

# Gene-specific rates
python run_simulation.py --gene-rate-groups "50:0.004,50:0.0045,50:0.005,50:0.0055" --years 100

# Error case (both specified)
python run_simulation.py --rate 0.005 --gene-rate-groups "50:0.004,50:0.0045"
# ERROR: Cannot specify both rate and gene-rate-groups
```

## Testing Plan
1. Test basic functionality with gene rate groups
2. Test validation (can't specify both, must specify one, counts must match)
3. Test edge cases (single group, many groups, high rates, mixed rates)
4. Test serialization/deserialization
5. Test with gene_JSD tracking
6. Test backward compatibility with uniform rate
7. Test phase2 pipeline compatibility