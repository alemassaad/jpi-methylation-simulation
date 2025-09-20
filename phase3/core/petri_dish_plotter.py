"""
PetriDish visualization module.
Handles all plotting for PetriDish objects from phase1 simulations.
Moved from phase1/cell.py to consolidate visualization in phase3.
"""
import statistics
import math
import os
from typing import Dict, Any, Optional, List
import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False

# Import phase1 modules
import sys
phase1_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'phase1')
sys.path.insert(0, phase1_path)
from cell import PetriDish, Cell


class PetriDishPlotter:
    """
    Handles all plotting for PetriDish objects.
    Reuses the logic from plot_history.py but in an object-oriented way.
    """
    
    def __init__(self, petri: PetriDish):
        """
        Initialize with a PetriDish that has history.
        
        Args:
            petri: PetriDish object with history to plot
            
        Raises:
            ValueError: If PetriDish has no history
        """
        self.petri = petri
        self.stats = None  # Cache calculated statistics
        
        # Validate that petri has history
        if not hasattr(petri, 'cell_history') or not petri.cell_history:
            raise ValueError("PetriDish has no cell history to plot")
    
    def calculate_statistics(self) -> Dict:
        """Calculate statistics from petri's history."""
        stats = {
            'years': [],
            'population_size': [],
            'cell_jsd': {
                'mean': [], 'median': [], 'p5': [], 'p25': [],
                'p75': [], 'p95': [], 'min': [], 'max': []
            },
            'methylation': {
                'mean': [], 'median': [], 'p5': [], 'p25': [],
                'p75': [], 'p95': []
            }
        }
        
        # Sort years numerically
        sorted_years = sorted([int(year) for year in self.petri.cell_history.keys()])
        
        for year in sorted_years:
            year_data = self.petri.cell_history[str(year)]
            
            # Skip if no cells (edge case)
            if not year_data:
                continue
            
            # Extract values (handle both naming conventions)
            jsd_values = []
            meth_values = []
            for cell in year_data:
                # Handle both 'cell_jsd' and 'cell_JSD' for compatibility
                if isinstance(cell, dict):
                    jsd = cell.get('cell_jsd', cell.get('cell_JSD', 0.0))
                    jsd_values.append(jsd)
                    # Handle both old and new keys
                    if 'cell_methylation_proportion' in cell:
                        meth_values.append(cell['cell_methylation_proportion'] * 100)
                    elif 'methylation_proportion' in cell:
                        # Legacy support
                        meth_values.append(cell['methylation_proportion'] * 100)
                    else:
                        # Calculate from methylated array if available
                        if 'methylated' in cell:
                            meth_prop = sum(cell['methylated']) / len(cell['methylated'])
                            meth_values.append(meth_prop * 100)
                        else:
                            meth_values.append(0.0)
                else:
                    # Handle Cell objects
                    jsd_values.append(getattr(cell, 'cell_jsd', 0.0))
                    meth_values.append(getattr(cell, 'cell_methylation_proportion', 
                                             getattr(cell, 'methylation_proportion', 0.0)) * 100)
            
            stats['years'].append(year)
            stats['population_size'].append(len(year_data))
            
            # JSD statistics
            if jsd_values:  # Ensure we have data
                stats['cell_jsd']['mean'].append(statistics.mean(jsd_values))
                stats['cell_jsd']['median'].append(statistics.median(jsd_values))
                
                # For single values, all percentiles are the same
                if len(jsd_values) == 1:
                    val = jsd_values[0]
                    stats['cell_jsd']['p5'].append(val)
                    stats['cell_jsd']['p25'].append(val)
                    stats['cell_jsd']['p75'].append(val)
                    stats['cell_jsd']['p95'].append(val)
                    stats['cell_jsd']['min'].append(val)
                    stats['cell_jsd']['max'].append(val)
                else:
                    # Use statistics.quantiles for percentiles
                    import numpy as np
                    stats['cell_jsd']['p5'].append(np.percentile(jsd_values, 5))
                    stats['cell_jsd']['p25'].append(np.percentile(jsd_values, 25))
                    stats['cell_jsd']['p75'].append(np.percentile(jsd_values, 75))
                    stats['cell_jsd']['p95'].append(np.percentile(jsd_values, 95))
                    stats['cell_jsd']['min'].append(min(jsd_values))
                    stats['cell_jsd']['max'].append(max(jsd_values))
            
            # Methylation statistics
            if meth_values:
                stats['methylation']['mean'].append(statistics.mean(meth_values))
                stats['methylation']['median'].append(statistics.median(meth_values))
                
                if len(meth_values) == 1:
                    val = meth_values[0]
                    stats['methylation']['p5'].append(val)
                    stats['methylation']['p25'].append(val)
                    stats['methylation']['p75'].append(val)
                    stats['methylation']['p95'].append(val)
                else:
                    import numpy as np
                    stats['methylation']['p5'].append(np.percentile(meth_values, 5))
                    stats['methylation']['p25'].append(np.percentile(meth_values, 25))
                    stats['methylation']['p75'].append(np.percentile(meth_values, 75))
                    stats['methylation']['p95'].append(np.percentile(meth_values, 95))
        
        self.stats = stats
        return stats
    
    def detect_growth_phase(self) -> int:
        """Auto-detect growth phase from population dynamics."""
        if not self.stats:
            self.calculate_statistics()
        
        pop_sizes = self.stats['population_size']
        
        for i in range(1, len(pop_sizes)):
            if i > 1:
                expected = 2 ** i
                actual = pop_sizes[i]
                if actual != expected:
                    return i - 1
        
        # Fallback to petri's growth_phase or default
        return getattr(self.petri, 'growth_phase', 13)
    
    def plot_jsd(self, title: str = None, output_path: str = None,
                 width: int = 1200, height: int = 500):
        """
        Create Cell JSD plot with secondary y-axis for cell count.
        
        Args:
            title: Optional title for the plot
            output_path: Optional path to save the plot
            width: Plot width in pixels
            height: Plot height in pixels
            
        Returns:
            Plotly figure object
        """
        # Import plotly here to avoid dependency issues
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not self.stats:
            self.calculate_statistics()
        
        # Create figure with secondary y-axis
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        years = self.stats['years']
        growth_phase = self.detect_growth_phase()
        
        # Add vertical line for growth phase
        if growth_phase > 0 and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text=f"End of growth phase",
                annotation_position="top"
            )
        
        # Add 5-95 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p95'] + self.stats['cell_jsd']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add 25-75 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p75'] + self.stats['cell_jsd']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add mean line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['cell_jsd']['mean'],
                mode='lines',
                name='Mean Cell JSD',
                line=dict(color='rgb(99, 110, 250)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Cell JSD: %{y:.4f}<br>Population: %{customdata}<extra></extra>',
                customdata=self.stats['population_size']
            ),
            secondary_y=False
        )
        
        # Add cell count on secondary y-axis
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>'
            ),
            secondary_y=True
        )
        
        # Update layout
        if title is None:
            title = f"Cell JSD Score vs Time"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            height=height,
            margin=dict(t=120),
            hovermode='x unified',
            showlegend=True,
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            ),
            template='plotly_white'
        )
        
        # Set axis titles
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell JSD Score", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
        
        # Add annotation
        final_idx = -1
        final_year = years[final_idx]
        final_pop = self.stats['population_size'][final_idx]
        
        annotation_text = (
            f"<b>Final Statistics (Year {final_year}):</b><br>"
            f"Population: {final_pop} cells<br>"
            f"Cell JSD Mean: {self.stats['cell_jsd']['mean'][final_idx]:.4f}<br>"
            f"Cell JSD 25-75%: [{self.stats['cell_jsd']['p25'][final_idx]:.4f}, {self.stats['cell_jsd']['p75'][final_idx]:.4f}]<br>"
            f"Cell JSD 5-95%: [{self.stats['cell_jsd']['p5'][final_idx]:.4f}, {self.stats['cell_jsd']['p95'][final_idx]:.4f}]"
        )
        
        if growth_phase > 0:
            annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
        
        fig.add_annotation(
            text=annotation_text,
            xref="paper", yref="paper",
            x=0.5, y=1.15,
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1,
            font=dict(size=10),
            align="center",
            xanchor="center",
            yanchor="top"
        )
        
        if output_path:
            fig.write_image(output_path, width=width, height=height, scale=2)
            
        return fig
    
    def plot_cell_methylation_proportion(self, title: str = None, output_path: str = None,
                        width: int = 1200, height: int = 500):
        """Create cell methylation proportion plot with secondary y-axis for cell count."""
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not self.stats:
            self.calculate_statistics()
        
        # Create figure with secondary y-axis
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        years = self.stats['years']
        growth_phase = self.detect_growth_phase()
        
        # Add vertical line for growth phase
        if growth_phase > 0 and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text=f"End of growth phase",
                annotation_position="top"
            )
        
        # Add 5-95 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p95'] + self.stats['methylation']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add 25-75 percentile band
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p75'] + self.stats['methylation']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add mean line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['methylation']['mean'],
                mode='lines',
                name='Mean Cell Methylation Proportion',
                line=dict(color='rgb(239, 85, 59)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Cell Methylation Proportion: %{y:.2f}%<br>Population: %{customdata}<extra></extra>',
                customdata=self.stats['population_size']
            ),
            secondary_y=False
        )
        
        # Add cell count on secondary y-axis
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>'
            ),
            secondary_y=True
        )
        
        # Update layout
        if title is None:
            title = f"Cell Methylation Proportion vs Time"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            height=height,
            margin=dict(t=120),
            hovermode='x unified',
            showlegend=True,
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            ),
            template='plotly_white'
        )
        
        # Set axis titles
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Methylation Proportion (%)", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
        
        # Add annotation
        final_idx = -1
        final_year = years[final_idx]
        final_pop = self.stats['population_size'][final_idx]
        
        annotation_text = (
            f"<b>Final Statistics (Year {final_year}):</b><br>"
            f"Population: {final_pop} cells<br>"
            f"Methylation Mean: {self.stats['methylation']['mean'][final_idx]:.2f}%<br>"
            f"Methylation 25-75%: [{self.stats['methylation']['p25'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p75'][final_idx]:.2f}%]<br>"
            f"Methylation 5-95%: [{self.stats['methylation']['p5'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p95'][final_idx]:.2f}%]"
        )
        
        if growth_phase > 0:
            annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
        
        fig.add_annotation(
            text=annotation_text,
            xref="paper", yref="paper",
            x=0.5, y=1.15,
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1,
            font=dict(size=10),
            align="center",
            xanchor="center",
            yanchor="top"
        )
        
        if output_path:
            fig.write_image(output_path, width=width, height=height, scale=2)
            
        return fig
    
    def plot_combined(self, title: str = None, output_path: str = None,
                     width: int = 1200, height: int = 800):
        """Create combined JSD and methylation plot."""
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            raise ImportError("Plotly is required for plotting. Install with: pip install plotly kaleido")
        
        if not self.stats:
            self.calculate_statistics()
        
        # Create subplots with secondary y-axes
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('Cell JSD Score vs Time', 'Methylation Proportion vs Time'),
            vertical_spacing=0.12,
            row_heights=[0.5, 0.5],
            specs=[[{"secondary_y": True}], [{"secondary_y": True}]]
        )
        
        years = self.stats['years']
        growth_phase = self.detect_growth_phase()
        
        # Add vertical lines
        if growth_phase > 0 and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text=f"End of growth phase",
                annotation_position="top",
                row=1, col=1
            )
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                row=2, col=1
            )
        
        # ===== JSD Plot (Row 1) =====
        # Add percentile bands
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p95'] + self.stats['cell_jsd']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            row=1, col=1, secondary_y=False
        )
        
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['cell_jsd']['p75'] + self.stats['cell_jsd']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(99, 110, 250, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            row=1, col=1, secondary_y=False
        )
        
        # Add mean JSD line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['cell_jsd']['mean'],
                mode='lines',
                name='Mean Cell JSD',
                line=dict(color='rgb(99, 110, 250)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Cell JSD: %{y:.4f}<extra></extra>'
            ),
            row=1, col=1, secondary_y=False
        )
        
        # Cell count for JSD plot
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
                showlegend=True
            ),
            row=1, col=1, secondary_y=True
        )
        
        # ===== Methylation Plot (Row 2) =====
        # Add percentile bands
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p95'] + self.stats['methylation']['p5'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.15)',
                line=dict(color='rgba(255,255,255,0)'),
                name='5-95 percentile',
                showlegend=False,
                hoverinfo='skip'
            ),
            row=2, col=1, secondary_y=False
        )
        
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=self.stats['methylation']['p75'] + self.stats['methylation']['p25'][::-1],
                fill='toself',
                fillcolor='rgba(239, 85, 59, 0.25)',
                line=dict(color='rgba(255,255,255,0)'),
                name='25-75 percentile',
                showlegend=False,
                hoverinfo='skip'
            ),
            row=2, col=1, secondary_y=False
        )
        
        # Add mean methylation line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['methylation']['mean'],
                mode='lines',
                name='Mean Methylation',
                line=dict(color='rgb(239, 85, 59)', width=2.5),
                hovertemplate='Year: %{x}<br>Mean Methylation: %{y:.2f}%<extra></extra>'
            ),
            row=2, col=1, secondary_y=False
        )
        
        # Cell count for methylation plot
        fig.add_trace(
            go.Scatter(
                x=years,
                y=self.stats['population_size'],
                mode='lines',
                name='Cell Count',
                line=dict(color='rgba(255, 127, 14, 0.7)', width=2, dash='dot'),
                hovertemplate='Year: %{x}<br>Cell Count: %{y}<extra></extra>',
                showlegend=False
            ),
            row=2, col=1, secondary_y=True
        )
        
        # Update layout
        if title is None:
            title = f"Simulation Results"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            height=height,
            margin=dict(t=120),
            hovermode='x unified',
            showlegend=True,
            legend=dict(
                yanchor="top", y=0.99, xanchor="left", x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            ),
            template='plotly_white'
        )
        
        # Update axes
        fig.update_xaxes(title_text="", row=1, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_xaxes(title_text="Age (years)", row=2, col=1, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell JSD Score", row=1, col=1, secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", row=1, col=1, secondary_y=True, showgrid=False)
        fig.update_yaxes(title_text="Methylation (%)", row=2, col=1, secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", row=2, col=1, secondary_y=True, showgrid=False)
        
        # Add annotation
        final_idx = -1
        final_year = years[final_idx]
        final_pop = self.stats['population_size'][final_idx]
        
        annotation_text = (
            f"<b>Final Statistics (Year {final_year}):</b><br>"
            f"Population: {final_pop} cells<br>"
            f"Cell JSD Mean: {self.stats['cell_jsd']['mean'][final_idx]:.4f}<br>"
            f"Cell JSD 25-75%: [{self.stats['cell_jsd']['p25'][final_idx]:.4f}, {self.stats['cell_jsd']['p75'][final_idx]:.4f}]<br>"
            f"Cell JSD 5-95%: [{self.stats['cell_jsd']['p5'][final_idx]:.4f}, {self.stats['cell_jsd']['p95'][final_idx]:.4f}]<br>"
            f"Methylation Mean: {self.stats['methylation']['mean'][final_idx]:.2f}%<br>"
            f"Methylation 25-75%: [{self.stats['methylation']['p25'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p75'][final_idx]:.2f}%]<br>"
            f"Methylation 5-95%: [{self.stats['methylation']['p5'][final_idx]:.2f}%, "
            f"{self.stats['methylation']['p95'][final_idx]:.2f}%]"
        )
        
        if growth_phase > 0:
            annotation_text += f"<br>Growth phase: Years 0-{growth_phase}"
        
        fig.add_annotation(
            text=annotation_text,
            xref="paper", yref="paper",
            x=0.5, y=1.15,
            showarrow=False,
            bgcolor="rgba(255, 255, 255, 0.9)",
            bordercolor="rgba(0, 0, 0, 0.2)",
            borderwidth=1,
            font=dict(size=10),
            align="center",
            xanchor="center",
            yanchor="top"
        )
        
        if output_path:
            fig.write_image(output_path, width=width, height=height, scale=2)
            
        return fig
    
    def plot_all(self, base_path: str, title_prefix: str = None):
        """
        Generate all three plot types.
        
        Args:
            base_path: Base path for output files (without extension)
            title_prefix: Optional prefix for plot titles
            
        Returns:
            Self for method chaining
        """
        title = title_prefix or "PetriDish Simulation"
        
        self.plot_jsd(title, f"{base_path}_jsd.png")
        self.plot_cell_methylation_proportion(title, f"{base_path}_cell_methylation_proportion.png")
        self.plot_combined(title, f"{base_path}_combined.png")
        
        # If gene JSD history exists, plot basic gene JSD (keep the existing simple plot)
        if hasattr(self.petri, 'gene_jsd_history') and self.petri.gene_jsd_history:
            self.plot_gene_jsd(title, f"{base_path}_gene_jsd.png")
            # Note: Advanced gene JSD plots (heatmap, rate comparison) are now generated in phase2
        
        return self
    
    def plot_gene_jsd(self, title: str = None, output_path: str = None,
                      width: int = 1200, height: int = 600):
        """
        Create Gene JSD plot showing per-gene JSD evolution over time.
        
        Args:
            title: Optional title for the plot
            output_path: Optional path to save the plot
            width: Plot width in pixels
            height: Plot height in pixels
        """
        try:
            import plotly.graph_objects as go
            from plotly.subplots import make_subplots
        except ImportError:
            print("Warning: Plotly not installed. Skipping gene JSD plot.")
            return
        
        # Check if gene JSD history exists
        if not hasattr(self.petri, 'gene_jsd_history') or not self.petri.gene_jsd_history:
            print("No gene JSD history available. Enable with track_gene_jsd=True")
            return
        
        # Extract data
        years = sorted([int(y) for y in self.petri.gene_jsd_history.keys()])
        n_genes = len(next(iter(self.petri.gene_jsd_history.values())))
        
        # Create figure with secondary y-axis for cell count
        fig = make_subplots(specs=[[{"secondary_y": True}]])
        
        # Calculate statistics for each year
        mean_gene_jsds = []
        median_gene_jsds = []
        p25_gene_jsds = []
        p75_gene_jsds = []
        p5_gene_jsds = []
        p95_gene_jsds = []
        
        for year in years:
            gene_jsds = self.petri.gene_jsd_history[str(year)]
            mean_gene_jsds.append(np.mean(gene_jsds))
            median_gene_jsds.append(np.median(gene_jsds))
            p25_gene_jsds.append(np.percentile(gene_jsds, 25))
            p75_gene_jsds.append(np.percentile(gene_jsds, 75))
            p5_gene_jsds.append(np.percentile(gene_jsds, 5))
            p95_gene_jsds.append(np.percentile(gene_jsds, 95))
        
        # Add percentile bands
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=p95_gene_jsds + p5_gene_jsds[::-1],
                fill='toself',
                fillcolor='rgba(31, 119, 180, 0.1)',
                line=dict(width=0),
                name='5-95 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        fig.add_trace(
            go.Scatter(
                x=years + years[::-1],
                y=p75_gene_jsds + p25_gene_jsds[::-1],
                fill='toself',
                fillcolor='rgba(31, 119, 180, 0.25)',
                line=dict(width=0),
                name='25-75 percentile',
                showlegend=True,
                hoverinfo='skip'
            ),
            secondary_y=False
        )
        
        # Add mean line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=mean_gene_jsds,
                mode='lines',
                name='Mean',
                line=dict(color='#1f77b4', width=2),
                hovertemplate='Year: %{x}<br>Mean Gene JSD: %{y:.4f}<extra></extra>'
            ),
            secondary_y=False
        )
        
        # Add median line
        fig.add_trace(
            go.Scatter(
                x=years,
                y=median_gene_jsds,
                mode='lines',
                name='Median',
                line=dict(color='#ff7f0e', width=2, dash='dash'),
                hovertemplate='Year: %{x}<br>Median Gene JSD: %{y:.4f}<extra></extra>'
            ),
            secondary_y=False
        )
        
        # Add cell count on secondary axis
        if self.stats and 'population_size' in self.stats:
            fig.add_trace(
                go.Scatter(
                    x=self.stats['years'],
                    y=self.stats['population_size'],
                    mode='lines',
                    name='Cell Count',
                    line=dict(color='gray', width=1, dash='dot'),
                    hovertemplate='Year: %{x}<br>Cells: %{y}<extra></extra>'
                ),
                secondary_y=True
            )
        
        # Add growth phase indicator
        growth_phase = self.detect_growth_phase()
        if growth_phase > 0 and years and growth_phase < years[-1]:
            fig.add_vline(
                x=growth_phase,
                line_dash="dash",
                line_color="gray",
                line_width=1,
                annotation_text="End of growth phase",
                annotation_position="top"
            )
        
        # Update layout
        if title is None:
            title = f"Gene JSD Score vs Time (Population-level, {n_genes} genes)"
        
        fig.update_layout(
            title=dict(text=title, font=dict(size=16)),
            xaxis_title="Age (years)",
            template='plotly_white',
            width=width,
            height=height,
            showlegend=True,
            legend=dict(
                yanchor="top",
                y=0.99,
                xanchor="left",
                x=0.01,
                bgcolor="rgba(255, 255, 255, 0.8)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            )
        )
        
        # Set axis titles
        fig.update_xaxes(showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Gene JSD Score", secondary_y=False, showgrid=True, gridcolor='rgba(0,0,0,0.1)')
        fig.update_yaxes(title_text="Cell Count", secondary_y=True, showgrid=False)
        
        # Add final statistics annotation
        if years:
            final_year = years[-1]
            final_mean = mean_gene_jsds[-1]
            final_median = median_gene_jsds[-1]
            
            annotation_text = (
                f"<b>Final Gene JSD Statistics (Year {final_year}):</b><br>"
                f"Mean: {final_mean:.4f}<br>"
                f"Median: {final_median:.4f}<br>"
                f"25-75%: [{p25_gene_jsds[-1]:.4f}, {p75_gene_jsds[-1]:.4f}]<br>"
                f"5-95%: [{p5_gene_jsds[-1]:.4f}, {p95_gene_jsds[-1]:.4f}]<br>"
                f"Genes tracked: {n_genes}"
            )
            
            fig.add_annotation(
                text=annotation_text,
                xref="paper", yref="paper",
                x=0.98, y=0.02,
                showarrow=False,
                font=dict(size=10, family="monospace"),
                align="right",
                xanchor="right",
                yanchor="bottom",
                bgcolor="rgba(255, 255, 255, 0.9)",
                bordercolor="rgba(0, 0, 0, 0.2)",
                borderwidth=1
            )
        
        # Save if path provided
        if output_path:
            os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
            fig.write_image(output_path, width=width, height=height, scale=2)
            print(f"Gene JSD plot saved to {output_path}")
        else:
            fig.show()
        
        return fig

    # Note: Additional plotting methods (heatmap, rate group comparison, trajectory, timeline)
    # have been omitted for brevity but would follow the same pattern