"""
	Adaptivity

Adaptive mesh refinement for the Pyrolysis.jl package.

This module provides:
- Error indicators (gradient, curvature, reaction rate, composite)
- h-refinement (cell splitting)
- h-coarsening (cell merging)
- r-refinement (node relocation via equidistribution)
- Mesh quality metrics (aspect ratio, smoothness, Jacobian)
- AMR controller for orchestrating adaptation
- Cell depletion detection and handling (both surface and interior)

## Usage

```julia
using Pyrolysis.Adaptivity

# Create error indicator
indicator = GradientIndicator(:T)

# Compute indicators
η = compute_indicator(indicator, mesh)

# Create AMR controller
controller = AMRController(
	error_indicator = indicator,
	θ_refine = 0.7,
	θ_coarsen = 0.1,
	h_min = 1e-6,
	max_level = 4
)

# Adapt mesh
result = adapt_mesh!(mesh, controller)

# Handle cell depletion (any cell, not just surface)
depleted = find_depleted_cells(mesh, material, 0.01)
for cell_id in depleted
	neighbor = find_merge_neighbor(mesh, cell_id, material)
	if neighbor > 0
		merge_adjacent_cells!(mesh, cell_id, neighbor, material)
	end
end
```

## References
Design Document Sections 4, 6.6, 7.4-7.6
Cell_Depletion_Implementation_Plan.md
"""
module Adaptivity

using ..Core
using ..Geometry
using ..Geometry: Mesh1D, Node1D, Cell1D, Face1D, CellState
using ..Geometry: cell_size, cell_center, update_geometry!, n_cells, n_nodes
using ..Geometry: domain_bounds

# Include subfiles
include("error_indicators.jl")
include("h_refinement.jl")
include("h_coarsening.jl")
include("r_refinement.jl")
include("mesh_quality.jl")
include("amr_controller.jl")
include("surface_depletion.jl")

# ============================================================================
# Exports
# ============================================================================

# Error indicators
export AbstractErrorIndicator
export GradientIndicator, CurvatureIndicator, ReactionRateIndicator
export InterfaceProximityIndicator, CompositeIndicator
export compute_indicator

# h-refinement
export refine_cell!, perform_h_refinement!
export mark_cells_for_refinement, enforce_2to1_balance!
export verify_refinement_conservation

# h-coarsening
export coarsen_cells!, perform_h_coarsening!
export mark_cells_for_coarsening
export can_coarsen, get_sibling
export verify_coarsening_conservation
# Extended h-coarsening for adjacent cell merging
export merge_adjacent_cells!, merge_surface_cell!
export compact_mesh_after_surface_merge!, compact_mesh_after_merge!
export find_shared_face, are_cells_adjacent
export merge_temperatures_conservatively, merge_temperatures_conservatively_from_values

# r-refinement
export MonitorFunction, GradientMonitor, InterfaceMonitor, CompositeMonitor
export relocate_nodes!, hr_adapt!
export evaluate_monitor, compute_target_positions

# Mesh quality
export AbstractQualityMetric
export AspectRatioMetric, SmoothnessMetric, JacobianMetric
export MinCellSizeMetric, CompositeQualityMetric
export check_quality, quality_violations
export compute_aspect_ratios, compute_smoothness
export enforce_quality!, mesh_quality_summary

# AMR controller
export AMRController
export adapt_mesh!, should_adapt, maybe_adapt!
export amr_statistics, reset_statistics!
export create_default_controller

# Cell depletion utilities
export find_surface_cell, find_bottom_cell
export find_interior_neighbor, find_exterior_neighbor
export find_merge_neighbor

end # module Adaptivity
