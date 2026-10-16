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
- Depleted-cell parking (surface and interior cells)

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

# Park a depleted cell: merge it into its neighbour and move the freed slot
# above the live surface
neighbor = find_merge_neighbor(mesh, cell_id)
park_cell!(mesh, cell_id, neighbor; material = material, ale = true)
```
"""
module Adaptivity

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Geometry: Mesh1D, Node1D, Cell1D, Face1D, CellState
using ..Geometry: cell_size, update_geometry!, n_cells, n_nodes
using ..Geometry: domain_bounds

# Include subfiles
include("error_indicators.jl")

# ============================================================================
# AMR invariants
# ============================================================================
#
# Refinement appends child cells but never compacts the backing arrays, so
# `n_active_cells < length(mesh.cells)` is expected after refinement. This
# assertion catches bookkeeping drift — e.g. an `n_active_cells` counter that
# disagrees with the `active_cell_mask`, or parallel arrays of mismatched
# length.
"""
	assert_amr_invariants(mesh::Mesh1D)

Cheap consistency check on AMR bookkeeping. Throws an `AssertionError` if
`n_active_cells` does not match the active mask, or if any of the parallel
per-cell arrays disagree in length. Intended to be called after every
refinement / coarsening step.
"""
function assert_amr_invariants(mesh::Mesh1D)
    n = length(mesh.cells)
    @assert length(mesh.active_cell_mask) == n "active_cell_mask length $(length(mesh.active_cell_mask)) ≠ cells $(n)"
    @assert length(mesh.cell_states) == n "cell_states length $(length(mesh.cell_states)) ≠ cells $(n)"
    @assert length(mesh.cell_volumes) == n "cell_volumes length $(length(mesh.cell_volumes)) ≠ cells $(n)"
    @assert length(mesh.cell_centers) == n "cell_centers length $(length(mesh.cell_centers)) ≠ cells $(n)"
    # The χ denominator is a fifth parallel per-cell array (empty only for
    # hand-built meshes that never sized it).
    @assert isempty(mesh.m_dry_solid_initial) || length(mesh.m_dry_solid_initial) == n "m_dry_solid_initial length $(length(mesh.m_dry_solid_initial)) ≠ cells $(n)"
    @assert mesh.n_active_cells == count(mesh.active_cell_mask) "n_active_cells $(mesh.n_active_cells) ≠ count(active_cell_mask) $(count(mesh.active_cell_mask))"
    return nothing
end

include("cell_merge.jl")
include("cell_park.jl")
include("r_refinement.jl")
include("mesh_quality.jl")
include("surface_depletion.jl")

# ============================================================================
# Exports
# ============================================================================

# Error indicators
export AbstractErrorIndicator
export GradientIndicator, CurvatureIndicator, ReactionRateIndicator
export InterfaceProximityIndicator, CompositeIndicator
export compute_indicator

# Depleted-cell parking (cell_park.jl) and the conservative merge it uses
export park_cell!
export merge_temperatures_conservatively

# r-refinement
export MonitorFunction, GradientMonitor, InterfaceMonitor, CompositeMonitor
export relocate_nodes!
export evaluate_monitor, compute_target_positions

# Mesh quality
export AbstractQualityMetric
export AspectRatioMetric, SmoothnessMetric, JacobianMetric
export MinCellSizeMetric, CompositeQualityMetric
export check_quality, quality_violations
export compute_aspect_ratios, compute_smoothness
export enforce_quality!, mesh_quality_summary

# Cell depletion utilities
export find_interior_neighbor, find_exterior_neighbor
export find_merge_neighbor

# AMR invariant checks
export assert_amr_invariants

end # module Adaptivity
