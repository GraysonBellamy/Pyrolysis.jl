"""
Combined h+r adaptation driver.

`hr_adapt!` runs h-adaptivity (topology change via the experimental h-AMR
refine/coarsen path) followed by r-adaptivity (node relocation). It lives in
Experimental because it depends on the h-AMR cluster. The pure r-refinement
machinery (`relocate_nodes!`, monitors, `extract_state_arrays`) stays in core
`Adaptivity/r_refinement.jl` and is imported here.
"""

using ..Geometry: Mesh1D
using ..Adaptivity:
    relocate_nodes!,
    compute_indicator,
    MonitorFunction,
    GradientMonitor,
    AbstractErrorIndicator,
    GradientIndicator,
    extract_state_arrays

"""
	hr_adapt!(mesh::Mesh1D{NC};
			  h_indicator::AbstractErrorIndicator=GradientIndicator(:T),
			  r_monitor::MonitorFunction=GradientMonitor(),
			  θ_refine::Float64=0.7, θ_coarsen::Float64=0.1,
			  h_min::Float64=1e-6, h_max::Float64=1e-2,
			  max_level::Int=4, r_iterations::Int=5,
			  material=nothing) where NC

Perform combined h-refinement and r-refinement.

`material` is forwarded to the coarsening pass for the energy-conserving
temperature merge; without it sibling merges fall back to volume-weighted-T
averaging (not energy-conserving for mixed compositions).

# Strategy (from Design Doc Section 4.6)
1. Perform h-adaptivity first (change topology)
2. Perform r-adaptivity second (smooth and optimize)
"""
function hr_adapt!(
    mesh::Mesh1D{NC};
    h_indicator::AbstractErrorIndicator = GradientIndicator(:T),
    r_monitor::MonitorFunction = GradientMonitor(),
    θ_refine::Float64 = 0.7,
    θ_coarsen::Float64 = 0.1,
    h_min::Float64 = 1e-6,
    h_max::Float64 = 1e-2,
    max_level::Int = 4,
    r_iterations::Int = 5,
    material = nothing,
) where {NC}

    # 1. h-adaptivity
    T, ξ = extract_state_arrays(mesh)
    η = compute_indicator(h_indicator, mesh, T, ξ)

    # Refine
    marked_refine = mark_cells_for_refinement(mesh, η, θ_refine, h_min, max_level)
    marked_refine = enforce_2to1_balance!(mesh, marked_refine, max_level)
    perform_h_refinement!(mesh, marked_refine)

    # Coarsen
    T, ξ = extract_state_arrays(mesh)  # Re-extract after refinement
    η = compute_indicator(h_indicator, mesh, T, ξ)
    marked_coarsen = mark_cells_for_coarsening(mesh, η, θ_coarsen, h_max)
    perform_h_coarsening!(mesh, marked_coarsen; material)

    # 2. r-adaptivity
    if r_iterations > 0
        relocate_nodes!(mesh; monitor = r_monitor, max_iter = r_iterations, h_min = h_min)
    end

    return nothing
end
