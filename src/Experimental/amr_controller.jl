"""
AMR controller.

The AMR controller orchestrates the full adaptation cycle:
1. Compute error indicators
2. Mark cells for refinement/coarsening
3. Enforce 2:1 balance
4. Perform h-refinement
5. Perform h-coarsening
6. Check and enforce mesh quality
"""

using ..Geometry: Mesh1D, update_geometry!, n_cells
using ..Adaptivity: extract_fields  # unexported core helper (error_indicators.jl)

# ============================================================================
# AMR Controller Type
# ============================================================================

"""
	AMRController{EI,QM}

Controls adaptive mesh refinement operations.

# Type Parameters
- `EI`: Error indicator type
- `QM`: Quality metric type

# Fields
- `error_indicator`: Error indicator for refinement decisions
- `quality_metric`: Mesh quality metric to enforce
- `θ_refine`: Refinement threshold (refine if η > θ_refine * max(η))
- `θ_coarsen`: Coarsening threshold (coarsen if η < θ_coarsen * max(η))
- `h_min`: Minimum allowed cell size
- `h_max`: Maximum allowed cell size
- `max_level`: Maximum refinement levels
- `check_every`: Timesteps between AMR checks
- `n_refinements`: Counter for total refinements performed
- `n_coarsenings`: Counter for total coarsenings performed
- `last_check_time`: Last simulation time when AMR was checked
"""
mutable struct AMRController{EI<:AbstractErrorIndicator,QM<:AbstractQualityMetric}
    error_indicator::EI
    quality_metric::QM
    θ_refine::Float64
    θ_coarsen::Float64
    h_min::Float64
    h_max::Float64
    max_level::Int
    check_every::Int
    n_refinements::Int
    n_coarsenings::Int
    last_check_time::Float64
end

"""
	AMRController(; kwargs...)

Construct an AMR controller with keyword arguments.

# Keyword Arguments
- `error_indicator`: Error indicator (default: `GradientIndicator(:T)`)
- `quality_metric`: Quality metric (default: `AspectRatioMetric()`)
- `θ_refine`: Refinement threshold (default: 0.7)
- `θ_coarsen`: Coarsening threshold (default: 0.1)
- `h_min`: Minimum cell size (default: 1e-6)
- `h_max`: Maximum cell size (default: 1e-2)
- `max_level`: Maximum refinement level (default: 4)
- `check_every`: Check frequency in timesteps (default: 10)
"""
function AMRController(;
    error_indicator::AbstractErrorIndicator = GradientIndicator(:T),
    quality_metric::AbstractQualityMetric = AspectRatioMetric(),
    θ_refine::Float64 = 0.7,
    θ_coarsen::Float64 = 0.1,
    h_min::Float64 = 1e-6,
    h_max::Float64 = 1e-2,
    max_level::Int = 4,
    check_every::Int = 10,
)
    AMRController(
        error_indicator,
        quality_metric,
        θ_refine,
        θ_coarsen,
        h_min,
        h_max,
        max_level,
        check_every,
        0,
        0,
        0.0,
    )
end

# ============================================================================
# Main Adaptation Interface
# ============================================================================

"""
	adapt_mesh!(mesh::Mesh1D{NC}, controller::AMRController;
				material=nothing) where NC -> NamedTuple

Perform full adaptive mesh refinement cycle.

`material` is forwarded to the coarsening pass for the energy-conserving
temperature merge; without it sibling merges fall back to volume-weighted-T
averaging (not energy-conserving for mixed compositions).

# Algorithm
1. Compute error indicators for all active cells
2. Mark cells for refinement (η > θ_refine * max(η), h > h_min, level < max_level)
3. Enforce 2:1 balance on marked cells
4. Perform h-refinement on marked cells
5. Recompute indicators on refined mesh
6. Mark sibling pairs for coarsening (both η < θ_coarsen * max(η), combined h < h_max)
7. Perform h-coarsening on marked pairs
8. Check mesh quality; apply smoothing if needed

# Returns
NamedTuple with:
- `n_refined`: Number of cells refined
- `n_coarsened`: Number of cell pairs coarsened
- `quality_ok`: Whether mesh passes quality check
"""
function adapt_mesh!(
    mesh::Mesh1D{NC},
    controller::AMRController;
    material = nothing,
) where {NC}
    # Extract current state
    T, ξ = extract_fields(mesh)

    # ========================================================================
    # 1. Compute error indicators
    # ========================================================================
    η = compute_indicator(controller.error_indicator, mesh, T, ξ)

    # ========================================================================
    # 2-4. h-Refinement
    # ========================================================================
    marked_refine = mark_cells_for_refinement(
        mesh,
        η,
        controller.θ_refine,
        controller.h_min,
        controller.max_level,
    )

    # Enforce 2:1 balance
    marked_refine = enforce_2to1_balance!(mesh, marked_refine, controller.max_level)

    # Perform refinement
    n_to_refine = length(marked_refine)
    perform_h_refinement!(mesh, marked_refine)

    controller.n_refinements += n_to_refine

    # ========================================================================
    # 5-7. h-Coarsening
    # ========================================================================
    # Recompute indicators after refinement
    T, ξ = extract_fields(mesh)
    η = compute_indicator(controller.error_indicator, mesh, T, ξ)

    marked_coarsen =
        mark_cells_for_coarsening(mesh, η, controller.θ_coarsen, controller.h_max)

    n_to_coarsen = length(marked_coarsen)
    perform_h_coarsening!(mesh, marked_coarsen; material)

    controller.n_coarsenings += n_to_coarsen

    # ========================================================================
    # 8. Quality enforcement
    # ========================================================================
    quality_ok = check_quality(mesh, controller.quality_metric)

    if !quality_ok
        remaining_violations = enforce_quality!(mesh, controller.quality_metric)
        quality_ok = remaining_violations == 0
    end

    return (n_refined = n_to_refine, n_coarsened = n_to_coarsen, quality_ok = quality_ok)
end

"""
	should_adapt(controller::AMRController, t::Float64, step::Int) -> Bool

Check if adaptation should be performed at current time/step.

# Arguments
- `controller`: AMR controller
- `t`: Current simulation time
- `step`: Current timestep number

# Returns
- `true` if adaptation should be performed
"""
function should_adapt(controller::AMRController, t::Float64, step::Int)
    return step % controller.check_every == 0
end

"""
	maybe_adapt!(mesh::Mesh1D{NC}, controller::AMRController, 
				 t::Float64, step::Int) where NC -> Union{Nothing,NamedTuple}

Conditionally perform adaptation based on check frequency.

# Returns
- `nothing` if adaptation was skipped
- Result of `adapt_mesh!` if adaptation was performed
"""
function maybe_adapt!(
    mesh::Mesh1D{NC},
    controller::AMRController,
    t::Float64,
    step::Int;
    material = nothing,
) where {NC}
    if !should_adapt(controller, t, step)
        return nothing
    end

    result = adapt_mesh!(mesh, controller; material)
    controller.last_check_time = t

    return result
end

# ============================================================================
# Statistics and Diagnostics
# ============================================================================

"""
	amr_statistics(controller::AMRController) -> NamedTuple

Get AMR statistics from controller.
"""
function amr_statistics(controller::AMRController)
    return (
        total_refinements = controller.n_refinements,
        total_coarsenings = controller.n_coarsenings,
        last_check_time = controller.last_check_time,
        θ_refine = controller.θ_refine,
        θ_coarsen = controller.θ_coarsen,
        h_min = controller.h_min,
        h_max = controller.h_max,
        max_level = controller.max_level,
    )
end

"""
	reset_statistics!(controller::AMRController)

Reset AMR statistics counters.
"""
function reset_statistics!(controller::AMRController)
    controller.n_refinements = 0
    controller.n_coarsenings = 0
    controller.last_check_time = 0.0
    return nothing
end

# ============================================================================
# Convenience Functions
# ============================================================================

"""
	create_default_controller(; h_target::Float64=1e-4) -> AMRController

Create an AMR controller with default settings tuned for pyrolysis.

The default indicator combines temperature gradient and curvature (curvature
targets the sharp fronts that reaction zones produce). A
`ReactionRateIndicator` cannot be used here: `adapt_mesh!` recomputes
indicators after topology changes, when any precomputed rates are stale, so
that indicator is standalone-only.

# Arguments
- `h_target`: Target cell size in reaction zones
"""
function create_default_controller(; h_target::Float64 = 1e-4)
    AMRController(
        error_indicator = CompositeIndicator(
            [GradientIndicator(:T, 1.0), CurvatureIndicator(:T, 0.5)],
            :max,
        ),
        quality_metric = CompositeQualityMetric([
            AspectRatioMetric(2.0),
            SmoothnessMetric(0.5),
            MinCellSizeMetric(h_target / 10),
        ]),
        θ_refine = 0.7,
        θ_coarsen = 0.1,
        h_min = h_target / 10,
        h_max = h_target * 100,
        max_level = 4,
        check_every = 10,
    )
end
