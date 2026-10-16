"""
Error indicators for adaptive refinement.

Error indicators estimate the local discretization error to guide
adaptive mesh refinement. High indicator values suggest refinement
is needed; low values suggest coarsening is safe.
"""

using ..Geometry: Mesh1D, cell_size, n_cells

# ============================================================================
# Error Indicator Types
# ============================================================================

"""
	GradientIndicator <: AbstractErrorIndicator

Gradient-based error indicator: η_i = h_i · |∇φ|_i

Refines where gradients are steep (temperature fronts, reaction zones).

# Fields
- `variable::Union{Symbol,Int}`: Variable to track — `:T` for temperature, or
  a component index `j` to track ξ[j, :]
- `weight::Float64`: Weight for combining with other indicators
"""
struct GradientIndicator <: AbstractErrorIndicator
    variable::Union{Symbol,Int}
    weight::Float64
end
GradientIndicator(var::Union{Symbol,Int}) = GradientIndicator(var, 1.0)

"""
	CurvatureIndicator <: AbstractErrorIndicator

Curvature-based error indicator: η_i = h_i² · |∇²φ|_i

More sensitive to sharp fronts than gradient indicator.

# Fields
- `variable::Union{Symbol,Int}`: Variable to track — `:T` for temperature, or
  a component index `j` to track ξ[j, :]
- `weight::Float64`: Weight for combining with other indicators
"""
struct CurvatureIndicator <: AbstractErrorIndicator
    variable::Union{Symbol,Int}
    weight::Float64
end
CurvatureIndicator(var::Union{Symbol,Int}) = CurvatureIndicator(var, 1.0)

"""
	ReactionRateIndicator <: AbstractErrorIndicator

Reaction rate error indicator: η_i = h_i · max_j(|r_j|_i / max(ξ_j, ε))

Refines where reactions are active.

Requires per-cell reaction rates: `compute_indicator(ind, mesh, T, ξ;
rates = r)` with `r` sized like ξ for the **current** mesh. It therefore
cannot run inside the `adapt_mesh!` cycle (which recomputes indicators after
topology changes, when any precomputed rates are stale) — use it standalone,
or supply gradient/curvature indicators to the controller instead.

# Fields
- `weight::Float64`: Weight for combining with other indicators
"""
struct ReactionRateIndicator <: AbstractErrorIndicator
    weight::Float64
end
ReactionRateIndicator() = ReactionRateIndicator(1.0)

"""
	InterfaceProximityIndicator <: AbstractErrorIndicator

Interface proximity indicator: η_i = 1 / (1 + d_interface / h_target)

Forces refinement near tracked interfaces.

# Fields
- `interface_position::Float64`: Current interface z-position
- `h_target::Float64`: Target cell size near interface
- `weight::Float64`: Weight for combining with other indicators
"""
struct InterfaceProximityIndicator <: AbstractErrorIndicator
    interface_position::Float64
    h_target::Float64
    weight::Float64
end
InterfaceProximityIndicator(pos::Float64, h::Float64) =
    InterfaceProximityIndicator(pos, h, 1.0)

"""
	CompositeIndicator <: AbstractErrorIndicator

Combines multiple error indicators.

Component indicators carry incommensurable dimensions (a temperature-gradient
indicator is in K, a reaction-rate indicator in 1/s, …), so each component is
**normalized by its own maximum over active cells** before combination. The
combined η is dimensionless; the component `weight`s act on these normalized
scales.

# Fields
- `indicators::Vector{AbstractErrorIndicator}`: Component indicators
- `combination::Symbol`: How to combine (:max, :l2, :sum)
"""
struct CompositeIndicator <: AbstractErrorIndicator
    indicators::Vector{AbstractErrorIndicator}
    combination::Symbol
end
CompositeIndicator(inds::Vector{<:AbstractErrorIndicator}) = CompositeIndicator(inds, :max)

# ============================================================================
# Indicator Computation - Main Interface
# ============================================================================

"""
	compute_indicator(ind::AbstractErrorIndicator, mesh::Mesh1D, 
					  T::Vector{Float64}, ξ::Matrix{Float64}) -> Vector{Float64}

Compute error indicator values for all active cells.

# Arguments
- `ind`: Error indicator type
- `mesh`: 1D mesh
- `T`: Temperature at cell centers [K]
- `ξ`: Mass concentrations [kg/m³], size (NC, n_cells)

# Returns
- Vector of indicator values for each cell
"""
function compute_indicator end

# Resolve the tracked field: `:T` → temperature; an Int j → ξ[j, :]. Any other
# Symbol is an error — component names are not resolvable here (no material in
# scope), and a silent ξ[1, :] fallback would track the wrong field.
function _indicator_field(
    variable::Union{Symbol,Int},
    T::Vector{Float64},
    ξ::Matrix{Float64},
    NC::Int,
)
    if variable === :T
        return T
    elseif variable isa Int
        1 <= variable <= NC ||
            throw(ArgumentError("indicator component index $variable out of range 1:$NC"))
        return @view ξ[variable, :]
    else
        throw(
            ArgumentError(
                "unknown indicator variable :$variable — use :T or a component index (Int)",
            ),
        )
    end
end

# ============================================================================
# Gradient Indicator
# ============================================================================

function compute_indicator(
    ind::GradientIndicator,
    mesh::Mesh1D{NC},
    T::Vector{Float64},
    ξ::Matrix{Float64};
    kwargs...,
) where {NC}
    n = n_cells(mesh)
    η = zeros(n)

    φ = _indicator_field(ind.variable, T, ξ, NC)

    @inbounds for i = 1:n
        if !mesh.active_cell_mask[i]
            continue
        end

        h = cell_size(mesh, i)

        # Compute gradient using finite differences
        grad = compute_gradient_1d(φ, mesh, i)

        η[i] = ind.weight * h * abs(grad)
    end

    return η
end

"""
	compute_gradient_1d(φ::Vector{Float64}, mesh::Mesh1D, i::Int) -> Float64

Compute gradient of field φ at cell i using central differences.
"""
function compute_gradient_1d(φ::AbstractVector{Float64}, mesh::Mesh1D, i::Int)
    cell = mesh.cells[i]
    left_face = mesh.faces[cell.face_ids[1]]
    right_face = mesh.faces[cell.face_ids[2]]

    left_id = left_face.cell_ids[1]
    right_id = right_face.cell_ids[2]

    # Get left and right values
    if left_id == 0  # Left boundary
        φ_L = φ[i]  # Use cell value (zero gradient at boundary approx)
    else
        φ_L = φ[left_id]
    end

    if right_id == 0  # Right boundary  
        φ_R = φ[i]  # Use cell value
    else
        φ_R = φ[right_id]
    end

    # Distance between neighbors
    if left_id == 0 && right_id == 0
        # Single cell - no gradient
        return 0.0
    elseif left_id == 0
        # One-sided (right)
        dz = mesh.cell_centers[right_id] - mesh.cell_centers[i]
        return (φ_R - φ[i]) / dz
    elseif right_id == 0
        # One-sided (left)
        dz = mesh.cell_centers[i] - mesh.cell_centers[left_id]
        return (φ[i] - φ_L) / dz
    else
        # Opposite-distance-weighted slopes: second order on nonuniform
        # meshes
        d_L = mesh.cell_centers[i] - mesh.cell_centers[left_id]
        d_R = mesh.cell_centers[right_id] - mesh.cell_centers[i]
        return opposite_distance_gradient(φ_L, φ[i], φ_R, d_L, d_R)
    end
end

# ============================================================================
# Curvature Indicator
# ============================================================================

function compute_indicator(
    ind::CurvatureIndicator,
    mesh::Mesh1D{NC},
    T::Vector{Float64},
    ξ::Matrix{Float64};
    kwargs...,
) where {NC}
    n = n_cells(mesh)
    η = zeros(n)

    φ = _indicator_field(ind.variable, T, ξ, NC)

    @inbounds for i = 1:n
        if !mesh.active_cell_mask[i]
            continue
        end

        h = cell_size(mesh, i)

        # Compute second derivative (curvature)
        curv = compute_curvature_1d(φ, mesh, i)

        η[i] = ind.weight * h^2 * abs(curv)
    end

    return η
end

"""
	compute_curvature_1d(φ::Vector{Float64}, mesh::Mesh1D, i::Int) -> Float64

Compute second derivative of field φ at cell i.
"""
function compute_curvature_1d(φ::AbstractVector{Float64}, mesh::Mesh1D, i::Int)
    cell = mesh.cells[i]
    left_face = mesh.faces[cell.face_ids[1]]
    right_face = mesh.faces[cell.face_ids[2]]

    left_id = left_face.cell_ids[1]
    right_id = right_face.cell_ids[2]

    # Need at least one neighbor on each side for curvature
    if left_id == 0 || right_id == 0
        return 0.0  # Can't compute curvature at boundaries
    end

    # Non-uniform second derivative. This reduces to the usual centered
    # stencil when h_L == h_R, while preserving zero curvature for linear
    # fields on graded meshes.
    h_L = mesh.cell_centers[i] - mesh.cell_centers[left_id]
    h_R = mesh.cell_centers[right_id] - mesh.cell_centers[i]

    if h_L <= 0.0 || h_R <= 0.0
        return 0.0
    end

    dφ_L = (φ[i] - φ[left_id]) / h_L
    dφ_R = (φ[right_id] - φ[i]) / h_R
    curv = 2.0 * (dφ_R - dφ_L) / (h_L + h_R)

    return curv
end

# ============================================================================
# Reaction Rate Indicator
# ============================================================================

function compute_indicator(
    ind::ReactionRateIndicator,
    mesh::Mesh1D{NC},
    T::Vector{Float64},
    ξ::Matrix{Float64};
    rates::Union{Nothing,Matrix{Float64}} = nothing,
    kwargs...,
) where {NC}
    n = n_cells(mesh)
    η = zeros(n)

    rates === nothing && throw(
        ArgumentError(
            "ReactionRateIndicator requires per-cell reaction rates: pass " *
            "`rates = r` (sized like ξ, for the current mesh). It cannot run " *
            "inside the adapt_mesh! cycle — see the indicator docstring.",
        ),
    )
    size(rates, 2) == n || throw(
        DimensionMismatch(
            "rates has $(size(rates, 2)) columns for $n cells — rates must be " *
            "computed on the current mesh",
        ),
    )

    ε = 1e-10  # Small value to avoid division by zero

    @inbounds for i = 1:n
        if !mesh.active_cell_mask[i]
            continue
        end

        h = cell_size(mesh, i)

        # Maximum normalized reaction rate
        max_rate = 0.0
        for j in axes(rates, 1)
            # Normalize by reactant concentration
            ξ_j = max(ξ[j, i], ε)
            normalized = abs(rates[j, i]) / ξ_j
            max_rate = max(max_rate, normalized)
        end

        η[i] = ind.weight * h * max_rate
    end

    return η
end

# ============================================================================
# Interface Proximity Indicator
# ============================================================================

function compute_indicator(
    ind::InterfaceProximityIndicator,
    mesh::Mesh1D{NC},
    T::Vector{Float64},
    ξ::Matrix{Float64};
    kwargs...,
) where {NC}
    n = n_cells(mesh)
    η = zeros(n)

    @inbounds for i = 1:n
        if !mesh.active_cell_mask[i]
            continue
        end

        z_c = mesh.cell_centers[i]
        d = abs(z_c - ind.interface_position)

        η[i] = ind.weight / (1.0 + d / ind.h_target)
    end

    return η
end

# ============================================================================
# Composite Indicator
# ============================================================================

function compute_indicator(
    ind::CompositeIndicator,
    mesh::Mesh1D{NC},
    T::Vector{Float64},
    ξ::Matrix{Float64};
    kwargs...,
) where {NC}
    n = n_cells(mesh)

    # Compute all component indicators
    component_ηs =
        [compute_indicator(sub_ind, mesh, T, ξ; kwargs...) for sub_ind in ind.indicators]

    # Normalize each component by its own max over active cells: the components
    # carry incommensurable dimensions (K vs 1/s vs …), so combining the raw
    # values would let whichever indicator has the largest unit-system numbers
    # dominate :max/:sum/:l2.
    for comp_η in component_ηs
        η_max = maximum(comp_η[i] for i = 1:n if mesh.active_cell_mask[i]; init = 0.0)
        if η_max > 0.0
            comp_η ./= η_max
        end
    end

    # Combine based on method
    η = zeros(n)

    if ind.combination == :max
        @inbounds for i = 1:n
            η[i] = maximum(comp_η[i] for comp_η in component_ηs)
        end
    elseif ind.combination == :l2
        @inbounds for i = 1:n
            η[i] = sqrt(sum(comp_η[i]^2 for comp_η in component_ηs))
        end
    elseif ind.combination == :sum
        @inbounds for i = 1:n
            η[i] = sum(comp_η[i] for comp_η in component_ηs)
        end
    else
        error("Unknown combination method: $(ind.combination)")
    end

    return η
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
	extract_fields(mesh::Mesh1D{NC}) -> (T::Vector{Float64}, ξ::Matrix{Float64})

Extract temperature and concentration fields from mesh cell states.
"""
function extract_fields(mesh::Mesh1D{NC}) where {NC}
    n = n_cells(mesh)
    T = zeros(n)
    ξ = zeros(NC, n)

    @inbounds for i = 1:n
        state = mesh.cell_states[i]
        T[i] = state.T
        for j = 1:NC
            ξ[j, i] = state.ξ[j]
        end
    end

    return T, ξ
end

"""
	compute_indicator(ind::AbstractErrorIndicator, mesh::Mesh1D; kwargs...) -> Vector{Float64}

Convenience method that extracts fields from mesh automatically.
"""
function compute_indicator(
    ind::AbstractErrorIndicator,
    mesh::Mesh1D{NC};
    kwargs...,
) where {NC}
    T, ξ = extract_fields(mesh)
    return compute_indicator(ind, mesh, T, ξ; kwargs...)
end
