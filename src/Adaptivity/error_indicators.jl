"""
Error indicators for adaptive refinement.

Reference: Design Document Section 4.1

Error indicators estimate the local discretization error to guide
adaptive mesh refinement. High indicator values suggest refinement
is needed; low values suggest coarsening is safe.
"""

using ..Geometry: Mesh1D, cell_size, cell_center, n_cells

# ============================================================================
# Error Indicator Types
# ============================================================================

"""
	GradientIndicator <: AbstractErrorIndicator

Gradient-based error indicator: η_i = h_i · |∇φ|_i

Refines where gradients are steep (temperature fronts, reaction zones).

# Fields
- `variable::Symbol`: Variable to track (:T for temperature, or component name)
- `weight::Float64`: Weight for combining with other indicators
"""
struct GradientIndicator <: AbstractErrorIndicator
	variable::Symbol
	weight::Float64
end
GradientIndicator(var::Symbol) = GradientIndicator(var, 1.0)

"""
	CurvatureIndicator <: AbstractErrorIndicator

Curvature-based error indicator: η_i = h_i² · |∇²φ|_i

More sensitive to sharp fronts than gradient indicator.

# Fields
- `variable::Symbol`: Variable to track
- `weight::Float64`: Weight for combining with other indicators
"""
struct CurvatureIndicator <: AbstractErrorIndicator
	variable::Symbol
	weight::Float64
end
CurvatureIndicator(var::Symbol) = CurvatureIndicator(var, 1.0)

"""
	ReactionRateIndicator <: AbstractErrorIndicator

Reaction rate error indicator: η_i = h_i · max_j(|r_j|_i / max(ξ_j, ε))

Refines where reactions are active.

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
InterfaceProximityIndicator(pos::Float64, h::Float64) = InterfaceProximityIndicator(pos, h, 1.0)

"""
	CompositeIndicator <: AbstractErrorIndicator

Combines multiple error indicators.

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

# ============================================================================
# Gradient Indicator
# ============================================================================

function compute_indicator(ind::GradientIndicator, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64}) where NC
	n = n_cells(mesh)
	η = zeros(n)

	# Select the variable to track
	φ = if ind.variable == :T
		T
	else
		# Find component index by name
		comp_idx = findfirst(i -> mesh.cell_states[1].ξ isa NTuple ? true : false, 1:NC)
		# Use column of ξ matrix corresponding to component
		# For now, assume variable is an integer index or :T
		comp_idx = if ind.variable isa Integer
			ind.variable
		else
			1  # Default to first component
		end
		ξ[comp_idx, :]
	end

	@inbounds for i in 1:n
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
		# Central difference
		dz = mesh.cell_centers[right_id] - mesh.cell_centers[left_id]
		return (φ_R - φ_L) / dz
	end
end

# ============================================================================
# Curvature Indicator
# ============================================================================

function compute_indicator(ind::CurvatureIndicator, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64}) where NC
	n = n_cells(mesh)
	η = zeros(n)

	# Select the variable to track
	φ = if ind.variable == :T
		T
	else
		comp_idx = ind.variable isa Integer ? ind.variable : 1
		ξ[comp_idx, :]
	end

	@inbounds for i in 1:n
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

	# Second derivative: (φ_R - 2φ_C + φ_L) / h²
	h_L = mesh.cell_centers[i] - mesh.cell_centers[left_id]
	h_R = mesh.cell_centers[right_id] - mesh.cell_centers[i]
	h_avg = 0.5 * (h_L + h_R)

	curv = (φ[right_id] - 2.0 * φ[i] + φ[left_id]) / (h_avg^2)

	return curv
end

# ============================================================================
# Reaction Rate Indicator
# ============================================================================

function compute_indicator(ind::ReactionRateIndicator, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64};
	material = nothing, rates::Union{Nothing, Matrix{Float64}} = nothing) where NC
	n = n_cells(mesh)
	η = zeros(n)

	# If rates not provided, cannot compute (return zeros)
	if rates === nothing
		return η
	end

	ε = 1e-10  # Small value to avoid division by zero

	@inbounds for i in 1:n
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

function compute_indicator(ind::InterfaceProximityIndicator, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64}) where NC
	n = n_cells(mesh)
	η = zeros(n)

	@inbounds for i in 1:n
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

function compute_indicator(ind::CompositeIndicator, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64}; kwargs...) where NC
	n = n_cells(mesh)

	# Compute all component indicators
	component_ηs = [compute_indicator(sub_ind, mesh, T, ξ; kwargs...)
					for sub_ind in ind.indicators]

	# Combine based on method
	η = zeros(n)

	if ind.combination == :max
		@inbounds for i in 1:n
			η[i] = maximum(comp_η[i] for comp_η in component_ηs)
		end
	elseif ind.combination == :l2
		@inbounds for i in 1:n
			η[i] = sqrt(sum(comp_η[i]^2 for comp_η in component_ηs))
		end
	elseif ind.combination == :sum
		@inbounds for i in 1:n
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
function extract_fields(mesh::Mesh1D{NC}) where NC
	n = n_cells(mesh)
	T = zeros(n)
	ξ = zeros(NC, n)

	@inbounds for i in 1:n
		state = mesh.cell_states[i]
		T[i] = state.T
		for j in 1:NC
			ξ[j, i] = state.ξ[j]
		end
	end

	return T, ξ
end

"""
	compute_indicator(ind::AbstractErrorIndicator, mesh::Mesh1D) -> Vector{Float64}

Convenience method that extracts fields from mesh automatically.
"""
function compute_indicator(ind::AbstractErrorIndicator, mesh::Mesh1D{NC}) where NC
	T, ξ = extract_fields(mesh)
	return compute_indicator(ind, mesh, T, ξ)
end
