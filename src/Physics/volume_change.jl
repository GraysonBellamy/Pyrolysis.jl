"""
Volume change computations.

Reference: Design Document Section 4.6

PERFORMANCE NOTE: All functions that iterate over material.components use
dispatch-based unrolling (explicit overloads for small tuple sizes) to avoid
allocations from indexing heterogeneous tuples with runtime indices.
"""

# ============================================================================
# Dispatch-based volume change rate helpers (allocation-free)
# ============================================================================

@inline function _volume_change_rate(dξ, T, comps::Tuple{C}, ::Val{1}) where C
	comp = comps[1]
	γ = comp.swelling_factor
	if γ > 0
		ρⱼ = comp.density(T)
		return γ * dξ[1] / ρⱼ
	end
	return 0.0
end

@inline function _volume_change_rate(dξ, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	θ = 0.0
	γ1 = comps[1].swelling_factor
	if γ1 > 0
		ρ1 = comps[1].density(T)
		θ += γ1 * dξ[1] / ρ1
	end
	γ2 = comps[2].swelling_factor
	if γ2 > 0
		ρ2 = comps[2].density(T)
		θ += γ2 * dξ[2] / ρ2
	end
	return θ
end

@inline function _volume_change_rate(dξ, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	θ = 0.0
	γ1 = comps[1].swelling_factor
	if γ1 > 0
		ρ1 = comps[1].density(T)
		θ += γ1 * dξ[1] / ρ1
	end
	γ2 = comps[2].swelling_factor
	if γ2 > 0
		ρ2 = comps[2].density(T)
		θ += γ2 * dξ[2] / ρ2
	end
	γ3 = comps[3].swelling_factor
	if γ3 > 0
		ρ3 = comps[3].density(T)
		θ += γ3 * dξ[3] / ρ3
	end
	return θ
end

# @generated version for NC > 3 components (allocation-free via compile-time unrolling)
@generated function _volume_change_rate(dξ, T, comps::Tuple{Vararg{Any, NC}}, ::Val{NC}) where NC
	# Build the unrolled expression at compile time
	exprs = Expr[]
	push!(exprs, :(θ = zero(typeof(T))))
	for i in 1:NC
		γ_sym = Symbol("γ_", i)
		ρ_sym = Symbol("ρ_", i)
		push!(exprs, quote
			@inbounds begin
				$γ_sym = comps[$i].swelling_factor
				if $γ_sym > 0
					$ρ_sym = comps[$i].density(T)
					θ += $γ_sym * dξ[$i] / $ρ_sym
				end
			end
		end)
	end
	push!(exprs, :(return θ))
	return Expr(:block, exprs...)
end

"""
	volume_change_rate(dξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute relative volume change rate [(m³/m³)/s].

dV/dt / V = Σⱼ (1/ρⱼ) (dξⱼ/dt)

Only includes contributions from components with swelling.

PERFORMANCE: Uses dispatch-based helpers to avoid allocations from heterogeneous tuple access.
"""
@inline function volume_change_rate(dξ::NTuple{NC, Float64}, T::Float64, material) where NC
	return _volume_change_rate(dξ, T, material.components, Val(NC))
end

"""
	volume_change_rate(dξ::Vector{Float64}, T::Float64, material::Material{NC}) -> Float64

Vector version for pre-allocated storage.

PERFORMANCE: Uses dispatch-based helpers to avoid allocations from heterogeneous tuple access.
"""
@inline function volume_change_rate(dξ::AbstractVector{Float64}, T::Float64, material::Material{NC}) where NC
	return _volume_change_rate(dξ, T, material.components, Val(NC))
end

"""
	compute_cell_volume_rates!(dV::Vector{Float64}, mesh::Mesh1D{NC}, 
							   dξ_matrix::Matrix{Float64}, material::Material{NC}) where NC

Compute volume change rate for all cells.

# Arguments
- `dV`: Output vector of volume change rates [m³/s]
- `mesh`: Mesh with current state
- `dξ_matrix`: Species source terms dξ_matrix[j,i] for component j in cell i
- `material`: Material definition
"""
function compute_cell_volume_rates!(dV::Vector{Float64}, mesh::Mesh1D{NC},
	dξ_matrix::Matrix{Float64}, material::Material{NC}) where NC
	n_cells_val = length(mesh.cells)

	@inbounds for i in 1:n_cells_val
		state = mesh.cell_states[i]
		V = state.volume
		T = state.T

		# Get dξ/dt for this cell
		dξ = view(dξ_matrix, :, i)

		# Relative volume change rate
		θ = volume_change_rate(dξ, T, material)

		# Absolute volume change
		dV[i] = V * θ
	end

	return nothing
end
