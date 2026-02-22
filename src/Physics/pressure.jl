"""
Pressure field computation for gas phases.

Uses ideal gas law: P = Σⱼ (ξⱼ·R·T / Mⱼ) for gas components

Reference: Darcy_Fick_Implementation_Plan.md
"""

using ..Core: R_GAS, GAS

# ============================================================================
# Ideal Gas Pressure Computation
# ============================================================================

"""
	_ideal_gas_pressure_impl(ξ, T, comps::Tuple{C}, ::Val{1}) where C

Single-component ideal gas pressure helper.
"""
@inline function _ideal_gas_pressure_impl(ξ, T, comps::Tuple{C}, ::Val{1}) where C
	comp = comps[1]
	if comp.phase == GAS && comp.molar_mass > 0.0
		return ξ[1] * R_GAS * T / comp.molar_mass
	end
	return 0.0
end

"""
	_ideal_gas_pressure_impl(ξ, T, comps::Tuple{C1,C2}, ::Val{2}) where {C1,C2}

Two-component ideal gas pressure helper.
"""
@inline function _ideal_gas_pressure_impl(ξ, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	P = 0.0
	c1 = comps[1]
	if c1.phase == GAS && c1.molar_mass > 0.0
		P += ξ[1] * R_GAS * T / c1.molar_mass
	end
	c2 = comps[2]
	if c2.phase == GAS && c2.molar_mass > 0.0
		P += ξ[2] * R_GAS * T / c2.molar_mass
	end
	return P
end

"""
	_ideal_gas_pressure_impl(ξ, T, comps::Tuple{C1,C2,C3}, ::Val{3}) where {C1,C2,C3}

Three-component ideal gas pressure helper.
"""
@inline function _ideal_gas_pressure_impl(ξ, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	P = 0.0
	c1 = comps[1]
	if c1.phase == GAS && c1.molar_mass > 0.0
		P += ξ[1] * R_GAS * T / c1.molar_mass
	end
	c2 = comps[2]
	if c2.phase == GAS && c2.molar_mass > 0.0
		P += ξ[2] * R_GAS * T / c2.molar_mass
	end
	c3 = comps[3]
	if c3.phase == GAS && c3.molar_mass > 0.0
		P += ξ[3] * R_GAS * T / c3.molar_mass
	end
	return P
end

"""
	_ideal_gas_pressure_impl(ξ, T, comps::Tuple, ::Val{NC}) where NC

Generic N-component ideal gas pressure helper using generated function.
"""
@generated function _ideal_gas_pressure_impl(ξ, T, comps::Tuple, ::Val{NC}) where NC
	expr = quote
		P = 0.0
	end
	for i in 1:NC
		push!(expr.args, quote
			c = comps[$i]
			if c.phase == GAS && c.molar_mass > 0.0
				P += ξ[$i] * R_GAS * T / c.molar_mass
			end
		end)
	end
	push!(expr.args, :(return P))
	return expr
end

"""
	ideal_gas_pressure(ξ::NTuple{NC,Float64}, T::Float64, material) -> Float64

Compute ideal gas pressure from mass concentrations of gas components.

P = Σⱼ (ξⱼ × R × T / Mⱼ)   for gaseous components j

# Arguments
- `ξ`: Mass concentrations of all components [kg/m³]
- `T`: Temperature [K]
- `material`: Material definition with components

# Returns
- `P`: Pressure [Pa]

# Notes
Only components with `phase == GAS` and `molar_mass > 0` contribute to pressure.
Solid and liquid components are ignored.
"""
@inline function ideal_gas_pressure(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	return _ideal_gas_pressure_impl(ξ, T, material.components, Val(NC))
end

"""
	compute_pressure!(P::Vector{Float64}, mesh::Mesh1D{NC}, 
					  material::Material{NC}) where NC

Compute pressure field from ideal gas law for all cells.

# Arguments
- `P`: Output pressure vector [Pa], length = n_cells
- `mesh`: Mesh with current cell states (T, ξ)
- `material`: Material definition with components
"""
function compute_pressure!(P::Vector{Float64}, mesh, material)
	@inbounds for i in eachindex(mesh.cells)
		state = mesh.cell_states[i]
		P[i] = ideal_gas_pressure(state.ξ, state.T, material)
	end
	return nothing
end

# ============================================================================
# Pressure Gradient Computation
# ============================================================================

"""
	compute_pressure_gradient!(∇P::Vector{Float64}, P::Vector{Float64}, mesh)

Compute pressure gradient at cell faces using central differencing.

∇P_f = (P_R - P_L) / d_LR

# Arguments
- `∇P`: Output gradient vector [Pa/m], length = n_faces
- `P`: Cell-centered pressure values [Pa]
- `mesh`: Mesh with geometry (cell_centers)

# Notes
- Boundary faces (left_id=0 or right_id=0) are set to zero; boundary 
  conditions should be applied separately.
- Gradient is positive when pressure increases in +z direction.
"""
function compute_pressure_gradient!(∇P::Vector{Float64}, P::Vector{Float64}, mesh)
	@inbounds for f in eachindex(mesh.faces)
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			# Boundary face: set to zero, BC will override
			∇P[f] = 0.0
		else
			d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])
			∇P[f] = (P[right_id] - P[left_id]) / d_LR
		end
	end
	return nothing
end

# ============================================================================
# Pressure Derivatives (for Jacobian computation)
# ============================================================================

"""
	pressure_derivatives(ξ::NTuple{NC,Float64}, T::Float64, material) where NC

Compute partial derivatives of pressure with respect to state variables.

For P = Σⱼ (ξⱼ × R × T / Mⱼ):
- ∂P/∂T = Σⱼ (ξⱼ × R / Mⱼ)
- ∂P/∂ξⱼ = R × T / Mⱼ   (for gas components)
- ∂P/∂ξⱼ = 0             (for solid/liquid components)

# Arguments
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material definition

# Returns
- `∂P_∂T::Float64`: Derivative with respect to temperature [Pa/K]
- `∂P_∂ξ::NTuple{NC,Float64}`: Derivatives with respect to each concentration [Pa/(kg/m³)]
"""
@inline function pressure_derivatives(ξ::NTuple{NC, Float64}, T::Float64, material) where NC
	return _pressure_derivatives_impl(ξ, T, material.components, Val(NC))
end

@inline function _pressure_derivatives_impl(ξ, T, comps::Tuple{C}, ::Val{1}) where C
	c = comps[1]
	if c.phase == GAS && c.molar_mass > 0.0
		∂P_∂T = ξ[1] * R_GAS / c.molar_mass
		∂P_∂ξ = (R_GAS * T / c.molar_mass,)
	else
		∂P_∂T = 0.0
		∂P_∂ξ = (0.0,)
	end
	return ∂P_∂T, ∂P_∂ξ
end

@inline function _pressure_derivatives_impl(ξ, T, comps::Tuple{C1, C2}, ::Val{2}) where {C1, C2}
	∂P_∂T = 0.0
	∂P_∂ξ1 = 0.0
	∂P_∂ξ2 = 0.0

	c1 = comps[1]
	if c1.phase == GAS && c1.molar_mass > 0.0
		∂P_∂T += ξ[1] * R_GAS / c1.molar_mass
		∂P_∂ξ1 = R_GAS * T / c1.molar_mass
	end

	c2 = comps[2]
	if c2.phase == GAS && c2.molar_mass > 0.0
		∂P_∂T += ξ[2] * R_GAS / c2.molar_mass
		∂P_∂ξ2 = R_GAS * T / c2.molar_mass
	end

	return ∂P_∂T, (∂P_∂ξ1, ∂P_∂ξ2)
end

@inline function _pressure_derivatives_impl(ξ, T, comps::Tuple{C1, C2, C3}, ::Val{3}) where {C1, C2, C3}
	∂P_∂T = 0.0
	∂P_∂ξ1 = 0.0
	∂P_∂ξ2 = 0.0
	∂P_∂ξ3 = 0.0

	c1 = comps[1]
	if c1.phase == GAS && c1.molar_mass > 0.0
		∂P_∂T += ξ[1] * R_GAS / c1.molar_mass
		∂P_∂ξ1 = R_GAS * T / c1.molar_mass
	end

	c2 = comps[2]
	if c2.phase == GAS && c2.molar_mass > 0.0
		∂P_∂T += ξ[2] * R_GAS / c2.molar_mass
		∂P_∂ξ2 = R_GAS * T / c2.molar_mass
	end

	c3 = comps[3]
	if c3.phase == GAS && c3.molar_mass > 0.0
		∂P_∂T += ξ[3] * R_GAS / c3.molar_mass
		∂P_∂ξ3 = R_GAS * T / c3.molar_mass
	end

	return ∂P_∂T, (∂P_∂ξ1, ∂P_∂ξ2, ∂P_∂ξ3)
end

@generated function _pressure_derivatives_impl(ξ, T, comps::Tuple, ::Val{NC}) where NC
	# Build the body of the function
	expr = quote
		∂P_∂T = 0.0
	end

	# Create variables for each derivative
	for i in 1:NC
		var = Symbol("∂P_∂ξ", i)
		push!(expr.args, :($var = 0.0))
	end

	# Compute derivatives for each component
	for i in 1:NC
		var = Symbol("∂P_∂ξ", i)
		push!(expr.args, quote
			c = comps[$i]
			if c.phase == GAS && c.molar_mass > 0.0
				∂P_∂T += ξ[$i] * R_GAS / c.molar_mass
				$var = R_GAS * T / c.molar_mass
			end
		end)
	end

	# Build the return tuple
	vars = [Symbol("∂P_∂ξ", i) for i in 1:NC]
	tuple_expr = Expr(:tuple, vars...)
	push!(expr.args, :(return ∂P_∂T, $tuple_expr))

	return expr
end

# ============================================================================
# In-place Pressure Derivatives (allocation-free for Jacobian computation)
# ============================================================================

"""
	pressure_derivatives!(∂P_∂ξ::AbstractVector, ξ::NTuple{NC,Float64}, T::Float64, material) where NC

Compute partial derivatives of pressure in-place (allocation-free).

Same as `pressure_derivatives` but writes to a pre-allocated vector instead of
returning a tuple. This eliminates heap allocations when called in a loop.

# Arguments
- `∂P_∂ξ`: Output vector of length NC for concentration derivatives [Pa/(kg/m³)]
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material definition

# Returns
- `∂P_∂T::Float64`: Derivative with respect to temperature [Pa/K]

# Notes
The output vector `∂P_∂ξ` is modified in-place.
"""
@inline function pressure_derivatives!(∂P_∂ξ::AbstractVector, ξ::NTuple{NC, Float64},
	T::Float64, material) where NC
	return _pressure_derivatives_inplace!(∂P_∂ξ, ξ, T, material.components, Val(NC))
end

@generated function _pressure_derivatives_inplace!(∂P_∂ξ::AbstractVector, ξ, T,
	comps::Tuple, ::Val{NC}) where NC
	# Build the body of the function
	expr = quote
		∂P_∂T = 0.0
	end

	# Compute derivatives for each component and write to output vector
	for i in 1:NC
		push!(expr.args, quote
			c = comps[$i]
			if c.phase == GAS && c.molar_mass > 0.0
				∂P_∂T += ξ[$i] * R_GAS / c.molar_mass
				@inbounds ∂P_∂ξ[$i] = R_GAS * T / c.molar_mass
			else
				@inbounds ∂P_∂ξ[$i] = 0.0
			end
		end)
	end

	push!(expr.args, :(return ∂P_∂T))

	return expr
end
