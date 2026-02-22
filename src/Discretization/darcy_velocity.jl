"""
Darcy velocity computation for gas flow through porous media.

This module computes the gas velocity field using Darcy's law:
	v_g = -(κ/μ) × ∇P

where:
- κ is the permeability [m²]
- μ is the gas viscosity [Pa·s]
- ∇P is the pressure gradient [Pa/m]

The velocity is computed at face centers using:
- Harmonic mean for permeability (conservative for heterogeneous media)
- Arithmetic mean for temperature → viscosity
- Central difference for pressure gradient
"""

using ..Core: GAS
using ..Properties: effective_gas_viscosity, harmonic_mean_permeability

# ============================================================================
# Single Face Darcy Velocity
# ============================================================================

"""
	darcy_velocity_face(κ_L::Float64, κ_R::Float64, 
						T_L::Float64, T_R::Float64,
						P_L::Float64, P_R::Float64, 
						d_LR::Float64) -> Float64

Compute Darcy velocity at a single face.

Uses harmonic mean for permeability and arithmetic mean for temperature
(for viscosity evaluation).

# Arguments
- `κ_L, κ_R`: Permeability in left and right cells [m²]
- `T_L, T_R`: Temperature in left and right cells [K]
- `P_L, P_R`: Pressure in left and right cells [Pa]
- `d_LR`: Distance between cell centers [m]

# Returns
Darcy velocity [m/s], positive = flow in +z direction (from L to R)
"""
@inline function darcy_velocity_face(κ_L::Float64, κ_R::Float64,
	T_L::Float64, T_R::Float64,
	P_L::Float64, P_R::Float64,
	d_LR::Float64)
	# Harmonic mean for permeability (handles heterogeneous interfaces)
	κ_face = harmonic_mean_permeability(κ_L, κ_R)

	# Arithmetic mean temperature for viscosity
	T_face = 0.5 * (T_L + T_R)
	μ_face = effective_gas_viscosity(T_face)

	# Pressure gradient (central difference)
	∇P = (P_R - P_L) / d_LR

	# Darcy velocity: v = -(κ/μ)∇P
	if κ_face > 0.0 && μ_face > 1e-20
		return -κ_face / μ_face * ∇P
	else
		return 0.0
	end
end

"""
	darcy_velocity_face(κ_face::Float64, μ_face::Float64,
						P_L::Float64, P_R::Float64, 
						d_LR::Float64) -> Float64

Compute Darcy velocity at a face given pre-computed face properties.

# Arguments
- `κ_face`: Face-averaged permeability [m²]
- `μ_face`: Face-averaged viscosity [Pa·s]
- `P_L, P_R`: Pressure in left and right cells [Pa]
- `d_LR`: Distance between cell centers [m]

# Returns
Darcy velocity [m/s], positive = flow in +z direction
"""
@inline function darcy_velocity_face(κ_face::Float64, μ_face::Float64,
	P_L::Float64, P_R::Float64,
	d_LR::Float64)
	# Pressure gradient
	∇P = (P_R - P_L) / d_LR

	# Darcy velocity
	if κ_face > 0.0 && μ_face > 1e-20
		return -κ_face / μ_face * ∇P
	else
		return 0.0
	end
end

# ============================================================================
# Full Mesh Darcy Velocity Computation
# ============================================================================

"""
	compute_darcy_velocity!(v_gas::Vector{Float64}, mesh::Mesh1D{NC},
							P::Vector{Float64},
							κ::Vector{Float64}) where NC

Compute Darcy gas velocity using pre-computed permeability field.

Simplified version that uses pre-computed permeability at cell centers.

# Arguments
- `v_gas`: Output array for face velocities [m/s]
- `mesh`: 1D mesh with cell states
- `P`: Pressure field at cell centers [Pa]
- `κ`: Permeability field at cell centers [m²]
"""
function compute_darcy_velocity!(v_gas::Vector{Float64}, mesh::Mesh1D{NC},
	P::Vector{Float64},
	κ::Vector{Float64}) where NC
	n_faces = length(mesh.faces)

	@inbounds for f in 1:n_faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		# Boundary faces: no flow
		if left_id == 0 || right_id == 0
			v_gas[f] = 0.0
			continue
		end

		# Get properties
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		κ_L = κ[left_id]
		κ_R = κ[right_id]

		# Distance between cell centers
		d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])

		# Compute Darcy velocity
		v_gas[f] = darcy_velocity_face(κ_L, κ_R, T_L, T_R,
			P[left_id], P[right_id], d_LR)
	end

	return nothing
end

# ============================================================================
# Darcy Mobility (κ/μ) Computation
# ============================================================================

"""
	compute_darcy_mobility!(mobility::Vector{Float64}, mesh::Mesh1D{NC},
							κ::Vector{Float64}) where NC

Compute Darcy mobility using pre-computed permeability field.

# Arguments
- `mobility`: Output array for face mobilities [m²/(Pa·s)]
- `mesh`: 1D mesh with cell states  
- `κ`: Permeability field at cell centers [m²]
"""
function compute_darcy_mobility!(mobility::Vector{Float64}, mesh::Mesh1D{NC},
	κ::Vector{Float64}) where NC
	n_faces = length(mesh.faces)

	@inbounds for f in 1:n_faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		if left_id == 0 || right_id == 0
			mobility[f] = 0.0
			continue
		end

		# Harmonic mean permeability
		κ_L = κ[left_id]
		κ_R = κ[right_id]
		κ_face = harmonic_mean_permeability(κ_L, κ_R)

		# Arithmetic mean temperature for viscosity
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		T_face = 0.5 * (T_L + T_R)
		μ_face = effective_gas_viscosity(T_face)

		if μ_face > 1e-20
			mobility[f] = κ_face / μ_face
		else
			mobility[f] = 0.0
		end
	end

	return nothing
end

# ============================================================================
# Velocity from Mobility (for efficient recomputation)
# ============================================================================

"""
	velocity_from_mobility!(v_gas::Vector{Float64}, mobility::Vector{Float64},
							∇P::Vector{Float64})

Compute Darcy velocity from pre-computed mobility and pressure gradient.

v_g = -mobility × ∇P

This is efficient when mobility is constant but pressure changes (implicit solve).

# Arguments
- `v_gas`: Output velocity array [m/s]
- `mobility`: Pre-computed κ/μ at faces [m²/(Pa·s)]
- `∇P`: Pressure gradient at faces [Pa/m]
"""
function velocity_from_mobility!(v_gas::Vector{Float64},
	mobility::Vector{Float64},
	∇P::Vector{Float64})
	n_faces = length(v_gas)

	@inbounds @simd for f in 1:n_faces
		v_gas[f] = -mobility[f] * ∇P[f]
	end

	return nothing
end

# ============================================================================
# Darcy Velocity Jacobian Coefficients
# ============================================================================

"""
	DarcyVelocityJacobianCoeffs

Struct holding Jacobian coefficients for Darcy velocity with respect to pressure.

For v_gas[f] = -(κ_face/μ_face) × (P_R - P_L) / d_LR:
- ∂v_gas[f]/∂P[left]  = +κ_face/(μ_face × d_LR)  (stored in coeff_L)
- ∂v_gas[f]/∂P[right] = -κ_face/(μ_face × d_LR)  (stored in coeff_R)

# Fields
- `coeff_L::Vector{Float64}`: ∂v_gas[f]/∂P[left_cell] for each face
- `coeff_R::Vector{Float64}`: ∂v_gas[f]/∂P[right_cell] for each face
- `left_ids::Vector{Int}`: Left cell index for each face (0 if boundary)
- `right_ids::Vector{Int}`: Right cell index for each face (0 if boundary)
"""
struct DarcyVelocityJacobianCoeffs
	coeff_L::Vector{Float64}
	coeff_R::Vector{Float64}
	left_ids::Vector{Int}
	right_ids::Vector{Int}
end

"""
	DarcyVelocityJacobianCoeffs(n_faces::Int)

Create empty Jacobian coefficient storage for n_faces faces.
"""
function DarcyVelocityJacobianCoeffs(n_faces::Int)
	DarcyVelocityJacobianCoeffs(
		zeros(n_faces),
		zeros(n_faces),
		zeros(Int, n_faces),
		zeros(Int, n_faces),
	)
end

"""
	compute_darcy_velocity_jacobian!(jac::DarcyVelocityJacobianCoeffs,
									 mesh::Mesh1D{NC},
									 κ::Vector{Float64}) where NC

Compute Jacobian coefficients using pre-computed permeability field.

# Arguments
- `jac`: Output coefficient structure
- `mesh`: 1D mesh with cell states
- `κ`: Permeability field at cell centers [m²]
"""
function compute_darcy_velocity_jacobian!(jac::DarcyVelocityJacobianCoeffs,
	mesh::Mesh1D{NC},
	κ::Vector{Float64}) where NC
	n_faces = length(mesh.faces)

	@inbounds for f in 1:n_faces
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids

		jac.left_ids[f] = left_id
		jac.right_ids[f] = right_id

		if left_id == 0 || right_id == 0
			jac.coeff_L[f] = 0.0
			jac.coeff_R[f] = 0.0
			continue
		end

		# Harmonic mean permeability
		κ_L = κ[left_id]
		κ_R = κ[right_id]
		κ_face = harmonic_mean_permeability(κ_L, κ_R)

		# Temperature and viscosity
		T_L = mesh.cell_states[left_id].T
		T_R = mesh.cell_states[right_id].T
		T_face = 0.5 * (T_L + T_R)
		μ_face = effective_gas_viscosity(T_face)

		# Distance
		d_LR = abs(mesh.cell_centers[right_id] - mesh.cell_centers[left_id])

		if μ_face > 1e-20 && d_LR > 0.0
			coeff = κ_face / (μ_face * d_LR)
		else
			coeff = 0.0
		end

		jac.coeff_L[f] = coeff
		jac.coeff_R[f] = -coeff
	end

	return nothing
end

# ============================================================================
# Pressure Boundary Condition Application
# ============================================================================

using ..Core: DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC
using ..Core: ImpermeableBC, AbstractMassBC
# ============================================================================
# Type-Stable BoundaryConditionSet Signatures
# ============================================================================

using ..Core: BoundaryConditionSet, get_mass_bc

"""
	apply_pressure_bc_to_darcy_velocity!(v_gas::Vector{Float64},
										 mesh::Mesh1D{NC},
										 P::Vector{Float64},
										 κ::Vector{Float64},
										 bc_set::BoundaryConditionSet,
										 top_face_ids::Vector{Int},
										 bottom_face_ids::Vector{Int}) where NC

Apply pressure boundary conditions using type-stable BoundaryConditionSet.

This is the preferred signature for maximum type stability and zero allocations.
Uses pre-cached face IDs computed once at problem setup.

# Arguments
"""


# Type-stable helper using pre-cached face IDs (avoids Dict lookup)
# boundary_id as Val for compile-time dispatch on position

function _apply_pressure_bc_faces!(v_gas, mesh, P, κ, face_ids::Vector{Int},
	bc::DirichletPressureBC, ::Val{B}) where B
	@inbounds for f in face_ids
		face = mesh.faces[f]

		# Find the interior cell (one of the cell_ids is 0 for boundary faces)
		cell_id = face.cell_ids[1] != 0 ? face.cell_ids[1] : face.cell_ids[2]

		# Get interior properties
		P_int = P[cell_id]
		κ_int = κ[cell_id]
		T_int = mesh.cell_states[cell_id].T

		# Distance from cell center to boundary face
		node_pos = mesh.nodes[face.node_id].z
		cell_center = mesh.cell_centers[cell_id]
		d = abs(node_pos - cell_center)

		if d < 1e-15
			v_gas[f] = 0.0
			continue
		end

		# Viscosity at interior temperature
		μ_int = effective_gas_viscosity(T_int)

		# Pressure gradient from interior to boundary
		# B=1 is top (node_pos > cell_center), B=2 is bottom (node_pos < cell_center)
		if B == 1  # Top boundary
			∇P = (bc.P - P_int) / d
		else  # Bottom boundary
			∇P = (P_int - bc.P) / d
		end

		# Darcy velocity
		if κ_int > 0.0 && μ_int > 1e-20
			v_gas[f] = -κ_int / μ_int * ∇P
		else
			v_gas[f] = 0.0
		end
	end
end

function _apply_pressure_bc_faces!(v_gas, mesh, P, κ, face_ids::Vector{Int},
	bc::NeumannPressureBC, ::Val{B}) where B
	@inbounds for f in face_ids
		if abs(bc.∂P∂n) < 1e-20
			v_gas[f] = 0.0
		else
			face = mesh.faces[f]
			cell_id = face.cell_ids[1] != 0 ? face.cell_ids[1] : face.cell_ids[2]

			κ_int = κ[cell_id]
			T_int = mesh.cell_states[cell_id].T
			μ_int = effective_gas_viscosity(T_int)

			if κ_int > 0.0 && μ_int > 1e-20
				v_gas[f] = -κ_int / μ_int * bc.∂P∂n
			else
				v_gas[f] = 0.0
			end
		end
	end
end

function _apply_pressure_bc_faces!(v_gas, mesh, P, κ, face_ids::Vector{Int},
	bc::ConvectivePressureBC, ::Val{B}) where B
	@inbounds for f in face_ids
		face = mesh.faces[f]
		cell_id = face.cell_ids[1] != 0 ? face.cell_ids[1] : face.cell_ids[2]
		P_int = P[cell_id]

		# Robin BC: v·n = h_P × (P_int - P_ext)
		if B == 1  # Top boundary: positive v = outflow
			v_gas[f] = bc.h_P * (P_int - bc.P_ext)
		else  # Bottom boundary: negative v = outflow
			v_gas[f] = -bc.h_P * (P_int - bc.P_ext)
		end
	end
end

# Fallback for other BC types - keep v_gas = 0 (default from compute_darcy_velocity!)
function _apply_pressure_bc_faces!(v_gas, mesh, P, κ, face_ids::Vector{Int},
	bc::AbstractMassBC, ::Val{B}) where B
	# No-op: ImpermeableBC and other mass BCs don't modify Darcy velocity
	return nothing
end

# CombinedMassBC: iterate through sub-BCs to find and apply pressure BCs
# This allows combining DiffusiveMassBC + DirichletPressureBC
# Use @generated to unroll at compile time and avoid allocation from heterogeneous tuple
@generated function _apply_pressure_bc_faces!(v_gas, mesh, P, κ, face_ids::Vector{Int},
	bc::CombinedMassBC{BCs}, boundary_val::Val{B}) where {BCs <: Tuple, B}
	N = length(BCs.parameters)
	exprs = [
		quote
			_apply_pressure_bc_faces!(v_gas, mesh, P, κ, face_ids, bc.bcs[$i], boundary_val)
		end for i in 1:N
	]
	return quote
		$(exprs...)
		return nothing
	end
end

using ..Core: CombinedMassBC

# Type-stable wrapper that applies mass/pressure BCs at the top and bottom
# faces using the concrete `BoundaryConditionSet` and precomputed face id lists.
function apply_pressure_bc_to_darcy_velocity!(v_gas::Vector{Float64}, mesh::Mesh1D{NC},
	P::Vector{Float64}, κ::Vector{Float64}, bc_set::BoundaryConditionSet,
	top_face_ids::Vector{Int}, bottom_face_ids::Vector{Int}) where NC
	# Top boundary (boundary id = 1)
	mass_top = get_mass_bc(bc_set, Val(1))
	_apply_pressure_bc_faces!(v_gas, mesh, P, κ, top_face_ids, mass_top, Val(1))

	# Bottom boundary (boundary id = 2)
	mass_bottom = get_mass_bc(bc_set, Val(2))
	_apply_pressure_bc_faces!(v_gas, mesh, P, κ, bottom_face_ids, mass_bottom, Val(2))

	return nothing
end