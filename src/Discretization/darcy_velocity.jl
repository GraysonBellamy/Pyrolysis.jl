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

using ..Materials: GAS
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
@inline function darcy_velocity_face(
    κ_L::Float64,
    κ_R::Float64,
    T_L::Float64,
    T_R::Float64,
    P_L::Float64,
    P_R::Float64,
    d_LR::Float64,
)
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
@inline function darcy_velocity_face(
    κ_face::Float64,
    μ_face::Float64,
    P_L::Float64,
    P_R::Float64,
    d_LR::Float64,
)
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
function compute_darcy_velocity!(
    v_gas::AbstractVector,
    state,
    mesh,
    P::AbstractVector,
    κ::AbstractVector,
)
    n_faces = length(state.A_face)

    @inbounds for f = 1:n_faces
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            v_gas[f] = zero(eltype(v_gas))
            continue
        end

        T_L = state.T[left_id]
        T_R = state.T[right_id]
        κ_L = κ[left_id]
        κ_R = κ[right_id]

        d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])

        v_gas[f] = darcy_velocity_face(κ_L, κ_R, T_L, T_R, P[left_id], P[right_id], d_LR)
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
function compute_darcy_mobility!(mobility::AbstractVector, state, mesh, κ::AbstractVector)
    n_faces = length(state.A_face)

    @inbounds for f = 1:n_faces
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        if left_id == 0 || right_id == 0
            mobility[f] = zero(eltype(mobility))
            continue
        end

        κ_L = κ[left_id]
        κ_R = κ[right_id]
        κ_face = harmonic_mean_permeability(κ_L, κ_R)

        T_L = state.T[left_id]
        T_R = state.T[right_id]
        T_face = (T_L + T_R) / 2
        μ_face = effective_gas_viscosity(T_face)

        if μ_face > 1e-20
            mobility[f] = κ_face / μ_face
        else
            mobility[f] = zero(eltype(mobility))
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
function velocity_from_mobility!(
    v_gas::Vector{Float64},
    mobility::Vector{Float64},
    ∇P::Vector{Float64},
)
    n_faces = length(v_gas)

    @inbounds @simd for f = 1:n_faces
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
function compute_darcy_velocity_jacobian!(
    jac::DarcyVelocityJacobianCoeffs,
    state,
    mesh,
    κ::AbstractVector,
)
    n_faces = length(state.A_face)

    @inbounds for f = 1:n_faces
        face = mesh.faces[f]
        left_id, right_id = face.cell_ids

        jac.left_ids[f] = left_id
        jac.right_ids[f] = right_id

        if left_id == 0 || right_id == 0
            jac.coeff_L[f] = 0.0
            jac.coeff_R[f] = 0.0
            continue
        end

        κ_L = κ[left_id]
        κ_R = κ[right_id]
        κ_face = harmonic_mean_permeability(κ_L, κ_R)

        T_L = state.T[left_id]
        T_R = state.T[right_id]
        T_face = (T_L + T_R) / 2
        μ_face = effective_gas_viscosity(T_face)

        d_LR = abs(state.cell_centers[right_id] - state.cell_centers[left_id])

        if μ_face > 1e-20 && d_LR > 0.0
            coeff = κ_face / (μ_face * d_LR)
        else
            coeff = zero(typeof(κ_face))
        end

        jac.coeff_L[f] = coeff
        jac.coeff_R[f] = -coeff
    end

    return nothing
end


# ============================================================================
# Pressure Boundary Condition Application
# ============================================================================

using ..BoundaryConditions: DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC
using ..Materials: AbstractMassBC
using ..BoundaryConditions: ImpermeableBC
# ============================================================================
# Type-Stable BoundaryConditionSet Signatures
# ============================================================================

using ..BoundaryConditions: BoundaryConditionSet, get_mass_bc

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


using ..BoundaryConditions: CombinedMassBC

"""
	boundary_darcy_velocity(bc, P_int, κ_int, μ_int, d, ::Val{tag}) -> v_z

Boundary-face Darcy velocity [m/s] in the +z convention, from a pressure BC
and the interior cell's state. This is the single shared kernel for the
residual (`apply_pressure_bc_to_darcy_velocity!`), the Jacobian
(`BCMassOperator`), and the saved-output diagnostics — keep them in lockstep
by changing only this function. AD-clean (arguments may be `Dual`s).

`tag` is `:top` (outward normal +z) or `:bottom` (outward normal −z).
Non-pressure mass BCs contribute 0; `CombinedMassBC` sums its children.

Methods:
- `DirichletPressureBC`: one-sided gradient v_z = −(κ/μ)·(P_bc − P_int)/d at
  the top (mirrored at the bottom); 0 for an impermeable interior cell.
- `NeumannPressureBC`: `∂P∂n` is the pressure gradient along the **outward
  normal**, so v_z = −(κ/μ)·∂P/∂z = −(κ/μ)·∂P∂n at the top and
  **+**(κ/μ)·∂P∂n at the bottom (n̂ = −ẑ there).
- `ConvectivePressureBC`: film law v_out = h_P·(P_int − P_ext) along the
  outward normal, gated on κ_int > 0 — an impermeable surface cell cannot
  vent regardless of h_P.
"""
function boundary_darcy_velocity end

@inline function boundary_darcy_velocity(
    bc::DirichletPressureBC,
    P_int::Real,
    κ_int::Real,
    μ_int::Real,
    d::Real,
    ::Val{tag},
) where {tag}
    Tv = promote_type(typeof(P_int), typeof(κ_int), typeof(μ_int), typeof(d))
    (d < 1e-15 || !(κ_int > 0) || !(μ_int > 1e-20)) && return zero(Tv)
    ∇P = tag === :top ? (bc.P - P_int) / d : (P_int - bc.P) / d
    return -κ_int / μ_int * ∇P
end

@inline function boundary_darcy_velocity(
    bc::NeumannPressureBC,
    P_int::Real,
    κ_int::Real,
    μ_int::Real,
    d::Real,
    ::Val{tag},
) where {tag}
    Tv = promote_type(typeof(P_int), typeof(κ_int), typeof(μ_int), typeof(d))
    (abs(bc.∂P∂n) < 1e-20 || !(κ_int > 0) || !(μ_int > 1e-20)) && return zero(Tv)
    # ∂P∂n is along the outward normal: ∂P/∂z = +∂P∂n at the top (n̂ = +ẑ),
    # −∂P∂n at the bottom (n̂ = −ẑ); v_z = −(κ/μ)·∂P/∂z.
    return tag === :top ? -κ_int / μ_int * bc.∂P∂n : κ_int / μ_int * bc.∂P∂n
end

@inline function boundary_darcy_velocity(
    bc::ConvectivePressureBC,
    P_int::Real,
    κ_int::Real,
    μ_int::Real,
    d::Real,
    ::Val{tag},
) where {tag}
    Tv = promote_type(typeof(P_int), typeof(κ_int), typeof(μ_int), typeof(d))
    κ_int > 0 || return zero(Tv)
    v_out = bc.h_P * (P_int - bc.P_ext)
    return tag === :top ? v_out : -v_out
end

# Non-pressure mass BCs contribute no boundary Darcy velocity.
@inline function boundary_darcy_velocity(
    ::AbstractMassBC,
    P_int::Real,
    κ_int::Real,
    μ_int::Real,
    d::Real,
    ::Val{tag},
) where {tag}
    return zero(promote_type(typeof(P_int), typeof(κ_int), typeof(μ_int), typeof(d)))
end

@generated function boundary_darcy_velocity(
    bc::CombinedMassBC{BCs},
    P_int::Real,
    κ_int::Real,
    μ_int::Real,
    d::Real,
    boundary_val::Val{tag},
) where {BCs<:Tuple,tag}
    N = length(BCs.parameters)
    exprs = [
        :(boundary_darcy_velocity(bc.bcs[$i], P_int, κ_int, μ_int, d, boundary_val)) for
        i = 1:N
    ]
    return :(+($(exprs...)))
end

# Per-boundary face loop: write the kernel's velocity into v_gas. Non-pressure
# BCs write 0 — same as the compute_darcy_velocity! default.
function _apply_pressure_bc_faces!(
    v_gas,
    state,
    mesh,
    P,
    κ,
    face_ids::Vector{Int},
    bc::AbstractMassBC,
    ::Val{tag},
) where {tag}
    @inbounds for f in face_ids
        face = mesh.faces[f]
        cell_id = face.cell_ids[1] != 0 ? face.cell_ids[1] : face.cell_ids[2]

        d = abs(state.face_positions[f] - state.cell_centers[cell_id])
        μ_int = effective_gas_viscosity(state.T[cell_id])

        v_gas[f] = boundary_darcy_velocity(bc, P[cell_id], κ[cell_id], μ_int, d, Val(tag))
    end
    return nothing
end

function apply_pressure_bc_to_darcy_velocity!(
    v_gas::AbstractVector,
    state,
    mesh,
    P::AbstractVector,
    κ::AbstractVector,
    bc_set::BoundaryConditionSet,
    top_face_ids::Vector{Int},
    bottom_face_ids::Vector{Int},
)
    mass_top = get_mass_bc(bc_set, Val(:top))
    _apply_pressure_bc_faces!(v_gas, state, mesh, P, κ, top_face_ids, mass_top, Val(:top))

    mass_bottom = get_mass_bc(bc_set, Val(:bottom))
    _apply_pressure_bc_faces!(
        v_gas,
        state,
        mesh,
        P,
        κ,
        bottom_face_ids,
        mass_bottom,
        Val(:bottom),
    )

    return nothing
end
