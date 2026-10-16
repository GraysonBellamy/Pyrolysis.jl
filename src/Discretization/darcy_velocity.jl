"""
Pressure boundary conditions for the Darcy gas velocity.

The interior face velocity `v_g = -(κ/μ) ∇P` is assembled in
`update_darcy_properties!`; this file supplies the boundary faces from the
pressure BCs (`apply_pressure_bc_to_darcy_velocity!`) and the shared
boundary-face kernel `boundary_darcy_velocity`.
"""

# ============================================================================
# Pressure Boundary Condition Application
# ============================================================================

using ..BoundaryConditions: DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC
using ..Materials: AbstractMassBC
using ..BoundaryConditions: ImpermeableBC
using ..BoundaryConditions: BoundaryConditionSet, get_mass_bc, CombinedMassBC

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
# BCs write 0.
function _apply_pressure_bc_faces!(
    v_gas,
    state,
    mesh,
    P,
    κ,
    μ,
    face_ids::Vector{Int},
    bc::AbstractMassBC,
    ::Val{tag},
) where {tag}
    @inbounds for f in face_ids
        face = mesh.faces[f]
        cell_id = face.cell_ids[1] != 0 ? face.cell_ids[1] : face.cell_ids[2]

        d = abs(state.face_positions[f] - state.cell_centers[cell_id])

        v_gas[f] =
            boundary_darcy_velocity(bc, P[cell_id], κ[cell_id], μ[cell_id], d, Val(tag))
    end
    return nothing
end

"""
	apply_pressure_bc_to_darcy_velocity!(v_gas, state, mesh, P, κ, μ, bc_set,
										 top_face_ids, bottom_face_ids)

Write the boundary-face Darcy velocities into `v_gas` from the pressure BCs
of `bc_set` (via [`boundary_darcy_velocity`](@ref)), over the face IDs cached
at problem setup. Type-stable and allocation-free.

`μ` is the per-cell pore-gas viscosity the caller already computed
(`darcy.μ`), not a value re-derived here: the boundary vent and the interior
faces must see the same viscosity, and with a material-supplied μ(T) curve
re-evaluating it would be both a second evaluation and a second chance to
disagree.
"""
function apply_pressure_bc_to_darcy_velocity!(
    v_gas::AbstractVector,
    state,
    mesh,
    P::AbstractVector,
    κ::AbstractVector,
    μ::AbstractVector,
    bc_set::BoundaryConditionSet,
    top_face_ids::Vector{Int},
    bottom_face_ids::Vector{Int},
)
    mass_top = get_mass_bc(bc_set, Val(:top))
    _apply_pressure_bc_faces!(
        v_gas,
        state,
        mesh,
        P,
        κ,
        μ,
        top_face_ids,
        mass_top,
        Val(:top),
    )

    mass_bottom = get_mass_bc(bc_set, Val(:bottom))
    _apply_pressure_bc_faces!(
        v_gas,
        state,
        mesh,
        P,
        κ,
        μ,
        bottom_face_ids,
        mass_bottom,
        Val(:bottom),
    )

    return nothing
end
