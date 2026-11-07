"""
Boundary temperature solver for cell-centered finite volume method.

This module implements Newton iteration to solve for the surface temperature T_s
that satisfies the surface energy balance at boundaries:

	q_BC(T_s) = q_cond(T_s, T_cell)

where:
- q_BC = q_ext - εσT_s⁴ - h(T_s - T_∞)  [applied heat flux, positive into domain]
- q_cond = k(T_s - T_cell) / Δz_half    [conduction into material, positive when T_s > T_cell]

This wall-heat balance adds no separate outgoing-gas enthalpy sink. The
boundary assembler combines signed donor enthalpy with heat at the same
control surface. For an LTE film it converts wall heat using
`q_r = q_s + H_r - H_s`; imposed heat and adiabatic conditions use incoming
reservoir enthalpy directly. Gas flux affects wall heat through the blowing
model. Coupled species/thermal roots are assembled in `Kernels/gas_boundary.jl`.
"""

using ..Materials: STEFAN_BOLTZMANN, T_PROPERTY_MIN, T_PROPERTY_MAX
using ..Materials: evaluate_param, evaluate_param_T
using ..Materials: RadiationModel, SURFACE_ABSORPTION, RADIATIVE_TRANSFER
using ForwardDiff
using StaticArrays: SVector
using ..Physics: primal, bounded_transport_root, implicit_transport_root, TransportFailure

# ============================================================================
# Main Solver
# ============================================================================

"""
	solve_boundary_temperature!(state::BoundaryState, T_cell, k_eff, Δz_half,
		bc, material, J_outward, t, ε_eff; kwargs...)

In-place Newton solve of the surface energy balance

	q_BC(T_s) = k_eff · (T_s − T_cell)/Δz_half

The condition is evaluated at `trial(T_s)`, the smooth clamp of the iterate
into the property bracket, so the balance has a root for any cell temperature
and the iterate itself is never clamped. The Newton starts from `T_cell`, so the result is a pure function of the
arguments. `state` is output only: it receives the converged surface
temperature and fluxes for diagnostics, and is left untouched when the final
iterate is not finite.

# Arguments
- `state::BoundaryState`: Diagnostic boundary state, written in place.
- `T_cell::Float64`: Adjacent cell-center temperature [K].
- `k_eff::Float64`: Effective conductivity at the boundary [W/(m·K)].
- `Δz_half::Float64`: Distance from boundary face to cell center [m].
- `bc::AbstractThermalBC`: Thermal BC.
- `material::Material{NC}`: Material (component cp_g, is_gas predicate).
- `J_outward::NTuple{NC,Float64}`: Per-component outward gas mass flux at this
  boundary face [kg/(m²·s)]. Outward is positive (gas leaving the slab).
- `t::Float64`: Current time [s].
- `ε_eff::Float64`: Effective surface emissivity (surface emission).
- `a_eff::Float64 = ε_eff` (keyword): Effective surface absorptivity for
  incident radiation, consumed by a `RadiativeFluxBC` with
  `absorptivity = nothing`. The default reproduces Kirchhoff's law.
- `trial` (keyword, required): the smooth clamp `T -> T_c` of a temperature
  into the property bracket at which the condition is evaluated.

# Blowing

Every `ConvectiveBC` inside `bc` applies its blowing model (`bc.blowing`, see
[`blowing_state`](@ref)) to a [`BlowingFace`](@ref) rebuilt at every Newton
iterate from `J_outward` and the current `T_s`.

# Non-convergence
The Newton works on the balance divided by its local scale, the largest of
the applied and conducted fluxes at the iterate. Where conduction is weak
against the boundary flux (a small `k_eff/Δz_half`), that scale is the flux
itself on either side of the root, the scaled residual stays near ±1 away
from it, and no step decreases it; the solve then converges first on the
balance divided by its scale at `T_cell` and finishes on the local scale,
raised to the round-off of the applied flux, a difference of larger absorbed
and emitted fluxes near the root. A balance on which neither solve converges
raises `TransportFailure`; no flux from an unconverged iterate is assembled.
Derivatives follow the converged residual through the implicit function
theorem.
"""
@inline function solve_boundary_temperature!(
    state::BoundaryState,
    T_cell::Real,
    k_eff::Real,
    Δz_half::Real,
    bc::AbstractThermalBC,
    material::Material{NC},
    J_outward::NTuple{NC,<:Real},
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    tol::Float64 = 1e-10,
    max_iter::Int = 50,
    trial,
) where {NC}
    # Solve a primal scalar balance, then differentiate the converged equation.
    # The balance at an iterate, its local scale, which includes the
    # resolvable conduction difference at T_s, and the resolvable difference
    # of the applied flux at T_s.
    balance =
        (z, T_cell, k_eff, Δz_half, ε_eff, a_eff, J_outward, t) -> begin
            temp_local=z[1]
            eval_local=trial(temp_local)
            face_local=BlowingFace(J_outward, eval_local, material)
            heat_local, slope_local=evaluate_bc_flux_and_derivative(
                bc,
                eval_local,
                t,
                ε_eff,
                radiation_model,
                face_local;
                a_eff,
            )
            conducted_local=k_eff*(temp_local-T_cell)/Δz_half
            scale_local=max(
                abs(primal(heat_local)),
                abs(primal(conducted_local)),
                1.0,
                64eps(Float64)*max(abs(primal(temp_local)), 1.0)*abs(primal(k_eff/Δz_half))/tol,
            )
            resolution_local=64eps(Float64)*max(abs(primal(temp_local)), 1.0)*abs(
                primal(slope_local),
            )/tol
            (heat_local-conducted_local, scale_local, resolution_local)
        end
    # The balance divided by its local scale, raised under weak conduction to
    # the applied flux's resolution.
    residual =
        (z, weak, T_cell, k_eff, Δz_half, ε_eff, a_eff, J_outward, t) -> begin
            b, s, r = balance(z, T_cell, k_eff, Δz_half, ε_eff, a_eff, J_outward, t)
            SVector(b/(weak ? max(s, r) : s))
        end
    plain = (
        primal(T_cell),
        primal(k_eff),
        primal(Δz_half),
        primal(ε_eff),
        primal(a_eff),
        map(primal, J_outward),
        primal(t),
    )
    z0=SVector(primal(T_cell))
    attempt=try
        root=bounded_transport_root(
            z->residual(z, false, plain...),
            z0;
            tolerance = tol,
            maxiters = max_iter,
            lower = SVector(-Inf),
        )
        (root = root, weak = false)
    catch failure
        failure isa TransportFailure || rethrow()
        (
            root = _weak_conduction_root(balance, residual, plain, z0, tol, max_iter),
            weak = true,
        )
    end
    root, weak = attempt.root, attempt.weak
    z=implicit_transport_root(
        z->residual(z, weak, T_cell, k_eff, Δz_half, ε_eff, a_eff, J_outward, t),
        root,
    )
    T_s=z[1]
    Te=trial(T_s)
    face=BlowingFace(J_outward, Te, material)
    q, _=evaluate_bc_flux_and_derivative(bc, Te, t, ε_eff, radiation_model, face; a_eff)
    qc=k_eff*(T_s-T_cell)/Δz_half
    if !(T_s isa ForwardDiff.Dual)
        state.T_surface=T_s
        state.q_applied=q
        state.q_conducted=qc
        state.T_cell=T_cell
        state.k_eff=k_eff
        state.Δz_half=Δz_half
        state.converged=true
        state.iterations=root.iterations
        _record_blowing!(state, bc, face, t)
    end
    return (T_s = T_s, q_applied = q, q_conducted = qc)
end

# The root of a balance whose conduction is weak against its boundary flux.
# The local scale is then the flux itself on either side of the root, so the
# scaled residual stays near ±1 away from it and no Newton step decreases it;
# and near the root the applied flux, a difference of larger absorbed and
# emitted ones, carries more round-off than the local scale resolves. Converge
# first on the balance divided once by its scale at the cell temperature, then
# finish on the local scale raised to the applied flux's resolution.
function _weak_conduction_root(balance, residual, plain, z0, tol, max_iter)
    _, s, r = balance(z0, plain...)
    scale = max(s, r)
    near = bounded_transport_root(
        z->SVector(first(balance(z, plain...))/scale),
        z0;
        tolerance = tol,
        maxiters = max_iter,
        lower = SVector(-Inf),
    )
    root = bounded_transport_root(
        z->residual(z, true, plain...),
        near.x;
        tolerance = tol,
        maxiters = max_iter,
        lower = SVector(-Inf),
    )
    return (
        x = root.x,
        linear = root.linear,
        iterations = near.iterations + root.iterations,
        residual = root.residual,
    )
end

# ============================================================================
# BC Flux and Derivative Evaluation
# ============================================================================

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff; a_eff = ε_eff)

Evaluate boundary heat flux and a frozen-coefficient temperature slope.
The slope holds heat-transfer, blowing and emissivity coefficients fixed.
Surface roots and their implicit derivatives differentiate the complete flux
residual with automatic differentiation; they do not use this approximate slope.

# Arguments
- `bc::AbstractThermalBC`: Boundary condition
- `T_s::Float64`: Surface temperature [K]
- `t::Float64`: Current time [s]
- `ε_eff::Float64`: Effective emissivity from material properties (emission)
- `a_eff::Float64 = ε_eff` (keyword): Effective surface absorptivity for
  incident radiation; only a `RadiativeFluxBC` with `absorptivity = nothing`
  reads it. Every overload accepts the keyword and forwards it.

# Returns
- `(q, dq_dT)`: Flux [W/m²] and frozen-coefficient slope [W/(m²·K)], positive into domain
"""
function evaluate_bc_flux_and_derivative end

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)

Radiation-model aware variant used by the boundary temperature solve.

It keeps `RadiativeFluxBC` consistent with the residual: under
`RADIATIVE_TRANSFER` the incident radiation is absorbed in depth and
contributes **zero** boundary flux here; under `SURFACE_ABSORPTION` it is
absorbed at the surface, at the surface absorptivity. Under
`RADIATIVE_TRANSFER` a `RadiativeBC` contributes zero as well: its exchange
with the surroundings crosses the face through the radiation rows.
"""
@inline function evaluate_bc_flux_and_derivative(
    bc::AbstractThermalBC,
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel;
    a_eff::Real = ε_eff,
)
    return evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff; a_eff = a_eff)
end

# Under `RADIATIVE_TRANSFER` the surroundings' radiation and the emission cross
# the exposed face through the radiation rows, so a `RadiativeBC` applies no
# flux at the face.
@inline function evaluate_bc_flux_and_derivative(
    bc::RadiativeBC,
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel;
    a_eff::Real = ε_eff,
)
    q, dq_dT = evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff; a_eff = a_eff)
    radiation_model == RADIATIVE_TRANSFER && return zero(q), zero(dq_dT)
    return q, dq_dT
end

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model, face)

Blowing-aware variant. `face` is the [`BlowingFace`](@ref) at this boundary
(outward species mass fluxes, surface-temperature iterate, material). Only
convective BCs consume it — each `ConvectiveBC` asks its own blowing model for
the ratio `χ = h_eff/h` against its own `h` and ambient temperature
([`blowing_ratio`](@ref)); every other BC type ignores `face` and falls
through to the radiation-model-aware method. Zero species fluxes reproduce the uncorrected flux exactly for every shipped model.
"""
@inline function evaluate_bc_flux_and_derivative(
    bc::AbstractThermalBC,
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel,
    face::BlowingFace;
    a_eff::Real = ε_eff,
)
    return evaluate_bc_flux_and_derivative(
        bc,
        T_s,
        t,
        ε_eff,
        radiation_model;
        a_eff = a_eff,
    )
end

# ConvectiveBC with blowing: q = h·χ·(T_inf − T_s), χ from the BC's model, with
# the frozen-χ slope dq/dT_s = −h·χ.
@inline function evaluate_bc_flux_and_derivative(
    bc::ConvectiveBC,
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel,
    face::BlowingFace;
    a_eff::Real = ε_eff,
)
    h = evaluate_param(bc.h, t)
    T_inf = evaluate_param(bc.T_inf, t)
    # h ≤ 0 leaves the blowing parameter undefined; the flux is then
    # ≤ 0·(T_inf − T_s) anyway, so apply no correction rather than divide by zero.
    χ = h > zero(h) ? blowing_ratio(bc.blowing, face, h, T_inf) : one(T_s)
    q = h * χ * (T_inf - T_s)
    dq_dT = -h * χ
    return q, dq_dT
end

@generated function evaluate_bc_flux_and_derivative(
    bc::CombinedThermalBC{BCs},
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel,
    face::BlowingFace;
    a_eff::Real = ε_eff,
) where {BCs}
    N = fieldcount(BCs)
    exprs = Expr[]
    push!(
        exprs,
        :(
            Tu = promote_type(
                typeof(T_s),
                typeof(t),
                typeof(ε_eff),
                typeof(a_eff),
                eltype(face.ṁ),
            )
        ),
    )
    push!(exprs, :(q_total = zero(Tu)))
    push!(exprs, :(dq_total = zero(Tu)))
    for i = 1:N
        push!(
            exprs,
            quote
                q_i, dq_i = evaluate_bc_flux_and_derivative(
                    bc.bcs[$i],
                    T_s,
                    t,
                    ε_eff,
                    radiation_model,
                    face;
                    a_eff = a_eff,
                )
                q_total += q_i
                dq_total += dq_i
            end,
        )
    end
    push!(exprs, :(return q_total, dq_total))
    return Expr(:block, exprs...)
end

# ----------------------------------------------------------------------------
# AdiabaticBC: q = 0
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    ::AdiabaticBC,
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    return 0.0, 0.0
end

# ----------------------------------------------------------------------------
# HeatFluxBC: q = prescribed flux (no T_s dependence)
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    bc::HeatFluxBC{<:Real},
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    return Float64(bc.flux), 0.0
end

@inline function evaluate_bc_flux_and_derivative(
    bc::HeatFluxBC{<:Function},
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    return bc.flux(t), 0.0
end

# ----------------------------------------------------------------------------
# RadiativeBC: q = ε × σ × (F × T_inf⁴ - T_s⁴). The view factor scales the
# irradiation from the surroundings only; the emission leaves into the whole
# hemisphere. Emission always uses ε_eff; `a_eff` is accepted for signature
# uniformity.
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    bc::RadiativeBC,
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    # Use effective emissivity if not specified
    ε = if bc.ε === nothing
        ε_eff
    else
        evaluate_param_T(bc.ε, T_s)
    end
    T_inf = evaluate_param(bc.T_inf, t)
    F = evaluate_param(bc.F, t)

    q = ε * STEFAN_BOLTZMANN * (F * T_inf^4 - T_s^4)
    dq_dT = -4 * ε * STEFAN_BOLTZMANN * T_s^3

    # This returned slope freezes ε. Surface solves differentiate q directly,
    # including dε_dT*σ*(F*T_inf^4-T_s^4) when ε depends on temperature.

    return q, dq_dT
end

# ----------------------------------------------------------------------------
# RadiativeFluxBC: q = α × incident_flux at the surface, or zero where the
# incident radiation is absorbed in depth. With `absorptivity = nothing` the BC
# absorbs at the material's mixture surface absorptivity `a_eff` (= ε_eff
# unless a component overrides it).
# - SURFACE_ABSORPTION: treat as boundary heat flux for correct T_surface
# - RADIATIVE_TRANSFER: handled volumetrically, return 0 to avoid
#   double-counting
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    bc::RadiativeFluxBC,
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel;
    a_eff::Real = ε_eff,
)
    if radiation_model == SURFACE_ABSORPTION
        # Radiation absorbed at boundary - include in boundary balance
        # This ensures T_surface is physically correct (hotter than T_cell when net heating)
        flux = evaluate_param(bc.flux, t)
        α = if bc.absorptivity === nothing
            a_eff
        else
            evaluate_param(bc.absorptivity, t)
        end
        return flux * α, 0.0
    end

    # In-depth absorption handles the incident radiation volumetrically, so
    # it must not enter the boundary balance.
    return 0.0, 0.0
end

# ----------------------------------------------------------------------------
# DirichletThermalBC: T = T_bc, as a stiff penalty flux toward T_bc
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    bc::DirichletThermalBC,
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    T_bc = evaluate_param(bc.T, t)

    # Use a very large heat transfer coefficient to enforce T_s ≈ T_bc
    # But not so large that it causes numerical issues
    h_large = 1e8  # W/(m²·K); larger values hurt the Newton conditioning
    q = h_large * (T_bc - T_s)
    dq_dT = -h_large

    return q, dq_dT
end
