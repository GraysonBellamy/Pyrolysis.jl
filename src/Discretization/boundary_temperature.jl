"""
Boundary temperature solver for cell-centered finite volume method.

This module implements Newton iteration to solve for the surface temperature T_s
that satisfies the surface energy balance at boundaries:

	q_BC(T_s) = q_cond(T_s, T_cell)

where:
- q_BC = q_ext - εσT_s⁴ - h(T_s - T_∞)  [applied heat flux, positive into domain]
- q_cond = k(T_s - T_cell) / Δz_half    [conduction into material, positive when T_s > T_cell]

The balance has no gas-enthalpy term. The sensible enthalpy of the gas crossing
the slab is carried by the volumetric advection source `−Σ_g cp_g J_g · ∇T`
(`S_conv`, see `Physics.convection`), the ThermaKin Eq.17 / FDS / Gpyro
convention. A surface term `Σ_g ṁ_g · Δh_g(T_s, T_∞)` on top of `S_conv` would
count the same enthalpy twice. The outward gas flux enters only through the
blowing models of the convective BCs.
"""

using ..Materials: STEFAN_BOLTZMANN, T_PROPERTY_MIN, T_PROPERTY_MAX
using ..Materials: evaluate_param, evaluate_param_T
using ..Materials:
    RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY
using ForwardDiff

# ============================================================================
# Main Solver
# ============================================================================

"""
	solve_boundary_temperature(T_cell, k_eff, Δz_half, bc, material, J_outward, t, ε_eff; kwargs...)

Solve for surface temperature satisfying energy balance using Newton iteration.

Finds T_s such that:
	q_BC(T_s) = k_eff × (T_s - T_cell) / Δz_half

Both sides are positive when heat flows INTO the material.

# Arguments
- `T_cell::Float64`: Temperature at adjacent cell center [K]
- `k_eff::Float64`: Effective thermal conductivity [W/(m·K)]
- `Δz_half::Float64`: Distance from boundary to cell center [m]
- `bc::AbstractThermalBC`: Thermal boundary condition
- `material::Material{NC}`: Material (gas heat capacities for blowing)
- `J_outward::NTuple{NC,Float64}`: Per-component outward gas mass flux at the
  face [kg/(m²·s)], positive leaving the slab; read by the blowing models
- `t::Float64`: Current time [s]
- `ε_eff::Float64`: Effective surface emissivity from material properties

# Keyword Arguments
- `a_eff::Float64`: Effective surface absorptivity for incident radiation
  (default: `ε_eff`, Kirchhoff). Used by a `RadiativeFluxBC` without an
  explicit `absorptivity`; emission (`RadiativeBC`) always uses `ε_eff`.
- `tol::Float64`: Convergence tolerance (default: 1e-10)
- `max_iter::Int`: Maximum iterations (default: 50)
- `relaxation::Float64`: Under-relaxation factor (default: 1.0)

# Returns
- `BoundaryState`: Solved boundary state with T_surface, fluxes, convergence info
"""
function solve_boundary_temperature(
    T_cell::Float64,
    k_eff::Float64,
    Δz_half::Float64,
    bc::AbstractThermalBC,
    material::Material{NC},
    J_outward::NTuple{NC,Float64},
    t::Float64,
    ε_eff::Float64;
    a_eff::Float64 = ε_eff,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
    tol::Float64 = 1e-10,
    max_iter::Int = 50,
    relaxation::Float64 = 1.0,
) where {NC}
    state = BoundaryState(T_cell, 0.0, 0.0, T_cell, k_eff, Δz_half, false, 0)
    solve_boundary_temperature!(
        state,
        T_cell,
        k_eff,
        Δz_half,
        bc,
        material,
        J_outward,
        t,
        ε_eff;
        a_eff = a_eff,
        radiation_model = radiation_model,
        tol = tol,
        max_iter = max_iter,
        relaxation = relaxation,
    )
    return state
end

"""
	solve_boundary_temperature!(state::BoundaryState, T_cell, k_eff, Δz_half,
		bc, material, J_outward, t, ε_eff; kwargs...)

In-place Newton solve of the surface energy balance

	q_BC(T_s) = k_eff · (T_s − T_cell)/Δz_half

The Newton starts from `T_cell` clamped to the property bracket, so the result
is a pure function of the arguments. `state` is output only: it receives the
converged surface temperature and fluxes for diagnostics, and is left untouched
when the final iterate is not finite.

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

# Blowing

Every `ConvectiveBC` inside `bc` applies its blowing model (`bc.blowing`, see
[`blowing_state`](@ref)) to a [`BlowingFace`](@ref) rebuilt at every Newton
iterate from `J_outward` and the current `T_s`. The Newton `df` uses the
frozen-χ derivative `−h·χ`: Newton still converges to the exact residual root,
marginally sub-quadratically, and the Dual fixed point is the exact
implicit-function-theorem derivative.

# Non-convergence and clamping
On max-iteration exhaustion or a singular Newton derivative the solve does NOT
abort: it warns (throttled, `maxlog = 10`) and evaluates the boundary flux at
the last iterate — the integrator's step-error control is the backstop. The
iterate is clamped to the physical bracket [200, 3000] K; a *final* iterate
pinned at the clamp also warns, because the balance then has no root in the
bracket and (in AD passes) the clamp reconstruction has zeroed the Dual
partials of `T_s` (the mathematically correct — but silent — clamp
derivative). The `state.converged` / `state.iterations` diagnostics record
the outcome either way.
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
    relaxation::Float64 = 1.0,
) where {NC}
    # Newton iterates carry whatever scalar type the inputs do (Float64 in
    # the canonical RHS, ForwardDiff.Dual under AD). Storage on `state` stays
    # Float64 — we strip with ForwardDiff.value at the end. The returned
    # (T_s, q_bc, q_cond) are typed and propagate derivatives correctly.
    Tu = promote_type(
        typeof(T_cell),
        typeof(k_eff),
        typeof(Δz_half),
        typeof(t),
        typeof(ε_eff),
        typeof(a_eff),
        eltype(J_outward),
    )
    # The Newton starts from the cell temperature clamped to the property bracket,
    # never from `state`: the result is a pure function of the arguments, so a wild
    # trial state of the outer Newton cannot leak into later residual evaluations.
    T_s::Tu = Tu(clamp(ForwardDiff.value(T_cell), T_PROPERTY_MIN, T_PROPERTY_MAX))
    dual_mode = Tu <: ForwardDiff.Dual
    converged = false
    iter = 0
    f_val = 0.0
    f_floor = 0.0
    scale = 1.0

    for i = 1:max_iter
        iter = i

        # The face at the current iterate: every convective member applies its
        # blowing model to it (`NoBlowing` returns the exact uncorrected flux).
        face = BlowingFace(J_outward, T_s, material)
        q_bc, dq_bc_dT = evaluate_bc_flux_and_derivative(
            bc,
            T_s,
            t,
            ε_eff,
            radiation_model,
            face;
            a_eff = a_eff,
        )

        q_cond = k_eff * (T_s - T_cell) / Δz_half
        dq_cond_dT = k_eff / Δz_half

        f = q_bc - q_cond
        df = dq_bc_dT - dq_cond_dT

        # Convergence on value (Newton residual; Dual derivatives ride along).
        f_val = ForwardDiff.value(f)
        scale = max(abs(ForwardDiff.value(q_bc)), abs(ForwardDiff.value(q_cond)), 1.0)
        # Round-off floor: the flux change from a few ulps of T_s. The conduction term is a
        # difference of two ~T numbers, so the residual cannot resolve below this; a face
        # with near-zero net flux would otherwise stall at round-off and report failure.
        f_floor =
            32 *
            eps(Float64) *
            max(abs(ForwardDiff.value(T_s)), 1.0) *
            (abs(ForwardDiff.value(dq_cond_dT)) + abs(ForwardDiff.value(dq_bc_dT)))
        if abs(f_val) < max(tol * scale, f_floor) && !(dual_mode && i == 1)
            converged = true
            break
        end

        if abs(ForwardDiff.value(df)) < 1e-30
            break
        end

        ΔT = -f / df
        T_s = T_s + relaxation * ΔT
        T_s_val = ForwardDiff.value(T_s)
        if T_s_val < T_PROPERTY_MIN || T_s_val > T_PROPERTY_MAX
            # Reconstructing through Tu(...) zeroes any accumulated Dual
            # partials. Self-healing: the next un-clamped Newton step rebuilds
            # them from (T_cell, f, df) — unless the LAST iterate is clamped,
            # which the post-loop check below reports.
            T_s = Tu(clamp(T_s_val, T_PROPERTY_MIN, T_PROPERTY_MAX))
        end
    end

    # Never fall through silently. The final flux evaluation
    # below proceeds at the unconverged T_s either way (the integrator's own
    # error control will reject the step if the resulting RHS is bad), but the
    # user must be able to see it happening.
    if !converged
        @warn "Surface-temperature Newton did not converge; using the last \
               iterate for the boundary flux" T_surface = ForwardDiff.value(T_s) T_cell =
            ForwardDiff.value(T_cell) iterations = iter residual_W_m2 = f_val tol_W_m2 =
            tol * scale roundoff_floor_W_m2 = f_floor t = ForwardDiff.value(t) maxlog = 10
    end

    # A final iterate pinned at the property-bracket clamp means
    # the surface energy balance has no root in the physical bracket — and in
    # AD passes the clamp has zeroed the Dual partials of T_s (zero is the
    # mathematically correct clamp derivative, but the silent kink is worth a
    # warning).
    T_s_pinned = ForwardDiff.value(T_s)
    if T_s_pinned == T_PROPERTY_MIN || T_s_pinned == T_PROPERTY_MAX
        @warn "Surface temperature pinned at the physical clamp; the energy \
               balance has no root in [$(T_PROPERTY_MIN), $(T_PROPERTY_MAX)] K \
               (AD sensitivities of T_s are zero at the clamp)" T_surface = T_s_pinned T_cell =
            ForwardDiff.value(T_cell) maxlog = 10
    end

    # Final flux evaluation (typed) — these are what the caller will read.
    face_final = BlowingFace(J_outward, T_s, material)
    q_bc_final, _ = evaluate_bc_flux_and_derivative(
        bc,
        T_s,
        t,
        ε_eff,
        radiation_model,
        face_final;
        a_eff = a_eff,
    )
    q_cond_final = k_eff * (T_s - T_cell) / Δz_half

    # Float64 storage on the diagnostic BoundaryState, skipped for a non-finite
    # iterate so the diagnostics keep the last finite solve.
    isfinite(ForwardDiff.value(T_s)) ||
        return (T_s = T_s, q_applied = q_bc_final, q_conducted = q_cond_final)
    state.T_surface = ForwardDiff.value(T_s)
    state.q_applied = ForwardDiff.value(q_bc_final)
    state.q_conducted = ForwardDiff.value(q_cond_final)
    _record_blowing!(state, bc, face_final, t)
    state.T_cell = ForwardDiff.value(T_cell)
    state.k_eff = ForwardDiff.value(k_eff)
    state.Δz_half = ForwardDiff.value(Δz_half)
    state.converged = converged
    state.iterations = iter

    return (T_s = T_s, q_applied = q_bc_final, q_conducted = q_cond_final)
end

# ============================================================================
# BC Flux and Derivative Evaluation
# ============================================================================

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff; a_eff = ε_eff)

Evaluate boundary condition flux and its temperature derivative.

# Arguments
- `bc::AbstractThermalBC`: Boundary condition
- `T_s::Float64`: Surface temperature [K]
- `t::Float64`: Current time [s]
- `ε_eff::Float64`: Effective emissivity from material properties (emission)
- `a_eff::Float64 = ε_eff` (keyword): Effective surface absorptivity for
  incident radiation; only a `RadiativeFluxBC` with `absorptivity = nothing`
  reads it. Every overload accepts the keyword and forwards it.

# Returns
- `(q, dq_dT)`: Flux [W/m²] and derivative [W/(m²·K)], positive into domain
"""
function evaluate_bc_flux_and_derivative end

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model)

Radiation-model aware variant used by the boundary temperature solve.

This is required to keep `RadiativeFluxBC` consistent with the main RHS, where
incident radiation is typically handled *volumetrically* (e.g. SURFACE_ABSORPTION)
and should therefore contribute **zero** boundary flux here.
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

"""
	evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff, radiation_model, face)

Blowing-aware variant. `face` is the [`BlowingFace`](@ref) at this boundary
(outward species mass fluxes, surface-temperature iterate, material). Only
convective BCs consume it — each `ConvectiveBC` asks its own blowing model for
the ratio `χ = h_eff/h` against its own `h` and ambient temperature
([`blowing_ratio`](@ref)); every other BC type ignores `face` and falls
through to the radiation-model-aware 5-argument method. Zero outflow
reproduces the uncorrected flux exactly for every shipped model.
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

# ConvectiveBC with blowing: q = h·χ·(T_inf − T_s), χ from the BC's model.
# The Newton derivative uses the frozen-χ form dq/dT_s = −h·χ, dropping the
# h·(∂χ/∂T_s)·(T_inf − T_s) term. The derivative is first-order consistent:
# Newton converges to the exact residual root, marginally sub-quadratically;
# the Dual fixed point is still the exact implicit-function-theorem derivative.
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
# ConvectiveBC: q = h(T_inf - T_s)
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    bc::ConvectiveBC,
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    h = evaluate_param(bc.h, t)
    T_inf = evaluate_param(bc.T_inf, t)

    q = h * (T_inf - T_s)
    dq_dT = -h
    return q, dq_dT
end

# ----------------------------------------------------------------------------
# RadiativeBC: q = F × ε × σ × (T_inf⁴ - T_s⁴)
# Emission always uses ε_eff; `a_eff` is accepted for signature uniformity.
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

    q = F * ε * STEFAN_BOLTZMANN * (T_inf^4 - T_s^4)
    dq_dT = -4 * F * ε * STEFAN_BOLTZMANN * T_s^3

    # Note: If ε depends on T, there's an additional term:
    # dq_dT += F * dε_dT * σ * (T_inf^4 - T_s^4)
    # This is usually small and we neglect it for simplicity

    return q, dq_dT
end

# Radiation-model aware handling for surface emission / ambient exchange.
# Under P1 the medium emits volumetrically (4σαT⁴) and exchanges with the
# ambient through the Marshak BC on G — a surface RadiativeBC on top of that
# double-counts both channels, so it contributes zero to the boundary balance.
# Under every other model (incl. BEER_LAMBERT, which deposits only the
# *incident* flux volumetrically) the surface emission stays at the boundary.
@inline function evaluate_bc_flux_and_derivative(
    bc::RadiativeBC,
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel;
    a_eff::Real = ε_eff,
)
    if radiation_model == P1_QUASI_STEADY
        return 0.0, 0.0
    end
    return evaluate_bc_flux_and_derivative(bc, T_s, t, ε_eff)
end

# ----------------------------------------------------------------------------
# RadiativeFluxBC: q = α × incident_flux (for NO_RADIATION mode)
# For in-depth radiation models, this is handled volumetrically.
# With `absorptivity = nothing` the BC absorbs at the material's mixture
# surface absorptivity `a_eff` (= ε_eff unless a component overrides it).
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    bc::RadiativeFluxBC,
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    # This BC represents incident radiation - absorptivity doesn't depend on T_s
    # If radiation is handled volumetrically, flux should be 0 here
    # For now, we return the absorbed flux (consistent with NO_RADIATION mode)
    flux = evaluate_param(bc.flux, t)
    α = if bc.absorptivity === nothing
        a_eff
    else
        evaluate_param(bc.absorptivity, t)
    end
    return flux * α, 0.0
end

# Radiation-model aware handling for incident radiation.
# - NO_RADIATION/SURFACE_ABSORPTION: treat as boundary heat flux for correct T_surface
# - BEER_LAMBERT/P1: handled volumetrically, return 0 to avoid double-counting
@inline function evaluate_bc_flux_and_derivative(
    bc::RadiativeFluxBC,
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel;
    a_eff::Real = ε_eff,
)
    if radiation_model == NO_RADIATION || radiation_model == SURFACE_ABSORPTION
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

    # Volumetric radiation models (BEER_LAMBERT, P1) handle incident
    # radiation separately, so we must not include it in the boundary balance.
    return 0.0, 0.0
end

# ----------------------------------------------------------------------------
# DirichletThermalBC: T = T_bc (special handling, not typical flux BC)
# For boundary temperature solve, we return the prescribed temperature directly
# since T_s = T_bc by definition
# ----------------------------------------------------------------------------
@inline function evaluate_bc_flux_and_derivative(
    bc::DirichletThermalBC,
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
)
    # Dirichlet BC prescribes temperature directly
    # The "flux" formulation: infinite conductance enforces T_s = T_bc
    # We handle this by setting large h with T_inf = T_bc
    T_bc = evaluate_param(bc.T, t)

    # Use a very large heat transfer coefficient to enforce T_s ≈ T_bc
    # But not so large that it causes numerical issues
    h_large = 1e8  # W/(m²·K) - reduced from 1e10 for stability
    q = h_large * (T_bc - T_s)
    dq_dT = -h_large

    return q, dq_dT
end

# ----------------------------------------------------------------------------
# CombinedThermalBC: Sum of component fluxes
# ----------------------------------------------------------------------------
# Use regular function instead of @generated to avoid excessive compilation time
@generated function evaluate_bc_flux_and_derivative(
    bc::CombinedThermalBC{BCs},
    T_s::Real,
    t::Real,
    ε_eff::Real;
    a_eff::Real = ε_eff,
) where {BCs}
    N = fieldcount(BCs)
    exprs = Expr[]
    push!(exprs, :(q_total = 0.0))
    push!(exprs, :(dq_total = 0.0))
    for i = 1:N
        push!(
            exprs,
            quote
                q_i, dq_i = evaluate_bc_flux_and_derivative(
                    bc.bcs[$i],
                    T_s,
                    t,
                    ε_eff;
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

@generated function evaluate_bc_flux_and_derivative(
    bc::CombinedThermalBC{BCs},
    T_s::Real,
    t::Real,
    ε_eff::Real,
    radiation_model::RadiationModel;
    a_eff::Real = ε_eff,
) where {BCs}
    N = fieldcount(BCs)
    exprs = Expr[]
    push!(exprs, :(Tu = promote_type(typeof(T_s), typeof(t), typeof(ε_eff), typeof(a_eff))))
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
                    radiation_model;
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

# ============================================================================
# Boundary State Update Function
# ============================================================================

"""
	update_boundary_state!(state_cache, mesh, boundary_id, bc, cache, material,
	                       J_outward, t, radiation_model)

Solve for `T_s` satisfying the surface energy balance
`q_BC(T_s) = k_eff·(T_s − T_cell)/Δz_half` at boundary `boundary_id`
(1 = top, 2 = bottom). `J_outward[j]` is the per-component outward gas mass
flux at the boundary face (positive = leaving the slab), read by the blowing
models of the convective BCs. Auxiliary callers without the live `J_gas` pass
zeros, the no-blowing limit.
"""
function update_boundary_state!(
    state_cache,
    mesh::Mesh1D,
    boundary_id::Int,
    bc::AbstractThermalBC,
    cache::PropertyCache,
    material::Material{NC},
    J_outward::NTuple{NC,<:Real},
    t::Real,
    radiation_model::RadiationModel = SURFACE_ABSORPTION,
) where {NC}
    # Cell adjacent to this boundary (1 = top/exposed at z = L; 2 = bottom).
    if boundary_id == 1
        cell_id = mesh.n_active_cells
        face_id = mesh.cells[cell_id].face_ids[2]
    else
        cell_id = 1
        face_id = mesh.cells[cell_id].face_ids[1]
    end

    T_cell = state_cache.T[cell_id]
    k_eff = cache.k[cell_id]
    ε_eff = cache.ε[cell_id]
    a_eff = cache.a[cell_id]
    Δz_half = state_cache.volume[cell_id] / state_cache.A_face[face_id] / 2

    state = mesh.boundary_states[boundary_id]
    # `state` storage is Float64 (diagnostics); the returned tuple carries
    # Tu typing so AD derivatives flow through `q_conducted` into `q[f]`.
    return solve_boundary_temperature!(
        state,
        T_cell,
        k_eff,
        Δz_half,
        bc,
        material,
        J_outward,
        t,
        ε_eff;
        a_eff = a_eff,
        radiation_model = radiation_model,
    )
end

# ============================================================================
# Diagnostic Functions
# ============================================================================

"""
	check_boundary_energy_balance(mesh)

Verify that boundary energy balance is satisfied.

# Returns
Dict{Int, Float64}: Boundary ID to relative imbalance (should be ~machine precision).
"""
function check_boundary_energy_balance(mesh::AbstractMesh)
    results = Dict{Int,Float64}()

    for (boundary_id, state) in mesh.boundary_states
        if abs(state.q_applied) > 1e-10
            imbalance = abs(state.q_applied - state.q_conducted) / abs(state.q_applied)
        else
            imbalance = abs(state.q_applied - state.q_conducted)
        end
        results[boundary_id] = imbalance
    end

    return results
end

"""
	get_surface_temperature(mesh)

Get the top surface temperature from the solved boundary state.
"""
function get_surface_temperature(mesh::AbstractMesh)
    return mesh.boundary_states[1].T_surface
end

"""
	get_back_face_temperature(mesh)

Get the back face (bottom) temperature from the solved boundary state.
"""
function get_back_face_temperature(mesh::AbstractMesh)
    return mesh.boundary_states[2].T_surface
end
