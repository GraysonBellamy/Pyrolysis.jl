"""
Boundary flux implementation.

This module provides flux evaluation functions for all boundary condition types:
- Thermal BCs: heat flux into domain [W/m²]
- Mass BCs: mass flux into domain [kg/(m²·s)]
"""

using ..Materials: STEFAN_BOLTZMANN
using ..Materials: evaluate_param, evaluate_param_T
using ..Materials:
    RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY

# ============================================================================
# Thermal Boundary Condition Flux Evaluation
# ============================================================================

"""
	thermal_bc_flux(bc::AbstractThermalBC, T_surface, t) -> Float64

Evaluate the heat flux from a thermal boundary condition.

# Arguments
- `bc`: The thermal boundary condition
- `T_surface`: Surface temperature [K]
- `t`: Current time [s]

# Returns
Heat flux [W/m²], positive into the domain.
"""
function thermal_bc_flux end

# ============================================================================
# Standard thermal_bc_flux with effective emissivity parameter
# ============================================================================
# All thermal BC flux functions now accept an optional ε_eff parameter.
# This is the effective surface emissivity computed from material properties.
# BCs that don't need it simply ignore the parameter.

# Adiabatic: zero flux
@inline thermal_bc_flux(::AdiabaticBC, T_surface, t, ε_eff) = 0.0

# Heat flux: prescribed value (constant or time-dependent)
@inline thermal_bc_flux(bc::HeatFluxBC{<:Real}, T_surface, t, ε_eff) = Float64(bc.flux)
@inline thermal_bc_flux(bc::HeatFluxBC{<:Function}, T_surface, t, ε_eff) = bc.flux(t)

# Convective: q" = h(T∞ - T_surface)
@inline function thermal_bc_flux(bc::ConvectiveBC, T_surface, t, ε_eff)
    h = evaluate_param(bc.h, t)
    T_inf = evaluate_param(bc.T_inf, t)
    return h * (T_inf - T_surface)
end

# Radiative: q" = F × ε × σ × (T∞⁴ - T_surface⁴)
# If bc.ε is nothing, use ε_eff from material properties; otherwise use specified value
@inline function thermal_bc_flux(bc::RadiativeBC, T_surface, t, ε_eff)
    # Use effective emissivity if not specified, otherwise use provided value
    ε = if bc.ε === nothing
        ε_eff
    else
        evaluate_param_T(bc.ε, T_surface)
    end
    T_inf = evaluate_param(bc.T_inf, t)
    F = evaluate_param(bc.F, t)  # View factor (default 1.0)
    return F * ε * STEFAN_BOLTZMANN * (T_inf^4 - T_surface^4)
end

# RadiativeFluxBC: External radiation source
# Behavior depends on radiation_model passed to solve():
# - NO_RADIATION/SURFACE_ABSORPTION: absorbed flux returned as boundary flux
#   This ensures the boundary temperature solve correctly accounts for radiation heating.
# - BEER_LAMBERT/P1: handled volumetrically, return 0 here
# If bc.absorptivity is nothing, use ε_eff (Kirchhoff's law: α = ε at equilibrium)
# NOTE: radiation_model is REQUIRED - no default behavior
@inline function thermal_bc_flux(
    bc::RadiativeFluxBC,
    T_surface,
    t,
    ε_eff,
    radiation_model::RadiationModel,
)
    if radiation_model == BEER_LAMBERT || radiation_model == P1_QUASI_STEADY
        # In-depth absorption handled by volumetric radiation solver
        # Don't add surface flux to avoid double-counting
        return 0.0
    else
        # NO_RADIATION or SURFACE_ABSORPTION: radiation absorbed at boundary
        # For SURFACE_ABSORPTION, this ensures physically correct T_surface
        # (surface hotter than cell when net heating)
        flux = bc.flux isa Function ? bc.flux(t) : Float64(bc.flux)
        α = if bc.absorptivity === nothing
            ε_eff
        else
            bc.absorptivity isa Function ? bc.absorptivity(t) : Float64(bc.absorptivity)
        end
        return flux * α  # Absorbed radiation at surface
    end
end

# Dirichlet: a prescribed-temperature BC has no closed-form flux at this call
# signature (the applied flux is the conducted flux k_eff·(T_bc − T_cell)/Δz_half,
# which needs the cell state). NaN remains the deliberate "not expressible
# here" sentinel — but every DIAGNOSTIC path now bypasses it: standalone
# Dirichlet dispatches to the dedicated `compute_boundary_heat_fluxes`
# method below, and the saved-solution extraction substitutes the solved
# `q_conducted` for `net`.
@inline thermal_bc_flux(::DirichletThermalBC, T_surface, t, ε_eff) = NaN
@inline thermal_bc_flux(::DirichletThermalBC, T_surface, t, ε_eff, ::RadiationModel) = NaN

# Combined: sum of all component fluxes using generated function for type stability
# NOTE: radiation_model is REQUIRED for CombinedThermalBC since it may contain RadiativeFluxBC
@generated function thermal_bc_flux(
    bc::CombinedThermalBC{BCs},
    T_surface,
    t,
    ε_eff,
    radiation_model::RadiationModel,
) where {BCs<:Tuple}
    N = length(BCs.parameters)
    exprs =
        [:(thermal_bc_flux(bc.bcs[$i], T_surface, t, ε_eff, radiation_model)) for i = 1:N]
    return :(+($(exprs...)))
end

# Fallback for BCs that don't care about radiation_model (just ignore it)
@inline thermal_bc_flux(::AdiabaticBC, T_surface, t, ε_eff, ::RadiationModel) = 0.0
@inline thermal_bc_flux(bc::HeatFluxBC{<:Real}, T_surface, t, ε_eff, ::RadiationModel) =
    Float64(bc.flux)
@inline thermal_bc_flux(bc::HeatFluxBC{<:Function}, T_surface, t, ε_eff, ::RadiationModel) =
    bc.flux(t)
@inline thermal_bc_flux(bc::ConvectiveBC, T_surface, t, ε_eff, ::RadiationModel) =
    thermal_bc_flux(bc, T_surface, t, ε_eff)
# RadiativeBC under P1 contributes zero: emission/ambient exchange is carried
# volumetrically (4σαT⁴ + Marshak BC on G) — see boundary_temperature.jl.
@inline function thermal_bc_flux(
    bc::RadiativeBC,
    T_surface,
    t,
    ε_eff,
    radiation_model::RadiationModel,
)
    radiation_model == P1_QUASI_STEADY && return 0.0
    return thermal_bc_flux(bc, T_surface, t, ε_eff)
end

# ============================================================================
# Mass Boundary Condition Flux Evaluation
# ============================================================================

"""
	mass_bc_flux(bc::AbstractMassBC, ξ_surface, t, j::Integer) -> Float64

Evaluate the mass flux from a mass boundary condition.

# Arguments
- `bc`: The mass boundary condition
- `ξ_surface`: Surface mass concentration [kg/m³]
- `t`: Current time [s]
- `j`: Component index of the gas species this flux applies to. Lets
  per-species ambients (`ConvectiveMassBC` with an `NTuple` ξ_inf) resolve
  the right ambient concentration.

# Returns
Mass flux [kg/(m²·s)], positive into the domain.
"""
function mass_bc_flux end

# Impermeable: zero flux
@inline mass_bc_flux(::ImpermeableBC, ξ_surface, t, j::Integer) = 0.0

# Prescribed mass flux
@inline mass_bc_flux(bc::MassFluxBC{<:Real}, ξ_surface, t, j::Integer) = Float64(bc.flux)
@inline mass_bc_flux(bc::MassFluxBC{<:Function}, ξ_surface, t, j::Integer) = bc.flux(t)

# Convective mass transfer: J = h_m(ξ∞ - ξ_surface)
@inline function mass_bc_flux(bc::ConvectiveMassBC, ξ_surface, t, j::Integer)
    h_m = evaluate_param(bc.h_m, t)
    ξ_inf = ambient_concentration(bc, t, j)
    return h_m * (ξ_inf - ξ_surface)
end

# Combined: sum of all component fluxes using generated function for type stability
@generated function mass_bc_flux(
    bc::CombinedMassBC{BCs},
    ξ_surface,
    t,
    j::Integer,
) where {BCs<:Tuple}
    N = length(BCs.parameters)
    exprs = [:(mass_bc_flux(bc.bcs[$i], ξ_surface, t, j)) for i = 1:N]
    return :(+($(exprs...)))
end

# ============================================================================
# Mass BC Flux with Diffusion Parameters (λρ, ρ, d)
# ============================================================================
# These overloads are used when the discretization provides the combined transport 
# coefficient λρ (= λ × ρ), density ρ, and distance d from cell center to boundary.
# This is consistent with the rigorous λρ formulation for gas transport.

"""
	mass_bc_flux(bc::AbstractMassBC, ξ_surface, λρ, ρ, d, t, j::Integer) -> Float64

Evaluate mass flux with diffusion parameters available.

# Arguments
- `bc`: The mass boundary condition
- `ξ_surface`: Surface mass concentration [kg/m³]
- `λρ`: Combined transport coefficient λ·ρ [kg/(m·s)]
- `ρ`: Mixture density at surface [kg/m³]
- `d`: Distance from cell center to boundary [m]
- `t`: Current time [s]
- `j`: Component index of the gas species this flux applies to (per-species
  ambients resolve through it)

# Returns
Mass flux [kg/(m²·s)], positive into the domain (negative = gas leaving).
"""

# Default fallback: ignore λρ, ρ, d for non-diffusive BCs
@inline mass_bc_flux(bc::ImpermeableBC, ξ_surface, λρ, ρ, d, t, j::Integer) = 0.0
@inline mass_bc_flux(bc::MassFluxBC{<:Real}, ξ_surface, λρ, ρ, d, t, j::Integer) =
    Float64(bc.flux)
@inline mass_bc_flux(bc::MassFluxBC{<:Function}, ξ_surface, λρ, ρ, d, t, j::Integer) =
    bc.flux(t)
@inline function mass_bc_flux(bc::ConvectiveMassBC, ξ_surface, λρ, ρ, d, t, j::Integer)
    h_m = evaluate_param(bc.h_m, t)
    ξ_inf = ambient_concentration(bc, t, j)
    return h_m * (ξ_inf - ξ_surface)
end

"""
	mass_bc_flux(bc::DiffusiveMassBC, ξ_surface, λρ, ρ, d, t, j::Integer) -> Float64

Diffusive mass BC with combined internal/external resistance model.

Using the rigorous λρ formulation where λρ = λ × ρ is the combined transport coefficient.

The flux is computed as:
	J = Y / (R_internal + R_external)

where Y = ξ/ρ is mass fraction and:
- R_internal = d / λρ = d / (λρ) (internal diffusion resistance)
- R_external = 1 / (h_m × ρ) (external mass transfer resistance)

When h_m = Inf (no external resistance):
	J = λρ × Y / d = λ × ξ / d

When h_m is finite:
	J = λ × h_m × ξ / (λ + h_m × d)

Note: We extract λ = λρ/ρ for the combined resistance formula since h_m
is defined in terms of concentration, not mass fraction.

# Sign Convention
Returns NEGATIVE flux for gas leaving the surface (consistent with the general
mass_bc_flux convention of "positive into domain"). When ξ_surface > 0 and gas
is diffusing out, the flux into the domain is negative.

This differs from `_compute_surface_flux()` in Physics/gas_transport.jl, which
returns POSITIVE values for gas escaping (used for ALE regression calculation).

The calling code in `_apply_mass_bc_to_flux!` (finite_volume.jl) handles the
face orientation by negating at the top boundary: J[f] = -J_bc. This results
in positive flux in the +z direction when gas escapes at the top surface,
which is correct for the divergence calculation.

Summary of sign conventions:
- mass_bc_flux():           Negative = gas escaping (flux-into-domain convention)
- _compute_surface_flux():  Positive = gas escaping (ALE regression convention)
- J[face] in finite volume: Positive = flux in +z direction
"""
@inline function mass_bc_flux(bc::DiffusiveMassBC, ξ_surface, λρ, ρ, d, t, j::Integer)
    h_m = evaluate_param(bc.h_m, t)

    if ρ <= 0.0
        return 0.0
    end

    # Extract λ for the combined resistance calculation. Do not floor ρ here:
    # callers pass λρ = λ*ρ, so exact division recovers λ for any positive
    # density, including tiny numerical edge cases.
    λ = λρ / ρ

    # Combined resistance model: J = -ξ / (d/λ + 1/h_m)
    # Negative sign: gas leaving = negative flux into domain
    # With h_m = Inf: J = -λ × ξ / d
    if isinf(h_m)
        # Pure internal diffusion (no external resistance)
        return -λ * ξ_surface / d
    else
        # Combined resistance: J = -λ × h_m × ξ / (λ + h_m × d)
        return -(λ * h_m * ξ_surface) / (λ + h_m * d)
    end
end

# Combined: sum of all component fluxes using generated function for type stability
@generated function mass_bc_flux(
    bc::CombinedMassBC{BCs},
    ξ_surface,
    λρ,
    ρ,
    d,
    t,
    j::Integer,
) where {BCs<:Tuple}
    N = length(BCs.parameters)
    exprs = [:(mass_bc_flux(bc.bcs[$i], ξ_surface, λρ, ρ, d, t, j)) for i = 1:N]
    return :(+($(exprs...)))
end

# ============================================================================
# Pressure BC Mass Flux (for CombinedMassBC compatibility)
# ============================================================================
# Pressure BCs (DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC)
# contribute no DIFFUSIVE film-law flux, so these methods return 0 and compose
# with DiffusiveMassBC inside CombinedMassBC. Their advective venting flux
# enters the species equations through the boundary-face Darcy velocity:
# `apply_pressure_bc_to_darcy_velocity!` (via `boundary_darcy_velocity`) sets
# v_b, and `compute_gas_fluxes_darcy_fick!` writes J = v_b·ξ_int/φ_int at the
# boundary face (outflow-upwinded against a zero-concentration ambient).

using ..BoundaryConditions: DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC

# Simple signature (no diffusion parameters)
@inline mass_bc_flux(::DirichletPressureBC, ξ_surface, t, j::Integer) = 0.0
@inline mass_bc_flux(::NeumannPressureBC, ξ_surface, t, j::Integer) = 0.0
@inline mass_bc_flux(::ConvectivePressureBC, ξ_surface, t, j::Integer) = 0.0

# Full signature with diffusion parameters (λρ, ρ, d)
@inline mass_bc_flux(::DirichletPressureBC, ξ_surface, λρ, ρ, d, t, j::Integer) = 0.0
@inline mass_bc_flux(::NeumannPressureBC, ξ_surface, λρ, ρ, d, t, j::Integer) = 0.0
@inline mass_bc_flux(::ConvectivePressureBC, ξ_surface, λρ, ρ, d, t, j::Integer) = 0.0

# ============================================================================
# Boundary Flux Decomposition for Diagnostics
# ============================================================================

"""
	compute_boundary_heat_fluxes(bc::AbstractThermalBC, T_surface, t, ε_eff,
	                             radiation_model)

Decompose the boundary heat flux into named components for diagnostics.

Returns a NamedTuple with:
- `absorbed`: incident radiation absorbed *at the boundary face* [W/m²]
- `reradiation`: re-radiation loss [W/m²] (positive = leaving)
- `convection`: convective loss [W/m²] (positive = leaving)
- `prescribed`: prescribed flux [W/m²]
- `net`: the flux the surface energy balance actually applies at the face
  [W/m²], positive into the domain

Convention: `net` is the **applied** face flux — per component it matches
`evaluate_bc_flux_and_derivative`, so the total satisfies the Newton balance
`net = q_cond + q_transp` at convergence. Under `BEER_LAMBERT` (and the
experimental P1) the incident radiation enters the domain *volumetrically*,
not through the face, so a `RadiativeFluxBC` contributes zero `absorbed` and
zero `net` here; the volumetric deposit is reported separately as
`extras.volumetric_absorption`. Under `NO_RADIATION` / `SURFACE_ABSORPTION`
the absorbed α·q_ext is applied at the face and included in `net`. For
combined BCs the components *and* `net` are summed over children — `net` is
never recomputed from the components.
"""
function compute_boundary_heat_fluxes(
    bc::RadiativeFluxBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel,
)
    if radiation_model == BEER_LAMBERT || radiation_model == P1_QUASI_STEADY
        # In-depth absorption: the incident radiation enters volumetrically
        # (compute_radiation_source! / the P1 cache), nothing is absorbed at
        # the face, and the Newton balance applies zero flux for this BC.
        # The volumetric deposit is reported via extras.volumetric_absorption.
        return (
            absorbed = 0.0,
            reradiation = 0.0,
            convection = 0.0,
            prescribed = 0.0,
            net = 0.0,
        )
    else
        # NO_RADIATION / SURFACE_ABSORPTION: α·q_ext is applied at the face by
        # the surface energy balance (evaluate_bc_flux_and_derivative).
        flux = bc.flux isa Function ? bc.flux(t) : Float64(bc.flux)
        α = if bc.absorptivity === nothing
            ε_eff
        else
            bc.absorptivity isa Function ? bc.absorptivity(t) : Float64(bc.absorptivity)
        end
        q_absorbed = flux * α
        return (
            absorbed = q_absorbed,
            reradiation = 0.0,
            convection = 0.0,
            prescribed = 0.0,
            net = q_absorbed,
        )
    end
end

function compute_boundary_heat_fluxes(
    bc::CombinedThermalBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel,
)
    # Sum contributions from all component BCs — including each child's `net`,
    # so the combined `net` stays equal to the flux the surface Newton applies
    # (recomputing it from the components would re-include radiation that the
    # model dispatch routes through the volume, e.g. α·q_ext under BEER_LAMBERT).
    absorbed = 0.0
    reradiation = 0.0
    convection = 0.0
    prescribed = 0.0
    net = 0.0

    for inner_bc in bc.bcs
        fluxes =
            compute_boundary_heat_fluxes(inner_bc, T_surface, t, ε_eff, radiation_model)
        absorbed += fluxes.absorbed
        reradiation += fluxes.reradiation
        convection += fluxes.convection
        prescribed += fluxes.prescribed
        net += fluxes.net
    end

    return (
        absorbed = absorbed,
        reradiation = reradiation,
        convection = convection,
        prescribed = prescribed,
        net = net,
    )
end

# Overloads with radiation_model for BCs that don't care about it
function compute_boundary_heat_fluxes(
    bc::AdiabaticBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel,
)
    return (
        absorbed = 0.0,
        reradiation = 0.0,
        convection = 0.0,
        prescribed = 0.0,
        net = 0.0,
    )
end

function compute_boundary_heat_fluxes(
    bc::HeatFluxBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel,
)
    q = thermal_bc_flux(bc, T_surface, t, ε_eff)
    return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = q, net = q)
end

function compute_boundary_heat_fluxes(
    bc::ConvectiveBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel,
)
    q = thermal_bc_flux(bc, T_surface, t, ε_eff)
    q_loss = -q
    return (
        absorbed = 0.0,
        reradiation = 0.0,
        convection = q_loss,
        prescribed = 0.0,
        net = q,
    )
end

function compute_boundary_heat_fluxes(
    bc::RadiativeBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel,
)
    # Under P1 the boundary balance applies zero for this BC (emission and
    # ambient exchange are volumetric) — `net` must match the applied flux.
    q = thermal_bc_flux(bc, T_surface, t, ε_eff, radiation_model)
    q_loss = -q
    return (
        absorbed = 0.0,
        reradiation = q_loss,
        convection = 0.0,
        prescribed = 0.0,
        net = q,
    )
end

function compute_boundary_heat_fluxes(
    bc::DirichletThermalBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel,
)
    # Prescribed temperature: the applied flux is whatever conduction demands
    # (k_eff·(T_bc − T_cell)/Δz_half), not computable from (T_surface, t)
    # alone. The decomposition is all-zero here; the saved-solution extraction
    # reports the solved q_conducted as the net applied flux (previously the
    # AbstractThermalBC fallback returned NaN and poisoned
    # every extras.surface_heat_flux_* field for Dirichlet runs).
    return (
        absorbed = 0.0,
        reradiation = 0.0,
        convection = 0.0,
        prescribed = 0.0,
        net = 0.0,
    )
end

function compute_boundary_heat_fluxes(
    bc::AbstractThermalBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel,
)
    # Fallback for other BC types
    q_total = thermal_bc_flux(bc, T_surface, t, ε_eff)
    return (
        absorbed = 0.0,
        reradiation = 0.0,
        convection = 0.0,
        prescribed = q_total,
        net = q_total,
    )
end

# ----------------------------------------------------------------------------
# Blowing-aware decomposition (trailing `mcp` argument)
# ----------------------------------------------------------------------------
# `mcp` is the net blowing conductance ṁ″·c_p,g at the converged T_s, read
# from `BoundaryState.blowing_mcp` (zero when the blowing correction is off).
# Mirrors `evaluate_bc_flux_and_derivative`: only convective children rescale
# their h by `blowing_factor(mcp/h)`; every other BC ignores `mcp`. This keeps
# the diagnostic decomposition consistent with the flux the surface Newton
# actually applied.

function compute_boundary_heat_fluxes(
    bc::AbstractThermalBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel,
    mcp::Float64,
)
    return compute_boundary_heat_fluxes(bc, T_surface, t, ε_eff, radiation_model)
end

function compute_boundary_heat_fluxes(
    bc::ConvectiveBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel,
    mcp::Float64,
)
    h = evaluate_param(bc.h, t)
    T_inf = evaluate_param(bc.T_inf, t)
    χ = h > 0.0 ? blowing_factor(mcp / h) : 1.0
    q = h * χ * (T_inf - T_surface)
    q_loss = -q
    return (
        absorbed = 0.0,
        reradiation = 0.0,
        convection = q_loss,
        prescribed = 0.0,
        net = q,
    )
end

function compute_boundary_heat_fluxes(
    bc::CombinedThermalBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel,
    mcp::Float64,
)
    absorbed = 0.0
    reradiation = 0.0
    convection = 0.0
    prescribed = 0.0
    net = 0.0

    for inner_bc in bc.bcs
        fluxes = compute_boundary_heat_fluxes(
            inner_bc,
            T_surface,
            t,
            ε_eff,
            radiation_model,
            mcp,
        )
        absorbed += fluxes.absorbed
        reradiation += fluxes.reradiation
        convection += fluxes.convection
        prescribed += fluxes.prescribed
        net += fluxes.net
    end

    return (
        absorbed = absorbed,
        reradiation = reradiation,
        convection = convection,
        prescribed = prescribed,
        net = net,
    )
end

"""
	compute_surface_mass_flux(bc::AbstractMassBC, ξ_surface, λρ, ρ, d, t, j::Integer) -> Float64

Compute the mass flux of gas species `j` leaving the surface (positive = leaving).

This differs from mass_bc_flux which uses flux-into-domain convention.
"""
function compute_surface_mass_flux(
    bc::AbstractMassBC,
    ξ_surface::Float64,
    λρ::Float64,
    ρ::Float64,
    d::Float64,
    t::Float64,
    j::Integer,
)
    # mass_bc_flux returns negative for gas escaping (flux into domain convention)
    # We want positive = gas leaving
    return -mass_bc_flux(bc, ξ_surface, λρ, ρ, d, t, j)
end

# Simple overload without diffusion parameters
function compute_surface_mass_flux(
    bc::AbstractMassBC,
    ξ_surface::Float64,
    t::Float64,
    j::Integer,
)
    return -mass_bc_flux(bc, ξ_surface, t, j)
end
