"""
Boundary flux implementation.

This module provides flux evaluation functions for all boundary condition types:
- Thermal BCs: heat flux into domain [W/m²]
- Mass BCs: mass flux into domain [kg/(m²·s)]
"""

using ..Materials: STEFAN_BOLTZMANN
using ..Materials: evaluate_param, evaluate_param_T
using ..Materials: RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT

# ============================================================================
# Boundary Flux Decomposition for Diagnostics
# ============================================================================

"""
	compute_boundary_heat_fluxes(bc::AbstractThermalBC, T_surface, t, ε_eff,
	                             radiation_model[, face]; a_eff = ε_eff)

Decompose the boundary heat flux into named components for diagnostics.
`a_eff` is the effective surface absorptivity for incident radiation
(Kirchhoff default `ε_eff`); it must match what the surface Newton used so
that `absorbed` reproduces the applied flux.

Returns a NamedTuple with:
- `absorbed`: incident radiation absorbed *at the boundary face* [W/m²]
- `reradiation`: re-radiation loss [W/m²] (positive = leaving)
- `convection`: convective loss [W/m²] (positive = leaving)
- `prescribed`: prescribed flux [W/m²]
- `net`: the flux the surface energy balance actually applies at the face
  [W/m²], positive into the domain

Convention: `net` is the **applied** face flux — per component it matches
`evaluate_bc_flux_and_derivative`, so the total satisfies the Newton balance
`net = q_cond` at convergence. Under `BEER_LAMBERT` the incident radiation
enters the domain *volumetrically*,
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
    radiation_model::RadiationModel;
    a_eff::Float64 = ε_eff,
)
    if radiation_model == BEER_LAMBERT
        # In-depth absorption: the incident radiation enters volumetrically,
        # nothing is absorbed at the face, and the Newton balance applies zero
        # flux for this BC.
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
            a_eff
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

# Overloads with radiation_model for BCs that don't care about it
function compute_boundary_heat_fluxes(
    bc::AdiabaticBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel;
    a_eff::Float64 = ε_eff,
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
    ::RadiationModel;
    a_eff::Float64 = ε_eff,
)
    q, _ = evaluate_bc_flux_and_derivative(bc, T_surface, t, ε_eff)
    return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = q, net = q)
end

function compute_boundary_heat_fluxes(
    bc::RadiativeBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel;
    a_eff::Float64 = ε_eff,
)
    q, _ = evaluate_bc_flux_and_derivative(bc, T_surface, t, ε_eff)
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
    ::RadiationModel;
    a_eff::Float64 = ε_eff,
)
    # Prescribed temperature: the applied flux is whatever conduction demands
    # (k_eff·(T_bc − T_cell)/Δz_half), not computable from (T_surface, t)
    # alone. The decomposition is all-zero here; the saved-solution extraction
    # reports the solved q_conducted as the net applied flux.
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
    ::RadiationModel;
    a_eff::Float64 = ε_eff,
)
    # Fallback for other BC types
    q_total, _ = evaluate_bc_flux_and_derivative(bc, T_surface, t, ε_eff)
    return (
        absorbed = 0.0,
        reradiation = 0.0,
        convection = 0.0,
        prescribed = q_total,
        net = q_total,
    )
end

# ----------------------------------------------------------------------------
# Blowing-aware decomposition (trailing `face` argument)
# ----------------------------------------------------------------------------
# `face` is the `BlowingFace` at the converged T_s (outward species fluxes,
# surface temperature, material). Mirrors `evaluate_bc_flux_and_derivative`:
# only convective children rescale their h by their blowing model's ratio;
# every other BC ignores `face`. This keeps the diagnostic decomposition
# consistent with the flux the surface Newton actually applied.

function compute_boundary_heat_fluxes(
    bc::AbstractThermalBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel,
    face::BlowingFace;
    a_eff::Float64 = ε_eff,
)
    return compute_boundary_heat_fluxes(
        bc,
        T_surface,
        t,
        ε_eff,
        radiation_model;
        a_eff = a_eff,
    )
end

function compute_boundary_heat_fluxes(
    bc::ConvectiveBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel,
    face::BlowingFace;
    a_eff::Float64 = ε_eff,
)
    h = evaluate_param(bc.h, t)
    T_inf = evaluate_param(bc.T_inf, t)
    χ = h > 0.0 ? blowing_ratio(bc.blowing, face, h, T_inf) : 1.0
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
    face::BlowingFace;
    a_eff::Float64 = ε_eff,
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
            face;
            a_eff = a_eff,
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
