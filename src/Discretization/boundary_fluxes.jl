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
# Mass Boundary Condition Flux Evaluation
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
function mass_bc_flux end

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

`apply_mass_bc_flux!` (finite_volume.jl) orients it to the face by negating
at the top boundary: J[f] = -J_bc. This gives a positive flux in the +z
direction when gas escapes at the top surface, as the divergence expects.

Summary of sign conventions:
- mass_bc_flux():               Negative = gas escaping (flux-into-domain convention)
- compute_surface_mass_flux():  Positive = gas escaping
- J[face] in finite volume:     Positive = flux in +z direction
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

@inline mass_bc_flux(::DirichletPressureBC, ξ_surface, λρ, ρ, d, t, j::Integer) = 0.0
@inline mass_bc_flux(::NeumannPressureBC, ξ_surface, λρ, ρ, d, t, j::Integer) = 0.0
@inline mass_bc_flux(::ConvectivePressureBC, ξ_surface, λρ, ρ, d, t, j::Integer) = 0.0

# ============================================================================
# Boundary-face gas mass flux (film law + Darcy venting)
# ============================================================================

"""
	_bc_boundary_face_fluxes(T_cell, ξ_cell, mass_bc, material,
	                         d_to_face, t, ::Val{tag}, darcy_active)
	    -> NTuple{NC}

Per-component boundary-face gas mass flux per area, **positive into the
domain**: the diffusive film law plus (when `darcy_active`) the pressure-BC
advective venting that mirrors the boundary-face branch of
`compute_gas_fluxes_darcy_fick!`. Non-gas components are zero. Generic over
Dual scalar types. Shared by the Jacobian's boundary operators
(`BCMassOperator`, `BCConvectionOperator`, `BCThermalOperator`) and the
solution extraction, so all of them see the flux the residual applies.
"""
function _bc_boundary_face_fluxes(
    T_cell::Real,
    ξ_cell::NTuple{NC,Tu},
    mass_bc,
    material::Material{NC},
    d_to_face::Real,
    t::Real,
    ::Val{tag},
    darcy_active::Bool,
) where {NC,Tu,tag}

    λ_eff = effective_gas_transfer(ξ_cell, T_cell, material)
    ρ_eff = mixture_density(ξ_cell, T_cell, material)
    λρ = λ_eff * ρ_eff

    Tv = promote_type(typeof(T_cell), Tu)
    F = ntuple(Val(NC)) do j
        if is_gas(material.components[j])
            Tv(mass_bc_flux(mass_bc, ξ_cell[j], λρ, ρ_eff, d_to_face, t, j))
        else
            zero(Tv)
        end
    end

    if darcy_active
        # Advective venting: outward Darcy speed advecting the intrinsic
        # density ξ/φ out of the boundary cell (inflow advects only the
        # zero-concentration ambient, hence the one-sided clamp).
        φ = porosity(ξ_cell, T_cell, material)
        P_int = ideal_gas_pressure(ξ_cell, T_cell, φ, material)
        κ_int = effective_permeability(ξ_cell, T_cell, material)
        μ_int = mixture_gas_viscosity(ξ_cell, T_cell, material)
        v_z = boundary_darcy_velocity(mass_bc, P_int, κ_int, μ_int, d_to_face, Val(tag))
        v_out = tag === :top ? v_z : -v_z
        if v_out > 0
            coef = v_out / max(φ, 1e-12)
            # Return a fresh tuple instead of rebinding `F` — a self-referential
            # rebind would box the variable (type-unstable, and inference can't
            # prove the boxed read defined).
            return ntuple(Val(NC)) do j
                if is_gas(material.components[j])
                    F[j] - coef * ξ_cell[j]
                else
                    F[j]
                end
            end
        end
    end

    return F
end

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
`net = q_cond` at convergence. Under `BEER_LAMBERT` (and the
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
    radiation_model::RadiationModel;
    a_eff::Float64 = ε_eff,
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

function compute_boundary_heat_fluxes(
    bc::CombinedThermalBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    radiation_model::RadiationModel;
    a_eff::Float64 = ε_eff,
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
        fluxes = compute_boundary_heat_fluxes(
            inner_bc,
            T_surface,
            t,
            ε_eff,
            radiation_model;
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
    bc::ConvectiveBC,
    T_surface::Float64,
    t::Float64,
    ε_eff::Float64,
    ::RadiationModel;
    a_eff::Float64 = ε_eff,
)
    q, _ = evaluate_bc_flux_and_derivative(bc, T_surface, t, ε_eff)
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
    radiation_model::RadiationModel;
    a_eff::Float64 = ε_eff,
)
    # Under P1 the boundary balance applies zero for this BC (emission and
    # ambient exchange are volumetric) — `net` must match the applied flux.
    q, _ = evaluate_bc_flux_and_derivative(bc, T_surface, t, ε_eff, radiation_model)
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
