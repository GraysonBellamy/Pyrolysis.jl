"""
Boundary flux implementation.

Reference: Design Document Section 2.5, 7.2

This module provides flux evaluation functions for all boundary condition types:
- Thermal BCs: heat flux into domain [W/m²]
- Mass BCs: mass flux into domain [kg/(m²·s)]
"""

using ..Core: STEFAN_BOLTZMANN
using ..Core: evaluate_param, evaluate_param_T  # Unified parameter evaluation from Core
using ..Core: RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT, P1_QUASI_STEADY

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
@inline function thermal_bc_flux(bc::RadiativeFluxBC, T_surface, t, ε_eff, radiation_model::RadiationModel)
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

# Dirichlet: no flux evaluation (handled differently in solver)
# Return NaN to indicate this BC type doesn't produce a flux
@inline thermal_bc_flux(::DirichletThermalBC, T_surface, t, ε_eff) = NaN
@inline thermal_bc_flux(::DirichletThermalBC, T_surface, t, ε_eff, ::RadiationModel) = NaN

# Combined: sum of all component fluxes using generated function for type stability
# NOTE: radiation_model is REQUIRED for CombinedThermalBC since it may contain RadiativeFluxBC
@generated function thermal_bc_flux(bc::CombinedThermalBC{BCs}, T_surface, t, ε_eff, radiation_model::RadiationModel) where {BCs <: Tuple}
	N = length(BCs.parameters)
	exprs = [:(thermal_bc_flux(bc.bcs[$i], T_surface, t, ε_eff, radiation_model)) for i in 1:N]
	return :(+($(exprs...)))
end

# Fallback for BCs that don't care about radiation_model (just ignore it)
@inline thermal_bc_flux(::AdiabaticBC, T_surface, t, ε_eff, ::RadiationModel) = 0.0
@inline thermal_bc_flux(bc::HeatFluxBC{<:Real}, T_surface, t, ε_eff, ::RadiationModel) = Float64(bc.flux)
@inline thermal_bc_flux(bc::HeatFluxBC{<:Function}, T_surface, t, ε_eff, ::RadiationModel) = bc.flux(t)
@inline thermal_bc_flux(bc::ConvectiveBC, T_surface, t, ε_eff, ::RadiationModel) = thermal_bc_flux(bc, T_surface, t, ε_eff)
@inline thermal_bc_flux(bc::RadiativeBC, T_surface, t, ε_eff, ::RadiationModel) = thermal_bc_flux(bc, T_surface, t, ε_eff)
@inline thermal_bc_flux(bc::FlameBC, T_surface, t, ε_eff::Float64, ::RadiationModel) = 0.0

"""
	thermal_bc_flux(bc::FlameBC, T_surface, t, mass_fluxes) -> Float64

Evaluate flame BC with critical mass flux ignition model.

# Arguments
- `bc`: FlameBC boundary condition
- `T_surface`: Surface temperature [K]
- `t`: Current time [s]
- `mass_fluxes`: Dict or vector of mass fluxes [kg/(m²·s)] for each gas component

# Returns
Flame heat flux [W/m²]. The value depends on the flame model:
- Binary mode (ramp_width=0): Full q_flame_max when ignited, 0 otherwise
- Gradual mode (ramp_width>0): Sigmoid ramp based on CI

The Combustion Index CI = Σⱼ(Jⱼ/Jⱼ,crit)
- Flame ignites when CI ≥ 1.0
- Flame extinguishes when CI < extinction_criterion (hysteresis)

# Note
FlameBC is stateless - flame state must be provided externally via the 
`flame_active` parameter. The solver maintains this state in 
`PyrolysisODEFunction.flame_state`.
"""
@inline function thermal_bc_flux(bc::FlameBC, T_surface, t, mass_fluxes::Dict, flame_active::Bool)
	# Compute combustion index
	CI = _compute_combustion_index(bc.critical_mass_fluxes, mass_fluxes)
	return _evaluate_flame_flux(bc, CI, t, flame_active)
end

# Vector version for mass fluxes indexed by component number
@inline function thermal_bc_flux(bc::FlameBC, T_surface, t, mass_fluxes::AbstractVector, flame_active::Bool)
	CI = _compute_combustion_index_vector(bc.critical_mass_fluxes, mass_fluxes)
	return _evaluate_flame_flux(bc, CI, t, flame_active)
end

# FlameBC with ε_eff but no mass fluxes - returns 0 (no ignition possible)
# This signature is used by CombinedThermalBC
@inline thermal_bc_flux(bc::FlameBC, T_surface, t, ε_eff::Float64) = 0.0

# Overloads without flame_active for simple cases
@inline function thermal_bc_flux(bc::FlameBC, T_surface, t, mass_fluxes::Dict)
	CI = _compute_combustion_index(bc.critical_mass_fluxes, mass_fluxes)
	return _evaluate_flame_flux_simple(bc, CI, t)
end

@inline function thermal_bc_flux(bc::FlameBC, T_surface, t, mass_fluxes::AbstractVector)
	CI = _compute_combustion_index_vector(bc.critical_mass_fluxes, mass_fluxes)
	return _evaluate_flame_flux_simple(bc, CI, t)
end

"""
Compute combustion index from Dict of mass fluxes.
"""
@inline function _compute_combustion_index(critical_mass_fluxes, mass_fluxes::Dict)
	CI = 0.0
	for (component_idx, J_crit) in critical_mass_fluxes
		J = get(mass_fluxes, component_idx, 0.0)
		if J_crit > 0
			CI += J / J_crit
		end
	end
	return CI
end

"""
Compute combustion index from vector of mass fluxes.
"""
@inline function _compute_combustion_index_vector(critical_mass_fluxes, mass_fluxes::AbstractVector)
	CI = 0.0
	for (component_idx, J_crit) in critical_mass_fluxes
		if component_idx <= length(mass_fluxes) && J_crit > 0
			CI += mass_fluxes[component_idx] / J_crit
		end
	end
	return CI
end

"""
Evaluate flame heat flux based on combustion index, with hysteresis and optional gradual onset.

# Hysteresis Logic
- If flame is off and CI ≥ 1.0: ignite (turn on)
- If flame is on and CI < extinction_criterion: extinguish (turn off)
- Otherwise: maintain current state

# Gradual Onset (ramp_width > 0)
When flame is active, the heat flux ramps smoothly:
	q = q_max × sigmoid((CI - 1.0) / ramp_width)
This prevents discontinuities at ignition/extinction.

# Arguments
- `bc`: FlameBC boundary condition
- `CI`: Combustion index
- `t`: Current time
- `flame_active`: Current flame state (from external state storage)
"""
@inline function _evaluate_flame_flux(bc::FlameBC, CI::Float64, t, flame_active::Bool)
	# If flame is off, no heat flux
	if !flame_active
		return 0.0
	end

	# Flame is on - compute heat flux
	q_max = evaluate_param(bc.q_flame_max, t)

	# Check for gradual onset mode
	if bc.ramp_width > 0.0
		# Sigmoid ramp: q = q_max × σ((CI - 1.0) / ramp_width)
		# σ(x) = 1 / (1 + exp(-x))
		x = (CI - 1.0) / bc.ramp_width
		# Clamp to prevent overflow in exp
		x_clamped = clamp(x, -20.0, 20.0)
		sigmoid = 1.0 / (1.0 + exp(-x_clamped))
		return q_max * sigmoid
	else
		# Binary mode: full heat flux when ignited
		return q_max
	end
end

"""
Simple flame flux evaluation without hysteresis (for backward compatibility).
Flame is active when CI >= 1.0, inactive otherwise.
"""
@inline function _evaluate_flame_flux_simple(bc::FlameBC, CI::Float64, t)
	if CI < 1.0
		return 0.0
	end

	q_max = evaluate_param(bc.q_flame_max, t)

	if bc.ramp_width > 0.0
		x = (CI - 1.0) / bc.ramp_width
		x_clamped = clamp(x, -20.0, 20.0)
		sigmoid = 1.0 / (1.0 + exp(-x_clamped))
		return q_max * sigmoid
	else
		return q_max
	end
end

"""
	update_flame_state!(flame_state, boundary_id, bc, CI)

Update flame state with hysteresis logic. Returns the new state.

# Arguments
- `flame_state`: Dict{Int, Bool} mapping boundary_id to flame active state
- `boundary_id`: The boundary ID
- `bc`: FlameBC boundary condition
- `CI`: Current combustion index

# Returns
The new flame state (true = active, false = inactive)
"""
function update_flame_state!(flame_state::Dict{Int, Bool}, boundary_id::Int,
	bc::FlameBC, CI::Float64)
	current = get(flame_state, boundary_id, false)

	if current
		# Flame is currently on - check for extinction
		if CI < bc.extinction_criterion
			flame_state[boundary_id] = false
		end
	else
		# Flame is currently off - check for ignition
		if CI >= 1.0
			flame_state[boundary_id] = true
		end
	end

	return get(flame_state, boundary_id, false)
end

# ============================================================================
# Mass Boundary Condition Flux Evaluation
# ============================================================================

"""
	mass_bc_flux(bc::AbstractMassBC, ξ_surface, t) -> Float64

Evaluate the mass flux from a mass boundary condition.

# Arguments
- `bc`: The mass boundary condition
- `ξ_surface`: Surface mass concentration [kg/m³]
- `t`: Current time [s]

# Returns
Mass flux [kg/(m²·s)], positive into the domain.
"""
function mass_bc_flux end

# Impermeable: zero flux
@inline mass_bc_flux(::ImpermeableBC, ξ_surface, t) = 0.0

# Prescribed mass flux
@inline mass_bc_flux(bc::MassFluxBC{<:Real}, ξ_surface, t) = Float64(bc.flux)
@inline mass_bc_flux(bc::MassFluxBC{<:Function}, ξ_surface, t) = bc.flux(t)

# Convective mass transfer: J = h_m(ξ∞ - ξ_surface)
@inline function mass_bc_flux(bc::ConvectiveMassBC, ξ_surface, t)
	h_m = evaluate_param(bc.h_m, t)
	ξ_inf = evaluate_param(bc.ξ_inf, t)
	return h_m * (ξ_inf - ξ_surface)
end

# Combined: sum of all component fluxes using generated function for type stability
@generated function mass_bc_flux(bc::CombinedMassBC{BCs}, ξ_surface, t) where {BCs <: Tuple}
	N = length(BCs.parameters)
	exprs = [:(mass_bc_flux(bc.bcs[$i], ξ_surface, t)) for i in 1:N]
	return :(+($(exprs...)))
end

# ============================================================================
# Mass BC Flux with Diffusion Parameters (λρ, ρ, d)
# ============================================================================
# These overloads are used when the discretization provides the combined transport 
# coefficient λρ (= λ × ρ), density ρ, and distance d from cell center to boundary.
# This is consistent with the rigorous λρ formulation for gas transport.

"""
	mass_bc_flux(bc::AbstractMassBC, ξ_surface, λρ, ρ, d, t) -> Float64

Evaluate mass flux with diffusion parameters available.

# Arguments
- `bc`: The mass boundary condition
- `ξ_surface`: Surface mass concentration [kg/m³]  
- `λρ`: Combined transport coefficient λ·ρ [kg/(m·s)]
- `ρ`: Mixture density at surface [kg/m³]
- `d`: Distance from cell center to boundary [m]
- `t`: Current time [s]

# Returns
Mass flux [kg/(m²·s)], positive into the domain (negative = gas leaving).
"""

# Default fallback: ignore λρ, ρ, d for non-diffusive BCs
@inline mass_bc_flux(bc::ImpermeableBC, ξ_surface, λρ, ρ, d, t) = 0.0
@inline mass_bc_flux(bc::MassFluxBC{<:Real}, ξ_surface, λρ, ρ, d, t) = Float64(bc.flux)
@inline mass_bc_flux(bc::MassFluxBC{<:Function}, ξ_surface, λρ, ρ, d, t) = bc.flux(t)
@inline function mass_bc_flux(bc::ConvectiveMassBC, ξ_surface, λρ, ρ, d, t)
	h_m = evaluate_param(bc.h_m, t)
	ξ_inf = evaluate_param(bc.ξ_inf, t)
	return h_m * (ξ_inf - ξ_surface)
end

"""
	mass_bc_flux(bc::DiffusiveMassBC, ξ_surface, λρ, ρ, d, t) -> Float64

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
@inline function mass_bc_flux(bc::DiffusiveMassBC, ξ_surface, λρ, ρ, d, t)
	h_m = evaluate_param(bc.h_m, t)

	# Extract λ for combined resistance calculation
	λ = λρ / max(ρ, 1e-20)

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
@generated function mass_bc_flux(bc::CombinedMassBC{BCs}, ξ_surface, λρ, ρ, d, t) where {BCs <: Tuple}
	N = length(BCs.parameters)
	exprs = [:(mass_bc_flux(bc.bcs[$i], ξ_surface, λρ, ρ, d, t)) for i in 1:N]
	return :(+($(exprs...)))
end

# ============================================================================
# Pressure BC Mass Flux (for CombinedMassBC compatibility)
# ============================================================================
# Pressure BCs (DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC)
# control Darcy velocity via apply_pressure_bc_to_darcy_velocity!, not diffusive flux.
# These methods return 0 so they can be combined with DiffusiveMassBC in CombinedMassBC.

using ..Core: DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC

# Simple signature (no diffusion parameters)
@inline mass_bc_flux(::DirichletPressureBC, ξ_surface, t) = 0.0
@inline mass_bc_flux(::NeumannPressureBC, ξ_surface, t) = 0.0
@inline mass_bc_flux(::ConvectivePressureBC, ξ_surface, t) = 0.0

# Full signature with diffusion parameters (λρ, ρ, d)
@inline mass_bc_flux(::DirichletPressureBC, ξ_surface, λρ, ρ, d, t) = 0.0
@inline mass_bc_flux(::NeumannPressureBC, ξ_surface, λρ, ρ, d, t) = 0.0
@inline mass_bc_flux(::ConvectivePressureBC, ξ_surface, λρ, ρ, d, t) = 0.0

# ============================================================================
# Mesh Boundary Condition Evaluation
# ============================================================================

"""
	mesh_bc_velocity(bc::AbstractMeshBC, t) -> Float64

Evaluate the mesh velocity at a boundary.

# Returns
Mesh velocity [m/s] at boundary.
"""
function mesh_bc_velocity end

@inline mesh_bc_velocity(::FixedMeshBC, t) = 0.0
@inline mesh_bc_velocity(bc::PrescribedVelocityBC{<:Real}, t) = Float64(bc.velocity)
@inline mesh_bc_velocity(bc::PrescribedVelocityBC{<:Function}, t) = bc.velocity(t)

# FreeSurfaceBC and SlidingBC need additional information (material velocity, etc.)
# These are handled in the ALE module with full state information

# ============================================================================
# Boundary Flux Decomposition for Diagnostics
# ============================================================================

"""
	compute_boundary_heat_fluxes(bc::AbstractThermalBC, T_surface, t, ε_eff)

Compute individual heat flux components for diagnostics.

Returns a NamedTuple with:
- `absorbed`: Absorbed radiation [W/m²] 
- `reradiation`: Re-radiation losses [W/m²] (positive = leaving)
- `convection`: Convective losses [W/m²] (positive = leaving)  
- `net`: Net heat entering domain [W/m²]

For combined BCs, iterates through each component and categorizes appropriately.
"""
function compute_boundary_heat_fluxes(bc::AbstractThermalBC, T_surface::Float64, t::Float64, ε_eff::Float64)
	# Default: compute total flux, can't decompose
	q_total = thermal_bc_flux(bc, T_surface, t, ε_eff)
	return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = q_total, net = q_total)
end

function compute_boundary_heat_fluxes(::AdiabaticBC, T_surface::Float64, t::Float64, ε_eff::Float64)
	return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = 0.0, net = 0.0)
end

function compute_boundary_heat_fluxes(bc::HeatFluxBC, T_surface::Float64, t::Float64, ε_eff::Float64)
	q = thermal_bc_flux(bc, T_surface, t, ε_eff)
	return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = q, net = q)
end

function compute_boundary_heat_fluxes(bc::ConvectiveBC, T_surface::Float64, t::Float64, ε_eff::Float64)
	# Convective: q = h(T∞ - T_surface)
	# Positive when T∞ > T_surface (heat entering)
	# We want to report loss (positive = leaving), so negate
	q = thermal_bc_flux(bc, T_surface, t, ε_eff)
	q_loss = -q  # Convert to "loss" convention
	return (absorbed = 0.0, reradiation = 0.0, convection = q_loss, prescribed = 0.0, net = q)
end

function compute_boundary_heat_fluxes(bc::RadiativeBC, T_surface::Float64, t::Float64, ε_eff::Float64)
	# Radiative: q = F × ε × σ × (T∞⁴ - T_surface⁴)
	# Positive when T∞ > T_surface (heat entering)
	# We report re-radiation as loss (positive = leaving), so negate
	q = thermal_bc_flux(bc, T_surface, t, ε_eff)
	q_loss = -q  # Convert to "loss" convention (positive = emitting more than receiving)
	return (absorbed = 0.0, reradiation = q_loss, convection = 0.0, prescribed = 0.0, net = q)
end

function compute_boundary_heat_fluxes(bc::RadiativeFluxBC, T_surface::Float64, t::Float64, ε_eff::Float64, radiation_model::RadiationModel)
	# RadiativeFluxBC: absorbed radiation (positive = entering)
	# NOTE: radiation_model is REQUIRED - determines surface vs in-depth absorption
	flux = bc.flux isa Function ? bc.flux(t) : Float64(bc.flux)
	α = if bc.absorptivity === nothing
		ε_eff
	else
		bc.absorptivity isa Function ? bc.absorptivity(t) : Float64(bc.absorptivity)
	end
	q_absorbed = flux * α

	if radiation_model == BEER_LAMBERT || radiation_model == P1_QUASI_STEADY || radiation_model == SURFACE_ABSORPTION
		# Volumetric absorption: radiation enters domain volumetrically (Beer-Lambert) 
		# or in surface cell (SURFACE_ABSORPTION). The absorbed flux is handled by 
		# compute_radiation_source!, so net=0 here to avoid double-counting in boundary flux.
		# We still report q_absorbed for diagnostics.
		return (absorbed = q_absorbed, reradiation = 0.0, convection = 0.0, prescribed = 0.0, net = 0.0)
	else
		# NO_RADIATION: radiation handled via boundary flux only
		return (absorbed = q_absorbed, reradiation = 0.0, convection = 0.0, prescribed = 0.0, net = q_absorbed)
	end
end

function compute_boundary_heat_fluxes(bc::CombinedThermalBC, T_surface::Float64, t::Float64, ε_eff::Float64, radiation_model::RadiationModel)
	# Sum contributions from all component BCs
	# NOTE: radiation_model is REQUIRED since CombinedThermalBC may contain RadiativeFluxBC
	absorbed = 0.0
	reradiation = 0.0
	convection = 0.0
	prescribed = 0.0

	for inner_bc in bc.bcs
		fluxes = compute_boundary_heat_fluxes(inner_bc, T_surface, t, ε_eff, radiation_model)
		absorbed += fluxes.absorbed
		reradiation += fluxes.reradiation
		convection += fluxes.convection
		prescribed += fluxes.prescribed
	end

	net = absorbed - reradiation - convection + prescribed
	return (absorbed = absorbed, reradiation = reradiation, convection = convection, prescribed = prescribed, net = net)
end

# Overloads with radiation_model for BCs that don't care about it
function compute_boundary_heat_fluxes(bc::AdiabaticBC, T_surface::Float64, t::Float64, ε_eff::Float64, ::RadiationModel)
	return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = 0.0, net = 0.0)
end

function compute_boundary_heat_fluxes(bc::HeatFluxBC, T_surface::Float64, t::Float64, ε_eff::Float64, ::RadiationModel)
	q = thermal_bc_flux(bc, T_surface, t, ε_eff)
	return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = q, net = q)
end

function compute_boundary_heat_fluxes(bc::ConvectiveBC, T_surface::Float64, t::Float64, ε_eff::Float64, ::RadiationModel)
	q = thermal_bc_flux(bc, T_surface, t, ε_eff)
	q_loss = -q
	return (absorbed = 0.0, reradiation = 0.0, convection = q_loss, prescribed = 0.0, net = q)
end

function compute_boundary_heat_fluxes(bc::RadiativeBC, T_surface::Float64, t::Float64, ε_eff::Float64, ::RadiationModel)
	q = thermal_bc_flux(bc, T_surface, t, ε_eff)
	q_loss = -q
	return (absorbed = 0.0, reradiation = q_loss, convection = 0.0, prescribed = 0.0, net = q)
end

function compute_boundary_heat_fluxes(bc::AbstractThermalBC, T_surface::Float64, t::Float64, ε_eff::Float64, ::RadiationModel)
	# Fallback for other BC types
	q_total = thermal_bc_flux(bc, T_surface, t, ε_eff)
	return (absorbed = 0.0, reradiation = 0.0, convection = 0.0, prescribed = q_total, net = q_total)
end

"""
	compute_surface_mass_flux(bc::AbstractMassBC, ξ_surface, λρ, ρ, d, t) -> Float64

Compute the mass flux leaving the surface (positive = leaving).

This differs from mass_bc_flux which uses flux-into-domain convention.
"""
function compute_surface_mass_flux(bc::AbstractMassBC, ξ_surface::Float64, λρ::Float64, ρ::Float64, d::Float64, t::Float64)
	# mass_bc_flux returns negative for gas escaping (flux into domain convention)
	# We want positive = gas leaving
	return -mass_bc_flux(bc, ξ_surface, λρ, ρ, d, t)
end

# Simple overload without diffusion parameters
function compute_surface_mass_flux(bc::AbstractMassBC, ξ_surface::Float64, t::Float64)
	return -mass_bc_flux(bc, ξ_surface, t)
end
