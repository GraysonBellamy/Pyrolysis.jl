"""
Material component definitions.

Components represent distinct phases (solid, liquid, gas) within a material,
each with their own thermophysical properties.

Reference: Design Document Section 6.2
"""

# ============================================================================
# Phase Enumeration
# ============================================================================

"""
	Phase

Enumeration of material phases.

- `SOLID`: Solid phase (does not transport)
- `LIQUID`: Liquid phase (may transport)
- `GAS`: Gas phase (transports through porous matrix)
"""
@enum Phase begin
	SOLID = 1
	LIQUID = 2
	GAS = 3
end

# ============================================================================
# Component Type
# ============================================================================

"""
	Component{Pρ,Pc,Pk,Pλ,Pκ} <: AbstractComponent

A material component with temperature-dependent thermophysical properties.

# Type Parameters
- `Pρ`: Type of density property function
- `Pc`: Type of heat capacity property function
- `Pk`: Type of thermal conductivity property function
- `Pλ`: Type of gas transfer coefficient property function
- `Pκ`: Type of permeability property function (Nothing for gases)

# Fields
- `name::Symbol`: Component identifier
- `phase::Phase`: Physical phase (SOLID, LIQUID, GAS)
- `molar_mass::Float64`: Molar mass [kg/mol] (for gases)
- `density::Pρ`: Density [kg/m³] or density function
- `heat_capacity::Pc`: Specific heat capacity [J/(kg·K)]
- `thermal_conductivity::Pk`: Thermal conductivity [W/(m·K)]
- `gas_transfer_coeff::Pλ`: Gas transfer coefficient [m²/s]
- `permeability::Pκ`: Permeability [m²] for Darcy flow (Nothing for gases)
- `emissivity::Float64`: Surface emissivity [-]
- `absorption_coeff::Float64`: Mass absorption coefficient [m²/kg]
- `swelling_factor::Float64`: Volume expansion factor [0-1]

# Example
```julia
virgin = Component(
	name = :virgin,
	phase = SOLID,
	molar_mass = 0.0,
	density = ConstantProperty(1190.0),
	heat_capacity = LinearProperty(1420.0, 2.5),
	thermal_conductivity = ConstantProperty(0.21),
	gas_transfer_coeff = ConstantProperty(0.0),
	permeability = ConstantProperty(1e-15),
	emissivity = 0.86,
	absorption_coeff = 1.5,
	swelling_factor = 0.0
)
```
"""
struct Component{
	Pρ <: AbstractPropertyFunction,
	Pc <: AbstractPropertyFunction,
	Pk <: AbstractPropertyFunction,
	Pλ <: AbstractPropertyFunction,
	Pκ <: Union{Nothing, AbstractPropertyFunction},
} <: AbstractComponent
	name::Symbol
	phase::Phase
	molar_mass::Float64           # kg/mol (for gases)

	# Thermophysical properties (all T-dependent)
	density::Pρ                   # kg/m³
	heat_capacity::Pc             # J/(kg·K)
	thermal_conductivity::Pk      # W/(m·K)
	gas_transfer_coeff::Pλ        # m²/s
	permeability::Pκ              # m² (Nothing for gases)

	# Radiative properties
	emissivity::Float64           # [-]
	absorption_coeff::Float64     # m²/kg

	# Volume change
	swelling_factor::Float64      # 0-1
end

# ============================================================================
# Convenience Constructors
# ============================================================================

"""
	SolidComponent(name; ρ, c, k, ε=0.9, α=0.0, γ=1.0, λ=0.0, κ=0.0)

Create a solid component with simplified interface.

# Arguments
- `name::Symbol`: Component identifier
- `ρ`: Density [kg/m³] (number or property function)
- `c`: Heat capacity [J/(kg·K)] (number or property function)
- `k`: Thermal conductivity [W/(m·K)] (number or property function)
- `ε=0.9`: Surface emissivity [-]
- `α=0.0`: Mass absorption coefficient [m²/kg]
- `γ=1.0`: Swelling factor [-] (default 1.0: solid contributes to volume change)
- `λ=0.0`: Gas transfer coefficient [m²/s] for mixture property calculation
- `κ=0.0`: Permeability [m²] for Darcy flow (0 = impermeable)

# Note on swelling_factor
The swelling factor controls whether this component contributes to volume change:
- γ=1.0: Full contribution (default for solids - when solid depletes, volume shrinks)
- γ=0.0: No contribution (component doesn't affect matrix volume)
- γ∈(0,1): Partial contribution (advanced use cases)

# Note on gas_transfer_coeff (λ)
Although solids don't transport themselves, they contribute to the effective gas
transfer coefficient of the mixture through volume-weighted mixing rules. Set to
a non-zero value (e.g., 2e-5 m²/s) to match ThermaKin conventions.

# Note on permeability (κ)
Permeability determines resistance to pressure-driven (Darcy) gas flow through
the solid matrix. Typical values:
- Virgin wood: ~10⁻¹⁵ m² (tight cellular structure)
- Char: ~10⁻¹² m² (porous with open channels)
- Impermeable: 0.0 m² (default)

# Example
```julia
virgin = SolidComponent(:virgin, ρ=1190, c=1420, k=0.21, ε=0.86)
char = SolidComponent(:char, ρ=200, c=1000, k=0.10, ε=0.95, α=5.0, κ=1e-12)
# Char that doesn't contribute to volume change:
rigid_char = SolidComponent(:char, ρ=200, c=1000, k=0.10, γ=0.0)
```
"""
function SolidComponent(name::Symbol;
	ρ, c, k,
	ε::Real = 0.9,
	α::Real = 0.0,
	γ::Real = 1.0,
	λ::Real = 0.0,
	κ::Real = 0.0)
	Component(
		name,
		SOLID,
		0.0,  # molar mass not used for solids
		to_property(ρ),
		to_property(c),
		to_property(k),
		to_property(λ),  # contributes to effective gas transfer coefficient
		to_property(κ),  # permeability for Darcy flow
		Float64(ε),
		Float64(α),
		Float64(γ),
	)
end

"""
	GaseousComponent(name; M, c, k, λ, ε=0.0, α=0.0, γ=0.0)

Create a gaseous component with simplified interface.

For ideal gases, density is computed from the ideal gas law.
Gases do NOT have permeability - they flow THROUGH the solid matrix.

# Arguments
- `name::Symbol`: Component identifier
- `M`: Molar mass [kg/mol]
- `c`: Heat capacity [J/(kg·K)] (number or property function)
- `k`: Thermal conductivity [W/(m·K)] (number or property function)
- `λ`: Gas transfer coefficient [m²/s] (number or property function)
- `ε=0.0`: Surface emissivity [-] (typically 0 for gases)
- `α=0.0`: Mass absorption coefficient [m²/kg]
- `γ=0.0`: Swelling factor [-] (default 0.0: gas escapes, doesn't affect matrix)

# Note on swelling_factor
The swelling factor controls whether this component contributes to volume change:
- γ=0.0: No contribution (default for gases - gas escapes through pores)
- γ=1.0: Full contribution (for intumescent materials where gas creates new volume)

# Note on permeability
Gases do not have permeability because they flow through the solid matrix.
Only solid and liquid components contribute to the effective permeability.

# Example
```julia
MMA = GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
water_vapor = GaseousComponent(:H2O, M=0.018, c=2000, k=0.025, λ=2e-5)
# Intumescent gas that creates new volume:
intumescent_gas = GaseousComponent(:gas, M=0.040, c=1000, k=0.02, λ=1e-5, γ=1.0)
```
"""
function GaseousComponent(name::Symbol;
	M::Real,
	c, k, λ,
	ε::Real = 0.0,
	α::Real = 0.0,
	γ::Real = 0.0)
	# Ideal gas density: ρ = PM/(RT)
	M_float = Float64(M)
	ρ_func = CallableProperty(T -> P_REF * M_float / (R_GAS * T))

	Component(
		name,
		GAS,
		M_float,
		ρ_func,
		to_property(c),
		to_property(k),
		to_property(λ),
		nothing,  # gases have no permeability - they flow THROUGH the matrix
		Float64(ε),
		Float64(α),
		Float64(γ),
	)
end

"""
	LiquidComponent(name; ρ, c, k, λ=2e-5, κ=0.0, ε=0.0, α=0.0, γ=0.0)

Create a liquid component with simplified interface.

# Arguments
- `name::Symbol`: Component identifier
- `ρ`: Density [kg/m³] (number or property function)
- `c`: Heat capacity [J/(kg·K)] (number or property function)
- `k`: Thermal conductivity [W/(m·K)] (number or property function)
- `λ=2e-5`: Transport coefficient [m²/s] for mixture property calculation
- `κ=0.0`: Permeability [m²] for Darcy flow (0 = blocks flow, typical for liquids)
- `ε=0.0`: Surface emissivity [-]
- `α=0.0`: Mass absorption coefficient [m²/kg]
- `γ=0.0`: Swelling factor [-]

# Note on gas_transfer_coeff (λ)
Like solids, liquids contribute to the effective gas transfer coefficient of the
mixture through volume-weighted mixing rules. The default value (2e-5 m²/s)
matches ThermaKin conventions.

# Note on permeability (κ)
Liquids typically block gas flow (κ=0), as the liquid fills pore space.
Set to non-zero only for special cases like partially saturated porous media.
"""
function LiquidComponent(name::Symbol;
	ρ, c, k,
	λ = 2e-5,
	κ::Real = 0.0,
	ε::Real = 0.0,
	α::Real = 0.0,
	γ::Real = 0.0)
	Component(
		name,
		LIQUID,
		0.0,  # molar mass not typically used for liquids
		to_property(ρ),
		to_property(c),
		to_property(k),
		to_property(λ),
		to_property(κ),  # permeability for Darcy flow
		Float64(ε),
		Float64(α),
		Float64(γ),
	)
end

# ============================================================================
# Property Accessors
# ============================================================================

"""
	get_density(comp::Component, T::Real)

Get component density at temperature T.
"""
@inline get_density(comp::Component, T::Real) = comp.density(T)

"""
	get_heat_capacity(comp::Component, T::Real)

Get component heat capacity at temperature T.
"""
@inline get_heat_capacity(comp::Component, T::Real) = comp.heat_capacity(T)

"""
	get_conductivity(comp::Component, T::Real)

Get component thermal conductivity at temperature T.
"""
@inline get_conductivity(comp::Component, T::Real) = comp.thermal_conductivity(T)

"""
	get_gas_transfer(comp::Component, T::Real)

Get component gas transfer coefficient at temperature T.
"""
@inline get_gas_transfer(comp::Component, T::Real) = comp.gas_transfer_coeff(T)

"""
	get_permeability(comp::Component, T::Real) -> Float64

Get component permeability at temperature T.
Returns 0.0 for gases (permeability = nothing).
"""
@inline function get_permeability(comp::Component, T::Real)
	if comp.permeability === nothing
		return 0.0
	else
		return comp.permeability(T)
	end
end

"""
	has_permeability(comp::Component) -> Bool

Check if component contributes to matrix permeability.
Returns false for gases (which flow through the matrix, not part of it).
"""
@inline has_permeability(comp::Component) = comp.permeability !== nothing

"""
	is_solid(comp::Component)

Check if component is in solid phase.
"""
@inline is_solid(comp::Component) = comp.phase == SOLID

"""
	is_liquid(comp::Component)

Check if component is in liquid phase.
"""
@inline is_liquid(comp::Component) = comp.phase == LIQUID

"""
	is_gas(comp::Component)

Check if component is in gas phase.
"""
@inline is_gas(comp::Component) = comp.phase == GAS

# ============================================================================
# Radiation Property Queries
# ============================================================================

"""
	has_absorption_coefficients(components::Tuple) -> Bool

Returns true if any solid component has a non-zero absorption coefficient.

This is used to warn users when SURFACE_ABSORPTION radiation model is selected
but material components have α values set (which will be ignored).

# Arguments
- `components`: Tuple of Component objects from a Material

# Returns
- `true` if any component has absorption_coeff > 0
- `false` otherwise

# Example
```julia
if has_absorption_coefficients(material.components)
	@warn "Absorption coefficients will be ignored with SURFACE_ABSORPTION model"
end
```
"""
function has_absorption_coefficients(components::Tuple)
	for comp in components
		if comp.absorption_coeff > 0.0
			return true
		end
	end
	return false
end
