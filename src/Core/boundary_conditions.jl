"""
Boundary condition type definitions.

Reference: Design Document Section 6.4
"""

# ============================================================================
# Validation Helpers
# ============================================================================

"""
	validate_bc_params(; kwargs...)

Validate boundary condition parameters. Called by BC constructors when
`validate=true` is passed.

# Supported validations
- `h`: Heat transfer coefficient must be non-negative
- `T`: Temperature must be positive (Kelvin)
- `ε`: Emissivity must be in [0, 1]
- `F`: View factor must be in [0, 1]
- `α`: Absorptivity must be in [0, 1]
- `h_m`: Mass transfer coefficient must be non-negative
- `P`: Pressure must be positive

Throws `ArgumentError` if validation fails.
"""
function validate_bc_params(;
	h::Union{Real, Nothing} = nothing,
	T::Union{Real, Nothing} = nothing,
	ε::Union{Real, Nothing} = nothing,
	F::Union{Real, Nothing} = nothing,
	α::Union{Real, Nothing} = nothing,
	h_m::Union{Real, Nothing} = nothing,
	P::Union{Real, Nothing} = nothing,
	flux::Union{Real, Nothing} = nothing)

	if h !== nothing && h < 0
		throw(ArgumentError("Heat transfer coefficient h must be non-negative, got $h"))
	end
	if T !== nothing && T <= 0
		throw(ArgumentError("Temperature must be positive (in Kelvin), got $T"))
	end
	if ε !== nothing && (ε < 0 || ε > 1)
		throw(ArgumentError("Emissivity ε must be in [0, 1], got $ε"))
	end
	if F !== nothing && (F < 0 || F > 1)
		throw(ArgumentError("View factor F must be in [0, 1], got $F"))
	end
	if α !== nothing && (α < 0 || α > 1)
		throw(ArgumentError("Absorptivity α must be in [0, 1], got $α"))
	end
	if h_m !== nothing && h_m < 0
		throw(ArgumentError("Mass transfer coefficient h_m must be non-negative, got $h_m"))
	end
	if P !== nothing && P <= 0
		throw(ArgumentError("Pressure must be positive, got $P"))
	end
	return nothing
end

"""
	_get_scalar_value(param)

Extract a scalar value from a parameter for validation purposes.
Returns nothing if the parameter is a function (can't validate at construction time)
or if the parameter is nothing (meaning use material property at runtime).
"""
_get_scalar_value(::Nothing) = nothing
_get_scalar_value(param::Real) = Float64(param)
_get_scalar_value(param::Function) = nothing
_get_scalar_value(param::AbstractPropertyFunction) = nothing

# ============================================================================
# Thermal Boundary Conditions
# ============================================================================

"""
	AdiabaticBC <: AbstractThermalBC

Zero heat flux (perfectly insulated) boundary condition.
"""
struct AdiabaticBC <: AbstractThermalBC end

"""
	HeatFluxBC{F} <: AbstractThermalBC

Prescribed heat flux boundary condition: q" = flux(t) or constant.

# Fields
- `flux::F`: Heat flux [W/m²], positive into domain

# Example
```julia
bc = HeatFluxBC(50e3)                    # Constant 50 kW/m²
bc = HeatFluxBC(t -> 50e3 * (1 - exp(-t/10)))  # Time-varying
```
"""
struct HeatFluxBC{F} <: AbstractThermalBC
	flux::F  # q"(t) or constant, positive = into domain
end

# Convenience constructor to convert integers to Float64
HeatFluxBC(flux::Integer) = HeatFluxBC{Float64}(Float64(flux))

"""
	ConvectiveBC{Fh,FT} <: AbstractThermalBC

Convective heat transfer boundary condition: q" = h(T∞ - T_surface)

# Fields
- `h::Fh`: Heat transfer coefficient [W/(m²·K)]
- `T_inf::FT`: Ambient temperature [K]

# Example
```julia
bc = ConvectiveBC(h=10.0, T_inf=300.0)
bc = ConvectiveBC(h=10.0, T_inf=t -> 300 + 0.1*t)  # Warming ambient
```
"""
struct ConvectiveBC{Fh, FT} <: AbstractThermalBC
	h::Fh           # h(t) or constant
	T_inf::FT       # T∞(t) or constant
end

function ConvectiveBC(; h, T_inf, validate::Bool = true)
	if validate
		h_val = _get_scalar_value(h)
		T_val = _get_scalar_value(T_inf)
		validate_bc_params(h = h_val, T = T_val)
	end
	ConvectiveBC(h, T_inf)
end

"""
	RadiativeBC{Fε,FT,FF} <: AbstractThermalBC

Radiative heat transfer boundary condition: q" = F × ε × σ × (T∞⁴ - T_surface⁴)

Models thermal radiation exchange between a surface and its surroundings.
The heat flux is positive when T∞ > T_surface (heating) and negative when
T_surface > T∞ (cooling/re-radiation loss).

# Fields
- `ε::Fε`: Surface emissivity [-], range 0-1. Can be:
  - `nothing` (default): Use effective emissivity from material properties (ε_eff = Σⱼ vⱼ εⱼ).
	This is the recommended setting as it automatically tracks composition changes.
  - A constant value (e.g., 0.9)
  - A temperature-dependent function (e.g., T -> 0.85 + 0.1*(T-300)/700)
  Typical values:
  - Polymers: 0.85 - 0.95
  - Char/carbon: 0.90 - 0.98
  - Polished metals: 0.02 - 0.10
  - Oxidized metals: 0.40 - 0.80
- `T_inf::FT`: Ambient/surroundings temperature [K]. Can be constant or time-dependent.
- `F::FF`: View factor [-], range 0-1. Fraction of radiation leaving the surface
  that reaches the surroundings at T_inf. Default is 1.0 (surface fully sees surroundings).
  Use values < 1 for:
  - Partially enclosed surfaces
  - Surfaces facing walls at different temperatures
  - Complex geometries with self-viewing

# Example
```julia
# Use effective emissivity from material (recommended, default)
bc = RadiativeBC(T_inf=300.0)

# Override with fixed emissivity
bc = RadiativeBC(ε=0.9, T_inf=300.0)

# With view factor for partially enclosed surface
bc = RadiativeBC(T_inf=300.0, F=0.7)

# Temperature-dependent emissivity (e.g., char formation)
bc = RadiativeBC(ε=T -> 0.85 + 0.1*(T-300)/700, T_inf=300.0)
```

# Note
When ε=nothing (default), the effective emissivity is computed as the volume-weighted
average of component emissivities: ε_eff = Σⱼ vⱼ εⱼ. This automatically tracks the
change in surface properties as material decomposes (e.g., polymer → char).
"""
struct RadiativeBC{Fε, FT, FF} <: AbstractThermalBC
	ε::Fε           # nothing = use effective, or ε(T), or constant
	T_inf::FT       # T∞(t) or constant
	F::FF           # View factor (0-1), default 1.0
end

# Constructors with optional emissivity and view factor
function RadiativeBC(; ε = nothing, T_inf, F = 1.0, validate::Bool = true)
	if validate
		ε_val = _get_scalar_value(ε)
		T_val = _get_scalar_value(T_inf)
		F_val = _get_scalar_value(F)
		validate_bc_params(ε = ε_val, T = T_val, F = F_val)
	end
	RadiativeBC(ε, T_inf, F)
end

"""
	DirichletThermalBC{FT} <: AbstractThermalBC

Prescribed temperature boundary condition: T = T_bc

# Fields
- `T::FT`: Boundary temperature [K]
"""
struct DirichletThermalBC{FT} <: AbstractThermalBC
	T::FT
end

"""
	CombinedThermalBC{BCs<:Tuple} <: AbstractThermalBC

Combined thermal boundary condition (sum of multiple BCs).

# Fields
- `bcs::BCs`: Tuple of boundary conditions

# Example
```julia
bc = HeatFluxBC(50e3) + RadiativeBC(ε=0.9, T_inf=300.0) + ConvectiveBC(h=10.0, T_inf=300.0)
```
"""
struct CombinedThermalBC{BCs <: Tuple} <: AbstractThermalBC
	bcs::BCs
end

# Combine BCs with + operator
Base.:+(bc1::AbstractThermalBC, bc2::AbstractThermalBC) =
	CombinedThermalBC((bc1, bc2))
Base.:+(cbc::CombinedThermalBC, bc::AbstractThermalBC) =
	CombinedThermalBC((cbc.bcs..., bc))
Base.:+(bc::AbstractThermalBC, cbc::CombinedThermalBC) =
	CombinedThermalBC((bc, cbc.bcs...))

"""
	RadiativeFluxBC{F,A} <: AbstractThermalBC

External radiative heat source (e.g., cone heater, furnace, radiant panel).

This BC represents incident radiation from an external source. Unlike `HeatFluxBC`
which applies flux directly at the surface, `RadiativeFluxBC` accounts for:
1. Surface absorptivity: fraction (absorptivity) enters material, rest is reflected
2. Radiation distribution: controlled by `radiation_model` parameter in `solve()`

The absorbed flux is: q_absorbed = absorptivity × incident_flux

The distribution of absorbed radiation depends on the `radiation_model` passed to `solve()`:
- `SURFACE_ABSORPTION`: All absorbed radiation deposited in surface cell (opaque materials)
- `BEER_LAMBERT`: In-depth volumetric absorption via Beer-Lambert law (semi-transparent)
- `P1_QUASI_STEADY`: P1 radiative transfer with volumetric emission

# Fields
- `flux::F`: Incident radiative flux [W/m²] (constant or time-dependent)
- `absorptivity::A`: Surface absorptivity (0-1). Can be:
  - `nothing` (default): Use effective emissivity from material properties (ε_eff = Σⱼ vⱼ εⱼ).
	By Kirchhoff's law, absorptivity = emissivity at thermal equilibrium.
  - A constant value (e.g., 0.9)
  - A time-dependent function (e.g., t -> 0.9 + 0.05*t/100)

# Example
```julia
# Use effective emissivity as absorptivity (recommended)
bc = RadiativeFluxBC(flux=25e3)

# Override with fixed absorptivity
bc = RadiativeFluxBC(flux=50e3, absorptivity=0.9)

# Time-varying flux with dynamic absorptivity
bc = RadiativeFluxBC(flux=t -> 25e3 * (1 - exp(-t/10)))

# Control radiation distribution via solve()
solve(problem; radiation_model=SURFACE_ABSORPTION)  # Opaque material
solve(problem; radiation_model=BEER_LAMBERT)        # Semi-transparent
```

# Note
When absorptivity=nothing (default), the effective emissivity is used as absorptivity
(Kirchhoff's law). This automatically tracks composition changes as material decomposes.

Do NOT combine `RadiativeFluxBC` with `HeatFluxBC` for the same radiation source -
this would double-count the energy input.
"""
struct RadiativeFluxBC{F, A} <: AbstractThermalBC
	flux::F           # Incident radiative flux [W/m²]
	absorptivity::A   # nothing = use effective, or constant, or function
end

# Constructors with optional absorptivity
function RadiativeFluxBC(; flux, absorptivity = nothing, validate::Bool = true)
	if validate && absorptivity !== nothing
		α_val = _get_scalar_value(absorptivity)
		validate_bc_params(α = α_val)
	end
	RadiativeFluxBC(flux, absorptivity)
end
function RadiativeFluxBC(flux; absorptivity = nothing, validate::Bool = true)
	if validate && absorptivity !== nothing
		α_val = _get_scalar_value(absorptivity)
		validate_bc_params(α = α_val)
	end
	RadiativeFluxBC(flux, absorptivity)
end

"""
	FlameBC{Fq,Fc} <: AbstractThermalBC

Flame heat flux boundary condition based on critical mass flux ignition model.

This BC models flame heat feedback to a pyrolyzing surface using the Combustion
Index (CI) approach. The flame turns on when sufficient fuel mass flux is generated
and extinguishes when the mass flux drops below a threshold.

# Fields
- `q_flame_max::Fq`: Maximum flame heat flux [W/m²] at full combustion
- `critical_mass_fluxes::Fc`: Dict mapping component index to critical mass flux [kg/(m²·s)]
  for ignition. Components not in this dict are ignored (non-combustible).
- `extinction_criterion::Float64`: CI threshold for extinction. Default 0.5.
  Use values < 1.0 to create hysteresis (harder to extinguish than ignite).
- `ramp_width::Float64`: Width of the sigmoid ramp for gradual flame onset. Default 0.0.
  - 0.0: Binary on/off at CI = 1.0 (original behavior)
  - 0.1-0.3: Gradual transition over CI range (more realistic)
  The ramp is centered at CI = 1.0 with width controlling steepness.

# Combustion Index Model
The Combustion Index is defined as:
	CI = Σⱼ (Jⱼ / Jⱼ,crit)
where Jⱼ is the mass flux of combustible component j and Jⱼ,crit is its
critical mass flux for ignition.

# Flame Heat Flux Model
With gradual onset (ramp_width > 0):
	q_flame = q_flame_max × σ((CI - 1.0) / ramp_width)
where σ(x) = 1/(1 + exp(-x)) is the sigmoid function.

With binary mode (ramp_width = 0):
	q_flame = q_flame_max if CI ≥ 1.0, else 0

# Hysteresis
When extinction_criterion < 1.0, the flame exhibits hysteresis:
- Flame ignites when CI first exceeds 1.0
- Once ignited, flame stays on until CI drops below extinction_criterion
- This prevents rapid on/off cycling near the ignition threshold

# Example
```julia
# Simple binary flame model
bc = FlameBC(
	q_flame_max = 25e3,  # 25 kW/m² flame heat flux
	critical_mass_fluxes = Dict(2 => 3e-3),  # Component 2 at 3 g/(m²·s)
)

# Gradual onset with hysteresis
bc = FlameBC(
	q_flame_max = 25e3,
	critical_mass_fluxes = Dict(2 => 3e-3, 3 => 5e-3),  # Multiple fuels
	extinction_criterion = 0.5,  # Extinguish at CI < 0.5
	ramp_width = 0.2,  # Gradual transition
)
```

# Note
The flame model assumes the gas-phase combustion is fast compared to the
pyrolysis rate, so the flame heat flux depends only on the instantaneous
mass flux, not on accumulated fuel or gas-phase transport.

# Thread Safety
As of v2.0, FlameBC is stateless - the flame on/off state is tracked externally
in the solver (PyrolysisODEFunction.flame_state). This enables thread-safe
operation for parallel ensemble solvers.
"""
struct FlameBC{Fq, Fc} <: AbstractThermalBC
	q_flame_max::Fq              # maximum flame heat flux [W/m²]
	critical_mass_fluxes::Fc     # Dict: component index => J_crit [kg/(m²·s)]
	extinction_criterion::Float64 # CI < this => extinction (default 0.5)
	ramp_width::Float64           # Sigmoid ramp width for gradual onset (default 0.0)
	# NOTE: flame_active state is now tracked externally in PyrolysisODEFunction
end

# Constructors
function FlameBC(; q_flame_max, critical_mass_fluxes,
	extinction_criterion::Float64 = 0.5,
	ramp_width::Float64 = 0.0)
	FlameBC(q_flame_max, critical_mass_fluxes, extinction_criterion, ramp_width)
end

# Backward compatible constructor (positional arguments)
FlameBC(q_flame_max, critical_mass_fluxes, extinction_criterion::Float64) =
	FlameBC(q_flame_max, critical_mass_fluxes, extinction_criterion, 0.0)

# ============================================================================
# Mass Boundary Conditions
# ============================================================================

"""
	ImpermeableBC <: AbstractMassBC

Zero mass flux (impermeable) boundary condition.
"""
struct ImpermeableBC <: AbstractMassBC end

"""
	MassFluxBC{F} <: AbstractMassBC

Prescribed mass flux boundary condition.

# Fields
- `flux::F`: Mass flux [kg/(m²·s)] for each gas component
"""
struct MassFluxBC{F} <: AbstractMassBC
	flux::F  # J(t) for each gas component
end

"""
	ConvectiveMassBC{Fh,Fξ} <: AbstractMassBC

Convective mass transfer boundary condition for gas escape at boundaries.

The mass flux is computed as:
	J = h_m × (ξ_inf - ξ_surface) [kg/(m²·s)]

where ξ is mass concentration [kg/m³]. This represents external convective
mass transfer where gas is carried away by the external flow.

# Fields
- `h_m::Fh`: External mass transfer coefficient [m/s]. This is a velocity-based
  coefficient representing the convective mass transfer rate. Typical values:
  - Natural convection: 0.01 - 0.05 m/s
  - Forced convection (low velocity): 0.05 - 0.1 m/s
  - Forced convection (high velocity): 0.1 - 0.5 m/s
- `ξ_inf::Fξ`: Ambient mass concentration [kg/m³]. Typically 0 for pyrolysis
  gases escaping to atmosphere.

# Physical Interpretation
The coefficient h_m can be estimated from heat-mass transfer analogy:
	h_m ≈ h / (ρ_air × cp_air)
where h is the convective heat transfer coefficient [W/(m²·K)].

For a typical cone calorimeter with h ≈ 10 W/(m²·K):
	h_m ≈ 10 / (1.2 × 1000) ≈ 0.008 m/s

# Example
```julia
bc = ConvectiveMassBC(h_m=0.05, ξ_inf=0.0)  # Gas escapes to atmosphere
bc = ConvectiveMassBC(h_m=t -> 0.05*(1 + 0.1*t), ξ_inf=0.0)  # Time-varying
```
"""
struct ConvectiveMassBC{Fh, Fξ} <: AbstractMassBC
	h_m::Fh         # mass transfer coefficient [m/s]
	ξ_inf::Fξ       # ambient concentration [kg/m³]
end

function ConvectiveMassBC(; h_m, ξ_inf = 0.0, validate::Bool = true)
	if validate
		h_m_val = _get_scalar_value(h_m)
		validate_bc_params(h_m = h_m_val)
	end
	ConvectiveMassBC(h_m, ξ_inf)
end

"""
	DiffusiveMassBC{Fh} <: AbstractMassBC

Diffusive mass transfer boundary condition using the internal diffusion model.

The flux is computed as: J = λ × ξ / d × h_m_factor

where:
- λ = gas transfer coefficient from material [m²/s]
- ξ = gas concentration at surface cell [kg/m³]
- d = distance from cell center to boundary [m]
- h_m_factor = external resistance factor (default Inf = no external resistance)

For combined internal + external resistance:
  J = ξ / (d/λ + 1/h_m) = λ × h_m × ξ / (λ + h_m × d)

When h_m → ∞: J = λ × ξ / d (pure diffusion, default)
When h_m → 0: J = 0 (impermeable)

This BC ensures consistency between species transport and ALE regression.

# Fields
- `h_m::Fh`: External mass transfer coefficient [m/s] (Inf = no external resistance)
"""
struct DiffusiveMassBC{Fh} <: AbstractMassBC
	h_m::Fh  # external mass transfer coefficient (Inf = no resistance)
end

# Convenience constructors
DiffusiveMassBC() = DiffusiveMassBC{Float64}(Inf)

"""
	CombinedMassBC{BCs<:Tuple} <: AbstractMassBC

Combined mass boundary condition (sum of multiple BCs).

This mirrors `CombinedThermalBC` and allows combining multiple mass BCs
at a single boundary. The total mass flux is the sum of individual BC fluxes.

# Fields
- `bcs::BCs`: Tuple of boundary conditions

# Example
```julia
# Combine diffusive and convective mass transfer
bc = DiffusiveMassBC(Inf) + ConvectiveMassBC(h_m=0.01, ξ_inf=0.0)
```

# Note
For most pyrolysis applications, a single mass BC (usually `DiffusiveMassBC`
or `ConvectiveMassBC`) is sufficient. Combined mass BCs are useful for:
- Modeling competing transport mechanisms
- Complex surface resistance models
- Research/experimental validation
"""
struct CombinedMassBC{BCs <: Tuple} <: AbstractMassBC
	bcs::BCs
end

# Combine BCs with + operator
Base.:+(bc1::AbstractMassBC, bc2::AbstractMassBC) =
	CombinedMassBC((bc1, bc2))
Base.:+(cbc::CombinedMassBC, bc::AbstractMassBC) =
	CombinedMassBC((cbc.bcs..., bc))
Base.:+(bc::AbstractMassBC, cbc::CombinedMassBC) =
	CombinedMassBC((bc, cbc.bcs...))

# ============================================================================
# Pressure Boundary Conditions (for Darcy Flow)
# ============================================================================

"""
	DirichletPressureBC <: AbstractMassBC

Fixed pressure at boundary (typically atmospheric at exposed surface).

When Darcy flow is enabled, this BC computes the gas velocity from the 
pressure gradient between the boundary and the interior cell:
	v_g = -(κ/μ) × (P_bc - P_interior) / d

where d is the distance from the cell center to the boundary.

# Fields
- `P::Float64`: Boundary pressure [Pa]

# Example
```julia
bc = DirichletPressureBC(101325.0)  # Atmospheric pressure
bc = DirichletPressureBC()          # Default: 101325 Pa
```
"""
struct DirichletPressureBC <: AbstractMassBC
	P::Float64  # Pressure [Pa]
end

# Constructors with optional validation
function DirichletPressureBC(P::Real; validate::Bool = true)
	if validate
		validate_bc_params(P = P)
	end
	DirichletPressureBC(Float64(P))
end

# Default constructor: standard atmospheric pressure
DirichletPressureBC() = DirichletPressureBC(101325.0)

"""
	NeumannPressureBC <: AbstractMassBC

Zero or prescribed pressure gradient at boundary (impermeable or specified flux).

For Darcy flow:
- ∂P/∂n = 0 (default): Zero gas velocity (impermeable)
- ∂P/∂n ≠ 0: Gas velocity v_g = -(κ/μ) × ∂P/∂n

# Fields
- `∂P∂n::Float64`: Pressure gradient normal to boundary [Pa/m] (default: 0 = impermeable)

# Example
```julia
bc = NeumannPressureBC()      # Impermeable (zero gradient)
bc = NeumannPressureBC(0.0)   # Same as above
bc = NeumannPressureBC(1000.0)  # Prescribed pressure gradient
```
"""
struct NeumannPressureBC <: AbstractMassBC
	∂P∂n::Float64  # Pressure gradient [Pa/m]
end

# Default constructor: zero gradient (impermeable)
NeumannPressureBC() = NeumannPressureBC(0.0)

"""
	ConvectivePressureBC <: AbstractMassBC

Pressure-driven outflow with resistance (Robin-type BC for pressure).

The gas velocity at the boundary is computed as:
	v_g·n = h_P × (P_internal - P_external)

This models a surface with finite resistance to gas flow, such as:
- Surface char layer with limited permeability
- Protective coating or membrane
- Blowing resistance at high mass flux

# Fields
- `h_P::Float64`: Pressure transfer coefficient [m/(Pa·s)]
- `P_ext::Float64`: External pressure [Pa]

# Example
```julia
bc = ConvectivePressureBC(1e-6, 101325.0)  # Typical resistance, atmospheric exterior
```

# Note
h_P → ∞ approaches DirichletPressureBC behavior
h_P → 0 approaches ImpermeableBC behavior
"""
struct ConvectivePressureBC <: AbstractMassBC
	h_P::Float64     # Pressure transfer coefficient [m/(Pa·s)]
	P_ext::Float64   # External pressure [Pa]
end

# ============================================================================
# Mesh Boundary Conditions (for ALE)
# ============================================================================

"""
	FixedMeshBC <: AbstractMeshBC

Fixed mesh position (no mesh motion at boundary).
"""
struct FixedMeshBC <: AbstractMeshBC end

"""
	FreeSurfaceBC <: AbstractMeshBC

Free surface that follows material motion (Lagrangian).
"""
struct FreeSurfaceBC <: AbstractMeshBC end

"""
	PrescribedVelocityBC{F} <: AbstractMeshBC

Prescribed mesh velocity at boundary.

# Fields
- `velocity::F`: Mesh velocity [m/s]
"""
struct PrescribedVelocityBC{F} <: AbstractMeshBC
	velocity::F  # w(t)
end

"""
	SlidingBC <: AbstractMeshBC

Node can slide along boundary (constrained motion).

# Fields
- `direction::Symbol`: Allowed motion direction (:r or :z)
"""
struct SlidingBC <: AbstractMeshBC
	direction::Symbol  # :r or :z
end

# ============================================================================
# Type-Stable Boundary Condition Set
# ============================================================================

"""
	BoundaryConditionSet{TB1, TB2, MB1, MB2}

Type-stable container for boundary conditions.

This struct provides compile-time type information for all boundary conditions,
enabling static dispatch and eliminating the overhead of Dict-based BC storage.

For 1D problems:
- Boundary 1 (top): z = L, exposed surface
- Boundary 2 (bottom): z = 0, substrate/back face

# Type Parameters
- `TB1`: Thermal BC type at top boundary
- `TB2`: Thermal BC type at bottom boundary
- `MB1`: Mass BC type at top boundary
- `MB2`: Mass BC type at bottom boundary

# Fields
- `thermal_top::TB1`: Thermal BC at top (exposed) surface
- `thermal_bottom::TB2`: Thermal BC at bottom (substrate) surface
- `mass_top::MB1`: Mass BC at top surface
- `mass_bottom::MB2`: Mass BC at bottom surface

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε=0.9, T_inf=300),
	thermal_bottom = AdiabaticBC(),
	mass_top = DiffusiveMassBC(Inf),
	mass_bottom = ImpermeableBC()
)
```

# Notes
- All type parameters are concrete, enabling full specialization
- Combined BCs (via `+` operator) are supported
- Default values are `AdiabaticBC()` for thermal and `ImpermeableBC()` for mass
"""
struct BoundaryConditionSet{TB1 <: AbstractThermalBC, TB2 <: AbstractThermalBC,
	MB1 <: AbstractMassBC, MB2 <: AbstractMassBC}
	thermal_top::TB1
	thermal_bottom::TB2
	mass_top::MB1
	mass_bottom::MB2
end

# Default constructor with keyword arguments
"""
	BoundaryConditionSet(; thermal_top, thermal_bottom, mass_top, mass_bottom)

Construct a BoundaryConditionSet with keyword arguments.

# Keyword Arguments
- `thermal_top::AbstractThermalBC = AdiabaticBC()`: Thermal BC at top (exposed) surface
- `thermal_bottom::AbstractThermalBC = AdiabaticBC()`: Thermal BC at bottom surface
- `mass_top::AbstractMassBC = ImpermeableBC()`: Mass BC at top surface
- `mass_bottom::AbstractMassBC = ImpermeableBC()`: Mass BC at bottom surface

# Example
```julia
# Minimal: only specify what you need
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3)
)

# Full specification
bc_set = BoundaryConditionSet(
	thermal_top = HeatFluxBC(50e3) + RadiativeBC(ε=0.9, T_inf=300),
	thermal_bottom = AdiabaticBC(),
	mass_top = DiffusiveMassBC(Inf),
	mass_bottom = ImpermeableBC()
)
```
"""
function BoundaryConditionSet(;
	thermal_top::AbstractThermalBC = AdiabaticBC(),
	thermal_bottom::AbstractThermalBC = AdiabaticBC(),
	mass_top::AbstractMassBC = ImpermeableBC(),
	mass_bottom::AbstractMassBC = ImpermeableBC())
	BoundaryConditionSet(thermal_top, thermal_bottom, mass_top, mass_bottom)
end

# Accessor functions using Val for compile-time dispatch
"""
	get_thermal_bc(bc_set, ::Val{1}) -> thermal_top
	get_thermal_bc(bc_set, ::Val{2}) -> thermal_bottom

Get thermal BC for a boundary using compile-time boundary ID.
"""
@inline get_thermal_bc(bc_set::BoundaryConditionSet, ::Val{1}) = bc_set.thermal_top
@inline get_thermal_bc(bc_set::BoundaryConditionSet, ::Val{2}) = bc_set.thermal_bottom

"""
	get_mass_bc(bc_set, ::Val{1}) -> mass_top
	get_mass_bc(bc_set, ::Val{2}) -> mass_bottom

Get mass BC for a boundary using compile-time boundary ID.
"""
@inline get_mass_bc(bc_set::BoundaryConditionSet, ::Val{1}) = bc_set.mass_top
@inline get_mass_bc(bc_set::BoundaryConditionSet, ::Val{2}) = bc_set.mass_bottom

# Runtime boundary ID versions (for compatibility, less efficient)
"""
	get_thermal_bc(bc_set, boundary_id::Int) -> AbstractThermalBC

Get thermal BC for a boundary using runtime boundary ID.
Less efficient than Val version - use Val{1} or Val{2} in hot paths.
"""
@inline function get_thermal_bc(bc_set::BoundaryConditionSet, boundary_id::Int)
	boundary_id == 1 && return bc_set.thermal_top
	boundary_id == 2 && return bc_set.thermal_bottom
	error("Invalid boundary_id: $boundary_id (expected 1 or 2)")
end

"""
	get_mass_bc(bc_set, boundary_id::Int) -> AbstractMassBC

Get mass BC for a boundary using runtime boundary ID.
Less efficient than Val version - use Val{1} or Val{2} in hot paths.
"""
@inline function get_mass_bc(bc_set::BoundaryConditionSet, boundary_id::Int)
	boundary_id == 1 && return bc_set.mass_top
	boundary_id == 2 && return bc_set.mass_bottom
	error("Invalid boundary_id: $boundary_id (expected 1 or 2)")
end

# ============================================================================
# BC Evaluation Helpers (stubs - full implementation in Discretization)
# ============================================================================

"""
	evaluate_bc(bc, t, T_surface, state...)

Evaluate a boundary condition to get the heat/mass flux.
Returns flux value (positive into domain).
"""
evaluate_bc(bc::AdiabaticBC, t, T_surface, args...) = 0.0
evaluate_bc(bc::HeatFluxBC{<:Real}, t, T_surface, args...) = bc.flux
evaluate_bc(bc::HeatFluxBC{<:Function}, t, T_surface, args...) = bc.flux(t)
evaluate_bc(bc::ImpermeableBC, t, args...) = 0.0

# ============================================================================
# RadiativeFluxBC Evaluation Helpers
# ============================================================================

"""
	get_radiative_flux_bc(bc::AbstractThermalBC) -> Union{RadiativeFluxBC, Nothing}

Extract RadiativeFluxBC from a (possibly combined) thermal BC.
Returns nothing if no RadiativeFluxBC is present.
"""
get_radiative_flux_bc(bc::RadiativeFluxBC) = bc
get_radiative_flux_bc(bc::AbstractThermalBC) = nothing

function get_radiative_flux_bc(bc::CombinedThermalBC)
	for inner_bc in bc.bcs
		result = get_radiative_flux_bc(inner_bc)
		if result !== nothing
			return result
		end
	end
	return nothing
end

# ============================================================================
# Type-Stable BoundaryConditionSet Signatures
# ============================================================================

"""
	extract_radiative_flux_from_bcs(bc_set::BoundaryConditionSet, ::Val{B}, t::Real) where B -> Tuple{Float64, Bool}

Extract the incident radiation intensity from thermal BCs using type-stable BoundaryConditionSet.

This is the preferred signature for maximum type stability. Uses compile-time
boundary ID for dispatch.

# Arguments
- `bc_set`: Type-stable BoundaryConditionSet
- `::Val{B}`: Compile-time boundary ID (Val(1) for top, Val(2) for bottom)
- `t`: Current time [s]

# Returns
(I_incident, in_depth) where:
- I_incident: incident radiation intensity at surface [W/m²]
  (absorptivity is applied later in compute_radiation_source! using ε_eff from prop_cache)
- in_depth: true if radiation should be distributed via Beer-Lambert

# Note
The absorptivity is NOT applied here. For SURFACE_ABSORPTION and BEER_LAMBERT modes,
the effective emissivity from the material properties is used in compute_radiation_source!.
"""
function extract_radiative_flux_from_bcs(bc_set::BoundaryConditionSet, ::Val{B}, t::Real) where B
	bc = get_thermal_bc(bc_set, Val(B))

	rad_bc = get_radiative_flux_bc(bc)
	if rad_bc === nothing
		return (0.0, true)
	end

	# Evaluate flux (may be time-dependent)
	flux = rad_bc.flux isa Function ? rad_bc.flux(t) : Float64(rad_bc.flux)

	# Return INCIDENT flux - absorptivity is applied in compute_radiation_source!
	# using the effective emissivity from the surface cell's properties.
	# If explicit absorptivity is specified in RadiativeFluxBC, apply it here.
	if rad_bc.absorptivity !== nothing
		α = rad_bc.absorptivity isa Function ? rad_bc.absorptivity(t) : Float64(rad_bc.absorptivity)
		return (flux * α, false)
	else
		# No explicit absorptivity - return incident flux, let ε_eff be applied later
		return (flux, false)
	end
end

"""
	extract_radiative_flux_from_bcs(bc_set::BoundaryConditionSet, boundary_id::Int, t::Real) -> Tuple{Float64, Bool}

Extract the incident radiation intensity using BoundaryConditionSet with runtime boundary_id.

Less efficient than Val version - use Val{1} or Val{2} in hot paths.
"""
function extract_radiative_flux_from_bcs(bc_set::BoundaryConditionSet, boundary_id::Int, t::Real)
	bc = get_thermal_bc(bc_set, boundary_id)

	rad_bc = get_radiative_flux_bc(bc)
	if rad_bc === nothing
		return (0.0, true)
	end

	flux = rad_bc.flux isa Function ? rad_bc.flux(t) : Float64(rad_bc.flux)

	# Return INCIDENT flux - absorptivity is applied in compute_radiation_source!
	# using the effective emissivity from the surface cell's properties.
	# If explicit absorptivity is specified in RadiativeFluxBC, apply it here.
	if rad_bc.absorptivity !== nothing
		α = rad_bc.absorptivity isa Function ? rad_bc.absorptivity(t) : Float64(rad_bc.absorptivity)
		return (flux * α, false)
	else
		# No explicit absorptivity - return incident flux, let ε_eff be applied later
		return (flux, false)
	end
end

"""
	has_radiative_flux_bc(bc_set::BoundaryConditionSet, ::Val{B}) where B -> Bool

Check if a boundary has a RadiativeFluxBC using type-stable BoundaryConditionSet.

Uses compile-time boundary ID for dispatch.
"""
function has_radiative_flux_bc(bc_set::BoundaryConditionSet, ::Val{B}) where B
	bc = get_thermal_bc(bc_set, Val(B))
	return get_radiative_flux_bc(bc) !== nothing
end

"""
	has_radiative_flux_bc(bc_set::BoundaryConditionSet, boundary_id::Int) -> Bool

Check if a boundary has a RadiativeFluxBC using BoundaryConditionSet with runtime boundary_id.

Less efficient than Val version - use Val{1} or Val{2} in hot paths.
"""
function has_radiative_flux_bc(bc_set::BoundaryConditionSet, boundary_id::Int)
	bc = get_thermal_bc(bc_set, boundary_id)
	return get_radiative_flux_bc(bc) !== nothing
end