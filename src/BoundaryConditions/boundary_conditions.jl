"""
Boundary condition type definitions.
"""

# ============================================================================
# Validation Helpers
# ============================================================================

"""
	validate_bc_params(; kwargs...)

Validate boundary condition parameters, as the BC constructors do.

# Supported validations
- `h`: Heat transfer coefficient must be non-negative
- `T`: Temperature must be positive (Kelvin)
- `ε`: Emissivity must be in [0, 1]
- `F`: View factor must be in [0, 1]
- `α`: Absorptivity must be in [0, 1]

Throws `ArgumentError` if validation fails.
"""
function validate_bc_params(;
    h::Union{Real,Nothing} = nothing,
    T::Union{Real,Nothing} = nothing,
    ε::Union{Real,Nothing} = nothing,
    F::Union{Real,Nothing} = nothing,
    α::Union{Real,Nothing} = nothing,
)

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

"""
	ConvectiveBC{Fh,FT,BM} <: AbstractThermalBC

Convective heat transfer boundary condition: q" = h(T∞ - T_surface)

# Fields
- `h::Fh`: Heat transfer coefficient without outflow [W/(m²·K)]
- `T_inf::FT`: Ambient temperature [K]
- `blowing::BM`: how gas leaving through this face changes `h`, an
  `AbstractBlowingModel`. `NoBlowing()` (default) leaves `h` alone;
  `CouetteFilm()` is the classical film law; `PropertyRatioFilm` and
  `PoolingFilm` take the gas and the free stream into account. The face's
  outward gas flux is read from the mass balance at every surface iterate.

# Example
```julia
bc = ConvectiveBC(h=10.0, T_inf=300.0)
bc = ConvectiveBC(h=10.0, T_inf=t -> 300 + 0.1*t)  # Warming ambient
bc = ConvectiveBC(h=11.2, T_inf=360.0, blowing=PoolingFilm(NITROGEN))
```
"""
struct ConvectiveBC{Fh,FT,BM<:AbstractBlowingModel} <: AbstractThermalBC
    h::Fh                    # h(t) or constant
    T_inf::FT                # T∞(t) or constant
    blowing::BM              # blowing model
end

function ConvectiveBC(; h, T_inf, blowing::AbstractBlowingModel = NoBlowing())
    validate_bc_params(h = _get_scalar_value(h), T = _get_scalar_value(T_inf))
    return ConvectiveBC(h, T_inf, blowing)
end

"""
	RadiativeBC{Fε,FT,FF} <: AbstractThermalBC

Radiative exchange with gray surroundings: q" = ε σ (F T∞⁴ − T_surface⁴).

The surface emits ε σ T_surface⁴ into its whole hemisphere whatever it faces,
and receives ε σ T∞⁴ from the fraction F of that hemisphere occupied by the
surroundings at T∞. The flux is positive into the material when the
irradiation exceeds the emission and negative (a re-radiation loss) once the
surface is hot.

# Fields
- `ε::Fε`: Surface emissivity [-], range 0-1. Can be:
  - `nothing` (default): Use effective emissivity from material properties (ε_eff = Σⱼ vⱼ εⱼ),
	which tracks composition changes.
  - A constant value (e.g., 0.9)
  - A temperature-dependent function (e.g., T -> 0.85 + 0.1*(T-300)/700)
  Typical values:
  - Polymers: 0.85 - 0.95
  - Char/carbon: 0.90 - 0.98
  - Polished metals: 0.02 - 0.10
  - Oxidized metals: 0.40 - 0.80
- `T_inf::FT`: Ambient/surroundings temperature [K]. Can be constant or time-dependent.
- `F::FF`: View factor [-] from the surface to the surroundings at `T_inf`,
  range 0-1, default 1.0. It scales the irradiation received from them only,
  never the surface's own emission. Use `F < 1` when the rest of the
  hemisphere is taken up by a source whose irradiation is supplied
  separately (a heater through `RadiativeFluxBC`, with its own view factor
  folded into that flux). A cone flux set with a water-cooled gauge already
  carries the heater's geometry, so a cone calorimeter is modelled with
  `F = 1`.

# Example
```julia
# Use effective emissivity from material (default)
bc = RadiativeBC(T_inf=300.0)

# Override with fixed emissivity
bc = RadiativeBC(ε=0.9, T_inf=300.0)

# Surroundings over 40 % of the hemisphere; the heater filling the other
# 60 % is supplied through a RadiativeFluxBC
bc = RadiativeBC(T_inf=300.0, F=0.4)

# Temperature-dependent emissivity (e.g., char formation)
bc = RadiativeBC(ε=T -> 0.85 + 0.1*(T-300)/700, T_inf=300.0)
```

# Note
When ε=nothing (default), the effective emissivity is computed as the volume-weighted
average of component emissivities: ε_eff = Σⱼ vⱼ εⱼ. This automatically tracks the
change in surface properties as material decomposes (e.g., polymer → char).
"""
struct RadiativeBC{Fε,FT,FF} <: AbstractThermalBC
    ε::Fε           # nothing = use effective, or ε(T), or constant
    T_inf::FT       # T∞(t) or constant
    F::FF           # View factor to the surroundings (0-1), default 1.0
end

function RadiativeBC(; ε = nothing, T_inf, F = 1.0)
    validate_bc_params(
        ε = _get_scalar_value(ε),
        T = _get_scalar_value(T_inf),
        F = _get_scalar_value(F),
    )
    return RadiativeBC(ε, T_inf, F)
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
struct CombinedThermalBC{BCs<:Tuple} <: AbstractThermalBC
    bcs::BCs
end

# Combine BCs with + operator
Base.:+(bc1::AbstractThermalBC, bc2::AbstractThermalBC) = CombinedThermalBC((bc1, bc2))
Base.:+(cbc::CombinedThermalBC, bc::AbstractThermalBC) = CombinedThermalBC((cbc.bcs..., bc))
Base.:+(bc::AbstractThermalBC, cbc::CombinedThermalBC) = CombinedThermalBC((bc, cbc.bcs...))
# Disambiguate cbc + cbc (both prepend and append methods would otherwise match).
Base.:+(cbc1::CombinedThermalBC, cbc2::CombinedThermalBC) =
    CombinedThermalBC((cbc1.bcs..., cbc2.bcs...))

"""
	AngularDistribution

How the radiation of a [`RadiativeFluxBC`](@ref) arrives at its face:
[`Diffuse`](@ref), [`Collimated`](@ref) or [`Cone`](@ref). It sets where the
radiation is absorbed under `RADIATIVE_TRANSFER` and has no effect under
`SURFACE_ABSORPTION`.
"""
abstract type AngularDistribution end

"""
	Diffuse()

Radiation arriving equally from every direction of the hemisphere above the
face. Under `RADIATIVE_TRANSFER` it crosses the face at the emissivity of the
cell beside it and travels on as hemispherical flux.
"""
struct Diffuse <: AngularDistribution end

"""
	Collimated()

Radiation arriving along the face normal, as one beam.
"""
struct Collimated <: AngularDistribution end

"""
	Cone(θ₁, θ₂; directions = 3)

Radiation of uniform intensity arriving between the polar angles `θ₁` and
`θ₂` from the face normal [rad], `0 ≤ θ₁ < θ₂ ≤ π/2`: the view of a heater
ring from the sample. It travels as `directions` beams, at the Gauss–Legendre
nodes of the direction cosine `μ` over `[cos θ₂, cos θ₁]`, each carrying the
share of the flux that its weight times its `μ` gives.

# Example
```julia
Cone(deg2rad(20), deg2rad(73))   # a cone heater seen from the sample
```
"""
struct Cone <: AngularDistribution
    θ₁::Float64
    θ₂::Float64
    directions::Int
    function Cone(θ₁::Real, θ₂::Real, directions::Integer)
        0 <= θ₁ < θ₂ <= π / 2 || throw(
            ArgumentError("a cone needs 0 ≤ θ₁ < θ₂ ≤ π/2 [rad], got θ₁ = $θ₁, θ₂ = $θ₂"),
        )
        directions >= 1 ||
            throw(ArgumentError("a cone needs at least one direction, got $directions"))
        return new(θ₁, θ₂, directions)
    end
end

Cone(θ₁::Real, θ₂::Real; directions::Integer = 3) = Cone(θ₁, θ₂, directions)

"""
	RadiativeFluxBC{F,A,D} <: AbstractThermalBC

External radiative heat source (e.g., cone heater, furnace, radiant panel).

This BC represents incident radiation from an external source. Unlike `HeatFluxBC`
which applies flux directly at the surface, `RadiativeFluxBC` accounts for:
1. Surface absorptivity: fraction (absorptivity) enters material, rest is reflected
2. Radiation distribution: controlled by `radiation_model` parameter in `solve()`

Where the absorbed radiation goes depends on the `radiation_model` passed to
`solve()`:
- `SURFACE_ABSORPTION`: `absorptivity × flux` is absorbed at the exposed face,
  through the surface energy balance (opaque materials), whatever the
  `distribution`
- `RADIATIVE_TRANSFER`: the radiation crosses the face and is absorbed in depth.
  [`Diffuse`](@ref) radiation crosses the exposed face at the emissivity of the
  exposed cell, so it takes no `absorptivity`; [`Collimated`](@ref) and
  [`Cone`](@ref) radiation enters as beams at the absorptivity. Only diffuse
  radiation may arrive at `:bottom`.

# Fields
- `flux::F`: Incident radiative flux [W/m²] (constant or time-dependent)
- `absorptivity::A`: Surface absorptivity (0-1). Can be:
  - `nothing` (default): Use the material's effective surface absorptivity
	`a_eff = Σⱼ vⱼ aⱼ`, where `aⱼ` is the component's `absorptivity` keyword
	when set and its emissivity `εⱼ(T)` otherwise. With no component
	absorptivity set this is Kirchhoff's law, `a_eff = ε_eff`; setting one
	decouples the absorbed heater flux from the temperature-dependent
	emission (see `SolidComponent`).
  - A constant value (e.g., 0.9)
  - A time-dependent function (e.g., t -> 0.9 + 0.05*t/100)
- `distribution::D`: how the radiation arrives, an [`AngularDistribution`](@ref):
  `Diffuse()` (default), `Collimated()` or `Cone(θ₁, θ₂)`.

# Example
```julia
# Use the material's surface absorptivity (ε_eff by default)
bc = RadiativeFluxBC(25e3)

# Override with fixed absorptivity
bc = RadiativeFluxBC(50e3, absorptivity=0.9)

# Time-varying flux with dynamic absorptivity
bc = RadiativeFluxBC(t -> 25e3 * (1 - exp(-t/10)))

# A cone heater, absorbed in depth under RADIATIVE_TRANSFER
bc = RadiativeFluxBC(25e3; distribution = Cone(deg2rad(20), deg2rad(73)))
```

# Note
When absorptivity=nothing (default), the material's mixture surface absorptivity
`a_eff` is used — identical to the effective emissivity (Kirchhoff's law) unless a
component carries its own `absorptivity`. Either way it automatically tracks
composition and temperature changes as material decomposes.

Do NOT combine `RadiativeFluxBC` with `HeatFluxBC` for the same radiation source -
this would double-count the energy input.
"""
struct RadiativeFluxBC{F,A,D<:AngularDistribution} <: AbstractThermalBC
    flux::F           # Incident radiative flux [W/m²]
    absorptivity::A   # nothing = use effective, or constant, or function
    distribution::D   # how the radiation arrives at the face
end

function RadiativeFluxBC(
    flux;
    absorptivity = nothing,
    distribution::AngularDistribution = Diffuse(),
)
    validate_bc_params(α = _get_scalar_value(absorptivity))
    return RadiativeFluxBC(flux, absorptivity, distribution)
end

# ============================================================================
# Mass Boundary Conditions
# ============================================================================

"""
	ImpermeableBC <: AbstractMassBC

Zero mass flux (impermeable) boundary condition.
"""
struct ImpermeableBC <: AbstractMassBC end

include("gas_boundaries.jl")

# ============================================================================
# Type-Stable Boundary Condition Set
# ============================================================================

"""
	BoundaryConditionSet{TT <: NamedTuple, MT <: NamedTuple}

Type-stable, **tag-keyed** container for boundary conditions. Each tag
(`:top` and `:bottom`, and the other tags of a mesh) maps to one
`AbstractThermalBC` and one `AbstractMassBC`. Internal storage is a typed
`NamedTuple` per BC kind, so symbol lookup compiles to a direct field read
and the residual specializes statically.

# Fields
- `thermal::TT`: `NamedTuple{tags, Tuple{<: AbstractThermalBC, ...}}`
- `mass::MT`: `NamedTuple{tags, Tuple{<: AbstractMassBC, ...}}`

# Tags
- `:top` ↔ exposed surface
- `:bottom` ↔ substrate

# Construction
Each BC kind takes a NamedTuple of `tag => BC` mappings, one per boundary tag of
the mesh; `:top` and `:bottom` default to `AdiabaticBC()` and `ImpermeableBC()`.

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300),
	           bottom = AdiabaticBC()),
	mass    = (top    = AmbientGasBC(GasReservoir(P=0., T=300.)),
	           bottom = ImpermeableBC()),
)
```

The conditions are read by tag: `bc_set.thermal.top`, `bc_set.mass.bottom`.
"""
struct BoundaryConditionSet{TT<:NamedTuple,MT<:NamedTuple}
    thermal::TT
    mass::MT
end

# The conditions of `:top` and `:bottom` unless given.
const _DEFAULT_THERMAL_BCS = (top = AdiabaticBC(), bottom = AdiabaticBC())
const _DEFAULT_MASS_BCS = (top = ImpermeableBC(), bottom = ImpermeableBC())

"""
	BoundaryConditionSet(; thermal::NamedTuple = (top = AdiabaticBC(),  bottom = AdiabaticBC()),
						   mass::NamedTuple   = (top = ImpermeableBC(), bottom = ImpermeableBC()))

The conditions by boundary tag: each kind (`thermal`, `mass`) takes a NamedTuple
of `tag = BC` entries, merged over the defaults of `:top` and `:bottom`.
"""
function BoundaryConditionSet(;
    thermal::NamedTuple = _DEFAULT_THERMAL_BCS,
    mass::NamedTuple = _DEFAULT_MASS_BCS,
)
    thermal_full = merge(_DEFAULT_THERMAL_BCS, thermal)
    mass_full = merge(_DEFAULT_MASS_BCS, mass)
    # Validate that thermal and mass payloads are AbstractThermalBC / AbstractMassBC
    # at every tag — fail loudly if a user passes the wrong BC kind to the wrong slot.
    for (tag, bc) in pairs(thermal_full)
        bc isa AbstractThermalBC || throw(
            ArgumentError("thermal[$tag] must be an AbstractThermalBC, got $(typeof(bc))"),
        )
    end
    for (tag, bc) in pairs(mass_full)
        bc isa AbstractMassBC ||
            throw(ArgumentError("mass[$tag] must be an AbstractMassBC, got $(typeof(bc))"))
    end
    BoundaryConditionSet(thermal_full, mass_full)
end

# ============================================================================
# Time Dependence
# ============================================================================

_any_function(parameters...) = any(x -> x isa Function, parameters)

"""
	reads_time(bc) -> Bool

Whether a boundary condition, or a [`BoundaryConditionSet`](@ref), can read the
time: whether a parameter the residual evaluates at the time is a function. A
function of the temperature, as `RadiativeBC`'s `ε` and a film's conductance are,
does not count. A condition with no method of its own is taken to read it.
"""
reads_time(::Union{AbstractThermalBC,AbstractMassBC}) = true
reads_time(::Union{AdiabaticBC,ImpermeableBC}) = false
reads_time(bc::HeatFluxBC) = _any_function(bc.flux)
reads_time(bc::ConvectiveBC) = _any_function(bc.h, bc.T_inf)
reads_time(bc::RadiativeBC) = _any_function(bc.T_inf, bc.F)
reads_time(bc::DirichletThermalBC) = _any_function(bc.T)
reads_time(bc::RadiativeFluxBC) = _any_function(bc.flux, bc.absorptivity)
reads_time(bc::CombinedThermalBC) = any(reads_time, bc.bcs)
reads_time(res::GasReservoir) = _any_function(res.P, res.T, res.mole_fractions)
reads_time(bc::AmbientGasBC) = reads_time(bc.reservoir)
reads_time(bc::ConvectivePressureBC) =
    reads_time(bc.reservoir) || _any_function(bc.resistance)
reads_time(bc::NeumannPressureBC) = _any_function(bc.gradient, bc.T, bc.mole_fractions)
reads_time(bc::MassFluxBC) =
    _any_function(bc.fluxes...) || _any_function(bc.inflow_temperature...)
reads_time(set::BoundaryConditionSet) =
    any(reads_time, values(set.thermal)) || any(reads_time, values(set.mass))

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
