"""
Boundary condition type definitions.
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

# Convenience constructor to convert integers to Float64
HeatFluxBC(flux::Integer) = HeatFluxBC{Float64}(Float64(flux))

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

ConvectiveBC(h, T_inf) = ConvectiveBC(h, T_inf, NoBlowing())

function ConvectiveBC(;
    h,
    T_inf,
    blowing::AbstractBlowingModel = NoBlowing(),
    validate::Bool = true,
)
    if validate
        h_val = _get_scalar_value(h)
        T_val = _get_scalar_value(T_inf)
        validate_bc_params(h = h_val, T = T_val)
    end
    ConvectiveBC(h, T_inf, blowing)
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
struct RadiativeBC{Fε,FT,FF} <: AbstractThermalBC
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
	RadiativeFluxBC{F,A} <: AbstractThermalBC

External radiative heat source (e.g., cone heater, furnace, radiant panel).

This BC represents incident radiation from an external source. Unlike `HeatFluxBC`
which applies flux directly at the surface, `RadiativeFluxBC` accounts for:
1. Surface absorptivity: fraction (absorptivity) enters material, rest is reflected
2. Radiation distribution: controlled by `radiation_model` parameter in `solve()`

The absorbed flux is: q_absorbed = absorptivity × incident_flux

The distribution of absorbed radiation depends on the `radiation_model` passed to `solve()`:
- `NO_RADIATION`, `SURFACE_ABSORPTION`: absorbed at the exposed face, through the surface energy balance (opaque materials)
- `BEER_LAMBERT`: In-depth volumetric absorption via Beer-Lambert law (semi-transparent)

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

# Example
```julia
# Use the material's surface absorptivity (ε_eff by default) — recommended
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
When absorptivity=nothing (default), the material's mixture surface absorptivity
`a_eff` is used — identical to the effective emissivity (Kirchhoff's law) unless a
component carries its own `absorptivity`. Either way it automatically tracks
composition and temperature changes as material decomposes.

Do NOT combine `RadiativeFluxBC` with `HeatFluxBC` for the same radiation source -
this would double-count the energy input.
"""
struct RadiativeFluxBC{F,A} <: AbstractThermalBC
    flux::F           # Incident radiative flux [W/m²]
    absorptivity::A   # nothing = use effective, or constant, or function
end

# Constructors with optional absorptivity (flux positional or keyword)
RadiativeFluxBC(; flux, absorptivity = nothing, validate::Bool = true) =
    RadiativeFluxBC(flux; absorptivity = absorptivity, validate = validate)
function RadiativeFluxBC(flux; absorptivity = nothing, validate::Bool = true)
    if validate && absorptivity !== nothing
        α_val = _get_scalar_value(absorptivity)
        validate_bc_params(α = α_val)
    end
    RadiativeFluxBC(flux, absorptivity)
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
(e.g. `:top`, `:bottom` for 1D, additional tags for 2D/3D) maps to one
`AbstractThermalBC` and one `AbstractMassBC`. Internal storage is a typed
`NamedTuple` per BC kind, so symbol lookup compiles to a direct field read
and the residual specializes statically.

# Fields
- `thermal::TT`: `NamedTuple{tags, Tuple{<: AbstractThermalBC, ...}}`
- `mass::MT`: `NamedTuple{tags, Tuple{<: AbstractMassBC, ...}}`

# 1D convention
- `:top` ↔ exposed surface (z = L)
- `:bottom` ↔ substrate (z = 0)

# Construction
The keyword constructor is the user-facing API. Each BC kind takes a
NamedTuple of `tag => BC` mappings. A future 2D constructor takes a
richer tag set without changing the field shape.

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top    = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300),
	           bottom = AdiabaticBC()),
	mass    = (top    = AmbientGasBC(GasReservoir(P=0., T=300.)),
	           bottom = ImpermeableBC()),
)
```

# Accessors
- `get_thermal_bc(bc_set, Val(:top))` / `Val(:bottom)`
- `bc_set.mass.top` / `bc_set.mass.bottom`
"""
struct BoundaryConditionSet{TT<:NamedTuple,MT<:NamedTuple}
    thermal::TT
    mass::MT
end

# Default tag sets for 1D
const _BC_DEFAULT_THERMAL_1D = (top = AdiabaticBC(), bottom = AdiabaticBC())
const _BC_DEFAULT_MASS_1D = (top = ImpermeableBC(), bottom = ImpermeableBC())

"""
	BoundaryConditionSet(; thermal::NamedTuple = (top = AdiabaticBC(),  bottom = AdiabaticBC()),
						   mass::NamedTuple   = (top = ImpermeableBC(), bottom = ImpermeableBC()))

Tag-keyed constructor. Each BC kind (`thermal`, `mass`) takes a NamedTuple
of `tag => BC` mappings. The 1D convention uses `:top` and `:bottom`; future
2D/3D extensions add tags without changing the runtime contract.

# Example

```julia
bc_set = BoundaryConditionSet(
    thermal = (top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300),
               bottom = AdiabaticBC()),
    mass    = (top = AmbientGasBC(GasReservoir(P=0., T=300.)),
               bottom = ImpermeableBC()),
)
```
"""
function BoundaryConditionSet(;
    thermal::NamedTuple = _BC_DEFAULT_THERMAL_1D,
    mass::NamedTuple = _BC_DEFAULT_MASS_1D,
)
    # Merge user-supplied partial NamedTuples with the 1D defaults so callers
    # can specify only the slots they care about.
    thermal_full = merge(_BC_DEFAULT_THERMAL_1D, thermal)
    mass_full = merge(_BC_DEFAULT_MASS_1D, mass)
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

# ----- Tag-keyed accessors -----

"""
	get_thermal_bc(bc_set, ::Val{tag::Symbol})

Look up the thermal BC at a tag (`:top`, `:bottom`, …). Compile-time dispatch
via `Val{tag}` resolves to a direct NamedTuple field read.
"""
@inline get_thermal_bc(bc_set::BoundaryConditionSet, ::Val{tag}) where {tag} =
    bc_set.thermal[tag]


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
