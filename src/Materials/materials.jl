"""
Material type definitions.

A Material combines components and reactions into a complete
description of a pyrolyzing substance.
"""

# ============================================================================
# Mixing Rule Enumeration
# ============================================================================

"""
	MixingRule

Enumeration of property mixing rules for multi-component materials.

- `PARALLEL`: Parallel (volume-weighted) mixing: k = Σᵢ vᵢkᵢ
- `SERIES`: Series (harmonic) mixing: 1/k = Σᵢ vᵢ/kᵢ
- `WEIGHTED`: Weighted combination: k = β·k_parallel + (1-β)·k_series
- `BRUGGEMAN`: Bruggeman effective medium: Σᵢ vᵢ(kᵢ - k_eff)/(kᵢ + 2k_eff) = 0.
  Defined only for scalar transport coefficients (thermal conductivity `k`);
  rejected for the gas-transfer coefficient `λ` and for permeability `κ`.
- `CARMAN_KOZENY`: Porosity-based permeability closure (κ channel only):
  κ = ℓ²·φ³ / (180·(1-φ)²), with ℓ the characteristic particle/pore diameter
  supplied as `MixingSpec(CARMAN_KOZENY; length_scale = ℓ)`. This is the
  EMT-replacement for permeability — `BRUGGEMAN` is not a valid κ closure.
"""
@enum MixingRule begin
    PARALLEL = 1
    SERIES = 2
    WEIGHTED = 3
    BRUGGEMAN = 4
    CARMAN_KOZENY = 5
end

# ============================================================================
# MixingSpec
# ============================================================================

"""
	MixingSpec(rule::MixingRule, β::Float64 = 0.5; length_scale::Real = NaN)

A mixing-rule + weighting pair. `β` is only meaningful when `rule == WEIGHTED`,
but is always carried so a single uniform shape can describe every channel.

`length_scale` (ℓ) is only meaningful when `rule == CARMAN_KOZENY` (permeability
channel), where it is the characteristic particle/pore diameter [m] entering
κ = ℓ²·φ³/(180·(1-φ)²). It is required (a positive, finite value) for that rule
and ignored otherwise; the constructor fails fast if `CARMAN_KOZENY` is requested
without it.

Used by `Material.mixing` to bundle the mixing rule and its parameters for each
of the three mixed properties: thermal conductivity (`k`), gas transfer
coefficient (`λ`), and permeability (`κ`).
"""
struct MixingSpec
    rule::MixingRule
    β::Float64
    ℓ::Float64   # Carman–Kozeny characteristic length [m]; NaN unless rule == CARMAN_KOZENY

    function MixingSpec(rule::MixingRule, β::Real = 0.5; length_scale::Real = NaN)
        0 ≤ β ≤ 1 || throw(ArgumentError("β must be in [0,1], got $β"))
        if rule == CARMAN_KOZENY
            (isfinite(length_scale) && length_scale > 0) || throw(
                ArgumentError(
                    "CARMAN_KOZENY requires a positive `length_scale` (characteristic " *
                    "particle/pore diameter in metres), got $length_scale",
                ),
            )
        end
        new(rule, Float64(β), Float64(length_scale))
    end
end

"""
	MaterialMixing(; k = MixingSpec(WEIGHTED, 0.5),
				    λ = MixingSpec(PARALLEL, 0.5),
				    κ = MixingSpec(WEIGHTED, 0.5))

Build the `mixing` NamedTuple expected by `Material`. `k`, `λ`, `κ` may be passed
either as `MixingSpec` instances or as bare `MixingRule` enum values (β defaults
to 0.5 in that case).
"""
function MaterialMixing(;
    k::Union{MixingSpec,MixingRule} = MixingSpec(WEIGHTED, 0.5),
    λ::Union{MixingSpec,MixingRule} = MixingSpec(PARALLEL, 0.5),
    κ::Union{MixingSpec,MixingRule} = MixingSpec(WEIGHTED, 0.5),
)
    k_spec = k isa MixingSpec ? k : MixingSpec(k, 0.5)
    λ_spec = λ isa MixingSpec ? λ : MixingSpec(λ, 0.5)
    κ_spec = κ isa MixingSpec ? κ : MixingSpec(κ, 0.5)
    return (k = k_spec, λ = λ_spec, κ = κ_spec)
end

# ============================================================================
# Material Type
# ============================================================================

"""
	Material{NC,NR,C,RT,LSL,Pμ,DL} <: AbstractMaterial

A complete material definition with components, reactions, mixing rules, and gas transport.

# Type Parameters
- `NC`: Number of components
- `NR`: Number of reactions
- `C`: Component tuple type (preserves heterogeneous component types)
- `RT`: Reaction tuple type (preserves heterogeneous reaction types)
- `LSL`: Lateral shrinkage law type (concrete `Function` subtype, or `Nothing`).
  Parameterizing here lets calls through `material.lateral_shrinkage_law`
  specialize statically — no dynamic dispatch through `Union{Function,Nothing}`.
- `Pμ`: Pore-gas viscosity property-function type. Parameterized for the same
  reason as `LSL` — `material.gas_viscosity(T)` stays monomorphic.

# Fields
- `name::Symbol`: Material identifier
- `components::C`: Material components tuple
- `reactions::RT`: Decomposition reactions tuple
- `mixing::NamedTuple{(:k,:λ,:κ),Tuple{MixingSpec,MixingSpec,MixingSpec}}`:
  Mixing rules (and per-rule β) for thermal conductivity (`k`), gas transfer
  coefficient (`λ`), and permeability (`κ`).
- `lateral_shrinkage_law::LSL`: Optional A(t)/A_0 = law(χ̄). `nothing` ⇒ identity.
- `gas_viscosity::Pμ`: Dynamic viscosity μ(T) of the pore gas [Pa·s]. Accepts
  anything `to_property` understands — a number, a `(p0, p1)` tuple, a
  `Function`, or a property function. Defaults to
  [`SutherlandProperty`](@ref)`()` (air).

# Pore-gas viscosity
μ enters the model *only* as the Darcy mobility κ/μ, so an error in μ is
equivalent to the reciprocal error in permeability. The default is air, which
is within ~2 % for an N₂-swept pore space but overstates the viscosity of
pyrolysis volatiles by roughly 1.3–1.7× over 300–1000 K — a bias in the same
direction everywhere, and one that grows as the pores fill with volatiles.
Supply a blend correlation for the gas the material actually makes:

```julia
# Cantera mixture-averaged fit for wood pyrolyzates, ascending coefficients
# (`PolynomialProperty` uses `evalpoly`), valid 280–1200 K.
wood = Material(
	name = :WhitePine,
	components = (...),
	reactions = (...),
	gas_viscosity = PolynomialProperty(mu_coeffs),
)
```

Because this is one curve per material, it is an average over the composition
the material is expected to produce; it cannot resolve the pore gas shifting
from background air to volatiles as pyrolysis proceeds.

# Gas Transport
Darcy-Fick gas transport is enabled when any solid/liquid component has κ > 0.
The effective permeability is computed from component values via mixing rules,
similar to thermal conductivity and gas transfer coefficient.

Specify κ on each component:
```julia
charring = Material(
	name = :CharringWood,
	components = (
		SolidComponent(:virgin, ρ=530, c=2000, k=0.12, κ=1e-15),  # tight structure
		SolidComponent(:char, ρ=150, c=1500, k=0.20, κ=1e-12),    # porous char
		GaseousComponent(:gas, M=0.044, c=1800, k=0.03, λ=1e-4),
	),
	reactions = (...),
	mixing = MaterialMixing(κ = MixingSpec(WEIGHTED, 0.5)),
)
```

# Example
```julia
# Simple PMMA (diffusion-only gas transport, κ=0 on all components)
pmma = Material(
	name = :PMMA,
	components = (
		SolidComponent(:virgin, ρ=1190, c=1420, k=0.21, ε=0.86),
		GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
	),
	reactions = (
		Reaction(:decomposition, 1 => 2, A=8.5e12, E=188e3, h=870e3),  # h > 0 = endothermic
	),
	mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5)),
)
```
"""
struct Material{NC,NR,C<:Tuple,RT<:Tuple,LSL,Pμ,DL<:AbstractDepletionLimiter} <:
       AbstractMaterial
    name::Symbol
    components::C
    reactions::RT  # Full tuple type to preserve heterogeneous reaction types

    # Bundled mixing specs for the three mixed channels (k, λ, κ). NamedTuple
    # field order is fixed so accessors compile down to direct field reads.
    mixing::NamedTuple{(:k, :λ, :κ),Tuple{MixingSpec,MixingSpec,MixingSpec}}

    # Variable cross-section: optional law A(t)/A_0 = law(χ̄). `nothing` ⇒ identity.
    # LSL is parameterized on the concrete law type so calls specialize (no
    # dynamic dispatch through Union{Function, Nothing}).
    lateral_shrinkage_law::LSL

    # Depletion limiter applied to every reaction-rate evaluation for this
    # material (see DepletionLimiter). DL is the concrete
    # limiter type, so rate evaluation specializes on it.
    depletion_limiter::DL

    # Dynamic viscosity of the pore gas, μ(T) [Pa·s]. Enters the model only as
    # the Darcy mobility κ/μ, so it is a single material-level curve rather
    # than a per-component property: it describes the gas mixture this material
    # actually produces. Defaults to `SutherlandProperty()` (air).
    gas_viscosity::Pμ

    # Inner constructor with validation
    function Material{NC,NR,C,RT,LSL,Pμ,DL}(
        name::Symbol,
        components::C,
        reactions::RT,
        mixing::NamedTuple{(:k, :λ, :κ),Tuple{MixingSpec,MixingSpec,MixingSpec}},
        lateral_shrinkage_law::LSL,
        depletion_limiter::DL,
        gas_viscosity::Pμ,
    ) where {NC,NR,C<:Tuple,RT<:Tuple,LSL,Pμ,DL<:AbstractDepletionLimiter}
        # Validate component count matches
        length(components) == NC || throw(ArgumentError("Component count mismatch"))
        length(reactions) == NR || throw(ArgumentError("Reaction count mismatch"))
        LSL === Nothing ||
            LSL <: Function ||
            throw(ArgumentError("lateral_shrinkage_law must be a Function or nothing"))
        _validate_gas_viscosity(gas_viscosity, name)
        new{NC,NR,C,RT,LSL,Pμ,DL}(
            name,
            components,
            reactions,
            mixing,
            lateral_shrinkage_law,
            depletion_limiter,
            gas_viscosity,
        )
    end
end

# μ enters only as the denominator of the Darcy mobility κ/μ, so a curve that is
# non-positive or non-finite anywhere in the modelled range does not degrade the
# answer — it reverses or blows up the gas velocity. Correlations fitted over a
# finite window (the Cantera surrogate's poly4 is fitted 280–1200 K) can do
# exactly that on extrapolation. Fail at construction, in the spirit of
# `MixingSpec`'s CARMAN_KOZENY/`length_scale` check, rather than at hour three of
# a solve. The probe spans the pyrolysis range plus generous margin; the runtime
# floor `Properties.MU_FLOOR`, which the state mapping applies, remains the
# backstop for Newton trial states that wander outside it.
function _validate_gas_viscosity(gas_viscosity, name::Symbol)
    for T in (200.0, 300.0, 500.0, 800.0, 1200.0, 1800.0)
        μ = gas_viscosity(T)
        (isfinite(μ) && μ > 0) || throw(
            ArgumentError(
                "Material :$name has a gas_viscosity that is not positive and finite " *
                "at T = $T K (got $μ). μ appears only as the Darcy mobility κ/μ, so a " *
                "non-positive value reverses the gas velocity rather than merely " *
                "degrading it. Check the correlation's validity range — a polynomial " *
                "fitted over a narrow window will do this on extrapolation.",
            ),
        )
    end
    return nothing
end

# ============================================================================
# Constructors
# ============================================================================

"""
	Material(; name, components, reactions=(),
			  mixing = MaterialMixing(),
			  lateral_shrinkage_law = nothing,
			  depletion_limiter = DepletionLimiter(),
			  gas_viscosity = SutherlandProperty())

Create a material with keyword arguments.

# Arguments
- `name::Symbol`: Material identifier
- `components`: Tuple of components
- `reactions`: Tuple of reactions (can be empty tuple `()`)
- `mixing`: A `NamedTuple{(:k,:λ,:κ),Tuple{MixingSpec,MixingSpec,MixingSpec}}`,
  typically built via [`MaterialMixing`](@ref). Defaults: k = `WEIGHTED, β=0.5`,
  λ = `PARALLEL, β=0.5` (gas transfer; `BRUGGEMAN` is rejected at evaluation
  time as physically unjustified for this channel), κ = `WEIGHTED, β=0.5`.
- `lateral_shrinkage_law`: optional A(t)/A_0 = law(χ̄). `nothing` ⇒ identity.
- `depletion_limiter`: the [`AbstractDepletionLimiter`](@ref) applied to every
  reaction-rate evaluation for this material. Default: a
  [`DepletionLimiter`](@ref) (a roll-off near reactant depletion that is zero
  at and below zero concentration) with `threshold = 1.0` kg/m³, enabled.
- `gas_viscosity`: dynamic viscosity μ(T) of the pore gas [Pa·s], a number or
  a property function of temperature, used only in the Darcy mobility κ/μ.
  Default: [`SutherlandProperty`](@ref) (air). It must be positive and finite
  from 200 K to 1800 K.

# Example
```julia
mat = Material(
	name = :PMMA,
	components = (virgin, gas),
	reactions = (rxn,),
	mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.6),
							λ = MixingSpec(PARALLEL),
							κ = MixingSpec(WEIGHTED, 0.5)),
)
```
"""
function Material(;
    name::Symbol,
    components::C,
    reactions::Tuple = (),
    mixing::NamedTuple{(:k, :λ, :κ),Tuple{MixingSpec,MixingSpec,MixingSpec}} = MaterialMixing(),
    lateral_shrinkage_law::Union{Function,Nothing} = nothing,
    depletion_limiter::AbstractDepletionLimiter = DepletionLimiter(),
    gas_viscosity = SutherlandProperty(),
) where {C<:Tuple}
    # Resolve any PendingReactions (symbolic component names → integer indices)
    # against the components tuple's name → index map.
    resolved = _resolve_reactions(reactions, components)
    # Function barrier: constructing through a helper whose `LSL` is bound by the
    # argument type makes the `Material{…,LSL}` apply_type monomorphic in each
    # branch, instead of a `Union{Function,Nothing}` split JET flags as a possible
    # no-matching-method on the inner constructor.
    return _construct_material(
        name,
        components,
        resolved,
        mixing,
        lateral_shrinkage_law,
        depletion_limiter,
        to_property(gas_viscosity),
    )
end

function _construct_material(
    name::Symbol,
    components::C,
    resolved::RT,
    mixing,
    lateral_shrinkage_law::LSL,
    depletion_limiter::DL,
    gas_viscosity::Pμ,
) where {C<:Tuple,RT<:Tuple,LSL,Pμ,DL<:AbstractDepletionLimiter}
    Material{length(components),length(resolved),C,RT,LSL,Pμ,DL}(
        name,
        components,
        resolved,
        mixing,
        lateral_shrinkage_law,
        depletion_limiter,
        gas_viscosity,
    )
end

# Walk the reactions tuple, resolve every PendingReaction to a Reaction. The
# `name → index` map is built once per Material(...) call.
function _resolve_reactions(reactions::Tuple, components::Tuple)
    # Detect whether resolution is needed before building the map (avoids the
    # Dict allocation in the common no-symbolic-reactions case).
    needs_resolution = false
    for rxn in reactions
        if rxn isa PendingReaction
            needs_resolution = true
            break
        end
    end
    needs_resolution || return reactions
    name_to_index = Dict{Symbol,Int}()
    for (i, c) in enumerate(components)
        haskey(name_to_index, c.name) &&
            error("Duplicate component name :$(c.name) in components tuple")
        name_to_index[c.name] = i
    end
    return map(rxn -> resolve_reaction(rxn, name_to_index), reactions)
end

# ============================================================================
# Accessor Functions
# ============================================================================

"""
	n_components(m::Material{NC}) where NC

Get the number of components in a material.
"""
n_components(::Material{NC}) where {NC} = NC

"""
	n_reactions(m::Material{NC,NR}) where {NC,NR}

Get the number of reactions in a material.
"""
n_reactions(::Material{NC,NR}) where {NC,NR} = NR

"""
	has_darcy_flow(m::Material) -> Bool

Check if material uses Darcy-Fick gas transport.

Returns true if the permeability law of any solid/liquid component is not
[`identically_zero`](@ref).

The component-based approach: define κ on each solid component,
and the effective permeability will be computed via mixing rules.
"""
has_darcy_flow(m::Material) =
    any(c -> has_permeability(c) && !identically_zero(c.permeability), m.components)

"""
	component_index(m::Material, name::Symbol)

Get the index of a component by name.
"""
function component_index(m::Material, name::Symbol)
    for (i, c) in enumerate(m.components)
        c.name == name && return i
    end
    error("Component :$name not found in material $(m.name)")
end

# ============================================================================
# Display
# ============================================================================

function Base.show(io::IO, m::Material{NC,NR}) where {NC,NR}
    print(io, "Material(:$(m.name), $NC components, $NR reactions)")
end

function Base.show(io::IO, ::MIME"text/plain", m::Material{NC,NR}) where {NC,NR}
    println(io, "Material: $(m.name)")
    println(io, "  Components ($NC):")
    for (i, c) in enumerate(m.components)
        println(io, "    $i. $(c.name) ($(phase(c)))")
    end
    if NR > 0
        println(io, "  Reactions ($NR):")
        for (i, r) in enumerate(m.reactions)
            println(io, "    $i. $(r.name): A=$(r.A), E=$(r.E) J/mol")
        end
    end
    println(io, "  Conductivity mixing: $(m.mixing.k.rule), β=$(m.mixing.k.β)")
    println(io, "  Gas transfer mixing: $(m.mixing.λ.rule), β=$(m.mixing.λ.β)")
    println(io, "  Permeability mixing: $(m.mixing.κ.rule), β=$(m.mixing.κ.β)")
    if !(m.gas_viscosity isa SutherlandProperty && m.gas_viscosity == SutherlandProperty())
        println(
            io,
            "  Pore-gas viscosity: $(typeof(m.gas_viscosity).name.name) " *
            "(μ(300 K) = $(round(m.gas_viscosity(300.0) * 1e6, digits = 2)) μPa·s, " *
            "μ(700 K) = $(round(m.gas_viscosity(700.0) * 1e6, digits = 2)) μPa·s)",
        )
    end
end

# ============================================================================
# Generated Functions for Type-Stable Iteration over Heterogeneous Components
# ============================================================================
#
# These functions use @generated to unroll loops at compile time, avoiding
# the type instability that occurs when indexing heterogeneous tuples at runtime.
# This is critical for performance in hot paths like RHS evaluation.

"""
	dry_pyrolysis_gas_indices(material::Material{NC}) where NC

Return a tuple of indices of gas components representing dry-pyrolysis gases —
i.e. all GAS-phase components except water vapor / moisture-derived gases.

Used by the variable cross-section model: only dry-pyrolysis-gas mass efflux
drives the radial shrinkage progress variable χ_i. Moisture loss does not
contribute to χ.

Recognized water-vapor names (case-insensitive substring match): `:water_vapor`,
`:WATER_VAPOR`, `:H2O`, anything containing `H2O` or `water`.
"""
@generated function dry_pyrolysis_gas_indices(material::Material{NC}) where {NC}
    return quote
        _collect_dry_pyrolysis_gas_indices(material.components)
    end
end

# Compile-time recursive collector for dry-pyrolysis-gas indices.
# Phase check is type-level via `is_gas`; the water-vapor name filter remains runtime.
@inline _collect_dry_pyrolysis_gas_indices(::Tuple{}) = ()
@inline function _collect_dry_pyrolysis_gas_indices(comps::Tuple)
    first_comp = comps[1]
    rest = Base.tail(comps)
    rest_indices = _collect_dry_pyrolysis_gas_indices(rest)
    if is_gas(first_comp) && !_is_water_vapor_name(first_comp.name)
        return (1, map(x -> x + 1, rest_indices)...)
    else
        return map(x -> x + 1, rest_indices)
    end
end

@inline function _is_water_vapor_name(name::Symbol)
    s = String(name)
    if occursin("H2O", s)
        return true
    end
    ls = lowercase(s)
    return occursin("water", ls) || ls == "h2o"
end
