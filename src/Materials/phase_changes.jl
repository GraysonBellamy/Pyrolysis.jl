"""
Phase changes of a liquid component and its vapour.
"""

# ============================================================================
# PhaseChange Type
# ============================================================================

"""
	PhaseChange{NS,A<:AbstractActivity,Ps<:AbstractPropertyFunction}

A two-way phase change between a liquid component and its vapour, held in
sorbing solids: per cell the liquid evaporates or desorbs at the rate

	ω_v = K (ξ_v^eq − ξ_v),   ξ_v^eq = φ M a(X, T) p_sat(T)/(R T),

toward the vapour concentration in equilibrium with the liquid, and takes the
heat `L_ref + h_v(T) − h_w(T) + Q_s(X, T)` per kilogram that becomes vapour
(negative `ω_v`, condensation or sorption, returns it). `X` is the moisture
content, the liquid's concentration over the weighted concentration of the
sorbing solids, `a` the activity law, `p_sat` the saturation curve, `φ` the
pore fraction, `M` the vapour's molar mass, `h_v` and `h_w` the components'
sensible enthalpies and `Q_s` the activity law's heat of sorption; see
[`phase_change_rate`](@ref) and [`phase_change_heat`](@ref).

[`Material`](@ref) builds it from the declaration of
[`PhaseChange`](@ref)`(name, liquid => vapour; ...)`, resolving the
components and copying in the vapour's molar mass.

# Fields
- `name::Symbol`: identifier
- `liquid::Int`: the liquid's component index
- `vapour::Int`: the vapour's component index
- `sorbents::NTuple{NS,Tuple{Int,Float64}}`: `(component_idx, weight)` per sorbing solid
- `activity::A`: the activity law, an [`AbstractActivity`](@ref)
- `saturation::Ps`: the liquid's saturation pressure [Pa] as a function of temperature
- `latent_heat::Float64`: the latent heat `L_ref` at `T_REF` [J/kg]
- `exchange_rate::Float64`: the exchange rate `K` [1/s]
- `molar_mass::Float64`: the vapour's molar mass [kg/mol]
"""
struct PhaseChange{NS,A<:AbstractActivity,Ps<:AbstractPropertyFunction}
    name::Symbol
    liquid::Int
    vapour::Int
    sorbents::NTuple{NS,Tuple{Int,Float64}}
    activity::A
    saturation::Ps
    latent_heat::Float64
    exchange_rate::Float64
    molar_mass::Float64
end

# A phase change as declared: its components as names or indices, and the
# sorbing solids with their weights. `Material` resolves it against its
# components.
struct PendingPhaseChange{L,V,S<:Tuple,A<:AbstractActivity,Ps<:AbstractPropertyFunction}
    name::Symbol
    liquid::L
    vapour::V
    sorbents::S
    activity::A
    saturation::Ps
    latent_heat::Float64
    exchange_rate::Float64
end

# ============================================================================
# Constructor
# ============================================================================

"""
	PhaseChange(name, liquid => vapour; activity, sorbents,
	            saturation = WATER_SATURATION, latent_heat = 2.4417e6,
	            exchange_rate = 1e3)

A phase change of the liquid component `liquid` into the gaseous component
`vapour`, each named by its component name or by its index in the material's
components, which [`Material`](@ref) resolves (`phase_changes = (...)`).

# Arguments
- `name::Symbol`: identifier
- `liquid => vapour`: the liquid and its vapour, names or indices
- `activity`: the activity law `a(X, T)`, an [`AbstractActivity`](@ref):
  [`CASTIsotherm`](@ref) (with [`CAST_SPRUCE`](@ref), [`CAST_BEECH`](@ref)),
  [`FreeLiquid`](@ref) or [`ActivityFunction`](@ref)
- `sorbents`: the sorbing solids, a tuple of component names or indices (weight
  one) or of `component => weight` pairs; the moisture content is the liquid's
  concentration over `Σₛ wₛ ξₛ`. Required: a char, say, does not hold water
  as the wood it comes from does, so no solid is a sorbent by default
- `saturation`: the liquid's saturation pressure [Pa] as a function of
  temperature, anything `to_property` takes (default [`WATER_SATURATION`](@ref))
- `latent_heat`: the latent heat at `T_REF` [J/kg], that of water at 25 °C by
  default
- `exchange_rate`: the exchange rate `K` [1/s], the inverse of the pore
  vapour's relaxation time

`Material` refuses, naming the phase change, a `liquid` that is not a
`LiquidComponent`, a `vapour` that is not a `GaseousComponent` with a positive
molar mass, a sorbent that is not a `SolidComponent` or is named twice, and
two phase changes of one liquid; this constructor refuses a weight that is
negative or not finite, sorbents none of which has a positive weight, an
`exchange_rate` that is not positive and finite and a `latent_heat` that is
not finite. A liquid may also be the reactant or a product of a
[`Reaction`](@ref).

# Example
```julia
PhaseChange(:moisture, :water => :vapour; activity = CAST_SPRUCE, sorbents = (:wood,))
PhaseChange(:moisture, 2 => 3; activity = FreeLiquid(), sorbents = (1 => 1.0, 4 => 0.5))
```
"""
function PhaseChange(
    name::Symbol,
    pair::Pair;
    activity::AbstractActivity,
    sorbents,
    saturation = WATER_SATURATION,
    latent_heat::Real = 2.4417e6,
    exchange_rate::Real = 1e3,
)
    _component_reference(pair.first, name)
    _component_reference(pair.second, name)
    entries = map(s -> _sorbent_entry(s, name), _sorbent_tuple(sorbents))
    isempty(entries) &&
        throw(ArgumentError("phase change :$name declares no sorbing solid"))
    any(e -> last(e) > 0, entries) ||
        throw(ArgumentError("phase change :$name: no sorbing solid has a positive weight"))
    (isfinite(exchange_rate) && exchange_rate > 0) || throw(
        ArgumentError(
            "phase change :$name has exchange_rate = $exchange_rate; the exchange " *
            "rate must be a finite positive number [1/s]",
        ),
    )
    isfinite(latent_heat) || throw(
        ArgumentError(
            "phase change :$name has latent_heat = $latent_heat; the latent heat " *
            "must be finite [J/kg]",
        ),
    )
    return PendingPhaseChange(
        name,
        pair.first,
        pair.second,
        entries,
        activity,
        to_property(saturation),
        Float64(latent_heat),
        Float64(exchange_rate),
    )
end

_component_reference(c::Union{Symbol,Integer}, name) = c
_component_reference(c, name) = throw(
    ArgumentError(
        "phase change :$name names the component $(repr(c)); give a component " *
        "name or index",
    ),
)

_sorbent_tuple(sorbents::Tuple) = sorbents
_sorbent_tuple(sorbent) = (sorbent,)

_sorbent_entry(s::Pair, name) =
    (_component_reference(s.first, name), _sorbent_weight(s.second, name))
_sorbent_entry(s, name) = (_component_reference(s, name), 1.0)

function _sorbent_weight(w, name)
    (w isa Real && isfinite(w) && w >= 0) || throw(
        ArgumentError(
            "phase change :$name has the sorbent weight $(repr(w)); a weight must " *
            "be finite and nonnegative",
        ),
    )
    return Float64(w)
end

# ============================================================================
# Resolution against the components
# ============================================================================

# The phase change with its components resolved against `components` (names
# `names`), checked, and the vapour's molar mass copied in.
function resolve_phase_change(pc::PendingPhaseChange, components::Tuple, names)
    liquid = _phase_component(pc.liquid, components, names, pc.name)
    vapour = _phase_component(pc.vapour, components, names, pc.name)
    is_liquid(components[liquid]) || throw(
        ArgumentError(
            "phase change :$(pc.name): component :$(names[liquid]) is not a " *
            "LiquidComponent",
        ),
    )
    M = molar_mass(components[vapour])
    (is_gas(components[vapour]) && isfinite(M) && M > 0) || throw(
        ArgumentError(
            "phase change :$(pc.name): component :$(names[vapour]) is not a " *
            "GaseousComponent with a positive molar mass",
        ),
    )
    sorbents = map(pc.sorbents) do (s, w)
        j = _phase_component(s, components, names, pc.name)
        is_solid(components[j]) || throw(
            ArgumentError(
                "phase change :$(pc.name): the sorbent :$(names[j]) is not a " *
                "SolidComponent",
            ),
        )
        (j, w)
    end
    allunique(map(first, sorbents)) ||
        throw(ArgumentError("phase change :$(pc.name) names a sorbing solid twice"))
    return PhaseChange(
        pc.name,
        liquid,
        vapour,
        sorbents,
        pc.activity,
        pc.saturation,
        pc.latent_heat,
        pc.exchange_rate,
        Float64(M),
    )
end

# A resolved phase change, checked again against `components`.
resolve_phase_change(pc::PhaseChange, components::Tuple, names) = resolve_phase_change(
    PendingPhaseChange(
        pc.name,
        pc.liquid,
        pc.vapour,
        pc.sorbents,
        pc.activity,
        pc.saturation,
        pc.latent_heat,
        pc.exchange_rate,
    ),
    components,
    names,
)

function _phase_component(c::Symbol, components, names, pc::Symbol)
    j = findfirst(==(c), names)
    j === nothing && throw(
        ArgumentError(
            "phase change :$pc references unknown component :$c; not found in the " *
            "material's components tuple",
        ),
    )
    return j
end

function _phase_component(j::Integer, components, names, pc::Symbol)
    1 <= j <= length(components) || throw(
        ArgumentError(
            "phase change :$pc references component $j, outside 1:$(length(components))",
        ),
    )
    return Int(j)
end

# The phase changes resolved against `components`, at most one per liquid.
function _resolve_phase_changes(phase_changes::Tuple, components::Tuple)
    isempty(phase_changes) && return ()
    all(pc -> pc isa Union{PendingPhaseChange,PhaseChange}, phase_changes) ||
        throw(ArgumentError("phase_changes must hold the declarations of PhaseChange"))
    names = map(c -> c.name, components)
    resolved = map(pc -> resolve_phase_change(pc, components, names), phase_changes)
    liquids = map(pc -> pc.liquid, resolved)
    allunique(liquids) || throw(
        ArgumentError(
            "the phase changes " *
            "$(join((repr(pc.name) for pc in resolved if count(==(pc.liquid), liquids) > 1), ", ")) " *
            "share a liquid; a liquid has at most one phase change",
        ),
    )
    return resolved
end
