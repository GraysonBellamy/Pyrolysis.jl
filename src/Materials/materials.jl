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
  permeability's effective-medium closure; `BRUGGEMAN` is not a valid one.
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

Bundled per mixed property by [`MaterialMixing`](@ref).
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
	MaterialMixing(; k = MixingSpec(WEIGHTED), λ = MixingSpec(PARALLEL),
	               κ = MixingSpec(WEIGHTED))

The mixing rules of a `Material`, one per mixed property: thermal conductivity
`k`, gas transfer coefficient `λ` and permeability `κ`. Each is a `MixingSpec`
or a bare `MixingRule`.
"""
struct MaterialMixing
    k::MixingSpec
    λ::MixingSpec
    κ::MixingSpec
end

MaterialMixing(;
    k::Union{MixingSpec,MixingRule} = MixingSpec(WEIGHTED),
    λ::Union{MixingSpec,MixingRule} = MixingSpec(PARALLEL),
    κ::Union{MixingSpec,MixingRule} = MixingSpec(WEIGHTED),
) = MaterialMixing(MixingSpec(k), MixingSpec(λ), MixingSpec(κ))

MixingSpec(spec::MixingSpec) = spec

# ============================================================================
# Material Type
# ============================================================================

"""
	Material{NC,NR,C,RT,LSL,Pμ,GT}

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
- `mixing::MaterialMixing`: Mixing rules (and per-rule β) for thermal
  conductivity (`k`), gas transfer coefficient (`λ`), and permeability (`κ`).
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
`gas_transport=LumpedTransfer()` selects the superficial partial-pressure law.
`DustyGas` selects coupled molecular/Knudsen friction. Both include independent
Darcy permeability and use the same physical pore EOS. Participating lumped
matrices require explicit condensed λ. Gas-bearing initial states require
resolved pore volume. See the model constructors for endpoint closures.

Specify κ on each component:
```julia
charring = Material(
	name = :CharringWood,
	components = (
		SolidComponent(:virgin, ρ=530, φ=.6, c=2000, k=0.12, λ=1e-5, κ=1e-15),  # tight structure
		SolidComponent(:char, ρ=150, φ=.85, c=1500, k=0.20, λ=1e-5, κ=1e-12),    # porous char
		GaseousComponent(:gas, M=0.044, c=1800, k=0.03,),
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
		GaseousComponent(:MMA, M=0.100, c=1500, k=0.02,)
	),
	reactions = (
		Reaction(:decomposition, 1 => 2, A=8.5e12, E=188e3, h=870e3),  # h > 0 = endothermic
	),
	mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5)),
)
```
"""
struct Material{NC,NR,C<:Tuple,RT<:Tuple,LSL,Pμ,GT<:GasTransport}
    name::Symbol
    components::C
    reactions::RT  # Full tuple type to preserve heterogeneous reaction types

    # Bundled mixing specs for the three mixed channels (k, λ, κ).
    mixing::MaterialMixing

    # Variable cross-section: optional law A(t)/A_0 = law(χ̄). `nothing` ⇒ identity.
    # LSL is parameterized on the concrete law type so calls specialize (no
    # dynamic dispatch through Union{Function, Nothing}).
    lateral_shrinkage_law::LSL

    # Depletion limiter applied to every reaction-rate evaluation for this
    # material (see DepletionLimiter).
    depletion_limiter::DepletionLimiter

    # Dynamic viscosity of the pore gas, μ(T) [Pa·s]. Enters the model only as
    # the Darcy mobility κ/μ, so it is a single material-level curve rather
    # than a per-component property: it describes the gas mixture this material
    # actually produces. Defaults to `SutherlandProperty()` (air).
    gas_viscosity::Pμ
    gas_transport::GT

    # Inner constructor with validation
    function Material{NC,NR,C,RT,LSL,Pμ,GT}(
        name::Symbol,
        components::C,
        reactions::RT,
        mixing::MaterialMixing,
        lateral_shrinkage_law::LSL,
        depletion_limiter::DepletionLimiter,
        gas_viscosity::Pμ,
        gas_transport::GT,
    ) where {NC,NR,C<:Tuple,RT<:Tuple,LSL,Pμ,GT<:GasTransport}
        # Validate component count matches
        length(components) == NC || throw(ArgumentError("Component count mismatch"))
        length(reactions) == NR || throw(ArgumentError("Reaction count mismatch"))
        LSL === Nothing ||
            LSL <: Function ||
            throw(ArgumentError("lateral_shrinkage_law must be a Function or nothing"))
        _validate_gas_viscosity(gas_viscosity, name)
        _validate_transport_inputs(components, gas_transport, mixing)
        new{NC,NR,C,RT,LSL,Pμ,GT}(
            name,
            components,
            reactions,
            mixing,
            lateral_shrinkage_law,
            depletion_limiter,
            gas_viscosity,
            gas_transport,
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
			  gas_viscosity = SutherlandProperty(),
			  gas_transport = LumpedTransfer())

Create a material with keyword arguments.

# Arguments
- `name::Symbol`: Material identifier
- `components`: Tuple of components
- `reactions`: Tuple of reactions (can be empty tuple `()`)
- `mixing`: A [`MaterialMixing`](@ref). Defaults: k = `WEIGHTED, β=0.5`,
  λ = `PARALLEL, β=0.5` (superficial matrix transfer), κ = `WEIGHTED, β=0.5`.
  `BRUGGEMAN` is rejected for both transport coefficients at construction.
- `lateral_shrinkage_law`: optional A(t)/A_0 = law(χ̄). `nothing` ⇒ identity.
- `depletion_limiter`: the [`DepletionLimiter`](@ref) applied to every
  reaction-rate evaluation for this material: a roll-off near reactant depletion
  that is zero at and below zero concentration, by default at
  `threshold = 1.0` kg/m³.
- `gas_transport`: [`LumpedTransfer`](@ref) or [`DustyGas`](@ref), including
  the explicit endpoint and its domain when needed. Gas-free materials bind
  the internal no-gas specialization.
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
    mixing::MaterialMixing = MaterialMixing(),
    lateral_shrinkage_law::Union{Function,Nothing} = nothing,
    depletion_limiter::DepletionLimiter = DepletionLimiter(),
    gas_viscosity = SutherlandProperty(),
    gas_transport::GasTransport = LumpedTransfer(),
) where {C<:Tuple}
    resolved = _resolve_reactions(reactions, components)
    # Function barrier: constructing through a helper whose `LSL` is bound by the
    # argument type makes the `Material{…,LSL}` apply_type monomorphic in each
    # branch, instead of a `Union{Function,Nothing}` split JET flags as a possible
    # no-matching-method on the inner constructor.
    material = _construct_material(
        name,
        components,
        resolved,
        mixing,
        lateral_shrinkage_law,
        depletion_limiter,
        to_property(gas_viscosity),
        gas_transport,
    )
    needs_lineage(material) && lineage_graph(material)
    return material
end

function _construct_material(
    name::Symbol,
    components::C,
    resolved::RT,
    mixing,
    lateral_shrinkage_law::LSL,
    depletion_limiter::DepletionLimiter,
    gas_viscosity::Pμ,
    gas_transport::GT,
) where {C<:Tuple,RT<:Tuple,LSL,Pμ,GT<:GasTransport}
    Material{length(components),length(resolved),C,RT,LSL,Pμ,GT}(
        name,
        components,
        resolved,
        mixing,
        lateral_shrinkage_law,
        depletion_limiter,
        gas_viscosity,
        gas_transport,
    )
end

# The reactions with their component names resolved against `components`.
function _resolve_reactions(reactions::Tuple, components::Tuple)
    any(r -> r isa PendingReaction, reactions) || return reactions
    names = map(c -> c.name, components)
    allunique(names) ||
        throw(ArgumentError("the component names of a material must differ"))
    return map(r -> resolve_reaction(r, names), reactions)
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
# Dry pyrolysis gases
# ============================================================================

"""
	dry_gas_flags(material::Material{NC}) -> NTuple{NC,Bool}

Flag the dry pyrolysis gases: the gases the solids give off through the
reaction graph. The lineage starts with every solid component and takes in,
until nothing changes, the products (yield > 0) of every reaction with a
reactant in it; its gases are the dry gases. Water released by decomposition
is dry gas, the vapour a liquid evaporates into is not, and a gas no reaction
produces is not.

A lineage gas that a reaction without a reactant in the lineage also produces
has no single origin, and the function throws an `ArgumentError` naming the
gas and the reactions; model it as two species, one per lineage.
"""
function dry_gas_flags(material::Material{NC}) where {NC}
    comps = material.components
    rxns = material.reactions
    lineage = Set(j for j = 1:NC if is_solid(comps[j]))
    from_lineage(r) = r.reactant in lineage
    produces(r, j) = any(p -> first(p) == j && last(p) > 0, r.products)
    grown = true
    while grown
        grown = false
        for r in rxns
            from_lineage(r) || continue
            for (j, y) in r.products
                if y > 0 && !(j in lineage)
                    push!(lineage, j)
                    grown = true
                end
            end
        end
    end
    for j = 1:NC
        (is_gas(comps[j]) && j in lineage) || continue
        outside = [r.name for r in rxns if !from_lineage(r) && produces(r, j)]
        isempty(outside) && continue
        inside = [r.name for r in rxns if from_lineage(r) && produces(r, j)]
        throw(
            ArgumentError(
                "gas :$(comps[j].name) is produced from the solids (reactions " *
                "$(join(repr.(inside), ", "))) and from outside their lineage " *
                "(reactions $(join(repr.(outside), ", "))), so whether it is dry " *
                "pyrolysis gas is ambiguous; split it into one species per lineage",
            ),
        )
    end
    return ntuple(j -> is_gas(comps[j]) && j in lineage, Val(NC))
end

# ============================================================================
# Reaction lineages
# ============================================================================

"""
	needs_lineage(material) -> Bool

Whether some reaction of `material` has an order other than one, which is
evaluated on the conversion of its reactant against the reactant's lineage
inventory. Every first-order reaction is evaluated on its concentration, and a
material with no other order carries no reference.
"""
needs_lineage(material::Material) = any(r -> r.order != 1, material.reactions)

"""
	LineageGraph

What the reaction network alone fixes of the reference of each reaction, per
reaction in the material's order: `reactant`, the reactant's component index;
`order`, its order; `basis`, the reaction's `basis` (`NaN` for the lineage
rules); `ancestors`, per component whether the reactant can be formed from it,
the reactant itself included (empty for a first-order reaction); `reach`, per
component the mass of the reactant one unit of its inventory can become (empty
for a first-order reaction and under a `basis`). `condensed` flags the
condensed components.
"""
struct LineageGraph
    names::Vector{Symbol}
    components::Vector{Symbol}
    condensed::Vector{Bool}
    reactant::Vector{Int}
    order::Vector{Float64}
    basis::Vector{Float64}
    ancestors::Vector{Vector{Bool}}
    reach::Vector{Vector{Float64}}
end

# The channels of every condensed component: per reaction it consumes, the
# condensed products other than itself at positive yield with their weights,
# the yields renormalised by what the reaction consumes. A reaction that
# returns its whole reactant consumes nothing and is no channel.
function _lineage_channels(material::Material{NC}, condensed) where {NC}
    channels = [Vector{Pair{Int,Float64}}[] for _ = 1:NC]
    for r in material.reactions
        a = r.reactant
        condensed[a] || continue
        consumed = 1.0 - sum(y for (p, y) in r.products if p == a; init = 0.0)
        consumed > 0 || continue
        push!(
            channels[a],
            [p => y / consumed for (p, y) in r.products if p != a && y > 0 && condensed[p]],
        )
    end
    return channels
end

# The nodes reachable from `seeds` along `next`, the seeds included.
function _reachable(next::Vector{Vector{Int}}, seeds)
    seen = falses(length(next))
    front = Int[seeds...]
    seen[front] .= true
    while !isempty(front)
        v = pop!(front)
        for p in next[v]
            seen[p] && continue
            seen[p] = true
            push!(front, p)
        end
    end
    return Vector{Bool}(seen)
end

# The nodes of `member` children first, and the cycle among them along
# `children`, `nothing` when there is none.
function _children_first(children::Vector{Vector{Int}}, member::Vector{Bool})
    state = zeros(Int8, length(children))
    order, path = Int[], Int[]
    function visit(v)
        state[v] = 1
        push!(path, v)
        for p in children[v]
            member[p] || continue
            state[p] == 1 && return path[findfirst(==(p), path):end]
            if state[p] == 0
                cycle = visit(p)
                cycle === nothing || return cycle
            end
        end
        pop!(path)
        state[v] = 2
        push!(order, v)
        return nothing
    end
    for v in eachindex(member)
        (member[v] && state[v] == 0) || continue
        cycle = visit(v)
        cycle === nothing || return order, cycle
    end
    return order, nothing
end

# The reach of the reactant `a` over its ancestors, visited children first:
# one unit of a node's inventory becomes the common value of its channels
# that reach `a`, a channel producing nothing toward `a` being ignored.
function _reach(channels, ancestors, order, a::Int, rxn::Symbol, names)
    Φ = zeros(length(ancestors))
    Φ[a] = 1.0
    for v in order
        v == a && continue
        value = NaN
        for channel in channels[v]
            any(e -> ancestors[first(e)], channel) || continue
            φ = sum(y * Φ[p] for (p, y) in channel; init = 0.0)
            if isnan(value)
                value = φ
            elseif abs(φ - value) > 1e-12 * max(abs(φ), abs(value))
                throw(
                    ArgumentError(
                        "reaction :$rxn: competing channels of :$(names[v]) reach " *
                        ":$(names[a]) with different yields ($value and $φ), so the " *
                        "reference depends on the rates: give `basis`",
                    ),
                )
            end
        end
        Φ[v] = value
    end
    return Φ
end

"""
	lineage_graph(material) -> LineageGraph

The graph rules of the reference of every reaction of order other than one, in
order: its reactant is condensed and not fed from outside the condensed
inventory (a condensed product of a gas-consuming reaction is fed, and so is
every node reachable from it along positive-yield edges), whatever its
`basis`; without a `basis`, the reactant's ancestors hold no cycle, and the
channels of each ancestor that reach it agree on the reach, the mass of the
reactant one unit of the ancestor's inventory can become. A reaction that
returns its whole reactant is no edge, and the yields of one that returns part
of it are renormalised by what it consumes. Throws an `ArgumentError` naming
the reaction otherwise, with `basis` as the remedy where there is one.
"""
function lineage_graph(material::Material{NC}) where {NC}
    for r in material.reactions
        all(j -> 1 <= j <= NC, (r.reactant, first.(r.products)...)) ||
            throw(ArgumentError("reaction :$(r.name) references a component outside 1:$NC"))
    end
    names = [c.name for c in material.components]
    condensed = [!is_gas(c) for c in material.components]
    channels = _lineage_channels(material, condensed)
    children = [unique(first(e) for channel in channels[v] for e in channel) for v = 1:NC]
    parents = [[v for v = 1:NC if p in children[v]] for p = 1:NC]
    gas_products = Int[
        p for r in material.reactions if !condensed[r.reactant] for
        (p, y) in r.products if y > 0 && condensed[p]
    ]
    fed = _reachable(children, gas_products)

    reactions = collect(material.reactions)
    ancestors = [Bool[] for _ in reactions]
    reach = [Float64[] for _ in reactions]
    for (k, r) in enumerate(reactions)
        a, n = r.reactant, r.order
        n == 1 && continue
        condensed[a] || throw(
            ArgumentError(
                "reaction :$(r.name) has the gas reactant :$(names[a]) at order $n: " *
                "an order other than one needs a condensed reactant",
            ),
        )
        fed[a] && throw(
            ArgumentError(
                "reaction :$(r.name): :$(names[a]) is fed by a gas-consuming reaction, " *
                "so its lineage has no initial inventory; no `basis` repairs this",
            ),
        )
        ancestors[k] = _reachable(parents, a)
        r.basis === nothing || continue
        order, cycle = _children_first(children, ancestors[k])
        cycle === nothing || throw(
            ArgumentError(
                "reaction :$(r.name): the ancestors of :$(names[a]) hold the cycle " *
                "$(join((names[v] for v in (cycle..., first(cycle))), " → ")): give `basis`",
            ),
        )
        reach[k] = _reach(channels, ancestors[k], order, a, r.name, names)
    end
    return LineageGraph(
        [r.name for r in reactions],
        names,
        condensed,
        [r.reactant for r in reactions],
        [r.order for r in reactions],
        [something(r.basis, NaN) for r in reactions],
        ancestors,
        reach,
    )
end
