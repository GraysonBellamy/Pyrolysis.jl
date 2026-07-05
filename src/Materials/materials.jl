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
	Material{NC,NR,C,RT,LSL} <: AbstractMaterial

A complete material definition with components, reactions, mixing rules, and gas transport.

# Type Parameters
- `NC`: Number of components
- `NR`: Number of reactions
- `C`: Component tuple type (preserves heterogeneous component types)
- `RT`: Reaction tuple type (preserves heterogeneous reaction types)
- `LSL`: Lateral shrinkage law type (concrete `Function` subtype, or `Nothing`).
  Parameterizing here lets calls through `material.lateral_shrinkage_law`
  specialize statically — no dynamic dispatch through `Union{Function,Nothing}`.

# Fields
- `name::Symbol`: Material identifier
- `components::C`: Material components tuple
- `reactions::RT`: Decomposition reactions tuple
- `mixing::NamedTuple{(:k,:λ,:κ),Tuple{MixingSpec,MixingSpec,MixingSpec}}`:
  Mixing rules (and per-rule β) for thermal conductivity (`k`), gas transfer
  coefficient (`λ`), and permeability (`κ`).
- `lateral_shrinkage_law::LSL`: Optional A(t)/A_0 = law(χ̄). `nothing` ⇒ identity.

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
struct Material{NC,NR,C<:Tuple,RT<:Tuple,LSL} <: AbstractMaterial
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
    # material (smooth tanh(ξ/threshold) roll-off; see DepletionLimiter).
    depletion_limiter::DepletionLimiter

    # Inner constructor with validation
    function Material{NC,NR,C,RT,LSL}(
        name::Symbol,
        components::C,
        reactions::RT,
        mixing::NamedTuple{(:k, :λ, :κ),Tuple{MixingSpec,MixingSpec,MixingSpec}},
        lateral_shrinkage_law::LSL,
        depletion_limiter::DepletionLimiter,
    ) where {NC,NR,C<:Tuple,RT<:Tuple,LSL}
        # Validate component count matches
        length(components) == NC || throw(ArgumentError("Component count mismatch"))
        length(reactions) == NR || throw(ArgumentError("Reaction count mismatch"))
        LSL === Nothing ||
            LSL <: Function ||
            throw(ArgumentError("lateral_shrinkage_law must be a Function or nothing"))
        new{NC,NR,C,RT,LSL}(
            name,
            components,
            reactions,
            mixing,
            lateral_shrinkage_law,
            depletion_limiter,
        )
    end
end

# ============================================================================
# Constructors
# ============================================================================

"""
	Material(; name, components, reactions=(),
			  mixing = MaterialMixing(),
			  lateral_shrinkage_law = nothing,
			  depletion_limiter = DepletionLimiter())

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
- `depletion_limiter`: [`DepletionLimiter`](@ref) applied to every
  reaction-rate evaluation for this material (smooth `tanh(ξ/threshold)`
  roll-off near reactant depletion). Default: `threshold = 1.0` kg/m³, enabled.

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
    depletion_limiter::DepletionLimiter = DepletionLimiter(),
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
    )
end

function _construct_material(
    name::Symbol,
    components::C,
    resolved::RT,
    mixing,
    lateral_shrinkage_law::LSL,
    depletion_limiter::DepletionLimiter,
) where {C<:Tuple,RT<:Tuple,LSL}
    Material{length(components),length(resolved),C,RT,LSL}(
        name,
        components,
        resolved,
        mixing,
        lateral_shrinkage_law,
        depletion_limiter,
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
	n_solids(m::Material)

Get the number of solid-phase components.
"""
function n_solids(m::Material)
    count(is_solid, m.components)
end

"""
	n_liquids(m::Material)

Get the number of liquid-phase components.
"""
function n_liquids(m::Material)
    count(is_liquid, m.components)
end

"""
	n_gases(m::Material)

Get the number of gas-phase components.
"""
function n_gases(m::Material)
    count(is_gas, m.components)
end

"""
	has_darcy_flow(m::Material) -> Bool

Check if material uses Darcy-Fick gas transport.

Returns true if any solid/liquid component has non-zero permeability (κ > 0).

The component-based approach: define κ on each solid component,
and the effective permeability will be computed via mixing rules.
"""
function has_darcy_flow(m::Material)
    # Check component-based permeability
    for comp in m.components
        if has_permeability(comp)
            # Check if κ > 0 at a reference temperature
            if get_permeability(comp, 500.0) > 0.0
                return true
            end
        end
    end
    return false
end

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
# Validation
# ============================================================================

"""
	validate_material(m::Material)

Validate material definition for consistency.

Checks:
- All reaction indices refer to valid components
- Stoichiometry is mass-balanced
- At least one component exists
"""
function validate_material(m::Material{NC,NR}) where {NC,NR}
    # Must have at least one component
    NC >= 1 || error("Material must have at least one component")

    # Validate reactions
    for (i, rxn) in enumerate(m.reactions)
        for (idx, _order) in rxn.reactants
            1 <= idx <= NC ||
                error("Reaction $i ($(rxn.name)): reactant index $idx out of range 1:$NC")
        end
        for (idx, _ν) in rxn.products
            1 <= idx <= NC ||
                error("Reaction $i ($(rxn.name)): product index $idx out of range 1:$NC")
        end

        # Check mass balance
        verify_stoichiometry(rxn) || @warn "Reaction $i ($(rxn.name)) fails mass balance"
    end

    return true
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
end

# ============================================================================
# Generated Functions for Type-Stable Iteration over Heterogeneous Components
# ============================================================================
#
# These functions use @generated to unroll loops at compile time, avoiding
# the type instability that occurs when indexing heterogeneous tuples at runtime.
# This is critical for performance in hot paths like RHS evaluation.

"""
	is_gas_phase(material::Material, j::Int) -> Bool

Check if component j is in gas phase. Uses dispatch for type stability.
"""
@generated function is_gas_phase(material::Material{NC}, ::Val{J}) where {NC,J}
    if J < 1 || J > NC
        return :(error("Component index $J out of range 1:$NC"))
    end
    return :(is_gas(material.components[$J]))
end

# Convenience method that takes Int (not type-stable, but used in non-critical paths)
is_gas_phase(material::Material, j::Int) = is_gas(material.components[j])

"""
	is_solid_phase(material::Material, j::Int) -> Bool

Check if component j is in solid phase. Uses dispatch for type stability.
"""
@generated function is_solid_phase(material::Material{NC}, ::Val{J}) where {NC,J}
    if J < 1 || J > NC
        return :(error("Component index $J out of range 1:$NC"))
    end
    return :(is_solid(material.components[$J]))
end

# Convenience method
is_solid_phase(material::Material, j::Int) = is_solid(material.components[j])

"""
	get_swelling_factor(material::Material, ::Val{J}) -> Float64

Get swelling factor for component J. Type-stable via dispatch.
"""
@generated function get_swelling_factor(material::Material{NC}, ::Val{J}) where {NC,J}
    if J < 1 || J > NC
        return :(error("Component index $J out of range 1:$NC"))
    end
    return :(material.components[$J].swelling_factor)
end

"""
	get_component_density(material::Material, ::Val{J}, T) -> Float64

Get density of component J at temperature T. Type-stable via dispatch.
"""
@generated function get_component_density(material::Material{NC}, ::Val{J}, T) where {NC,J}
    if J < 1 || J > NC
        return :(error("Component index $J out of range 1:$NC"))
    end
    return :(material.components[$J].density(T))
end

"""
	foreach_component(f, material::Material{NC}) where NC

Apply function f to each component index, using dispatch for type stability.
f should take (material, Val{j}) as arguments.

Example:
	foreach_component(material) do mat, j_val
		if is_gas_phase(mat, j_val)
			# do something
		end
	end
"""
@generated function foreach_component(f::F, material::Material{NC}) where {F,NC}
    exprs = [:(f(material, Val{$j}())) for j = 1:NC]
    return Expr(:block, exprs...)
end

"""
	foreach_gas_component(f, material::Material{NC}) where NC

Apply function f only to gas component indices. Type-stable.
f should take (material, component_index::Int) as arguments.
"""
@generated function foreach_gas_component(f::F, material::Material{NC}) where {F,NC}
    # Build list of gas components at compile time
    exprs = Expr[]
    for j = 1:NC
        push!(exprs, :(
            if is_gas(material.components[$j])
                ;f(material, $j);
            end
        ))
    end
    return Expr(:block, exprs...)
end

"""
	gas_component_indices(::Material{NC,NR,C,RT}) where {NC,NR,C,RT}

Return a tuple of gas component indices at compile time. Use this instead of
building an array at runtime.
"""
@generated function gas_component_indices(material::Material{NC}) where {NC}
    # This generates the expression to build a tuple of gas indices
    # We need to evaluate at generated function time, which is tricky
    # Instead, return code that builds the right tuple
    return quote
        _collect_gas_indices(material.components)
    end
end

# Helper that recursively collects gas indices - unrolled at compile time
# Phase test is now a type-level check via `is_gas`, so this resolves at compile time.
@inline _collect_gas_indices(::Tuple{}) = ()
@inline function _collect_gas_indices(comps::Tuple)
    first_comp = comps[1]
    rest = Base.tail(comps)
    rest_indices = _collect_gas_indices(rest)
    if is_gas(first_comp)
        return (1, map(x -> x + 1, rest_indices)...)
    else
        return map(x -> x + 1, rest_indices)
    end
end

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

# ============================================================================
# Unrolled Component Loops for Hot Paths
# ============================================================================

"""
	_unroll_gas_flux!(J::Matrix, j::Val{J}, ξ_L, ξ_R, f, flux_fn, material) where J

Unrolled gas flux computation for a single component.
If component J is GAS, compute flux and store in J[J, f].
Uses the volume fraction gradient formulation.
"""
@generated function _unroll_gas_flux!(
    J::Matrix,
    ::Val{IDX},
    ξ_L::NTuple{NC},
    ξ_R::NTuple{NC},
    f::Int,
    λ_face,
    T_L,
    T_R,
    d_LR,
    material::Material{NC},
) where {IDX,NC}
    if IDX < 1 || IDX > NC
        return :(error("Component index $IDX out of range 1:$NC"))
    end
    return quote
        @inbounds if is_gas(material.components[$IDX])
            ξ_L_j = ξ_L[$IDX]
            ξ_R_j = ξ_R[$IDX]
            # Volume fraction gradient: J = -(λ/T_face) × Δ(ξT) / d
            T_face = 0.5 * (T_L + T_R)
            J[$IDX, f] = -λ_face * (ξ_R_j * T_R - ξ_L_j * T_L) / (T_face * d_LR)
        end
    end
end

"""
	compute_gas_fluxes_unrolled!(J, f, ξ_L, ξ_R, λ_face, T_L, T_R, d_LR, material::Material{NC})

Compute gas fluxes for all components at face f, using compile-time unrolling.

Uses the volume fraction gradient formulation where λ_face is the harmonic mean of λ:
	J = -(λ_face / T_face) × (ξ_R × T_R - ξ_L × T_L) / d_LR
"""
@generated function compute_gas_fluxes_unrolled!(
    J::Matrix,
    f::Int,
    ξ_L::NTuple{NC},
    ξ_R::NTuple{NC},
    λ_face,
    T_L,
    T_R,
    d_LR,
    material::Material{NC},
) where {NC}

    exprs = Expr[]
    for j = 1:NC
        push!(
            exprs,
            quote
                @inbounds if is_gas(material.components[$j])
                    ξ_L_j = ξ_L[$j]
                    ξ_R_j = ξ_R[$j]
                    # Volume fraction gradient: J = -(λ/T_face) × Δ(ξT) / d
                    T_face = 0.5 * (T_L + T_R)
                    J[$j, f] = -λ_face * (ξ_R_j * T_R - ξ_L_j * T_L) / (T_face * d_LR)
                end
            end,
        )
    end
    return Expr(:block, exprs...)
end


"""
	compute_gas_fluxes_darcy_fick_unrolled!(J::Matrix, f::Int, v_g::Float64,
		ξ_L::NTuple{NC}, ξ_R::NTuple{NC}, φ_L, φ_R, λ_face, T_L, T_R, d_LR,
		material::Material{NC})

Compile-time unrolled Darcy-Fick flux computation for type stability.

Computes J = J_adv + J_diff for each gas component:
- J_adv = v_g × ξ_upwind / φ_face (first-order upwind; ξ is per bulk volume,
  v_g is superficial, so the advected scalar is the intrinsic density ξ/φ with
  φ_face = (φ_L + φ_R)/2 floored at 1e-12 like the pressure closure)
- J_diff = -(λ_face/T_face) × (ξ_R×T_R - ξ_L×T_L) / d_LR (volume fraction gradient)
"""
@generated function compute_gas_fluxes_darcy_fick_unrolled!(
    J::Matrix,
    f::Int,
    v_g,
    ξ_L::NTuple{NC},
    ξ_R::NTuple{NC},
    φ_L,
    φ_R,
    λ_face,
    T_L,
    T_R,
    d_LR,
    material::Material{NC},
) where {NC}

    exprs = Expr[:(inv_φ_face = 1.0 / max(0.5 * (φ_L + φ_R), 1e-12))]
    for j = 1:NC
        push!(
            exprs,
            quote
                @inbounds if is_gas(material.components[$j])
                    ξ_L_j = ξ_L[$j]
                    ξ_R_j = ξ_R[$j]

                    # Diffusive flux (volume fraction gradient)
                    T_face = 0.5 * (T_L + T_R)
                    J_diff = -λ_face * (ξ_R_j * T_R - ξ_L_j * T_L) / (T_face * d_LR)

                    # Advective flux: intrinsic density upwinded (first-order)
                    ξ_upwind = v_g >= 0.0 ? ξ_L_j : ξ_R_j
                    J_adv = v_g * ξ_upwind * inv_φ_face

                    # Combined Darcy-Fick flux
                    J[$j, f] = J_adv + J_diff
                end
            end,
        )
    end
    return Expr(:block, exprs...)
end

"""
	set_boundary_advective_flux_unrolled!(J::Matrix, f::Int, v_over_φ,
		ξ_int::NTuple{NC}, material::Material{NC})

Write the boundary-face advective venting flux `J[j, f] = v_over_φ · ξ_int[j]`
for every gas component (solids untouched). `v_over_φ = v_b/φ_int` is the
interstitial boundary Darcy velocity; the caller has already applied the
outflow-only upwind gate (the ambient ghost carries no tracked gas).
"""
@generated function set_boundary_advective_flux_unrolled!(
    J::Matrix,
    f::Int,
    v_over_φ,
    ξ_int::NTuple{NC},
    material::Material{NC},
) where {NC}
    exprs = Expr[]
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                J[$j, f] = v_over_φ * ξ_int[$j]
            end
        end)
    end
    return Expr(:block, exprs...)
end


"""
	apply_species_residual_unrolled!(dξ::Matrix, i::Int, div_J::Matrix, dξ_rxn::Matrix, material::Material{NC})

Apply species residual: for gases dξ = -div_J + dξ_rxn, for solids dξ = dξ_rxn.
Uses compile-time unrolling for type stability.
"""
@generated function apply_species_residual_unrolled!(
    dξ::Matrix,
    i::Int,
    div_J::Matrix,
    dξ_rxn::Matrix,
    material::Material{NC},
) where {NC}

    exprs = Expr[]
    for j = 1:NC
        push!(exprs, quote
            @inbounds if is_gas(material.components[$j])
                dξ[$j, i] = -div_J[$j, i] + dξ_rxn[$j, i]
            else
                dξ[$j, i] = dξ_rxn[$j, i]
            end
        end)
    end
    return Expr(:block, exprs...)
end
