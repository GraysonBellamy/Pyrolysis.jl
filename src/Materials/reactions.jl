"""
Reaction type definitions.

Reactions define the kinetics of thermal decomposition, including
stoichiometry, Arrhenius parameters, and heat of reaction.
"""

# ============================================================================
# Depletion Limiter
# ============================================================================

"""
	DepletionLimiter(; threshold = 1.0, join = 1.0, base_scale = 1e-8, enabled = true,
	                 conversion_threshold = 1e-3, conversion_base = 1e-8)

How a reaction rate vanishes as its reactant depletes; each `Material` carries
one. The reactant's rate factor is exactly zero at and below zero concentration
and C¹ through it, for every reaction order, so an absent
or negative trial concentration is never consumed and every rate has a finite
derivative at any real state.

A first-order reactant carries the factor on its concentration; a reactant of
any other order carries it on its conversion `x = ξ/ξ_ref`, with
`conversion_threshold` and `conversion_base` in place of `threshold` and
`base_scale`, both dimensionless and `join` shared. The conversion onset is
exactly one for `x ≥ 20 conversion_threshold` at `join ≤ 20`.

Per reactant the rate is multiplied by `g(ξ/threshold)`: the cubic onset
`(s − 2y) t³ + (3y − s) t²`, `t = x/join`, `y = tanh(join)`,
`s = join·sech²(join)`, which has zero value and slope at `x = 0` and matches
`tanh` and its slope at `x = join`, and `tanh(x)` above. The reaction-order
power is evaluated on the finite base `base_scale · softplus(ξ / base_scale)`,
which equals `ξ` to within `base_scale` and is finite for every real `ξ`. Far
above the threshold the factor rounds to one and the bare Arrhenius rate is
recovered (relative deficit `2e^(-2ξ/threshold)`: below 4 % at
`2·threshold`, below 5e-9 at `10·threshold`).

With `enabled = false` the factor is one and only the finite base remains (an
order-zero reaction then runs at its full rate at any concentration): an
absent reactant then reacts at the rate of a `base_scale` concentration, and
the solver's admissibility check is the guard.

Lower thresholds keep the kinetics closer to the bare rate law near depletion;
higher ones damp the stiffness of a vanishing reactant.

# Fields
- `threshold::Float64`: concentration scale of the roll-off [kg/m³]
- `join::Float64`: where the onset joins `tanh`, in units of `threshold`
- `base_scale::Float64`: scale of the finite power base [kg/m³]
- `enabled::Bool`: apply the factor
- `conversion_threshold::Float64`: conversion scale of the roll-off, orders other than one
- `conversion_base::Float64`: conversion scale of the finite power base, orders other than one

# Usage
```julia
Material(name = ..., components = ..., reactions = ...,
         depletion_limiter = DepletionLimiter(threshold = 0.1))
```
"""
struct DepletionLimiter
    threshold::Float64
    join::Float64
    base_scale::Float64
    enabled::Bool
    conversion_threshold::Float64
    conversion_base::Float64
    function DepletionLimiter(
        threshold::Real,
        join::Real,
        base_scale::Real,
        enabled::Bool,
        conversion_threshold::Real,
        conversion_base::Real,
    )
        threshold > 0 ||
            throw(ArgumentError("DepletionLimiter: threshold must be positive"))
        join > 0 || throw(ArgumentError("DepletionLimiter: join must be positive"))
        base_scale > 0 ||
            throw(ArgumentError("DepletionLimiter: base_scale must be positive"))
        conversion_threshold > 0 ||
            throw(ArgumentError("DepletionLimiter: conversion_threshold must be positive"))
        conversion_base > 0 ||
            throw(ArgumentError("DepletionLimiter: conversion_base must be positive"))
        return new(
            Float64(threshold),
            Float64(join),
            Float64(base_scale),
            enabled,
            Float64(conversion_threshold),
            Float64(conversion_base),
        )
    end
end

DepletionLimiter(;
    threshold::Real = 1.0,
    join::Real = 1.0,
    base_scale::Real = 1e-8,
    enabled::Bool = true,
    conversion_threshold::Real = 1e-3,
    conversion_base::Real = 1e-8,
) = DepletionLimiter(
    threshold,
    join,
    base_scale,
    enabled,
    conversion_threshold,
    conversion_base,
)

# ============================================================================
# Reaction Type
# ============================================================================

"""
	Reaction{NP,Ph<:AbstractPropertyFunction}

A single decomposition reaction with Arrhenius kinetics: one reactant and `NP`
products.

The reaction rate is computed as:
```
r = A * exp(-E/RT) * ξ_ref * x^n * H(T - T_min) * H(T_max - T),   x = ξ/ξ_ref
```
with `ξ` the reactant's concentration and `ξ_ref` the concentration of its
lineage inventory, the condensed mass its lineage started with at the material
point (first order: `ξ_ref` cancels and `r = k ξ`). `A` is in s⁻¹ for every
order. Per unit reaction extent the reactant is consumed with mass
stoichiometry **1.0** and each product is produced with its yield, so mass
conservation requires `Σ yields == 1`.

Reaction orders must be `n ≥ 0` (a negative order makes the rate singular at
reactant depletion). At runtime the reactant also contributes the material's
depletion factor to the rate, and the hard temperature gates `H(·)` are realized
as smooth tanh ramps of width `gate_width` [K]
(`0.5·(1 + tanh((T − T_min)/gate_width))` and mirrored) — see
`Pyrolysis.Physics.reaction_rate` and [`DepletionLimiter`](@ref). Smaller
widths approach the hard step at the cost of Jacobian stiffness near the gate.

# Type Parameters
- `NP`: Number of products
- `Ph`: Heat-of-reaction property type

# Fields
- `name::Symbol`: Reaction identifier
- `reactant::Int`: The reactant's component index
- `order::Float64`: The reaction order
- `products::NTuple{NP, Tuple{Int, Float64}}`: `(component_idx, yield)` per product
- `A::Float64`: Pre-exponential factor [s⁻¹]
- `E::Float64`: Activation energy [J/mol]
- `heat::Ph`: Heat of reaction [J/kg of reactant]
- `T_min::Float64`: Minimum temperature for reaction [K]
- `T_max::Float64`: Maximum temperature for reaction [K]
- `gate_width::Float64`: Width of the smooth tanh temperature gates [K]
- `basis::Float64`: The reference as a fraction of the material point's initial
  condensed mass, on the cells the reactant's lineage can reach, in place of
  the lineage rules; `NaN` for the lineage rules. Only for orders other than one.

# Heat of Reaction Convention
- `h > 0`: Endothermic (absorbs heat, cools material)
- `h < 0`: Exothermic (releases heat, heats material)
"""
struct Reaction{NP,Ph<:AbstractPropertyFunction}
    name::Symbol
    reactant::Int
    order::Float64
    products::NTuple{NP,Tuple{Int,Float64}}
    A::Float64
    E::Float64
    heat::Ph
    T_min::Float64
    T_max::Float64
    gate_width::Float64
    basis::Float64
end

# ============================================================================
# Validation Helpers
# ============================================================================

"""
	_validate_reaction_order(n::Real, name::Symbol)

Reject negative reaction orders at construction. In the rate a negative order
is singular as the reactant depletes (`x^n → ∞` as `x → 0`) — no depletion
limiter can repair that. Orders `n ≥ 0` (including 0 and fractional TGA fits) are supported; their
near-depletion behavior is regularized by the depletion limiter
([`DepletionLimiter`](@ref)).
"""
function _validate_reaction_order(n::Real, name::Symbol)
    n >= 0 || throw(
        ArgumentError(
            "Reaction :$name has reaction order n = $n; negative orders are " *
            "unphysical (the rate is singular as the reactant depletes). " *
            "Supported orders are n ≥ 0.",
        ),
    )
    return nothing
end

"""
	_validate_gate_width(w::Real, name::Symbol)

Reject non-positive or non-finite temperature-gate widths at construction.
The gates are evaluated as `0.5·(1 + tanh((T − T_min)/gate_width))` (mirrored
for `T_max`): `w ≤ 0` inverts or degenerates the ramp, and `w = Inf` freezes
both gates at 0.5 everywhere.
"""
function _validate_gate_width(w::Real, name::Symbol)
    (isfinite(w) && w > 0) || throw(
        ArgumentError(
            "Reaction :$name has gate_width = $w; the tanh temperature-gate " *
            "width must be a finite positive number of kelvin.",
        ),
    )
    return nothing
end

"""
	_validate_basis(basis, first_order::Bool, name::Symbol)

Reject a `basis` outside `(0, 1]` or on a first-order reaction, whose reference
cancels. `nothing` asks for the lineage rules and passes.
"""
function _validate_basis(basis, first_order::Bool, name::Symbol)
    basis === nothing && return nothing
    (isfinite(basis) && 0 < basis <= 1) || throw(
        ArgumentError(
            "Reaction :$name has basis = $basis; a basis is a fraction of the " *
            "initial condensed mass in (0, 1].",
        ),
    )
    first_order && throw(
        ArgumentError(
            "Reaction :$name is first order, whose reference cancels: give no basis.",
        ),
    )
    return nothing
end

"""
	_validate_yield(y::Real, name::Symbol)

Reject negative product mass stoichiometries (yields) at construction. A
negative yield consumes a product species in proportion to the reaction rate —
the rate law has no knowledge of that species' concentration, so its ξ is
driven negative without bound. Mass-balance checks alone cannot catch this:
yields (1.5, −0.5) sum to 1.
"""
function _validate_yield(y::Real, name::Symbol)
    y >= 0 || throw(
        ArgumentError(
            "Reaction :$name has a negative product yield ($y); negative " *
            "yields consume a species the rate law does not see and drive " *
            "its concentration negative without bound. Yields must be ≥ 0.",
        ),
    )
    return nothing
end

"""
	_validate_stoichiometry(products, name::Symbol)

Validate that the yields conserve mass: per unit reaction extent the reactant
contributes mass 1.0, so `Σ yields == 1` to a relative tolerance of 1e-10.
"""
function _validate_stoichiometry(products::NTuple{NP,Tuple{Int,Float64}}, name) where {NP}
    product_mass = sum(last, products; init = 0.0)
    isapprox(1.0, product_mass; atol = max(1e-10 * (1.0 + product_mass), 1e-14)) || throw(
        ArgumentError(
            "Reaction :$name violates the mass balance: the yields sum to " *
            "$product_mass, not 1",
        ),
    )
    return nothing
end

# ============================================================================
# Constructors
# ============================================================================

"""
	Reaction(name, reactant => products; A, E, h, yields = (1.0,), n = 1.0,
	         T_min = 0.0, T_max = Inf, gate_width = 1.0, basis = nothing)

A reaction of `reactant` into one product or a tuple of them, each named by its
index in the material's components or, on both sides of the pair, by its
component name, which [`Material`](@ref) resolves.

# Arguments
- `name::Symbol`: Reaction identifier
- `reactant => products`: Component indices or names
- `A`: Pre-exponential factor [s⁻¹]
- `E`: Activation energy [J/mol]
- `h`: Heat of reaction [J/kg] (number or property function), positive endothermic
- `yields`: Mass yield of each product, summing to one
- `n=1.0`: Reaction order
- `T_min=0.0`: Minimum temperature [K]
- `T_max=Inf`: Maximum temperature [K]
- `gate_width=1.0`: Width of the smooth tanh temperature gates [K] (> 0)
- `basis=nothing`: The reference as a fraction of the initial condensed mass (`n ≠ 1`)

# Example
```julia
Reaction(:decomposition, :virgin => :MMA, A = 8.5e12, E = 188e3, h = 870e3)
Reaction(:charring, :virgin => (:char, :gas), A = 1e15, E = 200e3, h = 500e3,
         yields = (0.2, 0.8))
```
"""
function Reaction(
    name::Symbol,
    pair::Pair{Int};
    A::Real,
    E::Real,
    h,
    yields::Tuple{Vararg{Real}} = (1.0,),
    n::Real = 1.0,
    T_min::Real = 0.0,
    T_max::Real = Inf,
    gate_width::Real = 1.0,
    basis::Union{Nothing,Real} = nothing,
)
    products = _products(pair.second)
    length(products) == length(yields) || throw(
        ArgumentError(
            "Reaction :$name has $(length(products)) products and $(length(yields)) yields",
        ),
    )
    _validate_reaction_order(n, name)
    foreach(y -> _validate_yield(y, name), yields)
    stoichiometry = map((p, y) -> (Int(p), Float64(y)), products, yields)
    _validate_stoichiometry(stoichiometry, name)
    _validate_gate_width(gate_width, name)
    _validate_basis(basis, n == 1, name)
    return Reaction(
        name,
        pair.first,
        Float64(n),
        stoichiometry,
        Float64(A),
        Float64(E),
        to_property(h),
        Float64(T_min),
        Float64(T_max),
        Float64(gate_width),
        basis === nothing ? NaN : Float64(basis),
    )
end

_products(products::Tuple) = products
_products(product) = (product,)

# A reaction whose reactant and products are still component names, as
# `Reaction` received it; `Material` resolves it against its components.
struct PendingReaction{P<:Tuple,K<:NamedTuple}
    name::Symbol
    reactant::Symbol
    products::P
    options::K
end

function Reaction(name::Symbol, pair::Pair{Symbol}; options...)
    products = _products(pair.second)
    # The options are checked here, before the components are known.
    Reaction(name, 1 => ntuple(identity, length(products)); options...)
    return PendingReaction(name, pair.first, products, NamedTuple(options))
end

# The reaction with its component names replaced by their places in `names`.
resolve_reaction(rxn::Reaction, names) = rxn
function resolve_reaction(rxn::PendingReaction, names)
    index(component) = _component_place(names, component, rxn.name)
    return Reaction(
        rxn.name,
        index(rxn.reactant) => map(index, rxn.products);
        rxn.options...,
    )
end

function _component_place(names, component::Symbol, reaction::Symbol)
    j = findfirst(==(component), names)
    j === nothing && throw(
        ArgumentError(
            "Reaction :$reaction references unknown component :$component; " *
            "not found in the material's components tuple",
        ),
    )
    return j
end
