"""
Reaction type definitions.

Reactions define the kinetics of thermal decomposition, including
stoichiometry, Arrhenius parameters, and heat of reaction.
"""

# ============================================================================
# Depletion Limiter
# ============================================================================

"""
	AbstractDepletionLimiter

How a reaction rate vanishes as a reactant depletes. Each `Material` carries
one; [`DepletionLimiter`](@ref) and [`OneSidedLimiter`](@ref) are the two
kinds.
"""
abstract type AbstractDepletionLimiter end

"""
Depletion limiter settings for reaction rate computation.

The depletion limiter removes the numerical stiffness of near-zero reactant
concentrations by multiplying the reaction rate with a smooth `tanh(ξ/threshold)`
factor per reactant.

# Fields
- `threshold::Float64`: Concentration scale of the rate roll-off [kg/m³]
- `enabled::Bool`: If false, no limiting is applied (may cause numerical issues)

# Behavior
The factor `tanh(ξ/threshold)` is applied to the *rate* (not the power-law
base) and is C^∞ in ξ — there is no threshold branch, so neither the rate nor
its derivative jumps anywhere:

- `ξ ≫ threshold`: factor → 1 double-exponentially (relative deficit
  `2e^(-2ξ/threshold)`: < 4% at `2·threshold`, < 5e-9 at `10·threshold`) — the
  bare Arrhenius rate is recovered.
- `ξ → 0`: the rate vanishes ∝ `ξ^(n+1)` with bounded `∂r/∂ξ` for every
  reaction order n ≥ 0 — including n = 0, where the unlimited rate would stay
  at k(T) through ξ = 0 and drive the concentration negative without bound.

# Usage
Each `Material` carries its own limiter; pass a custom one at construction:

```julia
Material(name = ..., components = ..., reactions = ...,
         depletion_limiter = DepletionLimiter(threshold = 0.1))
```

The solver's kinetics path (`Pyrolysis.Physics.evaluate_reactions!`) reads
`material.depletion_limiter` for every rate evaluation.

# Notes
Lower thresholds are more physically accurate but may cause numerical stiffness.
For validation studies, temporarily disable with `enabled=false`. Disabling
removes the only depletion guard: order n = 0 reactions then consume reactant
at full rate through ξ = 0 (unbounded negative ξ), and orders 0 < n < 1 have a
singular ∂r/∂ξ as ξ → 0.
"""
Base.@kwdef struct DepletionLimiter <: AbstractDepletionLimiter
    threshold::Float64 = 1.0    # kg/m³ - rate roll-off scale
    enabled::Bool = true        # Can disable for testing/validation
end

"""
	OneSidedLimiter(; threshold = 1.0, join = 1.0, base_scale = 1e-8, enabled = true)

Depletion limiter that is exactly zero at and below zero concentration and C¹
through it, for every reaction order.

Per reactant the rate is multiplied by `g(ξ/threshold)`: the cubic onset
`(s − 2y) t³ + (3y − s) t²`, `t = x/join`, `y = tanh(join)`,
`s = join·sech²(join)`, which has zero value and slope at `x = 0` and matches
`tanh` and its slope at `x = join`, and `tanh(x)` above. The reaction-order
power is evaluated on the finite base `base_scale · softplus(ξ / base_scale)`,
which equals `ξ` to within `base_scale` and is finite for every real `ξ`. Above
`join · threshold` the rate equals the [`DepletionLimiter`](@ref) rate with the
same `threshold` to round-off, so calibrated kinetics are unchanged outside the
onset interval.

With `enabled = false` the factor is one and only the finite base remains: an
absent reactant then reacts at the rate of a `base_scale` concentration, and
the solver's admissibility check is the guard.

# Fields
- `threshold::Float64`: concentration scale of the roll-off [kg/m³]
- `join::Float64`: where the onset joins `tanh`, in units of `threshold`
- `base_scale::Float64`: scale of the finite power base [kg/m³]
- `enabled::Bool`: apply the factor
"""
struct OneSidedLimiter <: AbstractDepletionLimiter
    threshold::Float64
    join::Float64
    base_scale::Float64
    enabled::Bool
    function OneSidedLimiter(threshold::Real, join::Real, base_scale::Real, enabled::Bool)
        threshold > 0 || throw(ArgumentError("OneSidedLimiter: threshold must be positive"))
        join > 0 || throw(ArgumentError("OneSidedLimiter: join must be positive"))
        base_scale > 0 ||
            throw(ArgumentError("OneSidedLimiter: base_scale must be positive"))
        return new(Float64(threshold), Float64(join), Float64(base_scale), enabled)
    end
end

OneSidedLimiter(;
    threshold::Real = 1.0,
    join::Real = 1.0,
    base_scale::Real = 1e-8,
    enabled::Bool = true,
) = OneSidedLimiter(threshold, join, base_scale, enabled)

"""
	OneSidedLimiter(l::AbstractDepletionLimiter) -> OneSidedLimiter

The one-sided limiter with the threshold and the enabled flag of `l`.
"""
OneSidedLimiter(l::OneSidedLimiter) = l
OneSidedLimiter(l::DepletionLimiter) =
    OneSidedLimiter(; threshold = l.threshold, enabled = l.enabled)

# ============================================================================
# Reaction Type
# ============================================================================

"""
	Reaction{NRin, NRout, Ph<:AbstractPropertyFunction} <: AbstractReaction

A single decomposition reaction with Arrhenius kinetics, n-ary in both
reactants and products.

The reaction rate is computed as:
```
r = A * exp(-E/RT) * ∏_i ξ_{r_i}^{n_i} * H(T - T_min) * H(T_max - T)
```
where the product runs over `reactants`. Per unit reaction extent, the first
reactant is consumed with mass stoichiometry **1.0**; products are produced
with their per-product `stoich` from the `products` tuple.

Reaction orders must be `n ≥ 0` (enforced at construction; a negative order
makes the rate singular at reactant depletion). At runtime each reactant also
contributes a smooth `tanh(ξ/ξ_th)` depletion factor to the rate, and the hard
temperature gates `H(·)` are realized as smooth tanh ramps of width
`gate_width` [K] (`0.5·(1 + tanh((T − T_min)/gate_width))` and mirrored) — see
`Pyrolysis.Physics.reaction_rate` and [`DepletionLimiter`](@ref). Smaller
widths approach the hard step at the cost of Jacobian stiffness near the gate.

# Type Parameters
- `NRin`:  Number of reactants
- `NRout`: Number of products
- `Ph`:    Heat-of-reaction property type

# Fields
- `name::Symbol`: Reaction identifier
- `reactants::NTuple{NRin, Tuple{Int, Float64}}`: `(component_idx, reaction_order)` per reactant.
- `products::NTuple{NRout, Tuple{Int, Float64}}`: `(component_idx, mass_stoich)` per product.
- `A::Float64`: Pre-exponential factor [s⁻¹] or [m³/(kg·s)]
- `E::Float64`: Activation energy [J/mol]
- `heat::Ph`: Heat of reaction [J/kg of reactant₁]
- `T_min::Float64`: Minimum temperature for reaction [K]
- `T_max::Float64`: Maximum temperature for reaction [K]
- `gate_width::Float64`: Width of the smooth tanh temperature gates [K] (> 0, default 1.0)

# Mass-balance convention

For multi-reactant reactions (NRin > 1), additional reactant mass stoichiometries
would be required; the current shape supports the single-reactant case used
throughout the existing models (`NRin == 1`, mass consumed = 1·r). Multi-reactant
support can extend the per-reactant tuple to 3-ary `(idx, order, mass_stoich)`
without changing this struct's outer shape.

For single-reactant reactions, mass conservation requires `Σ stoich(products) == 1.0`.

# Heat of Reaction Convention
- `h > 0`: Endothermic (absorbs heat, cools material)
- `h < 0`: Exothermic (releases heat, heats material)
"""
struct Reaction{NRin,NRout,Ph<:AbstractPropertyFunction} <: AbstractReaction
    name::Symbol
    reactants::NTuple{NRin,Tuple{Int,Float64}}    # (component_idx, reaction_order)
    products::NTuple{NRout,Tuple{Int,Float64}}    # (component_idx, mass_stoich)
    A::Float64
    E::Float64
    heat::Ph
    T_min::Float64
    T_max::Float64
    gate_width::Float64                           # tanh temperature-gate width [K]

    function Reaction(
        name::Symbol,
        reactants::NTuple{NRin,Tuple{Int,Float64}},
        products::NTuple{NRout,Tuple{Int,Float64}},
        A::Real,
        E::Real,
        heat::Ph,
        T_min::Real,
        T_max::Real;
        gate_width::Real = 1.0,
        validate_mass::Bool = true,
    ) where {NRin,NRout,Ph<:AbstractPropertyFunction}
        for (_, order) in reactants
            _validate_reaction_order(order, name)
        end
        for (_, y) in products
            _validate_yield(y, name)
        end
        if validate_mass && NRin == 1
            _validate_stoichiometry(products, name)
        end
        _validate_gate_width(gate_width, name)
        new{NRin,NRout,Ph}(name, reactants, products, A, E, heat, T_min, T_max, gate_width)
    end
end


# ============================================================================
# Validation Helpers
# ============================================================================

"""
	_validate_reaction_order(n::Real, name::Symbol)

Reject negative reaction orders at construction. In the ξ-power-law rate
`r = k(T)·ξ^n` a negative order makes the rate singular as the reactant
depletes (`ξ^n → ∞` as `ξ → 0`) — no depletion limiter can repair that.
Orders `n ≥ 0` (including 0 and fractional TGA fits) are supported; their
near-depletion behavior is regularized by the depletion limiter
([`DepletionLimiter`](@ref)).
"""
function _validate_reaction_order(n::Real, name::Symbol)
    n >= 0 || throw(
        ArgumentError(
            "Reaction :$name has reaction order n = $n; negative orders are " *
            "unphysical in the ξ-power-law rate (ξ^n → ∞ as the reactant " *
            "depletes). Supported orders are n ≥ 0.",
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

Validate that the product mass stoichiometries conserve mass for a
single-reactant reaction: per unit reaction extent the reactant contributes
mass 1.0, so mass conservation requires `Σ ν(products) == 1.0`.

Tolerance: 1e-10 (relative).

Called from the inner `Reaction` constructor for every `NRin == 1` reaction
unless `validate_mass=false` is passed (not recommended).
"""
function _validate_stoichiometry(
    products::NTuple{NRout,Tuple{Int,Float64}},
    name::Symbol,
) where {NRout}
    product_mass = 0.0
    for (_, ν) in products
        product_mass += ν
    end

    # Reactant mass is 1.0 by the single-reactant convention.
    tol = max(1e-10 * (1.0 + product_mass), 1e-14)

    if !isapprox(1.0, product_mass; atol = tol)
        error("""
           Reaction mass balance violated for reaction :$name
             Reactant mass: 1.0 (single-reactant convention)
             Product mass:  Σν = $product_mass
             Difference: $(1.0 - product_mass)

           Set `validate_mass=false` to disable this check (not recommended).
           """)
    end

    return nothing
end

# ============================================================================
# Convenience Constructors
# ============================================================================

"""
	Reaction(name, reactant => product; A, E, h, n=1.0, T_min=0.0, T_max=Inf, gate_width=1.0, validate_mass=true)

Create a simple single-reactant, single-product reaction.

# Arguments
- `name::Symbol`: Reaction identifier
- `reactant => product`: Pair of component indices
- `A`: Pre-exponential factor [s⁻¹]
- `E`: Activation energy [J/mol]
- `h`: Heat of reaction [J/kg] (number or property function)
- `n=1.0`: Reaction order
- `T_min=0.0`: Minimum temperature [K]
- `T_max=Inf`: Maximum temperature [K]
- `gate_width=1.0`: Width of the smooth tanh temperature gates [K] (> 0)
- `validate_mass=true`: If true, enforce mass balance in stoichiometry

# Example
```julia
rxn = Reaction(:decomp, 1 => 2, A=8.5e12, E=188e3, h=870e3)  # positive h = endothermic
```
"""
function Reaction(
    name::Symbol,
    pair::Pair{Int,Int};
    A::Real,
    E::Real,
    h,
    n::Real = 1.0,
    T_min::Real = 0.0,
    T_max::Real = Inf,
    gate_width::Real = 1.0,
    validate_mass::Bool = true,
)
    reactant, product = pair
    reactants = ((reactant, Float64(n)),)
    products = ((product, 1.0),)

    Reaction(
        name,
        reactants,
        products,
        Float64(A),
        Float64(E),
        to_property(h),
        Float64(T_min),
        Float64(T_max);
        gate_width = Float64(gate_width),
        validate_mass = validate_mass,
    )
end

"""
	Reaction(name, reactant => (product1, product2); A, E, h, yields, n=1.0, T_min=0.0, T_max=Inf, gate_width=1.0, validate_mass=true)

Create a reaction with two products.

# Arguments
- `name::Symbol`: Reaction identifier
- `reactant => (product1, product2)`: Reactant and product indices
- `A`: Pre-exponential factor [s⁻¹]
- `E`: Activation energy [J/mol]
- `h`: Heat of reaction [J/kg]
- `yields`: Tuple (y₁, y₂) of mass yields for each product
- `n=1.0`: Reaction order
- `T_min=0.0`: Minimum temperature [K]
- `T_max=Inf`: Maximum temperature [K]
- `gate_width=1.0`: Width of the smooth tanh temperature gates [K] (> 0)
- `validate_mass=true`: If true, enforce mass balance in stoichiometry

# Example
```julia
# Charring reaction: virgin -> 0.2 char + 0.8 gas (endothermic)
rxn = Reaction(:charring, 1 => (2, 3), A=1e15, E=200e3, h=500e3, yields=(0.2, 0.8))
```
"""
function Reaction(
    name::Symbol,
    pair::Pair{Int,Tuple{Int,Int}};
    A::Real,
    E::Real,
    h,
    yields::Tuple{Real,Real},
    n::Real = 1.0,
    T_min::Real = 0.0,
    T_max::Real = Inf,
    gate_width::Real = 1.0,
    validate_mass::Bool = true,
)
    reactant = pair.first
    product1, product2 = pair.second
    y1, y2 = Float64(yields[1]), Float64(yields[2])
    reactants = ((reactant, Float64(n)),)
    products = ((product1, y1), (product2, y2))

    Reaction(
        name,
        reactants,
        products,
        Float64(A),
        Float64(E),
        to_property(h),
        Float64(T_min),
        Float64(T_max);
        gate_width = Float64(gate_width),
        validate_mass = validate_mass,
    )
end

# ============================================================================
# Symbolic-name PendingReaction (resolved by Material constructor)
# ============================================================================

"""
	PendingReaction{R, P, Ph}

A reaction whose reactants/products are still expressed as component **symbols**
rather than integer indices. Produced by the symbolic constructors:

```julia
Reaction(:decomp, :virgin => :MMA;       A=…, E=…, h=…)
Reaction(:char,   :virgin => (:char, :gas); A=…, E=…, h=…, yields=(0.2, 0.8))
```

`PendingReaction`s are resolved into `Reaction`s by `Material(...)` once the
components tuple (and therefore the `name → index` map) is known. Users do not
construct `PendingReaction` directly; they call `Reaction(...)` with symbol
pairs.

Type parameters:
- `R`: reactant payload (`Symbol` or `Tuple{Vararg{Symbol}}`)
- `P`: product payload (`Symbol` or `Tuple{Vararg{Symbol}}`)
- `Ph`: heat-of-reaction property type
"""
struct PendingReaction{R,P,Ph<:AbstractPropertyFunction} <: AbstractReaction
    name::Symbol
    reactant::R
    product::P
    A::Float64
    E::Float64
    heat::Ph
    n::Float64
    yields::Union{Nothing,Tuple{Float64,Float64}}
    T_min::Float64
    T_max::Float64
    gate_width::Float64
    validate_mass::Bool
end

"""
	Reaction(name::Symbol, reactant::Symbol => product::Symbol; A, E, h, n=1.0, T_min=0.0, T_max=Inf, gate_width=1.0, validate_mass=true)

Symbolic single-reactant, single-product reaction constructor. Returns a
`PendingReaction` that is resolved to a numeric `Reaction` inside the
`Material(...)` constructor against the supplied components tuple.

```julia
Reaction(:decomposition, :virgin => :MMA, A = 8.5e12, E = 188e3, h = 870e3)
```
"""
function Reaction(
    name::Symbol,
    pair::Pair{Symbol,Symbol};
    A::Real,
    E::Real,
    h,
    n::Real = 1.0,
    T_min::Real = 0.0,
    T_max::Real = Inf,
    gate_width::Real = 1.0,
    validate_mass::Bool = true,
)
    _validate_reaction_order(n, name)
    _validate_gate_width(gate_width, name)
    PendingReaction(
        name,
        pair.first,
        pair.second,
        Float64(A),
        Float64(E),
        to_property(h),
        Float64(n),
        nothing,
        Float64(T_min),
        Float64(T_max),
        Float64(gate_width),
        validate_mass,
    )
end

"""
	Reaction(name::Symbol, reactant::Symbol => (product1::Symbol, product2::Symbol);
			 A, E, h, yields, n=1.0, T_min=0.0, T_max=Inf, gate_width=1.0, validate_mass=true)

Symbolic two-product reaction constructor. Returns a `PendingReaction` resolved
inside `Material(...)`.

```julia
Reaction(:charring, :virgin => (:char, :gas),
	A = 1e15, E = 200e3, h = 500e3, yields = (0.2, 0.8))
```
"""
function Reaction(
    name::Symbol,
    pair::Pair{Symbol,Tuple{Symbol,Symbol}};
    A::Real,
    E::Real,
    h,
    yields::Tuple{Real,Real},
    n::Real = 1.0,
    T_min::Real = 0.0,
    T_max::Real = Inf,
    gate_width::Real = 1.0,
    validate_mass::Bool = true,
)
    _validate_reaction_order(n, name)
    _validate_gate_width(gate_width, name)
    yields_f64 = (Float64(yields[1]), Float64(yields[2]))
    PendingReaction(
        name,
        pair.first,
        pair.second,
        Float64(A),
        Float64(E),
        to_property(h),
        Float64(n),
        yields_f64,
        Float64(T_min),
        Float64(T_max),
        Float64(gate_width),
        validate_mass,
    )
end

"""
	resolve_reaction(rxn::AbstractReaction, name_to_index::Dict{Symbol,Int}) -> Reaction

Resolve any `PendingReaction` (with `Symbol` payloads) to a numeric `Reaction`
by looking the symbols up in `name_to_index`. A concrete `Reaction` passes
through unchanged.

Used by `Material(...)` after the components tuple is fixed; not part of the
user-facing API.
"""
@inline resolve_reaction(rxn::Reaction, ::Dict{Symbol,Int}) = rxn

function resolve_reaction(
    rxn::PendingReaction{Symbol,Symbol},
    name_to_index::Dict{Symbol,Int},
)
    r_idx = _lookup_component(name_to_index, rxn.reactant, rxn.name)
    p_idx = _lookup_component(name_to_index, rxn.product, rxn.name)
    return Reaction(
        rxn.name,
        r_idx => p_idx;
        A = rxn.A,
        E = rxn.E,
        h = rxn.heat,
        n = rxn.n,
        T_min = rxn.T_min,
        T_max = rxn.T_max,
        gate_width = rxn.gate_width,
        validate_mass = rxn.validate_mass,
    )
end

function resolve_reaction(
    rxn::PendingReaction{Symbol,Tuple{Symbol,Symbol}},
    name_to_index::Dict{Symbol,Int},
)
    r_idx = _lookup_component(name_to_index, rxn.reactant, rxn.name)
    p1_idx = _lookup_component(name_to_index, rxn.product[1], rxn.name)
    p2_idx = _lookup_component(name_to_index, rxn.product[2], rxn.name)
    yields =
        rxn.yields === nothing ?
        throw(ArgumentError("Reaction :$(rxn.name) with two products requires `yields`")) :
        rxn.yields
    return Reaction(
        rxn.name,
        r_idx => (p1_idx, p2_idx);
        A = rxn.A,
        E = rxn.E,
        h = rxn.heat,
        n = rxn.n,
        yields = yields,
        T_min = rxn.T_min,
        T_max = rxn.T_max,
        gate_width = rxn.gate_width,
        validate_mass = rxn.validate_mass,
    )
end

@inline function _lookup_component(
    name_to_index::Dict{Symbol,Int},
    sym::Symbol,
    rxn_name::Symbol,
)
    idx = get(name_to_index, sym, 0)
    idx == 0 && error("Reaction :$rxn_name references unknown component :$sym; \
                        not found in the material's components tuple")
    return idx
end

# ============================================================================
# Stoichiometry Helpers
# ============================================================================

"""
	mass_balance_error(rxn::Reaction)

Compute the mass balance error for a reaction (should be ≈ 0).

Convention: per unit reaction extent, the first reactant contributes mass 1.0;
all products contribute their per-product mass stoich. Multi-reactant reactions
are not supported in this 1.5D shape; if `NRin > 1` this returns NaN.
"""
function mass_balance_error(rxn::Reaction{NRin}) where {NRin}
    NRin == 1 || return NaN
    reactant_mass = 1.0
    product_mass = 0.0
    @inbounds for (_, ν) in rxn.products
        product_mass += ν
    end
    return abs(reactant_mass - product_mass)
end

"""
	verify_stoichiometry(rxn::Reaction; tol=1e-10)

Verify that reaction stoichiometry satisfies mass balance.
"""
function verify_stoichiometry(rxn::Reaction; tol::Real = 1e-10)
    err = mass_balance_error(rxn)
    if err > tol
        @warn "Reaction $(rxn.name) has mass balance error: $err"
        return false
    end
    return true
end
