"""
Reaction kinetics.
"""

# Immutable default for the standalone `reaction_rate` entry point. The solver
# paths never use it: `evaluate_reactions!` reads `material.depletion_limiter`
# unless a caller passes its own limiter.
const DEFAULT_DEPLETION_LIMITER = DepletionLimiter()

# PERFORMANCE: Fast power function for common reaction orders (0, 1).
# Scalar-type-generic so AD passes carry Duals through.
@inline _fast_pow(x::Real, n::Real) = n == 1 ? x : (n == 0 ? one(x) : x^n)

# The factor of a first-order reactant, multiplied onto `conc`: the finite base
# times the onset, which is zero at and below zero.
@inline function _reactant_factor(conc, limiter::DepletionLimiter, ξ::Real)
    if limiter.enabled
        # The onset and its slope vanish at and below zero, so the rate and
        # every derivative of it are exactly zero there.
        ξ <= 0 && return zero(conc) * zero(ξ)
        conc *= _finite_base(ξ, limiter.base_scale)
        conc *= one_sided_factor(ξ / limiter.threshold, limiter.join)
        return conc
    end
    return conc * _finite_power(ξ, limiter.base_scale, 1.0)
end

# The factor of a reactant of order other than one, multiplied onto `conc`: the
# reference concentration times the power of the finite base of the conversion
# `x = ξ/ξ_ref`, times the onset on `x`, which is zero at and below zero mass: it
# gates order zero (`x⁰ ≡ 1`) and keeps a finite derivative at depletion for
# orders below one half.
@inline function _conversion_factor(
    conc,
    limiter::DepletionLimiter,
    ξ::Real,
    ξ_ref::Real,
    order::Real,
)
    if limiter.enabled
        ξ <= 0 && return zero(conc) * zero(ξ)
        x = ξ / ξ_ref
        conc *= ξ_ref * _fast_pow(_finite_base(x, limiter.conversion_base), order)
        conc *= one_sided_factor(x / limiter.conversion_threshold, limiter.join)
        return conc
    end
    return conc * ξ_ref * _finite_power(ξ / ξ_ref, limiter.conversion_base, order)
end

# δ·softplus(ξ/δ) for ξ > 0: equal to ξ to within δ·log 2, and exactly ξ once
# ξ/δ is large enough for the correction to round away.
@inline _finite_base(ξ::Real, δ::Float64) = ξ + δ * log1p(exp(-ξ / δ))

# (δ·softplus(ξ/δ))^order for every real ξ. Below zero the power is formed in
# log space, where log softplus(x) = x once exp(x) is below round-off, so the
# value and its slope underflow together instead of meeting as 0 · ∞.
@inline function _finite_power(ξ::Real, δ::Float64, order::Real)
    x = ξ / δ
    x > 0 && return _fast_pow(_finite_base(ξ, δ), order)
    order == 0 && return one(ξ)
    log_softplus = x < -36 ? x : log(log1p(exp(x)))
    return exp(order * (log(δ) + log_softplus))
end

"""
	one_sided_factor(x, join) -> Real

The onset of [`DepletionLimiter`](@ref): zero at and below zero, the cubic
`(s − 2y) t³ + (3y − s) t²` (`t = x/join`, `y = tanh(join)`,
`s = join·sech²(join)`) up to `join`, and `tanh(x)` above. C¹ everywhere.
"""
@inline function one_sided_factor(x::Real, join::Float64)
    x <= 0 && return zero(x)
    x >= join && return x < 20.0 ? tanh(x) : one(x)
    t = x / join
    y = tanh(join)
    s = join / cosh(join)^2
    return t * t * ((s - 2y) * t + (3y - s))
end

"""
	reaction_rate(rxn::Reaction, ξ::NTuple{NC,Float64}, T::Float64;
				  limiter=DEFAULT_DEPLETION_LIMITER) -> Float64

Compute the reaction rate [kg/(m³·s)] of a first-order reaction.

r = A * exp(-E/RT) * ξ * g(ξ/ξ_th) * H(T - T_min) * H(T_max - T)

Any other order needs a reference: `reaction_rate(rxn, ξ, T, ξ_ref)`; this method
throws an `ArgumentError` for it.

# Depletion Protection
The reactant contributes a rate factor set by the [`DepletionLimiter`](@ref)
(`ξ_th = limiter.threshold`, default 1.0 kg/m³): the power of a finite base
times a C¹ onset `g` ([`one_sided_factor`](@ref)) that is exactly zero at and
below ξ = 0 and equals `tanh(ξ/ξ_th)` above `join · ξ_th`, so the rate recovers
the bare Arrhenius law for ξ ≫ ξ_th.

Pass a limiter via the `limiter` keyword to override per-call.

# Temperature Gates
The hard steps `H(·)` are realized as smooth tanh ramps of width
`rxn.gate_width` [K]: `0.5·(1 + tanh((T − T_min)/gate_width))` and the
mirrored `T_max` form. The default width is 1.0 K; smaller widths approach the
hard step (all Jacobian couplings differentiate through this function with
ForwardDiff, so the width needs no derivative threading).
"""
function reaction_rate(
    rxn::Reaction,
    ξ::NTuple{NC,<:Real},
    T::Real,
    ξ_ref::Real = NaN;
    limiter::DepletionLimiter = DEFAULT_DEPLETION_LIMITER,
) where {NC}
    k = rxn.A * exp(-rxn.E / (R_GAS * T))
    conc = one(promote_type(eltype(ξ), typeof(T), typeof(ξ_ref)))
    @inbounds ξ_i = ξ[rxn.reactant]
    rxn.order == 1 &&
        return _temperature_gates(k * _reactant_factor(conc, limiter, ξ_i), rxn, T)
    (isfinite(ξ_ref) && ξ_ref > 0) || _refuse_reference(rxn.name, ξ_ref)
    return _temperature_gates(
        k * _conversion_factor(conc, limiter, ξ_i, ξ_ref, rxn.order),
        rxn,
        T,
    )
end

# Smooth replacement for H(T - T_min) * H(T_max - T): a hard step would
# break Newton/AD smoothness, so a tanh ramp of width rxn.gate_width [K]
# keeps the bounds effective (≈1 inside, ≈0 outside) while remaining
# differentiable.
@inline function _temperature_gates(r, rxn::Reaction, T::Real)
    if isfinite(rxn.T_min) && rxn.T_min > 0.0
        r *= 0.5 * (1 + tanh((T - rxn.T_min) / rxn.gate_width))
    end
    if isfinite(rxn.T_max)
        r *= 0.5 * (1 + tanh((rxn.T_max - T) / rxn.gate_width))
    end
    return r
end

@noinline _refuse_reference(name::Symbol, ξ_ref) = throw(
    ArgumentError(
        "reaction :$name needs a finite positive reference concentration, got $ξ_ref",
    ),
)

"""
	evaluate_reactions!(rates, heats, ξ, T, material,
	                    limiter = material.depletion_limiter, ξ_ref = NaNs)

Evaluate all reaction rates [kg/(m³·s)] and heats of reaction [J/kg] at the mass
concentrations `ξ` [kg/m³] and temperature `T` [K], allocation-free.

Rates are evaluated with `limiter`, by default the material's own; see
[`DepletionLimiter`](@ref). `ξ_ref` holds the reference concentrations, one per
reaction, not read at order one. A reaction of another order with a zero
reference has no lineage in the cell: its rate and heat are zero whatever `ξ` and
the limiter hold. With no reference given, such a reaction is refused.

Heats of reaction can depend on temperature only or on both temperature and
concentration state (e.g., a moisture-dependent heat of sorption), through a
[`StateDependentProperty`](@ref).
"""
function evaluate_reactions!(
    rates::AbstractVector,
    heats::AbstractVector,
    ξ::NTuple{NC,<:Real},
    T::Real,
    material::Material{NC,NR},
    limiter::DepletionLimiter = material.depletion_limiter,
    ξ_ref::NTuple{NR,<:Real} = ntuple(_ -> NaN, Val(NR)),
) where {NC,NR}
    _evaluate_reactions_gen!(
        rates,
        heats,
        ξ,
        T,
        material.reactions,
        limiter,
        ξ_ref,
        Val(NR),
    )
    return nothing
end

@generated function _evaluate_reactions_gen!(
    rates::AbstractVector,
    heats::AbstractVector,
    ξ::NTuple{NC,<:Real},
    T::Real,
    reactions::Tuple{Vararg{Any,NR}},
    limiter::DepletionLimiter,
    ξ_ref::NTuple{NR,<:Real},
    ::Val{NR},
) where {NC,NR}
    exprs = Expr[]
    for i = 1:NR
        push!(
            exprs,
            quote
                @inbounds begin
                    rxn = reactions[$i]
                    if iszero(ξ_ref[$i]) && rxn.order != 1
                        rates[$i] = zero(eltype(rates))
                        heats[$i] = zero(eltype(heats))
                    else
                        rates[$i] = reaction_rate(rxn, ξ, T, ξ_ref[$i]; limiter = limiter)
                        heats[$i] = rxn.heat(T, ξ)
                    end
                end
            end,
        )
    end
    push!(exprs, :(return nothing))
    return Expr(:block, exprs...)
end

"""
	species_source_terms!(sources, rates, material)

Compute species source terms from reaction rates into `sources`:

dξⱼ/dt = Σᵢ νᵢⱼ rᵢ

where νᵢⱼ is positive for products, negative for reactants. The reaction
loop is unrolled at compile time (allocation-free).

# Arguments
- `sources`: Output vector of species source terms [kg/(m³·s)]
- `rates`: Reaction rates [kg/(m³·s)]
- `material`: Material definition
"""
function species_source_terms!(
    sources::AbstractVector,
    rates::AbstractVector,
    material::Material{NC,NR},
) where {NC,NR}
    fill!(sources, zero(eltype(sources)))
    _species_source_terms_gen!(sources, rates, material.reactions, Val(NR))
    return nothing
end

# Per-reaction source application: per unit reaction extent the reactant
# contributes mass 1.0 and each product its yield.
@inline function _apply_reaction_source!(src, rxn::Reaction, r)
    @inbounds src[rxn.reactant] -= r
    @inbounds for (idx, ν) in rxn.products
        src[idx] += ν * r
    end
    return nothing
end

# @generated version for allocation-free iteration over reactions
@generated function _species_source_terms_gen!(
    sources::AbstractVector,
    rates::AbstractVector,
    reactions::Tuple{Vararg{Any,NR}},
    ::Val{NR},
) where {NR}
    exprs = Expr[]
    for i = 1:NR
        push!(exprs, quote
            @inbounds begin
                rxn = reactions[$i]
                _apply_reaction_source!(sources, rxn, rates[$i])
            end
        end)
    end
    push!(exprs, :(return nothing))
    return Expr(:block, exprs...)
end

"""
	reaction_heat_source(rates::Vector{Float64}, heats::Vector{Float64}) -> Float64

Compute total heat source from reactions [W/m³].

Q = -Σᵢ hᵢ rᵢ

The negative sign converts from enthalpy convention (h < 0 exothermic) to
heat source convention (Q > 0 heats the system).

For exothermic reactions: h < 0, r > 0, so Q = -h*r > 0 (heating)
For endothermic reactions: h > 0, r > 0, so Q = -h*r < 0 (cooling)
"""
@inline function reaction_heat_source(rates::AbstractVector, heats::AbstractVector)
    Tu = promote_type(eltype(rates), eltype(heats))
    Q = zero(Tu)
    @inbounds for i in eachindex(rates, heats)
        Q -= heats[i] * rates[i]  # Note the negative sign!
    end
    return Q
end
