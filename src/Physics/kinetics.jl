"""
Reaction kinetics.
"""

# The limiter types are defined in Materials/reactions.jl (each Material carries
# its own limiter in `material.depletion_limiter`); re-exported from Physics
# for the kinetics-facing API.

# Immutable default for the standalone `reaction_rate` entry point. The solver
# paths never use it: `evaluate_reactions!` reads `material.depletion_limiter`
# unless a caller passes its own limiter.
const DEFAULT_DEPLETION_LIMITER = DepletionLimiter()

# PERFORMANCE: Fast power function for common reaction orders (0, 1).
# Scalar-type-generic so AD passes carry Duals through.
@inline _fast_pow(x::Real, n::Real) = n == 1 ? x : (n == 0 ? one(x) : x^n)

# The concentration factor of one reactant, multiplied onto `conc`: the power
# of the finite base times the onset, which is zero at and below zero.
@inline function _reactant_factor(conc, limiter::DepletionLimiter, ξ::Real, order::Real)
    if limiter.enabled
        # The onset and its slope vanish at and below zero, so the rate and
        # every derivative of it are exactly zero there.
        ξ <= 0 && return zero(conc) * zero(ξ)
        conc *= _fast_pow(_finite_base(ξ, limiter.base_scale), order)
        conc *= one_sided_factor(ξ / limiter.threshold, limiter.join)
        return conc
    end
    return conc * _finite_power(ξ, limiter.base_scale, order)
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

Compute the reaction rate [kg/(m³·s)] for a single reaction.

r = A * exp(-E/RT) * ∏ᵢ ξᵢ^nᵢ * ∏ᵢ g(ξᵢ/ξ_th) * H(T - T_min) * H(T_max - T)

# Depletion Protection
Each reactant contributes a rate factor set by the [`DepletionLimiter`](@ref)
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
    rxn::Reaction{NRin},
    ξ::NTuple{NC,<:Real},
    T::Real;
    limiter::AbstractDepletionLimiter = DEFAULT_DEPLETION_LIMITER,
) where {NRin,NC}
    # Arrhenius rate constant
    k = rxn.A * exp(-rxn.E / (R_GAS * T))

    # Concentration product over reactants: ∏_i ξ_{r_i}^{n_i}, each with a
    # smooth tanh depletion factor on the *rate*, not the power-law base — a
    # base-applied limiter cannot gate order-0 reactants (ξ_eff^0 ≡ 1) and has
    # singular ∂r/∂ξ at depletion for orders n < 0.5. Tuple iteration is fully
    # unrolled by the compiler when NRin is a static type parameter.
    conc_factor = one(promote_type(eltype(ξ), typeof(T)))
    @inbounds for (idx, order) in rxn.reactants
        conc_factor = _reactant_factor(conc_factor, limiter, ξ[idx], order)
    end

    r = k * conc_factor

    # Smooth replacement for H(T - T_min) * H(T_max - T): a hard step would
    # break Newton/AD smoothness, so a tanh ramp of width rxn.gate_width [K]
    # keeps the bounds effective (≈1 inside, ≈0 outside) while remaining
    # differentiable.
    if isfinite(rxn.T_min) && rxn.T_min > 0.0
        r *= 0.5 * (1 + tanh((T - rxn.T_min) / rxn.gate_width))
    end
    if isfinite(rxn.T_max)
        r *= 0.5 * (1 + tanh((rxn.T_max - T) / rxn.gate_width))
    end

    return r
end

"""
	evaluate_reactions!(rates::Vector{Float64}, heats::Vector{Float64},
					   ξ::NTuple{NC,Float64}, T::Float64, material::Material,
					   limiter = material.depletion_limiter)

Evaluate all reaction rates and heats of reaction.

# Arguments
- `rates`: Output vector of reaction rates [kg/(m³·s)]
- `heats`: Output vector of heats of reaction [J/kg]
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material definition

This function is designed to be allocation-free using @generated for compile-time unrolling.

# Depletion Limiter
Rates are evaluated with `limiter`, by default the material's own
(`material.depletion_limiter`, settable via the `Material` constructor's
`depletion_limiter` keyword); see [`DepletionLimiter`](@ref).

# State-Dependent Heats
Heats of reaction can depend on temperature only (standard) or on both
temperature and concentration state (e.g., moisture-dependent heat of sorption).
The concentration tuple ξ is passed to allow StateDependentProperty evaluation.
"""
function evaluate_reactions!(
    rates::AbstractVector,
    heats::AbstractVector,
    ξ::NTuple{NC,<:Real},
    T::Real,
    material::Material{NC,NR},
    limiter::AbstractDepletionLimiter = material.depletion_limiter,
) where {NC,NR}
    _evaluate_reactions_gen!(rates, heats, ξ, T, material.reactions, limiter, Val(NR))
    return nothing
end

# @generated version for allocation-free iteration over reactions
# Note: heats are evaluated with (T, ξ) to support StateDependentProperty
@generated function _evaluate_reactions_gen!(
    rates::AbstractVector,
    heats::AbstractVector,
    ξ::NTuple{NC,<:Real},
    T::Real,
    reactions::Tuple{Vararg{Any,NR}},
    limiter::AbstractDepletionLimiter,
    ::Val{NR},
) where {NC,NR}
    exprs = Expr[]
    for i = 1:NR
        push!(exprs, quote
            @inbounds begin
                rxn = reactions[$i]
                rates[$i] = reaction_rate(rxn, ξ, T; limiter = limiter)
                heats[$i] = rxn.heat(T, ξ)  # Pass ξ for state-dependent heats
            end
        end)
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

# Per-reaction source application. Convention: per unit reaction extent the
# first (only) reactant contributes mass 1.0; products contribute their
# per-product `mass_stoich`. Multi-reactant (NRin > 1) is not supported in the
# current shape — the per-reactant tuple would need extending to (idx, order,
# mass_stoich).
@inline function _apply_reaction_source!(src, rxn::Reaction{NRin}, r) where {NRin}
    NRin == 1 || error(
        "Multi-reactant reaction $(rxn.name) (NRin=$NRin) requires per-reactant mass stoich; not yet supported",
    )
    @inbounds src[rxn.reactants[1][1]] -= r
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
