"""
Reaction kinetics.
"""

# DepletionLimiter is defined in Materials/reactions.jl (each Material carries
# its own limiter in `material.depletion_limiter`); re-exported from Physics
# for the kinetics-facing API.

# Immutable default for the standalone `reaction_rate` entry point. The solver
# path never uses it — `evaluate_reactions!` reads `material.depletion_limiter`.
const DEFAULT_DEPLETION_LIMITER = DepletionLimiter()

# PERFORMANCE: Fast power function for common reaction orders (0, 1).
# Scalar-type-generic so AD passes carry Duals through.
@inline _fast_pow(x::Real, n::Real) = n == 1 ? x : (n == 0 ? one(x) : x^n)

# tanh(x) rounds to exactly 1.0 in Float64 for x ≳ 19.1, so skipping the call
# above 20·threshold is bit-identical to evaluating it. PERFORMANCE: bulk
# concentrations sit orders of magnitude above the default 1 kg/m³ threshold,
# so the common path pays no transcendental.
@inline _depletion_factor(ξ::Real, threshold::Float64) =
    ξ < 20.0 * threshold ? tanh(ξ / threshold) : one(ξ)

"""
	reaction_rate(rxn::Reaction, ξ::NTuple{NC,Float64}, T::Float64;
				  limiter=DEFAULT_DEPLETION_LIMITER) -> Float64

Compute the reaction rate [kg/(m³·s)] for a single reaction.

r = A * exp(-E/RT) * ∏ᵢ ξᵢ^nᵢ * ∏ᵢ tanh(ξᵢ/ξ_th) * H(T - T_min) * H(T_max - T)

# Depletion Protection
Each reactant contributes a smooth `tanh(ξ/ξ_th)` rate factor
(`ξ_th = limiter.threshold`, default 1.0 kg/m³): the rate recovers the bare
Arrhenius law for ξ ≫ ξ_th and vanishes with bounded ∂r/∂ξ as ξ → 0 for every
reaction order n ≥ 0. See [`DepletionLimiter`](@ref).

Pass a custom `DepletionLimiter` via the `limiter` keyword to override per-call.

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
    limiter::DepletionLimiter = DEFAULT_DEPLETION_LIMITER,
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
        ξ_raw = max(ξ[idx], zero(eltype(ξ)))
        conc_factor *= _fast_pow(ξ_raw, order)
        if limiter.enabled
            conc_factor *= _depletion_factor(ξ_raw, limiter.threshold)
        end
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
					   ξ::NTuple{NC,Float64}, T::Float64, material::Material)

Evaluate all reaction rates and heats of reaction.

# Arguments
- `rates`: Output vector of reaction rates [kg/(m³·s)]
- `heats`: Output vector of heats of reaction [J/kg]
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material definition

This function is designed to be allocation-free using @generated for compile-time unrolling.

# Depletion Limiter
Rates are evaluated with the material's own [`DepletionLimiter`](@ref)
(`material.depletion_limiter`, settable via the `Material` constructor's
`depletion_limiter` keyword).

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
) where {NC,NR}
    _evaluate_reactions_gen!(
        rates,
        heats,
        ξ,
        T,
        material.reactions,
        material.depletion_limiter,
        Val(NR),
    )
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
    limiter::DepletionLimiter,
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
