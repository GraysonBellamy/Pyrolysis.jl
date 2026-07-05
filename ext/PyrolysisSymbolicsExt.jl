"""
PyrolysisSymbolicsExt — symbolic-derivative verification for mixing rules.

Loaded automatically when both `Pyrolysis` and `Symbolics` are present.
This extension is **test-only**: production code does not call it. Its
single deliverable is `verify_mixing_derivatives`, which derives
∂k_eff/∂ξ and ∂k_eff/∂T symbolically via `Symbolics.jl` for each mixing
rule, then compares against the package's AD-default path. Drift bugs in
mixing-rule derivatives surface as relative errors above tolerance.

Strategy:

1. Build the mixing-rule expression in `Symbolics` variables `ξ_i`, `T`.
2. Apply `Symbolics.jacobian` / `Symbolics.gradient` to get analytical
   derivatives.
3. Substitute concrete values, compare against `ForwardDiff.gradient`
   applied to the same closure backed by the package's mixing functions.
"""
module PyrolysisSymbolicsExt

using Pyrolysis
using Pyrolysis: verify_mixing_derivatives
using Pyrolysis.Properties: effective_conductivity
using Pyrolysis.Core: PARALLEL, SERIES, WEIGHTED, BRUGGEMAN
using Symbolics
using ForwardDiff

# ---------------------------------------------------------------------------
# Internal: maximum relative error between two equal-length vectors.
# ---------------------------------------------------------------------------
@inline function _max_rel_err(a::AbstractVector, b::AbstractVector; atol::Real = 1e-12)
    denom = max.(abs.(a), abs.(b), atol)
    return maximum(abs.((a .- b) ./ denom))
end

# ---------------------------------------------------------------------------
# Build a Symbolics expression for the (PARALLEL / SERIES / WEIGHTED) k_eff
# given concrete `k_j(T)` callables and a mixing β.
# ---------------------------------------------------------------------------
function _symbolic_k_eff(
    rule::Pyrolysis.Core.MixingRule,
    k_funcs::NTuple{NC,Any},
    β::Real,
) where {NC}

    @variables ξ[1:NC] T

    # Volume fractions assuming intrinsic densities = 1 (cancels in ratio
    # tests; for absolute checks we lift the user's intrinsic densities).
    v = [ξ[i] / sum(ξ) for i = 1:NC]
    k = [k_funcs[i](T) for i = 1:NC]

    if rule === PARALLEL
        k_eff = sum(v[i] * k[i] for i = 1:NC)
    elseif rule === SERIES
        k_eff = 1 / sum(v[i] / k[i] for i = 1:NC)
    elseif rule === WEIGHTED
        k_par = sum(v[i] * k[i] for i = 1:NC)
        k_ser = 1 / sum(v[i] / k[i] for i = 1:NC)
        k_eff = β * k_par + (1 - β) * k_ser
    else
        throw(
            ArgumentError(
                "BRUGGEMAN has no closed form; symbolic verification not supported. " *
                "Use forward-mode AD for that rule.",
            ),
        )
    end

    return ξ, T, k_eff
end

# ---------------------------------------------------------------------------
# Public method: verify mixing-rule derivatives against the AD-default path.
# ---------------------------------------------------------------------------
function Pyrolysis.verify_mixing_derivatives(
    material;
    rtol::Real = 1e-8,
    atol::Real = 1e-10,
    ξ_test = nothing,
    T_test::Real = 600.0,
)

    NC = length(material.components)
    ξ_vec =
        ξ_test === nothing ? [Float64(c.density(T_test)) for c in material.components] :
        Vector{Float64}(ξ_test)

    # AD reference: ForwardDiff over `effective_conductivity`.
    ad_grad_ξ = ForwardDiff.gradient(
        ξ -> effective_conductivity(NTuple{NC,Float64}(ξ), T_test, material),
        ξ_vec,
    )
    ad_grad_T = ForwardDiff.derivative(
        T -> effective_conductivity(NTuple{NC,Float64}(ξ_vec), T, material),
        T_test,
    )

    # Symbolic reference for the package's active rule. Only PARALLEL /
    # SERIES / WEIGHTED have closed forms here; BRUGGEMAN is skipped.
    rule = material.mixing.k.rule
    results = Dict{Symbol,Float64}()

    if rule === BRUGGEMAN
        results[:bruggeman] = NaN  # marker: AD-only
    else
        k_funcs = ntuple(i -> material.components[i].thermal_conductivity, NC)
        β = material.mixing.k.β
        ξsym, Tsym, k_eff = _symbolic_k_eff(rule, k_funcs, β)

        grad_ξ_sym = Symbolics.gradient(k_eff, collect(ξsym))
        grad_T_sym = Symbolics.derivative(k_eff, Tsym)

        subs = Dict([ξsym[i] => ξ_vec[i] for i = 1:NC]..., Tsym => T_test)
        sym_grad_ξ = Float64.(Symbolics.value.(substitute.(grad_ξ_sym, Ref(subs))))
        sym_grad_T = Float64(Symbolics.value(substitute(grad_T_sym, subs)))

        results[Symbol(rule)] = max(
            _max_rel_err(ad_grad_ξ, sym_grad_ξ; atol = atol),
            _max_rel_err([ad_grad_T], [sym_grad_T]; atol = atol),
        )
    end

    return NamedTuple{Tuple(keys(results))}(values(results))
end

end # module
