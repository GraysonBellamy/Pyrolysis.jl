"""
Volume change computations.

PERFORMANCE NOTE: All functions that iterate over material.components use
dispatch-based unrolling (explicit overloads for small tuple sizes) to avoid
allocations from indexing heterogeneous tuples with runtime indices.
"""

# ============================================================================
# Dispatch-based volume change rate helpers (allocation-free)
# ============================================================================

@inline function _volume_change_rate(dξ, T, comps::Tuple{C}, ::Val{1}) where {C}
    W = promote_type(typeof(T), eltype(dξ))
    comp = comps[1]
    γ = comp.swelling_factor
    if γ > 0
        ρⱼ = comp.density(T)
        return convert(W, γ * dξ[1] / ρⱼ)
    end
    return zero(W)
end

@inline function _volume_change_rate(dξ, T, comps::Tuple{C1,C2}, ::Val{2}) where {C1,C2}
    θ = zero(promote_type(typeof(T), eltype(dξ)))
    γ1 = comps[1].swelling_factor
    if γ1 > 0
        ρ1 = comps[1].density(T)
        θ += γ1 * dξ[1] / ρ1
    end
    γ2 = comps[2].swelling_factor
    if γ2 > 0
        ρ2 = comps[2].density(T)
        θ += γ2 * dξ[2] / ρ2
    end
    return θ
end

@inline function _volume_change_rate(
    dξ,
    T,
    comps::Tuple{C1,C2,C3},
    ::Val{3},
) where {C1,C2,C3}
    θ = zero(promote_type(typeof(T), eltype(dξ)))
    γ1 = comps[1].swelling_factor
    if γ1 > 0
        ρ1 = comps[1].density(T)
        θ += γ1 * dξ[1] / ρ1
    end
    γ2 = comps[2].swelling_factor
    if γ2 > 0
        ρ2 = comps[2].density(T)
        θ += γ2 * dξ[2] / ρ2
    end
    γ3 = comps[3].swelling_factor
    if γ3 > 0
        ρ3 = comps[3].density(T)
        θ += γ3 * dξ[3] / ρ3
    end
    return θ
end

# @generated version for NC > 3 components (allocation-free via compile-time unrolling)
@generated function _volume_change_rate(
    dξ,
    T,
    comps::Tuple{Vararg{Any,NC}},
    ::Val{NC},
) where {NC}
    # Build the unrolled expression at compile time
    exprs = Expr[]
    push!(exprs, :(θ = zero(promote_type(typeof(T), eltype(dξ)))))
    for i = 1:NC
        γ_sym = Symbol("γ_", i)
        ρ_sym = Symbol("ρ_", i)
        push!(exprs, quote
            @inbounds begin
                $γ_sym = comps[$i].swelling_factor
                if $γ_sym > 0
                    $ρ_sym = comps[$i].density(T)
                    θ += $γ_sym * dξ[$i] / $ρ_sym
                end
            end
        end)
    end
    push!(exprs, :(return θ))
    return Expr(:block, exprs...)
end

"""
	volume_change_rate(dξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute relative volume change rate [(m³/m³)/s].

dV/dt / V = Σⱼ (γⱼ/ρⱼ) (dξⱼ/dt)

Every component is weighted by its swelling factor γⱼ — the same universal
γ semantics used by `Properties.mixture_density` (1/ρ = Σγⱼ·Yⱼ/ρⱼ) and
`Properties.solid_volume_fractions`: γⱼ is the fraction of
component j's intrinsic volume that contributes to the bulk volume.
Canonical values are γ = 0 (gas, escapes the matrix) and γ = 1 (condensed).

PERFORMANCE: Uses dispatch-based helpers to avoid allocations from heterogeneous tuple access.
"""
@inline function volume_change_rate(dξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    return _volume_change_rate(dξ, T, material.components, Val(NC))
end
