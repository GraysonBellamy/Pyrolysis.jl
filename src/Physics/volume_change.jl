"""
Volume change computations.

PERFORMANCE NOTE: All functions that iterate over material.components are
generated, to avoid allocations from indexing heterogeneous tuples with runtime
indices.
"""

# ============================================================================
# Volume change rate helpers (allocation-free)
# ============================================================================
#
# `is_gas` resolves from the component type, so a gas term is folded away at
# compile time and its reference density is never evaluated.

# Unrolled at compile time, so the iteration does not allocate.
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
        ρ_sym = Symbol("ρ_", i)
        push!(exprs, quote
            @inbounds begin
                if !is_gas(comps[$i])
                    $ρ_sym = comps[$i].density(T)
                    θ += dξ[$i] / $ρ_sym
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

dV/dt / V = Σⱼ (dξⱼ/dt) / ρⱼ, summed over the condensed components

A condensed (solid or liquid) component holds the volume ξⱼ/ρⱼ of its bulk
density, so its source changes the cell volume by dξⱼ/ρⱼ; the same bulk
occupancy weights the property mixing (`Properties.condensed_volume_fractions`).
Gases never enter the strain.

PERFORMANCE: Uses dispatch-based helpers to avoid allocations from heterogeneous tuple access.
"""
@inline function volume_change_rate(dξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    return _volume_change_rate(dξ, T, material.components, Val(NC))
end
