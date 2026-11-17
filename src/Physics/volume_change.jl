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
# compile time and no gas density is read.

# The bulk density of a condensed component: a constant, since `setup` refuses
# a density law on a moving mesh.
@inline _bulk_density(p::ConstantProperty) = p.value

# Unrolled at compile time, so the iteration does not allocate.
@generated function _volume_change_rate(
    dξ,
    comps::Tuple{Vararg{Any,NC}},
    ::Val{NC},
) where {NC}
    exprs = Expr[]
    push!(exprs, :(θ = zero(eltype(dξ))))
    for i = 1:NC
        push!(exprs, quote
            @inbounds if !is_gas(comps[$i])
                θ += dξ[$i] / _bulk_density(comps[$i].density)
            end
        end)
    end
    push!(exprs, :(return θ))
    return Expr(:block, exprs...)
end

"""
	volume_change_rate(dξ::NTuple{NC,<:Real}, material::Material) -> Real

Compute relative volume change rate [(m³/m³)/s].

dV/dt / V = Σⱼ (dξⱼ/dt) / ρⱼ, summed over the condensed components

A condensed (solid or liquid) component holds the volume ξⱼ/ρⱼ of its bulk
density, so its source changes the cell volume by dξⱼ/ρⱼ; the same bulk
occupancy weights the property mixing (`Properties.condensed_volume_fractions`).
Gases never enter the strain. Each ρⱼ is a constant: a density that varied
with temperature would also change the volume through thermal dilation, a
term the strain rate does not carry, so `setup` refuses a condensed density
law on a moving mesh.

PERFORMANCE: Uses dispatch-based helpers to avoid allocations from heterogeneous tuple access.
"""
@inline function volume_change_rate(dξ::NTuple{NC,<:Real}, material) where {NC}
    return _volume_change_rate(dξ, material.components, Val(NC))
end
