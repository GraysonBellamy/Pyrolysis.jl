"""
Mixture property computations.

PERFORMANCE NOTE: All functions that iterate over material.components use
recursive dispatch to avoid allocations from indexing heterogeneous tuples
with runtime indices. Each component access uses a compile-time constant index.
"""

# Promoted working scalar type for property evaluation with a (possibly
# heterogeneous) tuple input and a temperature-like scalar. The colored
# Jacobian feeds `ForwardDiff.Dual`s through these helpers; computing
# accumulator inits and guarded zero branches in this type keeps every helper
# type-stable. Literal-`0.0` zeros instead produce `Union{Float64, Dual}`
# locals and returns, which box and allocate on the hot Jacobian path.
@inline _prop_type(x::Tuple) = promote_type(map(typeof, x)...)
@inline _prop_type(x::Tuple, T) = promote_type(typeof(T), map(typeof, x)...)

@inline _zero_ntuple(::Type{W}, ::Val{N}) where {W,N} = ntuple(_ -> zero(W), Val(N))

# ============================================================================
# Recursive helpers for heterogeneous tuple iteration
# These avoid runtime indexing by using compile-time dispatch
# ============================================================================

# Compute specific volume from components recursively
@inline function _sum_specific_vol(ξ, T, comps::Tuple{C}, ::Val{1}) where {C}
    W = _prop_type(ξ, T)
    ρ = comps[1].density(T)
    return ξ[1] > 0.0 && ρ > 0.0 ? convert(W, ξ[1] / ρ) : zero(W)
end

@inline function _sum_specific_vol(ξ, T, comps::Tuple{C1,C2}, ::Val{2}) where {C1,C2}
    v = zero(_prop_type(ξ, T))
    ρ1 = comps[1].density(T)
    if ξ[1] > 0.0 && ρ1 > 0.0
        v += ξ[1] / ρ1
    end
    ρ2 = comps[2].density(T)
    if ξ[2] > 0.0 && ρ2 > 0.0
        v += ξ[2] / ρ2
    end
    return v
end

@inline function _sum_specific_vol(ξ, T, comps::Tuple{C1,C2,C3}, ::Val{3}) where {C1,C2,C3}
    v = zero(_prop_type(ξ, T))
    ρ1 = comps[1].density(T)
    if ξ[1] > 0.0 && ρ1 > 0.0

        v += ξ[1] / ρ1
    end
    ρ2 = comps[2].density(T)
    if ξ[2] > 0.0 && ρ2 > 0.0

        v += ξ[2] / ρ2
    end
    ρ3 = comps[3].density(T)
    if ξ[3] > 0.0 && ρ3 > 0.0

        v += ξ[3] / ρ3
    end
    return v
end

# Generic fallback for any number of components (NC > 3) using generated function
@inline function _sum_specific_vol(ξ, T, comps::Tuple, ::Val{NC}) where {NC}
    return _sum_specific_vol_gen(ξ, T, comps)
end

# Compute volume fractions tuple
@inline function _vol_fracs(ξ, T, inv_v_total, comps::Tuple{C}, ::Val{1}) where {C}
    W = _prop_type(ξ, T)
    ρ = comps[1].density(T)
    v = ξ[1] > 0.0 && ρ > 0.0 ? convert(W, ξ[1] / ρ * inv_v_total) : zero(W)
    return (v,)
end

@inline function _vol_fracs(ξ, T, inv_v_total, comps::Tuple{C1,C2}, ::Val{2}) where {C1,C2}
    W = _prop_type(ξ, T)
    ρ1 = comps[1].density(T)
    v1 = ξ[1] > 0.0 && ρ1 > 0.0 ? convert(W, ξ[1] / ρ1 * inv_v_total) : zero(W)
    ρ2 = comps[2].density(T)
    v2 = ξ[2] > 0.0 && ρ2 > 0.0 ? convert(W, ξ[2] / ρ2 * inv_v_total) : zero(W)
    return (v1, v2)
end

@inline function _vol_fracs(
    ξ,
    T,
    inv_v_total,
    comps::Tuple{C1,C2,C3},
    ::Val{3},
) where {C1,C2,C3}
    W = _prop_type(ξ, T)
    ρ1 = comps[1].density(T)
    v1 = ξ[1] > 0.0 && ρ1 > 0.0 ? convert(W, ξ[1] / ρ1 * inv_v_total) : zero(W)
    ρ2 = comps[2].density(T)
    v2 = ξ[2] > 0.0 && ρ2 > 0.0 ? convert(W, ξ[2] / ρ2 * inv_v_total) : zero(W)
    ρ3 = comps[3].density(T)
    v3 = ξ[3] > 0.0 && ρ3 > 0.0 ? convert(W, ξ[3] / ρ3 * inv_v_total) : zero(W)
    return (v1, v2, v3)
end

# Generic fallback for any number of components (NC > 3) using generated function
@inline function _vol_fracs(ξ, T, inv_v_total, comps::Tuple, ::Val{NC}) where {NC}
    return _vol_fracs_gen(ξ, T, inv_v_total, comps)
end

"""
	volume_fractions(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> NTuple{NC,Float64}

Compute volume fractions from mass concentrations and component densities.

vⱼ = (ξⱼ/ρⱼ) / Σᵢ(ξᵢ/ρᵢ)

Uses dispatch-based unrolling to avoid allocation from heterogeneous tuple indexing.
"""
@inline function volume_fractions(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    comps = material.components
    v_total = _sum_specific_vol(ξ, T, comps, Val(NC))

    if v_total < 1e-20
        return _zero_ntuple(_prop_type(ξ, T), Val(NC))
    end

    inv_v_total = 1.0 / v_total
    return _vol_fracs(ξ, T, inv_v_total, comps, Val(NC))
end

# ============================================================================
# Swelling-weighted (matrix) volume fractions
# One universal γ semantics: every component contributes
# v_j = γ_j·ξ_j/ρ_j of bulk volume — the same weighting used by the dilation
# rate θ (Physics/volume_change.jl). For
# canonical materials (γ = 0 gas, γ = 1 condensed) this is a binary
# include/exclude; a partial-γ component contributes the same volume fraction
# here as it does to θ.
# ============================================================================

# Sum of swelling-weighted specific volumes Σ γ_j·ξ_j/ρ_j
@generated function _sum_solid_specific_vol(ξ, T, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        v = zero(_prop_type(ξ, T))
    end
    for i = 1:NC
        push!(expr.args, quote
            c = comps[$i]
            γ = c.swelling_factor
            if γ > 0.0
                ρ = c.density(T)
                if ξ[$i] > 0.0 && ρ > 0.0
                    v += γ * ξ[$i] / ρ
                end
            end
        end)
    end
    push!(expr.args, :(return v))
    return expr
end

# Swelling-weighted volume fractions tuple v_j = γ_j·ξ_j/ρ_j / v_total
@generated function _solid_vol_fracs(ξ, T, inv_v_total, comps::Tuple, ::Val{NC}) where {NC}
    var_exprs = []
    for i = 1:NC
        sym = Symbol("v$i")
        push!(
            var_exprs,
            quote
                c = comps[$i]
                γ = c.swelling_factor
                if γ > 0.0
                    ρ = c.density(T)
                    $sym =
                        ξ[$i] > 0.0 && ρ > 0.0 ? convert(W, γ * ξ[$i] / ρ * inv_v_total) :
                        zero(W)
                else
                    $sym = zero(W)
                end
            end,
        )
    end

    tuple_symbols = [Symbol("v$i") for i = 1:NC]
    tuple_expr = Expr(:tuple, tuple_symbols...)

    return quote
        W = _prop_type(ξ, T)
        $(var_exprs...)
        return $tuple_expr
    end
end

"""
	solid_volume_fractions(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> NTuple{NC,Float64}

Compute swelling-weighted (matrix) volume fractions: vⱼ = γⱼ·ξⱼ/ρⱼ,
normalized to sum to 1 over the contributing components.

The swelling factor γ has **one** meaning everywhere: the fraction of a
component's intrinsic volume that contributes to the bulk (matrix) volume.
This function shares its weights with the dilation rate θ = Σγⱼ(dξⱼ/dt)/ρⱼ.
Skeletal occupancy and pore volume are separate quantities. Consequences:

- Canonical gas (γ = 0): excluded — matches ThermaKin "when γ = 0, gaseous
  components do not contribute to the thermal conductivity" (ThermaKin2Ds
  Reference, Eqs. 8–10).
- Canonical condensed (γ = 1): full ξ/ρ weight.
- Partial γ: contributes the same γ-fraction of its volume here as it does
  to θ. A γ = 0 condensed phase has no dilation or matrix conductivity weight;
  it still occupies skeletal volume.

This gas-exclusion limit does not identify our per-component gamma with
ThermaKin's mixture gas swelling factor; in particular its liquid volume is
not removed by `L SWELLING: 0` (ThermaKin Eqs. 3–4).

Used for:
- Effective thermal conductivity
- Effective emissivity

Gas transport uses `transport_matrix_state`, which excludes gases even when
a gas component has nonzero γ, and has its own explicit endpoint closure.
"""
@inline function solid_volume_fractions(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    comps = material.components
    v_total = _sum_solid_specific_vol(ξ, T, comps, Val(NC))

    if v_total < 1e-20
        return _zero_ntuple(_prop_type(ξ, T), Val(NC))
    end

    inv_v_total = 1.0 / v_total
    return _solid_vol_fracs(ξ, T, inv_v_total, comps, Val(NC))
end

# Matrix volume fraction window over which the mixing weights blend from the
# matrix basis to the all-phase occupancy. Above `MIXING_BLEND_HI` the cell
# is an ordinary porous solid; below `MIXING_BLEND_LO` it is a void.
const MIXING_BLEND_LO = 1.0e-4
const MIXING_BLEND_HI = 1.0e-3

"""
	matrix_volume_fraction(ξ::NTuple{NC,<:Real}, T::Real, material) -> Real

Matrix volume fraction `φ_m = Σⱼ γⱼ ξⱼ/ρⱼ(T)` of a cell: the share of its
volume held by the swelling-weighted condensed matrix. A cell with
`φ_m < MIXING_BLEND_LO` is a void to the property model
([`mixing_volume_fractions`](@ref)).
"""
@inline matrix_volume_fraction(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC} =
    _sum_solid_specific_vol(ξ, T, material.components, Val(NC))

"""
	mixing_volume_fractions(ξ::NTuple{NC,<:Real}, T::Real, material) -> NTuple{NC}

Volume-fraction weights for conductivity/emissivity mixing: the
swelling-weighted matrix fractions ([`solid_volume_fractions`](@ref)) while
matrix volume remains, blending smoothly into the all-phase occupancy
([`volume_fractions`](@ref)) as the matrix vanishes.

The blend weight is a cubic smoothstep of the matrix volume fraction
`φ_m = Σⱼ γⱼ ξⱼ/ρⱼ` over `[MIXING_BLEND_LO, MIXING_BLEND_HI]`; above the
window the weights are exactly the matrix fractions, below it exactly the
occupancy. Without the occupancy limit a fully converted cell gets
k_eff = ε_eff = 0 — a perfect non-emitting insulator — which the harmonic
face mean then propagates to both adjacent faces, and Kirchhoff α = ε kills
absorption too; with it the void cell degrades to the gas-mixture
conductivity (heat still conducts through the pore gas). The blend keeps
∂k_eff/∂ξ bounded through the transition: the matrix fractions are
normalised by `φ_m`, so their derivatives grow as `1/φ_m²` and a hard switch
would hand the Jacobian entries of order `1/φ_m²` in a cell that stays in
the system on a fixed mesh.
"""
@inline function mixing_volume_fractions(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    φ_m = matrix_volume_fraction(ξ, T, material)
    φ_m >= MIXING_BLEND_HI && return solid_volume_fractions(ξ, T, material)
    v_full = volume_fractions(ξ, T, material)
    φ_m <= MIXING_BLEND_LO && return v_full
    v_solid = solid_volume_fractions(ξ, T, material)
    s = (φ_m - MIXING_BLEND_LO) / (MIXING_BLEND_HI - MIXING_BLEND_LO)
    w = s * s * (3 - 2 * s)
    return ntuple(j -> w * v_solid[j] + (1 - w) * v_full[j], Val(NC))
end

"""
	skeletal_volume_fraction(ξ::NTuple{NC,<:Real}, T::Real, material) -> Real

The bulk-volume fraction occupied by the condensed (solid and liquid)
skeleton, `Σⱼ ξⱼ / ρⱼ_intrinsic(T)` over the condensed components, without a
clamp: `1 − skeletal_volume_fraction` is the pore fraction (porosity), which a
caller floors smoothly where it divides by it and reads unclamped for
admissibility. The intrinsic density is that of the cell-wall material itself,
ignoring pore volume; a component constructed without `ρ_intrinsic` or `φ`
falls back to its bulk density and contributes no pore volume.
"""
@inline function skeletal_volume_fraction(
    ξ::NTuple{NC,<:Real},
    T::Real,
    material,
) where {NC}
    return _skeletal_volume(ξ, T, material.components, Val(NC))
end

@generated function _skeletal_volume(
    ξ::NTuple{NC,<:Real},
    T::Real,
    comps::Tuple,
    ::Val{NC},
) where {NC}
    expr = quote
        solid_vol = zero(_prop_type(ξ, T))
    end
    for i = 1:NC
        push!(expr.args, quote
            c = comps[$i]
            if is_solid(c) || is_liquid(c)
                ρ_intr = get_intrinsic_density(c, T)
                if ρ_intr > 0.0
                    solid_vol += ξ[$i] / ρ_intr
                end
            end
        end)
    end
    push!(expr.args, :(return solid_vol))
    return expr
end

# ============================================================================
# Generated functions for allocation-free computation with many components
# These generate fully unrolled code at compile time, avoiding runtime
# indexing into heterogeneous tuples which would allocate.
# ============================================================================

"""
Generated function to compute sum of specific volumes.
Generates fully unrolled code like:
	v = 0.0
	ρ1 = comps[1].density(T); if ξ[1] > 0.0 && ρ1 > 0.0; v += ξ[1] / ρ1; end
	ρ2 = comps[2].density(T); if ξ[2] > 0.0 && ρ2 > 0.0; v += ξ[2] / ρ2; end
	...
	return v
"""
@generated function _sum_specific_vol_gen(
    ξ::NTuple{NC,<:Real},
    T::Real,
    comps::Tuple,
) where {NC}
    expr = quote
        v = zero(_prop_type(ξ, T))
    end
    for i = 1:NC
        push!(expr.args, quote
            ρ = comps[$i].density(T)
            if ξ[$i] > 0.0 && ρ > 0.0
                v += ξ[$i] / ρ
            end
        end)
    end
    push!(expr.args, :(return v))
    return expr
end

"""
Generated function to compute volume fractions tuple.
Generates fully unrolled code returning a compile-time sized tuple.
"""
@generated function _vol_fracs_gen(
    ξ::NTuple{NC,<:Real},
    T::Real,
    inv_v_total::Real,
    comps::Tuple,
) where {NC}
    # Generate variable declarations for each volume fraction
    var_exprs = []
    for i = 1:NC
        sym = Symbol("v$i")
        push!(
            var_exprs,
            quote
                ρ = comps[$i].density(T)
                $sym =
                    ξ[$i] > 0.0 && ρ > 0.0 ? convert(W, ξ[$i] / ρ * inv_v_total) : zero(W)
            end,
        )
    end

    # Generate tuple return
    tuple_symbols = [Symbol("v$i") for i = 1:NC]
    tuple_expr = Expr(:tuple, tuple_symbols...)

    return quote
        W = _prop_type(ξ, T)
        $(var_exprs...)
        return $tuple_expr
    end
end
