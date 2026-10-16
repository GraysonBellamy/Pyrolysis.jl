"""
Mixture property computations.

PERFORMANCE NOTE: All functions that iterate over material.components use
recursive dispatch to avoid allocations from indexing heterogeneous tuples
with runtime indices. Each component access uses a compile-time constant index.
"""

# Promoted working scalar type for property evaluation with a (possibly
# heterogeneous) tuple input and a temperature-like scalar. The analytical
# Jacobian operators feed `ForwardDiff.Dual`s through these helpers; computing
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
# rate θ (Physics/volume_change.jl) and the mixture density below. For
# canonical materials (γ = 0 gas, γ = 1 condensed) this reproduces the old
# binary include/exclude exactly; a partial-γ component now contributes the
# same volume fraction here as it does to θ and ρ_bulk.
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
This function shares it with the dilation rate θ = Σγⱼ(dξⱼ/dt)/ρⱼ and the
mixture density 1/ρ = Σγⱼ·Yⱼ/ρⱼ. Consequences:

- Canonical gas (γ = 0): excluded — matches ThermaKin "when γ = 0, gaseous
  components do not contribute to the thermal conductivity" (ThermaKin2Ds
  Reference, Eqs. 8–10).
- Canonical condensed (γ = 1): full ξ/ρ weight.
- Partial γ: contributes the same γ-fraction of its volume here as it does
  to θ and ρ_bulk. A γ = 0 condensed phase is consistently *volumeless*
  (no strain, no density volume, no conductivity weight).

Used for:
- Effective thermal conductivity
- Effective emissivity

NOT used for the gas-transfer coefficient λ: gas diffusion happens *through*
the pore space, so [`effective_gas_transfer`](@ref) mixes on the all-phase
occupancy basis [`volume_fractions`](@ref) (gas included), never this
matrix basis.
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
would hand the Jacobian entries of order `1/φ_m²` in a cell that a fixed
mesh cannot park.
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

# ============================================================================
# mixture_density with swelling
# ============================================================================

# Swelling-weighted specific volume of the bulk, mass fractions computed
# inline: v = Σ γ_j·Y_j/ρ_j. Same universal γ semantics as
# solid_volume_fractions and the θ kernel: γ applies to every component, with
# a γ = 0 component contributing zero volume.
@generated function _specific_vol_with_swelling_direct(
    ξ,
    inv_ρ,
    T,
    comps::Tuple,
    ::Val{NC},
) where {NC}
    expr = quote
        v = zero(_prop_type(ξ, T))
    end
    for i = 1:NC
        push!(expr.args, quote
            c = comps[$i]
            γ = c.swelling_factor
            if γ > 0.0
                Y_i = ξ[$i] * inv_ρ  # Compute mass fraction inline
                v += γ * Y_i / c.density(T)
            end
        end)
    end
    push!(expr.args, :(return v))
    return expr
end

"""
	mixture_density(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute bulk mixture density from mass concentrations.

The bulk specific volume uses the universal swelling weighting

	1/ρ = Σⱼ γⱼ·Yⱼ/ρⱼ

— each component contributes the γ-fraction of its intrinsic volume to the
bulk, exactly as it does to the dilation rate θ and to
[`solid_volume_fractions`](@ref) (one γ semantics). For
canonical materials (γ = 0 gas, γ = 1 condensed) this is the familiar
condensed-only specific volume.

PERFORMANCE: Uses optimized version that computes mass fractions inline to avoid
allocating intermediate NTuple.
"""
@inline function mixture_density(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    W = _prop_type(ξ, T)
    ρ_total = sum(ξ)
    if ρ_total < 1e-20
        return zero(W)
    end

    inv_ρ = 1.0 / ρ_total
    # Use optimized version that computes mass fractions inline
    v_specific =
        _specific_vol_with_swelling_direct(ξ, inv_ρ, T, material.components, Val(NC))

    if v_specific < 1e-20
        return zero(W)
    end

    return convert(W, 1.0 / v_specific)
end

# ============================================================================
# Porosity
# ============================================================================

"""
	porosity(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute porosity (void fraction): the bulk-volume fraction not occupied by the
solid/liquid skeletal frame. With component **intrinsic** (skeletal) densities
ρⱼ_intrinsic:

	φ = clamp(1 - Σⱼ_solid_or_liquid (ξⱼ / ρⱼ_intrinsic), 0.0, 1.0)

The intrinsic density is the density of the cell-wall material itself,
ignoring pore volume — distinct from the bulk density used for thermal mixing
rules. When a component was constructed without an explicit `ρ_intrinsic` /
`φ` kwarg, `get_intrinsic_density` falls back to the bulk density, which
recovers the legacy "no-pores" definition.

The result enters the gas-pressure formula `P = N(ξ,T)/φ`. Derivatives are
taken by ForwardDiff through this function (the live structured Jacobian
differentiates the clamped expression directly — there is no separate
analytic-derivative path). The single canonical implementation used by
production property updates lives here.
"""
@inline function porosity(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    return _porosity_void_frac(ξ, T, material.components, Val(NC))
end

@generated function _porosity_void_frac(
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
    push!(expr.args, :(return clamp(1.0 - solid_vol, 0.0, 1.0)))
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
