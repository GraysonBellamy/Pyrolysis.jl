# ============================================================================
# Variable cross section
# ============================================================================

# Half-width of the band outside [0, 1] into which the column conversion is
# clamped before the lateral law: the clamp is the identity on every physical
# conversion and bounds the law's argument on trial states.
const LATERAL_CLAMP_WIDTH = 1e-3

"""
	LateralShrinkage{F,NC}

The variable cross section of a column whose area follows the material's
lateral shrinkage law, `A = A₀ law(χ̄)`, with `χ̄ = Σ Gᵢ / Σ Mᵢ` the column's
conversion: `Gᵢ` the cumulative dry pyrolysis gas cell `i` has produced, its
source the dry-gas reaction rates, and `Mᵢ` its initial dry solid mass. Two
algebraic rows close it:

	0 = A − A₀ law(χ̄_c),
	0 = θ_A − (d/dt) log law(χ̄_c),

the second the logarithmic rate of the first along `(Σ Ġ, Σ Ṁ)`, so the area
row's own rate is exactly `θ_A A`. `χ̄_c` is `χ̄` smoothly clamped into
`[−w, 1 + w]` (`w = LATERAL_CLAMP_WIDTH`), the identity on `[0, 1]`, and `Σ M`
is smoothly floored at `M_floor`. Both rows read every cell's history pair,
the lateral one also its rates: they are the system's two dense rows.

`dry_gas` flags the dry pyrolysis gases (every gas but water vapour).
"""
struct LateralShrinkage{F,NC}
    law::F
    A₀::Float64
    dry_gas::NTuple{NC,Bool}
    M_floor::Float64
end

"""
	LateralShrinkage(material, A₀, M₀)

The cross section of initial area `A₀` [m²] of a column holding the initial
dry solid mass `M₀` [kg]. The material's law must be positive and finite on
`[−w, 1 + w]` and equal one at zero conversion, and the column must hold dry
solid.
"""
function LateralShrinkage(material::Material{NC}, A₀::Real, M₀::Real) where {NC}
    law = material.lateral_shrinkage_law
    law === nothing && throw(ArgumentError("the material has no lateral shrinkage law"))
    M₀ > 0 || throw(
        ArgumentError(
            "a lateral shrinkage law needs dry solid mass: the column's conversion " *
            "is the dry gas produced over the initial dry solid mass, which is zero",
        ),
    )
    w = LATERAL_CLAMP_WIDTH
    for χ in range(-w, 1 + w; length = 201)
        a = law(χ)
        (isfinite(a) && a > 0) || throw(
            ArgumentError(
                "the lateral shrinkage law must be positive and finite on " *
                "[$(-w), $(1 + w)]; law($χ) = $a",
            ),
        )
    end
    abs(law(0.0) - 1) <= 1e-12 || throw(
        ArgumentError(
            "the lateral shrinkage law must equal one at zero conversion; law(0) = " *
            "$(law(0.0))",
        ),
    )
    dry = dry_pyrolysis_gas_indices(material)
    return LateralShrinkage(law, Float64(A₀), ntuple(j -> j in dry, Val(NC)), 1e-8 * M₀)
end

@inline _clamped_conversion(χ) =
    smooth_clamp(χ, -LATERAL_CLAMP_WIDTH, 1 + LATERAL_CLAMP_WIDTH, LATERAL_CLAMP_WIDTH)

"""
	chi_bar(lateral, ΣG, ΣM)

The column conversion `Σ G / Σ M`, the denominator smoothly floored.
"""
@inline chi_bar(lat::LateralShrinkage, ΣG, ΣM) = ΣG / smooth_floor(ΣM, lat.M_floor)

@inline _log_area_ratio(lat::LateralShrinkage, ΣG, ΣM) =
    log(lat.law(_clamped_conversion(chi_bar(lat, ΣG, ΣM))))

# The logarithmic rate of the area along the history rates.
@inline _lateral_rate(lat::LateralShrinkage, ΣG, ΣM, ΣĠ, ΣṀ) =
    ForwardDiff.derivative(τ -> _log_area_ratio(lat, ΣG + τ * ΣĠ, ΣM + τ * ΣṀ), 0.0)

"""
	lateral_rows(lateral, s) -> SVector{2}

The residuals of the area row and the lateral-rate row at the column sums
`s = (Σ G, Σ M, Σ Ġ, Σ Ṁ, θ_A, A)`. Every dependence of the two rows on the
state passes through these six numbers, which is how their Jacobian rows are
assembled.
"""
@inline function lateral_rows(lat::LateralShrinkage, s::SVector{6})
    ΣG, ΣM, ΣĠ, ΣṀ, θ_A, A = s
    χ = _clamped_conversion(chi_bar(lat, ΣG, ΣM))
    return SVector(A - lat.A₀ * lat.law(χ), θ_A - _lateral_rate(lat, ΣG, ΣM, ΣĠ, ΣṀ))
end

"""
	consistent_lateral(lateral, ΣG, ΣM, ΣĠ, ΣṀ) -> (A, θ_A)

The area and lateral rate that satisfy both rows.
"""
function consistent_lateral(lat::LateralShrinkage, ΣG, ΣM, ΣĠ, ΣṀ)
    χ = _clamped_conversion(chi_bar(lat, ΣG, ΣM))
    return lat.A₀ * lat.law(χ), _lateral_rate(lat, ΣG, ΣM, ΣĠ, ΣṀ)
end

# History sources: the dry gas each cell's reactions produce, and none for
# the initial dry solid mass.
@inline _history_sources!(c, ::Nothing, topo) = nothing

function _history_sources!(
    c::Contributions{NC},
    lat::LateralShrinkage{F,NC},
    topo::MeshTopology,
) where {F,NC}
    @inbounds for i = 1:n_cells(topo)
        S = zero(eltype(c.S_carried))
        for j = 1:NC
            lat.dry_gas[j] && (S += c.S_mass[j, i])
        end
        c.S_carried[1, i] = S
        c.S_carried[2, i] = zero(S)
    end
    return nothing
end
