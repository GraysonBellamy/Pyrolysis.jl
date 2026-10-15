"""
Energy-conserving temperature merge.

Shared by the depleted-cell park (`cell_park.jl`) and the AMR sibling
coarsening in `Experimental/h_coarsening.jl`.
"""

using ..Geometry: Mesh1D, Node1D, Cell1D, Face1D, CellState
using ..Materials: AbstractMaterial, Material, SOLID, T_REF, is_gas

# ============================================================================
# Energy-Conserving Temperature Merge
# ============================================================================

# 4-point Gauss–Legendre nodes/weights on [-1, 1]. Exact for polynomial cₚ(T)
# up to degree 7; merges are rare callback events, so the quadrature cost is
# irrelevant.
const _GL4_X =
    (-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526)
const _GL4_W =
    (0.34785484513745385, 0.6521451548625461, 0.6521451548625461, 0.34785484513745385)

# ∫_{T_lo}^{T_hi} cₚ(T′) dT′ by 4-point Gauss–Legendre quadrature.
@inline function _integrate_cp(c_p, T_lo::Float64, T_hi::Float64)
    half = 0.5 * (T_hi - T_lo)
    mid = 0.5 * (T_hi + T_lo)
    s = 0.0
    @inbounds for k = 1:4
        s += _GL4_W[k] * Float64(c_p(mid + half * _GL4_X[k]))
    end
    return half * s
end

"""
	merge_temperatures_conservatively(V_a::Float64, state_a::CellState{NC},
									   V_b::Float64, state_b::CellState{NC},
									   material::Material{NC}) where NC -> Float64

Compute the merged temperature by conserving the condensed-phase sensible
enthalpy.

# Principle
The conserved functional is

```
H(T) = Σ_{j ∉ gas} m_j ∫_{T_REF}^{T} cₚ,j(T′) dT′ ,   m_j = V·ξ_j
```

and `T_merged` solves `H_merged(T_merged) = H_a(T_a) + H_b(T_b)`. Two choices
here are deliberate:

- **∫cₚ dT, not cₚ(T)·T.** For temperature-dependent cₚ the product form is
  not an enthalpy: with cₚ = a + bT it carries the b-term with coefficient 1
  instead of ½, a few-K systematic bias per merge when a hot surface cell
  merges into a cooler interior one. The integral form is exact for any cₚ(T)
  (4-point Gauss–Legendre per component). The datum `T_REF` cancels because
  each component's total mass is identical before and after the merge.
- **Gas components are excluded.** Formulation A stores thermal energy in the
  condensed matrix only (the gas is thermally quasi-steady), so gas heat
  capacity does not belong in the merge weighting.

The solve is a bracketed Newton iteration on `[min(T_a,T_b), max(T_a,T_b)]`;
`H` is strictly monotone (cₚ > 0) so the root is unique and always inside the
bracket.

# Arguments
- `V_a`: Volume of cell a [m³]
- `state_a`: Cell state for cell a (T, ξ)
- `V_b`: Volume of cell b [m³]
- `state_b`: Cell state for cell b (T, ξ)
- `material`: Material definition with component heat capacity functions

# Returns
Merged temperature [K] conserving condensed-phase sensible enthalpy.

# Edge case
If neither cell carries condensed heat capacity (both essentially pure gas),
the condensed energy balance is degenerate and the volume-weighted average is
returned.
"""
function merge_temperatures_conservatively(
    V_a::Float64,
    state_a::CellState{NC},
    V_b::Float64,
    state_b::CellState{NC},
    material::Material{NC},
) where {NC}
    T_a = state_a.T
    T_b = state_b.T
    V_total = V_a + V_b

    # Condensed-phase volumetric heat capacities (degenerate-case test and
    # Newton initial guess) and total condensed sensible enthalpy.
    ρcₚ_a = 0.0
    ρcₚ_b = 0.0
    H_total = 0.0
    @inbounds for j = 1:NC
        comp = material.components[j]
        is_gas(comp) && continue
        c_j = comp.heat_capacity
        ρcₚ_a += state_a.ξ[j] * Float64(c_j(T_a))
        ρcₚ_b += state_b.ξ[j] * Float64(c_j(T_b))
        H_total +=
            V_a * state_a.ξ[j] * _integrate_cp(c_j, T_REF, T_a) +
            V_b * state_b.ξ[j] * _integrate_cp(c_j, T_REF, T_b)
    end

    if ρcₚ_a < 1e-20 && ρcₚ_b < 1e-20
        # No condensed matter in either cell: the condensed energy balance
        # cannot determine T. Fall back to volume-weighted averaging.
        return (V_a * T_a + V_b * T_b) / V_total
    end

    T_lo, T_hi = minmax(T_a, T_b)
    T_hi - T_lo < 1e-12 && return T_a

    # Merged condensed enthalpy functional and its derivative. Masses are
    # merged exactly (ξ is volume-weighted), so m_j = V_a ξ_a,j + V_b ξ_b,j.
    H_merged = T -> begin
        s = 0.0
        @inbounds for j = 1:NC
            comp = material.components[j]
            is_gas(comp) && continue
            m_j = V_a * state_a.ξ[j] + V_b * state_b.ξ[j]
            s += m_j * _integrate_cp(comp.heat_capacity, T_REF, T)
        end
        s
    end
    Hprime_merged = T -> begin
        s = 0.0
        @inbounds for j = 1:NC
            comp = material.components[j]
            is_gas(comp) && continue
            m_j = V_a * state_a.ξ[j] + V_b * state_b.ξ[j]
            s += m_j * Float64(comp.heat_capacity(T))
        end
        s
    end

    # Bracketed Newton: H is strictly increasing, and the root lies in
    # [T_lo, T_hi] because H_merged(T_lo) ≤ H_total ≤ H_merged(T_hi).
    T = (V_a * ρcₚ_a * T_a + V_b * ρcₚ_b * T_b) / (V_a * ρcₚ_a + V_b * ρcₚ_b)
    T = clamp(T, T_lo, T_hi)
    lo, hi = T_lo, T_hi
    for _ = 1:60
        f = H_merged(T) - H_total
        if f > 0.0
            hi = T
        else
            lo = T
        end
        fp = Hprime_merged(T)
        T_new = fp > 1e-20 ? T - f / fp : 0.5 * (lo + hi)
        if !(lo < T_new < hi)
            T_new = 0.5 * (lo + hi)  # Newton left the bracket: bisect
        end
        if abs(T_new - T) < 1e-9
            return T_new
        end
        T = T_new
    end
    return T
end

"""
	merge_temperatures_conservatively_from_values(V_a::Float64, T_a::Float64, ξ_a::NTuple{NC,Float64},
												   V_b::Float64, T_b::Float64, ξ_b::NTuple{NC,Float64},
												   material::Material{NC}) where NC -> Float64

Alternate interface for energy-conserving temperature merge using raw values instead of CellState.
Useful when working directly with state vectors.
"""
function merge_temperatures_conservatively_from_values(
    V_a::Float64,
    T_a::Float64,
    ξ_a::NTuple{NC,Float64},
    V_b::Float64,
    T_b::Float64,
    ξ_b::NTuple{NC,Float64},
    material::Material{NC},
) where {NC}

    # Create temporary cell states
    state_a = CellState{NC}(T_a, ξ_a, V_a, V_a)
    state_b = CellState{NC}(T_b, ξ_b, V_b, V_b)

    return merge_temperatures_conservatively(V_a, state_a, V_b, state_b, material)
end
