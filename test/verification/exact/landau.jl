# =============================================================================
# Landau ablation problem — exact steadily-regressing traveling wave
# =============================================================================
#
# H.G. Landau, "Heat conduction in a melting solid", Q. Appl. Math. 8:81 (1950).
#
# Semi-infinite solid, constant net surface flux q″, material removed at fixed
# ablation temperature T_p with heat of ablation ΔH_p. In the frame attached to
# the regressing surface (ζ = depth below the surface) the steady solution is a
# traveling wave:
#
#   v_ss = q″ / (ρ [c (T_p − T₀) + ΔH_p])          (regression rate)
#   T(ζ) = T₀ + (T_p − T₀) exp(−v_ss ζ / α)         (α = k/ρc)
#
# with steady MLR = ρ v_ss and thermal-wave thickness δ = α / v_ss.
#
# Derivation sketch (verified in selftests.jl): in the wave frame the energy
# equation is α T″ + v_ss T′ = 0 with T(0) = T_p, T(∞) = T₀ ⇒ the exponential
# profile; the surface balance q″ = k|dT/dζ|₀ + ρ v_ss ΔH_p then fixes v_ss,
# since k (T_p − T₀) v_ss/α = ρ c (T_p − T₀) v_ss.

"""
    landau_vss(; q, ρ, c, ΔH, T0, Tp) -> v_ss [m/s]

Exact steady regression rate of the Landau ablation problem.
`q` is the net absorbed surface flux [W/m²], `ΔH` the heat of ablation [J/kg]
(positive = endothermic), `Tp` the ablation temperature [K].
"""
landau_vss(; q::Real, ρ::Real, c::Real, ΔH::Real, T0::Real, Tp::Real) =
    q / (ρ * (c * (Tp - T0) + ΔH))

"""
    landau_T(ζ; T0, Tp, v, α) -> T [K]

Exact steady in-depth temperature at depth `ζ ≥ 0` below the regressing
surface, for regression rate `v` and thermal diffusivity `α`.
"""
landau_T(ζ::Real; T0::Real, Tp::Real, v::Real, α::Real) = T0 + (Tp - T0) * exp(-v * ζ / α)
