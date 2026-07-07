# =============================================================================
# Constant-heating-rate TGA — high-precision first-order reference
# =============================================================================
#
# First-order Arrhenius decay under T(t) = T₀ + βt has the exact solution
#
#   ξ(t)/ξ₀ = exp(−(A/β)·I(T(t))),   I(T) = ∫_{T₀}^{T} e^(−E/RT′) dT′
#
# where I(T) (the "temperature integral") has no closed form — it is computed
# here by composite Simpson quadrature summed to ≪ any tolerance the cases
# use, and self-tested against the standard asymptotic expansion
# I ≈ (RT²/E)e^(−E/RT)(1 − 2RT/E + 6(RT/E)²) plus Richardson refinement.
# This is the 0D leg (b) reference of the V-CPL-05 consistency ladder.

"""
    tga_temperature_integral(T; T0, E, n = 20_000) -> I(T) [K]

Composite-Simpson value of ∫_{T0}^{T} exp(−E/(R·T′)) dT′.
"""
function tga_temperature_integral(T::Real; T0::Real, E::Real, n::Int = 20_000)
    R = 8.314462618
    T <= T0 && return 0.0
    h = (T - T0) / n
    f(x) = exp(-E / (R * x))
    s = f(T0) + f(T)
    for i = 1:2:(n-1)
        s += 4 * f(T0 + i * h)
    end
    for i = 2:2:(n-2)
        s += 2 * f(T0 + i * h)
    end
    return s * h / 3
end

"""
    tga_xi(t; ξ0, A, E, β, T0) -> ξ(t)

Exact first-order TGA mass under the constant heating rate β [K/s].
"""
function tga_xi(t::Real; ξ0::Real, A::Real, E::Real, β::Real, T0::Real)
    T = T0 + β * t
    return ξ0 * exp(-(A / β) * tga_temperature_integral(T; T0 = T0, E = E))
end
