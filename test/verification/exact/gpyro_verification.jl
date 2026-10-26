# Problem provenance: Abdelhafiz et al., Gpyro Verification and Validation
# Manual, version 1.0 (15 January 2026), Gpyro V2.0, §§3.6, 4.3–4.4.
# Pin the guide revision so published case citations remain reproducible.
const GPYRO_VERIFICATION_GUIDE = "https://github.com/Gedeben/GPYRO-V2/blob/bff0126c68f2a35c084ef04490ae327185f281a3/documentation/gpyro_verification_manual.pdf"

"""
    robin_diffusion_uptake(t; D, b)

Integral of the semi-infinite Robin diffusion response over depth [m].
For Gpyro §4.4, `b = h_m/(φ ρ_g D)` and multiplying by
`φ ρ_g (Y∞ − Y₀)` gives the species inventory gained per unit area.
Integrating the surface flux `h_m (Y∞ − Y₀) erfcx(b√(Dt))` gives
`2√(Dt/π) + (erfcx(b√(Dt)) − 1)/b`.
"""
function robin_diffusion_uptake(t; D, b)
    t <= 0 && return 0.0
    return 2sqrt(D*t/π) + (erfcx(b*sqrt(D*t))-1)/b
end

"""
    darcy_imposed_flux_pressure(x; Pout, flux, μ, T, κ, M)

Steady isothermal pure-gas pressure at distance `x` inward from the outlet,
with positive imposed inflow magnitude `flux` [kg/(m² s)]. Gpyro §4.3,
Eq. (4.5), with `α = μRT/M`: `P² = Pout² + 2μRT flux x/(κM)`.
"""
darcy_imposed_flux_pressure(x; Pout, flux, μ, T, κ, M) =
    sqrt(Pout^2 + 2μ*Pyrolysis.Materials.R_GAS*T*flux*x/(κ*M))

"""
    adiabatic_reaction_temperature(t; T0, rate, heat, ca, cb)

Independent insulated-cell adaptation of Gpyro §3.6. A first-order,
temperature-independent reaction leaves fraction `a = exp(-rate*t)` of
the reactant; `cb` is the yield-weighted heat capacity of retained products.
The temperature equation is `[ca*a + cb*(1-a)] T' = -heat*rate*a`, hence
`T = T0 + heat/(ca-cb) log([ca*a + cb*(1-a)]/ca)`.
Heat is positive for an endothermic reaction. Equal heat capacities give
`T = T0 - heat*(1-a)/ca`. This follows the package's local-temperature
reaction-heat convention, including sensible-enthalpy formation in H form.
"""
function adiabatic_reaction_temperature(t; T0, rate, heat, ca, cb)
    converted = -expm1(-rate*t)
    ca == cb && return T0 - heat*converted/ca
    return T0 + heat/(ca-cb)*log1p((cb-ca)*converted/ca)
end
