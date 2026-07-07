# =============================================================================
# Frank-Kamenetskii critical slab — exact subcritical profile + criticality
# =============================================================================
#
# D.A. Frank-Kamenetskii, "Diffusion and Heat Exchange in Chemical Kinetics".
# Zero-order exothermic Arrhenius source in a slab with isothermal walls,
# high-activation-energy (FK) limit. With x̃ = x/L (L = HALF-width),
# θ = (E/RT₀²)(T − T₀), δ = (E/RT₀²)·S₀L²/k (S₀ = heat release rate at T₀):
#
#   θ″ + δ e^θ = 0,  θ(±1) = 0   ⇒   θ(x̃) = θ_m − 2 ln cosh(x̃ √(δ e^(θ_m)/2))
#
# where θ_m solves θ_m = 2 ln cosh(√(δ e^(θ_m)/2)) — two roots below the
# critical δ_c (lower = stable), merging at δ_c. Parametrically, with
# u = √(δ e^(θ_m)/2): δ(u) = 2u²/cosh²(u), maximized where u·tanh(u) = 1,
# giving δ_c = 0.8785 and θ_m,c = 1.1868 (textbook: 0.878, ≈1.19).
#
# The library also provides a full-Arrhenius steady reference (no FK
# truncation): k T″ + S₀ e^(E/R (1/T₀ − 1/T)) = 0 solved by RK4 shooting from
# the symmetric center with bisection on T_center — the high-precision pin
# the case compares against at finite E/RT₀ (dependency-free by design; the
# test manifest is pin-fragile, so no BoundaryValueDiffEq).

"""
    fk_ustar() -> u*

Root of `u·tanh(u) = 1` (parametric location of the FK slab turning point).
"""
fk_ustar() = bisect_root(u -> u * tanh(u) - 1.0, 0.5, 2.0)

"Critical FK parameter for the slab: δ_c = 2u*²/cosh²(u*) ≈ 0.8785."
function fk_delta_crit()
    u = fk_ustar()
    return 2.0 * u^2 / cosh(u)^2
end

"Critical center rise θ_m,c = 2 ln cosh(u*) ≈ 1.1868."
fk_theta_max_crit() = 2.0 * log(cosh(fk_ustar()))

"""
    fk_theta_max(δ) -> θ_m

Lower-branch (stable) root of `θ_m = 2 ln cosh(√(δ e^(θ_m)/2))` for
subcritical `δ < δ_c`; bracketed on [0, θ_m,c] where the residual changes
sign exactly once.
"""
function fk_theta_max(δ::Real)
    δc = fk_delta_crit()
    δ < δc || throw(ArgumentError("δ = $δ ≥ δ_c = $δc: no subcritical steady state"))
    g(θm) = θm - 2.0 * log(cosh(sqrt(δ * exp(θm) / 2.0)))
    return bisect_root(g, 0.0, fk_theta_max_crit())
end

"""
    fk_theta_profile(x̃; θm, δ) -> θ

Exact subcritical steady profile, `x̃ ∈ [−1, 1]` from the slab center.
"""
fk_theta_profile(x̃::Real; θm::Real, δ::Real) =
    θm - 2.0 * log(cosh(x̃ * sqrt(δ * exp(θm) / 2.0)))

"""
    fk_steady_shoot(; S0, E, T0, k, L; n_steps = 4000, θ_upper) -> (; T_center, T_at)

Full-Arrhenius steady reference: solve `k T″ + S₀ exp(E/R (1/T₀ − 1/T)) = 0`
on the half-slab x ∈ [0, L] with T′(0) = 0 (symmetry) and T(L) = T₀, by RK4
integration from the center and bisection on T_center to the STABLE (lower)
root. `θ_upper` must lie strictly between the two steady branches in θ units
(default: midpoint of the FK roots — adequate for E/RT₀ ≳ 20). Returns the
center temperature and an interpolating function `T_at(x)` on the RK4 grid.
"""
function fk_steady_shoot(;
    S0::Real,
    E::Real,
    T0::Real,
    k::Real,
    L::Real,
    n_steps::Int = 4000,
    θ_upper::Union{Nothing,Real} = nothing,
)
    R = 8.314462618
    εT = R * T0^2 / E   # 1 θ-unit in Kelvin
    S(T) = S0 * exp(E / R * (1.0 / T0 - 1.0 / T))

    # Integrate T from the center with T(0) = Tc, T′(0) = 0; return T(L)
    # and the grid values.
    function wall_T(Tc; keep::Bool = false)
        h = L / n_steps
        T, W = float(Tc), 0.0
        Ts = keep ? Vector{Float64}(undef, n_steps + 1) : Float64[]
        keep && (Ts[1] = T)
        f(T) = -S(T) / k
        for i = 1:n_steps
            k1T, k1W = W, f(T)
            k2T, k2W = W + h / 2 * k1W, f(T + h / 2 * k1T)
            k3T, k3W = W + h / 2 * k2W, f(T + h / 2 * k2T)
            k4T, k4W = W + h * k3W, f(T + h * k3T)
            T += h / 6 * (k1T + 2k2T + 2k3T + k4T)
            W += h / 6 * (k1W + 2k2W + 2k3W + k4W)
            keep && (Ts[i+1] = T)
        end
        return T, Ts
    end

    # Default bracket top: 60% of the way from the FK lower root toward the
    # critical θ_m,c — strictly between the stable and unstable branches for
    # the full-Arrhenius equation as long as its roots sit within O(RT₀/E)
    # of the FK ones (fine for E/RT₀ ≳ 20).
    δ = S0 * L^2 / (k * εT)
    θ_up = if θ_upper === nothing
        θ_lo = fk_theta_max(min(δ, 0.999 * fk_delta_crit()))
        θ_lo + 0.6 * (fk_theta_max_crit() - θ_lo)
    else
        float(θ_upper)
    end
    # Residual T(L) − T0 is negative at Tc = T0 (pure cooling) and positive
    # between the branches; bisect to the stable root.
    g(Tc) = wall_T(Tc)[1] - T0
    Tc = bisect_root(g, T0 + 1e-9, T0 + θ_up * εT)
    _, Ts = wall_T(Tc; keep = true)
    h = L / n_steps
    function T_at(x::Real)
        0.0 <= x <= L || throw(ArgumentError("x outside [0, L]"))
        j = min(floor(Int, x / h), n_steps - 1)
        w = (x - j * h) / h
        return (1 - w) * Ts[j+1] + w * Ts[j+2]
    end
    return (; T_center = Tc, T_at = T_at)
end
