# =============================================================================
# Exact solution: in-depth absorption with conduction (finite slab)
# =============================================================================
#
# Lautenberger dissertation §3.4.2 + Appendix A, re-derived here from scratch
# (formula-provenance rule). Slab 0 ≤ x ≤ δ, x = depth below the irradiated
# front face; volumetric source S(x) = q₀·κ·e^(−κx); convective loss
# k·θ_x(0) = h·θ(0) at the front (θ = T − T₀, T∞ = T₀); insulated back
# θ_x(δ) = 0; θ(x,0) = 0.
#
# Steady part (integrating θ'' = −(q₀κ/k)e^(−κx) with the two BCs):
#
#   θ_ss(x) = q₀[(1−e^(−κδ))/h + (1−e^(−κx))/(kκ)] − (q₀/k)e^(−κδ)·x
#
# (steady front loss h·θ_ss(0) = q₀(1−e^(−κδ)) = total deposited power; the
# transmitted fraction e^(−κδ) never enters the slab.)
#
# Transient: θ = θ_ss − Σₙ cₙ e^(−α λₙ² t) cos(λₙ(δ−x)). The eigenfunctions
# cos(λ(δ−x)) satisfy the back BC identically; the front BC gives
# cot(λδ) = kλ/h, which is uₙ·tan(uₙ) = Bi with u = λδ, Bi = hδ/k — the SAME
# eigencondition as the convective slab, so `biot_eigenvalues` is reused.
# Coefficients cₙ = ∫θ_ss·cos(λₙ(δ−x))dx / (δ/2 + sin(2uₙ)/(4λₙ)) are
# evaluated by composite Simpson (machine-accurate for these smooth,
# BC-compatible integrands); the self-tests verify the series against the IC
# and the PDE residual by finite differences.

"""
    make_indepth_absorption_exact(; δ, κ, k, h, α, q0, T0, n_roots = 60,
                                  quad_intervals = 4096)

Build the exact solution for the in-depth-absorption slab. Returns a
NamedTuple `(T_exact, θ_ss, λs, cs)` where `T_exact(x, t)` is the
temperature at depth `x` (from the irradiated face) and time `t`, and
`θ_ss(x)` the steady rise.
"""
function make_indepth_absorption_exact(;
    δ::Real,
    κ::Real,
    k::Real,
    h::Real,
    α::Real,
    q0::Real,
    T0::Real,
    n_roots::Int = 60,
    quad_intervals::Int = 4096,
)
    Bi = h * δ / k
    us = biot_eigenvalues(Float64(Bi), n_roots)
    λs = us ./ δ

    θ_ss(x) =
        q0 * ((1 - exp(-κ * δ)) / h + (1 - exp(-κ * x)) / (k * κ)) -
        (q0 / k) * exp(-κ * δ) * x

    # Composite-Simpson projection of θ_ss onto each eigenfunction.
    n = quad_intervals  # even
    Δ = δ / n
    xs = [i * Δ for i = 0:n]
    w = [i == 0 || i == n ? 1.0 : (isodd(i) ? 4.0 : 2.0) for i = 0:n]
    θvals = θ_ss.(xs)
    cs = map(eachindex(λs)) do m
        λ = λs[m]
        num = Δ / 3 * sum(w[i+1] * θvals[i+1] * cos(λ * (δ - xs[i+1])) for i = 0:n)
        den = δ / 2 + sin(2 * us[m]) / (4 * λ)
        num / den
    end

    function T_exact(x::Real, t::Real)
        t <= 0 && return float(T0)
        θ = θ_ss(x)
        for m in eachindex(λs)
            term = cs[m] * exp(-α * λs[m]^2 * t) * cos(λs[m] * (δ - x))
            θ -= term
            abs(term) < 1e-16 * max(1.0, abs(θ)) && break
        end
        return T0 + θ
    end

    return (T_exact = T_exact, θ_ss = θ_ss, λs = λs, cs = cs)
end
