# =============================================================================
# Exact solutions: finite-slab eigenfunction series
# =============================================================================
#
# Coordinate for both solutions: x̃ = x/L ∈ [0, 1] measured from the ADIABATIC
# back face (x̃ = 0), with the active boundary (convection or flux) at x̃ = 1.
# In mesh coordinates that is simply x̃ = z/L (bottom = adiabatic).

"""
    biot_eigenvalues(Bi, n_roots)

First `n_roots` positive roots of the slab eigencondition `λ tan λ = Bi`
(Carslaw & Jaeger §3.11; the symmetric-slab/Heisler problem). Root `n` lies
in `((n−1)π, (n−1)π + π/2)`; bracketed there and bisected on the pole-free
form `g(λ) = λ sin λ − Bi cos λ`.
"""
function biot_eigenvalues(Bi::Float64, n_roots::Int)
    Bi > 0 || throw(ArgumentError("Bi must be positive"))
    g(λ) = λ * sin(λ) - Bi * cos(λ)
    return [
        bisect_root(g, (n - 1) * π + 1e-12, (n - 1) * π + π / 2 - 1e-12) for n = 1:n_roots
    ]
end

"""
    slab_convective_theta(x̃, Fo, λs)

Dimensionless temperature `θ/θᵢ = (T − T∞)/(Tᵢ − T∞)` for a slab initially at
uniform `Tᵢ`, insulated at `x̃ = 0`, convecting to `T∞` at `x̃ = 1` with
`Bi = hL/k` (encoded in the precomputed eigenvalues `λs` from
[`biot_eigenvalues`](@ref)):

    θ/θᵢ = Σₙ Cₙ exp(−λₙ² Fo) cos(λₙ x̃),   Cₙ = 4 sin λₙ / (2λₙ + sin 2λₙ)

Truncation error is `O(exp(−λ_max² Fo))`: with 60 roots the series is
machine-accurate for `Fo ≳ 1.1e-3`.
"""
function slab_convective_theta(x̃::Real, Fo::Real, λs::Vector{Float64})
    Fo <= 0 && return 1.0
    acc = 0.0
    for λ in λs
        term = 4 * sin(λ) / (2 * λ + sin(2 * λ)) * exp(-λ^2 * Fo) * cos(λ * x̃)
        acc += term
        abs(term) < 1e-17 * max(1.0, abs(acc)) && break
    end
    return acc
end

"""
    slab_constflux_phi(x̃, Fo; max_terms = 5000)

Dimensionless temperature rise `φ = (T − T₀)·k/(q″L)` for a slab insulated at
`x̃ = 0` with constant flux `q″` applied at `x̃ = 1` (Carslaw & Jaeger §3.8):

    φ = Fo + (3x̃² − 1)/6 − (2/π²) Σₙ ((−1)ⁿ/n²) exp(−n²π²Fo) cos(nπx̃)

Exact identities used as self-tests: `φ(x̃, 0) = 0` (the Fourier series of
`(3x̃²−1)/6`), front–back plateau `φ(1) − φ(0) → 1/2` as `Fo → ∞`, and
spatial mean `⟨φ⟩ = Fo` (global energy balance).
"""
function slab_constflux_phi(x̃::Real, Fo::Real; max_terms::Int = 5000)
    Fo < 0 && return 0.0
    acc = Fo + (3 * x̃^2 - 1) / 6
    for n = 1:max_terms
        term = -(2 / π^2) * ((-1)^n / n^2) * exp(-n^2 * π^2 * Fo) * cos(n * π * x̃)
        acc += term
        # decay bound: remaining tail < term magnitude envelope
        abs(2 / π^2 / n^2 * exp(-n^2 * π^2 * Fo)) < 1e-17 * max(1.0, abs(acc)) && break
    end
    return acc
end
