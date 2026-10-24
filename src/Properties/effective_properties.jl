"""
Effective property computations.

PERFORMANCE NOTE: helpers that iterate over `material.components` are
`@generated` so the loop is unrolled at compile time per concrete `NC`. This
avoids the allocations that bare runtime indexing into a heterogeneous
`Tuple{...}` would otherwise incur.
"""

# ============================================================================
# Floors of the gas properties
# ============================================================================

# Robustness floor for the Chapman-Enskog diffusivity formula (T^{3/2} and a
# collision integral in T*). Property laws clamp their own argument to the
# property bracket; this formula is not a property law, so it carries its own
# floor. Newton trial states may present nonphysical temperatures, and the
# floor keeps the residual evaluable so the nonlinear solver can reject the
# step instead of the solve dying on a DomainError. 100 K is far below any
# physical pyrolysis state.
const T_GAS_PROPERTY_FLOOR = 100.0  # K

# Positive floor on the gas viscosity. μ appears in the model ONLY as the
# denominator of the Darcy mobility κ/μ, so a non-positive μ does not merely
# degrade the answer — it reverses the gas velocity or divides by zero. The
# built-in Sutherland form cannot go non-positive above `T_GAS_PROPERTY_FLOOR`,
# but a material may carry an arbitrary correlation (the Cantera blend fit is a
# quartic valid over 280–1200 K), and Newton trial states reach temperatures no
# fit was made for. `Material`'s constructor rejects a curve that is already
# unusable across the pyrolysis range; this is the runtime backstop for
# excursions outside it. 1e-12 Pa·s is seven orders below any real gas, so it
# never binds on a legitimate state; the state mapping applies it as a smooth
# floor.
const MU_FLOOR = 1e-12  # Pa·s

# ============================================================================
# Chapman-Enskog Diffusivity
# ============================================================================

# Physical constant for Chapman-Enskog formula (Bird, Stewart, Lightfoot)
# D_AB = 1.8583×10⁻³ × T^(3/2) × √(1/M_A + 1/M_B) / (P × σ_AB² × Ω_D)
# where T[K], M[g/mol], P[atm], σ[Å] → D[cm²/s]
const CHAPMAN_ENSKOG_PREFACTOR = 1.8583e-3  # cm²/s × atm × Å² / (K^1.5 × (g/mol)^0.5)

# Default Lennard-Jones parameters for air (background gas)
const SIGMA_AIR = 3.617   # Å - collision diameter for air
const EPSILON_K_AIR = 97.0  # K - well depth / k_B for air
const M_AIR = 0.029       # kg/mol - molar mass of air

"""
	collision_integral_diffusion(T_star::Float64) -> Float64

Compute the diffusion collision integral Ω_D using the Neufeld correlation.

This is an empirical fit to tabulated collision integrals for the 
Lennard-Jones (12-6) potential.

Ω_D(T*) = A/(T*)^B + C/exp(D·T*) + E/exp(F·T*) + G/exp(H·T*)

Valid for 0.3 ≤ T* ≤ 400 with < 0.1% error.

# Reference
Neufeld, Janzen, Aziz, J. Chem. Phys. 57, 1100 (1972)

# Arguments
- `T_star`: Reduced temperature T* = k_B T / ε_AB

# Returns
- `Ω_D`: Dimensionless diffusion collision integral
"""
@inline function collision_integral_diffusion(T_star::Real)
    # Neufeld et al. correlation coefficients
    A = 1.06036
    B = 0.15610
    C = 0.19300
    D = 0.47635
    E = 1.03587
    F = 1.52996
    G = 1.76474
    H = 3.89411

    return A / T_star^B + C / exp(D * T_star) + E / exp(F * T_star) + G / exp(H * T_star)
end

"""
	chapman_enskog_diffusivity(T, M; P=P_REF, σ=nothing, ε_k=nothing) -> Float64

Compute binary diffusion coefficient using Chapman-Enskog kinetic theory.

The Chapman-Enskog expression for binary diffusion of species A in background gas B:

D_AB = (5.954×10⁻⁴) × T^(3/2) × √(1/M_A + 1/M_B) / (P × σ_AB² × Ω_D)   [m²/s]

where:
- T: Temperature [K]
- M_A, M_B: Molar masses [kg/mol]
- P: Pressure [Pa]
- σ_AB: Average collision diameter [Å]
- Ω_D: Diffusion collision integral (function of T*)

(The Bird–Stewart–Lightfoot prefactor 1.8583×10⁻³ for cm²/s, g/mol and atm,
converted to these units.)

# Arguments
- `T`: Temperature [K]
- `M`: Molar mass of diffusing species [kg/mol]
- `P`: Total pressure [Pa] (default: 101325 Pa)
- `σ`: Collision diameter [Å] (default: estimated from Tc or molar mass)
- `ε_k`: Well depth / k_B [K] (default: estimated)

# Returns
- `D`: Binary diffusion coefficient [m²/s]

# Default Parameters
For pyrolysis gases, if σ and ε_k are not provided, uses:
- σ estimated from Tc-based correlation or default to 4 Å
- ε_k estimated as 0.77 × Tc or default to 300 K
- Background gas is air (σ_air = 3.617 Å, ε_k_air = 97 K, M_air = 0.029 kg/mol)

# Examples
```julia
# Water vapor at 500 K
D_H2O = chapman_enskog_diffusivity(500.0, 0.018)  # ≈ 5.41e-5 m²/s

# CO₂ at 600 K
D_CO2 = chapman_enskog_diffusivity(600.0, 0.044)  # ≈ 5.44e-5 m²/s

# With explicit Lennard-Jones parameters
D_CH4 = chapman_enskog_diffusivity(500.0, 0.016, σ=3.758, ε_k=148.6)
```

# Physical Behavior
- D ∝ T^1.5 to T^1.75 (depending on collision integral)
- D ∝ 1/P (inversely proportional to pressure)
- D ∝ √(1/M) (lighter species diffuse faster)

# Notes
- Valid for dilute gas limit (low density)
- Assumes ideal gas behavior
- Uses combining rules: σ_AB = (σ_A + σ_B)/2, ε_AB = √(ε_A × ε_B)
"""
@inline function chapman_enskog_diffusivity(
    T::Real,
    M::Real;
    P::Float64 = P_REF,
    σ::Union{Float64,Nothing} = nothing,
    ε_k::Union{Float64,Nothing} = nothing,
)

    # Estimate Lennard-Jones parameters if not provided
    # Using correlations for organic vapors
    σ_species = if σ !== nothing
        σ
    else
        # Rough estimate: σ ≈ 1.18 × (Vc)^(1/3) where Vc ~ M/ρ_c
        # For organic gases, σ typically 3-6 Å, estimate from M
        3.5 + 0.5 * log10(max(M / 0.016, 1.0))  # Å, anchored to CH4 ≈ 3.8 Å
    end

    ε_k_species = if ε_k !== nothing
        ε_k
    else
        # Rough estimate: ε/k ≈ 0.77 × Tc
        # For organics, ε/k typically 100-600 K
        150.0 + 100.0 * (M / 0.044)  # K, scaled roughly to molecular weight
    end

    # Combining rules for binary interaction parameters
    σ_AB = 0.5 * (σ_species + SIGMA_AIR)  # Å
    ε_k_AB = sqrt(ε_k_species * EPSILON_K_AIR)  # K

    # Clamp nonphysical Newton trial temperatures (see T_GAS_PROPERTY_FLOOR):
    # `T^1.5` below would throw a DomainError on a negative iterate.
    T_eff = max(T, T_GAS_PROPERTY_FLOOR)

    # Reduced temperature
    T_star = T_eff / ε_k_AB

    # Collision integral
    Ω_D = collision_integral_diffusion(T_star)

    # Chapman-Enskog diffusivity
    # Standard formula (Bird, Stewart, Lightfoot):
    # D_AB = 1.8583×10⁻³ × T^(3/2) × √(1/M_A + 1/M_B) / (P × σ_AB² × Ω_D)
    # where T in K, M in g/mol, P in atm, σ in Å → D in cm²/s
    #
    # Convert units:
    # - M: kg/mol → g/mol (multiply by 1000)
    # - P: Pa → atm (divide by 101325)
    # - D: cm²/s → m²/s (multiply by 1e-4)

    M_gmol = M * 1000.0        # kg/mol → g/mol
    M_air_gmol = M_AIR * 1000.0  # kg/mol → g/mol
    P_atm = P / 101325.0      # Pa → atm

    # Reduced mass factor with M in g/mol
    reduced_mass_factor = sqrt(1.0 / M_gmol + 1.0 / M_air_gmol)

    # D in cm²/s
    D_cm2s =
        CHAPMAN_ENSKOG_PREFACTOR * T_eff^1.5 * reduced_mass_factor / (P_atm * σ_AB^2 * Ω_D)

    # Convert to m²/s
    D = D_cm2s * 1e-4

    return D
end

# ============================================================================
# Conductivity helpers
# ============================================================================

# With `skip_absent`, a component whose volume fraction is exactly zero (a
# dual number's partials included) adds an exact zero to both sums for every
# finite conductivity, so its law is not evaluated.
@generated function _conductivity_sums(
    v,
    T,
    comps::Tuple,
    ::Val{NC},
    ::Val{skip_absent},
) where {NC,skip_absent}
    expr = quote
        W = _prop_type(v, T)
        k_parallel = zero(W)
        inv_k_series = zero(W)
    end
    for i = 1:NC
        term = quote
            k = comps[$i].thermal_conductivity(T)
            k_parallel += v[$i] * k
            if v[$i] > 1e-20 && k > 1e-20
                inv_k_series += v[$i] / k
            end
        end
        push!(expr.args, skip_absent ? :(iszero(v[$i]) || $term) : term)
    end
    push!(expr.args, :(return k_parallel, inv_k_series))
    return expr
end

# ============================================================================
# Bruggeman Effective Medium Theory
# ============================================================================

"""
	bruggeman_effective_property(v::NTuple{NC,Float64}, k::NTuple{NC,Float64};
								 tol=1e-10, maxiter=100) -> Float64

Solve for effective property using Bruggeman effective medium theory.

The Bruggeman (symmetric) equation for spherical inclusions is:
	Σⱼ vⱼ (kⱼ - k_eff) / (kⱼ + 2k_eff) = 0

This is solved iteratively using Newton-Raphson.

# Arguments
- `v`: Volume fractions of components (should sum to 1)
- `k`: Property values of components
- `tol`: Convergence tolerance (default: 1e-10)
- `maxiter`: Maximum iterations (default: 100)

# Returns
- `k_eff`: Effective property value

# Physical Interpretation
Bruggeman EMT assumes a self-consistent medium where each phase is surrounded
by the effective medium. It's particularly appropriate for:
- Random heterogeneous mixtures
- Self-similar microstructures (char formation)
- Media near percolation threshold

Unlike parallel/series bounds, Bruggeman can capture:
- Percolation transitions (property drop when continuous phase breaks)
- Non-linear property changes with composition

An insulating phase (k ≤ 1e-20) contributes its exact k → 0 limit −v/2 to the
self-consistency equation. That term is what produces the percolation
threshold: a conducting phase dispersed in an insulator yields k_eff = 0 until
its volume fraction exceeds 1/3 (e.g. v = (0.5, 0.5), k = (1, 0) → k_eff = 0.25).

# Notes
- For two phases, there's an analytical solution (not used here for generality)
- Initial guess is volume-weighted average (parallel mixing)
- Converges quadratically for well-posed problems

# Example
```julia
v = (0.7, 0.3)  # 70% phase 1, 30% phase 2
k = (1.0, 10.0)  # k₁ = 1, k₂ = 10
k_eff = bruggeman_effective_property(v, k)  # ≈ 2.261 (between series = 1.370 and parallel = 3.7)
```
"""
@inline function bruggeman_effective_property(
    v::NTuple{NC,<:Real},
    k::NTuple{NC,<:Real};
    tol::Float64 = 1e-10,
    maxiter::Int = 100,
) where {NC}

    # Working scalar type: keeps every early return and the Newton state in
    # one concrete type when `v`/`k` carry Duals (type-stable, no boxing).
    W = promote_type(_prop_type(v), _prop_type(k))

    # Handle degenerate cases
    total_v = sum(v)
    if total_v < 1e-20
        return zero(W)
    end

    # Count non-zero volume components and compute the parallel bound, the
    # initial guess.
    k_parallel = zero(W)
    n_active = 0
    for i = 1:NC
        if v[i] > 1e-20
            k_parallel += v[i] * k[i]
            n_active += 1
        end
    end
    k_parallel /= total_v

    k_eff = k_parallel

    # Single component - return its property
    if n_active == 1
        return k_parallel
    end

    # Newton-Raphson iteration to solve Σ vⱼ(kⱼ - k_eff)/(kⱼ + 2k_eff) = 0
    for _ = 1:maxiter
        f = zero(W)      # Function value
        df = zero(W)     # Derivative w.r.t. k_eff

        @inbounds for i = 1:NC
            if v[i] > 1e-20
                if k[i] > 1e-20
                    denom = k[i] + 2.0 * k_eff
                    f += v[i] * (k[i] - k_eff) / denom
                    # d/dk_eff of [v(k - k_eff)/(k + 2k_eff)]
                    # = v * [-(k + 2k_eff) - (k - k_eff)*2] / (k + 2k_eff)²
                    # = v * [-k - 2k_eff - 2k + 2k_eff] / (k + 2k_eff)²
                    # = v * [-3k] / (k + 2k_eff)²
                    df += v[i] * (-3.0 * k[i]) / (denom * denom)
                else
                    # Insulating phase: v(k - k_eff)/(k + 2k_eff) → -v/2 as
                    # k → 0 (df contribution → 0). This term is what produces
                    # the v = 1/3 percolation threshold; skipping it makes the
                    # conducting phase percolate at any volume fraction.
                    f -= 0.5 * v[i]
                end
            end
        end

        # Newton update with safeguard
        if abs(df) < 1e-20
            break  # Derivative too small, can't continue
        end

        dk = -f / df

        # Limit step size to prevent negative k_eff
        if k_eff + dk < 0.1 * k_eff
            dk = -0.9 * k_eff  # Don't reduce by more than 90%
        end

        k_eff += dk

        # Below the percolation threshold (conducting fraction ≤ 1/3) the root
        # is k_eff = 0 and f < 0 for every k_eff > 0, so the step-capped Newton
        # decays k_eff geometrically without ever satisfying |f| < tol. Stop
        # once it has collapsed instead of burning maxiter.
        if k_eff < 1e-30
            return zero(W)
        end

        # Convergence check
        if abs(f) < tol
            break
        end
    end

    return max(k_eff, zero(W))
end

# Helper to extract conductivities as tuple for Bruggeman
@generated function _conductivity_tuple(T, comps::Tuple, ::Val{NC}) where {NC}
    exprs = [:(comps[$i].thermal_conductivity(T)) for i = 1:NC]
    return :(($(exprs...),))
end

"""
	effective_conductivity(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute effective thermal conductivity [W/(m·K)].

Uses the mixing rule specified in the material:
- PARALLEL: k = Σⱼ vⱼ kⱼ
- SERIES: k = 1 / Σⱼ (vⱼ / kⱼ)
- WEIGHTED: k = β·k_parallel + (1-β)·k_series
- BRUGGEMAN: Σⱼ vⱼ(kⱼ - k_eff)/(kⱼ + 2k_eff) = 0

Note: Uses `mixing_volume_fractions` — the swelling-weighted (matrix) volume
fractions, excluding γ ≈ 0 components. This matches ThermaKin behavior: "when
γ = 0, gaseous components do not contribute to the thermal conductivity"
(ThermaKin2Ds Reference, Equations 8-10). For a fully depleted cell (no
matrix volume left) it falls back to the actual phase occupancy so k_eff
degrades to the gas-mixture conductivity instead of 0.

With `skip_absent = Val(true)` the parallel, series and weighted rules do not
evaluate the conductivity of a component whose volume fraction is exactly
zero; the result is the same for every finite conductivity.
"""
@inline function effective_conductivity(
    ξ::NTuple{NC,<:Real},
    T::Real,
    material;
    skip_absent::Val = Val(false),
) where {NC}
    v = mixing_volume_fractions(ξ, T, material)
    return _conductivity_from_v(v, T, material; skip_absent)
end

# The effective conductivity from precomputed mixing volume fractions. The
# tuple type is loose: AD-derived volume fractions can be a mixed-type tuple
# (some Float64 zeros, others Dual).
@inline function _conductivity_from_v(
    v::Tuple,
    T::Real,
    material::Material{NC};
    skip_absent::Val = Val(false),
) where {NC}
    W = _prop_type(v, T)
    mixing = material.mixing.k.rule
    if mixing == BRUGGEMAN
        k_tuple = _conductivity_tuple(T, material.components, Val(NC))
        return convert(W, bruggeman_effective_property(v, k_tuple))
    end

    k_parallel, inv_k_series =
        _conductivity_sums(v, T, material.components, Val(NC), skip_absent)
    k_series = inv_k_series > 1e-20 ? convert(W, 1.0 / inv_k_series) : zero(W)

    if mixing == PARALLEL
        return k_parallel
    elseif mixing == SERIES
        return k_series
    else  # WEIGHTED
        β = material.mixing.k.β
        return β * k_parallel + (1.0 - β) * k_series
    end
end

# ============================================================================
# Emissivity and surface-absorptivity helpers
# ============================================================================
#
# Both sums share the weights `v` (swelling-weighted matrix volume fractions)
# and the cell temperature. `get_absorptivity` returns `get_emissivity(comp, T)`
# for a component without an explicit surface absorptivity, so with none set
# the two sums are the same floating-point operations in the same order and
# `a_eff == ε_eff` bit-for-bit — the Kirchhoff default.

@generated function _sum_emissivity(v, T, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        result = zero(_prop_type(v, T))
    end
    for i = 1:NC
        push!(expr.args, :(result += v[$i] * get_emissivity(comps[$i], T)))
    end
    push!(expr.args, :(return result))
    return expr
end

@generated function _sum_absorptivity(v, T, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        result = zero(_prop_type(v, T))
    end
    for i = 1:NC
        push!(expr.args, :(result += v[$i] * get_absorptivity(comps[$i], T)))
    end
    push!(expr.args, :(return result))
    return expr
end

# ============================================================================
# Absorption helpers
# ============================================================================

@generated function _sum_absorption(ξ, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        result = zero(_prop_type(ξ))
    end
    for i = 1:NC
        push!(expr.args, :(result += ξ[$i] * comps[$i].absorption_coeff))
    end
    push!(expr.args, :(return result))
    return expr
end

"""
	mixture_emissivity(ξ::NTuple{NC,Real}, T::Real, material) -> Real
	mixture_absorptivity(ξ::NTuple{NC,Real}, T::Real, material) -> Real

The surface emissivity `Σⱼ vⱼ εⱼ(T)` and absorptivity for incident radiation
`Σⱼ vⱼ aⱼ`, `aⱼ` the component's `surface_absorptivity` or `εⱼ(T)` when unset
(Kirchhoff's law), weighted by the mixing volume fractions
(`mixing_volume_fractions`: the matrix at the surface, not the gas in its
pores). The weights are a convex combination and every component value lies in
`[0, 1]`, so the sums stay in `[0, 1]` to round-off, with no clamp and so no
kink in their derivatives. Each component's emissivity is a property function
of `T`, so a `ForwardDiff.Dual` `T` carries `dε/dT`.
"""
@inline function mixture_emissivity(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    v = mixing_volume_fractions(ξ, T, material)
    return _sum_emissivity(v, T, material.components, Val(NC))
end

@inline function mixture_absorptivity(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    v = mixing_volume_fractions(ξ, T, material)
    return _sum_absorptivity(v, T, material.components, Val(NC))
end

"""
	effective_absorption(ξ::NTuple{NC,Float64}, material::Material) -> Float64

Compute effective absorption coefficient [m⁻¹].

α_eff = Σⱼ ξⱼ αⱼ

where αⱼ is the mass absorption coefficient [m²/kg].
"""
@inline function effective_absorption(ξ::NTuple{NC,<:Real}, material) where {NC}
    return _sum_absorption(ξ, material.components, Val(NC))
end

# ============================================================================
# Gas transfer helpers
# ============================================================================

@generated function _gas_transfer_sums(v, T, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        W = _prop_type(v, T)
        λ_parallel = zero(W)
        inv_λ_series = zero(W)
    end
    for i = 1:NC
        push!(expr.args, quote
            λ = comps[$i].gas_transfer_coeff(T)
            λ_parallel += v[$i] * λ
            if v[$i] > 1e-20 && !isinf(inv_λ_series)
                if λ > 1e-20
                    inv_λ_series += v[$i] / λ
                else
                    inv_λ_series = oftype(inv_λ_series, Inf)
                end
            end
        end)
    end
    push!(expr.args, :(return λ_parallel, inv_λ_series))
    return expr
end

"""
	effective_gas_transfer(ξ::NTuple{NC,Float64}, T::Float64, material::Material) -> Float64

Compute effective gas transfer coefficient [m²/s].

Uses the same volume-fraction mixing *rules* as conductivity, EXCEPT BRUGGEMAN:
- PARALLEL: λ = Σⱼ vⱼ λⱼ
- SERIES: λ = 1 / Σⱼ (vⱼ / λⱼ)
- WEIGHTED: λ = β·λ_parallel + (1-β)·λ_series

but a different volume-fraction *basis*: the all-phase occupancy
[`volume_fractions`](@ref) (gas included), not the matrix basis
[`solid_volume_fractions`](@ref) used for k and ε — gas diffusion happens
through the pore space, which the matrix basis excludes.

BRUGGEMAN is NOT supported for gas transfer: the EMT closure has no physical
justification for diffusivities, so this function throws an `ArgumentError` if
`material.mixing.λ.rule == BRUGGEMAN`. Use PARALLEL (default), SERIES, or WEIGHTED.
"""
@inline function effective_gas_transfer(ξ::NTuple{NC,<:Real}, T::Real, material) where {NC}
    v = volume_fractions(ξ, T, material)
    return _gas_transfer_from_v(v, T, material)
end

# Cache-aware variant.
@inline function _gas_transfer_from_v(v::Tuple, T::Real, material::Material{NC}) where {NC}
    mixing = material.mixing.λ.rule
    if mixing == BRUGGEMAN
        throw(
            ArgumentError(
                "BRUGGEMAN mixing is not supported for gas-transfer (λ); the EMT closure has no physical justification for diffusivities. Use PARALLEL (default), SERIES, or WEIGHTED.",
            ),
        )
    end

    W = _prop_type(v, T)
    λ_parallel, inv_λ_series = _gas_transfer_sums(v, T, material.components, Val(NC))
    λ_series = inv_λ_series > 1e-20 ? convert(W, 1.0 / inv_λ_series) : zero(W)

    if mixing == PARALLEL
        return λ_parallel
    elseif mixing == SERIES
        return λ_series
    else  # WEIGHTED
        β = material.mixing.λ.β
        return β * λ_parallel + (1.0 - β) * λ_series
    end
end

# ============================================================================
# Effective Permeability (Component-Based)
# ============================================================================

# Helper: compute parallel and series sums for permeability on the
# renormalized solid basis (`solid_volume_fractions`). An impermeable phase
# (κ ≈ 0) with finite volume fraction makes the Wiener series (lower) bound
# exactly 0 — signalled by `inv_κ_series = Inf`, matching `_gas_transfer_sums`.
@generated function _permeability_sums(v, T, comps::Tuple, ::Val{NC}) where {NC}
    expr = quote
        W = _prop_type(v, T)
        κ_parallel = zero(W)
        inv_κ_series = zero(W)
    end
    for i = 1:NC
        push!(expr.args, quote
            κ = get_permeability(comps[$i], T)
            κ_parallel += v[$i] * κ
            if v[$i] > 1e-20 && !isinf(inv_κ_series)
                if κ > 1e-30
                    inv_κ_series += v[$i] / κ
                else
                    inv_κ_series = oftype(inv_κ_series, Inf)
                end
            end
        end)
    end
    push!(expr.args, :(return κ_parallel, inv_κ_series))
    return expr
end

"""
	effective_permeability(ξ::NTuple{NC,Float64}, T::Float64, material, φ) -> Float64

Compute effective permeability [m²] using component mixing rules.

Only solid and liquid components contribute to permeability (they form the matrix).
Gases flow THROUGH the matrix and do not contribute (has_permeability returns false
for gases, and get_permeability returns 0.0).

Mixing operates on the renormalized solid basis (`solid_volume_fractions`),
exactly like conductivity — gas volume is excluded and the surviving weights
sum to 1, so the Wiener bounds are the matrix's true bounds:
- PARALLEL: κ = Σⱼ vⱼ κⱼ (upper bound)
- SERIES: κ = 1 / Σⱼ (vⱼ / κⱼ) (lower bound; **exactly 0** when any phase
  with vⱼ > 0 is impermeable — gas cannot cross a series-stacked κ = 0 layer)
- WEIGHTED: κ = β·κ_parallel + (1-β)·κ_series
- CARMAN_KOZENY: porosity-based closure (see `_carman_kozeny_permeability`) on
  the porosity `φ` the caller gives, so a caller that floors the pore fraction
  `1 − skeletal_volume_fraction` smoothly gets a permeability smooth in the
  state; the per-component mixes do not read it

# Arguments
- `ξ`: Mass concentrations [kg/m³]
- `T`: Temperature [K]
- `material`: Material with components and mixing rules
- `φ`: Porosity, read by the Carman–Kozeny rule only

# Returns
Effective permeability [m²]. Returns 0.0 if no components have permeability.

# Physical Interpretation
Permeability determines resistance to pressure-driven (Darcy) gas flow:
	v_gas = -(κ/μ) ∇P

Typical values:
- Virgin wood: ~10⁻¹⁵ m² (tight cellular structure)
- Char: ~10⁻¹² m² (porous with open channels)
- Impermeable: 0.0 m²

As pyrolysis progresses and virgin material converts to char, the effective
permeability increases, allowing easier gas escape.
"""
@inline function effective_permeability(
    ξ::NTuple{NC,<:Real},
    T::Real,
    material,
    φ::Real,
) where {NC}
    material.mixing.κ.rule == CARMAN_KOZENY &&
        return _carman_kozeny_permeability(φ, material.mixing.κ.ℓ)
    v = solid_volume_fractions(ξ, T, material)
    return _permeability_from_v(v, T, material)
end

"""
	_carman_kozeny_permeability(φ, ℓ) -> κ [m²]

Kozeny–Carman permeability from porosity `φ` and characteristic particle/pore
diameter `ℓ` [m]:

	κ = min( ℓ² · φ³ / (180 · (1 - φ)²),  ℓ² )

The constant 180 is the Kozeny constant for a packed bed of spheres. As φ → 0
(fully solid) κ → 0 (impermeable). As φ → 1 the packed-bed closure loses
validity and diverges; κ is **capped at the physical free-flow scale ℓ²**
(an open channel of width ℓ has κ ~ ℓ²/32 — ℓ² is a generous upper bound).
The cap binds for φ ≳ 0.93 (where 180(1−φ)² < φ³). The 1e-30 denominator
guard only keeps φ = 1.0 finite (κ = ℓ²·1e30) before `min` caps it.

AD-clean: pure arithmetic, `φ` may be a `Dual`; `min` gives the consistent
(zero) dκ/dφ on the capped branch, so the ForwardDiff-extracted Jacobian
chain matches the residual everywhere — value and derivative saturate at the
same point.
"""
@inline function _carman_kozeny_permeability(φ::Real, ℓ::Real)
    one_minus = 1 - φ
    denom = 180 * one_minus * one_minus
    κ_ck = ℓ * ℓ * φ^3 / max(denom, 1e-30)
    return min(κ_ck, ℓ * ℓ)
end

# Cache-aware variant. `v` must be the renormalized solid basis
# (`solid_volume_fractions`), not the all-phase basis.
@inline function _permeability_from_v(v::Tuple, T::Real, material::Material{NC}) where {NC}
    mixing = material.mixing.κ.rule
    if mixing == BRUGGEMAN
        throw(
            ArgumentError(
                "BRUGGEMAN mixing is not supported for permeability (κ); EMT is only defined for scalar transport coefficients. Use PARALLEL/SERIES/WEIGHTED, or CARMAN_KOZENY (porosity-based) for permeability.",
            ),
        )
    elseif mixing == CARMAN_KOZENY
        throw(
            ArgumentError(
                "CARMAN_KOZENY permeability needs concentrations ξ to evaluate porosity; call effective_permeability(ξ, T, material, φ) rather than _permeability_from_v(v, ...).",
            ),
        )
    end

    W = _prop_type(v, T)
    κ_parallel, inv_κ_series = _permeability_sums(v, T, material.components, Val(NC))
    κ_series = inv_κ_series > 1e-30 ? convert(W, 1.0 / inv_κ_series) : zero(W)

    if mixing == PARALLEL
        return κ_parallel
    elseif mixing == SERIES
        return κ_series
    else  # WEIGHTED
        β = material.mixing.κ.β
        return β * κ_parallel + (1.0 - β) * κ_series
    end
end

# Guard of the smooth face mean [property units]: far below any physical
# conductivity or transfer coefficient.
const FACE_MEAN_GUARD = 1e-20

"""
	guarded_harmonic_mean(a, b, d_a, d_b, η = FACE_MEAN_GUARD)

Distance-weighted harmonic mean in the product form

	φ_face = (d_a + d_b) a b / (d_a b + d_b a + η (d_a + d_b))

smooth in every argument for non-negative values and positive distances. It
equals the resistance-in-series mean `(d_a + d_b)/(d_a/a + d_b/b)` to a
relative `η / φ`, and it goes to
zero continuously, with a finite derivative, as either value vanishes instead
of being cut off below a threshold. The guard `η` carries the property's units
and must lie far below its physical values: the default suits conductivities
and gas-transfer coefficients, a Darcy mobility needs a smaller one.
"""
@inline function guarded_harmonic_mean(
    a::Real,
    b::Real,
    d_a::Real,
    d_b::Real,
    η::Real = FACE_MEAN_GUARD,
)
    d = d_a + d_b
    return d * a * b / (d_a * b + d_b * a + η * d)
end
