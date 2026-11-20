# =============================================================================
# Phase change of water held in wood: the sealed cell
# =============================================================================
#
# Independent transcriptions of the closures of the phase change, and the
# references of a sealed, uniform cell built on them.
#
# - Saturation pressure of water: the vapour-pressure equation of Wagner and
#   Pruß (2002), J. Phys. Chem. Ref. Data 31, 387, eq. 2.5 (the auxiliary
#   equation of IAPWS-95), valid from the triple point to the critical point.
# - CAST isotherm, Dietenberger, Glass and Boardman (2024), eq. 18:
#   Δg(X, T) = [X (β₁ e^(−β₂/T) + β₃) + β₄]^(−β₅) − β₆ [kJ/mol], a = exp(−Δg/RT);
#   its heat of sorption, eq. 19, M Q_s = Δg + (Δg + β₆)/T · X β₁ β₂ β₅ e^(−β₂/T)/B
#   with B the bracket, which is (R T²/M) ∂ln a/∂T at constant X.
# - The model's two additions: the activity times X/(|X| + X₀), X₀ = 1e-4, with
#   the bracket evaluated at max(X, 0), and the soft cap
#   min_δ(x, 1) = 1 − δ ln(1 + e^((1 − x)/δ)), δ = 0.01, taken as x itself for
#   (1 − x)/δ > 40 and as one below −40; the heat carries the cap's
#   logarithmic slope x σ((1 − x)/δ)/min_δ(x, 1). Free water is the factor
#   alone, without the cap, and has no heat of sorption.
#
# The cell holds a sorbing solid at ξ_s, water at ξ_w, its vapour at ξ_v and
# air at ξ_a [kg/m³ bulk], every one uniform. Its pore fraction is
# φ = 1 − ξ_s/ρ_s − ξ_w/ρ_w, its moisture content X = ξ_w/(ξ_s + 10⁻³ kg/m³),
# and the vapour in equilibrium with the water is
# ξ_eq = φ M a(X, T) p_sat(T)/(R T). Sealed and uniform, nothing crosses its
# faces, and with ω = K (ξ_eq − ξ_v)
#
#   dξ_v/dt = ω,   ξ_w = ξ_w0 − (ξ_v − ξ_v0),   C dT/dt = −h ω,
#
# with C = Σ ξ_j c_j and h = L_ref + (c_v − c_w)(T − T_REF) + Q_s(X, T) for
# constant heat capacities. An athermal law (h ≡ 0) leaves T at its start and
# the first equation alone: its end state is the root of ξ_eq(ξ_w0 − x) = x,
# and its transient the quadrature t(ξ_v) = ∫ dx/(K (ξ_eq − x)), evaluated in
# s = −ln(1 − x/x∞), where the integrand stays finite up to the end state;
# past s = 20 the remaining gap, 2e-9 of the whole, relaxes linearly.
# The thermal case is integrated by classical Runge–Kutta at a step of
# 1/1000 of 1/K. For free water the total enthalpy
# Σ ξ_j c_j (T − T_REF) + L_ref ξ_v is conserved, so its end state is the root
# of one equation in ξ_v as well.

const PC_R = 8.31446261815324   # J/(mol K), CODATA 2018 (exact: N_A k_B)
const PC_T_REF = 298.15
const PC_X0 = 1e-4
const PC_DELTA = 0.01
const PC_SORBENT_FLOOR = 1e-3

# Table 3 of Dietenberger, Glass and Boardman (2024): β₁ … β₆ of the two sets.
const CAST_TABLE3 = (
    spruce = (304.2, 1661.0, 3.520, 0.3895, 2.886, 0.2006),
    beech = (46.33, 777.6, 2.028, 0.2210, 1.725, 0.3664),
)

"""Saturation pressure of water [Pa], Wagner and Pruß (2002) eq. 2.5."""
function pc_saturation(T)
    Tc, pc = 647.096, 22.064e6
    τ = 1 - T / Tc
    a = (-7.85951783, 1.84408259, -11.7866497, 22.6807411, -15.9618719, 1.80122502)
    s = a[1] * τ + a[2] * τ^1.5 + a[3] * τ^3 + a[4] * τ^3.5 + a[5] * τ^4 + a[6] * τ^7.5
    return pc * exp(Tc / T * s)
end

"""CAST eq. 18: Δg [J/mol] of the parameter set `β` at moisture `X` ≥ 0 and `T`."""
cast_gibbs(β, X, T) = 1e3 * ((X * (β[1] * exp(-β[2] / T) + β[3]) + β[4])^(-β[5]) - β[6])

"""CAST eq. 19: the heat of sorption [J/kg] of `β` at `X` ≥ 0 and `T`, vapour molar mass `M`."""
function cast_sorption_heat(β, X, T, M)
    e = exp(-β[2] / T)
    B = X * (β[1] * e + β[3]) + β[4]
    g = B^(-β[5]) - β[6]   # kJ/mol
    return 1e3 * (g + (g + β[6]) / T * X * β[1] * β[2] * β[5] * e / B) / M
end

"""The soft cap at one and its logarithmic slope."""
function pc_cap(x)
    z = (1 - x) / PC_DELTA
    z > 40 && return x
    z < -40 && return 1.0
    return 1 - PC_DELTA * log1p(exp(z))
end
function pc_cap_slope(x)
    z = (1 - x) / PC_DELTA
    z > 40 && return 1.0
    z < -40 && return 0.0
    return x / (1 + exp(-z)) / pc_cap(x)
end

pc_vanishing(X) = X / (abs(X) + PC_X0)

"""
    pc_activity(law, X, T)

The activity of `law`: `(:cast, β)`, the CAST isotherm of the set `β`;
`(:frozen, β, T_f)`, the CAST isotherm of `β` held at `T_f` whatever the
temperature; `(:free,)`, free water; `(:linear, X_s, X0)`, the linear law
X/X_s with the vanishing factor of scale `X0` and the cap.
"""
function pc_activity(law, X, T)
    kind = first(law)
    kind === :free && return pc_vanishing(X)
    kind === :linear && return pc_cap(max(X, 0.0) / law[2] * X / (abs(X) + law[3]))
    T_law = kind === :frozen ? law[3] : T
    return pc_cap(
        exp(-cast_gibbs(law[2], max(X, 0.0), T_law) / (PC_R * T_law)) * pc_vanishing(X),
    )
end

"""The heat of sorption [J/kg] of `law` (see [`pc_activity`](@ref))."""
function pc_sorption_heat(law, X, T, M)
    first(law) === :cast || return 0.0
    x = exp(-cast_gibbs(law[2], max(X, 0.0), T) / (PC_R * T)) * pc_vanishing(X)
    return cast_sorption_heat(law[2], max(X, 0.0), T, M) * pc_cap_slope(x)
end

"""
    SealedCell

A sealed, uniform cell: the sorbing solid at `ξ_s` of skeletal density `ρ_s`,
water at `ξ_w0` of density `ρ_w`, vapour at `ξ_v0` of molar mass `M`, air at
`ξ_a`, at `T0`; the activity `law`, exchange rate `K`, latent heat `L_ref`
at `T_REF`, heat capacities `c = (solid, water, vapour, air)`.
"""
Base.@kwdef struct SealedCell
    ξ_s::Float64
    ρ_s::Float64
    ξ_w0::Float64
    ρ_w::Float64 = 1000.0
    ξ_v0::Float64 = 0.0
    ξ_a::Float64
    M::Float64 = 0.018015
    T0::Float64
    law::Tuple
    K::Float64
    L_ref::Float64 = 2.4417e6
    c::NTuple{4,Float64}
end

pc_pore_fraction(cell::SealedCell, ξ_w) = 1 - cell.ξ_s / cell.ρ_s - ξ_w / cell.ρ_w
pc_moisture(cell::SealedCell, ξ_w) = ξ_w / (cell.ξ_s + PC_SORBENT_FLOOR)

"""The vapour [kg/m³] in equilibrium with water at `ξ_w` and `T`."""
pc_equilibrium_vapour(cell::SealedCell, ξ_w, T) =
    pc_pore_fraction(cell, ξ_w) *
    cell.M *
    pc_activity(cell.law, pc_moisture(cell, ξ_w), T) *
    pc_saturation(T) / (PC_R * T)

"""The vapour pressure [Pa] of vapour at `ξ_v` with water at `ξ_w` and `T`."""
pc_vapour_pressure(cell::SealedCell, ξ_w, ξ_v, T) =
    ξ_v / pc_pore_fraction(cell, ξ_w) * PC_R * T / cell.M

"""The equilibrium vapour pressure [Pa] of water at `ξ_w` and `T`."""
pc_equilibrium_pressure(cell::SealedCell, ξ_w, T) =
    pc_activity(cell.law, pc_moisture(cell, ξ_w), T) * pc_saturation(T)

_pc_water(cell::SealedCell, ξ_v) = cell.ξ_w0 - (ξ_v - cell.ξ_v0)

"""The end state of the athermal cell: the vapour [kg/m³] at equilibrium at `T0`."""
pc_athermal_end(cell::SealedCell) = bisect_root(
    x -> pc_equilibrium_vapour(cell, _pc_water(cell, x), cell.T0) - x,
    cell.ξ_v0,
    cell.ξ_v0 + cell.ξ_w0,
)

# The time [s] the athermal cell takes from its start to the vapour
# x∞ − (x∞ − ξ_v0) e^(−s), by Gauss–Legendre on 64 panels in s.
function _pc_athermal_time(cell::SealedCell, x∞, s)
    nodes, weights = _gauss_legendre(10)
    Δ = x∞ - cell.ξ_v0
    integrand(σ) = begin
        gap = Δ * exp(-σ)
        x = x∞ - gap
        gap / (cell.K * (pc_equilibrium_vapour(cell, _pc_water(cell, x), cell.T0) - x))
    end
    panels = 64
    h = s / panels
    return sum(
        h / 2 *
        sum(w * integrand((p - 0.5) * h + h / 2 * ξ) for (ξ, w) in zip(nodes, weights)) for
        p = 1:panels
    )
end

# Past s = 20 the gap to the end state is 2e-9 of the whole and the rate is
# linear in it to that order: the gap relaxes as e^(−λ t), λ = −K dg/dx at the
# end state, g(x) = ξ_eq(ξ_w0 − x) − x.
const PC_LINEAR_TAIL = 20.0

"""The vapour [kg/m³] of the athermal cell at time `t` [s], by quadrature."""
function pc_athermal_vapour(cell::SealedCell, t)
    t == 0 && return cell.ξ_v0
    x∞ = pc_athermal_end(cell)
    Δ = x∞ - cell.ξ_v0
    t_tail = _pc_athermal_time(cell, x∞, PC_LINEAR_TAIL)
    if t > t_tail
        g(x) = pc_equilibrium_vapour(cell, _pc_water(cell, x), cell.T0) - x
        h = 1e-6 * Δ
        λ = -cell.K * (g(x∞ + h) - g(x∞ - h)) / (2h)
        return x∞ - Δ * exp(-PC_LINEAR_TAIL) * exp(-λ * (t - t_tail))
    end
    s = bisect_root(σ -> _pc_athermal_time(cell, x∞, σ) - t, 0.0, PC_LINEAR_TAIL)
    return x∞ - Δ * exp(-s)
end

# The rates of the thermal cell at vapour `ξ_v` and temperature `T`.
function _pc_rates(cell::SealedCell, ξ_v, T)
    ξ_w = _pc_water(cell, ξ_v)
    ω = cell.K * (pc_equilibrium_vapour(cell, ξ_w, T) - ξ_v)
    cs, cw, cv, ca = cell.c
    C = cell.ξ_s * cs + ξ_w * cw + ξ_v * cv + cell.ξ_a * ca
    h =
        cell.L_ref +
        (cv - cw) * (T - PC_T_REF) +
        pc_sorption_heat(cell.law, pc_moisture(cell, ξ_w), T, cell.M)
    return ω, -h * ω / C
end

"""
    pc_thermal_history(cell, times; steps_per = 1000) -> (ξ_v, T)

The vapour [kg/m³] and temperature [K] of the thermal cell at `times`, each a
whole number of steps of 1/(steps_per K), by classical Runge–Kutta.
"""
function pc_thermal_history(cell::SealedCell, times; steps_per::Int = 1000)
    dt = 1 / (steps_per * cell.K)
    x, T, t = cell.ξ_v0, cell.T0, 0.0
    out_x, out_T = Float64[], Float64[]
    for target in times
        n = round(Int, (target - t) / dt)
        abs(n * dt - (target - t)) <= 1e-9 * max(target, dt) ||
            throw(ArgumentError("time $target is not a whole number of steps"))
        for _ = 1:n
            x, T = _pc_rk4_step(cell, x, T, dt)
        end
        t = target
        push!(out_x, x)
        push!(out_T, T)
    end
    return out_x, out_T
end

# One classical Runge–Kutta step of the thermal cell from (x, T).
function _pc_rk4_step(cell::SealedCell, x, T, dt)
    k1 = _pc_rates(cell, x, T)
    k2 = _pc_rates(cell, x + dt / 2 * k1[1], T + dt / 2 * k1[2])
    k3 = _pc_rates(cell, x + dt / 2 * k2[1], T + dt / 2 * k2[2])
    k4 = _pc_rates(cell, x + dt * k3[1], T + dt * k3[2])
    return x + dt / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1]),
    T + dt / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2])
end

"""
    pc_thermal_state(cell, t; steps_per = 1000) -> (ξ_v, T)

The vapour [kg/m³] and temperature [K] of the thermal cell at any time `t`,
by classical Runge–Kutta in the fewest equal steps no longer than
1/(steps_per K).
"""
function pc_thermal_state(cell::SealedCell, t; steps_per::Int = 1000)
    n = ceil(Int, t * steps_per * cell.K - 1e-9)
    n <= 0 && return cell.ξ_v0, cell.T0
    x, T = cell.ξ_v0, cell.T0
    for _ = 1:n
        x, T = _pc_rk4_step(cell, x, T, t / n)
    end
    return x, T
end

"""The total enthalpy [J/m³] of the cell, Σ ξ_j c_j (T − T_REF) + L_ref ξ_v."""
function pc_total_enthalpy(cell::SealedCell, ξ_v, T)
    ξ_w = _pc_water(cell, ξ_v)
    cs, cw, cv, ca = cell.c
    C = cell.ξ_s * cs + ξ_w * cw + ξ_v * cv + cell.ξ_a * ca
    return C * (T - PC_T_REF) + cell.L_ref * ξ_v
end

"""
    pc_free_water_end(cell) -> (ξ_v, T)

The end state of the adiabatic cell of free water, from the conservation of
the water and of the total enthalpy, with the vapour at equilibrium.
"""
function pc_free_water_end(cell::SealedCell)
    first(cell.law) === :free ||
        throw(ArgumentError("the end state from the enthalpy is free water's"))
    H0 = pc_total_enthalpy(cell, cell.ξ_v0, cell.T0)
    cs, cw, cv, ca = cell.c
    temperature(x) = begin
        C = cell.ξ_s * cs + _pc_water(cell, x) * cw + x * cv + cell.ξ_a * ca
        PC_T_REF + (H0 - cell.L_ref * x) / C
    end
    # Evaporation cools the cell, so its end lies below the vapour in
    # equilibrium at the start temperature.
    x = bisect_root(
        x -> pc_equilibrium_vapour(cell, _pc_water(cell, x), temperature(x)) - x,
        cell.ξ_v0,
        cell.ξ_v0 + 1.01 * pc_equilibrium_vapour(cell, cell.ξ_w0, cell.T0),
    )
    return x, temperature(x)
end

# =============================================================================
# The drying front of free water
# =============================================================================
#
# A slab of wood holding free water at moisture X, at T_i, has its face raised
# to T_s. The vapour leaves through the face at a total pressure P and carries
# no sensible heat. In the stiff limit the water leaves at a front held at the
# boiling temperature T_b, where a(X, T_b) p_sat(T_b) = P, and takes there the
# latent heat L(T_b) = L_ref + (c_v − c_w)(T_b − T_REF) per kilogram. Between
# the face and the front the slab is dry, of conductivity k_d and volumetric
# heat capacity C_d; beyond the front it is wet, k_w and C_w. This is the
# Neumann problem of exact/stefan.jl with the dry zone as its solid next to the
# wall at T∞ = T_s, the wet zone as its liquid at T0 = T_i, T_m = T_b, a unit
# density with the volumetric heat capacities, and ΔH_m = −ξ_w L(T_b): the
# front takes up, as it advances into the wet zone, the heat a solidification
# front gives out. Its eigencondition then reads
#
#   e^(−λ²)/erf λ − (k_w/k_d) ν (T_b − T_i)/(T_s − T_b) e^(−λ²ν²)/erfc(λν)
#       = √π λ ξ_w L(T_b) / (C_d (T_s − T_b)),   ν = √(α_d/α_w),
#
# with the front at s(t) = 2 λ √(α_d t) from the face.

"""
    pc_boiling_temperature(law, X, P) -> T_b [K]

The temperature at which water of moisture `X` under the activity `law` boils
at the total pressure `P` [Pa]: the root of a(X, T) p_sat(T) = P.
"""
pc_boiling_temperature(law, X, P) =
    bisect_root(T -> pc_activity(law, X, T) * pc_saturation(T) - P, 300.0, 500.0)

"""Kirchhoff's latent heat [J/kg] at `T` for the constant heat capacities `c_w` and `c_v`."""
pc_latent_heat(L_ref, c_w, c_v, T) = L_ref + (c_v - c_w) * (T - PC_T_REF)

"""
    pc_drying_front(; k_dry, k_wet, C_dry, C_wet, ξ_w, L, T_s, T_b, T_i) -> StefanProblem

The drying front as the Neumann problem of [`StefanProblem`](@ref): the dry
zone of conductivity `k_dry` and volumetric heat capacity `C_dry` between the
face at `T_s` and the front at `T_b`, the wet zone (`k_wet`, `C_wet`) at `T_i`
beyond it, and the water at `ξ_w` [kg/m³] taking the latent heat `L` [J/kg]
at the front. [`stefan_lambda`](@ref), [`stefan_front`](@ref) and
[`stefan_T`](@ref) then give the front and the temperatures, with depth
measured from the face.
"""
pc_drying_front(; k_dry, k_wet, C_dry, C_wet, ξ_w, L, T_s, T_b, T_i) = StefanProblem(
    k_s = k_dry,
    k_l = k_wet,
    ρ = 1.0,
    c_s = C_dry,
    c_l = C_wet,
    T_m = T_b,
    ΔH_m = -ξ_w * L,
    T0 = T_i,
    T∞ = T_s,
)

# =============================================================================
# Sorption into a slab through a face in contact with a humid reservoir
# =============================================================================
#
# A slab of thickness L, sealed at its back, holds water in equilibrium with
# its pore vapour at the activity a₀ and a uniform temperature. At t = 0 the
# face meets a reservoir at the activity a∞. The vapour moves by J = −Λ ∂ρ_v/∂z,
# ρ_v = a ρ_sat the vapour's density in the pores, and in the stiff limit the
# water follows it through the isotherm. Under the linear law a = X/X_s the
# water and its vapour, S X(a) + φ a ρ_sat per unit volume with S the sorbent
# concentration, are linear in ρ_v but for the pore fraction φ, which falls as
# the water's volume grows. With φ held, ρ_v obeys a diffusion equation of
#
#   D = Λ Δρ_v / (ΔW),   ΔW = S ΔX + Δ(φ ρ_v),
#
# the chord of the gain ΔW per unit volume against the step of the vapour,
# and the slab's gain per unit area is ΔW L times the uptake of a slab held at
# a fixed value at its face (Crank 1975, eq. 4.18),
#
#   M(t)/M∞ = 1 − Σ_{n≥0} 2/μₙ² exp(−μₙ² D t/L²),   μₙ = (2n + 1) π/2,
#
# the limit Bi → ∞ of the convective slab of exact/slab_series.jl. Under the
# dusty-gas closure the vapour, diffusing through air that does not move
# (Mason and Malinauskas 1983), has
#
#   1/Λ = 1/D_K + x_a/D_m,   D_K = (φ d/(3τ)) √(8 R T/(π M)),
#   D_m = (φ/τ) D_ref (T/T_ref)^1.75 P_ref/P,
#
# with x_a the air's mole fraction, d the pore diameter and τ the tortuosity.

"""
    pc_linear_moisture(a, X_s, X0) -> X

The moisture content at which the linear law `(:linear, X_s, X0)` reaches the
activity `a` below its cap: the positive root of X²/(X_s (X + X0)) = a.
"""
pc_linear_moisture(a, X_s, X0) = (a * X_s + sqrt((a * X_s)^2 + 4 * a * X_s * X0)) / 2

"""
    pc_slab_uptake(Fo) -> M/M∞

The fraction of its final gain that a slab sealed at its back has taken up
through a face held at a fixed value, at the Fourier number `Fo = D t/L²`
(Crank 1975, eq. 4.18). Summed until a term falls below 1e-17 with
μₙ² Fo > 50.
"""
function pc_slab_uptake(Fo)
    Fo <= 0 && return 0.0
    s = 0.0
    for n = 0:(10^7)
        μ = (2n + 1) * π / 2
        term = 2 / μ^2 * exp(-μ^2 * Fo)
        s += term
        term < 1e-17 && μ^2 * Fo > 50 && break
    end
    return 1 - s
end

"""
    pc_dusty_gas_transfer(; φ, τ, d, D_ref, T_ref, P_ref, T, P, M, x_a) -> Λ [m²/s]

The coefficient of the vapour's mass flux, J = −Λ ∂ρ_v/∂z, under the dusty-gas
model with one molecular diffusivity, the vapour (molar mass `M`) diffusing
through a background gas of mole fraction `x_a` that does not move:
1/Λ = 1/D_K + x_a/D_m, with the Knudsen coefficient of pores of diameter `d`
and the molecular diffusivity `D_ref` at `T_ref` and `P_ref`, both reduced by
`φ/τ`.
"""
function pc_dusty_gas_transfer(; φ, τ, d, D_ref, T_ref, P_ref, T, P, M, x_a)
    D_K = φ * d / (3τ) * sqrt(8 * PC_R * T / (π * M))
    D_m = φ / τ * D_ref * (T / T_ref)^1.75 * P_ref / P
    return 1 / (1 / D_K + x_a / D_m)
end
