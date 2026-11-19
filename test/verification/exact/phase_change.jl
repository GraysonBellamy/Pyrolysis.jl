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
temperature; `(:free,)`, free water.
"""
function pc_activity(law, X, T)
    kind = first(law)
    kind === :free && return pc_vanishing(X)
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
            k1 = _pc_rates(cell, x, T)
            k2 = _pc_rates(cell, x + dt / 2 * k1[1], T + dt / 2 * k1[2])
            k3 = _pc_rates(cell, x + dt / 2 * k2[1], T + dt / 2 * k2[2])
            k4 = _pc_rates(cell, x + dt * k3[1], T + dt * k3[2])
            x += dt / 6 * (k1[1] + 2k2[1] + 2k3[1] + k4[1])
            T += dt / 6 * (k1[2] + 2k2[2] + 2k3[2] + k4[2])
        end
        t = target
        push!(out_x, x)
        push!(out_T, T)
    end
    return out_x, out_T
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
