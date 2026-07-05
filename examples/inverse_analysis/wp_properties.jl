"""
White Pine baseline thermophysical and kinetic properties.

Source: examples/inverse_analysis/data/macfp_white_pine_property_file.json
        (MaCFP White Pine, Bellamy et al.)

Single source of truth for the inverse analysis. Anything that the optimizer is
NOT adjusting lives here. The thermal conductivities (`WP_K_VIRGIN_INIT`,
`WP_K_CHAR_INIT`) are present as initial guesses and bounds, but their actual
values during optimization come from the parameter vector θ.

Some properties (emissivities, gas-phase properties, transport coefficients,
convective h's) are estimated pending measured values.
"""

# ============================================================================
# Initial composition (mass fractions on total-solid basis)
# ============================================================================

const WP_MASS_FRAC_VIRGIN = 0.996
const WP_MASS_FRAC_ASH = 0.004

# ============================================================================
# Densities — kg/m³
# ============================================================================

const WP_RHO_VIRGIN = 357.5
const WP_RHO_INT = (317.4, 283.9, 134.5, 83.2)
const WP_RHO_CHAR = 89.2
const WP_RHO_ASH = WP_RHO_CHAR  # JSON note: ash ρ = char ρ for porous-media solver

# Fixed component densities from dilatometry (kg/m³), using the scaled column:
#
# component  T_peak_C  Y_peak  rho_raw_exp  rho_scaled
# virgin        27      1.00     357.5        357.5
# int1         283      0.71     317.4        317.4
# int2         337      0.58     283.9        283.9
# int3         386      0.20     134.5        134.5
# int4         476      0.16      83.2         83.2
# char         812      0.15      89.2         89.2
#
# The old q-based front-loaded interpolation for intermediate density was a
# closure fit to cone-test residue geometry. These dilatometry values are
# direct component inputs, so density is no longer inferred from progress.
const WP_RHO_DILATOMETRY_RAW =
    (virgin = 357.5, int = (317.4, 283.9, 134.5, 83.2), char = 89.2)

# Initial mass concentrations (used to fill ξ_initial)
const WP_C0_VIRGIN = WP_MASS_FRAC_VIRGIN * WP_RHO_VIRGIN   # 358.56 kg/m³
const WP_C0_ASH = WP_MASS_FRAC_ASH * WP_RHO_VIRGIN   # 1.44   kg/m³

# Solid yields ν_i (fraction of solid passing to next intermediate per step)
# Source: JSON kinetics.reactions[*].nu
const WP_NU = (0.914206, 0.785307, 0.386929, 0.778861, 0.697933)

# Cumulative reactive char yield (final residue on ash-free basis)
const WP_FINAL_CHAR_YIELD = WP_NU[1] * WP_NU[2] * WP_NU[3] * WP_NU[4] * WP_NU[5]   # 0.151004

# Progress to char as cumulative fraction of conversion completed.
# Used for cp interpolation and stage-1/stage-3 conductivity interpolation.
const WP_INTERP_PROGRESS = [
    (1 - WP_NU[1]) / (1 - WP_FINAL_CHAR_YIELD),
    (1 - WP_NU[1]*WP_NU[2]) / (1 - WP_FINAL_CHAR_YIELD),
    (1 - WP_NU[1]*WP_NU[2]*WP_NU[3]) / (1 - WP_FINAL_CHAR_YIELD),
    (1 - WP_NU[1]*WP_NU[2]*WP_NU[3]*WP_NU[4]) / (1 - WP_FINAL_CHAR_YIELD),
]

# Downstream-progress to char (used for emissivity interpolation, matching
# the DF example convention)
const WP_DOWNSTREAM_PROGRESS = [
    WP_NU[2]*WP_NU[3]*WP_NU[4]*WP_NU[5],   # int1 → char
    WP_NU[3] * WP_NU[4] * WP_NU[5],            # int2 → char
    WP_NU[4] * WP_NU[5],                     # int3 → char
    WP_NU[5],                              # int4 → char
]

# ============================================================================
# Heat capacity — J/(kg·K)
# Source: JSON heat_capacity (regression intercepts/slopes are in kJ/(kg·K),
# converted here to J/(kg·K) via ×1000).
# ============================================================================

# Virgin: 3 explicit linear segments — HFM (low-T), bridge (mid-T), STA (high-T).
# All four breakpoints come straight from the JSON.
const WP_T_HFM_LO = 298.15
const WP_T_HFM_HI = 318.15
const WP_T_BRIDGE_HI = 425.0
const WP_T_STA_HI = 460.0

# HFM segment (T ≤ 318.15 K)
const WP_CP_HFM_INTERCEPT = -0.45048819329636997 * 1e3   # J/(kg·K)
const WP_CP_HFM_SLOPE = 0.006411747341451224 * 1e3   # J/(kg·K²)

# Bridge segment (318.15 < T ≤ 425.0 K)
const WP_CP_BRIDGE_INTERCEPT = 1.4415826349220475 * 1e3
const WP_CP_BRIDGE_SLOPE = 0.000464644313890584 * 1e3

# STA segment (T > 425.0 K) — extrapolated past 460 K for high-temp evaluation
const WP_CP_STA_INTERCEPT = -1.4474065850950482 * 1e3
const WP_CP_STA_SLOPE = 0.007262266008048457 * 1e3

function wp_virgin_cp(T)
    if T <= WP_T_HFM_HI
        WP_CP_HFM_INTERCEPT + WP_CP_HFM_SLOPE * T
    elseif T <= WP_T_BRIDGE_HI
        WP_CP_BRIDGE_INTERCEPT + WP_CP_BRIDGE_SLOPE * T
    else
        WP_CP_STA_INTERCEPT + WP_CP_STA_SLOPE * T
    end
end

# Char (and ash, by JSON convention)
const WP_CP_CHAR_INTERCEPT = 0.6234 * 1e3
const WP_CP_CHAR_SLOPE = 0.001253 * 1e3

wp_char_cp(T) = WP_CP_CHAR_INTERCEPT + WP_CP_CHAR_SLOPE * T

# Ash cp = char cp per JSON note
wp_ash_cp(T) = wp_char_cp(T)

# Intermediate cp interpolated from virgin and char by cumulative-yield progress
function wp_int_cp(T, progress)
    (1 - progress) * wp_virgin_cp(T) + progress * wp_char_cp(T)
end

# Water vapor cp — fixed (NIST polynomial fit, 300–1200 K)
wp_water_vapor_cp(T) = 1716.7 + 0.4154 * T + 1.530e-4 * T^2

# ============================================================================
# Pyrolysis-gas specific heats — J/(kg·K)
#
# Single source of truth for WP_gas{1..5}. Consumed by the GaseousComponent
# constructors in wp_forward.jl AND by the Kirchhoff Δc_p closures below.
# Keeping one definition prevents the two consumers from drifting — when gas
# cp is upgraded to T-dependent NIST polynomials, change wp_gas_cp(i, T) and
# both consumers update automatically.
# ============================================================================

const WP_GAS_CP = (1400.0, 1400.0, 1300.0, 1500.0, 1600.0)

@inline wp_gas_cp(i::Integer, T) = WP_GAS_CP[i]

# Moisture cp — kept for schema consistency; with zero initial moisture this
# is unused but the component still exists.
function wp_moisture_cp(T)
    t = T / 1000.0
    M = 0.01801528
    A, B, C, D, E = -203.6060, 1523.290, -3196.413, 2474.455, 3.855326
    (A + B*t + C*t^2 + D*t^3 + E/t^2) / M
end

# ============================================================================
# Heat of sorption — kept for schema consistency; with zero initial moisture
# the moisture reaction has nothing to react and `wp_heat_of_sorption` is
# never called.
# ============================================================================

# Index of moisture in the 14-component schema (see wp_forward.jl for ordering)
const WP_MOISTURE_IDX = 8
const WP_RHO_DRY_BASIS = WP_RHO_VIRGIN

function wp_heat_of_sorption(T, ξ)
    MC = 100.0 * ξ[WP_MOISTURE_IDX] / WP_RHO_DRY_BASIS
    (2450.420909987441 + 875.6441221735336 * exp(-0.1121684273316275 * MC)) * 1e3
end

# ============================================================================
# Kinetic parameters — 5-step pyrolysis + moisture evaporation
# Source: JSON kinetics.reactions
#   - A from log10(1/s) → s⁻¹: A = 10^A_log10
#   - E from kJ/mol → J/mol:   E ×= 1000
#   - h from kJ/kg → J/kg:     h ×= 1000   (positive = endothermic)
# Component-index pairs assume the 14-component ordering in wp_forward.jl:
#   1=Virgin 2-5=int1..4 6=Char 7=Ash 8=Moisture 9-13=gas1..5 14=H2O_vapor
# ============================================================================

# ----------------------------------------------------------------------------
# Kirchhoff-corrected Δh_r(T) constructor
# ----------------------------------------------------------------------------
#
# For each step, the DSC-measured Δh_r is the value at (approximately) the
# DSC rate-weighted-mean peak temperature. The cone reactions occur at
# higher T (~50–150 K above the DSC peak), so applying the constant DSC
# value over-extracts energy from the matrix per Kirchhoff's law. This
# helper builds a Δh_r(T) closure using the midpoint rule for the Kirchhoff
# integral, matching the package's existing idiom for the transpiration-BC
# enthalpy integral (components.jl:493-496).
#
# The returned closure has signature `T -> Real` — single-argument, so that
# `to_property(::Function)` wraps it in `CallableProperty`, whose 1- and
# 2-arg dispatch both forward `func(T)`. A `(T, ξ) -> ...` closure would
# MethodError under the 2-arg path; see property_functions.jl:191.
# ----------------------------------------------------------------------------
function kirchhoff_h(h_dsc::Real, T_peak_dsc::Real, Δcp::Function)
    return T -> begin
        T_mid = (T + T_peak_dsc) * 0.5
        h_dsc + Δcp(T_mid) * (T - T_peak_dsc)
    end
end

const WP_KINETICS = (
    (
        name = :pyrolysis_step1,
        indices = 1 => (2, 9),
        A = 10^5.224581,
        E = 76.009062 * 1e3,
        h = -3.072917986831273 * 1e3,
        yields = (WP_NU[1], 1 - WP_NU[1]),
    ),
    (
        name = :pyrolysis_step2,
        indices = 2 => (3, 10),
        A = 10^9.623536,
        E = 131.985833 * 1e3,
        h = 9.962405812813428 * 1e3,
        yields = (WP_NU[2], 1 - WP_NU[2]),
    ),
    (
        name = :pyrolysis_step3,
        indices = 3 => (4, 11),
        A = 10^12.086253,
        E = 173.004437 * 1e3,
        h = 223.10369172664198 * 1e3,
        yields = (WP_NU[3], 1 - WP_NU[3]),
    ),
    (
        name = :pyrolysis_step4,
        indices = 4 => (5, 12),
        A = 10^2.59569,
        E = 66.703667 * 1e3,
        h = -211.701917462439 * 1e3,
        yields = (WP_NU[4], 1 - WP_NU[4]),
    ),
    (
        name = :pyrolysis_step5,
        indices = 5 => (6, 13),
        A = 10^-1.525875,
        E = 23.336201 * 1e3,
        h = -8.737698321457055 * 1e3,
        yields = (WP_NU[5], 1 - WP_NU[5]),
    ),
)

# ============================================================================
# Kirchhoff-lifted kinetics — derived from WP_KINETICS via dh_r/dT = Δc_p
#
# WP_KINETICS above is preserved exactly as the published constant-Δh_r DSC
# deconvolution result. The tuples below are *derived* T-dependent
# interpretations of those same DSC numbers — not new measurements. They
# replace the constant `h` with a `kirchhoff_h(h_dsc, T_peak_dsc, Δcp)`
# closure that respects Kirchhoff's law per cell.
#
# T_peak,DSC values (WP_T_PEAK_DSC below) are the full 5-step coupled-ODE rate
# maxima at β = 10 K/min (the DSC heating rate), extracted by integrating the
# constant-Δh_r kinetics at that ramp and reading off each step's peak rate. A
# single-reaction TGA peak would be wrong for our sequential scheme (steps 2..5
# are rate-limited by upstream production).
#
# Two variants:
#   • WP_KINETICS_KIRCHHOFF_123 — lift steps 1, 2, 3 only; 4, 5 stay
#     constant. Captures ~45% of the budget bias; ~3% residual from
#     steps 4, 5 (per §5.2). Useful for diagnostic A/B against _ALL.
#   • WP_KINETICS_KIRCHHOFF_ALL — lift all five steps. Eliminates the
#     residual; recommended for production runs.
# ============================================================================

# Δc_p(T) closures, J/(kg·K) per kg of reactant. For each step:
#   Δc_p = ν_solid · c_p,solid_product + (1 − ν_solid) · c_p,gas − c_p,reactant
# The reactant cp goes negative; products are weighted by their mass yields.
const _Δcp_step1 =
    T ->
        WP_NU[1] * wp_int_cp(T, WP_INTERP_PROGRESS[1]) + (1 - WP_NU[1]) * wp_gas_cp(1, T) -
        wp_virgin_cp(T)
const _Δcp_step2 =
    T ->
        WP_NU[2] * wp_int_cp(T, WP_INTERP_PROGRESS[2]) + (1 - WP_NU[2]) * wp_gas_cp(2, T) -
        wp_int_cp(T, WP_INTERP_PROGRESS[1])
const _Δcp_step3 =
    T ->
        WP_NU[3] * wp_int_cp(T, WP_INTERP_PROGRESS[3]) + (1 - WP_NU[3]) * wp_gas_cp(3, T) -
        wp_int_cp(T, WP_INTERP_PROGRESS[2])
const _Δcp_step4 =
    T ->
        WP_NU[4] * wp_int_cp(T, WP_INTERP_PROGRESS[4]) + (1 - WP_NU[4]) * wp_gas_cp(4, T) -
        wp_int_cp(T, WP_INTERP_PROGRESS[3])
const _Δcp_step5 =
    T ->
        WP_NU[5] * wp_char_cp(T) + (1 - WP_NU[5]) * wp_gas_cp(5, T) -
        wp_int_cp(T, WP_INTERP_PROGRESS[4])

# DSC rate-weighted-mean peak temperatures (K), from the coupled-ODE rate maxima
# at β = 10 K/min (see note above).
const WP_T_PEAK_DSC = (530.1, 587.8, 638.4, 680.0, 788.0)

# Lifted heat closures, one per step. Same Δh_r,DSC values as WP_KINETICS;
# only the per-cell evaluation changes from constant to T-dependent.
const _h_kirch_1 = kirchhoff_h(-3.072917986831273e3, WP_T_PEAK_DSC[1], _Δcp_step1)
const _h_kirch_2 = kirchhoff_h(9.962405812813428e3, WP_T_PEAK_DSC[2], _Δcp_step2)
const _h_kirch_3 = kirchhoff_h(223.10369172664198e3, WP_T_PEAK_DSC[3], _Δcp_step3)
const _h_kirch_4 = kirchhoff_h(-211.701917462439e3, WP_T_PEAK_DSC[4], _Δcp_step4)
const _h_kirch_5 = kirchhoff_h(-8.737698321457055e3, WP_T_PEAK_DSC[5], _Δcp_step5)

# Variant A: lift steps 1, 2, 3 only.
const WP_KINETICS_KIRCHHOFF_123 = (
    (
        name = :pyrolysis_step1,
        indices = 1 => (2, 9),
        A = 10^5.224581,
        E = 76.009062 * 1e3,
        h = _h_kirch_1,
        yields = (WP_NU[1], 1 - WP_NU[1]),
    ),
    (
        name = :pyrolysis_step2,
        indices = 2 => (3, 10),
        A = 10^9.623536,
        E = 131.985833 * 1e3,
        h = _h_kirch_2,
        yields = (WP_NU[2], 1 - WP_NU[2]),
    ),
    (
        name = :pyrolysis_step3,
        indices = 3 => (4, 11),
        A = 10^12.086253,
        E = 173.004437 * 1e3,
        h = _h_kirch_3,
        yields = (WP_NU[3], 1 - WP_NU[3]),
    ),
    # Steps 4, 5: kept constant for diagnostic A/B against _ALL.
    (
        name = :pyrolysis_step4,
        indices = 4 => (5, 12),
        A = 10^2.59569,
        E = 66.703667 * 1e3,
        h = -211.701917462439 * 1e3,
        yields = (WP_NU[4], 1 - WP_NU[4]),
    ),
    (
        name = :pyrolysis_step5,
        indices = 5 => (6, 13),
        A = 10^-1.525875,
        E = 23.336201 * 1e3,
        h = -8.737698321457055 * 1e3,
        yields = (WP_NU[5], 1 - WP_NU[5]),
    ),
)

# Variant B: lift all five steps.
const WP_KINETICS_KIRCHHOFF_ALL = (
    (
        name = :pyrolysis_step1,
        indices = 1 => (2, 9),
        A = 10^5.224581,
        E = 76.009062 * 1e3,
        h = _h_kirch_1,
        yields = (WP_NU[1], 1 - WP_NU[1]),
    ),
    (
        name = :pyrolysis_step2,
        indices = 2 => (3, 10),
        A = 10^9.623536,
        E = 131.985833 * 1e3,
        h = _h_kirch_2,
        yields = (WP_NU[2], 1 - WP_NU[2]),
    ),
    (
        name = :pyrolysis_step3,
        indices = 3 => (4, 11),
        A = 10^12.086253,
        E = 173.004437 * 1e3,
        h = _h_kirch_3,
        yields = (WP_NU[3], 1 - WP_NU[3]),
    ),
    (
        name = :pyrolysis_step4,
        indices = 4 => (5, 12),
        A = 10^2.59569,
        E = 66.703667 * 1e3,
        h = _h_kirch_4,
        yields = (WP_NU[4], 1 - WP_NU[4]),
    ),
    (
        name = :pyrolysis_step5,
        indices = 5 => (6, 13),
        A = 10^-1.525875,
        E = 23.336201 * 1e3,
        h = _h_kirch_5,
        yields = (WP_NU[5], 1 - WP_NU[5]),
    ),
)

# Moisture evaporation reaction (kept for schema completeness; inactive when
# initial moisture concentration is zero)
const WP_MOISTURE_RXN =
    (name = :moisture_evaporation, indices = 8 => 14, A = 2.997e5, E = 5.0e4)

# ============================================================================
# Component schema fields not in the JSON (estimated pending measurements)
# ============================================================================
# Emissivity — these are ABSORPTIVITIES for the cone spectrum: the FSRI
# integrating-sphere reflectivity ε(λ)=1−ρ(λ) Planck-weighted at the cone
# effective source temperature, so they vary with q̇". SPECTRUM-VALIDATED:
# `derive_surface_emissivity.py` reproduces every row below to <0.003 from the
# raw spectra (data/FSRI_Wood*_REFLECTIVITY_*.csv) at the recovered cone temps
# T_cone = 939 / 1070 / 1163 K for 30 / 50 / 70 kW/m².
#
#   group (char by formation flux)  n     30 kW/m²   50 kW/m²   70 kW/m²
#   Virgin                          3     0.838      0.807      0.783
#   Char 20 kW                      7     0.686      0.678      0.675
#   Char 40 kW                      6     0.819      0.829      0.837
#   Char 60 kW                      5     0.826      0.834      0.839
#   Char 40+60 pooled               11    0.822      0.831      0.838   ← used for char
#
# Char defaults to the pooled dataset (largest n). To use the per-source
# datasets instead, swap the values in `WP_EPS_CHAR_TABLE`.
#
# CAVEAT 1 — char depends on FORMATION flux. char40 ≈ char60 (blackness
# saturates above ~40 kW), but char20 is ~0.14 lower. The pooled 40+60 table is
# right for the 50/70 kW tests (char forms at ≥40 kW) but likely OVER-states the
# 30 kW char absorptivity (char forms at ~30 kW, between char20 and char40, in
# the steep transition). No 30 kW-formed char spectrum exists to fix this; the
# 30 kW absorbed flux ε·q_ext may be a few % high late-time. Flagged, not fixed.
#
# CAVEAT 2 — emission reuses these absorptivities (gray approx). The same ε_eff
# drives both ε·q_ext (absorption, correct) and the −ε·σ·T_s⁴ loss (should use
# the surface-temperature emissivity). Re-integrating the spectra at T_s gives
# emission ε: virgin ~0.90–0.94, char ~0.80–0.82. Since the dominant (hot char)
# emitter already sits at ~0.82 and the large virgin emissivity applies only
# when the surface is cool (small εσT⁴), the gray reuse is bounded at ≲1% of
# flux — left as-is. See derive_surface_emissivity.py for the ε_emit(T_s) table.

# Each entry is (q̇" [W/m²], ε)
const WP_EPS_VIRGIN_TABLE = ((30.0e3, 0.838), (50.0e3, 0.807), (70.0e3, 0.783))
const WP_EPS_CHAR_TABLE = ((30.0e3, 0.822), (50.0e3, 0.831), (70.0e3, 0.838))

# Linear interpolation in q̇"; flat extrapolation outside the tabulated range.
function _interp_eps(table, q_ext)
    q = float(q_ext)
    q <= table[1][1] && return table[1][2]
    q >= table[end][1] && return table[end][2]
    for i = 1:(length(table)-1)
        qa, qb = table[i][1], table[i+1][1]
        if qa <= q <= qb
            f = (q - qa) / (qb - qa)
            return table[i][2] + f * (table[i+1][2] - table[i][2])
        end
    end
    error("unreachable: q=$q outside emissivity table")
end

wp_eps_virgin(q_ext) = _interp_eps(WP_EPS_VIRGIN_TABLE, q_ext)
wp_eps_char(q_ext) = _interp_eps(WP_EPS_CHAR_TABLE, q_ext)
# Ash has no separate measurement; treated as char (JSON convention).
wp_eps_ash(q_ext) = wp_eps_char(q_ext)

# Intermediate emissivity: interpolated between virgin and char by
# downstream progress-to-char, evaluated at the operating heat flux.
function wp_int_emissivity(progress_to_char, q_ext)
    εv = wp_eps_virgin(q_ext)
    εc = wp_eps_char(q_ext)
    εv + (εc - εv) * progress_to_char
end

# ============================================================================
# EMISSION emissivity ε(T_s) — distinct from the cone-spectrum ABSORPTIVITY
# tables above. The tables above are α weighted at the hot cone spectrum
# (T_cone ≈ 939–1163 K) and drive ε·q_ext (absorption). These below are ε
# weighted at the SURFACE's own Planck spectrum (T_s ≈ 300–1000 K) and drive
# the −ε·σ·T_s⁴ emission loss. They differ because wood/char are non-gray
# (blacker in the mid/far-IR): virgin emits at ~0.88–0.94, well above its
# cone absorptivity ~0.78–0.84, while char emission (~0.80–0.82) ≈ its
# absorptivity. SPECTRUM-DERIVED: reproduced by derive_surface_emissivity.py
# (EMISSION EMISSIVITY ε(T_s) block) from the same FSRI reflectivity spectra.
#
# NOTE: the production solver carries ONE composition-weighted ε_eff and uses
# it for both absorption (RadiativeFluxBC) and emission (RadiativeBC). The
# white-pine example splits them at the example level (see wp_forward.jl):
# component ε fields are set to these EMISSION values (evaluated once at
# WP_T_EMIT_REP — a static approximation; true ε is T_s-resolved), and the
# absorption term pins its own cone-spectrum absorptivity. Promoting this to a
# T_s-resolved, two-property treatment is a (small) core change.
const WP_EPS_EMIT_VIRGIN_TABLE = (
    (300.0, 0.940),
    (400.0, 0.929),
    (500.0, 0.914),
    (600.0, 0.899),
    (700.0, 0.884),
    (800.0, 0.867),
    (900.0, 0.847),
    (1000.0, 0.824),
)
const WP_EPS_EMIT_CHAR_TABLE = (
    (300.0, 0.783),
    (400.0, 0.792),
    (500.0, 0.797),
    (600.0, 0.801),
    (700.0, 0.806),
    (800.0, 0.812),
    (900.0, 0.819),
    (1000.0, 0.826),
)

# Linear interpolation in T_s; flat extrapolation (reuses _interp_eps, which is
# a generic monotonic-table linear interpolant).
wp_eps_emit_virgin(T_s) = _interp_eps(WP_EPS_EMIT_VIRGIN_TABLE, T_s)
wp_eps_emit_char(T_s) = _interp_eps(WP_EPS_EMIT_CHAR_TABLE, T_s)
wp_eps_emit_ash(T_s) = wp_eps_emit_char(T_s)

function wp_int_emit(progress_to_char, T_s)
    εv = wp_eps_emit_virgin(T_s)
    εc = wp_eps_emit_char(T_s)
    εv + (εc - εv) * progress_to_char
end

# Representative surface temperature for the static example-only emission ε.
# ≈ the 30 kW surface-peak band where the virgin emission correction matters
# most; emission ε is weakly T-dependent so a single anchor is adequate here.
const WP_T_EMIT_REP = 700.0

# ============================================================================
# Cone heater geometry — sample→cone-heater view factor
# ============================================================================
# Fraction of the sample-surface hemisphere subtended by the cone heater. The
# Schmidt-Boelter flux gauge calibrates against the cone-only contribution
# (cold-wall background self-cancels via gauge emission), so the calibration
# identity is q̇"_ext = F_gauge · ε_c · σ · T_cone⁴. The same view factor then
# splits the top-BC hemisphere into the cone share (incoming = ε·q_ext, outgoing
# emission absorbed by cone) and the chamber-wall share (net exchange against
# T_top_surr). NOTE: in the assembled BC the view factor cancels out of the
# emission loss (always full-hemisphere −ε·σ·T_s⁴); its ONLY surviving role is
# scaling the surroundings input (1−F)·ε·σ·T_top_surr⁴.
#
# Derived (not assumed) with the helical-coil contour-integral model of
#   Pinto et al., Int. J. Heat Mass Transfer 245 (2025) 126976,
#   doi:10.1016/j.ijheatmasstransfer.2025.126976
# for THIS geometry: 70 mm circular sample (R = 35 mm) sitting 40 mm below the
# cone base-plate bottom (+5 mm plate ⇒ model gap d = 45 mm), standard cone
# R1=R2=40 mm, H=65 mm. See `derive_cone_view_factor.py` (reproducible, with the
# published-example regression assert). Results:
#     F_center (r=0)            = 0.586   ← used here (back-TC is at center, so
#                                            the 1D column is the centerline)
#     F_sample,avg (disk-avg)   = 0.565
#     F_gauge (Schmidt-Boelter) = 0.586
# The earlier 0.656 was the standard-25 mm-spacing value and too high for the
# 40 mm CAPA gap. F_sample,avg/F_gauge ≈ 0.96 means the sample averages ~4% less
# than the gauge-nominal flux (it overhangs the ~66 mm aperture); that radial
# non-uniformity is a CENTERLINE-vs-MLR modeling choice and is deliberately NOT
# folded into q_ext here (the centerline column sees ≈ nominal).
const WP_F_CONE = 0.586

# Gas-phase transport properties — using DF defaults
# (Chapman-Enskog scaling: λ_gas ∝ √(M_ref/M))
const WP_LAMBDA_WOOD = 5.0e-7
const WP_LAMBDA_CHAR = 1.0e-5
const WP_LAMBDA_GAS = 3.0e-5
const WP_KAPPA_VIRGIN = 1.0e-15
const WP_KAPPA_CHAR = 1.0e-12

# ============================================================================
# Porosity anchors (pure-phase) — required by the porous-medium gas-pressure
# model `P = N(ξ,T) / φ`, where φ is computed from intrinsic (skeletal) cell-
# wall density: φ = 1 − Σ_solid (ξⱼ / ρⱼ_intrinsic). Bulk densities above are
# preserved as-is for thermal mixing rules / dilatometry-tied volume change;
# intrinsic densities live alongside them solely to feed the porosity formula.
# ============================================================================
#
# Virgin: φ = 0.53 — connected (Darcy-relevant) macroporosity, lumen + tracheid
#                    space accessible to gas flow. (Lower than total mercury-
#                    intrusion porosity because cell-wall nanopores aren't seen
#                    by Darcy flow.)
# Char:   φ = 0.94 — measured open porosity of the char skeleton.
# Ash:    inert, sits in the char matrix → uses char's φ.
#
# Intermediates: no direct porosity measurement. Intrinsic density is
# interpolated linearly in yield-progress between virgin and char anchors
# (consistent with the cp / k interpolation already used for intermediates).

const WP_PHI_VIRGIN = 0.53
const WP_PHI_CHAR = 0.94
const WP_PHI_ASH = WP_PHI_CHAR

# Derived intrinsic (skeletal) densities. Held to the user-provided pair
# (bulk, φ) — so ρ_intr = ρ_bulk / (1 − φ). For char this is 89.2/0.06 ≈ 1486,
# which differs slightly from the directly measured ρ_true_char ≈ 1400 because
# we hold WP_RHO_CHAR = 89.2 fixed (within the 6% measurement uncertainty).
const WP_RHO_INTRINSIC_VIRGIN = WP_RHO_VIRGIN / (1.0 - WP_PHI_VIRGIN)   # ≈ 760.6
const WP_RHO_INTRINSIC_CHAR = WP_RHO_CHAR / (1.0 - WP_PHI_CHAR)    # ≈ 1486.7
const WP_RHO_INTRINSIC_ASH = WP_RHO_INTRINSIC_CHAR

# Yield-progress interpolation for intermediates (no direct measurement)
function wp_rho_intrinsic_intermediate(progress::Real)
    WP_RHO_INTRINSIC_VIRGIN + (WP_RHO_INTRINSIC_CHAR - WP_RHO_INTRINSIC_VIRGIN) * progress
end

# ============================================================================
# Optimization parameter bounds and initial guesses
# ============================================================================

# Stage 1: [k_virgin, k_char] — 2 parameters
const WP_K_VIRGIN_INIT = 0.08021  # W/(m·K), updated stage-1 init
const WP_K_CHAR_INIT = 0.1889
const WP_K_VIRGIN_BOUNDS = (0.07, 0.4)   # 0.07 floor: just below the lowest measured prior (LFA-∥ 0.086) and the ~0.087–0.093 fit optimum; keeps it interior to Fminbox while excluding unphysical (<still-air) values
const WP_K_CHAR_BOUNDS = (0.03, 2.5)

# Stage 2: [k_virgin, k_int1..4, k_char] — 6 parameters
function wp_stage2_init()
    [
        WP_K_VIRGIN_INIT;
        [
            WP_K_VIRGIN_INIT + (WP_K_CHAR_INIT - WP_K_VIRGIN_INIT) * p for
            p in WP_INTERP_PROGRESS
        ];
        WP_K_CHAR_INIT
    ]
end
const WP_K_INT_BOUNDS = (0.05, 2.5)

# Stage 3: [k_v_lo, dk_v, k_c_lo, dk_c] — 4 parameters, k_x(T) globally linear,
# parameterized phase-specifically over the temperature window where each
# phase actually contributes mass to the energy balance:
#
#   virgin: k_v(T) = k_v_lo + ((T − T_V_LO) / (T_V_HI − T_V_LO)) · dk_v
#                    anchored on [298, 500] K — virgin dominates before pyrolysis,
#                    HFM/LFA priors live in this window
#   char:   k_c(T) = k_c_lo + ((T − T_C_LO) / (T_C_HI − T_C_LO)) · dk_c
#                    anchored on [500, 1100] K — char appears around T_C_LO
#                    and the cone-front operates in this range
#
# Intermediate phases use the convex combination k_i(T) = (1−p_i)·k_v(T) +
# p_i·k_c(T) by yield-progress, evaluated at the cell's actual T (no separate
# intermediate anchors — virgin's line and char's line are extrapolated as
# needed and combined). Mass-fraction weighting in `k_eff` makes irrelevant
# extrapolations harmless: virgin contributes ~0 mass at high T (where its
# extrapolation lives), char contributes ~0 mass at low T.
#
# Each dk lower bound sits slightly negative so init dk = 0 (flat k(T),
# equivalent to Stage 1) is interior to the box, avoiding the Fminbox μ→∞
# degeneracy. Bounds chosen so the worst-case linear extrapolation stays
# positive, with the high-T cap set by where each phase actually carries
# mass: virgin = 900 K (essentially fully consumed above this, so high-T
# extrap gets mass-weighted out of k_eff), char = 1200 K (dominates k_eff
# at high T, must stay positive there). Formula:
#     dk_min = -k_lo_min · (T_HI − T_LO) / (T_check − T_LO)
const WP_T_V_LO = 298.0    # K — virgin anchor lo (HFM measurement temp)
const WP_T_V_HI = 500.0    # K — virgin anchor hi (top of LFA prior range + small extrap)
const WP_T_C_LO = 500.0    # K — char anchor lo (~pyrolysis onset)
const WP_T_C_HI = 1100.0   # K — char anchor hi (~peak surface T at 70 kW/m²)
const WP_STAGE3_INIT = [
    0.1016,
    -0.02229, # virgin: k_lo, Δk
    0.1325,
    0.09838,  # char:   k_lo, Δk
]
const WP_STAGE3_DK_V_BOUNDS = (-0.025, 0.4)   # virgin: -0.08·202/(900−298)
const WP_STAGE3_DK_C_BOUNDS = (-0.0429, 2.5)   # char:   -0.05·600/(1200−500)
const WP_STAGE3_BOUNDS =
    [WP_K_VIRGIN_BOUNDS, WP_STAGE3_DK_V_BOUNDS, WP_K_CHAR_BOUNDS, WP_STAGE3_DK_C_BOUNDS]

# ============================================================================
# Variable cross-section — radial-shrinkage constitutive law
# ============================================================================
#
# FRAMING (volume-grounded — see derive_shrinkage.py for the full reduction).
# The CAPA tests give mass, diameter AND thickness per flux, so shrinkage is a
# closed, checkable system rather than one fitted area law:
#   • VOLUME is NOT fitted — the solver derives it from component densities,
#     dV/V = Σⱼ(γⱼ/ρⱼ)dξⱼ (src/Physics/volume_change.jl). The measured residue
#     density m_f/V_f independently confirms the density set: pre-test bulk ρ₀
#     = 358.9 ≈ WP_RHO_VIRGIN (+0.4%), and the ~fully-charred 70 kW residue
#     ρ_f = 86.2 ≈ WP_RHO_CHAR (−3.3%, within the ±6% char-ρ uncertainty).
#   • The RADIAL share A_f/A₀ is the ONLY empirically-fitted partition (below).
#   • AXIAL is EMERGENT, not "out of scope": L_f/L₀ = (V_f/V₀)/(A_f/A₀). Char
#     that swells axially at high flux (L_f/L₀ = 1.09 at 70 kW) is simply the
#     regime where radial shrinkage outruns the small volume loss. The
#     measured linear radial strain is ~flux-independent (~7%); ALL the flux
#     dependence is axial (shrink at 30 kW → swell at 70 kW). So: fit radial
#     once, let density drive volume, read axial off the balance.
#
# Pre-test diameters (UMD CAPA N₂, multi-replicate averages):
#
#   flux   n   D_0 (mm)   A_0 (mm²)
#   30     3   69.51      3796.0
#   50     4   69.43      3787.3
#   70     3   69.59      3804.8
#   all   10   69.50      3794.9
#
# Three-anchor LSQ refit (post-test diameters) against the MODEL-FREE
# experimental conversion χ̄_exp = 1 − m_f/m_0 (dry UMD samples; mass table
# below). Slope-through-origin LSQ in A_f/A_0 space, p = 1:
#
#   flux   D_f (mm)   A_f/A_0   χ̄_exp(1−m_f/m_0)   pred (p=1, s_r²=0.8005)
#   30 kW    65.0     0.8744    0.523               0.8957  (Δ +0.0212)  ← outlier
#   50 kW    64.4     0.8604    0.727               0.8550  (Δ −0.0054)
#   70 kW    64.4     0.8564    0.767               0.8470  (Δ −0.0094)
#
# RMS residual 0.0138 ≈ the ±0.013 D-measurement noise — the law sits AT the
# noise floor, not under it. Per-anchor s_r² = 0.760 / 0.808 / 0.813: 50 and
# 70 kW agree, 30 kW is a clear outlier. Dropping 30 kW gives s_r² = 0.8105
# (≈ the pre-refit value). A 2-parameter exponential-saturation form
# A/A_0 = s_r² + (1−s_r²)·exp(−χ̄/χ_*) degenerates to the linear limit when
# fitted (χ_* → ∞), so p = 1 stands; 3 scattered points cannot support a
# second shape parameter.
#
# WHY THE 30 kW OUTLIER — it is a CONVERSION mismatch, not a shrinkage one.
# The forward model predicts χ̄ = 0.648 at 30 kW, but the sample only lost
# 52.3% of its dry mass (χ̄_exp = 0.523): the kinetics over-pyrolyze at low
# flux. The PREVIOUS calibration anchored on the *simulated* χ̄
# (0.648/0.723/0.746), which collapses all three anchors onto one line
# (s_r² = 0.807, RMS 0.0004, "30× under noise") — but that tightness was
# partly circular: s_r² silently absorbed the model's 30 kW conversion bias.
# Anchoring on measured mass instead makes s_r² independent of the model's own
# χ̄; the cost is the now-visible 30 kW misfit, which correctly flags the
# kinetics rather than hiding it. Because the runtime law is still driven by
# the model's χ̄, the change to forward-model geometry is tiny (≈0.3% in
# diameter at every flux). If/when the low-flux kinetics are corrected so
# χ̄_sim → χ̄_exp, re-examine this fit.
#
# Note on the original single-anchor fit (s_r² = 0.778): used χ̄_anchor ≈ 0.65
# from end-of-test residue 23.8% on the R1 sample only (D_0 = 69.62); the
# multi-replicate D_0 average is 69.50.
#
# Axial thickness is NOT separately modeled or "out of scope" — it is the
# emergent remainder L_f/L₀ = (V_f/V₀)/(A_f/A₀), with V from the density model
# (validated above) and A from the radial law. The residual 30 kW disagreement
# (sim L_f ≈ 9.3 vs measured ≈ 7.6 mm) traces to the KINETICS conversion bias,
# not a missing axial term: at 30 kW the model's χ̄_sim ≈ 0.65 overshoots the
# measured χ̄_exp = 0.523, so it over-removes volume and (via the χ̄-driven
# area law) over-shrinks the radius. Fix the low-flux kinetics and the
# emergent L_f follows; see derive_shrinkage.py §[5] and validate_thickness.jl.
#
# Record: the former docs/variable_cross_section_plan.md has been removed;
# this block is the authoritative derivation.

const WP_D0 = 0.06950                    # m, overall pre-test mean (n=10)
const WP_A0 = π * (WP_D0 / 2)^2          # m² ≈ 3.795e-3
const WP_L0 = 8.5e-3                     # m, nominal pre-test thickness (== wp_forward sample_thickness); confirmed by ρ₀ = m₀/(A₀·L₀) ≈ WP_RHO_VIRGIN
# Per-flux pre-test diameters (use when matching a specific flux's mean sample)
const WP_D0_30 = 0.06951                    # m, n=3
const WP_D0_50 = 0.06943                    # m, n=4
const WP_D0_70 = 0.06959                    # m, n=3
wp_A0_for_flux(q_kW::Real) =
    π * (q_kW < 40 ? WP_D0_30 : q_kW < 60 ? WP_D0_50 : WP_D0_70)^2 / 4

const WP_SR2 = 0.8005                      # asymptotic radial fraction (χ̄ → 1); LSQ vs experimental χ̄ = 1 − m_f/m_0 (3-anchor)
const WP_SHRINK_P = 1.0                         # revisit if D_f vs χ̄ ever shows curvature

wp_lateral_shrinkage_law(χ̄) = 1.0 - (1.0 - WP_SR2) * χ̄^WP_SHRINK_P

# ============================================================================
# Multi-replicate post-test measurements (UMD CAPA N₂)
# ============================================================================
#
# Per-flux means used as calibration / validation anchors. Source: weighed and
# measured in lab post-test. All replicates dry-conditioned (N₂ chamber).
#
#   flux  n   m_0 (g)   Δm (g)    m_f/m_0   D_f (mm)   L_f (mm)
#   30    3   11.59     6.07      0.477     65.0       7.60  (n=3 thickness)
#   50    4   11.93     8.42      0.273     64.4       8.80  (n=4 thickness)
#   70    3   11.21     8.60      0.233     64.4       9.30  (n=3 thickness)

const WP_MASS_RATIO_30 = 0.477   # m_f/m_0, mean of 3 replicates (±0.026)
const WP_MASS_RATIO_50 = 0.273   # m_f/m_0, mean of 4 replicates (±0.010)
const WP_MASS_RATIO_70 = 0.233   # m_f/m_0, mean of 3 replicates (±0.008)

# Absolute pre-test masses (g) — needed for the measured residue bulk density
# ρ_f = m_f/V_f used as a density-model validation anchor.
const WP_M0_30_G = 11.59
const WP_M0_50_G = 11.93
const WP_M0_70_G = 11.21

const WP_LF_30_MM = 7.60    # post-test L mean (mm), 3 replicates
const WP_LF_50_MM = 8.80    # post-test L mean (mm), 4 replicates
const WP_LF_70_MM = 9.30    # post-test L mean (mm), 3 replicates

const WP_DF_30_MM = 65.00   # post-test D (mm), 1 replicate
const WP_DF_50_MM = 64.40
const WP_DF_70_MM = 64.40

# Measured A_f/A_0 per flux (uses per-flux pre-test D_0)
wp_Af_over_A0(q_kW::Real) =
    q_kW < 40 ? (WP_DF_30_MM / (WP_D0_30 * 1000))^2 :
    q_kW < 60 ? (WP_DF_50_MM / (WP_D0_50 * 1000))^2 : (WP_DF_70_MM / (WP_D0_70 * 1000))^2

# Measured axial L_f/L_0 per flux (pre-test L_0 = WP_L0, post-test L_f means)
wp_Lf_over_L0(q_kW::Real) =
    (q_kW < 40 ? WP_LF_30_MM : q_kW < 60 ? WP_LF_50_MM : WP_LF_70_MM) / (WP_L0 * 1000)

# Measured volume ratio V_f/V_0 = (A_f/A_0)(L_f/L_0) — the density-grounded
# quantity the solver predicts from component densities (volume_change.jl).
wp_Vf_over_V0(q_kW::Real) = wp_Af_over_A0(q_kW) * wp_Lf_over_L0(q_kW)

# Measured residue bulk density ρ_f = m_f / V_f (kg/m³). The ~fully-charred
# 70 kW value (≈86) is an independent check on WP_RHO_CHAR; the 30/50 kW
# values sit between virgin and char (partial conversion). See derive_shrinkage.py.
function wp_rho_residue(q_kW::Real)
    m0_g, Df_mm, Lf_mm =
        q_kW < 40 ? (WP_M0_30_G, WP_DF_30_MM, WP_LF_30_MM) :
        q_kW < 60 ? (WP_M0_50_G, WP_DF_50_MM, WP_LF_50_MM) :
        (WP_M0_70_G, WP_DF_70_MM, WP_LF_70_MM)
    mf_g = wp_measured_mass_ratio(q_kW) * m0_g
    Vf_mm3 = π * (Df_mm / 2)^2 * Lf_mm
    return (mf_g / Vf_mm3) * 1e6   # g/mm³ → kg/m³
end

wp_measured_mass_ratio(q_kW::Real) =
    q_kW < 40 ? WP_MASS_RATIO_30 : q_kW < 60 ? WP_MASS_RATIO_50 : WP_MASS_RATIO_70

# Look-up by flux for use in objective / validation. Includes the radial
# (Af_over_A0), axial (Lf_over_L0), and volume (Vf_over_V0) anchors plus the
# residue density (rho_f) so validation can check volume, not just radius.
function wp_measured_anchors(q_kW::Real)
    Lf_mm = q_kW < 40 ? WP_LF_30_MM : q_kW < 60 ? WP_LF_50_MM : WP_LF_70_MM
    (
        m_ratio = wp_measured_mass_ratio(q_kW),
        Lf_mm = Lf_mm,
        Af_over_A0 = wp_Af_over_A0(q_kW),
        Lf_over_L0 = wp_Lf_over_L0(q_kW),
        Vf_over_V0 = wp_Vf_over_V0(q_kW),
        rho_f = wp_rho_residue(q_kW),
    )
end
