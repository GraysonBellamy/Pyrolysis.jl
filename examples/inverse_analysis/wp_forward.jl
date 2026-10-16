"""
Parameterized white pine forward model for inverse analysis.

`build_and_solve(θ, scenario; stage, kinetics_set)` builds a `Material` whose
conductivities are determined by θ (per the selected stage), wires up the
boundary conditions for the given experimental scenario, runs the simulation,
and returns the `PyrolysisSolution`. All other thermophysical properties and the
selectable white-pine kinetics sets come from `wp_properties.jl`.

Component schema (NC = 14):
    1 = WP (virgin)        7 = Ash (inert)        12 = WP_gas4
    2 = WP_int1            8 = Moisture           13 = WP_gas5
    3 = WP_int2            9 = WP_gas1            14 = WATER_VAPOR
    4 = WP_int3           10 = WP_gas2
    5 = WP_int4           11 = WP_gas3
    6 = Char

Ash (component 7) is inert (no reactions involve it). It sits in the char
matrix throughout: same ρ, cp, ε, and k as char (per the property JSON).
Wood and char-family component densities are fixed from the scaled
dilatometry table in `wp_properties.jl`.

`Scenario` encapsulates everything that varies between experiments (heat flux,
boundary-condition time-series interpolators, sample geometry, duration).
"""

using Pyrolysis
using OrdinaryDiffEq
# OrdinaryDiffEq v7 split the integrators into sublibraries; the umbrella package
# only re-exports the default solver set, so `KenCarp4` must be pulled from
# `OrdinaryDiffEqSDIRK` explicitly (see src/Solver/solve.jl docstring). Without
# this the `KenCarp4(...)` call in `build_and_solve` throws an UndefVarError that
# the surrounding try/catch silently turns into a failed solve.
using OrdinaryDiffEqSDIRK: KenCarp4
using LinearSolve
using Sparspak
using Logging
using LinearSolve: SparspakFactorization

include("wp_properties.jl")

const NC_WP = 14   # Number of components

# ============================================================================
# Scenario container
# ============================================================================

"""
	Scenario(; q_ext, T_top_gas, T_top_surr, T_bot_gas, T_bot_surr,
	         h_top, h_bot, sample_thickness, n_cells, t_final, label,
	         T_initial=nothing, c0_moisture=0.0)

Boundary-condition time series are callables `t -> K`. They should be flat
outside the data window (the loader in `wp_objective.jl` handles this).

`T_initial` (K) sets the uniform sample temperature at t=0. When `nothing`,
the solver falls back to `T_bot_surr(0.0)`.

`c0_moisture` (kg/m³) sets the initial uniform moisture concentration in
component slot 8. Defaults to 0 (dry sample). For e.g. 3% MC on dry-mass
basis, pass `c0_moisture = 0.03 * WP_RHO_DRY_BASIS` (= 10.8 kg/m³ for white
pine). The optimizer's `build_and_solve` reads this off the scenario, so
moisture can vary per scenario (e.g. only the FSRI sample is conditioned).
"""
function Scenario(;
    q_ext::Real,
    T_top_gas::Function,
    T_top_surr::Function,
    T_bot_gas::Function,
    T_bot_surr::Function,
    h_top::Real = 8.6,
    h_bot::Real = 4.0,
    ε_bottom::Real = 0.92,
    F_cone::Real = WP_F_CONE,
    sample_thickness::Real = 8.5e-3,
    # Pre-test cross-sectional area. Defaults to the multi-replicate
    # overall mean (D_0 = 69.50 mm, n=10 across 30/50/70 kW). Use
    # `wp_A0_for_flux(q_kW)` to override with the per-flux mean if needed.
    # NOTE: earlier versions defaulted to A = 1 m² via create_uniform_mesh;
    # D_0 = 73 mm (nominal) and D_0 = 69.62 mm (single R1 sample) are both
    # wrong defaults. See docs/variable_cross_section_plan.md §1.1.
    A_0::Real = WP_A0,
    n_cells::Int = 160,
    t_final::Real,
    label::String,
    T_initial::Union{Nothing,Real} = nothing,
    c0_moisture::Real = 0.0,
    # When false, the material is built with `lateral_shrinkage_law = nothing`
    # and the solver skips the χ-block tracking + dense ∂physics/∂χⱼ Jacobian
    # coupling — the inverse fit doesn't need post-test thickness/diameter and
    # pays a real per-evaluation cost for it. Validation drivers
    # (validate_thickness.jl, closure_check_70kW.jl) leave the default on.
    use_variable_cross_section::Bool = true,
)
    return (
        q_ext = float(q_ext),
        T_top_gas = T_top_gas,
        T_top_surr = T_top_surr,
        T_bot_gas = T_bot_gas,
        T_bot_surr = T_bot_surr,
        h_top = float(h_top),
        h_bot = float(h_bot),
        ε_bottom = float(ε_bottom),
        F_cone = float(F_cone),
        sample_thickness = float(sample_thickness),
        A_0 = float(A_0),
        n_cells = n_cells,
        t_final = float(t_final),
        label = label,
        T_initial = T_initial === nothing ? nothing : float(T_initial),
        c0_moisture = float(c0_moisture),
        use_variable_cross_section = use_variable_cross_section,
    )
end

# ============================================================================
# θ → conductivity functions (stage-dispatched)
# ============================================================================
#
# Returns a 6-tuple of conductivity property values
# (virgin, int1, int2, int3, int4, char). Ash uses the char value directly in
# `build_material`.

function _conductivities(θ::AbstractVector, ::Val{1})
    # Stage 1: [k_virgin, k_char] — intermediates by yield-progress interpolation
    kv, kc = θ[1], θ[2]
    (
        kv,
        kv + (kc - kv) * WP_INTERP_PROGRESS[1],
        kv + (kc - kv) * WP_INTERP_PROGRESS[2],
        kv + (kc - kv) * WP_INTERP_PROGRESS[3],
        kv + (kc - kv) * WP_INTERP_PROGRESS[4],
        kc,
    )
end

function _conductivities(θ::AbstractVector, ::Val{2})
    # Stage 2: [k_virgin, k_int1, k_int2, k_int3, k_int4, k_char] — independent
    (θ[1], θ[2], θ[3], θ[4], θ[5], θ[6])
end

function _conductivities(θ::AbstractVector, ::Val{3})
    # Stage 3: [k_v_lo, dk_v, k_c_lo, dk_c] with phase-specific anchor windows.
    # Each phase has its own slope-intercept (a, b); intermediates are the
    # yield-progress-weighted convex combination of the two phase lines, which
    # is itself linear in T (closed-form a_i, b_i below).
    k_v_lo, dk_v, k_c_lo, dk_c = θ[1], θ[2], θ[3], θ[4]

    b_v = dk_v / (WP_T_V_HI - WP_T_V_LO)
    a_v = k_v_lo - b_v * WP_T_V_LO
    b_c = dk_c / (WP_T_C_HI - WP_T_C_LO)
    a_c = k_c_lo - b_c * WP_T_C_LO

    int_ks = ntuple(4) do i
        p = WP_INTERP_PROGRESS[i]
        # k_i(T) = (1−p)·(a_v + b_v·T) + p·(a_c + b_c·T) — still linear in T
        a_i = (1 - p) * a_v + p * a_c
        b_i = (1 - p) * b_v + p * b_c
        CallableProperty(T -> a_i + b_i * T)
    end
    (CallableProperty(T -> a_v + b_v * T), int_ks..., CallableProperty(T -> a_c + b_c * T))
end

# Model k(T) evaluator used by priors / diagnostics
function model_k_at_T(
    θ::AbstractVector,
    T::Real,
    stage::Int,
    phase::Symbol;
    int_idx::Int = 1,
)
    if stage == 3
        # Phase-specific anchor windows. Each phase's k(T) is linear; each is
        # extrapolated as needed and combined for intermediates by yield-progress.
        k_v_lo, dk_v, k_c_lo, dk_c = θ[1], θ[2], θ[3], θ[4]
        f_v = (T - WP_T_V_LO) / (WP_T_V_HI - WP_T_V_LO)
        f_c = (T - WP_T_C_LO) / (WP_T_C_HI - WP_T_C_LO)
        k_v = k_v_lo + f_v * dk_v
        k_c = k_c_lo + f_c * dk_c
        if phase == :virgin
            return k_v
        elseif phase == :char
            return k_c
        elseif phase == :intermediate
            p = WP_INTERP_PROGRESS[int_idx]
            return (1 - p) * k_v + p * k_c
        else
            error("unknown phase: $phase")
        end
    elseif stage == 1
        kv, kc = θ[1], θ[2]
        if phase == :virgin
            return kv
        elseif phase == :char
            return kc
        elseif phase == :intermediate
            p = WP_INTERP_PROGRESS[int_idx]
            return kv + (kc - kv) * p
        else
            error("unknown phase: $phase")
        end
    elseif stage == 2
        # θ = [k_virgin, k_int1..4, k_char]
        if phase == :virgin
            return θ[1]
        elseif phase == :char
            return θ[6]
        elseif phase == :intermediate
            return θ[1+int_idx]
        else
            error("unknown phase: $phase")
        end
    else
        error("unsupported stage: $stage")
    end
end

# ============================================================================
# Material construction
# ============================================================================

function _wp_kinetics_set(kinetics_set::Symbol)
    if kinetics_set === :constant
        return WP_KINETICS
    elseif kinetics_set === :kirchhoff_123
        return WP_KINETICS_KIRCHHOFF_123
    elseif kinetics_set === :kirchhoff_all
        return WP_KINETICS_KIRCHHOFF_ALL
    else
        error("unsupported white-pine kinetics_set: $kinetics_set")
    end
end

function build_material(
    θ::AbstractVector,
    stage::Int,
    q_ext::Real;
    lateral_shrinkage_law = wp_lateral_shrinkage_law,
    kinetics_set::Symbol = :constant,
    use_emission_split::Bool = true,
    kinetics_dE = ntuple(_ -> 0.0, 5),
)
    ks = _conductivities(θ, Val(stage))
    k_virgin, k_int1, k_int2, k_int3, k_int4, k_char = ks
    kin = _wp_kinetics_set(kinetics_set)
    # Optional per-step activation-energy perturbation (J/mol), for sensitivity
    # studies of the early/low-T kinetics. Raising E throttles the low-T rate
    # more than the high-T rate, so it preferentially reduces conversion in the
    # long low-flux dwell. E-only (A is left alone); applied additively.
    if any(x -> x != 0.0, kinetics_dE)
        kin = ntuple(length(kin)) do i
            merge(kin[i], (E = kin[i].E + kinetics_dE[i],))
        end
    end
    # Ash is part of the char matrix; share k with char
    k_ash = k_char

    # Component surface emissivities. With the absorption/emission split on,
    # they carry the EMISSION values ε(T_s) (distinct from the cone-spectrum
    # ABSORPTIVITY, which is pinned on the top BC in build_bcs); ε_eff = Σⱼ vⱼ εⱼ
    # still tracks surface composition, and ε(T_s) is evaluated once at the
    # representative surface temp WP_T_EMIT_REP (static example-only split — see
    # the EMISSION block in wp_properties.jl). With the split OFF (baseline) the
    # components reuse the cone-spectrum absorptivity for emission, reproducing
    # the pre-split convention. Toggle via build_and_solve(use_emission_split=…).
    if use_emission_split
        ε_virgin = wp_eps_emit_virgin(WP_T_EMIT_REP)
        ε_char = wp_eps_emit_char(WP_T_EMIT_REP)
        ε_ash = wp_eps_emit_ash(WP_T_EMIT_REP)
        ε_int = ntuple(i -> wp_int_emit(WP_DOWNSTREAM_PROGRESS[i], WP_T_EMIT_REP), 4)
    else
        ε_virgin = wp_eps_virgin(q_ext)
        ε_char = wp_eps_char(q_ext)
        ε_ash = wp_eps_ash(q_ext)
        ε_int = ntuple(i -> wp_int_emissivity(WP_DOWNSTREAM_PROGRESS[i], q_ext), 4)
    end

    # cp closures — captured from wp_properties.jl
    virgin_cp_func = CallableProperty(wp_virgin_cp)
    char_cp_func = CallableProperty(wp_char_cp)
    ash_cp_func = CallableProperty(wp_ash_cp)
    int_cp_funcs = ntuple(4) do i
        p = WP_INTERP_PROGRESS[i]
        CallableProperty(T -> wp_int_cp(T, p))
    end
    moisture_cp_func = CallableProperty(wp_moisture_cp)
    water_vapor_cp_func = CallableProperty(wp_water_vapor_cp)

    hos_func = StateDependentProperty(wp_heat_of_sorption)

    # Intrinsic (skeletal) densities for the porosity formula. Bulk densities
    # above are kept for thermal mixing / dilatometry-tied volume change;
    # `ρ_intrinsic` only feeds φ = 1 − Σ ξⱼ/ρⱼ_intrinsic. Anchors are virgin and
    # char (measured); intermediates interpolate between them by yield-progress.
    ρ_intr_int = ntuple(i -> wp_rho_intrinsic_intermediate(WP_INTERP_PROGRESS[i]), 4)

    # The conductivity laws are linear in T with slopes drawn from the search
    # box, so positivity over the property bracket is the optimizer bounds'
    # job (WP_STAGE3_*_BOUNDS keep k > 0 over the anchor windows); a θ outside
    # that must reach the objective, not throw at construction.
    Material(
        name = :WhitePine,
        components = (
            # 1. Virgin white pine
            SolidComponent(
                :WP,
                ρ = WP_RHO_VIRGIN,
                c = virgin_cp_func,
                k = k_virgin,
                validate_properties = false,
                ε = ε_virgin,
                α = 0.0,
                λ = WP_LAMBDA_WOOD,
                κ = WP_KAPPA_VIRGIN,
                φ = WP_PHI_VIRGIN,
            ),
            # 2-5. Intermediates 1-4
            SolidComponent(
                :WP_int1,
                ρ = WP_RHO_INT[1],
                c = int_cp_funcs[1],
                k = k_int1,
                validate_properties = false,
                ε = ε_int[1],
                α = 0.0,
                λ = WP_LAMBDA_WOOD * 2.0,
                κ = WP_KAPPA_VIRGIN * 3.0,
                ρ_intrinsic = ρ_intr_int[1],
            ),
            SolidComponent(
                :WP_int2,
                ρ = WP_RHO_INT[2],
                c = int_cp_funcs[2],
                k = k_int2,
                validate_properties = false,
                ε = ε_int[2],
                α = 0.0,
                λ = WP_LAMBDA_WOOD * 5.0,
                κ = WP_KAPPA_VIRGIN * 10.0,
                ρ_intrinsic = ρ_intr_int[2],
            ),
            SolidComponent(
                :WP_int3,
                ρ = WP_RHO_INT[3],
                c = int_cp_funcs[3],
                k = k_int3,
                validate_properties = false,
                ε = ε_int[3],
                α = 0.0,
                λ = WP_LAMBDA_CHAR * 0.3,
                κ = WP_KAPPA_CHAR * 0.01,
                ρ_intrinsic = ρ_intr_int[3],
            ),
            SolidComponent(
                :WP_int4,
                ρ = WP_RHO_INT[4],
                c = int_cp_funcs[4],
                k = k_int4,
                validate_properties = false,
                ε = ε_int[4],
                α = 0.0,
                λ = WP_LAMBDA_CHAR * 0.7,
                κ = WP_KAPPA_CHAR * 0.1,
                ρ_intrinsic = ρ_intr_int[4],
            ),
            # 6. Char
            SolidComponent(
                :Char,
                ρ = WP_RHO_CHAR,
                c = char_cp_func,
                k = k_char,
                validate_properties = false,
                ε = ε_char,
                α = 0.0,
                λ = WP_LAMBDA_CHAR,
                κ = WP_KAPPA_CHAR,
                φ = WP_PHI_CHAR,
            ),
            # 7. Ash (inert, lives in char matrix — shares ρ/cp/ε/k with char)
            SolidComponent(
                :Ash,
                ρ = WP_RHO_ASH,
                c = ash_cp_func,
                k = k_ash,
                validate_properties = false,
                ε = ε_ash,
                α = 0.0,
                λ = WP_LAMBDA_CHAR,
                κ = WP_KAPPA_CHAR,
                φ = WP_PHI_ASH,
            ),
            # 8. Moisture (zero IC, but keep schema complete; pseudo-density
            # ρ=1e6 already gives ξ/ρ ≈ 0 contribution to porosity, no override needed)
            LiquidComponent(
                :MOISTURE,
                ρ = 1e6,
                c = moisture_cp_func,
                k = 0.6,
                ε = ε_virgin,
                α = 0.0,
                κ = 0.0,
            ),
            # 9-13. Pyrolysis gases (same speciation as DF; update if measured values exist)
            GaseousComponent(
                :WP_gas1,
                M = 0.060,
                c = WP_GAS_CP[1],
                k = 0.020,
                λ = WP_LAMBDA_GAS * sqrt(0.029/0.060),
            ),
            GaseousComponent(
                :WP_gas2,
                M = 0.070,
                c = WP_GAS_CP[2],
                k = 0.020,
                λ = WP_LAMBDA_GAS * sqrt(0.029/0.070),
            ),
            GaseousComponent(
                :WP_gas3,
                M = 0.110,
                c = WP_GAS_CP[3],
                k = 0.018,
                λ = WP_LAMBDA_GAS * sqrt(0.029/0.110),
            ),
            GaseousComponent(
                :WP_gas4,
                M = 0.070,
                c = WP_GAS_CP[4],
                k = 0.020,
                λ = WP_LAMBDA_GAS * sqrt(0.029/0.070),
            ),
            GaseousComponent(
                :WP_gas5,
                M = 0.050,
                c = WP_GAS_CP[5],
                k = 0.025,
                λ = WP_LAMBDA_GAS * sqrt(0.029/0.050),
            ),
            # 14. Water vapor
            GaseousComponent(
                :WATER_VAPOR,
                M = 0.018,
                c = water_vapor_cp_func,
                k = 0.025,
                λ = WP_LAMBDA_GAS * sqrt(0.029/0.018),
            ),
        ),
        reactions = (
            Reaction(
                kin[1].name,
                kin[1].indices,
                A = kin[1].A,
                E = kin[1].E,
                h = kin[1].h,
                n = 1.0,
                yields = kin[1].yields,
            ),
            Reaction(
                kin[2].name,
                kin[2].indices,
                A = kin[2].A,
                E = kin[2].E,
                h = kin[2].h,
                n = 1.0,
                yields = kin[2].yields,
            ),
            Reaction(
                kin[3].name,
                kin[3].indices,
                A = kin[3].A,
                E = kin[3].E,
                h = kin[3].h,
                n = 1.0,
                yields = kin[3].yields,
            ),
            Reaction(
                kin[4].name,
                kin[4].indices,
                A = kin[4].A,
                E = kin[4].E,
                h = kin[4].h,
                n = 1.0,
                yields = kin[4].yields,
            ),
            Reaction(
                kin[5].name,
                kin[5].indices,
                A = kin[5].A,
                E = kin[5].E,
                h = kin[5].h,
                n = 1.0,
                yields = kin[5].yields,
            ),
            Reaction(
                WP_MOISTURE_RXN.name,
                WP_MOISTURE_RXN.indices,
                A = WP_MOISTURE_RXN.A,
                E = WP_MOISTURE_RXN.E,
                h = hos_func,
                n = 1.0,
            ),
        ),
        mixing = MaterialMixing(
            k = MixingSpec(WEIGHTED, 0.5),
            λ = MixingSpec(WEIGHTED, 0.5),
            κ = MixingSpec(WEIGHTED, 0.5),
        ),
        lateral_shrinkage_law = lateral_shrinkage_law,
    )
end

# ============================================================================
# Boundary conditions
# ============================================================================

function build_bcs(
    scenario;
    use_emission_split::Bool = true,
    blowing::AbstractBlowingModel = CouetteFilm(),
)
    # Top BC — hemispherical decomposition. The cone subtends F_cone of the
    # sample's upper hemisphere; chamber walls fill the remaining (1−F_cone).
    # The Schmidt-Boelter flux gauge is calibrated cone-only (cold-wall
    # background cancels via gauge self-emission), so q_ext = F_cone·σ·T_cone⁴.
    #
    #   α · q_ext                                    absorbed from cone (cone-spectrum α)
    #   − ε · F_cone · σ · T_s⁴                      sample emission absorbed by cone
    #   ε · (1−F_cone) · σ · (T_top_surr⁴ − T_s⁴)    net exchange w/ chamber walls
    #   h · (T_top_gas − T_s)                        plume convection
    #
    # ABSORPTION vs EMISSION SPLIT (example-only). q_ext is absorbed with the
    # cone-spectrum absorptivity α (weighted at T_cone ≈ 939–1163 K), pinned
    # here to the char value wp_eps_char(q_ext) — the surface is char for the
    # bulk of every test, and at 70 kW virgin/char α differ by ~0.055, so a
    # char pin keeps the high-flux control clean (early-virgin α error is
    # ≤0.016 at 30 kW / brief at 70 kW). EMISSION uses ε_eff, which the
    # component ε fields now carry as the surface-temperature emission
    # emissivity (composition-tracked via mixing). The high-ε, cold,
    # water-cooled chamber walls make the (1−F_cone) surroundings a clean
    # blackbody sink, validating the bare −ε·σ·T_s⁴ emission form (no enclosure
    # radiosity factor). See wp_properties.jl EMISSION block + derive_shrinkage
    # discussion. The cone-direction RadiativeBC uses T_inf=0 to inject only
    # the sample's F_cone share of self-emission (cone re-emission already in
    # q_ext); validate=false because the validator rejects T_inf≤0.
    F_cone = scenario.F_cone
    # With the split on, pin the cone-spectrum absorptivity (char value);
    # otherwise fall back to ε_eff (baseline: absorptivity = emission = ε_eff).
    rad_flux =
        use_emission_split ?
        RadiativeFluxBC(flux = scenario.q_ext, absorptivity = wp_eps_char(scenario.q_ext)) :
        RadiativeFluxBC(flux = scenario.q_ext)
    top_thermal =
        rad_flux +
        RadiativeBC(F = F_cone, T_inf = 0.0, validate = false) +
        RadiativeBC(F = 1.0 - F_cone, T_inf = scenario.T_top_surr) +
        ConvectiveBC(h = scenario.h_top, T_inf = scenario.T_top_gas, blowing = blowing)
    top_mass = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)

    # Bottom: free to surroundings, impermeable.
    # Bottom emissivity is fixed at 0.92 (e.g. high-emissivity coating /
    # back surface measurement target), independent of decomposition state.
    bottom_thermal =
        RadiativeBC(ε = scenario.ε_bottom, T_inf = scenario.T_bot_surr) +
        ConvectiveBC(h = scenario.h_bot, T_inf = scenario.T_bot_gas, blowing = blowing)
    bottom_mass = ImpermeableBC()

    BoundaryConditionSet(
        thermal = (top = top_thermal, bottom = bottom_thermal),
        mass = (top = top_mass, bottom = bottom_mass),
    )
end

# ============================================================================
# Forward solve
# ============================================================================

# Dry-baseline initial concentration tuple — 14 components. Slot 8 (moisture)
# is overridden per scenario from `scenario.c0_moisture` inside `build_and_solve`.
const WP_XI0_DRY = (
    WP_C0_VIRGIN,           # 1. virgin (mass_frac × ρ_virgin)
    0.0,
    0.0,
    0.0,
    0.0,     # 2-5. intermediates
    0.0,                    # 6. char
    WP_C0_ASH,              # 7. ash
    0.0,                    # 8. moisture (overridden by scenario.c0_moisture)
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,# 9-13. pyrolysis gases
    0.0,                    # 14. water vapor
)

"""
	build_and_solve(θ, scenario; stage=1, n_output=601, verbose=false, kinetics_set=:constant, jacobian_strategy=nothing)

Build a white pine pyrolysis problem with conductivities determined by θ at the
given parameter stage, run it for `scenario.t_final` seconds, and return the
solution. Returns `nothing` if the solve fails (allowing the caller to apply a
penalty in the loss function).

Initial moisture concentration is read from `scenario.c0_moisture` (kg/m³).
For older Scenario tuples lacking that field, the slot defaults to 0.
`kinetics_set` selects `:constant`, `:kirchhoff_123`, or `:kirchhoff_all`.

`jacobian_strategy` overrides the Jacobian linear-solve strategy. `nothing`
(default) keeps the package default `Structured(strategy = DirectSolve())`.
Pass `ApproxSolve()` to opt into the `:ale_reduced` sparsity profile — it
drops the cumulative ALE mesh-velocity widening and the dense z-column
geometry block, paying ~14% of the prototype nnz and ~1% of the per-step
solve cost on MaCFP-shape problems. Newton absorbs the approximation by
taking a few extra iterations; for inverse analysis where the same problem
shape is solved hundreds of times this is a clear win.

`blowing` is the blowing model every `ConvectiveBC` of the surface BC set
carries, `CouetteFilm()` by default: the film coefficient is rescaled by the
Couette-film factor `B/(exp(B)−1)`, `B = ṁ″·c_p,g/h`, so outgassing pyrolysate
suppresses convective gain/loss at the hot face. The model sits on the BC, so
the residual and the Jacobian see the same physics. NOTE: changing the model
changes the forward model, so conductivities fit with another one are no
longer optimal — re-run the calibration.
"""
function build_and_solve(
    θ::AbstractVector,
    scenario;
    stage::Int = 1,
    n_output::Int = 601,
    verbose::Bool = false,
    log_level::Union{Nothing,Logging.LogLevel} = nothing,
    kinetics_set::Symbol = :constant,
    blowing::AbstractBlowingModel = CouetteFilm(),
    use_emission_split::Bool = true,
    kinetics_dE = ntuple(_ -> 0.0, 5),
    jacobian_strategy::Union{Nothing,LinearSolveStrategy} = nothing,
)

    use_var_cs = get(scenario, :use_variable_cross_section, true)
    law = use_var_cs ? wp_lateral_shrinkage_law : nothing
    material = build_material(
        θ,
        stage,
        scenario.q_ext;
        lateral_shrinkage_law = law,
        kinetics_set = kinetics_set,
        use_emission_split = use_emission_split,
        kinetics_dE = kinetics_dE,
    )
    bc_set = build_bcs(scenario; use_emission_split = use_emission_split, blowing = blowing)

    T0 = scenario.T_initial === nothing ? scenario.T_bot_surr(0.0) : scenario.T_initial

    T_initial_func = _ -> T0

    c0_moist = float(get(scenario, :c0_moisture, 0.0))

    ξ0 = (
        WP_C0_VIRGIN,           # 1
        0.0,
        0.0,
        0.0,
        0.0,     # 2-5
        0.0,                    # 6
        WP_C0_ASH,              # 7
        c0_moist,               # 8 — moisture (per-scenario)
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,# 9-13
        0.0,                    # 14
    )
    ξ_initial_funcs = [
        let v = ξ0[j]
            _ -> v
        end for j = 1:NC_WP
    ]

    mesh = create_uniform_mesh(
        scenario.sample_thickness,
        scenario.n_cells,
        Val(NC_WP);
        T_initial = T0,
        ξ_initial = ξ0,
        cross_section_area = get(scenario, :A_0, WP_A0),
    )

    problem = PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        T_initial = T_initial_func,
        ξ_initial = ξ_initial_funcs,
        tspan = (0.0, scenario.t_final),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, scenario.t_final, length = n_output)),
    )

    try
        # Build a JacobianSpec only when the caller asked for a non-default
        # strategy. With `jacobian = nothing`, `solve` picks
        # `Structured(strategy = DirectSolve())` itself; passing the spec
        # explicitly is required to switch to `ApproxSolve` (the
        # :ale_reduced profile).
        jac_spec =
            jacobian_strategy === nothing ? nothing :
            JacobianSpec(
                Structured(strategy = jacobian_strategy),
                problem;
                use_ale = true,
                radiation_model = SURFACE_ABSORPTION,
            )

        solver_kwargs = (
            abstol = 1e-7,
            reltol = 1e-7,
            # Cap dt and total step count so high-flux stiffness fails fast
            # instead of grinding for hours at microscopic dt. dtmax ~ saveat
            # interval; maxiters chosen with healthy naccept ~3k headroom.
            dtmin = 1e-8,
            dtmax = 2.0,
            maxiters = 200_000,
            progress = verbose,
            use_ale = true,
            radiation_model = SURFACE_ABSORPTION,
            min_thickness = 0.5e-3,
            handle_depletion = false,
            jacobian = jac_spec,
            integrator = KenCarp4(linsolve = SparspakFactorization()),
        )
        if log_level === nothing
            solution = solve(problem; solver_kwargs...)
        else
            logger = SimpleLogger(stderr, log_level)
            solution = with_logger(logger) do
                solve(problem; solver_kwargs...)
            end
        end
        return solution
    catch err
        verbose && @warn "Forward solve failed" exception=(err, catch_backtrace()) θ=θ
        return nothing
    end
end
