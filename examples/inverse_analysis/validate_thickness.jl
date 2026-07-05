"""
Validation driver — variable cross-section thickness anchor.

Runs the calibrated White Pine model (`wp_lateral_shrinkage_law`,
`s_r²=0.8005`, `p=1`, per-flux `D_0`) at 30 / 50 / 70 kW dry and reports:

  - Final thickness L_f                  (multi-replicate CAPA mean by flux)
  - Final radial ratio A_f/A_0           (diameter-derived CAPA mean by flux)
  - Mass ratio m_f/m_0                   (multi-replicate CAPA mean by flux)
  - MLR peaks (gauge + local)

The current radial law is an LSQ fit to the 30 / 50 / 70 kW diameter anchors.
Thickness and mass remain validation checks, not radial-law fit targets.

Note on θ: `BEST_THETA` was fit against the constant-A model. These
predictions therefore mix new geometry with stale conductivities; a later
inverse-fit re-run updates θ. Treat residuals in mass/thickness as
diagnostics for BCs, kinetics, fixed dilatometry densities, and missing void/swelling
physics rather than as radial-law fitting targets.
"""

using Printf
using Pyrolysis: StateLayout, to_problem_def, column_chi_bar_from_state

include("wp_forward.jl")
include("wp_objective.jl")

# Hardcoded fit params (stage 2, [k_virgin, k_int1..4, k_char]).
const BEST_THETA = [0.09643, 0.09449, 0.07935, 0.1176, 0.2464, 0.2717]
const BEST_STAGE = 2

# Multi-replicate measured means (UMD CAPA N₂, post-test).
# Canonical values live in wp_properties.jl; bind into a per-flux Dict here.
const MEAS = Dict(
    "30kW" => (
        Lf_mm = WP_LF_30_MM,
        Df_mm = WP_DF_30_MM,
        m_ratio = WP_MASS_RATIO_30,
        Af_over_A0 = wp_Af_over_A0(30.0),
    ),
    "50kW" => (
        Lf_mm = WP_LF_50_MM,
        Df_mm = WP_DF_50_MM,
        m_ratio = WP_MASS_RATIO_50,
        Af_over_A0 = wp_Af_over_A0(50.0),
    ),
    "70kW" => (
        Lf_mm = WP_LF_70_MM,
        Df_mm = WP_DF_70_MM,
        m_ratio = WP_MASS_RATIO_70,
        Af_over_A0 = wp_Af_over_A0(70.0),
    ),
)

# Acceptance window — keep wide; a residual gap is anticipated.
const LF_ACCEPT_MIN_MM = 8.0
const LF_ACCEPT_MAX_MM = 10.0

struct ScenarioResult
    label::String
    q_ext_kW::Float64
    t_final::Float64
    L0_mm::Float64
    Lf_mm::Float64
    A_f_over_A0::Float64
    m_ratio::Float64
    MLR_gauge_peak::Float64    # g/(m²·s)
    MLR_local_peak::Float64    # g/(m²·s)
    T_back_final::Float64      # K
    T_surface_peak::Float64    # K
    χ̄_final::Float64
    Vf_over_V0::Float64        # emergent volume ratio = (A_f/A_0)(L_f/L_0)
    rho_f::Float64             # residue bulk density ρ_0·(m_f/m_0)/(V_f/V_0), kg/m³
end

function compute_χ̄(sol)
    u_end = sol.raw_solution.u[end]
    mesh = sol.problem.mesh
    def = to_problem_def(sol.problem; use_ale = true)
    layout = StateLayout(NC_WP, length(mesh.cells), length(mesh.nodes), def.mode)
    return column_chi_bar_from_state(u_end, layout, mesh)
end

function run_scenario(label::AbstractString, q_ext::Real)
    sample_thickness = 8.5e-3
    n_cells = 160

    scenario_base, t_exp, T_back_exp, T_back_std, _ = load_scenario(
        label,
        q_ext;
        sample_thickness = sample_thickness,
        n_cells = n_cells,
        verbose = false,
        use_variable_cross_section = true,
    )
    # Override A_0 with the per-flux pre-test mean diameter so the sim
    # matches the population this flux's replicates were drawn from.
    A_0_flux = wp_A0_for_flux(q_ext / 1000)
    scenario = merge(scenario_base, (A_0 = A_0_flux,))
    @info "  scenario loaded" label q_ext_kW=q_ext/1000 t_final=scenario.t_final A_0_mm2=scenario.A_0*1e6
    @info "  D_0 from per-flux mean" D0_mm=2*sqrt(scenario.A_0/π)*1000

    t0 = time()
    sol = build_and_solve(
        BEST_THETA,
        scenario;
        stage = BEST_STAGE,
        n_output = 51,
        verbose = false,
    )
    dt = time() - t0
    @info "  solve complete" wall_time_s=round(dt; digits = 1)

    if sol === nothing
        @error "  forward solve failed for $label"
        return nothing
    end

    A_t = sol.extras.cross_section_area_history
    A_0 = A_t[1]
    A_f = A_t[end]

    mass = sol.extras.total_mass     # already per-A_0 (kg/m²)
    m_ratio = mass[end] / mass[1]

    # Use raw mesh state at end to read χ̄
    χ̄_final = try
        ;
        compute_χ̄(sol);
    catch
        ;
        NaN;
    end

    # Emergent volume ratio and residue bulk density. Volume is density-driven
    # in the solver (volume_change.jl); here we read the geometric outcome.
    # V_f/V_0 = (A_f/A_0)(L_f/L_0); ρ_f = ρ_0·(m_f/m_0)/(V_f/V_0) with the
    # initial bulk density ρ_0 = WP_RHO_VIRGIN.
    L_f_over_L0 = (sol.extras.thickness_history[end] * 1000) / (sample_thickness * 1000)
    Vf_over_V0 = (A_f / A_0) * L_f_over_L0
    rho_f = WP_RHO_VIRGIN * m_ratio / Vf_over_V0

    return ScenarioResult(
        String(label),
        q_ext / 1000,
        scenario.t_final,
        sample_thickness * 1000,
        sol.extras.thickness_history[end] * 1000,
        A_f / A_0,
        m_ratio,
        maximum(sol.extras.MLR_gauge) * 1000,
        maximum(sol.mass_loss_rate) * 1000,
        sol.extras.back_temperature[end],
        maximum(sol.surface_T),
        χ̄_final,
        Vf_over_V0,
        rho_f,
    )
end

# ------------------------------------------------------------------------
# Main
# ------------------------------------------------------------------------

function main()
    println("=" ^ 78)
    println(
        "Validation — variable-A White Pine (s_r²=$(WP_SR2), p=$(WP_SHRINK_P), per-flux D_0)",
    )
    println("=" ^ 78)
    @info "θ" θ=BEST_THETA stage=BEST_STAGE source="results/best.json (fit at constant A)"

    results = ScenarioResult[]
    for (label, q) in (("70kW", 70_000.0), ("50kW", 50_000.0), ("30kW", 30_000.0))
        @info "Running $label ..."
        try
            r = run_scenario(label, q)
            r === nothing || push!(results, r)
        catch err
            @warn "  scenario $label failed — skipping" exception=(err, catch_backtrace())
        end
    end

    isempty(results) && (println("No scenarios completed."); return nothing)

    println()
    println("RESULTS — sim vs measured (multi-replicate means)")
    println("-" ^ 100)
    @printf(
        "%-6s  %-15s  %-15s  %-15s  %-7s\n",
        "flux",
        "L_f sim/meas",
        "A_f/A_0 sim/m",
        "m_f/m_0 sim/m",
        "χ̄ sim"
    )
    @printf("%-6s  %-15s  %-15s  %-15s  %-7s\n", "", "(mm)", "(-)", "(-)", "(-)")
    for r in results
        m = get(MEAS, r.label, nothing)
        if m === nothing
            @printf(
                "%-6s  %5.2f /   —      %5.3f /   —      %5.3f /   —      %5.3f\n",
                r.label,
                r.Lf_mm,
                r.A_f_over_A0,
                r.m_ratio,
                r.χ̄_final
            )
        else
            @printf(
                "%-6s  %5.2f / %5.2f    %5.3f / %5.3f    %5.3f / %5.3f    %5.3f\n",
                r.label,
                r.Lf_mm,
                m.Lf_mm,
                r.A_f_over_A0,
                m.Af_over_A0,
                r.m_ratio,
                m.m_ratio,
                r.χ̄_final
            )
        end
    end
    println("-" ^ 100)

    # Per-flux residual table
    println()
    println("RESIDUALS (sim − measured)")
    @printf("  %-6s  %-12s  %-14s  %-14s\n", "flux", "ΔL_f (mm)", "ΔA_f/A_0", "Δm_f/m_0")
    for r in results
        m = get(MEAS, r.label, nothing)
        m === nothing && continue
        ΔL = r.Lf_mm - m.Lf_mm
        ΔA = r.A_f_over_A0 - m.Af_over_A0
        Δm = r.m_ratio - m.m_ratio
        @printf("  %-6s  %+7.2f       %+7.4f         %+7.4f\n", r.label, ΔL, ΔA, Δm)
    end

    # Volume / density decomposition — the density model drives V_f; the
    # radial law splits it; axial L_f is the remainder. Checking V_f/V_0 and
    # ρ_f directly catches density/kinetics errors that the separate radial
    # and thickness checks alone can hide. See derive_shrinkage.py.
    println()
    println("VOLUME / DENSITY (sim vs measured)")
    @printf("  %-6s  %-15s  %-17s\n", "flux", "V_f/V_0 sim/meas", "ρ_f sim/meas (kg/m³)")
    for r in results
        Vf_meas = wp_Vf_over_V0(r.q_ext_kW)
        rho_meas = wp_rho_residue(r.q_ext_kW)
        @printf(
            "  %-6s  %5.3f / %5.3f    %6.1f / %6.1f\n",
            r.label,
            r.Vf_over_V0,
            Vf_meas,
            r.rho_f,
            rho_meas,
        )
    end
    @printf(
        "  (density anchors: WP_RHO_VIRGIN=%.1f, WP_RHO_CHAR=%.1f; 70 kW ρ_f≈char)\n",
        WP_RHO_VIRGIN,
        WP_RHO_CHAR,
    )

    # Acceptance gate on the 70 kW point
    r70 = nothing
    for r in results
        r.label == "70kW" && (r70 = r; break)
    end
    if r70 !== nothing
        println()
        println(
            "ACCEPTANCE GATE  ($(LF_ACCEPT_MIN_MM) mm < L_f < $(LF_ACCEPT_MAX_MM) mm at 70 kW)",
        )
        passed = LF_ACCEPT_MIN_MM < r70.Lf_mm < LF_ACCEPT_MAX_MM
        println("  Result: $(passed ? "PASS" : "FAIL")")
    end

    println("=" ^ 92)
    return results
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
