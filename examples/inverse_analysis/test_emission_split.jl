"""
Before/after test of the absorption/emission emissivity split (example-only).

Runs the calibrated White Pine forward model at 30 kW (the over-conversion
case) and 70 kW (control) with the split OFF (baseline: cone absorptivity reused
for emission) and ON (emission ε(T_s) drives −εσT_s⁴; absorption pinned to the
char cone-α). Reports χ̄, mass ratio, MLR peaks, final back-T, peak surface T.

Hypothesis: ON should lower the 30 kW χ̄ toward the measured 0.523 (~half the
gap from ~0.65) and reduce the back-T overshoot, while leaving 70 kW ≈ unchanged.

Run:  julia --project=examples/inverse_analysis examples/inverse_analysis/test_emission_split.jl
"""

using Printf
using Pyrolysis: StateLayout, to_problem_def, column_chi_bar_from_state, ApproxSolve
include("wp_forward.jl")
include("wp_objective.jl")

const BEST_THETA = [0.09643, 0.09449, 0.07935, 0.1176, 0.2464, 0.2717]
const BEST_STAGE = 2

function chibar(sol)
    u_end = sol.raw_solution.u[end]
    mesh = sol.problem.mesh
    def = to_problem_def(sol.problem; use_ale = true)
    layout = StateLayout(NC_WP, length(mesh.cells), length(mesh.nodes), def.mode)
    return column_chi_bar_from_state(u_end, layout, mesh)
end

function run_one(label, q_ext, split::Bool)
    scenario_base, _, T_back_exp, _, _ = load_scenario(
        label,
        q_ext;
        n_cells = 160,
        verbose = false,
        use_variable_cross_section = true,
    )
    scenario = merge(scenario_base, (A_0 = wp_A0_for_flux(q_ext / 1000),))
    sol = build_and_solve(
        BEST_THETA,
        scenario;
        stage = BEST_STAGE,
        n_output = 51,
        verbose = false,
        use_emission_split = split,
        jacobian_strategy = ApproxSolve(),
    )
    sol === nothing && return nothing
    χ = try
        chibar(sol)
    catch
        NaN
    end
    return (
        χ = χ,
        mratio = sol.extras.total_mass[end] / sol.extras.total_mass[1],
        mlr_gauge = maximum(sol.extras.MLR_gauge) * 1000,
        mlr_local = maximum(sol.mass_loss_rate) * 1000,
        Tback = sol.extras.back_temperature[end],
        Tsurf = maximum(sol.surface_T),
        Texp_end = T_back_exp[end],
    )
end

function main()
    println("="^84)
    println("Absorption/emission split — baseline (OFF) vs split (ON)")
    println("="^84)
    for (label, q) in (("30kW", 30_000.0), ("70kW", 70_000.0))
        χ_exp = 1 - wp_measured_mass_ratio(q / 1000)
        @info "running $label ..."
        base = run_one(label, q, false)
        split = run_one(label, q, true)
        if base === nothing || split === nothing
            @warn "  $label failed";
            continue
        end
        println()
        println(
            "── $label   (measured χ̄_exp = $(round(χ_exp; digits=3)), " *
            "m_f/m_0 = $(round(wp_measured_mass_ratio(q/1000); digits=3)), " *
            "back-T_exp end = $(round(base.Texp_end; digits=0)) K)",
        )
        @printf("  %-16s %10s %10s %8s\n", "quantity", "OFF(base)", "ON(split)", "Δ")
        @printf("  %-16s %10.3f %10.3f %+8.3f\n", "χ̄", base.χ, split.χ, split.χ - base.χ)
        @printf(
            "  %-16s %10.3f %10.3f %+8.3f\n",
            "m_f/m_0",
            base.mratio,
            split.mratio,
            split.mratio - base.mratio
        )
        @printf(
            "  %-16s %10.2f %10.2f %+8.2f\n",
            "MLR_gauge pk",
            base.mlr_gauge,
            split.mlr_gauge,
            split.mlr_gauge - base.mlr_gauge
        )
        @printf(
            "  %-16s %10.2f %10.2f %+8.2f\n",
            "MLR_local pk",
            base.mlr_local,
            split.mlr_local,
            split.mlr_local - base.mlr_local
        )
        @printf(
            "  %-16s %10.1f %10.1f %+8.1f\n",
            "T_back end (K)",
            base.Tback,
            split.Tback,
            split.Tback - base.Tback
        )
        @printf(
            "  %-16s %10.1f %10.1f %+8.1f\n",
            "T_surf pk (K)",
            base.Tsurf,
            split.Tsurf,
            split.Tsurf - base.Tsurf
        )
    end
    println()
    println("="^84)
end

main()
