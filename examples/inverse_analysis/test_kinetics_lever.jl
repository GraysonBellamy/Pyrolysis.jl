"""
Quantify the early-kinetics lever on the 30 kW over-conversion.

The emission split moved χ̄ by ~0.001 — the 30 kW over-conversion is upstream,
in the early/low-T kinetics. This sweeps a per-step activation-energy bump on
the early cascade (steps 1–3) and measures χ̄ at 30/50/70 kW. Raising E throttles
the low-T rate more than the high-T rate, so the hypothesis is that it cuts the
long low-flux (30 kW) dwell conversion toward the measured χ̄_exp = 0.523 while
leaving 70 kW (χ̄_exp = 0.767) ≈ unchanged.

ApproxSolve Jacobian; emission split ON (default). Reports χ̄ and the residual
gap to the measured anchors for each (ΔE, flux).

Run:  julia --project=examples examples/inverse_analysis/test_kinetics_lever.jl
"""

using Printf
using Pyrolysis: StateLayout, to_problem_def, column_chi_bar_from_state, ApproxSolve
include("wp_forward.jl")
include("wp_objective.jl")

const BEST_THETA = [0.09643, 0.09449, 0.07935, 0.1176, 0.2464, 0.2717]
const BEST_STAGE = 2

# (label, per-step ΔE in J/mol). Steps 1–3 are the early/low-T cascade.
const SCENARIOS = [
    ("baseline      ", (0.0, 0.0, 0.0, 0.0, 0.0)),
    ("E123 +10 kJ/mol", (1.0e4, 1.0e4, 1.0e4, 0.0, 0.0)),
    ("E123 +20 kJ/mol", (2.0e4, 2.0e4, 2.0e4, 0.0, 0.0)),
]

const FLUXES = [("30kW", 30_000.0), ("50kW", 50_000.0), ("70kW", 70_000.0)]

function chibar(sol)
    u_end = sol.raw_solution.u[end]
    mesh = sol.problem.mesh
    def = to_problem_def(sol.problem; use_ale = true)
    layout = StateLayout(NC_WP, length(mesh.cells), length(mesh.nodes), def.mode)
    return column_chi_bar_from_state(u_end, layout, mesh)
end

function run_one(label, q_ext, dE)
    scenario_base, _, _, _, _ = load_scenario(
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
        kinetics_dE = dE,
        jacobian_strategy = ApproxSolve(),
    )
    sol === nothing && return (χ = NaN, mratio = NaN)
    χ = try
        chibar(sol)
    catch
        NaN
    end
    return (χ = χ, mratio = sol.extras.total_mass[end] / sol.extras.total_mass[1])
end

function main()
    println("="^88)
    println(
        "Early-kinetics lever — χ̄ vs ΔE(steps 1–3) at each flux  (ApproxSolve, split ON)",
    )
    println("="^88)
    χexp = Dict(q => 1 - wp_measured_mass_ratio(q / 1000) for (_, q) in FLUXES)
    @printf(
        "\n  measured χ̄_exp:  30kW=%.3f  50kW=%.3f  70kW=%.3f\n\n",
        χexp[30_000.0],
        χexp[50_000.0],
        χexp[70_000.0]
    )
    @printf(
        "  %-16s %18s %18s %18s\n",
        "ΔE scenario",
        "30kW χ̄ (gap)",
        "50kW χ̄ (gap)",
        "70kW χ̄ (gap)"
    )
    for (slabel, dE) in SCENARIOS
        cols = String[]
        for (flabel, q) in FLUXES
            @info "running $slabel @ $flabel ..."
            r = run_one(flabel, q, dE)
            gap = r.χ - χexp[q]
            push!(cols, @sprintf("%6.3f (%+5.3f)", r.χ, gap))
        end
        @printf("  %-16s %18s %18s %18s\n", slabel, cols[1], cols[2], cols[3])
    end
    println()
    println("  (gap = χ̄_sim − χ̄_exp; baseline 30 kW gap is the over-conversion to close)")
    println("="^88)
end

main()
