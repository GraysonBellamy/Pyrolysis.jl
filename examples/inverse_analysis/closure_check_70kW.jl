"""
Closure check for the variable-cross-section plan prerequisites (Step 0).

Runs the 70 kW dry forward simulation with the best-fit conductivities and
reports:
  - Final mass fraction m_f/m_0           (measurement target: 0.238)
  - Final residue envelope bulk density   (measurement target: 91.1 kg/m³)

These are independent of the geometry refactor — they're properties of the
kinetics + fixed dilatometry-density sequence alone. They must match the 70 kW R1
post-test data before any geometry plumbing is touched.

Reference: docs/variable_cross_section_plan.md §3 (prerequisites) and §8
(verification gates 1, 2).
"""

using Printf
using Pyrolysis: StateLayout, to_problem_def
using Pyrolysis.Problem: chi_range

include("wp_forward.jl")
include("wp_objective.jl")

const RESULT_DIR = joinpath(@__DIR__, "results")

# Hardcoded from results/best.json (stage 1, [k_virgin, k_char])
const BEST_THETA = [9.256987e-02, 1.345222e-01]
const BEST_STAGE = 1

# ============================================================================
# Targets (UMD CAPA N₂ 70 kW R1 post-test)
# ============================================================================

const TARGET_M_RATIO = 0.238       # m_f / m_0
const TARGET_RHO_RES = 91.1        # kg/m³ envelope bulk
const TARGET_CHI_BAR = 0.65        # approximate kinetic conversion fraction
const M0_DATA = 11.482      # g (data file)
const MF_DATA = 2.684       # g (data file at t=595)
const D0_MM = 69.62       # pre-test diameter
const LF_MM = 9.20        # post-test thickness

# ============================================================================
# Load best-fit θ
# ============================================================================

function load_best_theta()
    @info "Using hardcoded best.json θ" stage=BEST_STAGE θ=BEST_THETA
    return copy(BEST_THETA), BEST_STAGE
end

# ============================================================================
# Compute end-of-test diagnostics
# ============================================================================

"""
Sum the mass of solid components in cell i.
Solid components are slots 1..7 (virgin, int1..4, char, ash) and 8 (moisture).
This excludes gases 9..14 which have already left through the surface.
"""
function cell_solid_mass(state, V, NC)
    m = 0.0
    @inbounds for j = 1:NC
        # Components 1..8 are condensed (solid + liquid moisture).
        # 9..14 are pyrolysis gases + water vapor (transient, will escape).
        if j <= 8
            m += state.ξ[j] * V
        end
    end
    return m
end

"""
Sum the volumetric concentration (sum of ξ over condensed components) in cell i.
This is the envelope bulk density of the cell's solid+liquid contents.
"""
function cell_bulk_density(state, NC)
    ρ = 0.0
    @inbounds for j = 1:NC
        if j <= 8
            ρ += state.ξ[j]
        end
    end
    return ρ
end

"""
Sum the mass of dry pyrolysis gases (slots 9..13, excluding water vapor at 14).
"""
function cell_dry_gas_mass(state, V, NC)
    m = 0.0
    @inbounds for j = 9:13
        if j <= NC
            m += state.ξ[j] * V
        end
    end
    return m
end

"""
Extract the final column-averaged dry-pyrolysis progress χ̄.

The state vector stores intensive per-cell `χ_i = m_released_dry,i /
m_dry_solid_initial,i`. The column diagnostic for the shrinkage law is the
mass-weighted average

	χ̄ = Σᵢ χ_i m_dry_solid_initial,i / Σᵢ m_dry_solid_initial,i

which is equivalent to cumulative released dry-pyrolysis-gas mass divided by
the initial dry-solid mass.
"""
function final_chi_diagnostic(sol)
    u_end = sol.raw_solution.u[end]
    mesh = sol.problem.mesh
    def = to_problem_def(sol.problem; use_ale = true)
    layout = StateLayout(NC_WP, length(mesh.cells), length(mesh.nodes), def.mode)
    rng = chi_range(layout)
    m_init = mesh.m_dry_solid_initial
    if isempty(rng) || isempty(m_init)
        return (
            χ̄ = 0.0,
            m_released_dry = 0.0,
            m_dry_initial = 0.0,
            χ_min = 0.0,
            χ_max = 0.0,
        )
    end
    χ = @view u_end[rng]
    n = min(length(χ), length(m_init))

    m_dry_initial = sum(m_init[i] for i = 1:n)
    m_released_dry = sum(χ[i] * m_init[i] for i = 1:n)

    return (
        χ̄ = m_released_dry / max(m_dry_initial, eps()),
        m_released_dry = m_released_dry,
        m_dry_initial = m_dry_initial,
        χ_min = minimum(χ[1:n]),
        χ_max = maximum(χ[1:n]),
    )
end

function compute_diagnostics(sol)
    mesh = sol.problem.mesh   # mesh has been mutated to final ALE state
    NC = NC_WP

    # === Mass diagnostics ===
    # `sol.extras.total_mass` is mass per cross-section area [kg/m²], all components.
    # At t=0 it includes only the initial solid (no gas yet); at t_end most gases
    # have escaped, so m_total[end]/A is approximately the residue per unit area.
    A0 = mesh.cross_section_area
    m_total_initial = sol.extras.total_mass[1] * A0   # kg
    m_total_final = sol.extras.total_mass[end] * A0 # kg

    # Per-cell condensed mass at end (slots 1..8)
    ξ_final = sol.ξ[:, :, end]   # NC × n_cells
    n_cells_active = size(ξ_final, 2)

    m_condensed_final = 0.0
    m_dry_gas_in_domain = 0.0
    for i = 1:n_cells_active
        V = mesh.cell_volumes[i]
        mesh.active_cell_mask[i] || continue
        for j = 1:8
            m_condensed_final += ξ_final[j, i] * V
        end
        for j = 9:13   # dry pyrolysis gases (water vapor at 14 excluded)
            m_dry_gas_in_domain += ξ_final[j, i] * V
        end
    end

    # Initial condensed mass in the domain (same calculation but with WP_XI0_DRY)
    ξ0 = WP_XI0_DRY
    ρ0_condensed = sum(ξ0[j] for j = 1:8)
    # Use initial volume = A_0 × L_0 from the scenario
    V0_total = A0 * sol.problem.tspan[2] * 0.0 + A0 * sol.extras.thickness_history[1]
    m_condensed_initial = ρ0_condensed * V0_total

    # === Density diagnostics ===
    # Bulk density per cell = Σ_{j=1..8} ξ_j (condensed only)
    # Top-region average for residue density (the heated zone is most-converted).
    n_top = max(1, n_cells_active ÷ 4)
    # Bottom = cell 1, top = cell n. Top region = last n_top cells.
    m_top = 0.0
    V_top = 0.0
    for i = (n_cells_active-n_top+1):n_cells_active
        mesh.active_cell_mask[i] || continue
        V = mesh.cell_volumes[i]
        ρ_cell = sum(ξ_final[j, i] for j = 1:8)
        m_top += ρ_cell * V
        V_top += V
    end
    ρ_top_avg = m_top / V_top

    # Column-averaged bulk density (residue + unreacted base)
    V_total_final =
        sum(mesh.cell_volumes[i] for i = 1:n_cells_active if mesh.active_cell_mask[i])
    ρ_column_avg = m_condensed_final / V_total_final

    # Per-cell ξ for inspection (top, mid, bottom)
    function cell_state_summary(i)
        ρ = sum(ξ_final[j, i] for j = 1:8)
        # Find dominant component
        max_j = argmax([ξ_final[j, i] for j = 1:8])
        labels = ["WP", "int1", "int2", "int3", "int4", "Char", "Ash", "Moisture"]
        return (
            cell = i,
            ρ_bulk = ρ,
            dominant = labels[max_j],
            frac_dominant = ξ_final[max_j, i] / max(ρ, 1e-30),
        )
    end

    cell_summaries = [
        cell_state_summary(1),                              # bottom
        cell_state_summary(n_cells_active ÷ 4),             # 25%
        cell_state_summary(n_cells_active ÷ 2),             # mid
        cell_state_summary(3 * n_cells_active ÷ 4),         # 75%
        cell_state_summary(n_cells_active),                 # top
    ]

    # Final thickness from solution
    L_final = sol.extras.thickness_history[end]
    χ_diag = final_chi_diagnostic(sol)

    return (
        m_total_initial = m_total_initial,
        m_total_final = m_total_final,
        m_condensed_initial = m_condensed_initial,
        m_condensed_final = m_condensed_final,
        m_dry_gas_in_domain = m_dry_gas_in_domain,
        m_ratio = m_total_final / m_total_initial,
        m_ratio_condensed = m_condensed_final / m_condensed_initial,
        L_final = L_final,
        ρ_top_avg = ρ_top_avg,
        ρ_column_avg = ρ_column_avg,
        t_final = sol.t[end],
        cell_summaries = cell_summaries,
        surface_temperature_final = sol.surface_T[end],
        back_temperature_final = sol.extras.back_temperature[end],
        χ̄ = χ_diag.χ̄,
        m_released_dry = χ_diag.m_released_dry,
        m_dry_initial = χ_diag.m_dry_initial,
        χ_min = χ_diag.χ_min,
        χ_max = χ_diag.χ_max,
    )
end

# ============================================================================
# Main
# ============================================================================

function main()
    θ, stage = load_best_theta()

    # UMD 70 kW dry case. sample_thickness = 8.5 mm matches the experimental
    # pre-test measurement L_0 = 8.50 mm.
    scenario, t_exp, T_back_exp, T_back_std, n_reps = load_scenario(
        "70kW",
        70_000.0;
        sample_thickness = 8.5e-3,
        n_cells = 160,
        verbose = true,
    )

    @info "Scenario" q_ext=scenario.q_ext t_final=scenario.t_final L_0_mm=scenario.sample_thickness*1000

    @info "Running forward solve..."
    t0 = time()
    sol = build_and_solve(θ, scenario; stage = stage, n_output = 51, verbose = false)
    dt = time() - t0
    @info "Solve complete" wall_time_s=round(dt; digits = 1)

    if sol === nothing
        @error "Forward solve failed"
        return
    end

    d = compute_diagnostics(sol)

    # Format as a report
    println()
    println("="^70)
    println("CLOSURE CHECK — 70 kW DRY (UMD CAPA N₂ R1)")
    println("="^70)
    @printf("Solve duration              : %.1f s\n", d.t_final)
    @printf(
        "Final thickness L_f         : %.2f mm    (data: %.2f mm)\n",
        d.L_final*1000,
        LF_MM
    )
    @printf("Final surface T             : %.0f K\n", d.surface_temperature_final)
    @printf("Final back T                : %.0f K\n", d.back_temperature_final)
    println("-"^70)
    println("PREREQUISITE GATE 1: kinetic residue (m_f / m_0)")
    @printf("  Simulated total mass init  : %.3f g\n", d.m_total_initial*1000)
    @printf(
        "  Simulated total mass final : %.3f g  (incl. dry gas in domain: %.4f g)\n",
        d.m_total_final*1000,
        d.m_dry_gas_in_domain*1000
    )
    @printf("  Simulated condensed final  : %.3f g\n", d.m_condensed_final*1000)
    @printf("  Simulated m_f/m_0 (total)  : %.4f\n", d.m_ratio)
    @printf("  Simulated m_f/m_0 (cond.)  : %.4f\n", d.m_ratio_condensed)
    @printf(
        "  Measured  m_f/m_0          : %.4f      (data: %.3f g / %.3f g)\n",
        MF_DATA/M0_DATA,
        MF_DATA,
        M0_DATA
    )
    gate1 = abs(d.m_ratio - MF_DATA/M0_DATA) < 0.03
    @printf(
        "  Δ (sim total − data)       : %+.4f    (%.1f%% relative)\n",
        d.m_ratio - MF_DATA/M0_DATA,
        100*(d.m_ratio - MF_DATA/M0_DATA)/(MF_DATA/M0_DATA)
    )
    @printf("  GATE 1 (±0.03 absolute)    : %s\n", gate1 ? "PASS" : "FAIL")
    println("-"^70)
    println("PREREQUISITE GATE 2: residue envelope bulk density")
    @printf(
        "  Simulated column-average ρ : %.1f kg/m³  (full column, residue + base)\n",
        d.ρ_column_avg
    )
    @printf(
        "  Simulated top-25%% average ρ: %.1f kg/m³  (heated zone — residue)\n",
        d.ρ_top_avg
    )
    @printf("  Measured residue ρ         : %.1f kg/m³  (envelope bulk)\n", TARGET_RHO_RES)
    gate2 = abs(d.ρ_top_avg - TARGET_RHO_RES) < 5.0
    @printf(
        "  Δ (top sim − data)         : %+.1f kg/m³  (%.1f%% relative)\n",
        d.ρ_top_avg - TARGET_RHO_RES,
        100*(d.ρ_top_avg - TARGET_RHO_RES)/TARGET_RHO_RES
    )
    @printf("  GATE 2 (±5 kg/m³ on top region): %s\n", gate2 ? "PASS" : "FAIL")
    println("-"^70)
    println("DIAGNOSTIC: dry-pyrolysis progress χ̄")
    @printf(
        "  Final mass-weighted χ̄      : %.3f       (expected ≈ %.2f)\n",
        d.χ̄,
        TARGET_CHI_BAR
    )
    @printf("  Released dry-gas numerator : %.3f g\n", d.m_released_dry*1000)
    @printf("  Initial dry-solid denom.   : %.3f g\n", d.m_dry_initial*1000)
    @printf("  Cell χ range               : %.3f – %.3f\n", d.χ_min, d.χ_max)
    println("-"^70)
    println("CELL STATE SAMPLES (bottom → top)")
    @printf("  %-10s %-12s %-12s %-12s\n", "cell", "ρ_bulk", "dominant", "frac")
    for s in d.cell_summaries
        @printf(
            "  %-10d %-12.1f %-12s %-12.3f\n",
            s.cell,
            s.ρ_bulk,
            s.dominant,
            s.frac_dominant
        )
    end
    println("-"^70)
    println("INTERPRETATION")
    if gate1 && gate2
        println("  Both prerequisite gates pass. Geometry refactor is the next step.")
    else
        println("  At least one prerequisite gate fails. The geometry refactor cannot")
        println("  rescue thickness — fix the failing gate first:")
        gate1 || println(
            "    - GATE 1: kinetics over/under-converts. Retune ν or rate constants.",
        )
        if !gate2
            println(
                "    - GATE 2: residue density disagrees with fixed dilatometry densities;",
            )
            println(
                "      inspect kinetics, void/swelling, or the component-density table.",
            )
        end
    end
    println("="^70)

    return d
end

main()
