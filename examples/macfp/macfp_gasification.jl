"""
MaCFP-4 idealized 1D gasification target simulations for white pine.

Six runs covering all combinations of:
  - Incident radiant heat flux: 10, 30, 60 kW/m²
  - Sample thickness:           6, 18 mm

Setup (per MaCFP-4 §4.4.2 Table 2):
  - Initial sample / ambient temperature: 293 K
  - Sample cross-section:                 10 cm × 10 cm = 0.01 m²
  - Top BC:    incident radiant heat flux + grey-body re-radiation to ambient.
               No convection. Absorptivity = effective emissivity (Kirchhoff
               α=ε), where ε is the material's heat-flux-dependent surface
               emissivity evaluated at the case's incident q̇″. Re-radiation
               sink is −ε·σ·T_s⁴ against T_amb = 293 K (without it the
               problem has no equilibrium and surface T runs away).
  - Bottom BC: perfectly insulated (adiabatic, impermeable)
  - Initial moisture content: zero (dry sample)
  - Lateral shrinkage: disabled (idealized constant cross-section)

Outputs per run, in `results/<timestamp>/`:
  - macfp_<flux>kW_<thickness>mm.csv — per-case time-series CSV
        time_s, mass_g, mass_norm, T_top_K, T_back_K, thickness_mm,
        MLR_g_per_m2_s, gas_gen_g_per_m2_s,
        q_absorbed_W_per_m2, q_rerad_W_per_m2, q_conv_W_per_m2, q_cond_W_per_m2,
        eps_surf, k_surf_W_per_mK,
        wf_virgin, wf_int, wf_char, wf_ash, wf_gas, wf_water_vapor
  - macfp_<case>_detail.png   — 6-panel per-case diagnostic figure
  - macfp_summary.png         — cross-case overview (mass / MLR / T / thickness)

Conductivity parameters (Stage 2, [k_virgin, k_int1..4, k_char]) and all other
thermophysical/kinetic properties are loaded from `wp_properties.jl`.
"""

using Printf
using Dates
using Plots
using DelimitedFiles
using LinearAlgebra: BLAS

# Pin BLAS to a single thread so OpenBLAS doesn't spawn its own pool inside
# each parallel solve (standard practice for `Threads.@spawn`-based workloads).
BLAS.set_num_threads(1)

include(joinpath(@__DIR__, "..", "inverse_analysis", "wp_forward.jl"))

# ============================================================================
# Configuration
# ============================================================================

const THETA_STAGE2 = [
    0.09758105970837253,
    0.10016852221877959,
    0.07854572793279077,
    0.10422921145975779,
    0.2669146371686245,
    0.223603568269036,
]
const STAGE_USED = 2

const T_INITIAL_K = 293.0
const SAMPLE_AREA_M = 0.10            # 10 cm side
const A_SAMPLE_M2 = SAMPLE_AREA_M^2 # = 0.01 m²
const N_CELLS = 160

# (label, q_inc [kW/m²], thickness [mm], t_final [s])
# Durations sized generously per case; trim later from the trace if too long.
const CASES = [
    ("10kW_6mm", 10.0, 6.0, 2000.0),
    ("10kW_18mm", 10.0, 18.0, 5000.0),
    ("30kW_6mm", 30.0, 6.0, 1200.0),
    ("30kW_18mm", 30.0, 18.0, 2000.0),
    ("60kW_6mm", 60.0, 6.0, 800.0),
    ("60kW_18mm", 60.0, 18.0, 1500.0),
]

# ============================================================================
# MaCFP idealized boundary conditions
# ============================================================================

function build_macfp_bcs(q_inc::Real)
    # Top: incident radiant flux + grey-body re-radiation to ambient.
    # `RadiativeFluxBC(flux=q_inc)` gives absorbed flux = ε_eff·q_inc using
    # the material's effective emissivity (Kirchhoff α=ε). Per-component ε
    # is set to the q-dependent values inside `build_material(θ, stage, q)`,
    # so absorption tracks the case's incident flux.
    # `RadiativeBC(T_inf = T_amb)` is the −ε·σ·T_s⁴ re-emission sink with
    # F = 1 (full hemisphere to chamber walls at ambient). With ε=nothing
    # the BC also uses the material's effective ε, keeping absorption and
    # emission consistent. Without this term the problem has no thermal
    # equilibrium and the surface T runs away.
    # No convection (per MaCFP-4 idealized spec).
    top_thermal = RadiativeFluxBC(flux = float(q_inc)) + RadiativeBC(T_inf = T_INITIAL_K)
    top_mass = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)
    # Bottom: perfectly insulated, impermeable.
    bottom_thermal = AdiabaticBC()
    bottom_mass = ImpermeableBC()
    BoundaryConditionSet(
        thermal = (top = top_thermal, bottom = bottom_thermal),
        mass = (top = top_mass, bottom = bottom_mass),
    )
end

# ============================================================================
# Forward solve (MaCFP idealized variant of build_and_solve)
# ============================================================================

function build_and_solve_macfp(
    θ::AbstractVector,
    q_inc::Real,
    thickness_m::Real,
    t_final::Real;
    stage::Int = STAGE_USED,
    n_output::Int = 601,
    verbose::Bool = false,
)

    # Disable lateral shrinkage — MaCFP idealized scenario uses a constant
    # 10 cm × 10 cm cross-section.
    material = build_material(θ, stage, q_inc; lateral_shrinkage_law = nothing)
    bc_set = build_macfp_bcs(q_inc)

    T_initial_func = _ -> T_INITIAL_K

    ξ0 = (
        WP_C0_VIRGIN,           # 1. virgin
        0.0,
        0.0,
        0.0,
        0.0,     # 2-5. intermediates
        0.0,                    # 6. char
        WP_C0_ASH,              # 7. ash
        0.0,                    # 8. moisture (dry)
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,# 9-13. pyrolysis gases
        0.0,                    # 14. water vapor
    )
    ξ_initial_funcs = [
        let v = ξ0[j]
            _ -> v
        end for j = 1:NC_WP
    ]

    mesh = create_uniform_mesh(
        thickness_m,
        N_CELLS,
        Val(NC_WP);
        T_initial = T_INITIAL_K,
        ξ_initial = ξ0,
        cross_section_area = A_SAMPLE_M2,
    )

    problem = PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        T_initial = T_initial_func,
        ξ_initial = ξ_initial_funcs,
        tspan = (0.0, t_final),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, t_final, length = n_output)),
    )

    return solve(
        problem;
        abstol = 1e-9,
        reltol = 1e-7,
        dtmin = 1e-14,
        progress = verbose,
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        min_thickness = 0.5e-3,
        handle_depletion = false,
        # Default Jacobian is `Structured(strategy = DirectSolve())`.
        integrator = KenCarp4(linsolve = SparspakFactorization()),
    )
end

# ============================================================================
# Output writers
# ============================================================================

# Component groupings for diagnostic mass-fraction reporting.
# Schema matches wp_forward.jl (NC=14):
#   1=virgin  2-5=int1..4  6=char  7=ash  8=moisture
#   9-13=gas1..5  14=water_vapor
const VIRGIN_IDX = 1
const INT_IDXS = 2:5
const CHAR_IDX = 6
const ASH_IDX = 7
const GAS_IDXS = 9:13
const VAPOR_IDX = 14

# Integrate a per-cell concentration over the current cross section + cell
# volume contribution to give a per-component total mass [kg].
# We reconstruct cell volumes from the current thickness assuming uniform cells.
function _per_component_mass(sol, idx_set, k::Int)
    ξ = sol.ξ[:, :, k]   # [NC × n_cells]
    n_cells = size(ξ, 2)
    L = sol.extras.thickness_history[k]
    A =
        isempty(sol.extras.cross_section_area_history) ? A_SAMPLE_M2 :
        sol.extras.cross_section_area_history[k]
    V_cell = (L / n_cells) * A
    s = 0.0
    @inbounds for j = 1:n_cells, c in idx_set
        s += ξ[c, j] * V_cell
    end
    return s
end

function write_case_csv(path::AbstractString, sol)
    n = length(sol.t)
    mass_g = sol.extras.total_mass .* A_SAMPLE_M2 .* 1000   # g
    mass_norm = sol.extras.total_mass ./ sol.extras.total_mass[1]
    thick_mm = sol.extras.thickness_history .* 1000

    # Per-component totals → mass fractions of CURRENT total mass.
    wf_virgin = zeros(n)
    wf_int = zeros(n)
    wf_char = zeros(n)
    wf_ash = zeros(n)
    wf_gas = zeros(n)
    wf_vap = zeros(n)
    for k = 1:n
        m_v = _per_component_mass(sol, (VIRGIN_IDX,), k)
        m_i = _per_component_mass(sol, INT_IDXS, k)
        m_c = _per_component_mass(sol, (CHAR_IDX,), k)
        m_a = _per_component_mass(sol, (ASH_IDX,), k)
        m_g = _per_component_mass(sol, GAS_IDXS, k)
        m_w = _per_component_mass(sol, (VAPOR_IDX,), k)
        m_tot = m_v + m_i + m_c + m_a + m_g + m_w
        denom = m_tot > 0 ? m_tot : 1.0
        wf_virgin[k] = m_v / denom
        wf_int[k] = m_i / denom
        wf_char[k] = m_c / denom
        wf_ash[k] = m_a / denom
        wf_gas[k] = m_g / denom
        wf_vap[k] = m_w / denom
    end

    open(path, "w") do io
        println(
            io,
            join(
                [
                    "time_s",
                    "mass_g",
                    "mass_norm",
                    "T_top_K",
                    "T_back_K",
                    "thickness_mm",
                    "MLR_g_per_m2_s",
                    "gas_gen_g_per_m2_s",
                    "q_absorbed_W_per_m2",
                    "q_rerad_W_per_m2",
                    "q_conv_W_per_m2",
                    "q_cond_W_per_m2",
                    "eps_surf",
                    "k_surf_W_per_mK",
                    "wf_virgin",
                    "wf_int",
                    "wf_char",
                    "wf_ash",
                    "wf_gas",
                    "wf_water_vapor",
                ],
                ",",
            ),
        )
        for i = 1:n
            @printf(
                io,
                "%.4f,%.6f,%.6f,%.4f,%.4f,%.6f,%.6e,%.6e,%.6e,%.6e,%.6e,%.6e,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                sol.t[i],
                mass_g[i],
                mass_norm[i],
                sol.surface_T[i],
                sol.extras.back_temperature[i],
                thick_mm[i],
                sol.extras.surface_mass_flux[i] * 1000,
                sol.extras.gas_generation_rate[i] * 1000,
                sol.extras.surface_heat_flux_absorbed[i],
                sol.extras.surface_heat_flux_reradiation[i],
                sol.extras.surface_heat_flux_convection[i],
                sol.extras.surface_heat_flux_conduction[i],
                sol.extras.surface_emissivity[i],
                sol.extras.surface_conductivity[i],
                wf_virgin[i],
                wf_int[i],
                wf_char[i],
                wf_ash[i],
                wf_gas[i],
                wf_vap[i],
            )
        end
    end
end

function summary_plot(sols, out_dir::AbstractString)
    # Layout: 3 columns (10/30/60 kW), 4 rows (mass, MLR, T, thickness)
    fluxes_kW = [10.0, 30.0, 60.0]
    panels = Plots.Plot[]

    # Row 1: mass [g]
    for q in fluxes_kW
        p = plot(
            xlabel = "Time [s]",
            ylabel = "Mass [g]",
            title = @sprintf("%g kW/m² — Sample mass", q),
            legend = :topright,
        )
        for (lbl, qx, th, _) in CASES
            qx == q || continue
            sol = sols[lbl]
            sol === nothing && continue
            plot!(
                p,
                sol.t,
                sol.extras.total_mass .* A_SAMPLE_M2 .* 1000,
                lw = 2,
                label = @sprintf("%g mm", th)
            )
        end
        push!(panels, p)
    end

    # Row 2: surface MLR [g/(m²·s)]
    for q in fluxes_kW
        p = plot(
            xlabel = "Time [s]",
            ylabel = "MLR [g/(m²·s)]",
            title = @sprintf("%g kW/m² — Surface MLR", q),
            legend = :topright,
            ylims = (0, Inf),
        )
        for (lbl, qx, th, _) in CASES
            qx == q || continue
            sol = sols[lbl]
            sol === nothing && continue
            plot!(
                p,
                sol.t,
                sol.extras.surface_mass_flux .* 1000,
                lw = 2,
                label = @sprintf("%g mm", th)
            )
        end
        push!(panels, p)
    end

    # Row 3: top + back temperatures
    for q in fluxes_kW
        p = plot(
            xlabel = "Time [s]",
            ylabel = "T [K]",
            title = @sprintf("%g kW/m² — Surface & back T", q),
            legend = :bottomright,
        )
        for (lbl, qx, th, _) in CASES
            qx == q || continue
            sol = sols[lbl]
            sol === nothing && continue
            plot!(
                p,
                sol.t,
                sol.surface_T,
                lw = 2,
                label = @sprintf("%g mm — top", th),
                linestyle = :solid,
            )
            plot!(
                p,
                sol.t,
                sol.extras.back_temperature,
                lw = 2,
                label = @sprintf("%g mm — back", th),
                linestyle = :dash,
            )
        end
        push!(panels, p)
    end

    # Row 4: thickness [mm]
    for q in fluxes_kW
        p = plot(
            xlabel = "Time [s]",
            ylabel = "Thickness [mm]",
            title = @sprintf("%g kW/m² — Sample thickness", q),
            legend = :topright,
        )
        for (lbl, qx, th, _) in CASES
            qx == q || continue
            sol = sols[lbl]
            sol === nothing && continue
            plot!(
                p,
                sol.t,
                sol.extras.thickness_history .* 1000,
                lw = 2,
                label = @sprintf("%g mm", th)
            )
        end
        push!(panels, p)
    end

    fig = plot(
        panels...;
        layout = (4, 3),
        size = (1500, 1400),
        left_margin = 8Plots.mm,
        bottom_margin = 6Plots.mm,
        plot_title = "MaCFP-4 White Pine — Idealized 1D Gasification",
    )
    out = joinpath(out_dir, "macfp_summary.png")
    savefig(fig, out)
    return out
end

# Per-case 6-panel diagnostic figure
function detail_plot(
    case_label::AbstractString,
    q_kW::Real,
    th_mm::Real,
    sol,
    out_dir::AbstractString,
)

    title_root = @sprintf("White Pine — %g kW/m², %g mm", q_kW, th_mm)

    # Panel 1: mass [g] + normalized
    p1 = plot(
        xlabel = "Time [s]",
        ylabel = "Mass [g]",
        title = "Sample mass",
        legend = :bottomleft,
    )
    plot!(
        p1,
        sol.t,
        sol.extras.total_mass .* A_SAMPLE_M2 .* 1000,
        lw = 2.2,
        color = :navy,
        label = "mass",
    )
    p1b = twinx(p1)
    plot!(
        p1b,
        sol.t,
        sol.extras.total_mass ./ sol.extras.total_mass[1],
        lw = 1.6,
        color = :gray,
        linestyle = :dash,
        label = "m/m₀",
        ylabel = "m/m₀ [-]",
        ylims = (0, 1.05),
    )

    # Panel 2: MLR + gas-generation rate (per m²)
    p2 = plot(
        xlabel = "Time [s]",
        ylabel = "[g/(m²·s)]",
        title = "Mass loss & gas generation",
        legend = :topright,
        ylims = (0, Inf),
    )
    plot!(
        p2,
        sol.t,
        sol.extras.surface_mass_flux .* 1000,
        lw = 2.2,
        color = :crimson,
        label = "Surface MLR",
    )
    plot!(
        p2,
        sol.t,
        sol.extras.gas_generation_rate .* 1000,
        lw = 1.5,
        color = :darkorange,
        linestyle = :dash,
        label = "Gas generation",
    )

    # Panel 3: top + back temperatures
    p3 = plot(
        xlabel = "Time [s]",
        ylabel = "T [K]",
        title = "Surface & back temperature",
        legend = :bottomright,
    )
    plot!(p3, sol.t, sol.surface_T, lw = 2.2, color = :darkred, label = "Top")
    plot!(
        p3,
        sol.t,
        sol.extras.back_temperature,
        lw = 2.2,
        color = :steelblue,
        label = "Back",
    )

    # Panel 4: thickness + surface heat-flux components
    p4 = plot(
        xlabel = "Time [s]",
        ylabel = "q″ [kW/m²]",
        title = "Surface heat-flux balance",
        legend = :topright,
    )
    plot!(
        p4,
        sol.t,
        sol.extras.surface_heat_flux_absorbed ./ 1000,
        lw = 2.0,
        color = :crimson,
        label = "Absorbed (incident·ε)",
    )
    plot!(
        p4,
        sol.t,
        .-sol.extras.surface_heat_flux_reradiation ./ 1000,
        lw = 1.6,
        color = :gray60,
        linestyle = :dash,
        label = "Re-radiation (loss)",
    )
    plot!(
        p4,
        sol.t,
        .-sol.extras.surface_heat_flux_convection ./ 1000,
        lw = 1.6,
        color = :gray40,
        linestyle = :dot,
        label = "Convection (loss)",
    )
    plot!(
        p4,
        sol.t,
        sol.extras.surface_heat_flux_conduction ./ 1000,
        lw = 2.0,
        color = :navy,
        label = "Conduction in",
    )

    # Panel 5: composition mass fractions
    n = length(sol.t)
    wf_v = zeros(n)
    wf_i = zeros(n)
    wf_c = zeros(n)
    wf_a = zeros(n)
    wf_g = zeros(n)
    wf_w = zeros(n)
    for k = 1:n
        m_v = _per_component_mass(sol, (VIRGIN_IDX,), k)
        m_i = _per_component_mass(sol, INT_IDXS, k)
        m_c = _per_component_mass(sol, (CHAR_IDX,), k)
        m_a = _per_component_mass(sol, (ASH_IDX,), k)
        m_g = _per_component_mass(sol, GAS_IDXS, k)
        m_w = _per_component_mass(sol, (VAPOR_IDX,), k)
        denom = m_v + m_i + m_c + m_a + m_g + m_w
        denom = denom > 0 ? denom : 1.0
        wf_v[k] = m_v / denom
        wf_i[k] = m_i / denom
        wf_c[k] = m_c / denom
        wf_a[k] = m_a / denom
        wf_g[k] = m_g / denom
        wf_w[k] = m_w / denom
    end
    p5 = plot(
        xlabel = "Time [s]",
        ylabel = "Mass fraction [-]",
        title = "Composition (in-domain)",
        legend = :outerright,
        ylims = (0, 1.05),
    )
    plot!(p5, sol.t, wf_v, lw = 2.0, color = :forestgreen, label = "Virgin")
    plot!(p5, sol.t, wf_i, lw = 2.0, color = :goldenrod, label = "Intermediates")
    plot!(p5, sol.t, wf_c, lw = 2.0, color = :black, label = "Char")
    plot!(p5, sol.t, wf_a, lw = 1.4, color = :gray60, label = "Ash")
    plot!(p5, sol.t, wf_g, lw = 1.6, color = :crimson, label = "Pyro gas (in-domain)")
    if maximum(wf_w) > 1e-6
        plot!(p5, sol.t, wf_w, lw = 1.4, color = :steelblue, label = "Water vapor")
    end

    # Panel 6: temperature profiles at selected times
    t_targets = [0.05, 0.15, 0.30, 0.50, 0.70, 0.90, 1.00] .* sol.t[end]
    p6 = plot(
        xlabel = "Depth from heated surface [mm]",
        ylabel = "T [K]",
        title = "Temperature profiles",
        legend = :outerright,
    )
    L0 = sol.extras.thickness_history[1]
    for t_t in t_targets
        idx = argmin(abs.(sol.t .- t_t))
        T_prof = sol.T[:, idx]
        L_t = sol.extras.thickness_history[idx]
        # The state is stored bottom→top after ALE; re-orient depth from heated surface
        n_c = length(T_prof)
        depth = collect(range(0.0, L_t * 1000, length = n_c))         # mm
        depth = reverse(depth)                                         # 0 at top
        T_plot = reverse(T_prof)
        plot!(p6, depth, T_plot, lw = 1.6, label = @sprintf("t=%.0fs", sol.t[idx]))
    end

    fig = plot(
        p1,
        p2,
        p3,
        p4,
        p5,
        p6;
        layout = (3, 2),
        size = (1400, 1100),
        left_margin = 8Plots.mm,
        bottom_margin = 6Plots.mm,
        plot_title = title_root,
    )
    out = joinpath(out_dir, "macfp_$(case_label)_detail.png")
    savefig(fig, out)
    return out
end

# ============================================================================
# Driver
# ============================================================================

function main()
    out_dir = joinpath(
        @__DIR__,
        "results",
        "macfp_" * Dates.format(Dates.now(), "yyyymmdd_HHMMSS"),
    )
    mkpath(out_dir)
    @info "MaCFP gasification run" out_dir threads=Threads.nthreads() θ=THETA_STAGE2 stage=STAGE_USED

    # Spawn all six in parallel when threading is available.
    tasks = Dict{String,Task}()
    for (lbl, q_kW, th_mm, t_final) in CASES
        q_inc = q_kW * 1000.0
        th_m = th_mm / 1000.0
        @info "  scheduling" case=lbl q_inc_kW=q_kW thickness_mm=th_mm t_final_s=t_final
        tasks[lbl] =
            Threads.@spawn build_and_solve_macfp(THETA_STAGE2, q_inc, th_m, t_final)
    end

    sols = Dict{String,Any}()
    for (lbl, q_kW, th_mm, _) in CASES
        t0 = time()
        sol = try
            fetch(tasks[lbl])
        catch err
            @warn "  solve threw" case=lbl exception=err
            nothing
        end
        wall = round(time() - t0, digits = 1)
        if sol === nothing
            @warn "  failed" case=lbl wall_s=wall
        else
            ok = (
                sol.retcode == :Success ||
                sol.retcode == ReturnCode.Success ||
                sol.retcode == :Terminated
            )
            @info "  done" case=lbl retcode=sol.retcode wall_s=wall n_times=length(sol.t) ok=ok
            csv_path = joinpath(out_dir, "macfp_$(lbl).csv")
            write_case_csv(csv_path, sol)
            png_path = detail_plot(lbl, q_kW, th_mm, sol, out_dir)
            @info "  wrote" csv=basename(csv_path) detail_png=basename(png_path)
            sols[lbl] = sol
        end
    end

    if !isempty(sols)
        png_path = summary_plot(sols, out_dir)
        @info "Summary plot" png=png_path
    end

    @info "All done" out_dir
    return out_dir
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
