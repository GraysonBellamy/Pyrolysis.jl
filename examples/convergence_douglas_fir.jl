"""
Grid Convergence Study for Douglas Fir Cone Calorimeter Simulation

Runs the 07_douglas_fir_cone simulation at varying spatial resolutions to
determine how coarse the mesh can be without meaningfully impacting results.

Grid sizes tested (finest = reference solution):
  480, 240, 120, 60, 30, 15 cells  →  0.025 to 0.8 mm cell size

Outputs:
  - Convergence table: scalar QoIs vs grid size (peak MLR, time to ignition, etc.)
  - MLR overlay plot: all resolutions on one figure
  - Relative error plot: % error vs cell count for each QoI
  - Surface and back temperature overlay plots
  - Solve time vs grid size (speedup chart)

Usage:
  julia convergence_douglas_fir.jl                  # sequential (1 worker)
  julia -p auto convergence_douglas_fir.jl          # parallel (auto-detect cores)
  julia -p 4 convergence_douglas_fir.jl             # parallel (4 workers)
"""

using Distributed

# Load packages on all workers
@everywhere begin
    using Pyrolysis
    using LinearSolve
    using Sparspak
    using OrdinaryDiffEq
    using OrdinaryDiffEqSDIRK: KenCarp4   # v7 split: KenCarp4 not re-exported by the umbrella pkg
    using LinearAlgebra
    using LinearSolve: SparspakFactorization
end

# Main process only
using Statistics
using Printf
using Plots
using DelimitedFiles

# ============================================================================
# Grid sizes to test (finest first = reference)
# ============================================================================

const GRID_SIZES = [480, 240, 120, 60, 30, 15]
const RESULTS_DIR = joinpath(@__DIR__, "convergence_results")

# ============================================================================
# Build and solve at a given grid resolution (available on all workers)
# ============================================================================

@everywhere function build_and_solve(n_cells::Int)
    # --- Temperature-dependent heat capacities (nominal from sensitivity analysis) ---
    # Virgin DF: piecewise linear (3 segments) from HFM (low-T) + STA/DSC (high-T)
    df_cp_a1 = 4.133836600273531   # slope, T ≤ 323K (HFM)
    df_cp_b1 = -18.01630709338176  # intercept, T ≤ 323K
    df_cp_a3 = 9.136799235573225   # slope, T > 425K (STA/DSC)
    df_cp_b3 = -2139.3044552435705 # intercept, T > 425K
    T_HFM_MAX = 323.16
    T_STA_MIN = 425.02219679025956
    T_INTERP_WIDTH = T_STA_MIN - T_HFM_MAX
    df_cp = CallableProperty(
        T -> begin
            if T <= T_HFM_MAX
                df_cp_b1 + df_cp_a1 * T
            elseif T <= T_STA_MIN
                cp_low = df_cp_b1 + df_cp_a1 * T_HFM_MAX
                cp_high = df_cp_b3 + df_cp_a3 * T_STA_MIN
                cp_low + (cp_high - cp_low) * (T - T_HFM_MAX) / T_INTERP_WIDTH
            else
                df_cp_b3 + df_cp_a3 * T
            end
        end,
    )

    # Char cp (linear, from independent char STA measurement)
    char_cp_a = 1.8372509265710314
    char_cp_b = 312.4273785167445
    char_cp_func = CallableProperty(T -> char_cp_b + char_cp_a * T)

    # Intermediate cp: cumulative mass-loss interpolation between virgin and char
    # α = (1 - cumulative_solid_yield) / (1 - final_char_yield)
    # final_char_yield = 0.95 × 0.80 × 0.37 × 0.79 × 0.74
    final_char_yield = 0.95 * 0.80 * 0.37 * 0.79 * 0.74
    α_int1 = (1 - 0.95) / (1 - final_char_yield)                      # ≈ 0.060
    α_int2 = (1 - 0.95 * 0.80) / (1 - final_char_yield)              # ≈ 0.287
    α_int3 = (1 - 0.95 * 0.80 * 0.37) / (1 - final_char_yield)      # ≈ 0.860
    α_int4 = (1 - 0.95 * 0.80 * 0.37 * 0.79) / (1 - final_char_yield)  # ≈ 0.931
    df_int1_cp = CallableProperty(T -> (1 - α_int1) * df_cp(T) + α_int1 * char_cp_func(T))
    df_int2_cp = CallableProperty(T -> (1 - α_int2) * df_cp(T) + α_int2 * char_cp_func(T))
    df_int3_cp = CallableProperty(T -> (1 - α_int3) * df_cp(T) + α_int3 * char_cp_func(T))
    df_int4_cp = CallableProperty(T -> (1 - α_int4) * df_cp(T) + α_int4 * char_cp_func(T))

    # Moisture cp: Shomate equation (Chase, 1998)
    moisture_cp = CallableProperty(T -> begin
        t = T / 1000.0
        M = 0.01801528
        A = -203.6060
        B = 1523.290
        C = -3196.413
        D = 2474.455
        E = 3.855326
        (A + B * t + C * t^2 + D * t^3 + E / t^2) / M
    end)

    # Water vapor cp (NIST Shomate fit, matches sensitivity analysis)
    water_vapor_cp = CallableProperty(T -> 1716.7 + 0.4154 * T + 1.530e-4 * T^2)

    # Heat of sorption (state-dependent) — nominal CR2 fit from sensitivity analysis
    moisture_idx = 7
    ρ_dry_basis = 530.0
    hos_a = 2450.420909987441    # kJ/kg - asymptotic baseline
    hos_b = 875.6441221735336    # kJ/kg - exponential amplitude
    hos_c = 0.1121684273316275   # 1/% - exponential decay rate
    hos = StateDependentProperty(
        (T, ξ) -> begin
            MC = 100.0 * ξ[moisture_idx] / ρ_dry_basis
            return (hos_a + hos_b * exp(-hos_c * MC)) * 1e3
        end,
    )

    # Gas transfer and permeability (physically grounded values)
    # λ_solid: D_free × ε/τ; λ_gas: free molecular diffusion at ~500K
    # κ: intrinsic matrix permeability for Darcy flow
    λ_wood = 5.0e-7     # m²/s - virgin wood (ε~0.05, τ~5)
    λ_char = 1.0e-5     # m²/s - porous char (ε~0.7, τ~2)
    λ_gas = 3.0e-5     # m²/s - reference free diffusion (M~29, ~500K)
    κ_virgin = 1.0e-15  # m² - tight virgin wood
    κ_char = 1.0e-12    # m² - open char pore network

    # --- Material definition ---
    # Thermodynamic properties (cp, h_rxn) use nominal values from sensitivity analysis
    material = Material(
        name = :DouglasFir,
        components = (
            SolidComponent(
                :DF,
                ρ = 530.0,
                c = df_cp,
                k = 0.145,
                ε = 0.81,
                α = 0.0,
                λ = λ_wood,
                κ = κ_virgin,
            ),
            SolidComponent(
                :DF_int1,
                ρ = 507.2,
                c = df_int1_cp,
                k = 0.292951,
                ε = 0.790965,
                α = 0.0,
                λ = λ_wood * 2.0,
                κ = κ_virgin * 3.0,
            ),
            SolidComponent(
                :DF_int2,
                ρ = 419.8,
                c = df_int2_cp,
                k = 0.329938,
                ε = 0.786207,
                α = 0.0,
                λ = λ_wood * 5.0,
                κ = κ_virgin * 10.0,
            ),
            SolidComponent(
                :DF_int3,
                ρ = 203.2,
                c = df_int3_cp,
                k = 0.644833,
                ε = 0.745694,
                α = 0.0,
                λ = λ_char * 0.3,
                κ = κ_char * 0.01,
            ),
            SolidComponent(
                :DF_int4,
                ρ = 176.6,
                c = df_int4_cp,
                k = 0.777700,
                ε = 0.728600,
                α = 0.0,
                λ = λ_char * 0.7,
                κ = κ_char * 0.1,
            ),
            SolidComponent(
                :Char,
                ρ = 150.0,
                c = char_cp_func,
                k = 1.0,
                ε = 0.70,
                α = 0.0,
                λ = λ_char,
                κ = κ_char,
            ),
            LiquidComponent(
                :MOISTURE,
                ρ = 1e6,
                c = moisture_cp,
                k = 0.145,
                ε = 0.81,
                α = 0.0,
                κ = 0.0,
            ),
            GaseousComponent(
                :DF_gas1,
                M = 0.060,
                c = 1400.0,
                k = 0.020,
                λ = λ_gas * sqrt(0.029/0.060),
            ),
            GaseousComponent(
                :DF_gas2,
                M = 0.070,
                c = 1400.0,
                k = 0.020,
                λ = λ_gas * sqrt(0.029/0.070),
            ),
            GaseousComponent(
                :DF_gas3,
                M = 0.110,
                c = 1300.0,
                k = 0.018,
                λ = λ_gas * sqrt(0.029/0.110),
            ),
            GaseousComponent(
                :DF_gas4,
                M = 0.070,
                c = 1500.0,
                k = 0.020,
                λ = λ_gas * sqrt(0.029/0.070),
            ),
            GaseousComponent(
                :DF_gas5,
                M = 0.050,
                c = 1600.0,
                k = 0.025,
                λ = λ_gas * sqrt(0.029/0.050),
            ),
            GaseousComponent(
                :WATER_VAPOR,
                M = 0.018,
                c = water_vapor_cp,
                k = 0.025,
                λ = λ_gas * sqrt(0.029/0.018),
            ),
        ),
        reactions = (
            # Heats of reaction: nominal values from sensitivity analysis calibration
            Reaction(
                :pyrolysis_step1,
                1 => (2, 8),
                A = 74816.95,
                E = 71.59e3,
                h = -19200.66288032251,
                n = 1.0,
                yields = (0.95, 0.05),
            ),
            Reaction(
                :pyrolysis_step2,
                2 => (3, 9),
                A = 23550492839.0,
                E = 140.9e3,
                h = -5263.201466602739,
                n = 1.0,
                yields = (0.80, 0.20),
            ),
            Reaction(
                :pyrolysis_step3,
                3 => (4, 10),
                A = 7.744618e12,
                E = 181.0e3,
                h = 264424.79037373835,
                n = 1.0,
                yields = (0.37, 0.63),
            ),
            Reaction(
                :pyrolysis_step4,
                4 => (5, 11),
                A = 103276.14,
                E = 98.94e3,
                h = -217151.69043762444,
                n = 1.0,
                yields = (0.79, 0.21),
            ),
            Reaction(
                :pyrolysis_step5,
                5 => (6, 12),
                A = 6.397348,
                E = 56.01e3,
                h = -326.40885469229025,
                n = 1.0,
                yields = (0.74, 0.26),
            ),
            Reaction(
                :moisture_evaporation,
                7 => 13,
                A = 2.997e5,
                E = 5.0e4,
                h = hos,
                n = 1.0,
            ),
        ),
        mixing = MaterialMixing(
            k = MixingSpec(WEIGHTED, 0.5),
            λ = MixingSpec(WEIGHTED, 0.5),
            κ = MixingSpec(WEIGHTED, 0.5),
        ),
    )

    # --- Problem setup ---
    sample_thickness = 12.0e-3
    T_ambient = 300.0
    q_external = 50e3
    h_conv = 10.0
    moisture_content_db = 0.06
    ρ_DF_initial = 530.0
    ρ_moisture_initial = moisture_content_db * ρ_DF_initial

    mesh = create_uniform_mesh(
        sample_thickness,
        n_cells,
        Val(13);
        T_initial = T_ambient,
        ξ_initial = (
            ρ_DF_initial,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            ρ_moisture_initial,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
        ),
    )

    bc_set = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = q_external) +
                  RadiativeBC(T_inf = 300) +
                  ConvectiveBC(h = h_conv, T_inf = 300),
            bottom = AdiabaticBC(),
        ),
        mass = (
            top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
            bottom = ImpermeableBC(),
        ),
    )

    T_initial_func = z -> T_ambient
    ξ_initial_funcs = [
        z -> ρ_DF_initial,
        z -> 0.0,
        z -> 0.0,
        z -> 0.0,
        z -> 0.0,
        z -> 0.0,
        z -> ρ_moisture_initial,
        z -> 0.0,
        z -> 0.0,
        z -> 0.0,
        z -> 0.0,
        z -> 0.0,
        z -> 0.0,
    ]

    problem = PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        T_initial = T_initial_func,
        ξ_initial = ξ_initial_funcs,
        tspan = (0.0, 600.0),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, 600.0, length = 3601)),
    )

    # --- Solve ---
    solution = solve(
        problem;
        abstol = 1e-9,
        reltol = 1e-7,
        dtmin = 1e-14,
        progress = false,
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        min_thickness = 0.5e-3,
        handle_depletion = false,
        # Default Jacobian under ALE is `Structured(strategy = ApproxSolve())`.
        integrator = KenCarp4(linsolve = SparspakFactorization()),
    )

    return solution
end

# ============================================================================
# QoI extraction
# ============================================================================

const MLR_IGNITION_THRESHOLD = 3e-3  # kg/(m²·s)

struct ConvergenceResult
    n_cells::Int
    cell_size_mm::Float64
    peak_mlr::Float64
    time_to_peak_mlr::Float64
    time_to_ignition::Float64
    total_mass_loss_pct::Float64
    peak_surface_temp::Float64
    peak_back_temp::Float64
    final_thickness_mm::Float64
    solve_time::Float64
    retcode::Symbol
    times::Vector{Float64}
    mlr_curve::Vector{Float64}
    surface_temp_curve::Vector{Float64}
    back_temp_curve::Vector{Float64}
end

function extract_result(sol, n_cells)
    t = sol.t
    mlr = sol.extras.surface_mass_flux
    T_surf = sol.surface_T
    T_back = sol.extras.back_temperature

    peak_mlr_val = maximum(mlr)
    peak_idx = argmax(mlr)
    time_to_peak = t[peak_idx]

    tig_idx = findfirst(x -> x > MLR_IGNITION_THRESHOLD, mlr)
    time_to_ignition = tig_idx !== nothing ? t[tig_idx] : Inf

    total_mass_loss = (1.0 - sol.extras.total_mass[end] / sol.extras.total_mass[1]) * 100.0

    sample_thickness = 12.0e-3
    return ConvergenceResult(
        n_cells,
        sample_thickness / n_cells * 1000,
        peak_mlr_val,
        time_to_peak,
        time_to_ignition,
        total_mass_loss,
        maximum(T_surf),
        maximum(T_back),
        sol.extras.thickness_history[end] * 1000,
        sol.solve_time,
        sol.retcode,
        t,
        mlr,
        T_surf,
        T_back,
    )
end

# ============================================================================
# Save results to CSV
# ============================================================================

"""
Write two kinds of CSV files to `dir`:
  convergence_scalars.csv        – one row per grid size, all scalar QoIs
  convergence_curves_<N>cells.csv – time-series (time, mlr, T_surf, T_back) per run
"""
function save_results(results::Vector{ConvergenceResult}, dir::String)
    # --- Scalar QoI summary ---
    scalar_path = joinpath(dir, "convergence_scalars.csv")
    open(scalar_path, "w") do io
        println(
            io,
            "n_cells,cell_size_mm,peak_mlr_g_m2s,time_to_ignition_s," *
            "time_to_peak_mlr_s,total_mass_loss_pct,peak_surface_temp_K," *
            "peak_back_temp_K,final_thickness_mm,solve_time_s,retcode",
        )
        for r in results
            tig =
                isfinite(r.time_to_ignition) ? @sprintf("%.4f", r.time_to_ignition) : "Inf"
            println(
                io,
                @sprintf(
                    "%d,%.6f,%.6f,%s,%.4f,%.6f,%.4f,%.4f,%.6f,%.4f,%s",
                    r.n_cells,
                    r.cell_size_mm,
                    r.peak_mlr * 1000,
                    tig,
                    r.time_to_peak_mlr,
                    r.total_mass_loss_pct,
                    r.peak_surface_temp,
                    r.peak_back_temp,
                    r.final_thickness_mm,
                    r.solve_time,
                    string(r.retcode)
                )
            )
        end
    end
    println("  Saved convergence_scalars.csv")

    # --- Per-run time-series curves ---
    for r in results
        curve_path = joinpath(dir, "convergence_curves_$(r.n_cells)cells.csv")
        open(curve_path, "w") do io
            println(io, "time_s,mlr_g_m2s,surface_temp_K,back_temp_K")
            for j in eachindex(r.times)
                println(
                    io,
                    @sprintf(
                        "%.6f,%.8f,%.6f,%.6f",
                        r.times[j],
                        r.mlr_curve[j] * 1000,
                        r.surface_temp_curve[j],
                        r.back_temp_curve[j]
                    )
                )
            end
        end
    end
    println("  Saved convergence_curves_<N>cells.csv for each grid size")
end

# ============================================================================
# Main
# ============================================================================

function main()
    n_workers = nworkers()
    println("\n" * "="^70)
    println("  Grid Convergence Study — Douglas Fir Cone Calorimeter")
    println("="^70)
    println("  Grid sizes: ", join(GRID_SIZES, ", "), " cells")
    println("  Sample thickness: 12 mm")
    println(
        "  Cell sizes: ",
        join([@sprintf("%.3f", 12.0 / n) for n in GRID_SIZES], ", "),
        " mm",
    )
    println("  Workers: $n_workers (use `julia -p N` to run in parallel)")
    println("="^70 * "\n")

    mkpath(RESULTS_DIR)

    if n_workers > 1
        # --- Parallel execution via pmap ---
        println(
            "  Running $(length(GRID_SIZES)) grid sizes in parallel on $n_workers workers...",
        )
        t_wall_start = time()

        solutions = pmap(build_and_solve, GRID_SIZES)

        t_wall = time() - t_wall_start
        println(@sprintf("  All solves completed in %.1f s wall time\n", t_wall))

        # Extract results in original order
        results = [extract_result(sol, n) for (sol, n) in zip(solutions, GRID_SIZES)]
        for r in results
            println(
                @sprintf(
                    "    n_cells = %3d  [%s]  peak MLR = %.4f g/(m²·s)  solve = %.1f s",
                    r.n_cells,
                    r.retcode,
                    r.peak_mlr * 1000,
                    r.solve_time
                )
            )
        end
    else
        # --- Sequential execution ---
        results = ConvergenceResult[]
        for (i, n) in enumerate(GRID_SIZES)
            cell_mm = 12.0 / n
            println(
                @sprintf(
                    "  [%d/%d] Running n_cells = %d (Δx = %.3f mm)...",
                    i,
                    length(GRID_SIZES),
                    n,
                    cell_mm
                )
            )

            t_start = time()
            sol = build_and_solve(n)
            elapsed = time() - t_start

            r = extract_result(sol, n)
            push!(results, r)

            println(
                @sprintf(
                    "         Done in %.1f s  [%s]  peak MLR = %.4f g/(m²·s)",
                    elapsed,
                    r.retcode,
                    r.peak_mlr * 1000
                )
            )
        end
    end

    # Reference = finest grid
    ref = results[1]

    # ========================================================================
    # Save results to CSV
    # ========================================================================

    println("\nSaving results...")
    save_results(results, RESULTS_DIR)

    # ========================================================================
    # Convergence Table
    # ========================================================================

    println("\n" * "="^70)
    println("  Convergence Results (reference: $(ref.n_cells) cells)")
    println("="^70)

    # Header
    println(
        @sprintf(
            "\n  %-8s %-8s  %-12s %-12s %-12s %-12s %-12s %-10s",
            "n_cells",
            "Δx [mm]",
            "Peak MLR",
            "t_ign [s]",
            "t_peak [s]",
            "Mass Loss",
            "T_surf,max",
            "Solve [s]"
        )
    )
    println(
        @sprintf(
            "  %-8s %-8s  %-12s %-12s %-12s %-12s %-12s %-10s",
            "",
            "",
            "[g/(m²·s)]",
            "",
            "",
            "[%]",
            "[K]",
            ""
        )
    )
    println("  " * "-"^96)

    for r in results
        println(
            @sprintf(
                "  %-8d %-8.3f  %-12.4f %-12.2f %-12.2f %-12.2f %-12.1f %-10.1f",
                r.n_cells,
                r.cell_size_mm,
                r.peak_mlr * 1000,
                r.time_to_ignition,
                r.time_to_peak_mlr,
                r.total_mass_loss_pct,
                r.peak_surface_temp,
                r.solve_time
            )
        )
    end

    # Relative errors
    println(
        @sprintf(
            "\n  %-8s %-8s  %-12s %-12s %-12s %-12s %-12s %-10s",
            "n_cells",
            "Δx [mm]",
            "Δ Peak MLR",
            "Δ t_ign",
            "Δ t_peak",
            "Δ Mass Loss",
            "Δ T_surf",
            "Speedup"
        )
    )
    println("  " * "-"^96)

    for r in results
        Δ_mlr = (r.peak_mlr - ref.peak_mlr) / ref.peak_mlr * 100
        Δ_tig = if isfinite(r.time_to_ignition) && isfinite(ref.time_to_ignition)
            (r.time_to_ignition - ref.time_to_ignition) / ref.time_to_ignition * 100
        else
            NaN
        end
        Δ_tp = (r.time_to_peak_mlr - ref.time_to_peak_mlr) / ref.time_to_peak_mlr * 100
        Δ_ml =
            (r.total_mass_loss_pct - ref.total_mass_loss_pct) / ref.total_mass_loss_pct *
            100
        Δ_ts = (r.peak_surface_temp - ref.peak_surface_temp) / ref.peak_surface_temp * 100
        speedup = ref.solve_time / r.solve_time

        println(
            @sprintf(
                "  %-8d %-8.3f  %+11.2f%% %+11.2f%% %+11.2f%% %+11.2f%% %+11.2f%% %9.1fx",
                r.n_cells,
                r.cell_size_mm,
                Δ_mlr,
                Δ_tig,
                Δ_tp,
                Δ_ml,
                Δ_ts,
                speedup
            )
        )
    end

    # ========================================================================
    # Compute curve-level L2 errors
    # ========================================================================

    println("\n  Curve-level L2 relative errors vs $(ref.n_cells)-cell reference:")
    println("  " * "-"^60)
    println(
        @sprintf(
            "  %-8s  %-15s %-15s %-15s",
            "n_cells",
            "MLR L2 [%]",
            "T_surf L2 [%]",
            "T_back L2 [%]"
        )
    )
    println("  " * "-"^60)

    for r in results
        # L2 relative error: ||f - f_ref||₂ / ||f_ref||₂
        mlr_l2 = norm(r.mlr_curve .- ref.mlr_curve) / norm(ref.mlr_curve) * 100
        tsurf_l2 =
            norm(r.surface_temp_curve .- ref.surface_temp_curve) /
            norm(ref.surface_temp_curve) * 100
        tback_l2 =
            norm(r.back_temp_curve .- ref.back_temp_curve) / norm(ref.back_temp_curve) * 100

        println(
            @sprintf(
                "  %-8d  %-15.3f %-15.3f %-15.3f",
                r.n_cells,
                mlr_l2,
                tsurf_l2,
                tback_l2
            )
        )
    end

    # ========================================================================
    # Plots
    # ========================================================================

    println("\nGenerating plots...")

    colors = [:black, :blue, :green, :orange, :red, :purple]

    # --- 1. MLR overlay ---
    plt_mlr = plot(
        xlabel = "Time [s]",
        ylabel = "MLR [g/(m²·s)]",
        title = "Mass Loss Rate — Grid Convergence",
        legend = :topright,
        size = (900, 500),
    )
    for (i, r) in enumerate(results)
        lw = r.n_cells == ref.n_cells ? 2.5 : 1.5
        ls = r.n_cells == ref.n_cells ? :solid : :dash
        plot!(
            r.times,
            r.mlr_curve .* 1000,
            label = "$(r.n_cells) cells ($(r.cell_size_mm) mm)",
            lw = lw,
            ls = ls,
            color = colors[i],
        )
    end
    savefig(plt_mlr, joinpath(RESULTS_DIR, "convergence_mlr.png"))
    println("  Saved convergence_mlr.png")

    # --- 2. Surface temperature overlay ---
    plt_tsurf = plot(
        xlabel = "Time [s]",
        ylabel = "Surface Temperature [K]",
        title = "Surface Temperature — Grid Convergence",
        legend = :bottomright,
        size = (900, 500),
    )
    for (i, r) in enumerate(results)
        lw = r.n_cells == ref.n_cells ? 2.5 : 1.5
        ls = r.n_cells == ref.n_cells ? :solid : :dash
        plot!(
            r.times,
            r.surface_temp_curve,
            label = "$(r.n_cells) cells",
            lw = lw,
            ls = ls,
            color = colors[i],
        )
    end
    savefig(plt_tsurf, joinpath(RESULTS_DIR, "convergence_surface_temp.png"))
    println("  Saved convergence_surface_temp.png")

    # --- 3. Back temperature overlay ---
    plt_tback = plot(
        xlabel = "Time [s]",
        ylabel = "Back Temperature [K]",
        title = "Back Temperature — Grid Convergence",
        legend = :bottomright,
        size = (900, 500),
    )
    for (i, r) in enumerate(results)
        lw = r.n_cells == ref.n_cells ? 2.5 : 1.5
        ls = r.n_cells == ref.n_cells ? :solid : :dash
        plot!(
            r.times,
            r.back_temp_curve,
            label = "$(r.n_cells) cells",
            lw = lw,
            ls = ls,
            color = colors[i],
        )
    end
    savefig(plt_tback, joinpath(RESULTS_DIR, "convergence_back_temp.png"))
    println("  Saved convergence_back_temp.png")

    # --- 4. Relative error vs grid size ---
    # Helper: replace non-positive values with NaN so log10 scale doesn't blow up
    pos_or_nan(x) = x > 0.0 ? x : NaN

    n_vec = [r.n_cells for r in results]
    Δ_peak_mlr =
        [pos_or_nan(abs(r.peak_mlr - ref.peak_mlr) / ref.peak_mlr * 100) for r in results]
    Δ_tig = [
        isfinite(r.time_to_ignition) && isfinite(ref.time_to_ignition) ?
        pos_or_nan(
            abs(r.time_to_ignition - ref.time_to_ignition) / ref.time_to_ignition * 100,
        ) : NaN for r in results
    ]
    Δ_mass = [
        pos_or_nan(
            abs(r.total_mass_loss_pct - ref.total_mass_loss_pct) / ref.total_mass_loss_pct * 100,
        ) for r in results
    ]
    Δ_tsurf = [
        pos_or_nan(
            abs(r.peak_surface_temp - ref.peak_surface_temp) / ref.peak_surface_temp * 100,
        ) for r in results
    ]

    # Compute y-axis limits from finite data so log scale always has valid ticks
    all_errors = filter(isfinite, vcat(Δ_peak_mlr, Δ_tig, Δ_mass, Δ_tsurf))
    y_min = isempty(all_errors) ? 1e-4 : max(minimum(all_errors) / 10, 1e-6)
    y_max = isempty(all_errors) ? 100.0 : maximum(all_errors) * 10

    plt_err = plot(
        xlabel = "Number of Cells",
        ylabel = "Relative Error [%]",
        title = "Scalar QoI Convergence",
        xscale = :log10,
        yscale = :log10,
        legend = :topright,
        size = (900, 500),
        xlims = (10, 600),
        ylims = (y_min, y_max),
    )
    # Skip reference (error = 0 by definition)
    idx = 2:length(results)
    plot!(n_vec[idx], Δ_peak_mlr[idx], marker = :circle, lw = 2, label = "Peak MLR")
    if any(isfinite, Δ_tig[idx])
        plot!(n_vec[idx], Δ_tig[idx], marker = :square, lw = 2, label = "Time to Ignition")
    end
    plot!(n_vec[idx], Δ_mass[idx], marker = :diamond, lw = 2, label = "Total Mass Loss")
    plot!(
        n_vec[idx],
        Δ_tsurf[idx],
        marker = :utriangle,
        lw = 2,
        label = "Peak Surface Temp",
    )
    # 1% threshold line
    hline!([1.0], color = :gray, ls = :dash, lw = 1, label = "1% error")
    savefig(plt_err, joinpath(RESULTS_DIR, "convergence_error.png"))
    println("  Saved convergence_error.png")

    # --- 5. Solve time / speedup ---
    solve_times = [r.solve_time for r in results]
    speedups = [ref.solve_time / r.solve_time for r in results]

    plt_time = plot(
        layout = (1, 2),
        size = (1000, 400),
        title = ["Solve Time vs Grid Size" "Speedup vs Reference ($(ref.n_cells) cells)"],
    )
    bar!(
        plt_time[1],
        string.(n_vec),
        solve_times,
        xlabel = "Number of Cells",
        ylabel = "Solve Time [s]",
        legend = false,
        color = :steelblue,
        alpha = 0.8,
    )
    bar!(
        plt_time[2],
        string.(n_vec),
        speedups,
        xlabel = "Number of Cells",
        ylabel = "Speedup Factor",
        legend = false,
        color = :coral,
        alpha = 0.8,
    )
    hline!(plt_time[2], [1.0], color = :black, ls = :dash, lw = 1)
    savefig(plt_time, joinpath(RESULTS_DIR, "convergence_timing.png"))
    println("  Saved convergence_timing.png")

    # ========================================================================
    # Summary
    # ========================================================================

    println("\n" * "="^70)
    println("  Grid Convergence Study Complete!")
    println("  Results in: $RESULTS_DIR")
    println("="^70)
end

# Run
main()
