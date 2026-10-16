"""
Moisture Heat-of-Sorption Model Comparison

Compares two HoS models for moisture evaporation in Douglas Fir cone calorimeter
simulation across four initial moisture contents (3%, 6%, 9%, 12% dry basis):

1. Exponential decay (state-dependent):
   HoS(MC) = (a + b·exp(-c·MC)) × 1e3 J/kg
   where MC is instantaneous moisture content (% dry basis).
   Higher energy cost as moisture depletes — physically represents increasingly
   tightly bound water.

2. Constant average:
   HoS_avg = [∫₀^MC₀ (a + b·exp(-c·MC)) dMC] / MC₀ × 1e3 J/kg
   A single constant computed as the analytical integral average of the exponential
   model from 0 to MC_initial. Same total energy but uniform distribution in time.

The exponential model applies lower HoS at high MC (easy to release) and higher HoS
near depletion (tightly bound). The constant model applies a uniform penalty.

Output:
  - Per-run CSVs with full time-series diagnostics
  - Summary CSV with scalar metrics and percent differences
  - Comparison plots: HoS curves, MLR, temperatures, mass, differences
"""

using Pyrolysis
using Plots
using LinearSolve
using Sparspak
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK: KenCarp4   # v7 split: KenCarp4 not re-exported by the umbrella pkg
using Statistics
using LinearSolve: SparspakFactorization
using Printf
using DelimitedFiles

# ============================================================================
# Output Directory
# ============================================================================

const RESULTS_DIR = joinpath(@__DIR__, "moisture_hos_comparison_results")
mkpath(RESULTS_DIR)

# ============================================================================
# HoS Model Parameters (from sensitivity analysis nominal values)
# ============================================================================

const HOS_A = 2450.420909987441      # kJ/kg — asymptote
const HOS_B = 875.6441221735336      # kJ/kg — amplitude
const HOS_C = 0.1121684273316275     # 1/%   — decay rate

# MC levels to test (percent, dry basis)
const MC_LEVELS = [3.0, 6.0, 9.0, 12.0]

# ============================================================================
# Compute Constant Average HoS
# ============================================================================

"""
Analytical integral average of HoS(MC) = a + b·exp(-c·MC) over [0, MC₀]:

    avg = [a·MC₀ + (b/c)·(1 - exp(-c·MC₀))] / MC₀
"""
function compute_average_hos_kJ(a, b, c, MC₀)
    integral = a * MC₀ + (b / c) * (1.0 - exp(-c * MC₀))
    return integral / MC₀  # kJ/kg
end

# Print the HoS values for reference
println("Heat of Sorption Model Comparison")
println("="^60)
println("\nExponential model: HoS(MC) = $(HOS_A) + $(HOS_B)·exp(-$(HOS_C)·MC) [kJ/kg]")
println("\nConstant average values at each MC level:")
for mc in MC_LEVELS
    avg = compute_average_hos_kJ(HOS_A, HOS_B, HOS_C, mc)
    hos_at_mc = HOS_A + HOS_B * exp(-HOS_C * mc)
    hos_at_0 = HOS_A + HOS_B
    println(
        "  MC=$(Int(mc))%: avg = $(round(avg, digits=2)) kJ/kg " *
        "(range: $(round(hos_at_mc, digits=1)) → $(round(hos_at_0, digits=1)) kJ/kg)",
    )
end

# ============================================================================
# Heat Capacity Definitions (shared across all runs)
# ============================================================================

# Virgin DF: piecewise linear from HFM (low-T) and STA (high-T) regressions
const T_HFM_MAX = 323.16
const T_STA_MIN = 425.02219679025956
const T_INTERP_WIDTH = T_STA_MIN - T_HFM_MAX

const DF_CP_A1 = 4.133836600273531
const DF_CP_B1 = -18.01630709338176
const DF_CP_A3 = 9.136799235573225
const DF_CP_B3 = -2139.3044552435705

const CHAR_CP_A = 1.8372509265710314
const CHAR_CP_B = 312.4273785167445

function virgin_cp(T)
    if T <= T_HFM_MAX
        DF_CP_B1 + DF_CP_A1 * T
    elseif T <= T_STA_MIN
        cp_low = DF_CP_B1 + DF_CP_A1 * T_HFM_MAX
        cp_high = DF_CP_B3 + DF_CP_A3 * T_STA_MIN
        cp_low + (cp_high - cp_low) * (T - T_HFM_MAX) / T_INTERP_WIDTH
    else
        DF_CP_B3 + DF_CP_A3 * T
    end
end

# Cumulative-yield interpolation for intermediate cp's
const FINAL_CHAR_YIELD = 0.95 * 0.80 * 0.37 * 0.79 * 0.74
const INTERP_PROGRESS = [
    (1 - 0.95) / (1 - FINAL_CHAR_YIELD),
    (1 - 0.95 * 0.80) / (1 - FINAL_CHAR_YIELD),
    (1 - 0.95 * 0.80 * 0.37) / (1 - FINAL_CHAR_YIELD),
    (1 - 0.95 * 0.80 * 0.37 * 0.79) / (1 - FINAL_CHAR_YIELD),
]

# Gas transfer coefficients
const λ_WOOD = 5.0e-7
const λ_CHAR = 1.0e-5
const λ_GAS = 3.0e-5
const κ_VIRGIN = 1.0e-15
const κ_CHAR = 1.0e-12

# Component index for moisture (used in StateDependentProperty closure)
const MOISTURE_IDX = 7
const ρ_DRY_BASIS = 530.0

# ============================================================================
# Build-and-Solve Function
# ============================================================================

"""
	build_and_solve(mc_db, hos_model)

Build and solve a Douglas Fir cone calorimeter simulation.

Arguments:
- `mc_db::Float64`: moisture content as fraction on dry basis (e.g. 0.06 for 6%)
- `hos_model::Symbol`: `:exponential` or `:constant`

Returns the PyrolysisSolution.
"""
function build_and_solve(mc_db::Float64, hos_model::Symbol)
    # --- Heat capacities ---
    df_cp = CallableProperty(virgin_cp)
    char_cp_func = CallableProperty(T -> CHAR_CP_B + CHAR_CP_A * T)

    int1_cp = CallableProperty(
        T ->
            (1 - INTERP_PROGRESS[1]) * virgin_cp(T) +
            INTERP_PROGRESS[1] * (CHAR_CP_B + CHAR_CP_A * T),
    )
    int2_cp = CallableProperty(
        T ->
            (1 - INTERP_PROGRESS[2]) * virgin_cp(T) +
            INTERP_PROGRESS[2] * (CHAR_CP_B + CHAR_CP_A * T),
    )
    int3_cp = CallableProperty(
        T ->
            (1 - INTERP_PROGRESS[3]) * virgin_cp(T) +
            INTERP_PROGRESS[3] * (CHAR_CP_B + CHAR_CP_A * T),
    )
    int4_cp = CallableProperty(
        T ->
            (1 - INTERP_PROGRESS[4]) * virgin_cp(T) +
            INTERP_PROGRESS[4] * (CHAR_CP_B + CHAR_CP_A * T),
    )

    moisture_cp = CallableProperty(T -> begin
        t = T / 1000.0
        M = 0.01801528
        A = -203.6060
        B = 1523.290
        C = -3196.413
        D = 2474.455
        E = 3.855326
        (A + B*t + C*t^2 + D*t^3 + E/t^2) / M
    end)

    water_vapor_cp = CallableProperty(T -> 1716.7 + 0.4154 * T + 1.530e-4 * T^2)

    # --- Heat of sorption (the only difference between models) ---
    mc_pct = mc_db * 100.0
    if hos_model == :exponential
        h_evap = StateDependentProperty(
            (T, ξ) -> begin
                MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS
                return (HOS_A + HOS_B * exp(-HOS_C * MC)) * 1e3
            end,
        )
    elseif hos_model == :constant
        avg_kJ = compute_average_hos_kJ(HOS_A, HOS_B, HOS_C, mc_pct)
        h_evap = avg_kJ * 1e3  # Convert to J/kg — passed as Float64 → ConstantProperty
    else
        error("Unknown HoS model: $hos_model. Use :exponential or :constant.")
    end

    # --- Material ---
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
                λ = λ_WOOD,
                κ = κ_VIRGIN,
            ),
            SolidComponent(
                :DF_int1,
                ρ = 507.2,
                c = int1_cp,
                k = 0.292951,
                ε = 0.790965,
                α = 0.0,
                λ = λ_WOOD * 2.0,
                κ = κ_VIRGIN * 3.0,
            ),
            SolidComponent(
                :DF_int2,
                ρ = 419.8,
                c = int2_cp,
                k = 0.329938,
                ε = 0.786207,
                α = 0.0,
                λ = λ_WOOD * 5.0,
                κ = κ_VIRGIN * 10.0,
            ),
            SolidComponent(
                :DF_int3,
                ρ = 203.2,
                c = int3_cp,
                k = 0.644833,
                ε = 0.745694,
                α = 0.0,
                λ = λ_CHAR * 0.3,
                κ = κ_CHAR * 0.01,
            ),
            SolidComponent(
                :DF_int4,
                ρ = 176.6,
                c = int4_cp,
                k = 0.777700,
                ε = 0.728600,
                α = 0.0,
                λ = λ_CHAR * 0.7,
                κ = κ_CHAR * 0.1,
            ),
            SolidComponent(
                :Char,
                ρ = 150.0,
                c = char_cp_func,
                k = 1.0,
                ε = 0.70,
                α = 0.0,
                λ = λ_CHAR,
                κ = κ_CHAR,
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
                λ = λ_GAS * sqrt(0.029/0.060),
            ),
            GaseousComponent(
                :DF_gas2,
                M = 0.070,
                c = 1400.0,
                k = 0.020,
                λ = λ_GAS * sqrt(0.029/0.070),
            ),
            GaseousComponent(
                :DF_gas3,
                M = 0.110,
                c = 1300.0,
                k = 0.018,
                λ = λ_GAS * sqrt(0.029/0.110),
            ),
            GaseousComponent(
                :DF_gas4,
                M = 0.070,
                c = 1500.0,
                k = 0.020,
                λ = λ_GAS * sqrt(0.029/0.070),
            ),
            GaseousComponent(
                :DF_gas5,
                M = 0.050,
                c = 1600.0,
                k = 0.025,
                λ = λ_GAS * sqrt(0.029/0.050),
            ),
            GaseousComponent(
                :WATER_VAPOR,
                M = 0.018,
                c = water_vapor_cp,
                k = 0.025,
                λ = λ_GAS * sqrt(0.029/0.018),
            ),
        ),
        reactions = (
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
                h = h_evap,
                n = 1.0,
            ),
        ),
        mixing = MaterialMixing(
            k = MixingSpec(WEIGHTED, 0.5),
            λ = MixingSpec(WEIGHTED, 0.5),
            κ = MixingSpec(WEIGHTED, 0.5),
        ),
    )

    # --- Mesh and initial conditions ---
    sample_thickness = 12.0e-3
    n_cells = 240
    T_ambient = 300.0
    ρ_DF_initial = 530.0
    ρ_moisture_initial = mc_db * ρ_DF_initial

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
            top = RadiativeFluxBC(flux = 50000.0) +
                  RadiativeBC(T_inf = 300) +
                  ConvectiveBC(h = 10.0, T_inf = 300),
            bottom = AdiabaticBC(),
        ),
        mass = (
            top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
            bottom = ImpermeableBC(),
        ),
    )

    problem = PyrolysisProblem(
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        T_initial = z -> T_ambient,
        ξ_initial = [
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
        ],
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
# QoI Extraction
# ============================================================================

"""
Extract scalar metrics and time-series from a solution.
Returns a NamedTuple with all quantities needed for comparison.
"""
function extract_results(sol)
    times = sol.t
    mlr = sol.extras.surface_mass_flux           # kg/(m²·s)
    T_surface = sol.surface_T   # K
    T_back = sol.extras.back_temperature         # K
    total_mass = sol.extras.total_mass           # kg/m²
    thickness = sol.extras.thickness_history     # m

    # Scalar metrics
    peak_mlr = maximum(mlr) * 1000                        # g/(m²·s)
    peak_mlr_idx = argmax(mlr)
    time_to_peak_mlr = times[peak_mlr_idx]                # s

    mlr_critical = 3e-3  # kg/(m²·s)
    tig_idx = findfirst(mlr .> mlr_critical)
    time_to_ignition = tig_idx !== nothing ? times[tig_idx] : NaN  # s

    total_mass_loss_pct = (1 - total_mass[end] / total_mass[1]) * 100  # %
    final_thickness_mm = thickness[end] * 1000                          # mm
    peak_surface_temp = maximum(T_surface)                              # K

    # Back temperature at t=100s (indicative of drying front passage)
    idx_100 = argmin(abs.(times .- 100.0))
    T_back_100s = T_back[idx_100]

    # Average moisture concentration over time
    n_times = length(times)
    avg_moisture = zeros(n_times)
    for i = 1:n_times
        ξ_profile = sol.ξ[:, :, i]
        n_cells_i = size(ξ_profile, 2)
        if n_cells_i > 0
            avg_moisture[i] = sum(ξ_profile[MOISTURE_IDX, :]) / n_cells_i
        end
    end

    # Time when average moisture drops below 1% of initial
    moisture_initial = avg_moisture[1]
    if moisture_initial > 0
        depletion_idx = findfirst(avg_moisture .< 0.01 * moisture_initial)
        time_moisture_depleted = depletion_idx !== nothing ? times[depletion_idx] : NaN
    else
        time_moisture_depleted = 0.0
    end

    return (
        # Time series
        times = times,
        mlr = mlr,
        T_surface = T_surface,
        T_back = T_back,
        total_mass = total_mass,
        thickness = thickness,
        avg_moisture = avg_moisture,
        # Scalars
        peak_mlr = peak_mlr,
        time_to_peak_mlr = time_to_peak_mlr,
        time_to_ignition = time_to_ignition,
        total_mass_loss_pct = total_mass_loss_pct,
        final_thickness_mm = final_thickness_mm,
        peak_surface_temp = peak_surface_temp,
        T_back_100s = T_back_100s,
        time_moisture_depleted = time_moisture_depleted,
        solve_time = sol.solve_time,
        retcode = sol.retcode,
        # Full solution (for CSV output)
        gas_generation_rate = sol.extras.gas_generation_rate,
        q_absorbed = sol.extras.surface_heat_flux_absorbed,
        q_reradiation = sol.extras.surface_heat_flux_reradiation,
        q_convection = sol.extras.surface_heat_flux_convection,
        q_net = sol.extras.surface_heat_flux_net,
        surface_emissivity = sol.extras.surface_emissivity,
    )
end

# ============================================================================
# Write Per-Run CSV
# ============================================================================

function write_run_csv(filepath, res)
    open(filepath, "w") do io
        header = [
            "time",
            "mass",
            "mlr",
            "thickness",
            "front_temp",
            "back_temp",
            "avg_moisture",
            "gas_generation_rate",
            "q_absorbed",
            "q_reradiation",
            "q_convection",
            "q_net",
            "surface_emissivity",
        ]
        println(io, join(header, ","))

        for i in eachindex(res.times)
            data = [
                res.times[i],
                res.total_mass[i],
                res.mlr[i],
                res.thickness[i],
                res.T_surface[i],
                res.T_back[i],
                res.avg_moisture[i],
                res.gas_generation_rate[i],
                res.q_absorbed[i],
                res.q_reradiation[i],
                res.q_convection[i],
                res.q_net[i],
                res.surface_emissivity[i],
            ]
            println(io, join(data, ","))
        end
    end
end

# ============================================================================
# Run All Simulations
# ============================================================================

println("\n" * "="^60)
println("Running simulations (4 MC levels × 2 models = 8 total)")
println("="^60)

results = Dict{String,Any}()

for mc_pct in MC_LEVELS
    mc_db = mc_pct / 100.0
    for model in [:exponential, :constant]
        label = "$(model)_MC$(Int(mc_pct))"
        println("\n--- $label ---")
        @printf("  MC = %.0f%% dry basis, model = %s\n", mc_pct, model)

        t_start = time()
        sol = build_and_solve(mc_db, model)
        wall_time = time() - t_start

        res = extract_results(sol)
        results[label] = res

        # Save per-run CSV
        csv_path = joinpath(RESULTS_DIR, "timeseries_$(label).csv")
        write_run_csv(csv_path, res)

        @printf(
            "  Retcode: %s | Solve: %.1fs | Wall: %.1fs\n",
            res.retcode,
            res.solve_time,
            wall_time
        )
        @printf(
            "  Peak MLR: %.2f g/(m²·s) at t=%.1fs\n",
            res.peak_mlr,
            res.time_to_peak_mlr
        )
        @printf(
            "  Ignition (MLR>3): %.1fs | Moisture depleted: %.1fs\n",
            res.time_to_ignition,
            res.time_moisture_depleted
        )
    end
end

# ============================================================================
# Summary Table and CSV
# ============================================================================

println("\n" * "="^60)
println("Summary: Scalar Metrics")
println("="^60)

# Header
@printf(
    "\n%-6s  %-12s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
    "MC%",
    "Model",
    "PeakMLR",
    "t_peak",
    "t_ign",
    "ΔM%",
    "t_dry",
    "T_s_max",
    "T_b100",
    "thick_f"
)
@printf(
    "%-6s  %-12s  %8s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
    "",
    "",
    "g/m²s",
    "s",
    "s",
    "%",
    "s",
    "K",
    "K",
    "mm"
)
println("-"^100)

for mc_pct in MC_LEVELS
    for model in [:exponential, :constant]
        label = "$(model)_MC$(Int(mc_pct))"
        r = results[label]
        @printf(
            "%-6.0f  %-12s  %8.2f  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %8.1f  %8.2f\n",
            mc_pct,
            model,
            r.peak_mlr,
            r.time_to_peak_mlr,
            r.time_to_ignition,
            r.total_mass_loss_pct,
            r.time_moisture_depleted,
            r.peak_surface_temp,
            r.T_back_100s,
            r.final_thickness_mm
        )
    end
end

# Percent differences
println("\nPercent Difference (exponential vs constant):")
@printf(
    "%-6s  %8s  %8s  %8s  %8s  %8s  %8s  %8s\n",
    "MC%",
    "PeakMLR",
    "t_peak",
    "t_ign",
    "ΔM%",
    "t_dry",
    "T_s_max",
    "T_b100"
)
println("-"^74)

pct_diff(a, b) = abs(a) > 1e-10 ? (a - b) / abs(a) * 100.0 : 0.0

for mc_pct in MC_LEVELS
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]
    @printf(
        "%-6.0f  %+7.2f%%  %+7.2f%%  %+7.2f%%  %+7.2f%%  %+7.2f%%  %+7.2f%%  %+7.2f%%\n",
        mc_pct,
        pct_diff(re.peak_mlr, rc.peak_mlr),
        pct_diff(re.time_to_peak_mlr, rc.time_to_peak_mlr),
        pct_diff(re.time_to_ignition, rc.time_to_ignition),
        pct_diff(re.total_mass_loss_pct, rc.total_mass_loss_pct),
        pct_diff(re.time_moisture_depleted, rc.time_moisture_depleted),
        pct_diff(re.peak_surface_temp, rc.peak_surface_temp),
        pct_diff(re.T_back_100s, rc.T_back_100s)
    )
end

# Write summary CSV
summary_path = joinpath(RESULTS_DIR, "hos_comparison_summary.csv")
open(summary_path, "w") do io
    println(
        io,
        "mc_pct,model,peak_mlr_g_m2s,time_to_peak_mlr_s,time_to_ignition_s," *
        "total_mass_loss_pct,time_moisture_depleted_s,final_thickness_mm," *
        "peak_surface_temp_K,T_back_100s_K,solve_time_s",
    )
    for mc_pct in MC_LEVELS
        for model in [:exponential, :constant]
            label = "$(model)_MC$(Int(mc_pct))"
            r = results[label]
            println(
                io,
                join(
                    [
                        mc_pct,
                        model,
                        r.peak_mlr,
                        r.time_to_peak_mlr,
                        r.time_to_ignition,
                        r.total_mass_loss_pct,
                        r.time_moisture_depleted,
                        r.final_thickness_mm,
                        r.peak_surface_temp,
                        r.T_back_100s,
                        r.solve_time,
                    ],
                    ",",
                ),
            )
        end
    end
end
println("\nSummary saved to: $summary_path")

# ============================================================================
# Plot 1: HoS Model Curves
# ============================================================================

println("\nGenerating plots...")

mc_range = range(0.01, 15.0, length = 500)
hos_curve = [HOS_A + HOS_B * exp(-HOS_C * mc) for mc in mc_range]

p_hos = plot(
    mc_range,
    hos_curve,
    xlabel = "Moisture Content [% dry basis]",
    ylabel = "Heat of Sorption [kJ/kg]",
    title = "Heat of Sorption Models",
    label = "Exponential decay",
    linewidth = 2.5,
    color = :black,
    legend = :topright,
    ylims = (2400, 3400),
)

colors_mc = [:blue, :green, :orange, :red]
for (i, mc_pct) in enumerate(MC_LEVELS)
    avg = compute_average_hos_kJ(HOS_A, HOS_B, HOS_C, mc_pct)
    hline!(
        [avg],
        linestyle = :dash,
        color = colors_mc[i],
        linewidth = 1.5,
        label = "Avg, MC₀=$(Int(mc_pct))% ($(round(avg, digits=1)) kJ/kg)",
    )
    # Mark the MC range with a vertical line
    vline!([mc_pct], linestyle = :dot, color = colors_mc[i], linewidth = 1, label = false)
end

savefig(p_hos, joinpath(RESULTS_DIR, "hos_model_curves.png"))
println("  Saved: hos_model_curves.png")

# ============================================================================
# Plot 2: MLR Comparison (2×2 grid)
# ============================================================================

mlr_plots = []
for (i, mc_pct) in enumerate(MC_LEVELS)
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]

    p = plot(
        re.times,
        re.mlr .* 1000,
        xlabel = "Time [s]",
        ylabel = "MLR [g/(m²·s)]",
        title = "MC = $(Int(mc_pct))%",
        label = "Exponential",
        linewidth = 2,
        color = :black,
        legend = :topright,
    )
    plot!(
        rc.times,
        rc.mlr .* 1000,
        label = "Constant avg",
        linewidth = 2,
        color = colors_mc[i],
        linestyle = :dash,
    )
    hline!([3.0], linestyle = :dot, color = :gray, label = false, linewidth = 0.8)

    push!(mlr_plots, p)
end

p_mlr = plot(mlr_plots..., layout = (2, 2), size = (1200, 800))
plot!(plot_title = "MLR Comparison: Exponential vs Constant HoS")
savefig(p_mlr, joinpath(RESULTS_DIR, "hos_comparison_mlr.png"))
println("  Saved: hos_comparison_mlr.png")

# ============================================================================
# Plot 3: Temperature Comparison (2×2 grid)
# ============================================================================

temp_plots = []
for (i, mc_pct) in enumerate(MC_LEVELS)
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]

    p = plot(
        re.times,
        re.T_surface,
        xlabel = "Time [s]",
        ylabel = "Temperature [K]",
        title = "MC = $(Int(mc_pct))%",
        label = "Surface (exp)",
        linewidth = 2,
        color = :red,
        legend = :bottomright,
    )
    plot!(
        rc.times,
        rc.T_surface,
        label = "Surface (const)",
        linewidth = 2,
        color = :red,
        linestyle = :dash,
    )
    plot!(re.times, re.T_back, label = "Back (exp)", linewidth = 2, color = :blue)
    plot!(
        rc.times,
        rc.T_back,
        label = "Back (const)",
        linewidth = 2,
        color = :blue,
        linestyle = :dash,
    )

    push!(temp_plots, p)
end

p_temp = plot(temp_plots..., layout = (2, 2), size = (1200, 800))
plot!(plot_title = "Temperature Comparison: Exponential vs Constant HoS")
savefig(p_temp, joinpath(RESULTS_DIR, "hos_comparison_temperatures.png"))
println("  Saved: hos_comparison_temperatures.png")

# ============================================================================
# Plot 4: Mass Evolution (2×2 grid)
# ============================================================================

mass_plots = []
for (i, mc_pct) in enumerate(MC_LEVELS)
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]

    p = plot(
        re.times,
        re.total_mass .* 1000,
        xlabel = "Time [s]",
        ylabel = "Mass [g/m²]",
        title = "MC = $(Int(mc_pct))%",
        label = "Exponential",
        linewidth = 2,
        color = :black,
        legend = :topright,
        ylims = (0, :auto),
    )
    plot!(
        rc.times,
        rc.total_mass .* 1000,
        label = "Constant avg",
        linewidth = 2,
        color = colors_mc[i],
        linestyle = :dash,
    )

    push!(mass_plots, p)
end

p_mass = plot(mass_plots..., layout = (2, 2), size = (1200, 800))
plot!(plot_title = "Total Mass: Exponential vs Constant HoS")
savefig(p_mass, joinpath(RESULTS_DIR, "hos_comparison_mass.png"))
println("  Saved: hos_comparison_mass.png")

# ============================================================================
# Plot 5: Difference Plots (exponential − constant)
# ============================================================================

diff_plots = []

# MLR difference
p_dmlr = plot(
    xlabel = "Time [s]",
    ylabel = "ΔMLR [g/(m²·s)]",
    title = "MLR Difference (exp − const)",
    legend = :topright,
)
for (i, mc_pct) in enumerate(MC_LEVELS)
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]
    n = min(length(re.mlr), length(rc.mlr))
    plot!(
        re.times[1:n],
        (re.mlr[1:n] .- rc.mlr[1:n]) .* 1000,
        label = "MC=$(Int(mc_pct))%",
        linewidth = 1.5,
        color = colors_mc[i],
    )
end
hline!([0.0], linestyle = :dot, color = :gray, label = false)
push!(diff_plots, p_dmlr)

# Surface temperature difference
p_dts = plot(
    xlabel = "Time [s]",
    ylabel = "ΔT_surface [K]",
    title = "Surface Temp Difference (exp − const)",
    legend = :topright,
)
for (i, mc_pct) in enumerate(MC_LEVELS)
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]
    n = min(length(re.T_surface), length(rc.T_surface))
    plot!(
        re.times[1:n],
        re.T_surface[1:n] .- rc.T_surface[1:n],
        label = "MC=$(Int(mc_pct))%",
        linewidth = 1.5,
        color = colors_mc[i],
    )
end
hline!([0.0], linestyle = :dot, color = :gray, label = false)
push!(diff_plots, p_dts)

# Back temperature difference
p_dtb = plot(
    xlabel = "Time [s]",
    ylabel = "ΔT_back [K]",
    title = "Back Temp Difference (exp − const)",
    legend = :topright,
)
for (i, mc_pct) in enumerate(MC_LEVELS)
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]
    n = min(length(re.T_back), length(rc.T_back))
    plot!(
        re.times[1:n],
        re.T_back[1:n] .- rc.T_back[1:n],
        label = "MC=$(Int(mc_pct))%",
        linewidth = 1.5,
        color = colors_mc[i],
    )
end
hline!([0.0], linestyle = :dot, color = :gray, label = false)
push!(diff_plots, p_dtb)

# Average moisture concentration difference
p_dmoist = plot(
    xlabel = "Time [s]",
    ylabel = "Δξ_moisture [kg/m³]",
    title = "Avg Moisture Difference (exp − const)",
    legend = :topright,
)
for (i, mc_pct) in enumerate(MC_LEVELS)
    re = results["exponential_MC$(Int(mc_pct))"]
    rc = results["constant_MC$(Int(mc_pct))"]
    n = min(length(re.avg_moisture), length(rc.avg_moisture))
    plot!(
        re.times[1:n],
        re.avg_moisture[1:n] .- rc.avg_moisture[1:n],
        label = "MC=$(Int(mc_pct))%",
        linewidth = 1.5,
        color = colors_mc[i],
    )
end
hline!([0.0], linestyle = :dot, color = :gray, label = false)
push!(diff_plots, p_dmoist)

p_diff = plot(diff_plots..., layout = (2, 2), size = (1200, 800))
plot!(plot_title = "Model Differences: Exponential − Constant HoS")
savefig(p_diff, joinpath(RESULTS_DIR, "hos_comparison_delta.png"))
println("  Saved: hos_comparison_delta.png")

# ============================================================================
# Done
# ============================================================================

println("\n" * "="^60)
println("All results saved to: $RESULTS_DIR")
println("  - 8 time-series CSVs (timeseries_*.csv)")
println("  - 1 summary CSV (hos_comparison_summary.csv)")
println("  - 5 plots (hos_model_curves, hos_comparison_mlr/temperatures/mass/delta)")
println("="^60)
