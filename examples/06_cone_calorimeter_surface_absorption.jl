"""
Cone Calorimeter Simulation of PMMA with Surface-Only Radiation Absorption

This example demonstrates the use of the SURFACE_ABSORPTION radiation model,
which deposits all incident radiation in the surface cell for computationally
efficient simulation of opaque materials.

RADIATION MODEL:
Uses `radiation_model = SURFACE_ABSORPTION` which provides:
- O(1) computational cost (vs O(n) for Beer-Lambert)
- Tridiagonal Jacobian sparsity (vs upper triangular for Beer-Lambert)
- 5-15× speedup compared to high-α workaround
- Physically appropriate for opaque materials (optical thickness τ >> 1)

Physical Setup:
- Sample thickness: 5.8 mm (typical for cone calorimeter samples)
- Exposed surface (z=L): radiant heating + convective/radiative losses + mass release
- Back surface (z=0): re-radiation + convection, impermeable
- Surface regresses as material decomposes and gas escapes

Coordinate Convention:
- z = 0: Bottom (substrate/sample holder), boundary_id = 2, fixed
- z = L(t): Top (exposed surface), boundary_id = 1, moves with material
- Heat enters from top, material shrinks toward bottom
- Mesh follows the shrinking material (ALE formulation)
"""

using Pyrolysis
using Plots
using DelimitedFiles

# ============================================================================
# Material Definition: PMMA with Surface-Only Absorption
# ============================================================================

"""
Define PMMA material for surface absorption radiation model.

NOTE: With radiation_model = SURFACE_ABSORPTION, the absorption coefficient α
is NOT used - all radiation is deposited directly in the surface cell regardless
of material optical properties. We set α = 0 to make this explicit.

For comparison with in-depth absorption, see example 06_cone_calorimeter.jl
which uses α = 1.47 m²/kg with radiation_model = BEER_LAMBERT.
"""

# Temperature-dependent properties for PMMA (piecewise at T = 378 K); the fits
# are held constant outside their ~300–700 K measured range (see example 06).
pmma_cp = CallableProperty(
    T -> T < 378.0 ? (-2290.0 + 11.2 * max(T, 300.0)) : (1040.0 + 3.08 * T),
)
pmma_k = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * min(T, 700.0)))

pmma = Material(
    name = :PMMA_surface_abs,
    components = (
        # Virgin PMMA - solid phase
        # NOTE: α = 0 because SURFACE_ABSORPTION model ignores material absorption coefficients
        SolidComponent(
            :virgin,
            ρ = 1210.0,                    # kg/m³ - density
            c = pmma_cp,                   # J/(kg·K) - piecewise temperature-dependent
            k = pmma_k,                    # W/(m·K) - piecewise temperature-dependent
            ε = 0.96,                      # surface emissivity
            α = 0.0,                        # m²/kg - not used with SURFACE_ABSORPTION model
        ),
        # Char - solid residue
        SolidComponent(
            :char,
            ρ = 1210.0,                    # kg/m³ - same as virgin PMMA
            c = pmma_cp,                   # J/(kg·K) - same as virgin PMMA
            k = pmma_k,                    # W/(m·K) - same as virgin PMMA
            ε = 0.96,                      # same emissivity as virgin PMMA
            α = 0.0,                        # m²/kg - not used with SURFACE_ABSORPTION model
        ),
        # MMA gas - gaseous product
        GaseousComponent(
            :MMA,
            M = 0.10012,                   # kg/mol - molar mass of MMA
            c = 1500.0,                    # J/(kg·K) - gas heat capacity
            k = 0.02,                      # W/(m·K) - gas conductivity
            λ = 1e-5,                       # m²/s - gas transfer coefficient
        ),
    ),
    reactions = (
        # Single-step decomposition reaction with small char residue
        Reaction(
            :decomposition,
            1 => (2, 3),
            A = 1.5e14,                    # s⁻¹ - pre-exponential factor
            E = 2.03e5,                    # J/mol - activation energy
            h = 820e3,                     # J/kg - heat of reaction (positive = endothermic)
            n = 1.0,                       # reaction order
            yields = (0.002, 0.998),        # 0.2% char, 99.8% gas
        ),
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
)

# ============================================================================
# Simulation Parameters
# ============================================================================

# Sample geometry
sample_thickness = 5.8e-3    # 5.8 mm typical cone calorimeter sample
n_cells = 72                 # 72 cells → 0.08 mm resolution

# Time parameters  
t_final = 1600.0             # simulation duration
dt_initial = 1e-3            # initial timestep

# Environmental conditions
T_ambient = 300.0            # K - ambient temperature (27°C)
q_external = 0.9813 * 25e3   # W/m² - 25 kW/m² external heat flux

# Convective heat transfer coefficient
h_conv = 7.2                 # W/(m²·K)

# ============================================================================
# Mesh Creation
# ============================================================================

mesh = create_uniform_mesh(
    sample_thickness,
    n_cells,
    Val(3);
    T_initial = T_ambient,
    ξ_initial = (1210.0, 0.0, 0.0),
)

println("=" ^ 70)
println("SURFACE-ONLY RADIATION ABSORPTION TEST")
println("=" ^ 70)
println()
println("PURPOSE: Test if radiation absorption mechanism causes ThermaKin vs")
println("         Pyrolysis.jl differences by eliminating in-depth absorption")
println()
println("ABSORPTION COEFFICIENT:")
println("  Original: α = 1.47 m²/kg (semi-transparent)")
println("  Modified: α = $(α_surface_only) m²/kg (surface-only)")
println()
println("EXPECTED RESULT:")
println("  With surface-only absorption, Monte Carlo and continuous Beer-Lambert")
println("  should give IDENTICAL results since all energy goes to surface cell.")
println()
println("Mesh created:")
println("  - Thickness: $(sample_thickness*1000) mm")
println("  - Number of cells: $n_cells")
println("  - Cell size: $(sample_thickness/n_cells*1000) mm")

# ============================================================================
# Boundary Conditions
# ============================================================================

top_thermal_bc =
    RadiativeFluxBC(flux = q_external, absorptivity = 0.96) +
    RadiativeBC(ε = 0.96, T_inf = 330) +
    ConvectiveBC(h = h_conv, T_inf = 330)

top_mass_bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)

bottom_thermal_bc = RadiativeBC(ε = 0.95, T_inf = 307) + ConvectiveBC(h = 4, T_inf = 307)
bottom_mass_bc = ImpermeableBC()

bc_set = BoundaryConditionSet(
    thermal = (top = top_thermal_bc, bottom = bottom_thermal_bc),
    mass = (top = top_mass_bc, bottom = bottom_mass_bc),
)

println("\nBoundary conditions:")
println(
    "  - Top surface: RadiativeFluxBC ($(q_external/1000) kW/m², α=0.96) + re-radiation + convection",
)
println("  - Bottom surface: Re-radiation + convection, impermeable")

# ============================================================================
# Problem Setup
# ============================================================================

T_initial_func = z -> T_ambient
ξ_initial_funcs = [z -> 1210.0, z -> 0.0, z -> 0.0]

problem = PyrolysisProblem(
    mesh = mesh,
    material = pmma,
    bc_set = bc_set,
    T_initial = T_initial_func,
    ξ_initial = ξ_initial_funcs,
    tspan = (0.0, t_final),
    dt_initial = dt_initial,
    output_times = collect(range(0.0, t_final, length = 601)),
)

println("\nProblem created:")
println("  - Material: $(pmma.name)")
println("  - Duration: $(t_final) s")
println("  - External flux: $(q_external/1000) kW/m²")
println("  - Initial temperature: $(T_ambient) K")

# ============================================================================
# Solve with ALE Mesh Motion
# ============================================================================

println("\nSolving with ALE mesh motion...")
println("This may take a few minutes...\n")

# Solve using implicit solver with ALE enabled and SURFACE_ABSORPTION radiation model
solution = solve(
    problem;
    abstol = 1e-8,
    reltol = 1e-6,
    dtmin = 1e-14,
    progress = true,
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,  # All radiation deposited in surface cell
    min_thickness = 0.01e-3,
    handle_depletion = true,
    depletion_threshold = 0.05,
    min_cells = 1,
)

println("\nSimulation completed!")
println("  - Return code: $(solution.retcode)")
println("  - Solve time: $(round(solution.solve_time, digits=2)) s")
println(
    "  - Final thickness: $(round(solution.extras.thickness_history[end]*1000, digits=2)) mm",
)
println("  - Number of output times: $(length(solution.t))")

# ============================================================================
# Results Analysis
# ============================================================================

println("\n" * "=" ^ 70)
println("RESULTS SUMMARY (SURFACE_ABSORPTION MODEL)")
println("=" ^ 70)
println("Radiation model: SURFACE_ABSORPTION")
println("  - All radiation deposited in surface cell")
println("  - O(1) computational cost")
println("  - Tridiagonal Jacobian sparsity")
println("  - Ideal for opaque materials")
println("=" ^ 70)

# Extract key results at t = 600s (quasi-steady state)
idx_600 = findfirst(t -> t >= 600.0, solution.t)
if idx_600 !== nothing
    T_surface_600 = solution.surface_T[idx_600]
    mlr_600 = solution.extras.gas_generation_rate[idx_600] * 1000  # kg/(m²·s) → g/(m²·s)

    println("\nAt t = 600s (quasi-steady state):")
    println("  Surface temperature: $(round(T_surface_600, digits=1)) K")
    println("  Mass loss rate: $(round(mlr_600, digits=2)) g/(m²·s)")
    println()
    println("Compare with ThermaKin results using PMMA_cone_surface_absorption.cmp")
    println("If results match, radiation absorption algorithm is confirmed as root cause.")
end

# ============================================================================
# Plotting
# ============================================================================

println("\nGenerating plots...")

# Create output directory
output_dir = joinpath(@__DIR__, "..", "output")
mkpath(output_dir)

# Plot 1: Mass Loss Rate vs Time
p1 = plot(
    solution.t,
    solution.extras.gas_generation_rate .* 1000,
    xlabel = "Time [s]",
    ylabel = "MLR [g/(m²·s)]",
    title = "Mass Loss Rate (Surface-Only Absorption)",
    label = "Pyrolysis.jl (α=1e6)",
    linewidth = 2,
    legend = :topright,
)

# Load ThermaKin surface absorption reference data if available
thermakin_output_dir =
    joinpath(@__DIR__, "..", "docs", "thermakin_surface_absorption_output")
if isdir(thermakin_output_dir)
    println("Loading ThermaKin surface absorption reference data...")

    # MLR data: time_s, total_mass_loss_rate_kg_s
    thermakin_mlr_data =
        readdlm(joinpath(thermakin_output_dir, "pmma_cone_mlr.csv"), ',', skipstart = 1)
    thermakin_time = thermakin_mlr_data[:, 1]
    thermakin_mlr = thermakin_mlr_data[:, 2] .* 1000  # Convert kg/(m²·s) to g/(m²·s)

    # Front temperature data: time_s, from_back_m, front_temperature_K
    thermakin_temp_data = readdlm(
        joinpath(thermakin_output_dir, "pmma_cone_front_temp.csv"),
        ',',
        skipstart = 1,
    )
    thermakin_tc2 = thermakin_temp_data[:, 3]  # Front (surface) temperature

    # Summary data: time_s, mass_kg_m2, thickness_m
    thermakin_summary_data =
        readdlm(joinpath(thermakin_output_dir, "pmma_cone_summary.csv"), ',', skipstart = 1)
    thermakin_thickness = thermakin_summary_data[:, 3] .* 1000  # Convert m to mm

    # Back temperature - extract from spatial data
    # The parser converts all coordinates to from_back convention:
    #   - from_back = 0 is the back (substrate)
    #   - from_back = thickness is the front (heated surface)
    # Use efficient single-pass algorithm to find minimum position (back) for each time
    thermakin_spatial_path = joinpath(thermakin_output_dir, "pmma_cone_spatial.csv")
    println("  Loading spatial data (this may take a moment for large files)...")
    thermakin_spatial_data = readdlm(thermakin_spatial_path, ',', skipstart = 1)
    spatial_times = thermakin_spatial_data[:, 2]
    spatial_positions = thermakin_spatial_data[:, 3]  # from_back_m (always from_back convention)
    spatial_temps = thermakin_spatial_data[:, 4]
    println("  Processing $(length(spatial_times)) spatial data points...")

    # Build dictionary: time -> (min_from_back, temp_at_min) for back temperature
    back_temp_by_time = Dict{Float64,Tuple{Float64,Float64}}()
    for i in eachindex(spatial_times)
        t = spatial_times[i]
        pos = spatial_positions[i]
        temp = spatial_temps[i]
        if !haskey(back_temp_by_time, t) || pos < back_temp_by_time[t][1]
            back_temp_by_time[t] = (pos, temp)
        end
    end

    # Extract back temperatures for each MLR timestep - O(n)
    thermakin_tc1 = zeros(length(thermakin_time))
    for (i, t) in enumerate(thermakin_time)
        if haskey(back_temp_by_time, t)
            thermakin_tc1[i] = back_temp_by_time[t][2]
        else
            thermakin_tc1[i] = 300.0
        end
    end
    println("  Done loading ThermaKin data.")

    # Add ThermaKin to MLR plot
    plot!(
        p1,
        thermakin_time,
        thermakin_mlr,
        linewidth = 2,
        linestyle = :dash,
        label = "ThermaKin (α=1e6)",
    )

    has_thermakin_data = true
else
    println("ThermaKin surface absorption data not found at: $thermakin_output_dir")
    println("Run ThermaKin with PMMA_cone_surface_absorption.cmp first.")
    has_thermakin_data = false
end

savefig(p1, joinpath(output_dir, "mlr_surface_absorption.png"))

# Plot 2: Surface and Back Temperature vs Time (4 curves)
p2 = plot(
    solution.t,
    solution.surface_T,
    xlabel = "Time [s]",
    ylabel = "Temperature [K]",
    title = "Temperature Comparison (Surface-Only Absorption)",
    label = "Pyrolysis.jl Surface (z=L)",
    linewidth = 2,
    color = :red,
    legend = :bottomright,
)
plot!(
    p2,
    solution.t,
    solution.extras.back_temperature,
    linewidth = 2,
    color = :blue,
    label = "Pyrolysis.jl Back (z=0)",
)
if has_thermakin_data
    plot!(
        p2,
        thermakin_time,
        thermakin_tc2,
        linewidth = 2,
        linestyle = :dash,
        color = :red,
        label = "ThermaKin Surface (z=L)",
    )
    plot!(
        p2,
        thermakin_time,
        thermakin_tc1,
        linewidth = 2,
        linestyle = :dash,
        color = :blue,
        label = "ThermaKin Back (z=0)",
    )
end
savefig(p2, joinpath(output_dir, "temperature_comparison_surface_absorption.png"))

# Plot 3: Sample Thickness vs Time
p3 = plot(
    solution.t,
    solution.extras.thickness_history .* 1000,
    xlabel = "Time [s]",
    ylabel = "Thickness [mm]",
    title = "Sample Thickness (Surface-Only Absorption)",
    label = "Pyrolysis.jl (α=1e6)",
    linewidth = 2,
    legend = :topright,
)
if has_thermakin_data
    plot!(
        p3,
        thermakin_time,
        thermakin_thickness,
        linewidth = 2,
        linestyle = :dash,
        label = "ThermaKin (α=1e6)",
    )
end
savefig(p3, joinpath(output_dir, "thickness_surface_absorption.png"))

println("Plots saved to: $output_dir")

# ============================================================================
# Export Results for Comparison
# ============================================================================

# Export time series data
output_file = joinpath(output_dir, "pyrolysis_surface_absorption_results.csv")
open(output_file, "w") do io
    println(io, "# Pyrolysis.jl Surface-Only Absorption Results")
    println(io, "# α = $(α_surface_only) m²/kg (surface-only, was 1.47)")
    println(io, "time_s,mlr_g_m2_s,surface_temp_K,thickness_mm")
    for i in eachindex(solution.t)
        println(
            io,
            "$(solution.t[i]),$(solution.extras.gas_generation_rate[i]*1000),$(solution.surface_T[i]),$(solution.extras.thickness_history[i]*1000)",
        )
    end
end
println("Results exported to: $output_file")

# ============================================================================
# Quantitative Comparison with ThermaKin (if data available)
# ============================================================================

println("\n" * "=" ^ 70)
if has_thermakin_data
    println("COMPARISON: Pyrolysis.jl vs ThermaKin (Surface-Only Absorption)")
    println("=" ^ 70)

    # Find values at t = 600s
    idx_pj = findfirst(t -> t >= 600.0, solution.t)
    idx_tk = findfirst(t -> t >= 600.0, thermakin_time)

    if idx_pj !== nothing && idx_tk !== nothing
        T_surf_pj = solution.surface_T[idx_pj]
        T_surf_tk = thermakin_tc2[idx_tk]
        mlr_pj = solution.extras.gas_generation_rate[idx_pj] * 1000
        mlr_tk = thermakin_mlr[idx_tk]

        println("\nAt t = 600s (quasi-steady state):")
        println("─" ^ 70)
        println("                           Pyrolysis.jl    ThermaKin    Difference")
        println("─" ^ 70)
        println(
            "Surface temperature [K]    $(round(T_surf_pj, digits=1))           $(round(T_surf_tk, digits=1))         $(round(T_surf_pj - T_surf_tk, digits=1))",
        )
        println(
            "Mass loss rate [g/(m²·s)]  $(round(mlr_pj, digits=2))            $(round(mlr_tk, digits=2))          $(round(mlr_pj - mlr_tk, digits=2))",
        )
        println("─" ^ 70)

        # Interpretation
        T_diff = abs(T_surf_pj - T_surf_tk)
        mlr_diff = abs(mlr_pj - mlr_tk)

        println("\nINTERPRETATION:")
        if T_diff < 5.0 && mlr_diff < 0.5
            println("  ✓ Results MATCH within tolerance (ΔT < 5K, ΔMLR < 0.5 g/m²/s)")
            println("  ✓ This CONFIRMS radiation absorption algorithm is the root cause")
            println("    of differences seen with semi-transparent (α=1.47) case.")
        else
            println(
                "  ✗ Results DIFFER significantly (ΔT = $(round(T_diff, digits=1))K, ΔMLR = $(round(mlr_diff, digits=2)) g/m²/s)",
            )
            println("  → Other factors besides radiation absorption may be involved.")
            println("  → Check boundary conditions, material properties, mesh resolution.")
        end
    end
else
    println("NEXT STEPS:")
    println("=" ^ 70)
    println("1. Run ThermaKin with PMMA_cone_surface_absorption.cmp")
    println("2. Place output files in docs/thermakin_surface_absorption_output/")
    println("3. Re-run this script to compare results")
    println("4. If results match: Confirms radiation absorption algorithm is root cause")
    println("5. If results differ: Other factors may be involved")
end
println("=" ^ 70)
