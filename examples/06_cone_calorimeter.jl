"""
Cone Calorimeter Simulation of PMMA with ALE Mesh Motion

This example simulates the pyrolysis of PMMA (Poly(methyl methacrylate))
under standard cone calorimeter conditions with ALE (Arbitrary Lagrangian-Eulerian)
mesh motion to track the shrinking material as it decomposes:
- 1D simulation through sample thickness
- ALE mesh motion for accurate shrinkage modeling

Reference: Standard cone calorimeter test per ISO 5660-1

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

ALE Mesh Motion:
- Node velocities are the cumulative volumetric strain rate from the pinned bottom
  node upward; the mesh moves with the solid
- Cell volumes are recomputed from the node positions at every evaluation
- Simulation terminates when thickness falls below minimum threshold

PMMA Properties (from literature):
- Virgin PMMA density: 1190 kg/m³
- Specific heat capacity: ~1420 J/(kg·K) (with slight temperature dependence)
- Thermal conductivity: ~0.21 W/(m·K)
- Emissivity: 0.96
- Decomposition: single-step to MMA monomer + small char residue
- Products: 0.2% char, 99.8% MMA gas
- Activation energy: ~203 kJ/mol
- Pre-exponential factor: ~1.5×10¹⁴ s⁻¹
- Heat of decomposition: ~820 kJ/kg (endothermic)
"""

using Pyrolysis
using Plots
using DelimitedFiles

# ============================================================================
# Material Definition: PMMA
# ============================================================================

"""
Define PMMA material with temperature-dependent properties.

PMMA undergoes thermal depolymerization, releasing MMA gas with a small char residue:
	PMMA (solid) → 0.2% char (solid) + 99.8% MMA (gas)

The reaction is endothermic (absorbs heat).
"""

# Create PMMA material with three components:
# 1. Virgin PMMA (solid)
# 2. Char residue (solid)
# 3. MMA gas (gaseous product)

# Temperature-dependent properties for PMMA (piecewise at T = 378 K). The
# fits are measured on ~300–700 K and are held constant outside that range:
# the low-T cp line is negative below 204 K and the high-T k line crosses
# zero at 872 K, so extrapolating either would give an unphysical law.
# Heat capacity: cp = -2290 + 11.2*T for T < 378 K
#                cp = 1040 + 3.08*T  for T > 378 K
pmma_cp = CallableProperty(
    T -> T < 378.0 ? (-2290.0 + 11.2 * max(T, 300.0)) : (1040.0 + 3.08 * T),
)

# Thermal conductivity: k = 0.15 for T < 378 K
#                       k = 0.34 - 3.9e-4*T for T > 378 K
pmma_k = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * min(T, 700.0)))

pmma = Material(
    name = :PMMA,
    components = (
        # Virgin PMMA - solid phase
        SolidComponent(
            :virgin,
            ρ = 1210.0,                    # kg/m³ - density
            c = pmma_cp,                   # J/(kg·K) - piecewise temperature-dependent
            k = pmma_k,                    # W/(m·K) - piecewise temperature-dependent
            ε = 0.96,                      # surface emissivity
            α = 1.47,                        # m²/kg - absorption coefficient for in-depth radiation
        ),
        # Char - solid residue (same properties as virgin PMMA)
        SolidComponent(
            :char,
            ρ = 1210.0,                    # kg/m³ - same as virgin PMMA
            c = pmma_cp,                   # J/(kg·K) - same as virgin PMMA
            k = pmma_k,                    # W/(m·K) - same as virgin PMMA
            ε = 0.96,                      # same emissivity as virgin PMMA
            α = 1.47,                       # m²/kg - same absorption coefficient as virgin PMMA
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
        # virgin PMMA (index 1) → char (index 2) + MMA gas (index 3)
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
n_cells = 72               # 72 cells → 0.08 mm resolution

# Time parameters  
t_final = 1600.0            # 10 minutes - typical cone test duration
dt_initial = 1e-3          # initial timestep

# Environmental conditions
T_ambient = 300.0          # K - ambient temperature (27°C)
q_external = 0.9813 * 25e3          # W/m² - 25 kW/m² external heat flux, CAPA adjusted for radial dependence

# Convective heat transfer coefficient for cone calorimeter
# Natural convection + some forced convection from exhaust
h_conv = 7.2              # W/(m²·K)

# ============================================================================
# Mesh Creation
# ============================================================================

# Create uniform mesh; the initial temperature and concentrations (virgin PMMA
# + char + MMA gas) come from the problem's T_initial and ξ_initial below.
# Convention: z=0 is the back face, z=L the exposed surface
mesh = create_uniform_mesh(
    sample_thickness,      # domain length [m]
    n_cells,               # number of cells
)

println("Mesh created:")
println("  - Thickness: $(sample_thickness*1000) mm")
println("  - Number of cells: $n_cells")
println("  - Cell size: $(sample_thickness/n_cells*1000) mm")

# ============================================================================
# Boundary Conditions
# ============================================================================

# Top surface (boundary_id = 1): exposed surface
# - External radiant heat flux (25 kW/m²) with in-depth absorption controlled by radiation_model
# - Re-radiation losses (surface radiates to surroundings)
# - Convective losses (natural + forced convection)
#
# IMPORTANT: Radiation model comparison
# - RadiativeFluxBC: External radiative source (cone heater)
#   → Subject to absorptivity (fraction absorbed = α × incident)
#   → radiation_model=BEER_LAMBERT: absorbed radiation distributed via Beer-Lambert
#   → radiation_model=SURFACE_ABSORPTION: all absorbed in surface cell
# - HeatFluxBC: Direct contact heating or fully opaque surface absorption
#   → NOT subject to absorptivity - full flux applied at surface
# - RadiativeBC: Re-radiation losses to surroundings (σεT⁴)
#   → Loss term, not a heat input

top_thermal_bc =
    RadiativeFluxBC(flux = q_external, absorptivity = 0.96) +  # Cone heater (model-dependent)
    # RadiativeBC(ε = 0.96, T_inf = 360) +                 # Re-radiation losses (original)
    RadiativeBC(ε = 0.96, T_inf = 330) +                 # Re-radiation losses (ThermaKin match)
    ConvectiveBC(h = h_conv, T_inf = 330)                # Convective losses

# Top surface mass BC: diffusive gas escape
# DiffusiveMassBC uses J = λ × ξ / d (internal diffusion model) by default
# This is consistent with the surface regression calculation
# top_mass_bc = DiffusiveMassBC()  # Pure internal diffusion (h_m = Inf) (original)
top_mass_bc = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0)  # ThermaKin match: h_m = 0.05 m/s

# Bottom surface (boundary_id = 2): back face with heat losses
# - Re-radiation losses to ambient (back face radiates to surroundings)
# - Convective losses (natural convection from back face)
bottom_thermal_bc =
    RadiativeBC(ε = 0.95, T_inf = 307) +   # Re-radiation losses
    ConvectiveBC(h = 4, T_inf = 307)  # Convective losses
bottom_mass_bc = ImpermeableBC()   # Impermeable (no gas escape)

# Collect boundary conditions into BoundaryConditionSet
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

# Initial conditions
T_initial_func = z -> T_ambient           # Uniform initial temperature
ξ_initial_funcs = [
    z -> 1210.0,    # Virgin PMMA at full density
    z -> 0.0,       # No char initially
    z -> 0.0,       # No gas initially
]

# Create the pyrolysis problem
# Radiative heating is supplied via RadiativeFluxBC in `bc_set`.
problem = PyrolysisProblem(
    mesh = mesh,
    material = pmma,
    bc_set = bc_set,
    T_initial = T_initial_func,
    ξ_initial = ξ_initial_funcs,
    tspan = (0.0, t_final),
    dt_initial = dt_initial,
    # PERFORMANCE: Reduce output times from 6030 to 601 points (every 1s instead of 0.1s)
    # This significantly reduces interpolation overhead and memory usage
    output_times = collect(range(0.0, t_final, length = 601)),  # Output every 1s
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

# Solve using implicit solver with ALE enabled
# - use_ale=true: enables integrated ALE mesh motion (node positions in state vector)
# - min_thickness: stops simulation when material becomes too thin
# - deactivation: removes the surface cell once its condensed mass is below 1% of its
#   initial mass (the default `threshold`); the mesh correction thins the layer so that
#   it goes in one event
# - dtmin: allows very small timesteps during rapid reaction fronts
solution = solve(
    problem;
    abstol = 1e-8,
    reltol = 1e-6,
    dtmin = 1e-14,    # Allow very small timesteps during rapid reaction
    progress = true,
    use_ale = true,                        # Enable ALE mesh motion
    min_thickness = 0.01e-3,               # Stop when thickness < 0.01 mm
    deactivation = true,                   # Deactivate burned-out surface cells
)

println("\nSolution completed!")
println("  - Return code: $(solution.retcode)")
println("  - Solve time: $(round(solution.solve_time, digits=2)) s")

# Thickness history is now part of the solution (from state vector)
thickness_times = solution.t
thickness_values = solution.extras.thickness_history
final_thickness = isempty(thickness_values) ? sample_thickness : thickness_values[end]
println("  - Final thickness: $(round(final_thickness*1000, digits=2)) mm")

# ============================================================================
# Results Analysis
# ============================================================================

# Extract key results
times = solution.t
gas_generation_rate = solution.extras.gas_generation_rate  # kg/(m²·s) - volumetric gas production
T_surface = solution.surface_T  # K
T_back = solution.extras.back_temperature     # K
total_mass = solution.extras.total_mass       # kg/m²

println("\nKey Results:")
println(
    "  - Initial mass: $(round(total_mass[1]*1000, digits=2)) g/m² = $(round(total_mass[1], digits=4)) kg/m²",
)
println(
    "  - Final mass: $(round(total_mass[end]*1000, digits=2)) g/m² = $(round(total_mass[end], digits=4)) kg/m²",
)
println("  - Mass loss: $(round((1 - total_mass[end]/total_mass[1])*100, digits=1))%")
println(
    "  - Peak gas generation rate: $(round(maximum(gas_generation_rate)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak surface temperature: $(round(maximum(T_surface), digits=1)) K = $(round(maximum(T_surface)-273.15, digits=1)) °C",
)
println("  - Initial thickness: $(sample_thickness*1000) mm")
println("  - Final thickness: $(round(final_thickness*1000, digits=2)) mm")
println(
    "  - Thickness change: $(round((1 - final_thickness/sample_thickness)*100, digits=1))%",
)

# Diagnostic: compare mass and thickness evolution
println("\n" * "-" ^ 60)
println("Mass vs Thickness Diagnostic (separating solid and gas):")
println("-" ^ 60)
ρ_solid = 1210.0  # PMMA solid intrinsic density
# Focus on the active pyrolysis phase
for i = 1:length(times)
    t_i = times[i]
    L_i = thickness_values[i]
    m_total = total_mass[i]

    # Compute SOLID mass only (component 1 = virgin PMMA)
    ξ_profile = solution.ξ[:, :, i]
    n_cells = size(ξ_profile, 2)

    # Get cell volumes from node positions at this time
    # Note: we don't have direct access to volumes, so use total_mass calculation
    # For now, compute solid fraction as ξ_solid / (ξ_solid + ξ_gas)
    total_solid = 0.0
    total_gas = 0.0
    for j = 1:n_cells
        total_solid += ξ_profile[1, j]  # solid component (index 1)
        total_gas += ξ_profile[2, j]    # gas component (index 2)
    end

    # Average concentrations (per cell, which ~ per unit volume since cells are similar size)
    avg_solid_conc = total_solid / n_cells
    avg_gas_conc = total_gas / n_cells

    # Solid mass ≈ avg_solid_conc × L × A
    # The total_mass already integrates correctly, so solid_mass = total_mass × (solid_fraction)
    # Actually, let's be more careful: solid_mass = Σ ξ_solid[j] × V[j]
    # For now, assume uniform cell sizes: solid_mass ≈ (Σ ξ_solid / Σ(ξ_solid + ξ_gas)) × total_mass
    if (total_solid + total_gas) > 1e-10
        solid_fraction = total_solid / (total_solid + total_gas)
    else
        solid_fraction = 0.0
    end
    m_solid = m_total * solid_fraction

    # Expected thickness from solid mass
    L_expected = m_solid / ρ_solid
    ratio = L_expected > 1e-10 ? L_i / L_expected : Inf

    # Only print at key points
    should_print = (i <= 10) || (i >= length(times)-2) || (mod(i, 10) == 0)  # Print every 10th point
    if should_print
        println(
            "  t=$(round(t_i, digits=1))s: L=$(round(L_i*1000, digits=3))mm, m_solid=$(round(m_solid, digits=4))kg/m², m_gas=$(round(m_total-m_solid, digits=4))kg/m², L_exp=$(round(L_expected*1000, digits=3))mm, ratio=$(round(ratio, digits=3)), n=$n_cells",
        )
    end
end
println("-" ^ 60)

# Find time to ignition (when gas generation exceeds critical value ~3-6 g/(m²·s))
# Note: This uses gas_generation_rate which is a proxy for actual MLR
gen_rate_critical = 3e-3  # kg/(m²·s) ≈ 3 g/(m²·s) - approximate ignition threshold
tig_idx = findfirst(gas_generation_rate .> gen_rate_critical)
if tig_idx !== nothing
    t_ignition = times[tig_idx]
    println("  - Time to ignition (generation > 3 g/m²s): $(round(t_ignition, digits=1)) s")
end

# ============================================================================
# Plotting Results
# ============================================================================

println("\nGenerating plots...")

# ThermaKin reference outputs for this case (`docs/thermakin_output/`, from a
# ThermaKin run of `docs/PMMA_cone.cmp` + `docs/PMMA_cone.cnd`). The overlay is
# skipped when they are absent.
thermakin_output_dir = joinpath(@__DIR__, "..", "docs", "thermakin_output")
has_thermakin = isdir(thermakin_output_dir)
if has_thermakin

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
    # Spatial data columns: conc_PMMA_glass_kg_m3, time_s, from_back_m, temperature_K, ...
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
else
    println(
        "  No ThermaKin reference outputs in $(thermakin_output_dir); plotting Pyrolysis.jl only.",
    )
end

# Create figure with 4 subplots
p1 = plot(
    times,
    gas_generation_rate .* 1000,  # Convert to g/(m²·s)
    xlabel = "Time [s]",
    ylabel = "Gas Generation [g/(m²·s)]",
    title = "Gas Generation Rate (proxy for MLR)",
    legend = :topright,
    linewidth = 2,
    color = :red,
    label = "Pyrolysis.jl",
)
has_thermakin && plot!(
    thermakin_time,
    thermakin_mlr,
    linewidth = 2,
    color = :blue,
    linestyle = :dash,
    label = "ThermaKin MLR",
)
hline!(
    [gen_rate_critical*1000],
    linestyle = :dot,
    color = :gray,
    label = "Ignition threshold",
)

p2 = plot(
    times,
    T_surface,  # Already in Kelvin
    xlabel = "Time [s]",
    ylabel = "Temperature [K]",
    title = "Surface & Back Temperature",
    legend = :bottomright,
    linewidth = 2,
    color = :orange,
    label = "Surface (Pyrolysis.jl)",
)
has_thermakin && plot!(
    thermakin_time,
    thermakin_tc2,
    linewidth = 2,
    color = :red,
    linestyle = :dash,
    label = "Surface (ThermaKin)",
)
plot!(times, T_back, linewidth = 2, color = :blue, label = "Back (Pyrolysis.jl)")
has_thermakin && plot!(
    thermakin_time,
    thermakin_tc1,
    linewidth = 2,
    color = :green,
    linestyle = :dash,
    label = "Back (ThermaKin)",
)

p3 = plot(
    times,
    total_mass .* 1000,  # Convert to g/m²
    xlabel = "Time [s]",
    ylabel = "Mass [g/m²]",
    title = "Total Mass",
    legend = false,
    linewidth = 2,
    color = :green,
    ylims = (0, :auto),  # Start y-axis at zero for perspective on mass loss
)

# Sample thickness over time (from ALE data) - compare with ThermaKin
p4 = plot(
    thickness_times,
    thickness_values .* 1000,  # Convert to mm
    xlabel = "Time [s]",
    ylabel = "Thickness [mm]",
    title = "Sample Thickness",
    legend = :topright,
    linewidth = 2,
    color = :purple,
    label = "Pyrolysis.jl",
)
has_thermakin && plot!(
    thermakin_time,
    thermakin_thickness,
    linewidth = 2,
    color = :blue,
    linestyle = :dash,
    label = "ThermaKin",
)
hline!([0.1], linestyle = :dot, color = :gray, label = "Min thickness")

# Temperature profiles at selected times
# Note: For ALE, depth is relative to current (shrinking) domain
p5 = plot(xlabel = "Depth [mm]", ylabel = "Temperature [K]", title = "Temperature Profiles")
profile_times = [30, 60, 120, 300, 600, 1200, 1600]  # seconds
for t in profile_times
    # Find closest output time
    idx = argmin(abs.(times .- t))
    if idx <= size(solution.T, 2)
        T_profile = solution.T[:, idx]
        # Skip if only 1 cell remains (can't make a profile)
        if length(T_profile) <= 1
            continue
        end
        # Use current thickness at this time (approximate from nearest ALE data point)
        t_idx = isempty(thickness_times) ? 1 : argmin(abs.(thickness_times .- times[idx]))
        current_thickness =
            isempty(thickness_values) ? sample_thickness : thickness_values[t_idx]
        z_coords = range(0, current_thickness*1000, length = length(T_profile))
        plot!(
            z_coords,
            T_profile,
            label = "t=$(round(times[idx], digits=0))s",
            linewidth = 1.5,
        )
    end
end

# Combine all plots (now 5 plots in a 2x3 layout, with one empty)
final_plot = plot(p1, p2, p3, p4, p5, layout = (2, 3), size = (1200, 800))
plot!(plot_title = "Cone Calorimeter Simulation: PMMA at 25 kW/m² (ALE Mesh)")

# Save plot
savefig(final_plot, "cone_calorimeter_pmma_results.png")
println("Results saved to: cone_calorimeter_pmma_results.png")

# Display plot
display(final_plot)

# ============================================================================
# Summary Table
# ============================================================================

println("\n" * "="^60)
println("CONE CALORIMETER SIMULATION SUMMARY (ALE)")
println("="^60)
println("Material: PMMA (Poly(methyl methacrylate))")
println("Initial sample thickness: $(sample_thickness*1000) mm")
println("External heat flux: $(q_external/1000) kW/m²")
println("ALE mesh motion: enabled")
println("-"^60)
println("Simulation Results:")
println(
    "  Time to generation > 3 g/(m²·s): $(tig_idx !== nothing ? "$(round(times[tig_idx], digits=1)) s" : "N/A")",
)
println(
    "  Peak gas generation: $(round(maximum(gas_generation_rate)*1000, digits=2)) g/(m²·s)",
)
println("  Peak surface temperature: $(round(maximum(T_surface), digits=1)) K")
println("  Final mass loss: $(round((1 - total_mass[end]/total_mass[1])*100, digits=1))%")
println("-"^60)
println("ALE Mesh Results:")
println("  Final thickness: $(round(final_thickness*1000, digits=2)) mm")
println(
    "  Thickness reduction: $(round((1 - final_thickness/sample_thickness)*100, digits=1))%",
)
println(
    "  Simulation ended: $(solution.retcode == :Terminated ? "Material burned through" : string(solution.retcode))",
)
println("="^60)

# ============================================================================
# Detailed Energy Balance Analysis (for comparison with ThermaKin)
# ============================================================================

println("\n" * "="^70)
println("DETAILED ENERGY BALANCE ANALYSIS AT QUASI-STEADY STATE")
println("="^70)

# Find quasi-steady state period (around 400-600s)
for target_t in [200.0, 400.0, 600.0, 800.0]
    idx = argmin(abs.(times .- target_t))
    if idx > length(times)
        continue
    end
    t = times[idx]

    # Extract state at this time
    T_surf = T_surface[idx]
    T_back_val = T_back[idx]
    gen_rate_val = gas_generation_rate[idx]
    thick = thickness_values[idx]

    # Skip if not enough data
    if idx > size(solution.T, 2)
        continue
    end
    T_prof = solution.T[:, idx]
    n_cells = length(T_prof)

    # Energy balance components at top surface
    # External radiation
    q_ext = q_external

    # Re-radiation (using our BC parameters)
    T_inf_rad = 330.0  # From BC setup
    ε_surf = 0.96
    σ = 5.670374419e-8
    q_rerad = ε_surf * σ * (T_surf^4 - T_inf_rad^4)

    # Convection
    T_inf_conv = 330.0
    q_conv_loss = h_conv * (T_surf - T_inf_conv)

    # Net heat into material
    q_absorbed = q_ext * ε_surf  # Absorbed external radiation
    q_net = q_absorbed - q_rerad - q_conv_loss

    # Reaction heat absorption
    h_rxn = 820e3  # J/kg
    q_rxn = gen_rate_val * h_rxn

    # Conduction estimate (q_net - q_rxn goes into conduction + storage)
    q_remaining = q_net - q_rxn

    println("\n--- t = $(round(t, digits=1)) s ---")
    println("State:")
    println("  Surface T: $(round(T_surf, digits=1)) K")
    println("  Back T: $(round(T_back_val, digits=1)) K")
    println("  MLR: $(round(gen_rate_val*1000, digits=3)) g/(m²·s)")
    println("  Thickness: $(round(thick*1000, digits=3)) mm")
    println("  Cells: $n_cells")

    println("\nEnergy flux at top surface [kW/m²]:")
    println("  External radiation (absorbed): +$(round(q_absorbed/1000, digits=2))")
    println("  Re-radiation loss:             -$(round(q_rerad/1000, digits=2))")
    println("  Convection loss:               -$(round(q_conv_loss/1000, digits=2))")
    println("  ----------------------------------------")
    println("  NET into material:             +$(round(q_net/1000, digits=2))")
    println("  Reaction absorption:           -$(round(q_rxn/1000, digits=2))")
    println("  Conduction + storage:          +$(round(q_remaining/1000, digits=2))")

    # Surface cell concentrations
    if idx <= size(solution.ξ, 3)
        ξ_prof = solution.ξ[:, :, idx]
        if size(ξ_prof, 2) >= n_cells
            ξ_virgin = ξ_prof[1, n_cells]  # Surface cell
            ξ_char = ξ_prof[2, n_cells]
            ξ_gas = ξ_prof[3, n_cells]

            println("\nSurface cell composition [kg/m³]:")
            println("  Virgin: $(round(ξ_virgin, digits=2))")
            println("  Char: $(round(ξ_char, digits=4))")
            println("  Gas: $(round(ξ_gas, digits=4))")

            # Estimate reaction rate at surface
            A = 1.5e14
            E = 2.03e5
            R_gas = 8.314
            k_rate = A * exp(-E/(R_gas*T_surf))
            r_vol = k_rate * ξ_virgin

            println("\nReaction kinetics at surface:")
            println("  k(T_surf) = $(round(k_rate, sigdigits=4)) 1/s")
            println("  r = k × ξ_virgin = $(round(r_vol, digits=3)) kg/(m³·s)")
        end
    end
end

println("\n" * "="^70)
println("COMPARISON WITH THERMAKIN AT t ≈ 600s")
println("="^70)

# Find our values at t=600s
idx_600 = argmin(abs.(times .- 600.0))
if idx_600 <= length(times)
    println("\n                          Pyrolysis.jl    ThermaKin    Difference")
    println("-"^70)

    T_surf_PJ = T_surface[idx_600]
    T_surf_TK = 649.0
    println(
        "Surface T [K]             $(round(T_surf_PJ, digits=1))          $T_surf_TK        $(round(T_surf_TK - T_surf_PJ, digits=1))",
    )

    T_back_PJ = T_back[idx_600]
    T_back_TK = 507.2
    println(
        "Back T [K]                $(round(T_back_PJ, digits=1))          $T_back_TK        $(round(T_back_TK - T_back_PJ, digits=1))",
    )

    gen_rate_PJ = gas_generation_rate[idx_600] * 1000  # g/m²/s
    mlr_TK = 6.47  # ThermaKin reports actual MLR
    println(
        "Gen rate [g/(m²·s)]       $(round(gen_rate_PJ, digits=2))           $mlr_TK         $(round(mlr_TK - gen_rate_PJ, digits=2))",
    )

    thick_PJ = thickness_values[idx_600] * 1000  # mm
    thick_TK = 3.89
    println(
        "Thickness [mm]            $(round(thick_PJ, digits=2))           $thick_TK         $(round(thick_TK - thick_PJ, digits=2))",
    )

    println("-"^70)
    println("\nAnalysis:")
    println(
        "  ThermaKin shows HIGHER surface T (+$(round(T_surf_TK - T_surf_PJ, digits=0)) K)",
    )
    println(
        "  ThermaKin shows LOWER MLR (-$(round(gen_rate_PJ - mlr_TK, digits=2)) g/(m²·s))",
    )
    println()
    println("  This is physically consistent:")
    println("  - Higher T → more re-radiation loss → less heat for pyrolysis")
    println("  - OR: concentrated reaction zone → higher local T, lower total MLR")
    println()

    # Estimate reaction zone width
    A = 1.5e14
    E = 2.03e5
    R_gas = 8.314

    k_TK = A * exp(-E/(R_gas*T_surf_TK))
    k_PJ = A * exp(-E/(R_gas*T_surf_PJ))

    ξ_TK = 1192.8  # ThermaKin PMMA_glass concentration at surface
    ξ_PJ = idx_600 <= size(solution.ξ, 3) ? solution.ξ[:, :, idx_600][1, end] : 1200.0

    r_TK = k_TK * ξ_TK
    r_PJ = k_PJ * ξ_PJ

    # Reaction zone width (if uniform): δ = MLR / r
    δ_TK = mlr_TK/1000 / r_TK * 1000  # mm
    δ_PJ = gen_rate_PJ/1000 / r_PJ * 1000  # mm

    println("  Estimated reaction zone width:")
    println("    Pyrolysis.jl: $(round(δ_PJ, digits=2)) mm")
    println("    ThermaKin:    $(round(δ_TK, digits=2)) mm")
    println()
    println("  Pyrolysis.jl has WIDER reaction zone spreading heat absorption")
    println("  → Lower peak temperature, higher total MLR")
end

println("\n" * "="^70)
