"""
Cone Calorimeter Simulation of Pressure Treated Wood (PTW) with ALE Mesh Motion and Surface Absorption

This example simulates the pyrolysis of Pressure Treated Wood under standard cone
calorimeter conditions with ALE (Arbitrary Lagrangian-Eulerian) mesh motion
to track the material as it decomposes:
- 1D simulation through sample thickness
- ALE mesh motion for accurate shrinkage modeling
- Surface absorption radiation model for efficient opaque material simulation
- Multi-step pyrolysis kinetics with char formation
- Structured Jacobian (`Structured(strategy=DirectSolve())`) for fast assembly
- Bruggeman mixing rules for effective properties
- Temperature-limited reactions (piecewise activation)

Reference: Standard cone calorimeter test per ISO 5660-1

Physical Setup:
- Sample thickness: 12 mm (typical for wood samples)
- Exposed surface (z=L): radiant heating + convective/radiative losses + mass release
- Back surface (z=0): adiabatic + impermeable
- Surface regresses as material decomposes and gas escapes

Surface Absorption Radiation Model:
- All incident radiation deposited in the surface cell (z=L)
- Computationally efficient O(1) with tridiagonal Jacobian
- Appropriate for opaque materials like wood (optical thickness τ >> 1)
- Absorption coefficients (α) not used - radiation does not penetrate into depth
- Provides ~5-15× speedup over Beer-Lambert for opaque materials

Coordinate Convention:
- z = 0: Bottom (substrate/sample holder), boundary_id = 2, fixed
- z = L(t): Top (exposed surface), boundary_id = 1, moves with material
- Heat enters from top, material shrinks toward bottom
- Mesh follows the shrinking material (ALE formulation)

PTW Properties (from PTW.cmp):
- Virgin wood density: 498.46 kg/m³
- Multi-step decomposition with 4 intermediate states
- Final char density: 98.28 kg/m³
- Includes moisture content with evaporation
- Products: char + various gases + water vapor
- Temperature-limited reactions: reactions activate above T_limit

ThermaKin Temperature Limits and Piecewise Properties:
- ThermaKin uses "TEMP LIMIT: L T_limit" to prevent reactions below T_limit
- This creates discontinuities in effective properties at reaction onset
- In our model, we use temperature-dependent reaction rates (Arrhenius) naturally
- The temperature limits are implemented as piecewise activation thresholds
- Properties remain continuous but reaction contribution changes at T_limit

New Features Demonstrated:
- Surface absorption radiation model (radiation_model=SURFACE_ABSORPTION)
- Temperature-limited reactions from ThermaKin format
- Darcy-Fick gas transport (diffusion + pressure-driven flow)
- Temperature-dependent permeability evolution (virgin wood → porous char)
- ALE mesh motion with analytical Jacobian for performance
"""

using Pyrolysis
using Plots
using LinearSolve
using Sparspak
using OrdinaryDiffEq
using OrdinaryDiffEqSDIRK: KenCarp4   # v7 split: KenCarp4 not re-exported by the umbrella pkg
using Statistics
using LinearSolve: SparspakFactorization

# ============================================================================
# Material Definition: Pressure Treated Wood (PTW)
# ============================================================================

"""
Define PTW material with temperature-dependent properties.

PTW undergoes multi-step thermal decomposition:
	PTW → PTW1 → PTW2 → PTW3 → Char + gases

The material also contains moisture that evaporates:
	WATER → WATER_VAPOR

Heat capacity format from ThermaKin: cp = c0 + c1*T + c2*T^c3
For PTW components: cp = c0 + c1*T (linear, c2=0)

ThermaKin Temperature Limits:
- TEMP LIMIT: L 378 means reaction is inactive below 378 K
- TEMP LIMIT: L 290 for water evaporation (below 290 K)
- These create piecewise behavior as reactions turn on/off
- In our implementation, reactions are always active but naturally slow below T_limit
"""

# Temperature-dependent heat capacity functions
# Format from PTW.cmp: c0, c1, c2, c3 where cp = c0 + c1*T + c2*T^c3
# For PTW, all have c2=0, so cp = c0 + c1*T (linear)
ptw_cp = CallableProperty(T -> -302.0 + 4.93 * T)
ptw1_cp = CallableProperty(T -> -70.6 + 4.002 * T)
ptw2_cp = CallableProperty(T -> 160.8 + 3.0735 * T)
ptw3_cp = CallableProperty(T -> 392.2 + 2.145 * T)
char_cp = CallableProperty(T -> 623.6 + 1.217 * T)

# Water heat capacity: 5.2e3 - 6.7*T + 1.1e-2*T^2
water_cp = CallableProperty(T -> 5.2e3 - 6.7 * T + 1.1e-2 * T^2)

# Water vapor heat capacity: fitted to NIST Shomate data (300-1200K)
# NIST: 1865 (300K), 1955 (500K), 2150 (800K), 2291 (1000K) J/(kg·K)
water_vapor_cp = CallableProperty(T -> 1716.7 + 0.4154 * T + 1.530e-4 * T^2)

# ============================================================================
# Gas Transfer Coefficients and Permeability
# ============================================================================
#
# ThermaKin uses TRANSPORT parameter which combines both Fickian diffusion
# and Darcy flow. For PTW, TRANSPORT = 2e-5 m²/s for most components.
#
# We separate these into:
# - λ: Fickian gas transfer coefficient [m²/s]
# - κ: Darcy permeability [m²]
#
# Physical guidance for λ (from literature on wood pyrolysis):
# - Virgin wood: λ ~ 1e-6 m²/s (tight cellular structure)
# - Porous char: λ ~ 1e-5 m²/s (open pore network after volatiles escape)

# Base gas transfer coefficients for pyrolysis gases (Fickian diffusion)
λ_wood = 2.0e-5     # m²/s - matches ThermaKin TRANSPORT for PTW
λ_char = 2.0e-5     # m²/s - same for PTW components
λ_gas = 2.0e-5      # m²/s - gas phase diffusion

# Permeability parameters for Darcy flow (pressure-driven transport)
# PTW appears to have more uniform transport properties than Douglas Fir
κ_virgin = 1.0e-15  # m² - virgin PTW
κ_char = 1.0e-15    # m² - char (PTW has uniform transport)

# ============================================================================
# Component Indices for Reference
# ============================================================================
const WATER_IDX = 1            # Liquid water (moisture)
const WATER_VAPOR_IDX = 2      # Water vapor
const PTW_IDX = 3              # Virgin PTW
const PTW1_IDX = 4             # PTW intermediate 1
const PTW2_IDX = 5             # PTW intermediate 2
const PTW3_IDX = 6             # PTW intermediate 3
const CHAR_IDX = 7             # Final char
const PTW_GAS1_IDX = 8         # Gas product from reaction 2
const PTW_GAS2_IDX = 9         # Gas product from reaction 3
const PTW_GAS3_IDX = 10        # Gas product from reaction 4
const PTW_GAS4_IDX = 11        # Gas product from reaction 5

# Create PTW material with multiple components:
# Liquid: WATER
# Gas: WATER_VAPOR, PTW_gas1-4
# Solid: PTW, PTW1, PTW2, PTW3, Char

ptw_material = Material(
    name = :PressureTreatedWood,
    components = (
        # 1. Water (liquid moisture) - SWELLING: 0 means negligible volume
        # ThermaKin uses ρ=10000 kg/m³ as a "pseudo-density" for liquid water
        # This matches L SWELLING: 0 behavior where bound liquid doesn't contribute
        # to the material's bulk volume (absorbed into wood cell walls)
        LiquidComponent(
            :WATER,
            ρ = 1e5,                       # kg/m³ - pseudo-density (high to minimize volume)
            c = water_cp,                  # J/(kg·K) - temperature-dependent
            k = 0.145,                      # W/(m·K) - thermal conductivity
            ε = 0.81,                      # surface emissivity
            α = 1000.0,                    # m²/kg - absorption coefficient (not used with SURFACE_ABSORPTION)
            κ = 0.0,                        # m² - liquid blocks Darcy flow
        ),
        # 2. Water vapor - gaseous product from moisture evaporation
        GaseousComponent(
            :WATER_VAPOR,
            M = 0.018,                     # kg/mol - H₂O (18 g/mol)
            c = water_vapor_cp,            # J/(kg·K) - NIST polynomial fit
            k = 0.025,                     # W/(m·K) - H₂O vapor at ~400K
            λ = λ_gas * 1.6,               # m²/s - lightest molecule, fastest diffusion
        ),
        # 3. Virgin Pressure Treated Wood - solid phase
        SolidComponent(
            :PTW,
            ρ = 498.46,                    # kg/m³ - density
            c = ptw_cp,                    # J/(kg·K) - temperature-dependent
            k = 0.145,                      # W/(m·K) - thermal conductivity (virgin)
            ε = 0.81,                      # surface emissivity
            α = 1000.0,                    # m²/kg - not used with SURFACE_ABSORPTION model
            λ = λ_wood,                    # m²/s - gas transfer coefficient
            κ = κ_virgin,                   # m² - permeability
            # Cell-wall (skeletal) density of wood substance ~1500 kg/m³:
            # porosity φ = 1 - 498.46/1500 ≈ 0.67. Without it the pore space
            # has zero volume and cannot hold background air.
            ρ_intrinsic = 1500.0,
        ),
        # 4. PTW1 (first intermediate)
        SolidComponent(
            :PTW1,
            ρ = 364.25,
            c = ptw1_cp,
            k = 0.145,                      # Gradually increasing toward char
            ε = 0.783,                     # Interpolated toward char
            α = 1000.0,
            λ = λ_wood,
            κ = κ_virgin,
            # Skeletal density interpolated toward char (1950), same progress as ε
            ρ_intrinsic = 1610.0,
        ),
        # 5. PTW2 (second intermediate)
        SolidComponent(
            :PTW2,
            ρ = 156.37,
            c = ptw2_cp,
            k = 0.26,                       # Increasing conductivity
            ε = 0.755,                     # Continuing toward char
            α = 1000.0,
            λ = λ_wood,
            κ = κ_virgin,
            # Skeletal density interpolated toward char (1950), same progress as ε
            ρ_intrinsic = 1725.0,
        ),
        # 6. PTW3 (third intermediate, nearly char)
        SolidComponent(
            :PTW3,
            ρ = 127.83,
            c = ptw3_cp,
            k = 0.8,                        # Approaching char conductivity
            ε = 0.728,                     # Nearly char-like
            α = 1000.0,
            λ = λ_wood,
            κ = κ_virgin,
            # Skeletal density interpolated toward char (1950), same progress as ε
            ρ_intrinsic = 1835.0,
        ),
        # 7. Char - final solid residue
        SolidComponent(
            :Char,
            ρ = 98.28,
            c = char_cp,
            k = 1.0,                        # Char endpoint conductivity
            ε = 0.7,                       # Char emissivity
            α = 1000.0,
            λ = λ_char,
            κ = κ_char,
            # Carbon skeletal density ~1950 kg/m³ → char-only φ ≈ 0.95
            ρ_intrinsic = 1950.0,
        ),
        # 8. PTW_gas1 - gaseous product from reaction 2 (PTW → PTW1)
        # First decomposition step: hemicellulose-like breakdown
        # Products: acetic acid (60), furfural (96), formic acid (46), light tars
        GaseousComponent(
            :PTW_gas1,
            M = 0.060,                     # kg/mol - mix of acids and light organics
            c = 1500.0,                    # J/(kg·K) - organic vapor cp range
            k = 0.020,                     # W/(m·K) - organic vapors
            λ = λ_gas,
        ),
        # 9. PTW_gas2 - gaseous product from reaction 3 (PTW1 → PTW2)
        # Main cellulose decomposition (57% gas yield — largest gas-producing step)
        # Products: levoglucosan (162), anhydrosugars, heavy tars
        GaseousComponent(
            :PTW_gas2,
            M = 0.110,                     # kg/mol - dominated by levoglucosan/tars
            c = 1300.0,                    # J/(kg·K) - heavy organics have lower mass-specific cp
            k = 0.020,                     # W/(m·K)
            λ = λ_gas,
        ),
        # 10. PTW_gas3 - gaseous product from reaction 4 (PTW2 → PTW3)
        # Secondary decomposition (18% gas yield)
        # Products: lighter organics, phenolics from lignin
        GaseousComponent(
            :PTW_gas3,
            M = 0.080,                     # kg/mol - medium organics
            c = 1500.0,                    # J/(kg·K)
            k = 0.020,                     # W/(m·K)
            λ = λ_gas,
        ),
        # 11. PTW_gas4 - gaseous product from reaction 5 (PTW3 → Char)
        # Final charring stage (23% gas yield)
        # Products: lightest remaining volatiles, some CH₄, phenol fragments
        GaseousComponent(
            :PTW_gas4,
            M = 0.050,                     # kg/mol - lighter end of pyrolysis products
            c = 1600.0,                    # J/(kg·K) - lighter molecules, higher mass-specific cp
            k = 0.020,                     # W/(m·K)
            λ = λ_gas,
        ),
        # 12. Background pore air, initialized at ~1 atm partial pressure:
        # without it the pore space starts as a vacuum against the 101325 Pa
        # pressure boundary and early product gas is spuriously advected
        # inward down the absolute-pressure slope. Inert (no reactions).
        GaseousComponent(
            :AIR,
            M = 0.02897,                   # kg/mol - air
            c = 1006.0,                    # J/(kg·K) - air near 300-500 K
            k = 0.026,                     # W/(m·K) - air at ~300 K
            λ = λ_gas,                     # m²/s
        ),
    ),
    reactions = (
        # Reaction 1: WATER → WATER_VAPOR
        # ARRHENIUS: 15532.2434 43528.868 → A=15532.2434 s⁻¹, E=43528.868 J/mol
        # HEAT: -2.78E6 J/kg (negative = exothermic in ThermaKin)
        # TEMP LIMIT: L 290 → reaction active above 290 K
        # In our convention: h > 0 for endothermic, so h = +2.78e6
        Reaction(
            :water_evaporation,
            1 => 2,
            A = 15532.2434,                # s⁻¹
            E = 43528.868,                 # J/mol
            h = 2.78e6,                    # J/kg (endothermic in our convention)
            n = 1.0,
        ),
        # Reaction 2: PTW → 0.73075 PTW1 + 0.26925 PTW_gas1
        # ARRHENIUS: 19175814.8751 106116.6114
        # HEAT: -6.82E3 J/kg
        # TEMP LIMIT: L 378
        Reaction(
            :ptw_step1,
            3 => (4, 8),
            A = 19175814.8751,             # s⁻¹
            E = 106116.6114,               # J/mol
            h = 6.82e3,                    # J/kg (endothermic in our convention)
            n = 1.0,
            yields = (0.73075, 0.26925),
        ),
        # Reaction 3: PTW1 → 0.4293 PTW2 + 0.5707 PTW_gas2
        # ARRHENIUS: 186156937169.7363 161502.9426
        # HEAT: -1.67E5 J/kg
        # TEMP LIMIT: L 378
        Reaction(
            :ptw_step2,
            4 => (5, 9),
            A = 186156937169.7363,         # s⁻¹
            E = 161502.9426,               # J/mol
            h = 1.67e5,                    # J/kg (endothermic in our convention)
            n = 1.0,
            yields = (0.4293, 0.5707),
        ),
        # Reaction 4: PTW2 → 0.81744 PTW3 + 0.18256 PTW_gas3
        # ARRHENIUS: 179027.4062 100963.197
        # HEAT: 1.82E5 J/kg (positive in ThermaKin = endothermic)
        # TEMP LIMIT: L 378
        Reaction(
            :ptw_step3,
            5 => (6, 10),
            A = 179027.4062,               # s⁻¹
            E = 100963.197,                # J/mol
            h = -1.82e5,                   # J/kg (exothermic in our convention)
            n = 1.0,
            yields = (0.81744, 0.18256),
        ),
        # Reaction 5: PTW3 → 0.76882 Char + 0.23118 PTW_gas4
        # ARRHENIUS: 0.21061 34426.0053
        # HEAT: 2.02E5 J/kg (positive in ThermaKin = endothermic)
        # TEMP LIMIT: L 378
        Reaction(
            :ptw_step4,
            6 => (7, 11),
            A = 0.21061,                   # s⁻¹
            E = 34426.0053,                # J/mol
            h = -2.02e5,                   # J/kg (exothermic in our convention)
            n = 1.0,
            yields = (0.76882, 0.23118),
        ),
    ),
    # Mixing rules for effective properties (matches ThermaKin2Ds defaults)
    mixing = MaterialMixing(
        k = MixingSpec(WEIGHTED, 0.5),  # PARALL CONDUCTIVITY: 0.5
        λ = MixingSpec(WEIGHTED, 0.5),  # PARALL TRANSPORT: 0.5
        κ = MixingSpec(WEIGHTED, 0.5),  # permeability for Darcy flow
    ),
)

# ============================================================================
# Simulation Parameters
# ============================================================================

# Sample geometry
sample_thickness = 12.0e-3   # 12 mm typical wood sample
n_cells = 240             # 240 cells → 0.05 mm resolution

# Time parameters
t_final = 600.0             # 10 minutes - typical cone test duration for wood
dt_initial = 1e-3           # initial timestep

# Environmental conditions
T_ambient = 300.0           # K - ambient temperature (27°C)
q_external = 50e3  # W/m² - 50 kW/m² external heat flux

# Convective heat transfer coefficient for cone calorimeter
h_conv = 10.0                # W/(m²·K)

# ============================================================================
# Mesh Creation
# ============================================================================

# Number of components (NC = 12):
# 1 liquid (WATER) + 1 gas vapor (WATER_VAPOR) + 5 solids (PTW, PTW1-3, Char) + 5 gases (PTW_gas1-4, AIR)

# Initial concentrations:
# - PTW at full density (498.46 kg/m³)
# - Water moisture content (adjustable)
# - Background pore air at ~1 atm partial pressure
# - All others initially zero
moisture_content_db = 0.10  # 10% moisture content on dry basis
ρ_PTW_initial = 498.46      # Initial PTW density (dry)
# For dry basis, moisture density = MC_db * ρ_dry
ρ_water_initial = moisture_content_db * ρ_PTW_initial

# Pore air concentration giving P = 101325 Pa in the initial pore space
ξ_condensed_initial =
    (ρ_water_initial, 0.0, ρ_PTW_initial, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
ξ_air_initial =
    background_gas_concentration(ptw_material, :AIR, ξ_condensed_initial, T_ambient)

mesh = create_uniform_mesh(
    sample_thickness,      # domain length [m]
    n_cells,               # number of cells
    Val(12);               # number of components (NC)
    T_initial = T_ambient, # initial temperature
    ξ_initial = (
        ρ_water_initial,
        0.0,
        ρ_PTW_initial,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        0.0,
        ξ_air_initial,
    ),  # WATER + PTW + pore air, rest zero
)

println("Mesh created:")
println("  - Thickness: $(sample_thickness*1000) mm")
println("  - Number of cells: $n_cells")
println("  - Cell size: $(sample_thickness/n_cells*1000) mm")
println("  - Initial PTW density: $ρ_PTW_initial kg/m³")
println(
    "  - Initial water: $(round(ρ_water_initial, digits=1)) kg/m³ ($(round(moisture_content_db*100, digits=1))% dry basis)",
)

# Diagnostic: Check mass concentrations and fractions
total_initial_concentration = ρ_PTW_initial + ρ_water_initial
ptw_mass_fraction = ρ_PTW_initial / total_initial_concentration
water_mass_fraction = ρ_water_initial / total_initial_concentration

println(
    "  - Initial total concentration: $(round(total_initial_concentration, digits=1)) kg/m³",
)
println("  - PTW mass fraction: $(round(ptw_mass_fraction, digits=3))")
println("  - Water mass fraction: $(round(water_mass_fraction, digits=3))")
println(
    "  - Sum of mass fractions: $(round(ptw_mass_fraction + water_mass_fraction, digits=6)) (should be 1.0)",
)

# ============================================================================
# Boundary Conditions
# ============================================================================

# Top surface (boundary_id = 1): exposed surface
# - External radiant heat flux with surface absorption
# - Re-radiation losses (surface radiates to surroundings)
# - Convective losses (natural + forced convection)
top_thermal_bc =
    RadiativeFluxBC(flux = q_external) +                      # Cone heater (uses material ε)
    RadiativeBC(T_inf = 300) +                                # Re-radiation losses (uses material ε)
    ConvectiveBC(h = h_conv, T_inf = 300)                     # Convective losses

# Top surface mass BC: diffusive gas escape with Darcy flow to atmosphere
# h_m = 0.05 m/s matches ThermaKin typical boundary condition
top_mass_bc = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)

# Bottom surface (boundary_id = 2): back face
bottom_thermal_bc = AdiabaticBC()    # Adiabatic (no heat loss through substrate)
bottom_mass_bc = ImpermeableBC()   # Impermeable (no gas escape through substrate)

# Collect boundary conditions into BoundaryConditionSet
bc_set = BoundaryConditionSet(
    thermal = (top = top_thermal_bc, bottom = bottom_thermal_bc),
    mass = (top = top_mass_bc, bottom = bottom_mass_bc),
)

println("\nBoundary conditions:")
println(
    "  - Top surface: RadiativeFluxBC ($(q_external/1000) kW/m²) + re-radiation + convection",
)
println("  - Top mass: DiffusiveMassBC + DirichletPressureBC (101325 Pa)")
println("  - Bottom surface: Adiabatic, impermeable")

# ============================================================================
# Problem Setup
# ============================================================================

# Initial conditions
T_initial_func = z -> T_ambient           # Uniform initial temperature
ξ_initial_funcs = [
    z -> ρ_water_initial,  # 1. WATER (moisture)
    z -> 0.0,              # 2. WATER_VAPOR
    z -> ρ_PTW_initial,    # 3. PTW (virgin wood)
    z -> 0.0,              # 4. PTW1
    z -> 0.0,              # 5. PTW2
    z -> 0.0,              # 6. PTW3
    z -> 0.0,              # 7. Char
    z -> 0.0,              # 8. PTW_gas1
    z -> 0.0,              # 9. PTW_gas2
    z -> 0.0,              # 10. PTW_gas3
    z -> 0.0,              # 11. PTW_gas4
    z -> ξ_air_initial,    # 12. AIR (background pore air at ~1 atm)
]

println("Creating PyrolysisProblem...")
flush(stdout)

# Create the pyrolysis problem
problem = PyrolysisProblem(
    mesh = mesh,
    material = ptw_material,
    bc_set = bc_set,
    T_initial = T_initial_func,
    ξ_initial = ξ_initial_funcs,
    tspan = (0.0, t_final),
    dt_initial = dt_initial,
    output_times = collect(range(0.0, t_final, length = 3601)),  # Output every 1s
)

println("Problem created successfully!")
flush(stdout)

println("\nProblem created:")
println("  - Material: $(ptw_material.name)")
println("  - Duration: $(t_final) s")
println("  - External flux: $(q_external/1000) kW/m²")
println("  - Initial temperature: $(T_ambient) K")
println(
    "  - Conductivity mixing: $(ptw_material.mixing.k.rule) (β = $(ptw_material.mixing.k.β))",
)
println(
    "  - Gas transfer mixing: $(ptw_material.mixing.λ.rule) (β = $(ptw_material.mixing.λ.β))",
)
println("  - Darcy flow: ENABLED (κ_virgin = $(κ_virgin) m², κ_char = $(κ_char) m²)")

# ============================================================================
# Solve with ALE Mesh Motion and Surface Absorption
# ============================================================================

println(
    "\nSolving with ALE mesh motion, Darcy-Fick transport, surface absorption, and analytical Jacobian...",
)
println("Surface absorption model deposits all radiation in the surface cell.")
println("Darcy-Fick combines diffusion (λ-driven) + pressure-driven flow (κ-driven).")
println("This may take several minutes for wood pyrolysis...\n")

# Solve using implicit solver with ALE enabled and surface absorption radiation
solution = solve(
    problem;
    abstol = 1e-9,
    reltol = 1e-7,
    dtmin = 1e-14,                         # Allow very small timesteps during rapid reaction
    progress = true,
    use_ale = true,                        # Enable ALE mesh motion
    radiation_model = SURFACE_ABSORPTION,  # All radiation absorbed at surface (opaque material)
    min_thickness = 0.5e-3,                # Stop when sample is < 0.5 mm
    handle_depletion = false,              # Disable surface cell depletion (uses ALE instead)
    # Default Jacobian is Structured(strategy = DirectSolve()). Pass
    # `jacobian = JacobianSpec(Structured(strategy = ApproxSolve()), problem)`
    # to opt into the Schur-style approximation that drops the dense ALE
    # mesh-velocity coupling.
    integrator = KenCarp4(linsolve = SparspakFactorization()),
)

println("Solve() returned at $(time())!")
flush(stdout)

println("\nSolution completed!")
println("  - Return code: $(solution.retcode)")
println("  - Solve time: $(round(solution.solve_time, digits=2)) s")

# Thickness history
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

# Extract boundary fluxes (computed during solve for accuracy)
mlr = solution.extras.surface_mass_flux       # kg/(m²·s) - actual mass leaving at surface
q_absorbed = solution.extras.surface_heat_flux_absorbed      # W/m² - absorbed radiation
q_reradiation = solution.extras.surface_heat_flux_reradiation  # W/m² - re-radiation loss
q_convection = solution.extras.surface_heat_flux_convection    # W/m² - convective loss
q_conduction_surface = solution.extras.surface_heat_flux_conduction  # W/m² - conduction into material
q_net_surface = solution.extras.surface_heat_flux_net         # W/m² - net heat entering
ε_eff_front = solution.extras.surface_emissivity              # [-] - effective emissivity
k_eff_front = solution.extras.surface_conductivity            # W/(m·K) - effective conductivity

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
println("  - Peak MLR (surface flux): $(round(maximum(mlr)*1000, digits=2)) g/(m²·s)")
println(
    "  - Peak surface temperature: $(round(maximum(T_surface), digits=1)) K = $(round(maximum(T_surface)-273.15, digits=1)) °C",
)
println("  - Initial thickness: $(sample_thickness*1000) mm")
println("  - Final thickness: $(round(final_thickness*1000, digits=2)) mm")
println(
    "  - Thickness change: $(round((1 - final_thickness/sample_thickness)*100, digits=1))%",
)

println("\nBoundary Heat Fluxes (from solver):")
println("  - Peak absorbed radiation: $(round(maximum(q_absorbed)/1000, digits=2)) kW/m²")
println("  - Peak re-radiation loss: $(round(maximum(q_reradiation)/1000, digits=2)) kW/m²")
println("  - Peak convective loss: $(round(maximum(q_convection)/1000, digits=2)) kW/m²")
println(
    "  - Peak conductive flux: $(round(maximum(abs.(q_conduction_surface))/1000, digits=2)) kW/m²",
)

println("\nSurface Properties (from solver):")
println("  - Final emissivity: $(round(ε_eff_front[end], digits=3))")
println("  - Final conductivity: $(round(k_eff_front[end], digits=4)) W/(m·K)")

# ============================================================================
# Compute Component-Specific Gas Generation Rates (Post-Processing)
# ============================================================================

println("\nComputing component-specific gas generation rates...")

# Component indices for gas products
const GAS_COMPONENTS = [8, 9, 10, 11]  # PTW_gas1-4
const GAS_NAMES = ["PTW_gas1", "PTW_gas2", "PTW_gas3", "PTW_gas4"]

# Initialize component rate arrays
n_times = length(times)

# Gas generation rates (volumetric source terms)
mlr_gen_gas1 = zeros(n_times)      # From reaction 2
mlr_gen_gas2 = zeros(n_times)      # From reaction 3
mlr_gen_gas3 = zeros(n_times)      # From reaction 4
mlr_gen_gas4 = zeros(n_times)      # From reaction 5
mlr_gen_water = zeros(n_times)     # From reaction 1

# Pre-allocate arrays for reaction rate computation
NC = 12  # Number of components
NR = 5   # Number of reactions
rates = zeros(NR)
heats = zeros(NR)

# Compute component gas generation rates at each time
for i = 1:n_times
    ξ_profile = solution.ξ[:, :, i]
    T_profile = solution.T[:, i]
    n_cells_i = size(ξ_profile, 2)

    # Get volumes at this time (from thickness history for ALE)
    current_thickness = thickness_values[i]
    # Assume uniform cells for simplicity (good approximation for this case)
    cell_vol = current_thickness / n_cells_i * 1.0  # cross-section area = 1 m²

    # === Gas Generation Rates (Volumetric) ===
    for j = 1:n_cells_i
        T = T_profile[j]
        ξ = ntuple(k -> ξ_profile[k, j], NC)

        # Evaluate reaction rates
        Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, ptw_material)

        # Accumulate gas production from each reaction
        # Reaction 1: WATER → WATER_VAPOR (single product reaction)
        mlr_gen_water[i] += 1.0 * rates[1] * cell_vol  # yield = 1.0 for single product

        # Reaction 2: PTW → 0.73075 PTW1 + 0.26925 PTW_gas1
        mlr_gen_gas1[i] += ptw_material.reactions[2].stoich[4] * rates[2] * cell_vol

        # Reaction 3: PTW1 → 0.4293 PTW2 + 0.5707 PTW_gas2
        mlr_gen_gas2[i] += ptw_material.reactions[3].stoich[4] * rates[3] * cell_vol

        # Reaction 4: PTW2 → 0.81744 PTW3 + 0.18256 PTW_gas3
        mlr_gen_gas3[i] += ptw_material.reactions[4].stoich[4] * rates[4] * cell_vol

        # Reaction 5: PTW3 → 0.76882 Char + 0.23118 PTW_gas4
        mlr_gen_gas4[i] += ptw_material.reactions[5].stoich[4] * rates[5] * cell_vol
    end
end

println("\nGas Generation Rates (volumetric source):")
println("  - Peak H₂O generation: $(round(maximum(mlr_gen_water)*1000, digits=2)) g/(m²·s)")
println(
    "  - Peak PTW_gas1 generation: $(round(maximum(mlr_gen_gas1)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak PTW_gas2 generation: $(round(maximum(mlr_gen_gas2)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak PTW_gas3 generation: $(round(maximum(mlr_gen_gas3)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak PTW_gas4 generation: $(round(maximum(mlr_gen_gas4)*1000, digits=2)) g/(m²·s)",
)

# Find time to ignition (when MLR exceeds critical value)
mlr_critical = 3e-3  # kg/(m²·s) ≈ 3 g/(m²·s) - approximate ignition threshold
tig_idx = findfirst(mlr .> mlr_critical)
if tig_idx !== nothing
    t_ignition = times[tig_idx]
    println("\n  - Time to ignition (MLR > 3 g/m²s): $(round(t_ignition, digits=1)) s")
end

# ============================================================================
# Plotting Results
# ============================================================================

println("\nGenerating plots...")

# Create figure with 6 subplots
# Plot 1: Mass Loss Rate (MLR) - from solver (accurate)
p1 = plot(
    times,
    mlr .* 1000,  # Convert to g/(m²·s)
    xlabel = "Time [s]",
    ylabel = "MLR [g/(m²·s)]",
    title = "Mass Loss Rate (Surface Flux)",
    legend = :topright,
    linewidth = 2.5,
    color = :black,
    label = "Total MLR (solver)",
)
hline!([mlr_critical*1000], linestyle = :dot, color = :gray, label = "Ignition threshold")

# Plot 1b: Gas Generation Rates (volumetric source)
p1b = plot(
    times,
    gas_generation_rate .* 1000,
    xlabel = "Time [s]",
    ylabel = "Generation Rate [g/(m²·s)]",
    title = "Gas Generation Rate (Volumetric)",
    legend = :topright,
    linewidth = 2.5,
    color = :black,
    label = "Total (solver)",
)
# Add individual component generation rates
plot!(
    times,
    mlr_gen_water .* 1000,
    linewidth = 1.5,
    linestyle = :dot,
    label = "H₂O",
    color = :cyan,
)
plot!(
    times,
    mlr_gen_gas1 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "PTW_gas1",
    color = :blue,
)
plot!(
    times,
    mlr_gen_gas2 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "PTW_gas2",
    color = :green,
)
plot!(
    times,
    mlr_gen_gas3 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "PTW_gas3",
    color = :orange,
)
plot!(
    times,
    mlr_gen_gas4 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "PTW_gas4",
    color = :purple,
)

p2 = plot(
    times,
    T_surface,
    xlabel = "Time [s]",
    ylabel = "Temperature [K]",
    title = "Surface & Back Temperature",
    legend = :bottomright,
    linewidth = 2,
    color = :orange,
    label = "Surface",
)
plot!(times, T_back, linewidth = 2, color = :blue, label = "Back")

p3 = plot(
    times,
    total_mass .* 1000,  # Convert to g/m²
    xlabel = "Time [s]",
    ylabel = "Mass [g/m²]",
    title = "Total Mass",
    legend = false,
    linewidth = 2,
    color = :green,
    ylims = (0, :auto),
)

# Sample thickness over time
p4 = plot(
    thickness_times,
    thickness_values .* 1000,  # Convert to mm
    xlabel = "Time [s]",
    ylabel = "Thickness [mm]",
    title = "Sample Thickness",
    legend = false,
    linewidth = 2,
    color = :purple,
)

# Temperature profiles at selected times
p5 = plot(xlabel = "Depth [mm]", ylabel = "Temperature [K]", title = "Temperature Profiles")
profile_times = [30, 60, 120, 180, 300, 450, 600, 750, 900, 1000]  # seconds
for t in profile_times
    # Find closest output time
    idx = argmin(abs.(times .- t))
    if idx <= size(solution.T, 2)
        T_profile = solution.T[:, idx]
        # Skip if only 1 cell remains
        if length(T_profile) <= 1
            continue
        end
        # Use current thickness at this time
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

# Combine all plots
final_plot = plot(p1, p1b, p2, p3, p4, p5, layout = (2, 3), size = (1400, 800))
plot!(plot_title = "Cone Calorimeter Simulation: PTW at 50 kW/m² (ALE Mesh)")

# Save plot
savefig(final_plot, "cone_calorimeter_ptw_results.png")
println("Results saved to: cone_calorimeter_ptw_results.png")

# ============================================================================
# Component Evolution Plot
# ============================================================================

println("\nGenerating component evolution plot...")

# Extract component concentrations over time (averaged across thickness)
n_times = size(solution.ξ, 3)
avg_water = zeros(n_times)
avg_PTW = zeros(n_times)
avg_PTW1 = zeros(n_times)
avg_PTW2 = zeros(n_times)
avg_PTW3 = zeros(n_times)
avg_char = zeros(n_times)

for i = 1:n_times
    ξ_profile = solution.ξ[:, :, i]
    n_cells_i = size(ξ_profile, 2)
    if n_cells_i > 0
        avg_water[i] = sum(ξ_profile[1, :]) / n_cells_i
        avg_PTW[i] = sum(ξ_profile[3, :]) / n_cells_i
        avg_PTW1[i] = sum(ξ_profile[4, :]) / n_cells_i
        avg_PTW2[i] = sum(ξ_profile[5, :]) / n_cells_i
        avg_PTW3[i] = sum(ξ_profile[6, :]) / n_cells_i
        avg_char[i] = sum(ξ_profile[7, :]) / n_cells_i
    end
end

p_components = plot(
    times,
    avg_PTW,
    xlabel = "Time [s]",
    ylabel = "Concentration [kg/m³]",
    title = "Solid Component Evolution",
    label = "PTW (virgin)",
    linewidth = 2,
    legend = :right,
)
plot!(times, avg_PTW1, label = "PTW1", linewidth = 1.5)
plot!(times, avg_PTW2, label = "PTW2", linewidth = 1.5)
plot!(times, avg_PTW3, label = "PTW3", linewidth = 1.5)
plot!(times, avg_char, label = "Char", linewidth = 2, color = :black)
plot!(times, avg_water, label = "Water", linewidth = 1.5, linestyle = :dash)

savefig(p_components, "ptw_component_evolution.png")
println("Component evolution saved to: ptw_component_evolution.png")

# ============================================================================
# CSV Output
# ============================================================================

println("\nWriting CSV output...")

# Write CSV manually with key diagnostic quantities
open("ptw_results_50kWm2.csv", "w") do io
    # Header
    header = [
        "time",
        "mass",
        "mlr",
        "thickness",
        "front_temp",
        "back_temp",
        "PTW",
        "PTW1",
        "PTW2",
        "PTW3",
        "Char",
        "Water",
        "k_eff_front",
        "eps_eff_front",
        "q_absorbed",
        "q_reradiation",
        "q_convection",
        "q_conduction_surface",
        "q_net_surface",
        "mlr_gen_water",
        "mlr_gen_gas1",
        "mlr_gen_gas2",
        "mlr_gen_gas3",
        "mlr_gen_gas4",
    ]
    println(io, join(header, ","))

    # Data rows
    for i = 1:length(times)
        data = [
            times[i],
            total_mass[i],
            mlr[i],
            thickness_values[i],
            T_surface[i],
            T_back[i],
            avg_PTW[i],
            avg_PTW1[i],
            avg_PTW2[i],
            avg_PTW3[i],
            avg_char[i],
            avg_water[i],
            k_eff_front[i],
            ε_eff_front[i],
            q_absorbed[i],
            q_reradiation[i],
            q_convection[i],
            q_conduction_surface[i],
            q_net_surface[i],
            mlr_gen_water[i],
            mlr_gen_gas1[i],
            mlr_gen_gas2[i],
            mlr_gen_gas3[i],
            mlr_gen_gas4[i],
        ]
        println(io, join(data, ","))
    end
end

println("CSV output saved to: ptw_results_50kWm2.csv")

# ============================================================================
# Spatial Data Output
# ============================================================================

println("\nWriting spatial data CSV...")

# Select times for spatial output
spatial_output_times = [
    0.0,
    5.0,
    10.0,
    15.0,
    20.0,
    30.0,
    40.0,
    50.0,
    60.0,
    80.0,
    100.0,
    150.0,
    200.0,
    300.0,
    400.0,
    500.0,
    600.0,
    700.0,
    800.0,
    900.0,
    1000.0,
]

open("ptw_spatial_50kWm2.csv", "w") do io
    # Header
    header = [
        "time_s",
        "from_back_m",
        "from_front_m",
        "temperature_K",
        "conc_WATER_kg_m3",
        "conc_WATER_VAPOR_kg_m3",
        "conc_PTW_kg_m3",
        "conc_PTW1_kg_m3",
        "conc_PTW2_kg_m3",
        "conc_PTW3_kg_m3",
        "conc_Char_kg_m3",
        "conc_PTW_gas1_kg_m3",
        "conc_PTW_gas2_kg_m3",
        "conc_PTW_gas3_kg_m3",
        "conc_PTW_gas4_kg_m3",
        "conc_AIR_kg_m3",
        "total_density_kg_m3",
    ]
    println(io, join(header, ","))

    for target_t in spatial_output_times
        # Find closest time index
        idx = argmin(abs.(times .- target_t))
        actual_t = times[idx]

        # Skip if too far from target
        if abs(actual_t - target_t) > 1.0 && target_t > 0
            continue
        end

        # Get spatial data at this time
        T_profile = solution.T[:, idx]
        ξ_profile = solution.ξ[:, :, idx]
        n_cells_i = length(T_profile)

        # Get current thickness
        t_idx = argmin(abs.(thickness_times .- actual_t))
        current_thickness = thickness_values[t_idx]

        # Cell positions (from back, z=0 at back)
        dz = current_thickness / n_cells_i

        for cell_idx = 1:n_cells_i
            from_back = (cell_idx - 0.5) * dz  # Cell center position
            from_front = current_thickness - from_back

            T_cell = T_profile[cell_idx]
            ξ_cell = [ξ_profile[k, cell_idx] for k = 1:12]
            total_density = sum(ξ_cell)

            data = [
                actual_t,
                from_back,
                from_front,
                T_cell,
                ξ_cell[1],
                ξ_cell[2],
                ξ_cell[3],
                ξ_cell[4],
                ξ_cell[5],
                ξ_cell[6],
                ξ_cell[7],
                ξ_cell[8],
                ξ_cell[9],
                ξ_cell[10],
                ξ_cell[11],
                ξ_cell[12],
                total_density,
            ]
            println(io, join(data, ","))
        end
    end
end

println("Spatial data saved to: ptw_spatial_50kWm2.csv")

println("\nSimulation complete!")
println("\nFeatures used in this example:")
println("  ✓ WEIGHTED mixing rule for thermal conductivity (β=0.5, matches ThermaKin)")
println("  ✓ WEIGHTED mixing rule for gas transfer coefficient (β_λ=0.5)")
println("  ✓ Darcy-Fick gas transport (diffusion + pressure-driven flow)")
println("  ✓ Temperature-dependent permeability (κ = $(κ_virgin) m²)")
println("  ✓ ALE mesh motion with analytical Jacobian")
println("  ✓ Solver-computed boundary fluxes (MLR, heat flux decomposition)")
println("  ✓ Piecewise property handling from ThermaKin TEMP LIMIT parameters")
