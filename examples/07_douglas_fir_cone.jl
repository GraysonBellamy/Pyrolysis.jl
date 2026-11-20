"""
Cone Calorimeter Simulation of Douglas Fir with ALE Mesh Motion and Surface Absorption

This example simulates the pyrolysis of Douglas Fir wood under standard cone 
calorimeter conditions with ALE (Arbitrary Lagrangian-Eulerian) mesh motion 
to track the material as it decomposes:
- 1D simulation through sample thickness
- ALE mesh motion for accurate shrinkage modeling
- Surface absorption radiation model for efficient opaque material simulation
- Multi-step pyrolysis kinetics with char formation
- Bound water in two-way exchange with its pore vapour (a phase change)
- Sparse Jacobian by colored forward-mode differentiation, with the default
  integrator (Rodas5P with a sparse KLU solve)
- Blended parallel/series mixing rules for effective properties
- Component gas-transfer coefficients mixed into one shared face coefficient

Reference: Standard cone calorimeter test per ISO 5660-1

Physical Setup:
- Sample thickness: 12 mm (typical for wood samples)
- Exposed surface (z=L): radiant heating + convective/radiative losses + mass release
- Back surface (z=0): adiabatic + impermeable
- Surface regresses as material decomposes and gas escapes

Surface Absorption Radiation Model:
- The absorbed portion of incident radiation is deposited at the exposed face (z=L)
- Appropriate for opaque materials like wood (optical thickness τ >> 1)
- Absorption coefficients (α) not used - radiation does not penetrate into depth

Coordinate Convention:
- z = 0: Bottom (substrate/sample holder), tag :bottom, fixed
- z = L(t): Top (exposed surface), tag :top, moves with material
- Heat enters from top, material shrinks toward bottom
- Mesh follows the shrinking material (ALE formulation)

Douglas Fir Properties (from DF.cmp):
- Virgin wood density: 530 kg/m³
- Five decomposition steps through four intermediate solids to char
- Final char density: 150 kg/m³
- 6% moisture (dry basis), bound in the wood
- Products: char + various gases + water vapor

Moisture:
- The bound water and its pore vapour exchange mass toward sorption
  equilibrium, p_v → a(X, T) p_sat(T), with the activity a of the CAST
  isotherm for spruce (CAST_SPRUCE). The water evaporates where the wood
  heats and condenses where the vapour driven inward meets cooler wood, so
  the moisture rises ahead of the drying front and the back face holds near
  the temperature at which the bound water's equilibrium pressure meets the
  pore pressure.
- The virgin wood and its four intermediate solids hold water; char holds none.
- The heat of the step is Kirchhoff's latent heat from the water's and the
  vapour's heat capacities plus the isotherm's heat of sorption.
- The pore vapour starts in equilibrium with the initial moisture, and the
  ambient air carries the humidity in equilibrium with it at 300 K (30.7 %),
  so the conditioned sample is at rest before the heater acts.

New Features Demonstrated:
- Two-way phase change of the moisture and its vapour (PhaseChange)
- Surface absorption radiation model (radiation_model=SURFACE_ABSORPTION)
- Gas transfer coefficients of 5e-7 m²/s in the virgin wood, where λ stands for
  the bound water's transfer across the grain (Technical Reference §9.2.1), and
  1e-5 m²/s in the char, an estimate
- lumped wall and Darcy gas transport (diffusion + pressure-driven flow)
- Composition-dependent permeability evolution (virgin wood → porous char)
- Species molar masses in the ideal-gas pore pressure (Darcy)
- ALE mesh motion
"""

using Pyrolysis
using Plots
using Statistics

# ============================================================================
# Material Definition: Douglas Fir
# ============================================================================

"""
Define Douglas Fir material with temperature-dependent properties.

Douglas Fir undergoes multi-step thermal decomposition:
	DF → DF_int1 → DF_int2 → DF_int3 → DF_int4 → Char + gases

The material also holds moisture, which evaporates and condenses:
	MOISTURE ⇌ WATER_VAPOR

Heat capacity format from ThermaKin: cp = c0 + c1*T + c2*T^c3
For most components: cp = c0 + c1*T (linear)
"""

# Temperature-dependent heat capacity functions
# Virgin DF: piecewise linear from HFM (low-T) and STA (high-T) regressions
# Segment 1 (T ≤ 323.16 K): HFM regression
# Segment 2 (323.16–425.02 K): linearly interpolated between seg1 and seg3 endpoints
# Segment 3 (T > 425.02 K): STA per-curve regression
const T_HFM_MAX = 323.16
const T_STA_MIN = 425.02219679025956
const T_INTERP_WIDTH = T_STA_MIN - T_HFM_MAX

# Virgin cp regression coefficients [J/(kg·K)]
const DF_CP_A1 = 4.133836600273531     # slope, low-T (HFM)
const DF_CP_B1 = -18.01630709338176    # intercept, low-T (HFM)
const DF_CP_A3 = 9.136799235573225     # slope, high-T (STA)
const DF_CP_B3 = -2139.3044552435705   # intercept, high-T (STA)

# Char cp regression coefficients [J/(kg·K)] (independent char STA)
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

df_cp = CallableProperty(virgin_cp)
char_cp = CallableProperty(T -> CHAR_CP_B + CHAR_CP_A * T)

# Intermediate cp's derived via cumulative-yield interpolation between virgin and char
# Reaction yields: 0.95, 0.80, 0.37, 0.79, 0.74
# Final char yield: 0.95 * 0.80 * 0.37 * 0.79 * 0.74 = 0.164388...
const FINAL_CHAR_YIELD = 0.95 * 0.80 * 0.37 * 0.79 * 0.74
const INTERP_PROGRESS = [
    (1 - 0.95) / (1 - FINAL_CHAR_YIELD),                        # int1: α ≈ 0.0598
    (1 - 0.95 * 0.80) / (1 - FINAL_CHAR_YIELD),                 # int2: α ≈ 0.2872
    (1 - 0.95 * 0.80 * 0.37) / (1 - FINAL_CHAR_YIELD),          # int3: α ≈ 0.8602
    (1 - 0.95 * 0.80 * 0.37 * 0.79) / (1 - FINAL_CHAR_YIELD),   # int4: α ≈ 0.9309
]

df_int1_cp = CallableProperty(
    T ->
        (1 - INTERP_PROGRESS[1]) * virgin_cp(T) +
        INTERP_PROGRESS[1] * (CHAR_CP_B + CHAR_CP_A * T),
)
df_int2_cp = CallableProperty(
    T ->
        (1 - INTERP_PROGRESS[2]) * virgin_cp(T) +
        INTERP_PROGRESS[2] * (CHAR_CP_B + CHAR_CP_A * T),
)
df_int3_cp = CallableProperty(
    T ->
        (1 - INTERP_PROGRESS[3]) * virgin_cp(T) +
        INTERP_PROGRESS[3] * (CHAR_CP_B + CHAR_CP_A * T),
)
df_int4_cp = CallableProperty(
    T ->
        (1 - INTERP_PROGRESS[4]) * virgin_cp(T) +
        INTERP_PROGRESS[4] * (CHAR_CP_B + CHAR_CP_A * T),
)
# Moisture cp: Shomate equation (Chase, 1998), Cp = (A + B*t + C*t^2 + D*t^3 + E/t^2) / M [J/(kg·K)]
moisture_cp = CallableProperty(T -> begin
    t = T / 1000.0
    M = 0.01801528  # kg/mol (H2O)
    A = -203.6060
    B = 1523.290
    C = -3196.413
    D = 2474.455
    E = 3.855326
    (A + B*t + C*t^2 + D*t^3 + E/t^2) / M
end)
# Water vapor heat capacity: fitted to NIST Shomate data (300-1200K)
# NIST: 1865 (300K), 1955 (500K), 2150 (800K), 2291 (1000K) J/(kg·K)
water_vapor_cp = CallableProperty(T -> 1716.7 + 0.4154 * T + 1.530e-4 * T^2)

# ============================================================================
# Gas Transfer Coefficients and Permeability (Physically Realistic)
# ============================================================================
# 
# Superficial wall transfer is J_j = -M_j Λ/(R T) ∂p_j/∂z.
# The bulk occupancy of the transport matrix (the wood and char components)
# weights the effective Λ. These coefficients specify the pore matrix; gas
# components have caloric properties and molar masses but no transfer
# coefficient. Moisture declares neither λ nor κ: it occupies skeletal volume
# but is no part of the transport matrix.
λ_wood = 5.0e-7     # m²/s - superficial transfer through virgin wood
λ_char = 1.0e-5     # m²/s - superficial transfer through char

# Permeability for Darcy flow (pressure-driven bulk transport through solid matrix)
# κ determines resistance to pressure-driven gas flow, separate from lumped wall diffusion.
# Virgin softwood (transverse): κ ~ 1e-16 to 1e-14 m²
# Char (open pore network): κ ~ 1e-12 to 1e-11 m²
κ_virgin = 1.0e-15  # m² - tight virgin wood cellular structure
κ_char = 1.0e-12    # m² - open char pore network

# Create Douglas Fir material with multiple components:
# Condensed components: DF, DF_int1, DF_int2, DF_int3, DF_int4, Char, MOISTURE
# Gas components: DF_gas1-5, WATER_VAPOR

douglas_fir = Material(
    name = :DouglasFir,
    components = (
        # 1. Virgin Douglas Fir - solid phase
        SolidComponent(
            :DF,
            ρ = 530.0,                     # kg/m³ - density
            c = df_cp,                     # J/(kg·K) - temperature-dependent
            k = 0.145,                      # W/(m·K) - thermal conductivity (virgin)
            ε = 0.81,                      # surface emissivity
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            λ = λ_wood,                    # m²/s - tight virgin structure
            κ = κ_virgin,                   # m² - low permeability
            # Cell-wall (skeletal) density of wood substance ~1500 kg/m³
            # (Siau): porosity φ = 1 - 530/1500 ≈ 0.65. Without it the pore
            # space has zero volume and cannot hold background air.
            ρ_intrinsic = 1500.0,
        ),
        # 2. Intermediate 1 (partially decomposed, structure beginning to open)
        SolidComponent(
            :DF_int1,
            ρ = 507.2,
            c = df_int1_cp,
            # Thermal conductivity interpolated toward char using downstream solid-yield product
            # Progress to char = 0.80*0.37*0.79*0.74 ≈ 0.1730416
            # k = k_virgin + (k_char - k_virgin) * progress = 0.145 + (1.0-0.145)*progress
            # ≈ 0.145 + 0.855*0.1730416 ≈ 0.292951
            k = 0.292951,
            # Emissivity interpolated toward char based on downstream solid-yield chain
            # DF_int1 → DF_int2(0.80) → DF_int3(0.37) → DF_int4(0.79) → Char(0.74)
            # Progress to char = 0.80*0.37*0.79*0.74 ≈ 0.1730416
            ε = 0.790965,
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            λ = λ_wood * 2.0,             # structure beginning to open (ε/τ increasing)
            κ = κ_virgin * 3.0,            # m² - slightly more permeable
            # Skeletal density interpolated toward char (1950) by progress 0.173
            ρ_intrinsic = 1578.0,
        ),
        # 3. Intermediate 2 (significant decomposition, pore network developing)
        SolidComponent(
            :DF_int2,
            ρ = 419.8,
            c = df_int2_cp,
            # Progress to char = 0.37*0.79*0.74 ≈ 0.216302
            # k = 0.145 + (1.0-0.145)*0.216302 ≈ 0.329938
            k = 0.329938,
            # DF_int2 → DF_int3(0.37) → DF_int4(0.79) → Char(0.74)
            # Progress to char = 0.37*0.79*0.74 ≈ 0.216302
            ε = 0.786207,
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            λ = λ_wood * 5.0,             # pore network developing (ε~0.15, τ~3)
            κ = κ_virgin * 10.0,           # m² - pores opening up
            # Skeletal density interpolated toward char (1950) by progress 0.216
            ρ_intrinsic = 1597.0,
        ),
        # 4. Intermediate 3 (advanced decomposition, transitioning to char)
        SolidComponent(
            :DF_int3,
            ρ = 203.2,
            c = df_int3_cp,
            # Progress to char = 0.79*0.74 ≈ 0.5846
            # k = 0.145 + (1.0-0.145)*0.5846 ≈ 0.644833
            k = 0.644833,
            # DF_int3 → DF_int4(0.79) → Char(0.74)
            # Progress to char = 0.79*0.74 ≈ 0.5846
            ε = 0.745694,
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            λ = λ_char * 0.3,             # transitioning to char (ε~0.4, τ~2.5)
            κ = κ_char * 0.01,             # m² - significant porosity
            # Skeletal density interpolated toward char (1950) by progress 0.585
            ρ_intrinsic = 1763.0,
        ),
        # 5. Intermediate 4 (nearly char, highly porous)
        SolidComponent(
            :DF_int4,
            ρ = 176.6,
            c = df_int4_cp,
            # Progress to char = 0.74
            # k = 0.145 + (1.0-0.145)*0.74 ≈ 0.777700
            k = 0.777700,
            # DF_int4 → Char(0.74)
            # Progress to char = 0.74
            ε = 0.728600,
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            λ = λ_char * 0.7,             # nearly char-like (ε~0.55, τ~2)
            κ = κ_char * 0.1,              # m² - approaching char permeability
            # Skeletal density interpolated toward char (1950) by progress 0.74
            ρ_intrinsic = 1833.0,
        ),
        # 6. Char - final solid residue (fully porous). Its bulk density sets
        # the residue thickness (Technical Reference §10.2): with the chain's
        # solid yield of 0.164, the slab chars to 0.58 of its dry thickness,
        # the measured residue thickness.
        SolidComponent(
            :Char,
            ρ = 150.0,
            c = char_cp,
            k = 1.0,                      # Char has higher conductivity (char endpoint)
            ε = 0.70,                     # Char emissivity
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            λ = λ_char,                   # m²/s - fully porous char
            κ = κ_char,                     # m² - open pore network
            # Carbon skeletal density ~1950 kg/m³ → char-only φ ≈ 0.92
            ρ_intrinsic = 1950.0,
        ),
        # 7. Moisture: water bound in the cell walls, a liquid at the density
        # of water. It holds its own volume, so drying shrinks the wood and
        # sorption swells it, and it adds to k_eff and ε_eff by its volume
        # share.
        LiquidComponent(
            :MOISTURE,
            ρ = 1000.0,                    # kg/m³ - liquid water
            c = moisture_cp,               # Moisture contributes to heat capacity
            k = 0.6,                       # W/(m·K) - liquid water
            ε = 0.81,
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            κ = 0.0,                        # m² - liquid blocks Darcy flow
        ),
        # 8. DF_gas1 - gaseous product from reaction 1 (5% yield, low-T onset)
        # Early hemicellulose decomposition: acetic acid (60), furfural (96), formic acid (46)
        GaseousComponent(
            :DF_gas1,
            M = 0.060,                     # kg/mol - acetic acid-like (CH₃COOH = 60 g/mol)
            c = 1400.0,                    # J/(kg·K) - organic acids/aldehydes
            k = 0.020,                     # W/(m·K) - organic vapors
            # m²/s - ~2.1e-5 (moderate molecule)
        ),
        # 9. DF_gas2 - gaseous product from reaction 2 (20% yield)
        # Continued hemicellulose / early cellulose: mix of organic fragments, some CO₂
        GaseousComponent(
            :DF_gas2,
            M = 0.070,                     # kg/mol - mix of organics heavier than pure CO₂
            c = 1400.0,                    # J/(kg·K) - organic fragments
            k = 0.020,                     # W/(m·K)
            # m²/s - ~1.9e-5 (heavier organics)
        ),
        # 10. DF_gas3 - gaseous product from reaction 3 (63% yield — dominant gas-producing step)
        # Main cellulose pyrolysis: levoglucosan (162), anhydrosugars, heavy tars
        GaseousComponent(
            :DF_gas3,
            M = 0.110,                     # kg/mol - dominated by levoglucosan/anhydrosugars
            c = 1300.0,                    # J/(kg·K) - heavy organics, lower mass-specific cp
            k = 0.018,                     # W/(m·K) - heavy organic vapors
            # m²/s - ~1.5e-5 (heavy tars, slowest)
        ),
        # 11. DF_gas4 - gaseous product from reaction 4 (21% yield)
        # Secondary decomposition: lighter organics, phenol fragments from lignin
        GaseousComponent(
            :DF_gas4,
            M = 0.070,                     # kg/mol - phenolics (94), guaiacol (124), lighter organics
            c = 1500.0,                    # J/(kg·K) - medium organics
            k = 0.020,                     # W/(m·K)
            # m²/s - ~1.9e-5
        ),
        # 12. DF_gas5 - gaseous product from reaction 5 (26% yield, final charring)
        # Last volatiles released during charring: CH₄, light hydrocarbons, residual phenolics
        GaseousComponent(
            :DF_gas5,
            M = 0.050,                     # kg/mol - mix of CH₄ (16), phenol fragments, light organics
            c = 1600.0,                    # J/(kg·K) - lighter molecules, higher mass-specific cp
            k = 0.025,                     # W/(m·K)
            # m²/s - ~2.3e-5 (lighter mix)
        ),
        # 13. Water vapor - the moisture's vapour, made and taken back by the
        # phase change below (lightest, fastest)
        GaseousComponent(
            :WATER_VAPOR,
            M = 0.018,                     # kg/mol - H₂O (18 g/mol, lightest species)
            c = water_vapor_cp,            # J/(kg·K) - NIST polynomial fit
            k = 0.025,                     # W/(m·K) - H₂O vapor at ~400K
            # m²/s - ~3.8e-5 (lightest molecule, fastest)
        ),
        # 14. Background pore air, initialized at ~1 atm partial pressure:
        # without it the pore space starts as a vacuum against the 101325 Pa
        # pressure boundary and early product gas is spuriously advected
        # inward down the absolute-pressure slope. Inert (no reactions).
        GaseousComponent(
            :AIR,
            M = 0.02897,                   # kg/mol - air
            c = 1006.0,                    # J/(kg·K) - air near 300-500 K
            k = 0.026,                     # W/(m·K) - air at ~300 K
            # m²/s - reference M~29 diffusivity
        ),
    ),
    reactions = (
        # Reaction 1: DF → 0.95 DF_int1 + 0.05 DF_gas1
        Reaction(
            :pyrolysis_step1,
            1 => (2, 8),
            A = 74816.95,                  # s⁻¹
            E = 71.59e3,                   # J/mol
            h = -19200.66288032251,        # J/kg (exothermic)
            n = 1.0,
            yields = (0.95, 0.05),
        ),
        # Reaction 2: DF_int1 → 0.80 DF_int2 + 0.20 DF_gas2
        Reaction(
            :pyrolysis_step2,
            2 => (3, 9),
            A = 23550492839.0,             # s⁻¹
            E = 140.9e3,                   # J/mol
            h = -5263.201466602739,        # J/kg (exothermic)
            n = 1.0,
            yields = (0.80, 0.20),
        ),
        # Reaction 3: DF_int2 → 0.37 DF_int3 + 0.63 DF_gas3
        Reaction(
            :pyrolysis_step3,
            3 => (4, 10),
            A = 7.744618e12,               # s⁻¹
            E = 181.0e3,                   # J/mol
            h = 264424.79037373835,        # J/kg (endothermic)
            n = 1.0,
            yields = (0.37, 0.63),
        ),
        # Reaction 4: DF_int3 → 0.79 DF_int4 + 0.21 DF_gas4
        Reaction(
            :pyrolysis_step4,
            4 => (5, 11),
            A = 103276.14,                 # s⁻¹
            E = 98.94e3,                   # J/mol
            h = -217151.69043762444,       # J/kg (exothermic)
            n = 1.0,
            yields = (0.79, 0.21),
        ),
        # Reaction 5: DF_int4 → 0.74 Char + 0.26 DF_gas5
        Reaction(
            :pyrolysis_step5,
            5 => (6, 12),
            A = 6.397348,                  # s⁻¹
            E = 56.01e3,                   # J/mol
            h = -326.40885469229025,       # J/kg (exothermic)
            n = 1.0,
            yields = (0.74, 0.26),
        ),
    ),
    # MOISTURE ⇌ WATER_VAPOR: the pore vapour relaxes toward the equilibrium
    # pressure a(X, T) p_sat(T) of the bound water at the rate ω = K (ξ_eq − ξ_v),
    # K = 10³ s⁻¹ by default. X is the water over the sorbing solids: the wood
    # and its intermediates, not the char. The heat per kg evaporated is
    # 2.4417 MJ/kg at 298.15 K plus the vapour's sensible enthalpy less the
    # water's, plus the isotherm's heat of sorption.
    phase_changes = (
        PhaseChange(
            :moisture,
            :MOISTURE => :WATER_VAPOR;
            activity = CAST_SPRUCE,
            sorbents = (:DF, :DF_int1, :DF_int2, :DF_int3, :DF_int4),
        ),
    ),
    # Example mixing choices; TK has similar blends but different volume weights
    mixing = MaterialMixing(
        k = MixingSpec(WEIGHTED, 0.5),  # 50/50 parallel/series
        λ = MixingSpec(WEIGHTED, 0.5),  # weighted average for gas transfer
        κ = MixingSpec(WEIGHTED, 0.5),  # permeability for Darcy flow
    ),
)

# ============================================================================
# Simulation Parameters
# ============================================================================

# Sample geometry
sample_thickness = 12.0e-3   # 12 mm typical wood sample
n_cells = get(ENV, "PYROLYSIS_SMOKE", "0") == "1" ? 24 : 240  # Full run: 0.05 mm

# Time parameters
const SMOKE = get(ENV, "PYROLYSIS_SMOKE", "0") == "1"   # PYROLYSIS_SMOKE=1: a short run that checks the script works
t_final = SMOKE ? 20.0 : 600.0  # 10 minutes - typical cone test duration for wood
dt_initial = 1e-3           # initial timestep

# Environmental conditions
T_ambient = 300.0           # K - ambient temperature (27°C)
q_external = 75e3  # W/m² - external heat flux
# Output files carry the flux, e.g. douglas_fir_results_75kWm2_fit.csv.
results_csv = "douglas_fir_results_$(round(Int, q_external / 1000))kWm2_fit.csv"
spatial_csv = "douglas_fir_spatial_$(round(Int, q_external / 1000))kWm2_fit.csv"

# Convective heat transfer coefficient for cone calorimeter
h_conv = 10.0                # W/(m²·K)

# ============================================================================
# Mesh Creation
# ============================================================================

# Number of components (NC = 14):
# 6 solids (DF, DF_int1-4, Char) + 1 liquid (MOISTURE) + 7 gases (DF_gas1-5, WATER_VAPOR, AIR)

# Initial concentrations:
# - Bound water swells the wood by its own volume, so the dry wood (oven-dry
#   density 530 kg/m³) and the water share the cell: ξ_DF/530 + ξ_w/1000 = 1
#   with ξ_w = MC ξ_DF, giving ξ_DF = 513.7 and ξ_w = 30.8 kg/m³ at 6% MC
# - Pore vapour in equilibrium with the bound water, and pore air making up
#   the pressure to 101325 Pa
# - All others initially zero
moisture_content_db = 0.06  # 6% moisture content on dry basis
ρ_DF_oven_dry = 530.0       # Oven-dry density: oven-dry mass over oven-dry volume
ρ_water = 1000.0            # Density of the bound water (MOISTURE component)
ρ_DF_initial = ρ_DF_oven_dry / (1 + moisture_content_db * ρ_DF_oven_dry / ρ_water)
ρ_moisture_initial = moisture_content_db * ρ_DF_initial

# The pore gas at the start: the vapour in equilibrium with the bound water at
# T_ambient, and the air at the rest of 101325 Pa. Both helpers read the
# condensed entries of ξ_condensed_initial only (the porosity follows from the
# condensed composition and intrinsic densities).
ξ_condensed_initial = (
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
    0.0,
)
p_vapour_initial = equilibrium_vapour_pressure(douglas_fir, ξ_condensed_initial, T_ambient)
ξ_vapour_initial =
    equilibrium_vapour_concentration(douglas_fir, ξ_condensed_initial, T_ambient)
ξ_air_initial = background_gas_concentration(
    douglas_fir,
    :AIR,
    ξ_condensed_initial,
    T_ambient;
    P = 101325.0 - p_vapour_initial,
)

# The initial temperature and concentrations (DF + moisture + pore air, rest
# zero) come from the problem's T_initial and ξ_initial below.
mesh = uniform_mesh(
    sample_thickness,      # domain length [m]
    n_cells,               # number of cells
)

println("Mesh created:")
println("  - Thickness: $(sample_thickness*1000) mm")
println("  - Number of cells: $n_cells")
println("  - Cell size: $(sample_thickness/n_cells*1000) mm")
println(
    "  - Initial dry wood: $(round(ρ_DF_initial, digits=1)) kg/m³ (oven-dry density $ρ_DF_oven_dry kg/m³)",
)
println(
    "  - Initial moisture: $(round(ρ_moisture_initial, digits=1)) kg/m³ ($(round(moisture_content_db*100, digits=1))% dry basis)",
)
println(
    "  - Initial pore vapour: $(round(p_vapour_initial, digits=1)) Pa, in equilibrium with the moisture",
)

# Diagnostic: Check mass concentrations and fractions
total_initial_concentration = ρ_DF_initial + ρ_moisture_initial
wood_mass_fraction = ρ_DF_initial / total_initial_concentration
moisture_mass_fraction = ρ_moisture_initial / total_initial_concentration

println(
    "  - Initial total concentration: $(round(total_initial_concentration, digits=1)) kg/m³",
)
println("  - Wood mass fraction: $(round(wood_mass_fraction, digits=3))")
println("  - Moisture mass fraction: $(round(moisture_mass_fraction, digits=3))")
println(
    "  - Sum of mass fractions: $(round(wood_mass_fraction + moisture_mass_fraction, digits=6)) (should be 1.0)",
)

# ============================================================================
# Boundary Conditions
# ============================================================================

# Top surface (:top): exposed surface
# - External radiant heat flux, absorbed at the surface
# - Re-radiation losses (surface radiates to surroundings)
# - Convective losses (natural + forced convection)
# Note: RadiativeFluxBC and RadiativeBC use the effective emissivity from material
# properties by default (ε = nothing), which tracks the surface emissivity
# as the material decomposes (e.g., virgin wood → char transition). To override with
# a fixed value, pass ε = 0.9 or absorptivity = 0.9 explicitly.
# Note: radiation_model=SURFACE_ABSORPTION controls surface vs in-depth absorption
top_thermal_bc =
    RadiativeFluxBC(q_external) +                      # Cone heater (uses material ε)
    RadiativeBC(T_inf = 300) +                                # Re-radiation losses (uses material ε)
    ConvectiveBC(h = h_conv, T_inf = 300)                     # Convective losses

# Top surface mass BC: diffusive gas escape with Darcy flow to atmosphere
# The reservoir supplies ingress composition and temperature. Its humidity is
# the one that leaves the conditioned sample at rest: the vapour mole fraction
# in equilibrium with the initial moisture at the ambient temperature,
# x_v = p_eq(X₀, T∞)/P (30.7 % relative humidity at 300 K). The mass film is
# the twin of the convective film above, with its conductance from h_conv by
# the heat–mass transfer analogy (the default AnalogyFilm), and carries the
# vapour both ways. Interior wall and Darcy transport share the species-flux
# continuity closure.
x_vapour_ambient = p_vapour_initial / 101325.0
top_mass_bc = AmbientGasBC(
    GasReservoir(
        P = 101325.0,
        T = 300.0,
        mole_fractions = (AIR = 1.0 - x_vapour_ambient, WATER_VAPOR = x_vapour_ambient),
    ),
)

# Bottom surface (:bottom): back face
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
println(
    "  - Top mass: humid air at 101325 Pa (water vapour mole fraction $(round(x_vapour_ambient, digits=5))) with a species film",
)
println("  - Bottom surface: Adiabatic, impermeable")

# ============================================================================
# Problem Setup
# ============================================================================

# Initial conditions
T_initial_func = z -> T_ambient           # Uniform initial temperature
ξ_initial_funcs = [
    z -> ρ_DF_initial,     # 1. DF (virgin wood)
    z -> 0.0,              # 2. DF_int1
    z -> 0.0,              # 3. DF_int2
    z -> 0.0,              # 4. DF_int3
    z -> 0.0,              # 5. DF_int4
    z -> 0.0,              # 6. Char
    z -> ρ_moisture_initial, # 7. MOISTURE
    z -> 0.0,              # 8. DF_gas1
    z -> 0.0,              # 9. DF_gas2
    z -> 0.0,              # 10. DF_gas3
    z -> 0.0,              # 11. DF_gas4
    z -> 0.0,              # 12. DF_gas5
    z -> ξ_vapour_initial, # 13. WATER_VAPOR (in equilibrium with the moisture)
    z -> ξ_air_initial,    # 14. AIR (background pore air, the rest of 1 atm)
]

println("Creating PyrolysisProblem...")
flush(stdout)

# Create the pyrolysis problem
problem = PyrolysisProblem(
    mesh = mesh,
    material = douglas_fir,
    bc_set = bc_set,
    T_initial = T_initial_func,
    ξ_initial = ξ_initial_funcs,
    tspan = (0.0, t_final),
    dt_initial = dt_initial,
    output_times = collect(range(0.0, t_final, length = SMOKE ? 21 : 3601)),  # Output every 1s
)

println("Problem created successfully!")
flush(stdout)

println("\nProblem created:")
println("  - Material: $(douglas_fir.name)")
println("  - Duration: $(t_final) s")
println("  - External flux: $(q_external/1000) kW/m²")
println("  - Initial temperature: $(T_ambient) K")
println(
    "  - Conductivity mixing: $(douglas_fir.mixing.k.rule) (β = $(douglas_fir.mixing.k.β))",
)
println(
    "  - Gas transfer mixing: $(douglas_fir.mixing.λ.rule) (β = $(douglas_fir.mixing.λ.β))",
)
println("  - Darcy flow: ENABLED (κ_virgin = $(κ_virgin) m², κ_char = $(κ_char) m²)")
println(
    "  - Moisture: phase change MOISTURE ⇌ WATER_VAPOR, CAST spruce isotherm, K = $(douglas_fir.phase_changes[1].exchange_rate) s⁻¹",
)

# ============================================================================
# Solve with ALE Mesh Motion and Surface Absorption
# ============================================================================

println(
    "\nSolving with ALE mesh motion, lumped wall and Darcy transport and surface absorption...",
)
println("Surface absorption model absorbs all radiation at the exposed face.")
println(
    "lumped wall and Darcy combines diffusion (λ-driven) + pressure-driven flow (κ-driven).",
)
println("This may take several minutes for wood pyrolysis...\n")

# Solve using implicit solver with ALE enabled and surface absorption radiation.
# The default integrator (Rodas5P with a sparse KLU solve) takes the Jacobian
# from colored forward-mode differentiation over the mesh sparsity pattern.
solution = solve(
    problem;
    abstol = 1e-9,
    reltol = 1e-7,
    dtmin = 1e-14,                         # Allow very small timesteps during rapid reaction
    progress = true,
    use_ale = true,                        # Enable ALE mesh motion
    radiation_model = SURFACE_ABSORPTION,  # All radiation absorbed at surface (opaque material)
    min_thickness = 0.5e-3,                # Stop when sample is < 0.5 mm
)

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
#
# The solution provides:
# - solution.extras.surface_mass_flux: Total mass leaving at surface (computed during solve)
# - solution.extras.gas_generation_rate: Total volumetric gas production (computed during solve)
# - solution.extras.evaporation_rate, condensation_rate: the phase change's
#   gross exchange; their difference is the net production of water vapour
#
# The rates of the pyrolysis gases still require post-processing to separate
# them by species.
#
# ============================================================================

println("\nComputing component-specific gas generation rates...")

# Component indices for gas products
const GAS_COMPONENTS = [8, 9, 10, 11, 12, 13]  # DF_gas1-5, WATER_VAPOR
const GAS_NAMES = ["DF_gas1", "DF_gas2", "DF_gas3", "DF_gas4", "DF_gas5", "H₂O"]

# Initialize component rate arrays
n_times = length(times)

# Gas generation rates (volumetric source terms)
mlr_gen_gas1 = zeros(n_times)      # From reaction 1
mlr_gen_gas2 = zeros(n_times)      # From reaction 2
mlr_gen_gas3 = zeros(n_times)      # From reaction 3
mlr_gen_gas4 = zeros(n_times)      # From reaction 4
mlr_gen_gas5 = zeros(n_times)      # From reaction 5

# Water vapour from the phase change, net: evaporation less condensation
evaporation = solution.extras.evaporation_rate      # kg/(m²·s)
condensation = solution.extras.condensation_rate    # kg/(m²·s)
mlr_gen_water = evaporation .- condensation

# Pre-allocate arrays for reaction rate computation
NC = 14  # Number of components
NR = 5   # Number of reactions
rates = zeros(NR)
heats = zeros(NR)

# Compute component gas generation rates at each time
for i = 1:n_times
    ξ_profile = solution.ξ[:, :, i]
    T_profile = solution.T[:, i]
    n_cells_i = size(ξ_profile, 2)

    # Integrate over the actual moving cells, per unit initial exposed area.
    nodes = get_profile(solution, :z, i)
    area_ratio =
        solution.extras.cross_section_area_history[i] /
        solution.extras.cross_section_area_history[1]
    cell_depths = diff(nodes) .* area_ratio

    # === Gas Generation Rates (Volumetric) ===
    for j = 1:n_cells_i
        isfinite(T_profile[j]) || continue
        cell_vol = cell_depths[j]
        T = T_profile[j]
        ξ = ntuple(k -> ξ_profile[k, j], NC)

        # Evaluate reaction rates
        Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, douglas_fir)

        # Accumulate gas production from each reaction
        # Reaction 1: DF → 0.95 DF_int1 + 0.05 DF_gas1
        mlr_gen_gas1[i] += douglas_fir.reactions[1].products[2][2] * rates[1] * cell_vol

        # Reaction 2: DF_int1 → 0.80 DF_int2 + 0.20 DF_gas2
        mlr_gen_gas2[i] += douglas_fir.reactions[2].products[2][2] * rates[2] * cell_vol

        # Reaction 3: DF_int2 → 0.37 DF_int3 + 0.63 DF_gas3
        mlr_gen_gas3[i] += douglas_fir.reactions[3].products[2][2] * rates[3] * cell_vol

        # Reaction 4: DF_int3 → 0.79 DF_int4 + 0.21 DF_gas4
        mlr_gen_gas4[i] += douglas_fir.reactions[4].products[2][2] * rates[4] * cell_vol

        # Reaction 5: DF_int4 → 0.74 Char + 0.26 DF_gas5
        mlr_gen_gas5[i] += douglas_fir.reactions[5].products[2][2] * rates[5] * cell_vol
    end
end

# Water evaporated and condensed over the run [kg/m²], by the trapezoid rule
# over the saved times
trapezoid(y) = sum((times[k+1] - times[k]) * (y[k] + y[k+1]) / 2 for k = 1:(n_times-1))
water_evaporated = trapezoid(evaporation)
water_condensed = trapezoid(condensation)

println("\nGas Generation Rates (volumetric source):")
println(
    "  - Peak DF_gas1 generation: $(round(maximum(mlr_gen_gas1)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak DF_gas2 generation: $(round(maximum(mlr_gen_gas2)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak DF_gas3 generation: $(round(maximum(mlr_gen_gas3)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak DF_gas4 generation: $(round(maximum(mlr_gen_gas4)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak DF_gas5 generation: $(round(maximum(mlr_gen_gas5)*1000, digits=2)) g/(m²·s)",
)
println(
    "  - Peak H₂O generation (net): $(round(maximum(mlr_gen_water)*1000, digits=2)) g/(m²·s)",
)

println("\nMoisture (phase change):")
println(
    "  - Initial water: $(round(ρ_moisture_initial * sample_thickness, digits=4)) kg/m²",
)
println("  - Water evaporated: $(round(water_evaporated, digits=4)) kg/m²")
println(
    "  - Water condensed: $(round(water_condensed, digits=4)) kg/m² ($(round(100 * water_condensed / water_evaporated, digits=1))% of the evaporated)",
)
println("  - Peak evaporation: $(round(maximum(evaporation)*1000, digits=2)) g/(m²·s)")

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
    mlr_gen_gas1 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "DF_gas1",
    color = :blue,
)
plot!(
    times,
    mlr_gen_gas2 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "DF_gas2",
    color = :green,
)
plot!(
    times,
    mlr_gen_gas3 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "DF_gas3",
    color = :orange,
)
plot!(
    times,
    mlr_gen_gas4 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "DF_gas4",
    color = :purple,
)
plot!(
    times,
    mlr_gen_gas5 .* 1000,
    linewidth = 1.5,
    linestyle = :dash,
    label = "DF_gas5",
    color = :brown,
)
plot!(
    times,
    mlr_gen_water .* 1000,
    linewidth = 1.5,
    linestyle = :dot,
    label = "H₂O (net)",
    color = :cyan,
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
        # Cell centres of the moving mesh at this time
        z_coords = solution.cell_centers[:, idx] .* 1000
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
plot!(
    plot_title = "Cone Calorimeter Simulation: Douglas Fir at $(round(Int, q_external / 1000)) kW/m² (ALE Mesh)",
)

# Save plot
savefig(final_plot, "cone_calorimeter_douglas_fir_results.png")
println("Results saved to: cone_calorimeter_douglas_fir_results.png")

# ============================================================================
# Component Evolution Plot
# ============================================================================

println("\nGenerating component evolution plot...")

# Thickness average of component k at saved time i. The cells of the moving
# mesh have unequal widths, so each concentration is weighted by its cell's.
function thickness_average(solution, k, i)
    widths = diff(get_profile(solution, :z, i))
    live = findall(isfinite, solution.T[:, i])
    return sum(solution.ξ[k, j, i] * widths[j] for j in live) / sum(widths[live])
end

# Extract component concentrations over time (averaged across thickness)
n_times = size(solution.ξ, 3)
avg_DF = zeros(n_times)
avg_int1 = zeros(n_times)
avg_int2 = zeros(n_times)
avg_int3 = zeros(n_times)
avg_int4 = zeros(n_times)
avg_char = zeros(n_times)
avg_moisture = zeros(n_times)

for i = 1:n_times
    ξ_profile = solution.ξ[:, :, i]
    n_cells_i = size(ξ_profile, 2)
    if n_cells_i > 0
        avg_DF[i] = thickness_average(solution, 1, i)
        avg_int1[i] = thickness_average(solution, 2, i)
        avg_int2[i] = thickness_average(solution, 3, i)
        avg_int3[i] = thickness_average(solution, 4, i)
        avg_int4[i] = thickness_average(solution, 5, i)
        avg_char[i] = thickness_average(solution, 6, i)
        avg_moisture[i] = thickness_average(solution, 7, i)
    end
end

p_components = plot(
    times,
    avg_DF,
    xlabel = "Time [s]",
    ylabel = "Concentration [kg/m³]",
    title = "Solid Component Evolution",
    label = "DF (virgin)",
    linewidth = 2,
    legend = :right,
)
plot!(times, avg_int1, label = "DF_int1", linewidth = 1.5)
plot!(times, avg_int2, label = "DF_int2", linewidth = 1.5)
plot!(times, avg_int3, label = "DF_int3", linewidth = 1.5)
plot!(times, avg_int4, label = "DF_int4", linewidth = 1.5)
plot!(times, avg_char, label = "Char", linewidth = 2, color = :black)
plot!(times, avg_moisture, label = "Moisture", linewidth = 1.5, linestyle = :dash)

savefig(p_components, "douglas_fir_component_evolution.png")
println("Component evolution saved to: douglas_fir_component_evolution.png")

# ============================================================================
# Compute Additional Diagnostics (Beyond Solver-Computed Quantities)
# ============================================================================
#
# The solver provides these quantities computed during the solve:
# - solution.extras.surface_mass_flux: Total MLR (already extracted as `mlr`)
# - solution.extras.surface_heat_flux_absorbed: Absorbed radiation (already `q_absorbed`)
# - solution.extras.surface_heat_flux_reradiation: Re-radiation loss (already `q_reradiation`)
# - solution.extras.surface_heat_flux_convection: Convective loss (already `q_convection`)
# - solution.extras.surface_heat_flux_conduction: Conduction into material (already `q_conduction_surface`)
# - solution.extras.surface_heat_flux_net: Net heat entering (already `q_net_surface`)
# - solution.extras.surface_emissivity: Effective emissivity (already `ε_eff_front`)
# - solution.extras.surface_conductivity: Effective conductivity (already `k_eff_front`)
#
# This section computes ADDITIONAL diagnostics that require post-processing:
# - Mass balance checks
# - Component concentrations and averages
# - Temperature at specific depths
# - Energy balance and reaction heat rates
#
# ============================================================================

println("\nComputing additional diagnostics...")

# Pre-allocate diagnostic arrays for quantities NOT in solver output
n_times = length(times)

# Mass balance diagnostics
total_solid_conc = zeros(n_times)
total_liquid_conc = zeros(n_times)
total_gas_conc = zeros(n_times)
sum_all_components = zeros(n_times)
bulk_density = zeros(n_times)
mass_balance_error = zeros(n_times)

# Gas concentrations at front surface
gas1_conc_front = zeros(n_times)
gas2_conc_front = zeros(n_times)
gas3_conc_front = zeros(n_times)
gas4_conc_front = zeros(n_times)
gas5_conc_front = zeros(n_times)
water_vapor_conc_front = zeros(n_times)

# Average gas concentrations
gas1_conc_avg = zeros(n_times)
gas2_conc_avg = zeros(n_times)
gas3_conc_avg = zeros(n_times)
gas4_conc_avg = zeros(n_times)
gas5_conc_avg = zeros(n_times)
water_vapor_conc_avg = zeros(n_times)

# Additional effective properties (not in solver output)
cp_eff_front = zeros(n_times)
rho_cp_front = zeros(n_times)
alpha_thermal = zeros(n_times)

# Temperature at specific depths (for spatial comparison)
T_25pct = zeros(n_times)  # Temperature at 25% depth from surface
T_50pct = zeros(n_times)  # Temperature at 50% depth (midplane)
T_75pct = zeros(n_times)  # Temperature at 75% depth from surface

# Energy balance
energy_stored = zeros(n_times)         # Total thermal energy in sample [J/m²]
energy_rate = zeros(n_times)           # dE/dt [W/m²]
total_rxn_heat_domain = zeros(n_times) # Total reaction and phase-change heat (integrated over domain) [W/m²]

# Reaction rate diagnostics at front surface
rxn1_rate_front = zeros(n_times)
rxn2_rate_front = zeros(n_times)
rxn3_rate_front = zeros(n_times)
rxn4_rate_front = zeros(n_times)
rxn5_rate_front = zeros(n_times)
phase_change_rate_front = zeros(n_times)   # net evaporation [kg/(m³·s)]
total_rxn_heat_front = zeros(n_times)

# The phase change's net rate ω [kg/(m³·s)] and its heat h_pc [J/kg] in a cell
# of concentrations ξ at temperature T and pore fraction φ, as the solver
# evaluates them: the sorbing solids at their positive part, and the heat from
# the sensible enthalpies the solver builds from the water's and the vapour's
# heat capacities.
moisture_change = only(douglas_fir.phase_changes)
water_enthalpy = Pyrolysis.States.enthalpy_law(
    douglas_fir.components[moisture_change.liquid].heat_capacity,
)
vapour_enthalpy = Pyrolysis.States.enthalpy_law(
    douglas_fir.components[moisture_change.vapour].heat_capacity,
)
function moisture_exchange(ξ, T, φ)
    pc = moisture_change
    X = Pyrolysis.Physics.moisture_content(pc, ξ, map(x -> max(x, 0.0), ξ))
    ω = Pyrolysis.Physics.phase_change_rate(pc, X, T, φ, ξ[pc.vapour])
    h_v = Pyrolysis.States.sensible_enthalpy(vapour_enthalpy, T)
    h_w = Pyrolysis.States.sensible_enthalpy(water_enthalpy, T)
    return ω, Pyrolysis.Physics.phase_change_heat(pc, X, T, h_v, h_w)
end

# Compute diagnostics for each timestep
for i = 1:n_times
    ξ_profile = solution.ξ[:, :, i]
    T_profile = solution.T[:, i]
    n_cells_i = size(ξ_profile, 2)

    if n_cells_i > 0
        # === Mass Balance Diagnostics (averaged over domain) ===
        solids =
            avg_DF[i] + avg_int1[i] + avg_int2[i] + avg_int3[i] + avg_int4[i] + avg_char[i]
        liquids = avg_moisture[i]

        # Average gas concentrations
        gas1_conc_avg[i] = thickness_average(solution, 8, i)
        gas2_conc_avg[i] = thickness_average(solution, 9, i)
        gas3_conc_avg[i] = thickness_average(solution, 10, i)
        gas4_conc_avg[i] = thickness_average(solution, 11, i)
        gas5_conc_avg[i] = thickness_average(solution, 12, i)
        water_vapor_conc_avg[i] = thickness_average(solution, 13, i)
        air_conc_avg = thickness_average(solution, 14, i)

        gases =
            gas1_conc_avg[i] +
            gas2_conc_avg[i] +
            gas3_conc_avg[i] +
            gas4_conc_avg[i] +
            gas5_conc_avg[i] +
            water_vapor_conc_avg[i] +
            air_conc_avg

        total_solid_conc[i] = solids
        total_liquid_conc[i] = liquids
        total_gas_conc[i] = gases
        sum_all_components[i] = solids + liquids + gases
        bulk_density[i] = total_mass[i] / thickness_values[i]
        mass_balance_error[i] = (sum_all_components[i] - bulk_density[i]) / bulk_density[i]

        # === Front Surface Diagnostics ===
        # Front surface is the last cell (z = L, tag :top)
        T_front = T_profile[end]
        ξ_front = ntuple(k -> ξ_profile[k, end], 14)

        # Gas concentrations at front
        gas1_conc_front[i] = ξ_front[8]
        gas2_conc_front[i] = ξ_front[9]
        gas3_conc_front[i] = ξ_front[10]
        gas4_conc_front[i] = ξ_front[11]
        gas5_conc_front[i] = ξ_front[12]
        water_vapor_conc_front[i] = ξ_front[13]

        # Compute additional effective properties at front surface
        # Total density at front
        ρ_front = sum(ξ_front)

        # Effective heat capacity: cp_eff = Σ(ξⱼ·cpⱼ) / ρ
        cp_numerator = 0.0
        for j = 1:14
            cp_j = douglas_fir.components[j].heat_capacity(T_front)
            cp_numerator += ξ_front[j] * cp_j
        end
        cp_eff_front[i] = cp_numerator / ρ_front
        rho_cp_front[i] = ρ_front * cp_eff_front[i]

        # Thermal diffusivity (using solver-computed k_eff_front)
        alpha_thermal[i] = k_eff_front[i] / rho_cp_front[i]

        # === Temperature at Specific Depths ===
        # Depth measured from surface (z=L): the cell whose centre is nearest
        # to z = 0.75 L (25% depth), 0.50 L (midplane) and 0.25 L (75% depth)
        centres = solution.cell_centers[:, i]
        nearest(z) = argmin(abs.(replace(centres, NaN => Inf) .- z))
        idx_25pct = nearest(0.75 * thickness_values[i])
        idx_50pct = nearest(0.50 * thickness_values[i])
        idx_75pct = nearest(0.25 * thickness_values[i])
        T_25pct[i] = T_profile[idx_25pct]
        T_50pct[i] = T_profile[idx_50pct]
        T_75pct[i] = T_profile[idx_75pct]

        # === Energy Balance ===
        # Total thermal energy: E = ∫ ρ cp T dV (per unit area)
        # E = Σ (ρ_i × cp_i × T_i × Δz_i)
        widths = diff(get_profile(solution, :z, i))   # cell widths of the moving mesh
        energy_stored[i] = 0.0
        total_rxn_heat_domain[i] = 0.0
        for cell_idx = 1:n_cells_i
            T_cell = T_profile[cell_idx]
            ξ_cell = ntuple(k -> ξ_profile[k, cell_idx], 14)

            # ρcp for this cell
            ρ_cell = sum(ξ_cell)
            cp_cell = 0.0
            for j = 1:14
                cp_j = douglas_fir.components[j].heat_capacity(T_cell)
                cp_cell += ξ_cell[j] * cp_j / ρ_cell
            end

            # Contribution to thermal energy (relative to 298 K reference)
            energy_stored[i] += ρ_cell * cp_cell * (T_cell - 298.15) * widths[cell_idx]

            # Reaction heat at this cell: Q = -Σ(hᵢ × rᵢ) [W/m³]
            # heats[i] is in J/kg, rates[i] is in kg/(m³·s)
            # So heats[i] * rates[i] = W/m³ (volumetric heating rate)
            # Note: In our convention, h < 0 means exothermic (releases heat)
            # So Q = -Σ(h × r) gives positive Q for exothermic reactions
            local_rates = zeros(NR)
            local_heats = zeros(NR)
            Pyrolysis.Physics.evaluate_reactions!(
                local_rates,
                local_heats,
                ξ_cell,
                T_cell,
                douglas_fir,
            )
            # Sum of (-h × r) × dz gives W/m² for each cell contribution
            for rxn_idx = 1:NR
                total_rxn_heat_domain[i] +=
                    -local_heats[rxn_idx] * local_rates[rxn_idx] * widths[cell_idx]
            end
            # The phase change takes h_pc per kg evaporated and returns it on
            # condensation: -h_pc × ω
            ω_cell, h_pc_cell =
                moisture_exchange(ξ_cell, T_cell, solution.gas.porosity[cell_idx, i])
            total_rxn_heat_domain[i] += -h_pc_cell * ω_cell * widths[cell_idx]
        end

        # Energy rate: dE/dt (finite difference)
        if i > 1 && times[i] > times[i-1]
            dt = times[i] - times[i-1]
            energy_rate[i] = (energy_stored[i] - energy_stored[i-1]) / dt
        end

        # === Reaction Rates at Front Surface ===
        global rates = zeros(NR)
        global heats = zeros(NR)
        Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ_front, T_front, douglas_fir)

        rxn1_rate_front[i] = rates[1]
        rxn2_rate_front[i] = rates[2]
        rxn3_rate_front[i] = rates[3]
        rxn4_rate_front[i] = rates[4]
        rxn5_rate_front[i] = rates[5]
        ω_front, h_pc_front =
            moisture_exchange(ξ_front, T_front, solution.gas.porosity[end, i])
        phase_change_rate_front[i] = ω_front

        # Total reaction and phase-change heat at front (W/m³)
        # Q = -Σ(hᵢ × rᵢ) - h_pc × ω, where h < 0 = exothermic gives Q > 0
        total_rxn_heat_front[i] = -h_pc_front * ω_front
        for rxn_idx = 1:NR
            total_rxn_heat_front[i] += -heats[rxn_idx] * rates[rxn_idx]
        end
    end
end

println(
    "  - Maximum mass balance error: $(round(maximum(abs.(mass_balance_error))*100, digits=2))%",
)
println(
    "  - Maximum gas concentration at surface: $(round(maximum([maximum(gas1_conc_front), maximum(gas2_conc_front), maximum(gas3_conc_front), maximum(gas4_conc_front), maximum(gas5_conc_front), maximum(water_vapor_conc_front)]), digits=4)) kg/m³",
)
println("  - Average k_eff at surface: $(round(mean(k_eff_front), digits=4)) W/(m·K)")
println(
    "  - Average net heat flux at surface: $(round(mean(q_net_surface)/1000, digits=2)) kW/m²",
)

# ============================================================================
# CSV Output with Comprehensive Diagnostics
# ============================================================================

println("\nWriting comprehensive CSV output...")

# Write CSV manually with all diagnostic quantities
open(results_csv, "w") do io
    # Header - organized in logical groups
    header = [
        # Time and basic state
        "time",
        "mass",
        "mlr",
        "thickness",
        "front_temp",
        "back_temp",

        # Temperature at specific depths (from surface)
        "T_25pct_depth",
        "T_50pct_depth",
        "T_75pct_depth",

        # Solid components (averaged)
        "DF",
        "DF_int1",
        "DF_int2",
        "DF_int3",
        "DF_int4",
        "Char",
        "Moisture",

        # Mass balance diagnostics
        "total_solid_conc",
        "total_liquid_conc",
        "total_gas_conc",
        "sum_all_components",
        "bulk_density",
        "mass_balance_error",

        # Gas concentrations at front surface
        "gas1_conc_front",
        "gas2_conc_front",
        "gas3_conc_front",
        "gas4_conc_front",
        "gas5_conc_front",
        "water_vapor_conc_front",

        # Gas concentrations (averaged)
        "gas1_conc_avg",
        "gas2_conc_avg",
        "gas3_conc_avg",
        "gas4_conc_avg",
        "gas5_conc_avg",
        "water_vapor_conc_avg",

        # Effective properties at front
        "k_eff_front",
        "cp_eff_front",
        "rho_cp_front",
        "alpha_thermal",
        "eps_eff_front",

        # Heat flux components (all in W/m²)
        "q_absorbed",
        "q_reradiation",
        "q_convection",
        "q_conduction_surface",
        "q_net_surface",

        # Energy balance (W/m² and J/m²)
        "energy_stored",
        "energy_rate",
        "total_rxn_heat_domain",

        # Reaction and phase-change rates at front
        "rxn1_rate_front",
        "rxn2_rate_front",
        "rxn3_rate_front",
        "rxn4_rate_front",
        "rxn5_rate_front",
        "phase_change_rate_front",
        "total_rxn_heat_front",

        # Component-specific gas generation rates (post-processed; water net
        # of condensation, from the solver's phase-change extras)
        "mlr_gen_gas1",
        "mlr_gen_gas2",
        "mlr_gen_gas3",
        "mlr_gen_gas4",
        "mlr_gen_gas5",
        "mlr_gen_water",
    ]
    println(io, join(header, ","))

    # Data rows
    for i = 1:length(times)
        data = [
            # Time and basic state
            times[i],
            total_mass[i],
            mlr[i],
            thickness_values[i],
            T_surface[i],
            T_back[i],

            # Temperature at specific depths
            T_25pct[i],
            T_50pct[i],
            T_75pct[i],

            # Solid components
            avg_DF[i],
            avg_int1[i],
            avg_int2[i],
            avg_int3[i],
            avg_int4[i],
            avg_char[i],
            avg_moisture[i],

            # Mass balance
            total_solid_conc[i],
            total_liquid_conc[i],
            total_gas_conc[i],
            sum_all_components[i],
            bulk_density[i],
            mass_balance_error[i],

            # Gas concentrations at front
            gas1_conc_front[i],
            gas2_conc_front[i],
            gas3_conc_front[i],
            gas4_conc_front[i],
            gas5_conc_front[i],
            water_vapor_conc_front[i],

            # Gas concentrations (averaged)
            gas1_conc_avg[i],
            gas2_conc_avg[i],
            gas3_conc_avg[i],
            gas4_conc_avg[i],
            gas5_conc_avg[i],
            water_vapor_conc_avg[i],

            # Effective properties (k_eff, ε_eff from solver; cp, rho_cp, alpha from post-process)
            k_eff_front[i],
            cp_eff_front[i],
            rho_cp_front[i],
            alpha_thermal[i],
            ε_eff_front[i],

            # Heat fluxes (from solver)
            q_absorbed[i],
            q_reradiation[i],
            q_convection[i],
            q_conduction_surface[i],
            q_net_surface[i],

            # Energy balance
            energy_stored[i],
            energy_rate[i],
            total_rxn_heat_domain[i],

            # Reaction and phase-change rates
            rxn1_rate_front[i],
            rxn2_rate_front[i],
            rxn3_rate_front[i],
            rxn4_rate_front[i],
            rxn5_rate_front[i],
            phase_change_rate_front[i],
            total_rxn_heat_front[i],

            # Component gas generation rates
            mlr_gen_gas1[i],
            mlr_gen_gas2[i],
            mlr_gen_gas3[i],
            mlr_gen_gas4[i],
            mlr_gen_gas5[i],
            mlr_gen_water[i],
        ]
        println(io, join(data, ","))
    end
end

println("Comprehensive CSV output saved to: $(results_csv)")
println(
    "  - Total columns: 60 (6 basic + 3 temps + 7 solids + 6 mass + 12 gas + 5 props + 5 heat + 3 energy + 7 reaction and phase-change rates + 6 gen)",
)

# ============================================================================
# Spatial Data Output (for comparison with ThermaKin)
# ============================================================================

println("\nWriting spatial data CSV...")

# Select times for spatial output (similar to ThermaKin output intervals)
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

open(spatial_csv, "w") do io
    # Header
    header = [
        "time_s",
        "from_back_m",
        "from_front_m",
        "temperature_K",
        "conc_DF_kg_m3",
        "conc_DF_int1_kg_m3",
        "conc_DF_int2_kg_m3",
        "conc_DF_int3_kg_m3",
        "conc_DF_int4_kg_m3",
        "conc_Char_kg_m3",
        "conc_Moisture_kg_m3",
        "conc_gas1_kg_m3",
        "conc_gas2_kg_m3",
        "conc_gas3_kg_m3",
        "conc_gas4_kg_m3",
        "conc_gas5_kg_m3",
        "conc_water_vapor_kg_m3",
        "conc_air_kg_m3",
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

        nodes = get_profile(solution, :z, idx)
        centers = solution.cell_centers[:, idx]
        live_front = maximum(filter(isfinite, nodes))

        for cell_idx = 1:n_cells_i
            isfinite(T_profile[cell_idx]) || continue
            from_back = centers[cell_idx] - nodes[1]
            from_front = live_front - centers[cell_idx]

            T_cell = T_profile[cell_idx]
            ξ_cell = [ξ_profile[k, cell_idx] for k = 1:14]
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
                ξ_cell[13],
                ξ_cell[14],
                total_density,
            ]
            println(io, join(data, ","))
        end
    end
end

println("Spatial data saved to: $(spatial_csv)")

# ============================================================================
# Boundary Temperature Information
# ============================================================================
# 
# In cell-centered finite volume, the temperature at the surface is not directly
# computed - we only know cell-center temperatures. For boundary conditions with
# strong temperature dependence (especially radiation with T⁴ scaling), using
# the cell-center temperature causes significant error.
#
# Pyrolysis.jl solves for the true surface temperature T_s that satisfies
# energy balance:
#     q_applied(T_s) = q_conducted(T_s, T_cell)
#
# The solution holds it at every saved time (`solution.surface_T`,
# `solution.boundary_T`), together with the face heat fluxes it balances
# (`solution.extras.surface_heat_flux_*`):
# - Correct boundary flux computation
# - Accurate surface temperature reporting
# - Energy balance verification

println("\n" * "="^70)
println("Boundary Temperature Analysis")
println("="^70)

# Final boundary temperatures (solved face temperatures)
T_surface_final = solution.surface_T[end]
T_back_final = solution.boundary_T[:bottom][end]

# Cell-center temperatures of the outermost and innermost live cells
T_cells_final = solution.T[:, end]
live_final = findall(isfinite, T_cells_final)
T_cell_top = T_cells_final[last(live_final)]
T_cell_bottom = T_cells_final[first(live_final)]

# Compute temperature jumps at boundaries
ΔT_surface = T_surface_final - T_cell_top
ΔT_back = T_back_final - T_cell_bottom

println("\nSurface (top) boundary:")
println("  T_surface (solved)  = $(round(T_surface_final, digits=2)) K")
println("  T_cell (top cell)   = $(round(T_cell_top, digits=2)) K")
println("  ΔT = T_surface - T_cell = $(round(ΔT_surface, digits=2)) K")
println("  (Positive ΔT means surface is hotter than cell center - expected with heating)")

println("\nBack (bottom) boundary:")
println("  T_back (solved)     = $(round(T_back_final, digits=2)) K")
println("  T_cell (bottom)     = $(round(T_cell_bottom, digits=2)) K")
println("  ΔT = T_back - T_cell = $(round(ΔT_back, digits=2)) K")

# Check energy balance at the surface: the applied flux against the flux
# conducted into the surface cell
q_applied_final = solution.extras.surface_heat_flux_net[end]
q_conducted_final = solution.extras.surface_heat_flux_conduction[end]
imbalance = abs(q_applied_final - q_conducted_final) / max(abs(q_applied_final), 1e-12)
println("\nEnergy balance at the top surface (should be < 1e-6):")
status = imbalance < 1e-6 ? "✓" : "✗"
println("  Top surface: relative imbalance = $(round(imbalance, sigdigits=3)) $status")

# Show boundary state details
println("\nBoundary state details:")
println("  Top surface:")
println("    q_applied   = $(round(q_applied_final/1000, digits=2)) kW/m²")
println("    q_conducted = $(round(q_conducted_final/1000, digits=2)) kW/m²")

println("\nSimulation complete!")
println("\nFeatures used in this example:")
println("  ✓ WEIGHTED conductivity mixing (β=0.5; volume weights differ from ThermaKin)")
println("  ✓ WEIGHTED mixing rule for gas transfer coefficient (β_λ=0.5)")
println(
    "  ✓ Gas transfer coefficients (virgin wood: $(λ_wood) m²/s, the bound water's transfer across the grain; char: $(λ_char) m²/s)",
)
println("  ✓ lumped wall and Darcy gas transport (diffusion + pressure-driven flow)")
println(
    "  ✓ Temperature-dependent permeability (virgin: $(κ_virgin) m², char: $(κ_char) m²)",
)
println("  ✓ Species molar masses in the ideal-gas pore pressure (Darcy)")
println(
    "  ✓ Two-way phase change of the moisture (CAST spruce isotherm), its vapour starting in equilibrium, humid ambient",
)
println("  ✓ ALE mesh motion")
println("  ✓ Solver-computed boundary fluxes (MLR, heat flux decomposition)")
