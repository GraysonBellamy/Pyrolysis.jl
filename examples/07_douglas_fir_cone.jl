"""
Cone Calorimeter Simulation of Douglas Fir with ALE Mesh Motion and Surface Absorption

This example simulates the pyrolysis of Douglas Fir wood under standard cone 
calorimeter conditions with ALE (Arbitrary Lagrangian-Eulerian) mesh motion 
to track the material as it decomposes:
- 1D simulation through sample thickness
- ALE mesh motion for accurate shrinkage modeling
- Surface absorption radiation model for efficient opaque material simulation
- Multi-step pyrolysis kinetics with char formation
- Structured Jacobian (`Structured(strategy=DirectSolve())`) for fast assembly
- Bruggeman mixing rules for effective properties
- Species-specific transport (light gases diffuse faster)

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

Douglas Fir Properties (from DF.cmp):
- Virgin wood density: 530 kg/m³
- Multi-step decomposition with 5 intermediate states
- Final char density: 150 kg/m³
- Includes moisture content with evaporation
- Products: char + various gases + water vapor

New Features Demonstrated:
- Surface absorption radiation model (radiation_model=SURFACE_ABSORPTION)
- Physically-based gas transfer coefficients (λ ~ 1e-5 to 1e-4 m²/s)
- Darcy-Fick gas transport (diffusion + pressure-driven flow)
- Temperature-dependent permeability evolution (virgin wood → porous char)
- Species-specific molar masses for transport scaling
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
# Material Definition: Douglas Fir
# ============================================================================

"""
Define Douglas Fir material with temperature-dependent properties.

Douglas Fir undergoes multi-step thermal decomposition:
	DF → DF_int1 → DF_int2 → DF_int3 → DF_int4 → Char + gases

The material also contains moisture that evaporates:
	MOISTURE → WATER_VAPOR

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

# Moisture-dependent heat of sorption
# HoS(MC) = 2450.4 + 875.6 * exp(-0.1122 * MC) [kJ/kg]
# where MC is moisture content (%) on dry basis: MC = 100 * ξ_moisture / ρ_dry
# As moisture decreases, more energy is required to release the remaining bound water
const MOISTURE_IDX = 7           # Component index for moisture
const ρ_DRY_BASIS = 530.0        # Dry wood density for MC calculation [kg/m³]

heat_of_sorption = StateDependentProperty(
    (T, ξ) -> begin
        # Moisture content on dry basis (%)
        MC = 100.0 * ξ[MOISTURE_IDX] / ρ_DRY_BASIS
        # # Heat of sorption [J/kg] - increases as moisture decreases
        return (2450.420909987441 + 875.6441221735336 * exp(-0.1121684273316275 * MC)) * 1e3
    end,
)

# ============================================================================
# Gas Transfer Coefficients and Permeability (Physically Realistic)
# ============================================================================
# 
# Gas transport uses the volume-fraction-gradient (partial pressure) formulation:
#   J = -(λ/T) × ∂(ξT)/∂x = -λ × [∂ξ/∂x + (ξ/T)∂T/∂x]
#
# Each component has a λ [m²/s] value representing its intrinsic gas diffusivity.
# The effective mixture λ is computed via WEIGHTED mixing (β=0.5, parallel/series).
#
# For SOLID components, λ represents gas diffusion through that phase's pore structure:
#   λ_solid ≈ D_free × (ε/τ)
#   where D_free is free molecular diffusion (~3e-5 m²/s at 500K for typical gas),
#   ε is the porosity, and τ is the tortuosity.
#
# For GAS components, λ represents free molecular diffusion of that species:
#   λ_gas ≈ D_free(T, M) from Chapman-Enskog theory
#   D_free scales as T^1.5 / √M. At 500K, air-like gas: ~3e-5 m²/s
#   Lighter molecules (H₂O) diffuse ~1.5-2× faster than heavier ones (tars).
#
# For LIQUID (moisture): λ represents moisture migration rate through the wood
#   structure (capillary + vapor pressure driven). Much faster than gas diffusion.

# Solid phase gas transfer coefficients (D_free × ε/τ estimates)
# D_free at ~500-600K for typical pyrolysis gas ≈ 3-5e-5 m²/s
λ_wood = 5.0e-7     # m²/s - virgin wood (ε~0.05, τ~5 → ε/τ~0.01)
λ_char = 1.0e-5     # m²/s - porous char (ε~0.7, τ~2 → ε/τ~0.35)

# Gas phase free molecular diffusion (at representative ~500K)
# Reference: D(M=29 gas in air, 500K) ≈ 3e-5 m²/s (Chapman-Enskog)
λ_gas = 3.0e-5      # m²/s - reference for M~30 gas at ~500K

# Permeability for Darcy flow (pressure-driven bulk transport through solid matrix)
# κ determines resistance to pressure-driven gas flow, separate from Fickian diffusion.
# Virgin softwood (transverse): κ ~ 1e-16 to 1e-14 m²
# Char (open pore network): κ ~ 1e-12 to 1e-11 m²
κ_virgin = 1.0e-15  # m² - tight virgin wood cellular structure
κ_char = 1.0e-12    # m² - open char pore network

# Create Douglas Fir material with multiple components:
# Solid components: DF, DF_int1, DF_int2, DF_int3, DF_int4, Char, MOISTURE
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
        # 6. Char - final solid residue (fully porous)
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
        # 7. Moisture (liquid water bound in wood structure)
        # ThermaKin uses ρ=1e6 kg/m³ as a "pseudo-density" to give moisture
        # negligible volume fraction. This matches L SWELLING: 0 behavior where
        # bound liquid doesn't contribute to the material's bulk volume.
        # Physical interpretation: moisture is absorbed into wood cell walls,
        # not occupying separate volume.
        LiquidComponent(
            :MOISTURE,
            ρ = 1e6,                       # kg/m³ - pseudo-density (matches ThermaKin)
            c = moisture_cp,               # Moisture contributes to heat capacity
            k = 0.145,                       # W/(m·K) - but not significantly to k_eff
            ε = 0.81,
            α = 0.0,                       # m²/kg - not used with SURFACE_ABSORPTION model
            κ = 0.0,                        # m² - liquid blocks Darcy flow
        ),
        # 8. DF_gas1 - gaseous product from reaction 1 (5% yield, low-T onset)
        # Early hemicellulose decomposition: acetic acid (60), furfural (96), formic acid (46)
        # Gas component λ values: free molecular diffusion scaled by √(M_ref/M)
        # D_free ∝ 1/√M (kinetic theory). Reference: λ_gas at M_ref=0.029 (air).
        # At 500K: D(M=18)≈5.4e-5, D(M=60)≈2.6e-5, D(M=110)≈1.7e-5 m²/s
        GaseousComponent(
            :DF_gas1,
            M = 0.060,                     # kg/mol - acetic acid-like (CH₃COOH = 60 g/mol)
            c = 1400.0,                    # J/(kg·K) - organic acids/aldehydes
            k = 0.020,                     # W/(m·K) - organic vapors
            λ = λ_gas * sqrt(0.029/0.060), # m²/s - ~2.1e-5 (moderate molecule)
        ),
        # 9. DF_gas2 - gaseous product from reaction 2 (20% yield)
        # Continued hemicellulose / early cellulose: mix of organic fragments, some CO₂
        GaseousComponent(
            :DF_gas2,
            M = 0.070,                     # kg/mol - mix of organics heavier than pure CO₂
            c = 1400.0,                    # J/(kg·K) - organic fragments
            k = 0.020,                     # W/(m·K)
            λ = λ_gas * sqrt(0.029/0.070), # m²/s - ~1.9e-5 (heavier organics)
        ),
        # 10. DF_gas3 - gaseous product from reaction 3 (63% yield — dominant gas-producing step)
        # Main cellulose pyrolysis: levoglucosan (162), anhydrosugars, heavy tars
        GaseousComponent(
            :DF_gas3,
            M = 0.110,                     # kg/mol - dominated by levoglucosan/anhydrosugars
            c = 1300.0,                    # J/(kg·K) - heavy organics, lower mass-specific cp
            k = 0.018,                     # W/(m·K) - heavy organic vapors
            λ = λ_gas * sqrt(0.029/0.110), # m²/s - ~1.5e-5 (heavy tars, slowest)
        ),
        # 11. DF_gas4 - gaseous product from reaction 4 (21% yield)
        # Secondary decomposition: lighter organics, phenol fragments from lignin
        GaseousComponent(
            :DF_gas4,
            M = 0.070,                     # kg/mol - phenolics (94), guaiacol (124), lighter organics
            c = 1500.0,                    # J/(kg·K) - medium organics
            k = 0.020,                     # W/(m·K)
            λ = λ_gas * sqrt(0.029/0.070), # m²/s - ~1.9e-5
        ),
        # 12. DF_gas5 - gaseous product from reaction 5 (26% yield, final charring)
        # Last volatiles released during charring: CH₄, light hydrocarbons, residual phenolics
        GaseousComponent(
            :DF_gas5,
            M = 0.050,                     # kg/mol - mix of CH₄ (16), phenol fragments, light organics
            c = 1600.0,                    # J/(kg·K) - lighter molecules, higher mass-specific cp
            k = 0.025,                     # W/(m·K)
            λ = λ_gas * sqrt(0.029/0.050), # m²/s - ~2.3e-5 (lighter mix)
        ),
        # 13. Water vapor - gaseous product from moisture evaporation (lightest, fastest)
        GaseousComponent(
            :WATER_VAPOR,
            M = 0.018,                     # kg/mol - H₂O (18 g/mol, lightest species)
            c = water_vapor_cp,            # J/(kg·K) - NIST polynomial fit
            k = 0.025,                     # W/(m·K) - H₂O vapor at ~400K
            λ = λ_gas * sqrt(0.029/0.018), # m²/s - ~3.8e-5 (lightest molecule, fastest)
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
            λ = λ_gas,                     # m²/s - reference M~29 diffusivity
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
        # Reaction 6: MOISTURE → WATER_VAPOR
        # Heat of sorption: HoS(MC) = 2450.4 + 875.6 * exp(-0.1122 * MC) [kJ/kg]
        # Moisture-dependent: increases as moisture content decreases (bound water)
        Reaction(
            :moisture_evaporation,
            7 => 13,
            A = 2.997e5,                   # s⁻¹
            E = 5.0e4,                     # J/mol
            h = heat_of_sorption,          # State-dependent heat [J/kg]
            n = 1.0,
        ),
    ),
    # Mixing rules for effective properties (matches ThermaKin2Ds default behavior)
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
n_cells = 240             # 240 cells → 0.05 mm resolution

# Time parameters  
t_final = 600.0             # 10 minutes - typical cone test duration for wood
dt_initial = 1e-3           # initial timestep

# Environmental conditions
T_ambient = 300.0           # K - ambient temperature (27°C)
q_external = 75e3  # W/m² - 50 kW/m² external heat flux

# Convective heat transfer coefficient for cone calorimeter
h_conv = 10.0                # W/(m²·K)

# ============================================================================
# Mesh Creation
# ============================================================================

# Number of components (NC = 14):
# 6 solids (DF, DF_int1-4, Char) + 1 liquid (MOISTURE) + 7 gases (DF_gas1-5, WATER_VAPOR, AIR)

# Initial concentrations:
# - DF at full density (530 kg/m³)
# - Moisture content 15% on dry basis: 530 * 0.15 ≈ 79.5 kg/m³
# - Background pore air at ~1 atm partial pressure
# - All others initially zero
moisture_content_db = 0.06  # 6% moisture content on dry basis
ρ_DF_initial = 530.0        # Initial wood density (dry)
# For dry basis, moisture density = MC_db * ρ_dry
ρ_moisture_initial = moisture_content_db * ρ_DF_initial

# Pore air concentration giving P = 101325 Pa in the initial pore space
# (porosity follows from the condensed composition and intrinsic densities)
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
ξ_air_initial =
    background_gas_concentration(douglas_fir, :AIR, ξ_condensed_initial, T_ambient)

mesh = create_uniform_mesh(
    sample_thickness,      # domain length [m]
    n_cells,               # number of cells
    Val(14);               # number of components (NC)
    T_initial = T_ambient, # initial temperature
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
        ξ_air_initial,
    ),  # DF + moisture + pore air, rest zero
)

println("Mesh created:")
println("  - Thickness: $(sample_thickness*1000) mm")
println("  - Number of cells: $n_cells")
println("  - Cell size: $(sample_thickness/n_cells*1000) mm")
println("  - Initial wood density: $ρ_DF_initial kg/m³")
println(
    "  - Initial moisture: $(round(ρ_moisture_initial, digits=1)) kg/m³ ($(round(moisture_content_db*100, digits=1))% dry basis)",
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

# Top surface (boundary_id = 1): exposed surface
# - External radiant heat flux with Beer-Lambert in-depth absorption
# - Re-radiation losses (surface radiates to surroundings)
# - Convective losses (natural + forced convection)
# Note: RadiativeFluxBC and RadiativeBC now use the effective emissivity from material
# properties by default (ε = nothing). This automatically tracks the surface emissivity
# as the material decomposes (e.g., virgin wood → char transition). To override with
# a fixed value, pass ε = 0.9 or absorptivity = 0.9 explicitly.
# Note: radiation_model=SURFACE_ABSORPTION controls surface vs in-depth absorption
top_thermal_bc =
    RadiativeFluxBC(flux = q_external) +                      # Cone heater (uses material ε)
    RadiativeBC(T_inf = 300) +                                # Re-radiation losses (uses material ε)
    ConvectiveBC(h = h_conv, T_inf = 300)                     # Convective losses

# Top surface mass BC: diffusive gas escape with Darcy flow to atmosphere
# For cone calorimeter with continuous exhaust in nitrogen:
# - DiffusiveMassBC(h_m): External mass transfer coefficient [m/s] (matches ThermaKin "a" parameter)
# - DirichletPressureBC: Surface at atmospheric pressure (exhaust maintains ~1 atm)
# This combines Fickian diffusion (driven by λ) with Darcy flow (driven by pressure)
# h_m = 0.05 m/s matches ThermaKin boundary condition: "DF_gas1 LIN 0.05 0"
top_mass_bc = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0)

# Bottom surface (boundary_id = 2): back face
# bottom_thermal_bc = RadiativeBC(T_inf = 300) +                                # Re-radiation losses (uses material ε)
# 					ConvectiveBC(h = h_conv, T_inf = 300)                     # Convective losses
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
    z -> 0.0,              # 13. WATER_VAPOR
    z -> ξ_air_initial,    # 14. AIR (background pore air at ~1 atm)
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
    output_times = collect(range(0.0, t_final, length = 3601)),  # Output every 1s
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
# The reduced Jacobian omits the dense ALE mesh velocity coupling, reducing the
# Jacobian from 588K to 141K non-zeros and restoring O(N) factorization scaling.
# Combined with the analytical Jacobian, this gives ~400× speedup over dense FD.
# The solution is identical to the full Jacobian within solver tolerance.
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
#
# The solution now provides:
# - solution.extras.surface_mass_flux: Total mass leaving at surface (computed during solve)
# - solution.extras.gas_generation_rate: Total volumetric gas production (computed during solve)
#
# Component-specific rates still require post-processing to separate by species.
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
mlr_gen_water = zeros(n_times)     # From reaction 6

# Pre-allocate arrays for reaction rate computation
NC = 14  # Number of components
NR = 6   # Number of reactions
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
        Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ, T, douglas_fir)

        # Accumulate gas production from each reaction
        # Reaction 1: DF → 0.95 DF_int1 + 0.05 DF_gas1
        mlr_gen_gas1[i] += douglas_fir.reactions[1].stoich[4] * rates[1] * cell_vol

        # Reaction 2: DF_int1 → 0.80 DF_int2 + 0.20 DF_gas2
        mlr_gen_gas2[i] += douglas_fir.reactions[2].stoich[4] * rates[2] * cell_vol

        # Reaction 3: DF_int2 → 0.37 DF_int3 + 0.63 DF_gas3
        mlr_gen_gas3[i] += douglas_fir.reactions[3].stoich[4] * rates[3] * cell_vol

        # Reaction 4: DF_int3 → 0.79 DF_int4 + 0.21 DF_gas4
        mlr_gen_gas4[i] += douglas_fir.reactions[4].stoich[4] * rates[4] * cell_vol

        # Reaction 5: DF_int4 → 0.74 Char + 0.26 DF_gas5
        mlr_gen_gas5[i] += douglas_fir.reactions[5].stoich[4] * rates[5] * cell_vol

        # Reaction 6: MOISTURE → WATER_VAPOR (single product reaction)
        mlr_gen_water[i] += 1.0 * rates[6] * cell_vol  # yield = 1.0 for single product
    end
end

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
println("  - Peak H₂O generation: $(round(maximum(mlr_gen_water)*1000, digits=2)) g/(m²·s)")

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
    label = "H₂O",
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
plot!(plot_title = "Cone Calorimeter Simulation: Douglas Fir at 50 kW/m² (ALE Mesh)")

# Save plot
savefig(final_plot, "cone_calorimeter_douglas_fir_results.png")
println("Results saved to: cone_calorimeter_douglas_fir_results.png")

# ============================================================================
# Component Evolution Plot
# ============================================================================

println("\nGenerating component evolution plot...")

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
        avg_DF[i] = sum(ξ_profile[1, :]) / n_cells_i
        avg_int1[i] = sum(ξ_profile[2, :]) / n_cells_i
        avg_int2[i] = sum(ξ_profile[3, :]) / n_cells_i
        avg_int3[i] = sum(ξ_profile[4, :]) / n_cells_i
        avg_int4[i] = sum(ξ_profile[5, :]) / n_cells_i
        avg_char[i] = sum(ξ_profile[6, :]) / n_cells_i
        avg_moisture[i] = sum(ξ_profile[7, :]) / n_cells_i
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
# The solver now provides these quantities computed during the solve:
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
total_rxn_heat_domain = zeros(n_times) # Total reaction heat (integrated over domain) [W/m²]

# Reaction rate diagnostics at front surface
rxn1_rate_front = zeros(n_times)
rxn2_rate_front = zeros(n_times)
rxn3_rate_front = zeros(n_times)
rxn4_rate_front = zeros(n_times)
rxn5_rate_front = zeros(n_times)
rxn6_rate_front = zeros(n_times)
total_rxn_heat_front = zeros(n_times)

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
        gas1_conc_avg[i] = sum(ξ_profile[8, :]) / n_cells_i
        gas2_conc_avg[i] = sum(ξ_profile[9, :]) / n_cells_i
        gas3_conc_avg[i] = sum(ξ_profile[10, :]) / n_cells_i
        gas4_conc_avg[i] = sum(ξ_profile[11, :]) / n_cells_i
        gas5_conc_avg[i] = sum(ξ_profile[12, :]) / n_cells_i
        water_vapor_conc_avg[i] = sum(ξ_profile[13, :]) / n_cells_i

        gases =
            gas1_conc_avg[i] +
            gas2_conc_avg[i] +
            gas3_conc_avg[i] +
            gas4_conc_avg[i] +
            gas5_conc_avg[i] +
            water_vapor_conc_avg[i]

        total_solid_conc[i] = solids
        total_liquid_conc[i] = liquids
        total_gas_conc[i] = gases
        sum_all_components[i] = solids + liquids + gases
        bulk_density[i] = total_mass[i] / thickness_values[i]
        mass_balance_error[i] = (sum_all_components[i] - bulk_density[i]) / bulk_density[i]

        # === Front Surface Diagnostics ===
        # Front surface is the last cell (z = L, boundary_id = 1)
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
        # Depth measured from surface (z=L), so:
        # 25% depth = cell at index ~0.75*n_cells from bottom
        # 50% depth = midplane
        # 75% depth = cell at index ~0.25*n_cells from bottom
        idx_25pct = max(1, round(Int, 0.75 * n_cells_i))
        idx_50pct = max(1, round(Int, 0.50 * n_cells_i))
        idx_75pct = max(1, round(Int, 0.25 * n_cells_i))
        T_25pct[i] = T_profile[idx_25pct]
        T_50pct[i] = T_profile[idx_50pct]
        T_75pct[i] = T_profile[idx_75pct]

        # === Energy Balance ===
        # Total thermal energy: E = ∫ ρ cp T dV (per unit area)
        # E = Σ (ρ_i × cp_i × T_i × Δz_i)
        dz_cell = thickness_values[i] / n_cells_i
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
            energy_stored[i] += ρ_cell * cp_cell * (T_cell - 298.15) * dz_cell

            # Reaction heat at this cell: Q = -Σ(hᵢ × rᵢ) [W/m³]
            # heats[i] is in J/kg, rates[i] is in kg/(m³·s)
            # So heats[i] * rates[i] = W/m³ (volumetric heating rate)
            # Note: In our convention, h < 0 means exothermic (releases heat)
            # So Q = -Σ(h × r) gives positive Q for exothermic reactions
            local_rates = zeros(6)
            local_heats = zeros(6)
            Pyrolysis.Physics.evaluate_reactions!(
                local_rates,
                local_heats,
                ξ_cell,
                T_cell,
                douglas_fir,
            )
            # Sum of (-h × r) × dz gives W/m² for each cell contribution
            for rxn_idx = 1:6
                total_rxn_heat_domain[i] +=
                    -local_heats[rxn_idx] * local_rates[rxn_idx] * dz_cell
            end
        end

        # Energy rate: dE/dt (finite difference)
        if i > 1 && times[i] > times[i-1]
            dt = times[i] - times[i-1]
            energy_rate[i] = (energy_stored[i] - energy_stored[i-1]) / dt
        end

        # === Reaction Rates at Front Surface ===
        global rates = zeros(6)
        global heats = zeros(6)
        Pyrolysis.Physics.evaluate_reactions!(rates, heats, ξ_front, T_front, douglas_fir)

        rxn1_rate_front[i] = rates[1]
        rxn2_rate_front[i] = rates[2]
        rxn3_rate_front[i] = rates[3]
        rxn4_rate_front[i] = rates[4]
        rxn5_rate_front[i] = rates[5]
        rxn6_rate_front[i] = rates[6]

        # Total reaction heat at front (W/m³)
        # Q = -Σ(hᵢ × rᵢ) where h < 0 = exothermic gives Q > 0
        total_rxn_heat_front[i] = 0.0
        for rxn_idx = 1:6
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
open("douglas_fir_results_75kWm2_fit.csv", "w") do io
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

        # Reaction rates at front
        "rxn1_rate_front",
        "rxn2_rate_front",
        "rxn3_rate_front",
        "rxn4_rate_front",
        "rxn5_rate_front",
        "rxn6_rate_front",
        "total_rxn_heat_front",

        # Component-specific gas generation rates (post-processed)
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

            # Reaction rates
            rxn1_rate_front[i],
            rxn2_rate_front[i],
            rxn3_rate_front[i],
            rxn4_rate_front[i],
            rxn5_rate_front[i],
            rxn6_rate_front[i],
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

println("Comprehensive CSV output saved to: douglas_fir_results_75kWm2_fit.csv")
println(
    "  - Total columns: 59 (6 basic + 3 temps + 7 solids + 6 mass + 12 gas + 5 props + 6 heat + 3 energy + 7 rxn + 6 gen)",
)
println("  - Note: mlr and heat flux components now from solver (more accurate)")

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

open("douglas_fir_spatial_75kWm2_fit.csv", "w") do io
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

        # Get current thickness
        t_idx = argmin(abs.(thickness_times .- actual_t))
        current_thickness = thickness_values[t_idx]

        # Cell positions (from back, z=0 at back)
        dz = current_thickness / n_cells_i

        for cell_idx = 1:n_cells_i
            from_back = (cell_idx - 0.5) * dz  # Cell center position
            from_front = current_thickness - from_back

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

println("Spatial data saved to: douglas_fir_spatial_75kWm2_fit.csv")

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
# This is stored in mesh.boundary_states and used for:
# - Correct boundary flux computation
# - Accurate surface temperature reporting
# - Energy balance verification

using Pyrolysis.Discretization:
    get_surface_temperature, get_back_face_temperature, check_boundary_energy_balance

println("\n" * "="^70)
println("Boundary Temperature Analysis")
println("="^70)

# Get final boundary state information. `solve` integrates its own copy of the
# mesh; `solution.problem.mesh` holds the final state.
final_mesh = solution.problem.mesh
T_surface_final = get_surface_temperature(final_mesh)
T_back_final = get_back_face_temperature(final_mesh)

# Get cell-center temperatures for comparison
n_cells_final = length(final_mesh.cell_states)
T_cell_top = final_mesh.cell_states[n_cells_final].T
T_cell_bottom = final_mesh.cell_states[1].T

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

# Check energy balance at boundaries
balance_results = check_boundary_energy_balance(final_mesh)
println("\nEnergy balance at boundaries (should be < 1e-6):")
for (boundary_id, imbalance) in balance_results
    boundary_name = boundary_id == 1 ? "Top surface" : "Back surface"
    status = imbalance < 1e-6 ? "✓" : "✗"
    println(
        "  $(boundary_name): relative imbalance = $(round(imbalance, sigdigits=3)) $status",
    )
end

# Show boundary state details
println("\nBoundary state details:")
top_state = final_mesh.boundary_states[1]
println("  Top surface:")
println("    q_applied   = $(round(top_state.q_applied/1000, digits=2)) kW/m²")
println("    q_conducted = $(round(top_state.q_conducted/1000, digits=2)) kW/m²")
println("    converged   = $(top_state.converged) ($(top_state.iterations) iterations)")

println("\nSimulation complete!")
println("\nFeatures used in this example:")
println("  ✓ WEIGHTED mixing rule for thermal conductivity (β=0.5, matches ThermaKin)")
println("  ✓ WEIGHTED mixing rule for gas transfer coefficient (β_λ=0.5)")
println("  ✓ Physically-based gas transfer coefficients (λ ~ 1e-5 to 1e-4 m²/s)")
println("  ✓ Darcy-Fick gas transport (diffusion + pressure-driven flow)")
println(
    "  ✓ Temperature-dependent permeability (virgin: $(κ_virgin) m², char: $(κ_char) m²)",
)
println("  ✓ Species-specific molar masses for transport scaling")
println("  ✓ ALE mesh motion with analytical Jacobian")
println("  ✓ Solver-computed boundary fluxes (MLR, heat flux decomposition)")
