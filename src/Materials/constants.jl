"""
Physical constants used throughout the Pyrolysis.jl package.
"""

# Fundamental constants
const R_GAS = 8.31446261815324           # J/(mol·K) - universal gas constant
const STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴) - Stefan-Boltzmann constant

# Reference conditions
const T_REF = 298.15                     # K - standard temperature (25°C)
const P_REF = 101325.0                   # Pa - standard pressure (1 atm)

# Floor on the matrix volumetric heat capacity ρcₚ_eff [J/(m³·K)]. Formulation A
# excludes gas storage, so a fully-converted (100%-gas-yield) Eulerian cell
# drives the matrix-only ρcₚ → 0 and 1/ρcₚ in dT/dt explodes (a stiffness
# blow-up, not literal Inf). The floor is the thermal mass actually present in
# such a cell — the pore gas (ρ ≈ 1.2 kg/m³ × cₚ ≈ 1000 J/(kg·K) ≈ 1.2e3 at
# ambient; Gpyro floors the same quantity). It binds only when the condensed
# phase is essentially gone: virgin solids have ρcₚ ~ 1e5–1e6 J/(m³·K).
const RHO_CP_FLOOR = 1.0e3               # J/(m³·K) - matrix ρcₚ floor
