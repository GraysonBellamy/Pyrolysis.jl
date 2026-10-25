"""
Physical constants used throughout the Pyrolysis.jl package.
"""

# Fundamental constants
const R_GAS = 8.31446261815324           # J/(mol·K) - universal gas constant
const STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴) - Stefan-Boltzmann constant

# Reference conditions
const T_REF = 298.15                     # K - standard temperature (25°C)
const P_REF = 101325.0                   # Pa - standard pressure (1 atm)

# Temperature bracket on which component property laws are defined. Every
# temperature-dependent property function clamps its argument to this
# interval before evaluating the user's law, the component constructors
# check positivity of ρ, c, k (and 0 ≤ ε ≤ 1) on it, and the surface-energy
# Newton evaluates its condition at a smooth clamp of its iterate into it (the
# iterate itself is never clamped). Together these make every property law
# total on ℝ: a diverging Newton trial state with a nonphysical temperature
# is a failed Newton step, not a `DomainError` that kills the solve, and a
# law that would turn a cell into a perfect insulator (k ≤ 0) is rejected
# when the component is built rather than discovered as a stalled burnout.
const T_PROPERTY_MIN = 200.0             # K
const T_PROPERTY_MAX = 3000.0            # K

# Floor on the mixture volumetric heat capacity ρcₚ [J/(m³·K)]. ρcₚ sums every
# component, gas included, so it reaches zero only in a cell emptied of all
# mass (a burned-out cell whose pore gas has also diffused out); there 1/ρcₚ
# in dT/dt would explode (a stiffness blow-up, not literal Inf). The floor is
# the thermal mass of pore gas at ambient (ρ ≈ 1.2 kg/m³ × cₚ ≈ 1000 J/(kg·K);
# see the heat-capacity regularization in the Technical Reference). It binds
# only when the cell holds
# essentially nothing: virgin solids have ρcₚ ~ 1e5–1e6 J/(m³·K).
const RHO_CP_FLOOR = 1.0e3               # J/(m³·K) - mixture ρcₚ floor
