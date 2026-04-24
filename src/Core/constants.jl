"""
Physical constants used throughout the Pyrolysis.jl package.

Reference: Design Document Appendix A
"""

# Fundamental constants
const R_GAS = 8.31446261815324           # J/(mol·K) - universal gas constant
const STEFAN_BOLTZMANN = 5.670374419e-8  # W/(m²·K⁴) - Stefan-Boltzmann constant

# Reference conditions
const T_REF = 298.15                     # K - standard temperature (25°C)
const P_REF = 101325.0                   # Pa - standard pressure (1 atm)

# Gas transport reference properties
const M_REF_DIFFUSION = 0.029            # kg/mol - reference molar mass for species diffusion scaling
                                         # Default is air/N₂ (≈ 0.029 kg/mol)
                                         # Diffusion coefficient scales as D ∝ 1/√M from kinetic theory
                                         # Species scaling: Dⱼ = D_ref × √(M_ref / Mⱼ)

# Typical material property ranges (for validation)
const RHO_POLYMER_RANGE = (800.0, 2000.0)      # kg/m³
const CP_POLYMER_RANGE = (1000.0, 3000.0)      # J/(kg·K)
const K_POLYMER_RANGE = (0.1, 0.5)             # W/(m·K)
const EMISSIVITY_RANGE = (0.7, 1.0)            # -

# Typical kinetic parameter ranges
const A_RANGE = (1e8, 1e20)                    # s⁻¹
const E_RANGE = (1e5, 3e5)                     # J/mol
