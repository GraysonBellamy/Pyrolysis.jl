"""
	BoundaryConditions

Concrete boundary-condition types for the Pyrolysis.jl package.

Holds the thermal / mass / pressure boundary-condition types, the type-stable
`BoundaryConditionSet` container, and the per-species ambient-concentration
helper.

It also holds the face-local, mesh-free evaluation of those conditions: the
surface energy balance of a boundary face (`solve_boundary_temperature!`) with
the flux and derivative of every thermal condition, the blowing models' effect
on a convective face (`BlowingFace`, `blowing_state`), the mass-condition film
fluxes (`mass_bc_flux`), the heat-flux decomposition for diagnostics, and the
`BoundaryState` that records a converged balance.

Depends on `Materials` for the abstract BC supertypes (`AbstractThermalBC`,
`AbstractMassBC`), the property-function
hierarchy, the component accessors, and physical constants, and on
`Properties` for the gas diffusivity of the blowing models.
"""
module BoundaryConditions

using ..Materials

include("boundary_conditions.jl")
include("boundary_state.jl")
include("blowing.jl")
include("boundary_fluxes.jl")
include("surface_balance.jl")

# Thermal boundary conditions
export AdiabaticBC, HeatFluxBC, ConvectiveBC, RadiativeBC, DirichletThermalBC
export CombinedThermalBC, RadiativeFluxBC

# Mass boundary conditions
export ImpermeableBC, MassFluxBC, ConvectiveMassBC, DiffusiveMassBC, CombinedMassBC

# Pressure boundary conditions (Darcy flow)
export DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC

# Radiative-flux helpers
export get_radiative_flux_bc

# Type-stable boundary-condition container
export BoundaryConditionSet, get_thermal_bc, get_mass_bc

# Ambient-temperature evaluation (evaluate_param / evaluate_param_T live in Materials)
export ambient_concentration

# Surface energy balance of a boundary face and its diagnostics
export BoundaryState
export solve_boundary_temperature!
export evaluate_bc_flux_and_derivative
export compute_boundary_heat_fluxes, compute_surface_mass_flux

# Blowing models evaluated at a face
export blowing_factor, BlowingFace, blowing_state, blowing_ratio

# Mass-condition film fluxes
export mass_bc_flux

end # module BoundaryConditions
