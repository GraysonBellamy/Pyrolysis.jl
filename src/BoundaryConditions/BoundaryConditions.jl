"""
	BoundaryConditions

Concrete boundary-condition types for the Pyrolysis.jl package.

Promoted out of the former `Core` module. Holds the thermal / mass / pressure /
mesh-motion boundary-condition types, the type-stable `BoundaryConditionSet`
container, and the ambient-temperature evaluation helper. The concrete
radiation-BC types (Marshak / Reflective / Dirichlet) are P1-only and live in
`Pyrolysis.Experimental`.

Depends on `Materials` for the abstract BC supertypes (`AbstractThermalBC`,
`AbstractMassBC`, `AbstractRadiationBC`), the property-function
hierarchy, the component accessors, and physical constants.
"""
module BoundaryConditions

using StaticArrays
using ..Materials

include("boundary_conditions.jl")

# Thermal boundary conditions
export AdiabaticBC, HeatFluxBC, ConvectiveBC, RadiativeBC, DirichletThermalBC
export CombinedThermalBC, RadiativeFluxBC

# Mass boundary conditions
export ImpermeableBC, MassFluxBC, ConvectiveMassBC, DiffusiveMassBC, CombinedMassBC

# Pressure boundary conditions (Darcy flow)
export DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC

# Mesh-motion boundary conditions

# Radiative-flux helpers
export get_radiative_flux_bc, extract_radiative_flux_from_bcs, has_radiative_flux_bc

# Type-stable boundary-condition container
export BoundaryConditionSet, get_thermal_bc, get_mass_bc

# Ambient-temperature evaluation (evaluate_param / evaluate_param_T live in Materials)
export ambient_temperature, ambient_concentration

# Concrete radiation boundary conditions (Marshak / Reflective / Dirichlet,
# RadiationBCSet, get_radiation_bc, default_radiation_bcs) live in
# Pyrolysis.Experimental — they are P1-only.

end # module BoundaryConditions
