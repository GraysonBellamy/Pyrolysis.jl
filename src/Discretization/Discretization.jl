"""
	Discretization

Spatial discretization (finite volume) for the Pyrolysis.jl package.

This module provides:
- Finite volume operators (divergence, RHS assembly)
- Flux limiters (minmod, vanleer, superbee)
- ALE flux computation for moving meshes
- Boundary flux implementation
- **Boundary temperature solve** for correct surface temperature in cell-centered FVM

## Boundary Temperature Enhancement

In cell-centered finite volume methods, temperature is stored at cell centers, not boundaries.
For temperature-dependent boundary conditions (especially radiation with T⁴ scaling),
using the cell-center temperature directly causes significant error (~20% energy imbalance).

The `solve_boundary_temperature` function uses Newton iteration to find the true
surface temperature T_s satisfying energy balance:

	q_applied(T_s) = q_conducted(T_s, T_cell)

where:
- `q_applied = q_ext - εσT_s⁴ - h(T_s - T∞)` is the net heat flux from the boundary condition
- `q_conducted = k(T_s - T_cell)/Δz` is the heat flux into the material

This ensures energy balance at boundaries is satisfied to machine precision.
The solved surface temperature is stored in `mesh.boundary_states` and used for:
- Correct flux computation in the RHS
- Accurate surface temperature output
- Energy balance diagnostics

See `docs/boundary_temperature_enhancement.md` for detailed physics and implementation.
"""
module Discretization

using ..Materials
using ..BoundaryConditions
using ..Materials: GAS
using ..BoundaryConditions: DirichletThermalBC
using ..Geometry
using ..Properties
using ..Physics

# Include subfiles
# Note: `darcy_velocity.jl` must be included before `finite_volume.jl`
# because `finite_volume.jl` calls `apply_pressure_bc_to_darcy_velocity!`.
include("blowing.jl")  # blowing models: BlowingFace, blowing_factor, blowing_state
include("boundary_fluxes.jl")
include("boundary_temperature.jl")  # Boundary temperature solve for correct surface T
include("darcy_velocity.jl")
include("finite_volume.jl")
include("flux_limiters.jl")
include("ale_fluxes.jl")

# Exports - Finite Volume
export divergence_1d!, divergence_species_1d!
export compute_rhs!
export RHSCache
export update_darcy_properties!, compute_rhs_darcy_fick!
# Type-stable BC application functions
export apply_thermal_bc_flux!, apply_mass_bc_flux!

# Exports - Flux Limiters
export Limiter, MINMOD, VANLEER, SUPERBEE, NONE, CENTRAL
export minmod, vanleer, superbee
export apply_limiter

# Exports - ALE mesh-relative temperature advection
export compute_ale_advection_term!

# Exports - Darcy Velocity
export apply_pressure_bc_to_darcy_velocity!  # Pressure BC application
export boundary_darcy_velocity  # Shared boundary-face Darcy velocity kernel

# Exports - Boundary Fluxes
export thermal_bc_flux, mass_bc_flux
export compute_boundary_heat_fluxes, compute_surface_mass_flux  # Diagnostic decomposition

# Exports - Boundary Temperature Solve
export solve_boundary_temperature, update_boundary_state!
export evaluate_bc_flux_and_derivative
export blowing_factor, BlowingFace, blowing_state, blowing_ratio
export check_boundary_energy_balance
export get_surface_temperature, get_back_face_temperature

# P1 radiation flux assembly (solve_quasi_steady_radiation!,
# initialize_g_from_temperature!, update_radiation_properties!,
# compute_radiation_flux_divergence!) lives in Pyrolysis.Experimental.

end # module Discretization
