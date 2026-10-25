"""
	Pyrolysis.jl

A Julia package for simulating thermal decomposition (pyrolysis) of materials
under high heat flux with a finite-volume method on a mesh that moves with the
solid (Arbitrary Lagrangian–Eulerian).

# Main Features
- Multi-component materials with temperature-dependent properties
- Arrhenius-type decomposition kinetics (single and multi-step)
- Gas generation and transport through porous char
- In-depth radiation absorption (Beer-Lambert)
- Material volume changes (shrinkage, swelling, intumescence)
- Surface regression and internal interface motion
- Surface burnout by cell deactivation

# Quick Start
```julia
using Pyrolysis

# Define a simple material
pmma = Material(
	name = :porous_polymer,
	gas_transport=LumpedTransfer(endpoint=LumpedEndpoint(λ=1e-6, κ=0., max_skeletal_occupancy=.001)),
	components = (
		SolidComponent(:virgin, ρ=1190, φ=.3, c=1420, k=0.21, ε=0.86, λ=1e-6),
		GaseousComponent(:MMA, M=0.100, c=1500, k=0.02,)
	),
	reactions = (
		Reaction(:decomposition, :virgin => :MMA, A=8.5e12, E=188e3, h=870e3),  # h > 0 = endothermic
	)
)

# Create mesh and solve
mesh = create_uniform_mesh(0.01, 100)  # 1 cm, 100 cells
bcs = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
	mass = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(
	mesh = mesh,
	material = pmma,
	bc_set = bcs,
	T_initial = z -> 300.0,
	ξ_initial = [z -> 1190.0, z -> 0.0],
	tspan = (0.0, 100.0),
)
# A mesh that moves with the solid; the cells the surface burns away deactivate
solution = solve(problem; use_ale = true, deactivation = true)
```
"""
module Pyrolysis

using LinearAlgebra
using SparseArrays
using StaticArrays
using ForwardDiff
using LinearSolve
using OrdinaryDiffEq


# ============================================================================
# Submodule includes
# ============================================================================

# Foundational types and constants (no dependencies on other submodules)
include("Materials/Materials.jl")
using .Materials

# Material property evaluation (depends only on Materials)
include("Properties/Properties.jl")
using .Properties

# Constitutive physics and local numerical solves (Materials, Properties).
include("Physics/Physics.jl")
using .Physics

# Boundary conditions use the same qualified local solve as gas interfaces.
include("BoundaryConditions/BoundaryConditions.jl")
using .BoundaryConditions

# The face-based mesh, the state and its mapping, the physics kernels, the
# mesh closure, the residual, cell deactivation and the colored Jacobian.
# Reached by module path (`Pyrolysis.Meshes`, …); the user-facing names are
# exported below.
include("Meshes/Meshes.jl")
include("States/States.jl")
include("Kernels/Kernels.jl")
include("Closures/Closures.jl")
include("Assembly/Assembly.jl")
include("Deactivation/Deactivation.jl")
include("ColoredJacobian/ColoredJacobian.jl")

# The problem definition and its validation
include("Problem/Problem.jl")
using .Problem

# The solution and its accessors, and the time integration of a
# `PyrolysisProblem`
include("Outputs/Outputs.jl")
include("Integration/Integration.jl")

# Umbrella namespace re-exporting every submodule's public surface for advanced
# users (`using Pyrolysis.Internal`). Must load after every other submodule.
include("Internal.jl")

using .Meshes: create_uniform_mesh, create_stretched_mesh, create_mesh_from_nodes
using .Meshes: Planar, Cylindrical, Spherical, MeshTopology, MeshGeometry
using .Closures: NoCorrection, EquidistributionCorrection, ActivationTime
using .Outputs: PyrolysisSolution, SolverStatistics
using .Outputs:
    GasDiagnostics,
    get_porosity,
    get_pressure,
    get_partial_pressures,
    get_intrinsic_gas_densities,
    get_gas_composition,
    get_boundary_exchange,
    gas_species_audit
using .Outputs: get_times, get_temperatures, get_concentrations, get_gas_generation_rate
using .Outputs: n_timesteps, succeeded, completed
using .Outputs: get_profile, get_cell_positions, interpolate_solution
using .Outputs: ConservationRecorder, conservation_audit
using .Integration: solve, setup, integrate, CoreProblem, parameter_sweep
using .Integration: default_policy, rodas5p_policy, fbdf_policy, kencarp4_policy
using .Integration: StructuredKLU

# ============================================================================
# Public API surface.
# ============================================================================
#
# Everything not listed here is reachable as `Pyrolysis.<Submodule>.X` or
# through `Pyrolysis.Internal`.

# Materials
export Material, SolidComponent, GaseousComponent, LiquidComponent
export LumpedTransfer, DustyGas, LumpedEndpoint, DustyEndpoint
export phase, molar_mass, is_gas, is_solid, is_liquid
export Reaction, DepletionLimiter
export MixingRule, PARALLEL, SERIES, WEIGHTED, BRUGGEMAN, CARMAN_KOZENY
export MixingSpec, MaterialMixing

# Property functions
export ConstantProperty, LinearProperty, PolynomialProperty, SutherlandProperty
export ArrheniusProperty, TableProperty, CallableProperty, StateDependentProperty

# Boundary conditions
export HeatFluxBC,
    ConvectiveBC,
    RadiativeBC,
    RadiativeFluxBC,
    CombinedThermalBC,
    DirichletThermalBC,
    AdiabaticBC
export ImpermeableBC, MassFluxBC, GasReservoir, AmbientGasBC, MassTransferFilm
# Blowing models for ConvectiveBC and the gas data they read
export AbstractBlowingModel,
    NoBlowing, CouetteFilm, PropertyRatioFilm, PoolingFilm, FunctionBlowing
export GasSpecies, NITROGEN, nasa7_heat_capacity, eucken_viscosity, fuller_diffusivity
export BlowingFace, blowing_state, blowing_ratio, blowing_factor
export NeumannPressureBC, ConvectivePressureBC
export BoundaryConditionSet
# Background pore-air initialization for gas-bearing problems
export background_gas_concentration

# Radiation models
export RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT

# The mesh
export create_uniform_mesh, create_stretched_mesh, create_mesh_from_nodes
export Planar, Cylindrical, Spherical
export MeshTopology, MeshGeometry

# Problem & solve
export PyrolysisProblem, PyrolysisSolution, SolverStatistics
export solve, setup, integrate, CoreProblem
export parameter_sweep

# The mesh correction and the integrator
export NoCorrection, EquidistributionCorrection, ActivationTime
export default_policy, rodas5p_policy, fbdf_policy, kencarp4_policy
export StructuredKLU

# The solution's accessors and the conservation audit
export GasDiagnostics,
    get_porosity,
    get_pressure,
    get_partial_pressures,
    get_intrinsic_gas_densities,
    get_gas_composition,
    get_boundary_exchange,
    gas_species_audit
export get_times, get_temperatures, get_concentrations, get_gas_generation_rate
export n_timesteps, succeeded, completed
export get_profile, get_cell_positions, interpolate_solution
export ConservationRecorder, conservation_audit

# ============================================================================
# Precompile workload — solves the canonical configurations so that a
# session's first solve of such a problem reuses cached native code.
# ============================================================================

include("precompile.jl")

end # module Pyrolysis
