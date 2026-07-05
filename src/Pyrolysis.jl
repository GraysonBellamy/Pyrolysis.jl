"""
	Pyrolysis.jl

A Julia package for simulating thermal decomposition (pyrolysis) of materials
under high heat flux using an Arbitrary Lagrangian-Eulerian (ALE) adaptive
meshing framework.

# Main Features
- Multi-component materials with temperature-dependent properties
- Arrhenius-type decomposition kinetics (single and multi-step)
- Gas generation and transport through porous char
- In-depth radiation absorption (Beer-Lambert)
- Material volume changes (shrinkage, swelling, intumescence)
- Surface regression and internal interface motion
- Adaptive mesh refinement (h- and r-adaptivity)

# Quick Start
```julia
using Pyrolysis

# Define a simple material
pmma = Material(
	name = :PMMA,
	components = (
		SolidComponent(:virgin, ρ=1190, c=1420, k=0.21, ε=0.86),
		GaseousComponent(:MMA, M=0.100, c=1500, k=0.02, λ=1e-5)
	),
	reactions = (
		Reaction(:decomposition, :virgin => :MMA, A=8.5e12, E=188e3, h=870e3),  # h > 0 = endothermic
	)
)

# Create mesh and solve
mesh = create_uniform_mesh(0.01, 100, Val(2))  # 1 cm, 100 cells, 2 components
bcs = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
)

problem = PyrolysisProblem(
	mesh = mesh,
	material = pmma,
	bc_set = bcs,
	T_initial = z -> 300.0,
	ξ_initial = [z -> 1190.0, z -> 0.0],
	tspan = (0.0, 100.0),
)
solution = solve(problem; progress = false)
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

# Foundational types and constants (was Core; no dependencies on other submodules)
include("Materials/Materials.jl")
using .Materials

# Boundary-condition types (promoted out of Core; depends only on Materials)
include("BoundaryConditions/BoundaryConditions.jl")
using .BoundaryConditions

# Geometry and mesh infrastructure
include("Geometry/Geometry.jl")
using .Geometry

# Material property evaluation
include("Properties/Properties.jl")
using .Properties

# Physical models
include("Physics/Physics.jl")
using .Physics

# ALE mesh motion / GCL / remapping was folded into Physics during the reorg
# (the former `ALE` submodule no longer exists; its kernels live in Physics).

# Adaptive mesh refinement
include("Adaptivity/Adaptivity.jl")
using .Adaptivity

# Spatial discretization (finite volume)
include("Discretization/Discretization.jl")
using .Discretization

# Problem definition and setup (ProblemDef + Workspace + traits).
include("Problem/Problem.jl")
using .Problem

# Unified pure residual! (loaded after Problem; traits stay in Problem to avoid
# an upstream/downstream cycle).
include("Residual/Residual.jl")
using .Residual

# Operator-based Jacobian backend and solve strategies.
include("Jacobian/Jacobian.jl")
using .Jacobian

# Experimental namespace (research-grade P1 radiation, h-AMR) — load before Solver.
include("Experimental/Experimental.jl")

# First-class runtime conservation diagnostics (discrete + physical
# energy residuals; mass and reaction-balance ledgers).
include("Diagnostics/Diagnostics.jl")
using .Diagnostics

# Solution algorithms (consume Problem, Jacobian, Experimental, Diagnostics)
include("Solver/Solver.jl")
using .Solver

# Umbrella namespace re-exporting every submodule's public surface for advanced
# users (`using Pyrolysis.Internal`). Must load after every other submodule.
include("Internal.jl")

# ============================================================================
# Public API surface.
# ============================================================================
#
# Everything not listed here is reachable as `Pyrolysis.<Submodule>.X` or
# `Pyrolysis.Experimental.X`.

# Materials
export Material, SolidComponent, GaseousComponent, LiquidComponent
export phase, molar_mass, is_gas, is_solid, is_liquid
export Reaction, ReactionSet, DepletionLimiter
export MixingRule, PARALLEL, SERIES, WEIGHTED, BRUGGEMAN, CARMAN_KOZENY
export MixingSpec, MaterialMixing

# Property functions
export ConstantProperty, LinearProperty, PolynomialProperty
export ArrheniusProperty, TableProperty, CallableProperty, StateDependentProperty

# Boundary conditions
export HeatFluxBC,
    ConvectiveBC,
    RadiativeBC,
    RadiativeFluxBC,
    CombinedThermalBC,
    DirichletThermalBC,
    AdiabaticBC
export ImpermeableBC, MassFluxBC, ConvectiveMassBC, DiffusiveMassBC, CombinedMassBC
export DirichletPressureBC, NeumannPressureBC, ConvectivePressureBC
export BoundaryConditionSet
# Background pore-air initialization for DarcyFick-mode problems
export background_gas_concentration

# Radiation models (P1_QUASI_STEADY is exported from Pyrolysis.Experimental)
export RadiationModel, NO_RADIATION, SURFACE_ABSORPTION, BEER_LAMBERT

# Geometry
export create_uniform_mesh, create_stretched_mesh

# Problem & solve
export PyrolysisProblem, PyrolysisSolution
export solve
export parameter_sweep

# Pure residual + Workspace surface.
export SimulationMode
export StateLayout
export ProblemDef, Workspace
export to_problem_def, build_workspace, build_ad_workspace
export residual!
export state_length, pack_state!, update_mesh_from_state!, create_state_vector
export block_view
export get_temperature,
    get_concentration, get_node_position, get_node_positions, get_thickness_from_state
export column_chi_bar_from_state, cross_section_area_from_state
export MeshTopology, MeshGeometry

# Jacobian backends and harness. Re-exported at the top level so
# `using Pyrolysis` exposes them without forcing users to reach into
# `Pyrolysis.Jacobian`.
export JacobianBackend, JacobianSpec, compute_jacobian!
export Structured,
    LinearSolveStrategy, DirectSolve, StructuredSolve, ApproxSolve, MatrixFreeSolve
export DenseAD_ForwardDiff   # verifier reference only
export benchmark_jacobian_backends, verify_jacobian_at, is_approximate
export Operator,
    ResidualContext,
    StructuredJacobian,
    StructuredCoupling,
    PrefixSumCoupling,
    Rank1Coupling,
    SuffixSumCoupling,
    EmbeddedScaledPrefixSumCoupling,
    EmbeddedScaledSuffixSumCoupling,
    EmbeddedRank1Coupling,
    ImplicitTerm

# Linear-solver extensions (PyrolysisKLUExt, etc.) attach concrete
# factorizations to these stubs. Calling without the matching `using ...`
# errors via `MethodError`, prompting the user to load the dep.

"""
	klu_factorization()

KLU sparse-direct factorization for the integrator's linear solves, suitable
as the `linear_solver` keyword of `solve`. The concrete method is attached by
the `PyrolysisKLUExt` extension: load `KLU` first (`using KLU`); calling the
stub without it raises a `MethodError`.
"""
function klu_factorization end

"""
	pardiso_factorization()

MKL Pardiso sparse-direct factorization for the integrator's linear solves,
suitable as the `linear_solver` keyword of `solve`. The concrete method is
attached by the `PyrolysisPardisoExt` extension: load `Pardiso` first
(`using Pardiso`); calling the stub without it raises a `MethodError`.
"""
function pardiso_factorization end

"""
	amg_preconditioner(A)

Algebraic-multigrid preconditioner built from the sparse matrix `A`, for use
with the ILU/Krylov iterative linear solvers. The concrete method is attached
by the `PyrolysisAMGExt` extension: load `AlgebraicMultigrid` first
(`using AlgebraicMultigrid`); calling the stub without it raises a
`MethodError`.
"""
function amg_preconditioner end

export klu_factorization, pardiso_factorization, amg_preconditioner

# ============================================================================
# Sensitivity-analysis entry points — concrete methods are attached by
# `ext/PyrolysisSciMLSensitivityExt.jl`.
# ============================================================================

include("sensitivity_stubs.jl")

export forward_sensitivity, adjoint_sensitivity

# ============================================================================
# Symbolic-verification entry points (test-only) — concrete methods attached
# by `ext/PyrolysisSymbolicsExt.jl`.
# ============================================================================

include("symbolics_stubs.jl")

export verify_mixing_derivatives

# ============================================================================
# Precompile workload — warms up the canonical Material × SimulationMode
# combinations so TTFX drops to seconds for standard cases.
# ============================================================================

include("precompile.jl")

end # module Pyrolysis
