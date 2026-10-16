"""
	Jacobian

Operator-based Jacobian assembly and solve strategies.

This module provides:
- [`JacobianBackend`](@ref) — abstract supertype of every assembly strategy
- [`JacobianSpec`](@ref) — the single user-facing handle, carried as a
  keyword argument to `solve(...)`
- [`compute_jacobian!`](@ref) — the dispatcher every backend hooks into via
  `_compute_jacobian!`
- [`FastSparsityMap`](@ref) — O(1) sparse-write helper, used by every
  backend that writes into a `SparseMatrixCSC` directly

The single production backend is `Structured` (`src/Jacobian/structured_backend.jl`);
`DenseAD_ForwardDiff` (`dense_ad_reference.jl`) is preserved exclusively
as the reference path used by `verify_jacobian_at`. Sparsity construction
(`build_sparsity` and the four `_emit_*_pattern!` helpers) lives in
`sparsity.jl`; verification (`verify_jacobian_at`) in `verify.jl`.
"""
module Jacobian

using SparseArrays
using ForwardDiff
using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Properties
using ..Physics
using ..Discretization
using ..Problem
using ..Problem:
    ProblemDef, Workspace, StateLayout, state_length, z_offset, chi_offset, chi_range
using ..Residual: residual!, _compute_theta_A, _compute_ale_mesh_velocity!

include("interface.jl")
include("fast_sparsity_map.jl")
include("sparsity.jl")

# ============================================================================
# Layered architecture — Layer 1 / 2 / 3.
#
# These files provide the operator interface, structured-Jacobian operator,
# structured-coupling types, and linear-solve strategies that the Jacobian
# path is built on.
# ============================================================================
include("strategies.jl")
include("structured_couplings.jl")
include("operators/Operator.jl")
include("structured.jl")
include("operators/reaction_source.jl")
include("operators/conduction.jl")
include("operators/gas_transport.jl")
include("operators/darcy_transport.jl")
include("operators/radiation.jl")
include("operators/bc_thermal.jl")
include("operators/bc_mass.jl")
include("operators/bc_convection.jl")
include("operators/ale_mesh_velocity.jl")
include("operators/chi.jl")
# Geometry-block AD is wrapped in the `GeometryOperator` interface; the
# `Structured` backend dispatches via
# `contribute_to_jacobian!(GeometryOperator(), …)`.
include("operators/geometry.jl")
include("structured_backend.jl")

# Verifier-only reference path. `DenseAD_ForwardDiff` runs ForwardDiff over
# the entire residual and returns a dense Jacobian; production callers use
# `Structured(...)`, the verifier uses this as a ground-truth comparison.
include("dense_ad_reference.jl")

include("verify.jl")

# ============================================================================
# Exports — public surface
# ============================================================================

# Interface
export JacobianBackend, JacobianSpec, compute_jacobian!

# Sparsity construction
export build_sparsity

# Verifier reference. Internal use only — production callers use
# `Structured(...)` directly. Surfaced here so `verify_jacobian_at`
# can construct a reference spec.
export DenseAD_ForwardDiff

# Rebuilt-architecture backend and its strategies — the production path.
export Structured
export LinearSolveStrategy, DirectSolve, ApproxSolve

# Operator interface (Layer 1) — visible so tests / downstream phases can
# define their own operators.
export Operator
export ReactionSourceOperator, build_reaction_source_workspace
export ConductionOperator, build_conduction_workspace
export GasTransportOperator, build_gas_transport_workspace
export RadiationOperator, build_beer_lambert_workspace
export DarcyTransportOperator, build_darcy_transport_workspace
export BCThermalOperator, build_bc_thermal_workspace
export BCMassOperator, build_bc_mass_workspace
export BCConvectionOperator, build_bc_convection_workspace
export ALEMeshVelocityOperator,
    build_ale_mesh_velocity_workspace, compute_theta_A_derivatives!
export ChiOperator, build_chi_workspace
export StructuredJacobian,
    StructuredCoupling,
    PrefixSumCoupling,
    EmbeddedScaledSuffixSumCoupling,
    EmbeddedNodePrefixSumCoupling,
    EmbeddedRank1Coupling

export verify_jacobian_at, is_approximate
export GeometryOperator, GeometryADCache, apply_geometry_jacobian!

# Sparse-write helper.
export FastSparsityMap, fast_sparse_add!

end # module Jacobian
