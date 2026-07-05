"""
Operator type declarations carved out at the framework level but whose
implementations live in dedicated files. Keeping the type + dependency
declarations (`reads` / `writes` / `stencil`) here lets the topology
sort and the operator-list builder see every operator as a single
include, and keeps the dependency graph readable in one place.

These are live operator declarations, not stubs — the `apply!` and
`contribute_to_jacobian!` implementations live in the files below.

# Status table

The table maps each operator tag type to its implementation file:

| Operator                  | Real impl file |
|---------------------------|----------------|
| `DarcyTransportOperator`  | `darcy_transport.jl` |
| `BCThermalOperator`       | `bc_thermal.jl`     |
| `BCMassOperator`          | `bc_mass.jl`        |
| `ALEAdvectionOperator`    | `ale_advection.jl`  |
| `ALEMeshVelocityOperator` | `ale_mesh_velocity.jl` |

The full `Structured()` backend walks this operator list directly.
"""

using ..Materials: Material

# ============================================================================
# ALEAdvectionOperator (real impl in `ale_advection.jl`)
# ============================================================================

"""
	ALEAdvectionOperator <: Operator

ALE advective contribution `(w · ∂T/∂x, w · ∂ξ/∂x)` via the limited
gradient. Reads `:w_cell` from `ALEMeshVelocityOperator` so the framework
composes the local face Jacobian with the cumulative-mesh-velocity
[`PrefixSumCoupling`](@ref).

`apply!`, `contribute_to_jacobian!`, and the AD workspace
[`build_ale_advection_workspace`](@ref) live in
[`operators/ale_advection.jl`](ale_advection.jl).
"""
struct ALEAdvectionOperator <: Operator end

reads(::ALEAdvectionOperator) = (:w_cell, :theta)
writes(::ALEAdvectionOperator) = ()
stencil(::ALEAdvectionOperator) = (cells = -2:2,)

# ============================================================================
# ALEMeshVelocityOperator (real impl in `ale_mesh_velocity.jl`)
# ============================================================================
#
# The type, workspace, apply!, and bare-primitive coupling builders live
# in `operators/ale_mesh_velocity.jl`. This operator's `contribute_to_jacobian!`
# is a no-op — every contribution is a structured coupling: the T/ξ-row chain
# rule is composed by the consumer ALEAdvectionOperator from the provider's
# bare primitives, while the z-row × (T, ξ) prefix-sum blocks (and the WithChi
# θ_A rank-1 channels) are embedded by the provider's own
# `contribute_to_structured!`.

struct ALEMeshVelocityOperator <: Operator end

reads(::ALEMeshVelocityOperator) = (:reaction_rate,)
writes(::ALEMeshVelocityOperator) = (:w_node, :w_cell, :theta, :theta_A)
stencil(::ALEMeshVelocityOperator) = (cells = 0:0, cumulative = true)
