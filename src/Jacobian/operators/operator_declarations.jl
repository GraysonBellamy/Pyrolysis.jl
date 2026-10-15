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
| `ALEMeshVelocityOperator` | `ale_mesh_velocity.jl` |

The full `Structured()` backend walks this operator list directly.
"""

using ..Materials: Material

# ============================================================================
# ALEMeshVelocityOperator (real impl in `ale_mesh_velocity.jl`)
# ============================================================================
#
# The type, workspace, apply!, and bare-primitive coupling builders live
# in `operators/ale_mesh_velocity.jl`. This operator's `contribute_to_jacobian!`
# is a no-op — every contribution is a structured coupling: the z-row ×
# (T, ξ) prefix-sum blocks (and the WithChi θ_A rank-1 channels) are
# embedded by the provider's own `contribute_to_structured!`. The T and ξ
# rows carry no mesh-motion coupling: every component rides the mesh, and
# the dilation `−ξⱼθᵢ` is cell-local (assembled from this operator's cached
# strain-rate derivatives by the backend).

struct ALEMeshVelocityOperator <: Operator end

reads(::ALEMeshVelocityOperator) = (:reaction_rate,)
writes(::ALEMeshVelocityOperator) = (:w_node, :w_cell, :theta, :theta_A)
stencil(::ALEMeshVelocityOperator) = (cells = 0:0, cumulative = true)
