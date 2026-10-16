"""
`ALEMeshVelocityOperator` correctness test.

**`apply!` equivalence.** Calling
`apply!(ws, ALEMeshVelocityOperator(), state, mesh, material, t)` populates
`(w_node, theta, theta_A)` exactly equal to the mesh velocity and strain rate
`residual!` computes — the operator calls the same physics functions.

The z-row prefix-sum blocks the operator embeds are gated against DenseAD by
the structured-backend tests.
"""

using Test
using Pyrolysis
using Pyrolysis.Jacobian: ALEMeshVelocityOperator, build_ale_mesh_velocity_workspace, apply!

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# Set a uniform (T, ξ) state with mild perturbations so reactions are
# active and θ is non-trivial.
function _set_ale_state!(u, layout, T_val, ξ_vals)
    n = layout.n_cells
    @inbounds for i = 1:n
        u[i] = T_val
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n + (j - 1) * n
        for i = 1:n
            u[off+i] = ξj
        end
    end
    return u
end

@testset "ALEMeshVelocityOperator" begin

    @testset "apply! matches the residual's mesh velocity (NoChi)" begin
        def, ws, u, _ = build_ale_no_radiation()
        material = def.material
        layout = ws.layout

        # Hot-but-realistic state: 700 K with virgin fraction 0.7,
        # gas fraction 0.05 — reactions active, θ non-zero.
        _set_ale_state!(u, layout, 700.0, (0.7, 0.05))
        residual!(similar(u), u, def, 0.0, ws)

        op_ws = build_ale_mesh_velocity_workspace(material, layout.n_cells, layout.n_nodes)
        apply!(op_ws, ALEMeshVelocityOperator(), ws.state_cache, ws.mesh, material, 0.0)

        @test any(!iszero, op_ws.theta)
        @test op_ws.w_node == ws.w_work
        @test op_ws.theta == ws.θ_work
        @test op_ws.theta_A[] == 0.0    # NoChi
    end
end
