"""
`ALEMeshVelocityOperator` correctness test.

**`apply!` equivalence.** Calling
`apply!(ws, ALEMeshVelocityOperator(), state, mesh, material, t)` populates
`(w_node, theta, theta_A)` exactly equal to the `_compute_ale_mesh_velocity!`
path inside `residual!` — the operator wraps the shared physics function.

The z-row prefix-sum blocks the operator embeds are gated against DenseAD by
the structured-backend tests.
"""

using Test
using Pyrolysis
using Pyrolysis.Jacobian: ALEMeshVelocityOperator, build_ale_mesh_velocity_workspace, apply!
using Pyrolysis.Residual: _compute_ale_mesh_velocity!

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

    @testset "apply! matches _compute_ale_mesh_velocity! (NoChi)" begin
        def, ws, u, _ = build_ale_no_radiation()
        material = def.material
        layout = ws.layout
        mesh = ws.mesh

        # Hot-but-realistic state: 700 K with virgin fraction 0.7,
        # gas fraction 0.05 — reactions active, θ non-zero.
        _set_ale_state!(u, layout, 700.0, (0.7, 0.05))
        Pyrolysis.Problem.populate!(ws.state_cache, u, layout)
        Pyrolysis.update_cross_section_area!(
            ws.state_cache,
            material,
            mesh,
            u,
            layout,
            Pyrolysis.Problem.NoChi(),
        )
        Pyrolysis.update_state_geometry!(ws.state_cache, mesh, Pyrolysis.Problem.ALE())
        Pyrolysis.Properties.update_properties!(
            ws.prop_cache,
            ws.state_cache,
            mesh,
            material,
        )

        # Reference: existing physics path used by residual!.
        w_ref = zeros(Float64, layout.n_nodes)
        θ_ref = zeros(Float64, layout.n_cells)
        _compute_ale_mesh_velocity!(
            w_ref,
            ws.state_cache,
            mesh,
            material;
            θ_out = θ_ref,
            θ_A = 0.0,
        )

        # Operator path: build workspace, apply!, read out.
        op_ws = build_ale_mesh_velocity_workspace(material, layout.n_cells, layout.n_nodes)
        apply!(op_ws, ALEMeshVelocityOperator(), ws.state_cache, mesh, material, 0.0)

        @test op_ws.w_node ≈ w_ref rtol = 1e-12
        @test op_ws.theta ≈ θ_ref rtol = 1e-12
        @test op_ws.theta_A[] == 0.0    # NoChi
    end
end
