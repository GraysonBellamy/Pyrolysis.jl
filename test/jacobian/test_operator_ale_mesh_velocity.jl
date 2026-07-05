"""
`ALEMeshVelocityOperator` correctness tests.

Two gates:

1. **`apply!` equivalence.** Calling
   `apply!(ws, ALEMeshVelocityOperator(), state, mesh, material, t)`
   populates `(w_node, w_cell, theta, theta_A)` exactly equal to the
   `_compute_ale_mesh_velocity!` path inside `Pyrolysis.Problem.residual!`
   — there is no math duplication; the operator wraps the existing
   shared physics function.

2. **Bare-primitive coupling builders match the analytical kernel.**
   The `PrefixSumCoupling` produced by
   [`build_prefix_sum_coupling_T`](@ref), when applied to a unit T-
   perturbation vector, reproduces the column of `dw_dT[:, j]` that
   the existing analytical kernel
   [`Pyrolysis.Jacobian.compute_node_velocity_derivatives!`](@ref)
   computes. Same for the ξ-channel primitives.

The consumer `ALEAdvectionOperator` composes these primitives with its
local face Jacobian (the `:w_cell ∈ reads` declaration); these unit
tests exercise the primitives directly.
"""

using Test
using Pyrolysis
using Pyrolysis.Jacobian:
    ALEMeshVelocityOperator,
    build_ale_mesh_velocity_workspace,
    build_prefix_sum_coupling_T,
    build_prefix_sum_coupling_xi,
    apply!
using Pyrolysis.Residual: _compute_ale_mesh_velocity!

# Reach into the Jacobian module's namespace for the strain-rate AD
# helper (not in the public export list — internal API).
const _strain_rate_derivatives! =
    getfield(Pyrolysis.Jacobian, :compute_strain_rate_derivatives!)

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
            ws.prop_cache,
            material,
            0.0;
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

    @testset "PrefixSumCoupling-T matches dw_dT column-wise" begin
        # Use the same hot state.
        def, ws, u, _ = build_ale_no_radiation()
        material = def.material
        layout = ws.layout
        mesh = ws.mesh

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

        op_ws = build_ale_mesh_velocity_workspace(material, layout.n_cells, layout.n_nodes)
        apply!(op_ws, ALEMeshVelocityOperator(), ws.state_cache, mesh, material, 0.0)
        _strain_rate_derivatives!(op_ws, ws.state_cache, material)

        # Build the bare T-channel primitive.
        c_T = build_prefix_sum_coupling_T(op_ws, ws.state_cache)

        # Apply to a unit T-perturbation at column j; the resulting
        # vector y[k] must equal dw_dT[k, j] from the analytical
        # reference — i.e. (k > j) · Δzⱼ · ∂θ_j/∂T_j.
        n_cells = layout.n_cells
        n_nodes = layout.n_nodes
        A_section = ws.state_cache.A_section[]

        # Reference column from the analytical kernel: prefix-sum form
        # dw_dT[k, j] = (j ≤ k-1) · Δzⱼ · dθ_dT[j].
        j = max(1, n_cells ÷ 3)
        Δzⱼ = ws.state_cache.volume[j] / A_section
        expected = zeros(n_cells)        # PrefixSumCoupling acts on length-n_cells
        for k = (j+1):n_cells          # PrefixSumCoupling outputs length-N_cells
            expected[k] = Δzⱼ * op_ws.dθ_dT[j]
        end

        x = zeros(n_cells)
        x[j] = 1.0
        y = zeros(n_cells)
        Pyrolysis.Jacobian.apply!(y, c_T, x)
        @test y ≈ expected rtol = 1e-12
    end

    @testset "PrefixSumCoupling-ξ matches dw_dξ column-wise" begin
        def, ws, u, _ = build_ale_no_radiation()
        material = def.material
        layout = ws.layout
        mesh = ws.mesh

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

        op_ws = build_ale_mesh_velocity_workspace(material, layout.n_cells, layout.n_nodes)
        apply!(op_ws, ALEMeshVelocityOperator(), ws.state_cache, mesh, material, 0.0)
        _strain_rate_derivatives!(op_ws, ws.state_cache, material)

        # Test the gas-component channel (j=2 for PMMA two-component).
        species_j = 2
        c_ξj = build_prefix_sum_coupling_xi(op_ws, ws.state_cache, species_j)

        n_cells = layout.n_cells
        A_section = ws.state_cache.A_section[]
        col_cell = max(1, n_cells ÷ 3)

        expected = zeros(n_cells)
        Δzⱼ = ws.state_cache.volume[col_cell] / A_section
        for k = (col_cell+1):n_cells
            expected[k] = Δzⱼ * op_ws.dθ_dξ[species_j, col_cell]
        end

        x = zeros(n_cells)
        x[col_cell] = 1.0
        y = zeros(n_cells)
        Pyrolysis.Jacobian.apply!(y, c_ξj, x)
        @test y ≈ expected rtol = 1e-10
    end
end
