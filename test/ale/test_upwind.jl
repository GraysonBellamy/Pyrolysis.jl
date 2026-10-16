"""
Manufactured solution test for ALE upwind gradient on non-uniform meshes.

Verifies that the upwind gradient uses center-to-center distance (not cell width).
With cell width, the gradient has a zeroth-order bias on non-uniform meshes:
  (T_i - T_{i-1})/Δz_i → (d/Δz_i)×T'_i ≠ T'_i  as h→0
Using center-to-center distance gives correct first-order upwind.
"""

using Test
using Pyrolysis
using Pyrolysis.Problem: StateCache
using Pyrolysis.Discretization: compute_ale_advection_term!
using Pyrolysis.Geometry: update_geometry!, CellState

@testset "ALE Upwind Gradient on Non-Uniform Mesh" begin

    NC = 2

    # Create a non-uniform mesh by modifying node positions of a uniform mesh
    # Pattern: alternating small (h) and large (2h) cells
    function make_nonuniform_mesh(n_pairs::Int, L::Float64 = 0.03)
        n_cells = 2 * n_pairs
        h = L / (3 * n_pairs)  # 3h per pair (h + 2h)

        # Start with a uniform mesh (provides correct topology)
        mesh = create_uniform_mesh(L, n_cells, Val(NC))

        # Reposition nodes to create alternating h, 2h cells
        z = 0.0
        mesh.nodes[1].z = z
        for pair = 1:n_pairs
            # After small cell (width h)
            z += h
            mesh.nodes[2*(pair-1)+2].z = z
            # After large cell (width 2h)
            z += 2h
            mesh.nodes[2*(pair-1)+3].z = z
        end

        # Recompute geometry from new node positions
        update_geometry!(mesh)

        # Set T(z) = z² at cell centers
        for i = 1:n_cells
            z_c = mesh.cell_centers[i]
            mesh.cell_states[i] = CellState{NC}(
                z_c^2,
                (500.0, 0.0),
                mesh.cell_volumes[i],
                mesh.cell_volumes[i],
            )
        end

        return mesh, h
    end

    @testset "Gradient accuracy with T(z) = z²" begin
        n_pairs = 10
        mesh, h = make_nonuniform_mesh(n_pairs)
        n_cells = length(mesh.cells)
        n_nodes = length(mesh.nodes)

        # Verify mesh is actually non-uniform
        widths = [mesh.cell_volumes[i] / mesh.cross_section_area for i = 1:n_cells]
        @test maximum(widths) / minimum(widths) > 1.9  # Should be ~2:1

        # Set uniform mesh velocity w < 0 (contracting)
        w_nodes = fill(-1e-5, n_nodes)

        # Compute ALE advection term w·∇T
        advection = zeros(n_cells)
        compute_ale_advection_term!(advection, StateCache(mesh), mesh, w_nodes)

        # Exact: T(z) = z², T'(z) = 2z, w·T'(z) = w × 2z
        w = -1e-5
        max_rel_error = 0.0
        checked = 0
        for i = 2:(n_cells-1)
            z_c = mesh.cell_centers[i]
            exact = w * 2.0 * z_c
            if abs(exact) > 1e-15
                rel_err = abs(advection[i] - exact) / abs(exact)
                max_rel_error = max(max_rel_error, rel_err)
                checked += 1
            end
        end

        # First-order upwind error = d/(2z) × T''/T', which is large near z≈0.
        # Check that the error DECREASES away from z=0 (not a constant bias).
        # A cell-width gradient would show ~50% error in ALL cells regardless
        # of position; here only cells near z≈0 have large error (inherent to
        # upwind).
        @test checked > 5

        # Check that cells far from z=0 have small error (< 10%)
        far_from_zero_error = 0.0
        for i = (n_cells÷2):(n_cells-1)
            z_c = mesh.cell_centers[i]
            exact = w * 2.0 * z_c
            if abs(exact) > 1e-15
                far_from_zero_error =
                    max(far_from_zero_error, abs(advection[i] - exact) / abs(exact))
            end
        end
        @test far_from_zero_error < 0.10
    end

    @testset "Convergence rate ≥ 1 (first-order upwind)" begin
        errors = Float64[]
        hs = Float64[]

        for n_pairs in [5, 10, 20, 40]
            mesh, h = make_nonuniform_mesh(n_pairs)
            n_cells = length(mesh.cells)
            n_nodes = length(mesh.nodes)

            w_nodes = fill(-1e-5, n_nodes)
            advection = zeros(n_cells)
            compute_ale_advection_term!(advection, StateCache(mesh), mesh, w_nodes)

            # Measure error in mid-domain (away from z≈0 where T'≈0)
            w = -1e-5
            max_err = 0.0
            for i = (n_cells÷4):(3*n_cells÷4)
                z_c = mesh.cell_centers[i]
                exact = w * 2.0 * z_c
                if abs(exact) > 1e-12
                    max_err = max(max_err, abs(advection[i] - exact) / abs(exact))
                end
            end

            push!(errors, max_err)
            push!(hs, h)
        end

        # Compute convergence rates between successive refinements
        rates = Float64[]
        for k = 2:length(errors)
            if errors[k] > 1e-15 && errors[k-1] > 1e-15
                rate = log(errors[k-1] / errors[k]) / log(hs[k-1] / hs[k])
                push!(rates, rate)
            end
        end

        # First-order upwind should give rate ≥ 1
        # A cell-width gradient gives rate ≈ 0 (the error doesn't decrease on a
        # fixed stretch ratio)
        @test !isempty(rates)
        avg_rate = sum(rates) / length(rates)
        @test avg_rate > 0.8
    end

    @testset "Uniform mesh regression" begin
        # On uniform mesh, cell width = center-to-center, so no change
        n_cells = 20
        L = 0.01
        mesh = create_uniform_mesh(L, n_cells, Val(NC))

        for i = 1:n_cells
            z_c = mesh.cell_centers[i]
            mesh.cell_states[i] = CellState{NC}(
                z_c^2,
                (500.0, 0.0),
                mesh.cell_volumes[i],
                mesh.cell_volumes[i],
            )
        end

        n_nodes = length(mesh.nodes)
        w_nodes = fill(-1e-5, n_nodes)
        advection = zeros(n_cells)
        compute_ale_advection_term!(advection, StateCache(mesh), mesh, w_nodes)

        w = -1e-5
        for i = (n_cells÷4):(n_cells-1)  # Skip cells near z=0 where T'≈0
            z_c = mesh.cell_centers[i]
            exact = w * 2.0 * z_c
            if abs(exact) > 1e-12
                rel_err = abs(advection[i] - exact) / abs(exact)
                @test rel_err < 0.15
            end
        end
    end
end
