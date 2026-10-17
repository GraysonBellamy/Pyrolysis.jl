@testset "weighted closure and geometric volume derivative" begin
    @testset "two-square negative control" begin
        mesh = structured_mesh((2, 1); lengths = [2.0, 1.0])
        D, W = sweep_operator(mesh)
        free = [2n for n in axes(mesh.x, 2) if mesh.x[2, n] > 0]
        b = ones(2)
        w_bad = closure_velocity(D, W, b, free; weighted = false)
        w = closure_velocity(D, W, b, free)
        @test w_bad[free] ≈ [2 / 3, 4 / 3, 2 / 3] atol = 1e-12
        @test D * w_bad ≈ b atol = 1e-12
        @test maximum(abs.(w_bad[free] .- 1)) > 0.3
        @test w[free] ≈ ones(3) atol = 1e-12
        metric("closure: unweighted top velocities", w_bad[free])
        metric("closure: weighted top velocities", w[free])
    end

    for counts in ((2, 1), (3, 2), (2, 1, 2), (2, 2, 2))
        @testset "connected $(length(counts))D mesh $counts" begin
            d = length(counts)
            mesh = structured_mesh(counts; lengths = Float64.(collect(counts)))
            D, W = sweep_operator(mesh)
            volumes = [cell_volume(mesh.x, cell, mesh.corners) for cell in mesh.cells]
            @test all(>(0), W)
            @test volumes ≈ ones(prod(counts)) atol = 1e-12
            @test sum(W) / d ≈ sum(volumes) atol = 1e-12
            free = [d * n for n in axes(mesh.x, 2) if mesh.x[d, n] > 0]
            @test rank(D[:, free]) == length(mesh.cells)
            w = closure_velocity(D, W, 0.7volumes, free)
            @test D * w ≈ 0.7volumes atol = 2e-12
            @test w[free] ≈ 0.7vec(mesh.x)[free] atol = 2e-12
            # A fully pinned body cannot satisfy nonzero expansion.
            @test_throws ArgumentError closure_velocity(D, W, volumes, Int[])

            # Skew the connected geometry. The independent volume derivative
            # must match face sweeps, not merely a self-consistent D solve.
            transform = Matrix{Float64}(I, d, d)
            transform[1, d] = 0.3
            skew = (; mesh..., x = transform * mesh.x)
            Ds, Ws = sweep_operator(skew)
            volume_map(x) =
                [cell_volume(reshape(x, d, :), c, mesh.corners) for c in mesh.cells]
            gradient = ForwardDiff.jacobian(volume_map, vec(skew.x))
            @test Ds ≈ gradient atol = 2e-12 rtol = 2e-12
            @test all(>(0), Ws)
            # Also deform a single vertex: this is not an affine image of the
            # Cartesian grid, and a hexahedron may have a nonplanar face.
            distorted_x = copy(skew.x)
            distorted_x[1, end] += 0.15
            distorted = (; mesh..., x = distorted_x)
            Dd, Wd = sweep_operator(distorted)
            distorted_gradient = ForwardDiff.jacobian(volume_map, vec(distorted_x))
            @test Dd ≈ distorted_gradient atol = 2e-12 rtol = 2e-12
            @test all(>(0), Wd)
            metric("closure: rank on $counts with column boundary set", rank(D[:, free]))
        end
    end
end
