@testset "solid-carried fields and ordered material frame" begin
    @testset "reaction on, transport off" begin
        function rhs!(du, u, p, t)
            du[1] = -0.8u[1]
            du[2] = (u[2] / u[1]) * du[1]
        end
        sol = solve(
            ODEProblem(rhs!, [1.0, 2.0], (0.0, 1.0)),
            KenCarp4();
            reltol = 1e-10,
            abstol = 1e-12,
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ exp(-0.8) rtol = 1e-8
        @test all(u -> isapprox(u[2] / u[1], 2.0; atol = 1e-12), sol.u)
        # Negative control: omit the reaction source from the packed field.
        @test abs(2 / sol.u[end][1] - 2) > 2
    end

    @testset "reaction off, transport on" begin
        # Two periodic upwind cells. The same solid-carried face flux advects
        # carrier and packed field, with an independent matrix-exponential answer.
        function rhs!(du, u, p, t)
            du[1] = u[2] - u[1]
            du[2] = -du[1]
            du[3] = u[4] - u[3]
            du[4] = -du[3]
        end
        sol = solve(
            ODEProblem(rhs!, [1.0, 2.0, 2.0, 8.0], (0.0, 1.0)),
            KenCarp4();
            reltol = 1e-10,
            abstol = 1e-12,
        )
        @test SciMLBase.successful_retcode(sol)
        expected = [1.5 - 0.5exp(-2), 1.5 + 0.5exp(-2), 5 - 3exp(-2), 5 + 3exp(-2)]
        @test sol.u[end] ≈ expected rtol = 1e-8
        @test sum(sol.u[end][1:2]) ≈ 3 atol = 1e-12
        @test sum(sol.u[end][3:4]) ≈ 10 atol = 1e-12
    end

    @testset "grain/radial/tangential axis priority" begin
        reference = [1 -1 0; 1 1 0; 0 0 sqrt(2)] / sqrt(2)
        F = Diagonal([2.0, 1.0, 1.0])
        axes = material_axes(F, reference)
        expected = [2 -1 0; 1 2 0; 0 0 sqrt(5)] / sqrt(5)
        @test axes ≈ expected atol = 1e-14
        @test axes' * axes ≈ Matrix{Float64}(I, 3, 3) atol = 1e-14
        @test det(axes) ≈ 1 atol = 1e-14
        @test atan(axes[2, 1], axes[1, 1]) ≈ atan(1 / 2) atol = 1e-14
        tensor = axes * Diagonal([4.0, 1.0, 1.0]) * axes'
        @test tensor[1:2, 1:2] ≈ [3.4 1.2; 1.2 1.6] atol = 1e-14
        # The former polar-only choice leaves the fibre at 45 degrees.
        @test abs(atan(axes[2, 1], axes[1, 1]) - pi / 4) > 0.3
        # Reversing the declared grain/radial priority is a different model.
        swapped = material_axes(F, reference[:, [2, 1, 3]])
        alternate = swapped * Diagonal([1.0, 4.0, 1.0]) * swapped'
        @test alternate[1:2, 1:2] ≈ [1.6 1.2; 1.2 3.4] atol = 1e-14
        @test norm(alternate - tensor) > 2

        rotation = [0.0 -1 0; 1 0 0; 0 0 1]
        rotated = material_axes(rotation, Matrix{Float64}(I, 3, 3))
        @test rotated * Diagonal([4.0, 1.0, 1.0]) * rotated' ≈ Diagonal([1.0, 4.0, 1.0]) atol =
            1e-14
        @test material_axes(rotation * F, reference) ≈ rotation * axes atol = 1e-14
        @test all(
            isfinite,
            ForwardDiff.jacobian(
                x -> vec(material_axes(reshape(x, 3, 3), reference)),
                vec(Matrix(F)),
            ),
        )
    end
end
