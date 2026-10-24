@testset "prototype solver integration" begin
    @testset "weighted moving-mesh closure through AD and a DAE restart" begin
        base = structured_mesh((2, 1); lengths = [2.0, 1.0])
        top_nodes = [n for n in axes(base.x, 2) if base.x[2, n] > 0]
        free = 2 .* top_nodes
        theta = -0.2
        function geometry(u)
            x = Matrix{eltype(u)}(base.x)
            x[2, top_nodes] = u[1:3]
            mesh = (; base..., x)
            D, W = sweep_operator(mesh)
            volumes = [cell_volume(x, c, mesh.corners) for c in mesh.cells]
            return D[:, free], W[free], volumes
        end
        function consistent!(u)
            A, W, volumes = geometry(u)
            u[4:5] = (A * (A' ./ W)) \ (theta * volumes)
        end
        function rhs!(du, u, p, t)
            A, W, volumes = geometry(u)
            velocity = (A' * u[4:5]) ./ W
            du[1:3] = velocity
            du[4:5] = A * velocity - theta * volumes
        end
        residual(u) = (du = similar(u); rhs!(du, u, nothing, 0.0); du)
        initial = [1.0, 1.0, 1.0, 0.0, 0.0]
        consistent!(initial)
        @test norm(residual(initial)[4:5], Inf) < 1e-14
        jacobian = ForwardDiff.jacobian(residual, initial)
        h = 1e-6
        basis = Matrix{Float64}(I, 5, 5)
        finite_difference = hcat(
            [
                (
                    residual(initial + h * basis[:, j]) -
                    residual(initial - h * basis[:, j])
                ) / (2h) for j = 1:5
            ]...,
        )
        @test jacobian ≈ finite_difference atol = 1e-9 rtol = 1e-8
        function jump!(integrator)
            integrator.u[1:3] *= 0.8
            consistent!(integrator.u)
            derivative_discontinuity!(integrator, true)
        end
        callback = DiscreteCallback(
            (u, t, integ) -> t == 0.5,
            jump!;
            save_positions = (true, true),
        )
        f = ODEFunction(rhs!; mass_matrix = Diagonal([1.0, 1.0, 1.0, 0.0, 0.0]))
        sol = solve(
            ODEProblem(f, initial, (0.0, 1.0)),
            KenCarp4();
            callback,
            tstops = [0.5],
            initializealg = SciMLBase.NoInit(),
            reltol = 1e-9,
            abstol = 1e-11,
        )
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1:3] ≈ fill(0.8exp(theta), 3) rtol = 1e-7
        constraint_error = maximum(norm(residual(u)[4:5], Inf) for u in sol.u)
        @test constraint_error < 1e-8
        @test maximum(maximum(u[1:3]) - minimum(u[1:3]) for u in sol.u) < 1e-10
        metric("closure DAE: maximum accepted-state balance residual", constraint_error)
    end

    @testset "DAE consistent initialization and state-jump restart" begin
        # Differential y'= -y; algebraic z=y^2. The independent exact solution
        # is y=exp(-t), multiplied by 1/2 at t=0.5. This smoke test does not
        # replace Phase 0's stiff ladder or solver/bordered-solve comparison.
        function rhs!(du, u, p, t)
            du[1] = -u[1]
            du[2] = u[2] - u[1]^2
        end
        f = ODEFunction(rhs!; mass_matrix = Diagonal([1.0, 0.0]))
        function reinitialize!(integrator)
            integrator.u[1] *= 0.5
            integrator.u[2] = integrator.u[1]^2
            derivative_discontinuity!(integrator, true)
        end
        cb = DiscreteCallback(
            (u, t, integ) -> t == 0.5,
            reinitialize!;
            save_positions = (true, true),
        )
        sol = solve(
            ODEProblem(f, [1.0, 1.0], (0.0, 1.0)),
            KenCarp4();
            callback = cb,
            tstops = [0.5],
            initializealg = SciMLBase.NoInit(),
            reltol = 1e-10,
            abstol = 1e-12,
        )
        @test SciMLBase.successful_retcode(sol)
        residual = maximum(abs(u[2] - u[1]^2) for u in sol.u)
        @test residual < 1e-8
        @test sol.u[end][1] ≈ exp(-1) / 2 rtol = 1e-8
        @test sol.u[end][2] ≈ exp(-2) / 4 rtol = 1e-8
        metric("DAE: maximum accepted-state constraint residual", residual)
        # A state jump without its algebraic sweep is observably inconsistent.
        y, z = exp(-0.5), exp(-1)
        @test abs(z - (y / 2)^2) > 0.2
    end

    @testset "isothermal composition transport with common inventory fluxes" begin
        # Two periodic cells, two species with different cp. First-order face
        # states isolate flux compatibility; van-Albada reconstruction remains
        # a separate full-Phase-0 obligation.
        cp, T0 = [1.0, 3.0], 400.0
        initial_masses = [1.0 0.0; 0.0 1.0]
        H0 = vec(cp' * initial_masses) * T0
        initial = [vec(initial_masses); H0]
        function rhs!(du, u, p, t)
            masses = reshape(view(u, 1:4), 2, 2)
            T = [u[4+i] / dot(cp, masses[:, i]) for i = 1:2]
            for j = 1:2
                flux = masses[j, 2] - masses[j, 1]
                du[j] = flux
                du[j+2] = -flux
            end
            energy_flux =
                sum(cp[j] * (masses[j, 2] * T[2] - masses[j, 1] * T[1]) for j = 1:2)
            du[5] = energy_flux
            du[6] = -energy_flux
        end
        sol = solve(
            ODEProblem(rhs!, initial, (0.0, 1.0)),
            KenCarp4();
            reltol = 1e-10,
            abstol = 1e-12,
        )
        @test SciMLBase.successful_retcode(sol)
        max_temperature_error =
            maximum(abs(u[4+i] / dot(cp, u[(2i-1):(2i)]) - T0) for u in sol.u, i = 1:2)
        @test max_temperature_error < 1e-10
        @test all(u -> isapprox(sum(u[5:6]), sum(H0); atol = 1e-10), sol.u)
        metric("transport: maximum isothermal error", max_temperature_error)
        # Wrong energy flux (no energy transport) conserves total H but creates
        # a temperature change as the species composition changes.
        final_m1 = [(1 + exp(-2)) / 2, (1 - exp(-2)) / 2]
        @test abs(H0[1] / dot(cp, final_m1) - T0) > 100
    end
end
