include("nonlinear_caloric.jl")
using .NonlinearCaloric

@testset "nonlinear caloric inversion" begin
    model = CaloricModel((PolynomialLaw(700, 900, 1000), PeakLaw(800, 6000, 800, 12)))
    # Independent high-precision bisection: no Newton derivative or stopping rule.
    function oracle(H, m, model)
        setprecision(256) do
            lo, hi = big"-1e12", big"1e12"
            for _ = 1:240
                mid = (lo + hi) / 2
                if extended_energy(mid, BigFloat.(m), model) < BigFloat(H)
                    lo = mid
                else
                    hi = mid
                end
            end
            Float64((lo + hi) / 2)
        end
    end
    @testset "cp is the derivative of enthalpy" begin
        for law in model.laws, T in (-1000.0, 200.0, 300.0, 799.0, 800.0, 812.0, 3000.0)
            @test ForwardDiff.derivative(t -> enthalpy(law, t, model.Tref), T) ≈
                  heat_capacity(law, T) rtol = 1e-13
            @test heat_capacity(law, T) >= law.a
        end
    end
    max_iterations, max_error = 0, 0.0
    @testset "physical, absent, depleted, and negative trial inventories" begin
        masses = (
            [1.0, 0.0],
            [0.0, 1.0],
            [0.7, 0.3],
            [model.delta, 0.0],
            [0.5model.delta, -1.0],
            [0.0, 0.0],
            [-1.0, -2.0],
            [1e-18, 2e-18],
        )
        for m in masses, target in (-10000.0, -1.0, 0.0, 200.0, 300.0, 800.0, 3000.0, 1e6)
            H = extended_energy(target, m, model)
            report = inverse_report(H, m, model)
            reference = oracle(H, m, model)
            max_iterations = max(max_iterations, report.iterations)
            max_error = max(max_error, abs(report.T - reference) / max(1, abs(reference)))
            @test report.T ≈ reference atol = 2e-10 rtol = 2e-13
            @test report.T ≈ target atol = 2e-10 rtol = 2e-13
            @test extended_capacity(report.T, m, model) >= model.floor
            grad =
                ForwardDiff.gradient(v -> invert_temperature(v[1], v[2:3], model), [H; m])
            @test all(isfinite, grad)
            @test grad[1] ≈ 1 / extended_capacity(report.T, m, model) rtol = 1e-12
            @test all(
                isfinite,
                ForwardDiff.hessian(v -> invert_temperature(v[1], v[2:3], model), [H; m]),
            )
        end
        # This is an independent *physical* energy, without the extension.
        for m in ([1.0, 0.0], [0.0, model.delta], [model.delta, 100model.delta]),
            T in (200.0, 800.0, 3000.0)

            H = sum(m[i] * enthalpy(l, T, model.Tref) for (i, l) in enumerate(model.laws))
            @test regularization(m, model) == 0
            @test invert_temperature(H, m, model) ≈ T atol = 1e-10
        end
        @test_throws ErrorException inverse_report(1e7, [1.0, 0.0], model; maxiters = 1)
        @test_throws ArgumentError invert_temperature(Inf, [1.0, 0.0], model)
    end
    @testset "implicit first and second derivatives, including floor activation" begin
        # Nondimensional model makes finite differences well scaled near both
        # positive-part joints and the capacity-floor transition.
        scaled = CaloricModel(
            (PolynomialLaw(1, 2, 1), PeakLaw(2, 4, 1, 0.3));
            delta = 0.1,
            floor = 0.2,
            Tref = 0.0,
        )
        for m in ([-0.1, 0.0], [0.0, 0.15], [0.03, 0.08], [0.1, 0.2], [0.4, 0.7])
            H = extended_energy(1.3, m, scaled)
            f(v) = invert_temperature(v[1], v[2:3], scaled)
            v = [H; m]
            g = ForwardDiff.gradient(f, v)
            # At C1 joints a central difference has O(h), not O(h^2), error.
            step = 1e-8
            E = Matrix{Float64}(I, 3, 3)
            fd = [(f(v + step * E[:, j]) - f(v - step * E[:, j])) / (2step) for j = 1:3]
            @test g ≈ fd rtol = 2e-7 atol = 2e-6
            # The inverse identity must have zero first and second derivatives.
            identity(v) = extended_energy(f(v), v[2:3], scaled) - v[1]
            @test norm(ForwardDiff.gradient(identity, v), Inf) < 2e-12
            @test norm(ForwardDiff.hessian(identity, v), Inf) < 1e-10
        end
        # Identity masses alone do NOT guarantee physical agreement when the
        # capacity floor is active. This qualifies the document's exact region.
        floor_model =
            CaloricModel((PolynomialLaw(1, 0, 1),); delta = 1e-10, floor = 1e-6, Tref = 0.0)
        @test invert_temperature(1e-10, [1e-10], floor_model) ≈
              1e-10 / (1e-10 + regularization([1e-10], floor_model))
        @test regularization([1e-10], floor_model) > 0
        for mass in range(0, 3e-6; length = 31)
            @test regularization([mass], floor_model) >= 0
            @test extended_capacity(1.0, [mass], floor_model) >= floor_model.floor
        end
    end
    metric("nonlinear caloric: largest inverse iteration count", max_iterations)
    metric("nonlinear caloric: maximum relative error against 256-bit bisection", max_error)

    @testset "nonlinear inventory map inside an implicit, AD-driven reacting solve" begin
        T0, rate = 800.0, 20.0
        h = [enthalpy(l, T0, model.Tref) for l in model.laws]
        function rhs!(du, u, p, t)
            T = invert_temperature(u[3], view(u, 1:2), model)
            r = rate * u[1]
            du[1], du[2] = -r, r
            # Deliberate source gives an exact isothermal reacting trajectory.
            du[3] =
                r * (
                    enthalpy(model.laws[2], T, model.Tref) -
                    enthalpy(model.laws[1], T, model.Tref)
                )
        end
        prob = ODEProblem(rhs!, [1.0, 0.0, h[1]], (0.0, 0.2))
        run() = solve(prob, KenCarp4(); reltol = 1e-11, abstol = [1e-14, 1e-14, 1e-8])
        sol = run()
        @test SciMLBase.successful_retcode(sol)
        @test sol.u[end][1] ≈ exp(-4) rtol = 2e-7
        T_error = maximum(abs(invert_temperature(u[3], u[1:2], model) - T0) for u in sol.u)
        invariant_error = maximum(abs(u[3] - dot(h, u[1:2])) for u in sol.u)
        @test T_error < 1e-7
        @test invariant_error < 1e-6
        warmed = @timed run()
        metric("nonlinear caloric: reacting-solve maximum temperature error", T_error)
        metric(
            "nonlinear caloric: warmed solve seconds / allocated bytes",
            (warmed.time, warmed.bytes),
        )
    end
end
