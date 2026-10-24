import Pyrolysis

@testset "depletion and trial-temperature contracts" begin
    @testset "finite rates and derivatives at absent reactants" begin
        for order in (0.0, 0.3, 0.5, 1.0, 2.0), x in (-1e6, -1.0, 0.0, 1e-12, 1.0)
            f(y) = reaction_rate(y, order, 700.0)
            @test isfinite(f(x))
            @test isfinite(ForwardDiff.derivative(f, x))
            x <= 0 && @test f(x) == 0
        end
        # Negative control: a positive base alone consumes absent reactants.
        for order in (0.0, 0.5, 1.0)
            @test log(2.0)^order > 0
        end
    end

    @testset "C1 onset and unchanged factor above the join" begin
        @test depletion_factor(0.0) == 0
        @test ForwardDiff.derivative(depletion_factor, 0.0) == 0
        for x in (1.0, 2.0, 10.0, 20.0)
            @test depletion_factor(x) == tanh(x)
        end
        @test ForwardDiff.derivative(depletion_factor, 1.0 - 1e-8) ≈ 1 / cosh(1.0)^2 atol =
            1e-7
        for x in range(0.0, 1.0; length = 101)
            @test 0 <= depletion_factor(x) <= 1
            @test ForwardDiff.derivative(depletion_factor, x) >= -1e-14
        end
        # Independent central difference, away from the piecewise boundaries.
        for x in (0.2, 0.7, 1.5)
            h = 1e-6
            fd = (depletion_factor(x + h) - depletion_factor(x - h)) / (2h)
            @test ForwardDiff.derivative(depletion_factor, x) ≈ fd rtol = 1e-8
        end
    end

    @testset "rates above onset agree with the current production kernel" begin
        for order in (0.0, 0.3, 0.5, 1.0)
            rxn = Pyrolysis.Reaction(
                :phase0_reference,
                1 => 2;
                A = 1.0,
                E = 1000 * Pyrolysis.Materials.R_GAS,
                h = 0.0,
                n = order,
            )
            for x in (1.0, 2.0, 20.0), T in (300.0, 700.0, 1500.0)
                expected = Pyrolysis.Physics.reaction_rate(rxn, (x, 0.0), T)
                @test reaction_rate(x, order, T) ≈ expected rtol = 1e-14
            end
        end
    end

    @testset "trial temperature is safe in representative kernels" begin
        kernels(T) = begin
            t = trial_temperature(T)
            [exp(-1000 / t), sqrt(t), log(t), t^1.75, t^4, t / 0.05]
        end
        for T in (-1e6, -1.0, 0.0, 199.0, 200.0, 200.5, 300.0, 2999.5, 3000.0, 1e6)
            @test 200 <= trial_temperature(T) <= 3000
            @test all(isfinite, kernels(T))
            @test all(isfinite, ForwardDiff.jacobian(x -> kernels(x[1]), [T]))
        end
        @test trial_temperature(700.0) == 700.0
        @test ForwardDiff.derivative(trial_temperature, 700.0) == 1
        # Negative controls: a capacity floor cannot repair either expression.
        @test !isfinite(ForwardDiff.derivative(T -> exp(-1000 / T), 0.0))
        @test !isfinite(exp(-1000 / -1.0))
    end

    @testset "absent reactant through an actual reacting step" begin
        for order in (0.0, 0.5, 1.0)
            # An inert carrier keeps the cell alive while A is absent; another
            # reactant C is present and reacts, so this is not a zero RHS test.
            function rhs!(du, u, p, t)
                a = reaction_rate(u[2], order, 700.0)
                c = reaction_rate(u[3], 1.0, 700.0)
                du .= (0.0, -a, -c, a + c)
            end
            sol = solve(
                ODEProblem(rhs!, [1.0, 0.0, 1.0, 0.0], (0.0, 1.0)),
                KenCarp4();
                reltol = 1e-9,
                abstol = 1e-12,
                isoutofdomain = (u, p, t) -> any(<(0), u),
            )
            @test SciMLBase.successful_retcode(sol)
            @test all(u -> u[2] == 0, sol.u)
            @test sol.u[end][3] < 1
            @test sol.u[end][4] > 0
            @test sum(sol.u[end]) ≈ 2 atol = 1e-12
        end
    end

    @testset "per-order approach to depletion" begin
        for order in (0.0, 0.5, 1.0)
            rhs!(du, u, p, t) = (du[1] = -reaction_rate(u[1], order, 700.0))
            sol = solve(
                ODEProblem(rhs!, [0.1], (0.0, 100.0)),
                KenCarp4();
                reltol = 1e-9,
                abstol = 1e-12,
                isoutofdomain = (u, p, t) -> any(<(0), u),
            )
            @test SciMLBase.successful_retcode(sol)
            @test all(u -> 0 <= u[1] <= 0.1, sol.u)
            @test sol.u[end][1] < 0.1
        end
    end
end
