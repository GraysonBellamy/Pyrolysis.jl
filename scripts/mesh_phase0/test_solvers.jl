include("solver_probe.jl")
using .SolverProbe
using SparseArrays
import LinearSolve

@testset "solver and exact bordered comparison" begin
    @testset "lifted and eliminated Jacobians against dense AD" begin
        for form in (:lifted, :bordered)
            prob = make_problem(6, form)
            p, u = prob.p, prob.u0
            residual(u) = (du = similar(u); ladder_rhs!(du, u, p, 0.0, form == :lifted); du)
            analytical = zeros(length(u), length(u))
            ladder_jac!(analytical, u, p, 0.0, form == :lifted)
            @test analytical ≈ ForwardDiff.jacobian(residual, u) rtol = 1e-13 atol = 1e-13
            trial = u .+ 0.1 .* sin.(eachindex(u))
            ladder_jac!(analytical, trial, p, 0.0, form == :lifted)
            @test analytical ≈ ForwardDiff.jacobian(residual, trial) rtol = 1e-13 atol =
                1e-13
        end
    end
    @testset "bordered Newton corrections equal a dense exact solve" begin
        p = Ladder(8)
        y = 0.5 .+ (1:8) ./ 8
        J = reduced_jacobian(y, p)
        p.jac_derivative .= 2 .* y ./ p.n
        p.jac_revision += 1
        cache = BorderCache(p.n)
        b = cos.(1.0:p.n)
        for (mu, factor) in ((1.0, -1e-5), (1.0, -0.1), (-100.0, 1.0), (1.0, -10.0))
            A = mu * I + factor * J
            answer = border_solve!(cache, A, b, zeros(p.n), p)
            @test answer ≈ A \ b rtol = 1e-11 atol = 1e-12
            @test norm(A * answer - b, Inf) < 1e-11
            count = cache.factorizations
            @test border_solve!(cache, A, 2b, zeros(p.n), p) ≈ 2answer rtol = 1e-12
            @test cache.factorizations == count
        end
        # A missed recurrence derivative silently changes the Newton correction.
        A = I - 0.1J
        missing = Matrix(A)
        missing[end, 1] = 0
        @test norm(A * (missing \ b) - b, Inf) > 1e-5
        # A new state/Jacobian at an unchanged shift must invalidate the factor.
        border_solve!(cache, A, b, zeros(p.n), p)
        count = cache.factorizations
        y = 0.6 .* y
        p.jac_derivative .= 2 .* y ./ p.n
        p.jac_revision += 1
        updated = I - 0.1reduced_jacobian(y, p)
        correction = border_solve!(cache, updated, b, zeros(p.n), p)
        @test cache.factorizations == count + 1
        @test norm(updated * correction - b, Inf) < 1e-11
    end
    @testset "masked-norm negative control" begin
        full, masked = solve_mask_counterexample(false), solve_mask_counterexample(true)
        error(sol) = maximum(abs(u[2]^2 - (1 + t)) for (u, t) in zip(sol.u, sol.t))
        @test SciMLBase.successful_retcode(full)
        @test SciMLBase.successful_retcode(masked)
        @test error(full) < 1e-9
        @test error(masked) > 1e-5
        metric(
            "constraint counterexample: full / masked norm residual",
            (error(full), error(masked)),
        )
    end

    # Timings include problem construction, Jacobians, solves, and callback
    # restart, exclude first-call compilation, and retain all accepted states.
    # All forms use analytic Jacobians checked against AD above. This is NOT a
    # colored-AD cost gate. The dense reference uses LinearSolve's exact LHL
    # reduction, not an approximation. The reduced prototype materializes J;
    # its bordered factorization is sparse but its total cost is not an optimal
    # matrix-free bordered implementation. Do not extrapolate a winner to 2D/3D.
    rows = NamedTuple[]
    @testset "stiff ladder and consistent restart" begin
        for n in (8, 32, 96),
            tol in (1e-6, 1e-8),
            method in (:KenCarp4, :KenCarp4_raw, :FBDF, :Rodas5P),
            form in (:lifted, :bordered, :dense)

            result = run_ladder(n, method, form; tolerance = tol, verify = true)
            @test SciMLBase.successful_retcode(result.sol)
            @test result.p.jumps == 1
            err = solution_error(result.sol, result.p, form == :lifted)
            constraint = form == :lifted ? constraint_error(result.sol, result.p) : 0.0
            @test all(u -> all(isfinite, u), result.sol.u)
            @test last(result.sol.t) == 1.0
            # Local step tolerances are not a bound on global solution error.
            # Rank only runs meeting the SAME physical and constraint budgets.
            # In particular, default KenCarp4 can report success while missing
            # this budget. The unsmoothed variant tests that sensitivity.
            qualified = err < 1e-5 && constraint < 1e-7
            form == :bordered && (@test result.cache.calls > 0)
            form == :bordered && (@test result.cache.max_residual < 1e-10)
            # Repeat only for measurement after this configuration has warmed.
            samples = [@timed(run_ladder(n, method, form; tolerance = tol)) for _ = 1:3]
            order = sortperm([s.time for s in samples])
            sample = samples[order[2]]
            stats = result.sol.stats
            push!(
                rows,
                (;
                    n,
                    method,
                    form,
                    tolerance = tol,
                    masked = false,
                    qualified,
                    max_error = err,
                    max_scaled_constraint = constraint,
                    seconds = sample.time,
                    bytes = sample.bytes,
                    accepted = stats.naccept,
                    rejected = stats.nreject,
                    rhs_calls = stats.nf,
                    jacobians = stats.njacs,
                    linear_solves = stats.nsolve,
                    border_factors = result.cache.factorizations,
                    max_linear_residual = result.cache.max_residual,
                    jumps = result.p.jumps,
                ),
            )
            metric(
                "ladder n=$n $method $form tol=$tol: error / constraint / seconds",
                (err, constraint, sample.time),
            )
        end
        # A benign ladder need not expose the norm bug; keep the decisive scalar
        # counterexample above as well. Record masked ladder runs, never rank them.
        for method in (:KenCarp4, :FBDF, :Rodas5P)
            run_ladder(8, method, :lifted; tolerance = 1e-8, masked = true)
            sample = @timed run_ladder(8, method, :lifted; tolerance = 1e-8, masked = true)
            result = sample.value
            @test SciMLBase.successful_retcode(result.sol)
            @test result.p.jumps == 1
            stats = result.sol.stats
            push!(
                rows,
                (;
                    n = 8,
                    method,
                    form = :lifted,
                    tolerance = 1e-8,
                    masked = true,
                    qualified = false,
                    max_error = solution_error(result.sol, result.p, true),
                    max_scaled_constraint = constraint_error(result.sol, result.p),
                    seconds = sample.time,
                    bytes = sample.bytes,
                    accepted = stats.naccept,
                    rejected = stats.nreject,
                    rhs_calls = stats.nf,
                    jacobians = stats.njacs,
                    linear_solves = stats.nsolve,
                    border_factors = 0,
                    max_linear_residual = 0.0,
                    jumps = result.p.jumps,
                ),
            )
        end
        for n in (8, 32, 96), form in (:lifted, :bordered, :dense)
            @test any(r -> r.n == n && r.form == form && r.qualified, rows)
        end
        for n in (8, 32, 96),
            method in (:KenCarp4, :KenCarp4_raw, :FBDF, :Rodas5P),
            form in (:lifted, :bordered, :dense)

            pair = filter(
                r -> r.n == n && r.method == method && r.form == form && !r.masked,
                rows,
            )
            @test pair[2].max_error < pair[1].max_error
        end
    end
    length(rows) == 75 || error("solver comparison incomplete; refusing a partial CSV")
    destination = joinpath(@__DIR__, "results", "solver_comparison.csv")
    mkpath(dirname(destination))
    open(destination, "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(string.(values(row)), ','))
        end
    end
    metric("solver comparison CSV", destination)
    open(joinpath(dirname(destination), "solver_environment.txt"), "w") do io
        println(io, "Julia=", VERSION)
        println(io, "OrdinaryDiffEq=", pkgversion(OrdinaryDiffEq))
        println(io, "OrdinaryDiffEqSDIRK=", pkgversion(OrdinaryDiffEqSDIRK))
        println(io, "LinearSolve=", pkgversion(LinearSolve))
        println(io, "CPU=", Sys.CPU_NAME)
        println(io, "Julia_threads=", Threads.nthreads())
        println(io, "BLAS_threads=", BLAS.get_num_threads())
        println(io, "physical_error_budget=1e-5")
        println(io, "scaled_constraint_budget=1e-7")
        println(io, "normal_rows=72; masked_negative_controls=3")
        println(io, "time_statistic=median_of_three_warmed_runs; compilation_excluded")
        println(io, "lifted=UMFPACK; bordered=augmented_UMFPACK; dense=LHL")
        println(io, "KenCarp4_raw=smooth_est_false")
        println(
            io,
            "Jacobians=analytic_checked_against_dense_AD; production_coloring_not_measured",
        )
    end
end
