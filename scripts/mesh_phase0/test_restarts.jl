isdefined(@__MODULE__, :SolverProbe) || include("solver_probe.jl")
isdefined(@__MODULE__, :RestartProbe) || include("restart_probe.jl")
using .RestartProbe
using SparseArrays
import LinearSolve

@testset "frequent jumps and measured colored AD" begin
    @testset "independent manufactured trajectory and one-sided jumps" begin
        for events in (0, 1, 10, 100), schedule in (:uniform, :clustered)
            p = JumpLadder(8, events; schedule)
            q = events == 0 ? 1.0 : 0.6^(1/events)
            for k = 0:events
                left = k == 0 ? 0.0 : p.times[k]
                right = k == events ? 1.0 : p.times[k+1]
                t = (left+right)/2
                # Direct sum of all decayed kicks, independent of the
                # precomputed segment-amplitude recurrence used by the RHS.
                fast =
                    0.1*q^k*exp(-200t) + sum(
                        0.05*(-1)^(j+1)*q^(k-j)*exp(-200(t-p.times[j])) for j = 1:k;
                        init = 0.0,
                    )
                exact = q^k .* (1 .+ 0.2 .* p.slow_mode .* exp(-t)) .+ fast .* p.fast_mode
                y, dy = event_state(p, t, k)
                @test y ≈ exact rtol = 2e-14 atol = 2e-14
                @test minimum(y) > 0.5
                p.segment = k
                u = vec(permutedims(hcat(y, SolverProbe.recurrence(y))))
                rhs = similar(u)
                SolverProbe.ladder_rhs!(rhs, u, p, t, true)
                @test rhs[1:2:end] ≈ dy atol = 2e-13
                @test norm(rhs[2:2:end], Inf) < 3e-16
            end
            for (k, t) in enumerate(p.times)
                before, after = event_state(p, t, k-1)[1], event_state(p, t, k)[1]
                @test after ≈ q .* before .+ (0.05*(-1)^(k+1)) .* p.fast_mode
            end
        end
        p, old = JumpLadder(8, 1), SolverProbe.Ladder(8)
        for (t, k) in ((0.0, 0), (0.49, 0), (0.5, 1), (0.8, 1))
            @test event_state(p, t, k)[1] ≈ SolverProbe.exact_state(old, t; after = k==1)[1]
        end
    end

    @testset "structural coloring and derivatives" begin
        for n in (3, 8, 32, 96), lifted in (false, true)
            cache = ColorCache(n, lifted)
            pattern, colors = cache.pattern, cache.colors
            for i in axes(pattern, 1)
                columns = findall(!iszero, Array(pattern[i, :]))
                @test length(unique(colors[columns])) == length(columns)
            end
            @test length(cache.columns) == (lifted ? 4 : n)
            # Check dense AD at small sizes, including zeros that must not
            # remove structurally possible recurrence couplings.
            n > 8 && continue
            p = JumpLadder(n, 10)
            for k in (0, 1, 10), trial in (:physical, :perturbed, :zero)
                p.segment = k
                t = k == 0 ? 0.0 : p.times[k]
                y = event_state(p, t, k)[1]
                u = lifted ? vec(permutedims(hcat(y, SolverProbe.recurrence(y)))) : y
                trial == :perturbed && (u .+= 0.2 .* sin.(eachindex(u)))
                trial == :zero && fill!(u, 0)
                residual(u) =
                    (du = similar(u); SolverProbe.ladder_rhs!(du, u, p, t, lifted); du)
                reference = ForwardDiff.jacobian(residual, u)
                colored = copy(pattern)
                colored_jac!(colored, u, p, t, lifted, cache)
                analytic = zeros(length(u), length(u))
                SolverProbe.ladder_jac!(analytic, u, p, t, lifted)
                @test Matrix(colored) ≈ reference rtol = 1e-13 atol = 1e-13
                @test analytic ≈ reference rtol = 1e-13 atol = 1e-13
                @test all(
                    iszero(reference[i, j]) || !iszero(pattern[i, j]) for
                    i in axes(reference, 1), j in axes(reference, 2)
                )
            end
            @test cache.passes == cache.calls * length(cache.columns)
        end
    end

    destination = joinpath(@__DIR__, "results")
    mkpath(destination)
    function write_csv(name, rows)
        path = joinpath(destination, name)
        open(path, "w") do io
            println(io, join(string.(keys(first(rows))), ','))
            for row in rows
                println(io, join(string.(values(row)), ','))
            end
        end
        metric("CSV", path)
    end

    costs = [jacobian_cost(n, form) for n in (8, 32, 96), form in (:lifted, :dense)]
    write_csv("colored_jacobian_cost.csv", vec(costs))
    for row in costs
        metric(
            "Jacobian n=$(row.n) $(row.form): colors / RHS equivalents",
            (row.colors, row.colored_rhs_equivalents),
        )
    end

    configurations = NamedTuple[]
    function add(n, form, events, tol; schedule = :uniform, reset_dt = false)
        for backend in (:analytic, :colored), method in (:KenCarp4_raw, :FBDF, :Rodas5P)
            push!(
                configurations,
                (; n, form, events, tol, schedule, reset_dt, backend, method),
            )
        end
    end
    # The primary comparison covers event density at the actual V-BO-02
    # tolerance pair and at 100x tighter tolerances. The 32-cell runs check a
    # second size; reduced runs measure the coloring penalty after elimination.
    for events in (0, 10, 30, 100, 120), tol in (1e-6, 1e-8)
        add(96, :lifted, events, tol)
    end
    for events in (0, 100)
        add(32, :lifted, events, 1e-6)
        for form in (:bordered, :dense)
            add(96, form, events, 1e-6)
        end
    end
    add(96, :lifted, 100, 1e-6; schedule = :clustered)
    add(96, :lifted, 100, 1e-6; reset_dt = true)
    rows = NamedTuple[]
    @testset "work precision with actual AD and consistent restarts" begin
        for c in configurations
            (; n, method, form, backend, events, tol, schedule, reset_dt) = c
            run(; verify = false) = run_restarts(
                n,
                method,
                form,
                backend,
                events;
                tolerance = tol,
                schedule,
                reset_dt,
                verify,
            )
            checked = run(; verify = true)
            @test SciMLBase.successful_retcode(checked.sol)
            @test checked.p.jumps == events
            @test checked.event_times == checked.p.times
            @test last(checked.sol.t) == 1.0
            @test all(u -> all(isfinite, u), checked.sol.u)
            @test maximum(checked.event_constraints; init = 0.0) < 1e-14
            @test backend != :colored || checked.colors.calls == checked.sol.stats.njacs
            @test checked.colors.passes ==
                  checked.colors.calls * length(checked.colors.columns)
            @test form != :bordered || checked.border.max_residual < 1e-10
            if method == :FBDF
                @test length(checked.restart_orders) == events
                @test all(==(1), checked.restart_orders)
            end
            # Compare a common achieved-error budget, and also report the
            # global error / requested rtol. Local tolerance is not a global
            # bound; retcode alone is never a qualification criterion.
            err = trajectory_error(checked, form)
            jerr = jump_error(checked, form)
            map_error = jump_error(checked, form; map_only = true)
            @test map_error < 2e-14
            constraint =
                form == :lifted ? SolverProbe.constraint_error(checked.sol, checked.p) : 0.0
            qualified = err < 1e-5 && jerr < 1e-5 && constraint < min(tol, 1e-7)
            # Whole-solve timings include preparation, callbacks and output;
            # exclude compilation. Error computation is outside the timing.
            samples = map(1:3) do _
                # Exclude collection of earlier configurations' saved
                # solutions. Collections incurred by this solve stay timed.
                GC.gc()
                @timed run()
            end
            sample = samples[sortperm([s.time for s in samples])[2]]
            result, stats = sample.value, sample.value.sol.stats
            push!(
                rows,
                (;
                    n,
                    method,
                    form,
                    backend,
                    events,
                    schedule,
                    reset_dt,
                    reltol = tol,
                    abstol = tol*0.01,
                    qualified,
                    retcode = string(result.sol.retcode),
                    max_error = err,
                    error_over_rtol = err/tol,
                    max_jump_error = jerr,
                    max_jump_map_error = map_error,
                    max_scaled_constraint = constraint,
                    seconds = sample.time,
                    seconds_min = minimum(s.time for s in samples),
                    seconds_max = maximum(s.time for s in samples),
                    gc_seconds = sample.gctime,
                    bytes = sample.bytes,
                    accepted = stats.naccept,
                    rejected = stats.nreject,
                    rhs_calls = stats.nf,
                    jacobians = stats.njacs,
                    colors = length(result.colors.columns),
                    dual_rhs_passes = result.colors.passes,
                    linear_solves = stats.nsolve,
                    border_factors = result.border.factorizations,
                    max_linear_residual = checked.border.max_residual,
                    min_bdf_order = isempty(result.order_samples) ? 0 :
                                    minimum(result.order_samples),
                    max_bdf_order = maximum(result.order_samples; init = 0),
                    first_order_steps = count(==(1), result.order_samples),
                    restart_count = length(result.restart_orders),
                ),
            )
            metric(
                "restarts n=$n $method $form $backend J=$events tol=$tol $schedule reset_dt=$reset_dt: error / seconds / accepted",
                (err, sample.time, stats.naccept),
            )
        end
        for events in (0, 10, 30, 100, 120),
            backend in (:analytic, :colored),
            method in (:KenCarp4_raw, :FBDF, :Rodas5P)

            pair = filter(
                r ->
                    r.n==96 &&
                    r.form==:lifted &&
                    r.events==events &&
                    r.backend==backend &&
                    r.method==method &&
                    r.schedule==:uniform &&
                    !r.reset_dt,
                rows,
            )
            @test length(pair) == 2
            @test pair[2].max_error < pair[1].max_error
        end
        for events in (0, 10, 30, 100, 120), tol in (1e-6, 1e-8)
            @test any(
                r ->
                    r.n==96 &&
                    r.form==:lifted &&
                    r.events==events &&
                    r.reltol==tol &&
                    r.backend==:colored &&
                    r.qualified,
                rows,
            )
        end
    end
    length(rows) == length(configurations) || error("incomplete restart comparison")
    write_csv("restart_comparison.csv", rows)
    open(joinpath(destination, "restart_environment.txt"), "w") do io
        println(io, "Julia=", VERSION, "; CPU=", Sys.CPU_NAME)
        for mod in
            (OrdinaryDiffEq, OrdinaryDiffEqSDIRK, LinearSolve, ForwardDiff, SciMLBase)
            println(io, nameof(mod), "=", pkgversion(mod))
        end
        for (_, mod) in Base.loaded_modules
            if nameof(mod) in (:DiffEqBase, :OrdinaryDiffEqBDF, :OrdinaryDiffEqCore)
                println(io, nameof(mod), "=", pkgversion(mod), "; source=", pathof(mod))
            end
        end
        println(
            io,
            "Julia_threads=",
            Threads.nthreads(),
            "; BLAS_threads=",
            BLAS.get_num_threads(),
        )
        println(
            io,
            "rows=",
            length(rows),
            "; time_statistic=median_of_three_warmed_whole_solves",
        )
        println(
            io,
            "GC before each timed solve; within-solve GC included; timing range reported",
        )
        println(
            io,
            "production_tolerances=V-BO-02 reltol 1e-6 abstol 1e-8; refinement=100x tighter",
        )
        println(
            io,
            "qualification=max_error and both one-sided jump errors <1e-5; scaled_constraint<min(reltol,1e-7)",
        )
        println(
            io,
            "errors sampled at all accepted states and both sides of every jump; dense interpolation not gated",
        )
        println(
            io,
            "lifted=UMFPACK; bordered=exact augmented UMFPACK with dense J; dense=LHL",
        )
        println(
            io,
            "AD=greedy structural coloring; one single-direction ForwardDiff RHS pass per color; included in jac!",
        )
        println(
            io,
            "rhs_calls=solver statistics excluding dual passes; dual_rhs_passes reported separately",
        )
        println(
            io,
            "event_map=q*y+alternating_0.05_fast_mode; q^J=0.6; fast_decay=200; nondimensional time interval [0,1]",
        )
        println(io, "uniform_events=k/(J+1); clustered_events=0.45+0.1*k/(J+1)")
        println(
            io,
            "reset_dt=true is a toy sensitivity, not a dimensional reproduction of DT_AFTER_PARK",
        )
        println(
            io,
            "observer is inert (tested against none); BDF orders read before jump; event initialization is explicit recurrence plus NoInit",
        )
        println(
            io,
            "observer comparison removes only exactly identical saved (t,u) pairs; both event limits retained",
        )
        println(
            io,
            "limits=manufactured equations; prescribed events; no root-location cost, topology/pattern changes or production physics",
        )
        println(
            io,
            "bordered retains dense J preparation; sparse augmented factorization is exact but not an optimized matrix-free implementation",
        )
        println(
            io,
            "colors are this prototype's counts, not production or multi-D counts; no production solver edits",
        )
    end
end
