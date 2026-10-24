isdefined(@__MODULE__, :NonlinearCaloric) || include("nonlinear_caloric.jl")
using .NonlinearCaloric

@testset "nonlinear temperature versus inventory integration" begin
    law = PolynomialLaw(700, 900, 1000)
    model = CaloricModel((law,))
    masses, T0, bath, conductance = [1.0], 1400.0, 300.0, 7000.0
    H0 = enthalpy(law, T0, model.Tref)
    # cp(T)*T' = -k*(T-bath) has an independent separation-of-variables
    # solution. No time integrator or caloric inversion supplies this oracle.
    primitive(T) =
        (law.a + law.b * (bath / law.scale)^2) * log(T - bath) +
        law.b / law.scale^2 * (T^2 / 2 + bath * T)
    function exact_temperature(t)
        target = primitive(T0) - conductance * t
        lo, hi = bath, T0
        for _ = 1:60
            mid = (lo + hi) / 2
            if primitive(mid) > target
                hi = mid
            else
                lo = mid
            end
        end
        return (lo + hi) / 2
    end
    function run(form; dt = nothing, tolerance = 1e-10)
        function rhs!(du, u, p, t)
            T = form == :temperature ? u[1] : invert_temperature(u[1], masses, model)
            q = conductance * (T - bath)
            du[1] = form == :temperature ? -q / heat_capacity(law, T) : -q
            du[2] = q # stage-level energy delivered to the bath
        end
        initial = [form == :temperature ? T0 : H0, 0.0]
        prob = ODEProblem(rhs!, initial, (0.0, 1.0))
        atol = [form == :temperature ? 1e-10 : 1e-7, 1e-7]
        if isnothing(dt)
            return solve(prob, KenCarp4(); reltol = tolerance, abstol = atol)
        end
        return solve(prob, KenCarp4(); reltol = 1e-12, abstol = atol, adaptive = false, dt)
    end
    function errors(sol, form)
        temperatures = [
            form == :temperature ? u[1] : invert_temperature(u[1], masses, model) for
            u in sol.u
        ]
        error =
            maximum(abs(T - exact_temperature(t)) for (T, t) in zip(temperatures, sol.t))
        ledger =
            maximum(
                abs(
                    (form == :temperature ? enthalpy(law, u[1], model.Tref) : u[1]) + u[2] -
                    H0,
                ) for u in sol.u
            ) / H0
        return (; error, ledger)
    end
    rows = NamedTuple[]
    for form in (:temperature, :inventory)
        refinement = Float64[]
        ledger_refinement = Float64[]
        for dt in (0.1, 0.05, 0.025)
            sol = run(form; dt)
            @test SciMLBase.successful_retcode(sol)
            e = errors(sol, form)
            push!(refinement, e.error)
            push!(ledger_refinement, e.ledger)
            form == :inventory && (@test e.ledger < 1e-12)
        end
        orders = log2.(refinement[1:2] ./ refinement[2:3])
        @test all(>(3.3), orders)
        if form == :temperature
            @test ledger_refinement[1] > 1e-8
            @test all(>(3.3), log2.(ledger_refinement[1:2] ./ ledger_refinement[2:3]))
        end
        metric("nonlinear $form: temperature errors / orders", (refinement, orders))
        metric("nonlinear $form: relative energy-ledger errors", ledger_refinement)
        for tolerance in (1e-6, 1e-8, 1e-10)
            sol = run(form; tolerance)
            @test SciMLBase.successful_retcode(sol)
            e = errors(sol, form)
            samples = [@timed(run(form; tolerance)) for _ = 1:3]
            sample = samples[sortperm([s.time for s in samples])[2]]
            push!(
                rows,
                (;
                    form,
                    tolerance,
                    temperature_error = e.error,
                    relative_ledger_error = e.ledger,
                    seconds = sample.time,
                    bytes = sample.bytes,
                    accepted = sol.stats.naccept,
                    rejected = sol.stats.nreject,
                    rhs_calls = sol.stats.nf,
                    jacobians = sol.stats.njacs,
                ),
            )
        end
    end
    # Both forms must offer a measured temperature-accuracy match. No speed gate:
    # timing this manufactured one-cell problem cannot settle production cost.
    for form in (:temperature, :inventory)
        @test any(r -> r.form == form && r.temperature_error < 1e-4, rows)
    end
    destination = joinpath(@__DIR__, "results", "caloric_comparison.csv")
    mkpath(dirname(destination))
    open(destination, "w") do io
        println(io, join(string.(keys(first(rows))), ','))
        for row in rows
            println(io, join(string.(values(row)), ','))
        end
    end
    metric("nonlinear state comparison CSV", destination)
end
