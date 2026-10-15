"""
Newton-failure step cap: after an accepted step whose Newton needed retries,
the proposed step stays at or below `growth` times that step for `steps`
accepted steps, then is free. A run without Newton failures is untouched.
"""

using Test
using Pyrolysis
using Pyrolysis.Solver:
    create_newton_failure_cap_callback,
    FAILURE_MEMORY_GROWTH,
    FAILURE_MEMORY_STEPS,
    DiscreteCallback,
    KenCarp4,
    get_proposed_dt,
    derivative_discontinuity!
using OrdinaryDiffEq: ODEProblem, CallbackSet, init, step!, solve

# Stiff scalar relaxation towards 1: after the transient the controller grows
# the step by the largest factor it allows. The first RHS call past `t_fail`
# returns NaN, which fails that stage's Newton once.
function _relaxation_problem(; t_fail = Inf)
    armed = Ref(isfinite(t_fail))
    f(u, p, t) =
        if armed[] && t > t_fail
            armed[] = false
            oftype(u, NaN)
        else
            -50.0 * (u - 1.0) * (1.0 + u^2)
        end
    return ODEProblem(f, 0.0, (0.0, 20.0))
end

# Mirrors `solve`, where the mesh-update callback fires after every step and
# marks a derivative discontinuity.
_noop_callback() = DiscreteCallback(
    (u, t, integ) -> true,
    integ -> derivative_discontinuity!(integ, true);
    save_positions = (false, false),
)

@testset "Newton-failure step cap" begin
    g, M = FAILURE_MEMORY_GROWTH, FAILURE_MEMORY_STEPS
    @test g == 1.25
    @test M == 20

    @testset "cap holds for M steps after a failure, then releases" begin
        integ = init(
            _relaxation_problem(t_fail = 0.2),
            KenCarp4();
            callback = CallbackSet(_noop_callback(), create_newton_failure_cap_callback()),
            abstol = 1e-8,
            reltol = 1e-6,
        )
        rows = NamedTuple{(:dt, :proposed, :nlfail),Tuple{Float64,Float64,Int}}[]
        while integ.t < 20.0 && length(rows) < 400
            step!(integ)
            push!(
                rows,
                (
                    dt = integ.t - integ.tprev,
                    proposed = get_proposed_dt(integ),
                    nlfail = integ.stats.nnonlinconvfail,
                ),
            )
        end
        i_fail = findfirst(r -> r.nlfail > 0, rows)
        @test i_fail !== nothing
        @test rows[end].nlfail == 1
        cap = g * rows[i_fail].dt
        # The failing step itself and the next M − 1 accepted steps propose at
        # most g × the failing step's accepted size.
        for i = i_fail:(i_fail+M-1)
            @test rows[i].proposed <= cap * (1 + 1e-12)
        end
        # The cap binds: without it the controller would have grown the step.
        @test any(i -> rows[i].proposed == cap, i_fail:(i_fail+M-1))
        # From the M-th step after the failure the step is free again.
        @test rows[i_fail+M].proposed > cap
    end

    @testset "no Newton failure: identical trajectory" begin
        kw = (; abstol = 1e-8, reltol = 1e-6)
        prob = _relaxation_problem()
        free = solve(prob, KenCarp4(); callback = CallbackSet(_noop_callback()), kw...)
        capped = solve(
            prob,
            KenCarp4();
            callback = CallbackSet(_noop_callback(), create_newton_failure_cap_callback()),
            kw...,
        )
        @test string(free.retcode) == "Success"
        @test free.stats.nnonlinconvfail == 0
        @test capped.stats.naccept == free.stats.naccept
        @test capped.t == free.t
        @test capped.u == free.u
    end
end
