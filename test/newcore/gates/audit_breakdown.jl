# Where the new core's time goes on the audit decks of gate item 2
# (`audit_decks.jl`), at matched accuracy.
#
# For each deck (the well-posed variants of `audit_decks.jl`) and each run at
# the tolerance the matched comparison selected: the warm time, the steps,
# right-hand sides, Jacobians, factorizations and linear solves, and the time
# of one residual, one colored Jacobian and one factorization of the
# iteration matrix with one solve (KLU, the integrator's linear solver), each
# a warm median in a function. The decomposition counts × unit costs against
# the warm time attributes the solve's time to its parts; the remainder is
# the integrator's own work, the domain test and the callbacks. Beside the
# default (`Rodas5P`), the two other candidates on the new core
# (`KenCarp4` with the raw and with the smoothed error estimate) and the old
# core's default, each at the tolerance that matches the old default's error
# on the deck (the cheapest of its ladder).
#
# Run from the repository root after `audit_decks.jl`:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/audit_breakdown.jl
# Writes test/newcore/gates/measurements/audit_breakdown.md.

ENV["PYROLYSIS_CORE"] = "old"
const _HARNESS =
    normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "mesh_phase0", "baselines"))
include(joinpath(_HARNESS, "jacobian_common.jl"))

using LinearSolve: LinearSolve, KLUFactorization, LinearProblem
using SparseArrays: SparseMatrixCSC
using LinearAlgebra: I
using OrdinaryDiffEq: AutoFiniteDiff
using Pyrolysis.Solver: KenCarp4
using Pyrolysis.Integration: IntegratorPolicy, kencarp4_policy, DT_AFTER_EVENT
import Pyrolysis.Assembly
import Pyrolysis.ColoredJacobian
import Pyrolysis.States

const OUT = joinpath(@__DIR__, "measurements")
const REPEATS = 3
const DECKS = (
    (:eulerian, 100, 1e-6),
    (:ale, 25, 1e-6),
    (:ale, 100, 1e-6),
    (:ale_darcy, 25, nothing),
    (:ale_darcy, 100, nothing),
)
const LADDER = (1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8)
const REFERENCE = 1e-11
_tol(reltol) = (reltol = reltol, abstol = reltol / 100)

smoothed_kencarp4_policy(linsolve) = IntegratorPolicy(
    KenCarp4(autodiff = AutoFiniteDiff(), linsolve = linsolve),
    :solver_chosen,
    :keep,
    DT_AFTER_EVENT,
)

function timed(f)
    f()
    times = Float64[]
    for _ = 1:REPEATS
        GC.gc()
        push!(times, @elapsed f())
    end
    return median(times)
end

# Median of many timings of `f()`, in a function so nothing is global.
function unit_time(f; samples = 51)
    f()
    return median([@elapsed(f()) for _ = 1:samples])
end

# One residual, one colored Jacobian, one numeric factorization of `M − γ J`
# with a solve (as the integrator refactors, the symbolic analysis kept) and
# one solve alone, at the deck's initial state.
function unit_costs(problem, kw)
    cp = Pyrolysis.Integration.setup(problem; use_ale = kw.use_ale, kw.radiation_model)
    u, du = copy(cp.u0), similar(cp.u0)
    J = ColoredJacobian.jac_prototype(cp.jacobian)
    t_rhs = unit_time(() -> Assembly.residual!(du, u, cp.def, 0.0, cp.ws))
    t_jac = unit_time(() -> ColoredJacobian.jacobian!(J, cp.jacobian, u, 0.0))
    M = all(States.is_differential(cp.def.layout)) ? I : States.mass_matrix(cp.def.layout)
    W = SparseMatrixCSC(M - 1e-3 * J)
    b = rand(length(u))
    cache = LinearSolve.init(LinearProblem(W, b), KLUFactorization())
    LinearSolve.solve!(cache)
    t_factor = unit_time(() -> (cache.A = W; LinearSolve.solve!(cache)))
    t_solve = unit_time(() -> (cache.b = b; LinearSolve.solve!(cache)))
    return (; t_rhs, t_jac, t_factor, t_solve, colors = ColoredJacobian.color_count(cp.def))
end

time_error(sol, ref) = maximum(abs, sol.T .- ref.T)

# The cheapest ladder tolerance whose error is at or below `target`.
function matched(solve_at, ref, target)
    for reltol in LADDER
        sol = solve_at(_tol(reltol))
        sol.retcode === :Success && time_error(sol, ref) <= target && return reltol, sol
    end
    return NaN, nothing
end

function deck_rows(kind, n, solid_λ)
    problem, kw = problem_for(kind, n; solid_λ)
    deck = solid_λ === nothing ? "$(kind) $(n)" : "$(kind) $(n), solid λ $(solid_λ)"
    old_at(tol) = solve(problem; kw..., tol...)
    new_at(policy) =
        tol -> Pyrolysis.Integration.solve(
            problem;
            use_ale = kw.use_ale,
            kw.radiation_model,
            integrator_policy = policy,
            tol...,
        )
    old_ref = old_at(_tol(REFERENCE))
    old = old_at(_tol(1e-6))
    target = time_error(old, old_ref)
    unit = unit_costs(problem, kw)
    rows = []
    s = old.raw_solution.stats
    push!(
        rows,
        (;
            deck,
            run = "old core default",
            reltol = 1e-6,
            error = target,
            time = timed(() -> old_at(_tol(1e-6))),
            steps = s.naccept,
            rejected = s.nreject,
            rhs = s.nf,
            jacobians = s.njacs,
            factorizations = s.nw,
            solves = s.nsolve,
            unit...,
        ),
    )
    new_ref = new_at(nothing)(_tol(REFERENCE))
    for (run, policy) in (
        ("new Rodas5P (default)", nothing),
        ("new KenCarp4 raw", kencarp4_policy(KLUFactorization())),
        ("new KenCarp4 smoothed", smoothed_kencarp4_policy(KLUFactorization())),
    )
        f = new_at(policy)
        reltol, sol = matched(f, new_ref, target)
        if sol === nothing
            push!(
                rows,
                (;
                    deck,
                    run,
                    reltol,
                    error = NaN,
                    time = NaN,
                    steps = 0,
                    rejected = 0,
                    rhs = 0,
                    jacobians = 0,
                    factorizations = 0,
                    solves = 0,
                    unit...,
                ),
            )
            continue
        end
        s = sol.raw_solution.stats
        push!(
            rows,
            (;
                deck,
                run,
                reltol,
                error = time_error(sol, new_ref),
                time = timed(() -> f(_tol(reltol))),
                steps = s.naccept,
                rejected = s.nreject,
                rhs = s.nf,
                jacobians = s.njacs,
                factorizations = s.nw,
                solves = s.nsolve,
                unit...,
            ),
        )
    end
    for r in rows
        println(join(("$k=$v" for (k, v) in pairs(r)), " "))
    end
    flush(stdout)
    return rows
end

include(joinpath(@__DIR__, "stamp.jl"))
load_start = Sys.loadavg()
started = time()
rows = []
for (kind, n, solid_λ) in DECKS
    append!(rows, deck_rows(kind, n, solid_λ))
end
load_end = Sys.loadavg()

_g(x; d = 3) =
    x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
open(joinpath(OUT, "audit_breakdown.md"), "w") do io
    println(io, "# Gate item 2: where the time goes at matched accuracy\n")
    println(
        io,
        "Produced by `test/newcore/gates/audit_breakdown.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. Each run at the cheapest tolerance of its ladder whose error is at ",
        "or below the old default's (error against the same core at reltol $(REFERENCE)). ",
        "Unit costs of the new core at the initial state [µs]: one residual, one colored ",
        "Jacobian, one KLU refactorization of the iteration matrix with a solve, one ",
        "solve; the old ",
        "core's row repeats them. Parts: Jacobians × Jacobian, right-hand sides × ",
        "residual, factorizations × factorization, solves × solve, each as a share of ",
        "the warm time; the remainder is the rest.\n",
    )
    println(
        io,
        "| deck | run | reltol | error [K] | warm [s] | steps | rejected | RHS | Jacobians | factorizations | solves | residual [µs] | Jacobian [µs] | factorization [µs] | solve [µs] | Jacobians share | RHS share | factorizations share | solves share | rest |",
    )
    println(io, "|" * repeat("---|", 20))
    for r in rows
        new = startswith(r.run, "new")
        parts =
            new && isfinite(r.time) ?
            (
                r.jacobians * r.t_jac,
                r.rhs * r.t_rhs,
                r.factorizations * r.t_factor,
                r.solves * r.t_solve,
            ) ./ r.time : (NaN, NaN, NaN, NaN)
        rest = new && isfinite(r.time) ? 1 - sum(parts) : NaN
        println(
            io,
            "| $(r.deck) | $(r.run) | $(r.reltol) | $(_g(r.error)) | $(_g(r.time)) | ",
            "$(r.steps) | $(r.rejected) | $(r.rhs) | $(r.jacobians) | $(r.factorizations) | ",
            "$(r.solves) | $(_g(1e6r.t_rhs)) | $(_g(1e6r.t_jac)) | $(_g(1e6r.t_factor)) | ",
            "$(_g(1e6r.t_solve)) | $(join(_g.(parts; d = 2), " | ")) | $(_g(rest; d = 2)) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
