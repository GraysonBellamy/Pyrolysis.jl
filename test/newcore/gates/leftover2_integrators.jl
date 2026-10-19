# Leftover 2 of the Phase 1 plan (§7.2): the integrator on the new residual.
#
# The candidates, each under its own policy (`kencarp4_policy`: raw estimate,
# flag kept, solver-chosen restart; `rodas5p_policy` with the solver-chosen and
# the forced restart; `fbdf_policy`: flag cleared, forced restart), on the
# finite-volume core with its colored Jacobian:
#
#   the audit decks `eulerian`, `ale`, `ale_darcy` at 25 and 100 cells (the
#   `problem_for` decks of the Phase 0 harness: 80 s, three outputs), at the
#   solve's default tolerances (abstol 1e-8 K, reltol 1e-6);
#   the event-dense case V-CPL-02 (Landau ablation, about sixty
#   deactivations), at its own tolerances (abstol 1e-7, reltol 1e-5).
#
# Reported per run: return code, the achieved error against a tight
# reference of the same deck (the largest temperature difference over the
# cells at the final time; for V-CPL-02 the final thickness and surface
# temperature, and its QoIs), the largest scaled constraint residual of the
# accepted states, steps, rejections, Jacobians, nonlinear failures and the
# warm median of three solves. The reference is KenCarp4 at abstol 1e-10,
# reltol 1e-8; Rodas5P at the same tolerances is the second reference, and
# their difference bounds what an error below it means. The old core's
# default solve of each deck is timed beside them, in the same session.
#
# The decision is the owner's (plan §10). `qualifies` applies a provisional
# rule for the table only: success, every constraint residual at or below
# its tolerance, and the error within ten times the requested relative
# tolerance of the reference's largest temperature (V-CPL-02: every QoI).
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/leftover2_integrators.jl
# Writes test/newcore/gates/measurements/leftover2_integrators.{csv,md}.

const _HARNESS =
    normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "mesh_phase0", "baselines"))
include(joinpath(_HARNESS, "jacobian_common.jl"))
include(joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))

using Dates
using LinearSolve: KLUFactorization
using Pyrolysis.Integration: kencarp4_policy, rodas5p_policy, fbdf_policy

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const REPEATS = 3
const DECKS = (
    (:eulerian, 25),
    (:eulerian, 100),
    (:ale, 25),
    (:ale, 100),
    (:ale_darcy, 25),
    (:ale_darcy, 100),
)
const CANDIDATES = (
    ("KenCarp4", () -> kencarp4_policy(KLUFactorization())),
    ("Rodas5P", () -> rodas5p_policy(KLUFactorization())),
    (
        "Rodas5P forced restart",
        () -> rodas5p_policy(KLUFactorization(); restart = :forced_dt),
    ),
    ("FBDF", () -> fbdf_policy(KLUFactorization())),
)
const TIGHT = (abstol = 1e-10, reltol = 1e-8)

# One warm-up solve, then the median of `REPEATS` timed solves.
function timed(f)
    sol = f()
    times = Float64[]
    for _ = 1:REPEATS
        GC.gc()
        push!(times, @elapsed sol = f())
    end
    return sol, median(times)
end

_final_T(sol) = sol.T[:, end]

function deck_rows(kind, n)
    problem, kw = problem_for(kind, n)
    new_solve(policy; tol...) = Pyrolysis.Integration.solve(
        problem;
        use_ale = kw.use_ale,
        kw.radiation_model,
        integrator_policy = policy,
        tol...,
    )
    reference = new_solve(kencarp4_policy(KLUFactorization()); TIGHT...)
    second = new_solve(rodas5p_policy(KLUFactorization()); TIGHT...)
    T_ref = _final_T(reference)
    spread = maximum(abs, _final_T(second) .- T_ref)
    budget = 10 * 1e-6 * maximum(abs, T_ref)
    rows = []
    for (name, make) in CANDIDATES
        sol, t = timed(() -> new_solve(make()))
        st = sol.statistics
        ok = sol.retcode === :Success
        err = ok ? maximum(abs, _final_T(sol) .- T_ref) : NaN
        push!(
            rows,
            (;
                deck = "$(kind) $(n)",
                candidate = name,
                retcode = sol.retcode,
                error = err,
                reference_spread = spread,
                constraint = st.largest_constraint_residual,
                steps = st.accepted,
                rejected = st.rejected,
                jacobians = st.jacobians,
                nonlinear_failures = st.nonlinear_failures,
                time = t,
                qualifies = ok && st.largest_constraint_residual <= 1 && err <= budget,
            ),
        )
    end
    old, t_old = timed(() -> solve(problem; kw...))
    push!(
        rows,
        (;
            deck = "$(kind) $(n)",
            candidate = "old core default",
            retcode = Symbol(old.retcode),
            error = NaN,
            reference_spread = spread,
            constraint = NaN,
            steps = old.raw_solution.stats.naccept,
            rejected = old.raw_solution.stats.nreject,
            jacobians = old.raw_solution.stats.njacs,
            nonlinear_failures = old.raw_solution.stats.nnonlinconvfail,
            time = t_old,
            qualifies = false,
        ),
    )
    return rows
end

function event_rows(case)
    new_solve(policy; tol...) = Pyrolysis.Integration.solve(
        case.build(case.ci.n_cells);
        use_ale = true,
        min_thickness = 1e-4,
        deactivation = true,
        integrator_policy = policy,
        abstol = 1e-7,
        reltol = 1e-5,
        tol...,
    )
    reference = new_solve(kencarp4_policy(KLUFactorization()); TIGHT...)
    second = new_solve(rodas5p_policy(KLUFactorization()); TIGHT...)
    L_ref, T_ref = reference.extras.thickness_history[end], reference.surface_T[end]
    spread = abs(second.surface_T[end] - T_ref)
    rows = []
    for (name, make) in CANDIDATES
        sol, t = timed(() -> new_solve(make()))
        st = sol.statistics
        ok = sol.retcode === :Success
        qois = ok ? Base.invokelatest(case.qois, sol) : []
        failed = ok ? count(q -> !(q.error < q.tol), qois) : -1
        err = ok ? abs(sol.surface_T[end] - T_ref) : NaN
        push!(
            rows,
            (;
                deck = "V-CPL-02 ($(length(sol.events)) events)",
                candidate = name,
                retcode = sol.retcode,
                error = err,
                reference_spread = spread,
                constraint = st.largest_constraint_residual,
                steps = st.accepted,
                rejected = st.rejected,
                jacobians = st.jacobians,
                nonlinear_failures = st.nonlinear_failures,
                time = t,
                qualifies = ok && failed == 0 && st.largest_constraint_residual <= 1,
                thickness_error = ok ? abs(sol.extras.thickness_history[end] - L_ref) : NaN,
                qoi_failures = failed,
            ),
        )
    end
    old, t_old = timed(
        () ->
            solve(case.build(case.ci.n_cells); progress = false, case.solve_kwargs...),
    )
    push!(
        rows,
        (;
            deck = "V-CPL-02",
            candidate = "old core default",
            retcode = Symbol(old.retcode),
            error = NaN,
            reference_spread = spread,
            constraint = NaN,
            steps = old.raw_solution.stats.naccept,
            rejected = old.raw_solution.stats.nreject,
            jacobians = old.raw_solution.stats.njacs,
            nonlinear_failures = old.raw_solution.stats.nnonlinconvfail,
            time = t_old,
            qualifies = false,
            thickness_error = NaN,
            qoi_failures = -1,
        ),
    )
    return rows
end

load_start = Sys.loadavg()
started = now()
rows = []
for (kind, n) in DECKS
    for r in deck_rows(kind, n)
        push!(rows, merge(r, (thickness_error = NaN, qoi_failures = -1)))
        println(join(("$k=$v" for (k, v) in pairs(rows[end])), " "))
        flush(stdout)
    end
end
# The case is loaded at top level, so the closures its file defines are
# callable from `event_rows`.
event_case = only(filter(c -> c.id == "V-CPL-02", load_verification_cases()))
for r in event_rows(event_case)
    push!(rows, r)
    println(join(("$k=$v" for (k, v) in pairs(r)), " "))
    flush(stdout)
end
load_end = Sys.loadavg()

open(joinpath(OUT, "leftover2_integrators.csv"), "w") do io
    csv_row(io, string.(keys(first(rows)))...)
    for row in rows
        csv_row(io, values(row)...)
    end
end

_g(x) = x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = 3))) : string(x)
open(joinpath(OUT, "leftover2_integrators.md"), "w") do io
    println(io, "# Leftover 2: the integrator candidates on the finite-volume core\n")
    println(
        io,
        "Produced by `test/newcore/gates/leftover2_integrators.jl`. Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. Error: the largest final-time temperature difference against ",
        "KenCarp4 at abstol 1e-10, reltol 1e-8 [K] (V-CPL-02: the final surface ",
        "temperature); spread: Rodas5P's at the same tight tolerances against that ",
        "reference. Constraint: the largest scaled constraint residual of the accepted ",
        "states (at most one passes). Warm medians of $(REPEATS) solves. `qualifies` is a ",
        "provisional rule for this table, not the decision (plan §10).\n",
    )
    println(
        io,
        "| deck | candidate | retcode | error | spread | constraint | steps | rejected | Jacobians | nonlinear failures | warm [s] | QoI failures | qualifies |",
    )
    println(io, "|" * repeat("---|", 13))
    for r in rows
        println(
            io,
            "| $(r.deck) | $(r.candidate) | $(r.retcode) | $(_g(r.error)) | $(_g(r.reference_spread)) | ",
            "$(_g(r.constraint)) | $(r.steps) | $(r.rejected) | $(r.jacobians) | $(r.nonlinear_failures) | ",
            "$(_g(r.time)) | $(r.qoi_failures < 0 ? "" : r.qoi_failures) | $(r.qualifies) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(now() - started, Second))",
)
