# Shared helpers of the mesh-phase-0 baseline measurements.
#
# Run from the repository root with
#     julia --project=test --startup-file=no --threads=1 scripts/mesh_phase0/baselines/<driver>.jl
#
# Both drivers replace one solver callback in-process (method redefinition in
# `Pyrolysis.Solver`, never a source edit), solve verification cases with a
# wall-clock guard, and append one CSV row per configuration as soon as it is
# measured.

using Pyrolysis
using Pyrolysis.Solver: get_cell_positions
using Statistics: median
using Logging
using OrdinaryDiffEq: FBDF, Rodas5P, AutoFiniteDiff
using LinearSolve: KLUFactorization

# KenCarp4 lives in the SDIRK sublibrary, which the test environment reaches
# only through the package's own import.
const KenCarp4 = Pyrolysis.Solver.KenCarp4

const ROOT = normpath(joinpath(@__DIR__, "..", "..", ".."))
const RESULTS_DIR = joinpath(ROOT, "scripts", "mesh_phase0", "results")
mkpath(RESULTS_DIR)

# Verification-case definitions (Test-free) and the exact-solution library.
include(joinpath(ROOT, "test", "verification", "casedefs.jl"))

# ----------------------------------------------------------------------------
# Wall-clock guard
#
# Every solve runs under a deadline. The step-update callback (overridden by
# each driver) fires after every accepted step and terminates the integrator
# once the deadline has passed; the driver records the run as timed out. The
# iteration budget `maxiters` is the second guard.
# ----------------------------------------------------------------------------
const DEADLINE = Ref(Inf)
const TIMED_OUT = Ref(false)
# Integrator time at the last step-update fire: how far a run that stopped
# early got (the saved output times alone do not show it).
const LAST_T = Ref(NaN)
# Called right before every solve (drivers reset their per-solve counters).
const PRE_SOLVE_HOOK = Ref{Function}(() -> nothing)
const SOLVE_TIMEOUT_S = parse(Float64, get(ENV, "MESH_PHASE0_SOLVE_TIMEOUT_S", "300"))

# Discards every message below Warn without formatting it; warnings reach the
# `_VBOWarnCounter` wrappers of the burnout cases (they count before they
# forward) and are then discarded here.
struct QuietLogger <: Logging.AbstractLogger end
Logging.min_enabled_level(::QuietLogger) = Logging.Warn
Logging.shouldlog(::QuietLogger, args...) = true
Logging.catch_exceptions(::QuietLogger) = false
Logging.handle_message(::QuietLogger, args...; kwargs...) = nothing

"""
    guarded_solve(case, n_cells; timeout_s, kwargs...) -> NamedTuple

One solve of `case.build(n_cells)` through the case's own solve path under a
wall-clock deadline. Never throws: a failure is returned as `retcode =
:Error` with the message. `kwargs` are merged over `case.solve_kwargs`.
"""
function guarded_solve(case, n_cells::Int; timeout_s::Float64 = SOLVE_TIMEOUT_S, kwargs...)
    problem = case.build(n_cells)
    kw = merge(case.solve_kwargs, (; kwargs...))
    DEADLINE[] = time() + timeout_s
    TIMED_OUT[] = false
    LAST_T[] = NaN
    PRE_SOLVE_HOOK[]()
    GC.gc()
    timed = @timed try
        Logging.with_logger(QuietLogger()) do
            Base.invokelatest(_case_solve, case, problem; kw...)
        end
    catch err
        err
    end
    sol = timed.value
    if sol isa Exception
        return (;
            sol = nothing,
            retcode = :Error,
            timed_out = TIMED_OUT[],
            message = first(split(sprint(showerror, sol), '\n')),
            wall_s = timed.time,
        )
    end
    return (;
        sol,
        retcode = sol.retcode,
        timed_out = TIMED_OUT[],
        message = "",
        wall_s = timed.time,
    )
end

# Integrator counters of a finished solve (NaN-free integers; -1 when absent).
function solve_stats(sol)
    sol === nothing && return (;
        steps = -1,
        rejects = -1,
        nf = -1,
        njacs = -1,
        nw = -1,
        nsolve = -1,
        nonliniter = -1,
        nonlinfail = -1,
    )
    st = sol.raw_solution.stats
    return (;
        steps = st.naccept,
        rejects = st.nreject,
        nf = st.nf,
        njacs = st.njacs,
        nw = st.nw,
        nsolve = st.nsolve,
        nonliniter = st.nnonliniter,
        nonlinfail = st.nnonlinconvfail,
    )
end

"""
    timed_runs(case, n_cells; repeats = 3, kwargs...) -> NamedTuple

One warm-up solve (compilation) followed by `repeats` timed solves. The
repeats are skipped when the warm-up did not finish with the expected return
code. Returns the last solution, the first-call time, the warm median / min /
max, and whether every repeat reproduced the warm-up's step count.
"""
function timed_runs(case, n_cells::Int; repeats::Int = 3, kwargs...)
    warm = guarded_solve(case, n_cells; kwargs...)
    ok = warm.retcode == case.expected_retcode && !warm.timed_out
    if !ok
        return (;
            sol = warm.sol,
            retcode = warm.retcode,
            timed_out = warm.timed_out,
            message = warm.message,
            first_s = warm.wall_s,
            warm_median_s = NaN,
            warm_min_s = NaN,
            warm_max_s = NaN,
            repeats_consistent = false,
        )
    end
    ref_steps = solve_stats(warm.sol).steps
    times = Float64[]
    consistent = true
    last = warm
    for _ = 1:repeats
        r = guarded_solve(case, n_cells; kwargs...)
        if r.retcode != case.expected_retcode || r.timed_out
            return (;
                sol = r.sol,
                retcode = r.retcode,
                timed_out = r.timed_out,
                message = "repeat: " * r.message,
                first_s = warm.wall_s,
                warm_median_s = NaN,
                warm_min_s = NaN,
                warm_max_s = NaN,
                repeats_consistent = false,
            )
        end
        push!(times, r.wall_s)
        consistent &= solve_stats(r.sol).steps == ref_steps
        last = r
    end
    return (;
        sol = last.sol,
        retcode = last.retcode,
        timed_out = false,
        message = "",
        first_s = warm.wall_s,
        warm_median_s = median(times),
        warm_min_s = minimum(times),
        warm_max_s = maximum(times),
        repeats_consistent = consistent,
    )
end

# ----------------------------------------------------------------------------
# CSV output: one row appended and flushed per configuration.
# ----------------------------------------------------------------------------
csv_field(x::AbstractString) = '"' * replace(x, '"' => "\"\"") * '"'
csv_field(x::Bool) = string(x)
csv_field(x::Integer) = string(x)
csv_field(x::AbstractFloat) = isnan(x) ? "NaN" : repr(x)
csv_field(x) = csv_field(string(x))

function csv_row(io, vals...)
    println(io, join((csv_field(v) for v in vals), ","))
    flush(io)
end

function say(args...)
    println(args...)
    flush(stdout)
end
