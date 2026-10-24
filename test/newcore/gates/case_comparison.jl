# The verification cases on both cores: every case of `NEW_CORE_CASES` at
# its CI resolution and solve keywords, the finite-volume core in the
# temperature form. Per case and core: each quantity of interest as the case
# reports it (its value, the exact or reference value, the error and the
# case's tolerance), the case's convergence spot check (its error norm at each
# resolution of the check and the observed orders), the integrator's counters,
# and the warm time of the case's solve.
#
# The cases come from `casedefs.jl` loaded into two modules, one routed to
# each core, so each core's quantities are its own. The quantities are paired
# by name; one a core does not report is left blank. The timing part times the
# case's own solve (its builder outside the timed region): the median of five
# solves after a warm-up, of three when the warm-up exceeds 5 s. The V-BO
# cases are not timed here: their solve evaluates the residual at every step
# on the finite-volume core, and `cutover_comparison.jl` times them.
#
# The accuracy parts do not depend on the load and run as parallel processes;
# the timing part runs alone, pinned to one core, at a one-minute load below
# 2. Run from the repository root, one thread each:
#   julia --project=ENV case_comparison.jl accuracy CORE DIR   (CORE = old, new)
#   taskset -c 1 julia --project=ENV case_comparison.jl timing DIR
#   julia --project=ENV case_comparison.jl merge DIR [OUT]
# `merge` writes test/newcore/gates/measurements/case_comparison.{csv,md}, or
# into OUT when given.

using Pyrolysis
using Logging
using Printf
using Serialization
using Statistics: median

include(joinpath(@__DIR__, "stamp.jl"))

const REPO = normpath(joinpath(@__DIR__, "..", "..", ".."))
const OUT = joinpath(@__DIR__, "measurements")
const CORES = ("old", "new")
const REPEATS = 5
const SLOW_WARMUP = 5.0

quiet(f) = Logging.with_logger(f, Logging.NullLogger())
latest(m, name, args...; kw...) =
    Base.invokelatest(Base.invokelatest(getglobal, m, name), args...; kw...)

# `casedefs.jl` loaded into a module routed to `core`.
const CASE_MODULES = Dict{String,Any}()
function cases_on(core)
    return get!(CASE_MODULES, core) do
        ENV["PYROLYSIS_CORE"] = core
        ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
        m = Module(Symbol("cases_", core))
        Core.eval(m, :(include(p) = Base.include($m, p)))
        Base.include(m, joinpath(REPO, "test", "verification", "casedefs.jl"))
        prefixes = getglobal(m, :NEW_CORE_CASES)
        list = filter(
            c -> any(p -> startswith(c.id, p), prefixes),
            latest(m, :load_verification_cases),
        )
        (; m, cases = Dict(c.id => c for c in list), ids = sort([c.id for c in list]))
    end
end

timed_case(id) = !startswith(id, "V-BO-")

function counters(sol)
    if sol isa Pyrolysis.Outputs.CoreSolution
        s = sol.statistics
        return (;
            steps = s.accepted,
            rejected = s.rejected,
            rhs = s.rhs_evaluations,
            jacobians = s.jacobians,
        )
    end
    raw = hasproperty(sol, :raw_solution) ? sol.raw_solution : nothing
    raw === nothing && return nothing
    s = raw.stats
    return (; steps = s.naccept, rejected = s.nreject, rhs = s.nf, jacobians = s.njacs)
end

# ----------------------------------------------------------------------------
# The accuracy part
# ----------------------------------------------------------------------------

function case_accuracy(cm, id)
    case = cm.cases[id]
    run_case(problem; kw...) = quiet(() -> latest(cm.m, :_case_solve, case, problem; kw...))
    problem = Base.invokelatest(case.build, case.ci.n_cells)
    sol = run_case(problem; case.solve_kwargs...)
    qois = [
        (; q.name, q.value, q.exact, q.error, q.tol, pass = q.error < q.tol) for
        q in Base.invokelatest(case.qois, sol)
    ]
    conv = case.ci.convergence
    convergence = nothing
    if conv !== nothing
        build = conv.build === nothing ? case.build : conv.build
        errors = Float64[]
        for n in conv.n_cells
            s = run_case(Base.invokelatest(build, n); conv.solve_kwargs...)
            push!(errors, Base.invokelatest(conv.error_fn, s))
        end
        orders = latest(cm.m, :observed_orders, errors, 1.0 ./ conv.n_cells)
        convergence = (; conv.n_cells, errors, orders, conv.expected_order)
    end
    return (;
        id,
        title = case.title,
        cells = case.ci.n_cells,
        retcode = Symbol(sol.retcode),
        expected_retcode = case.expected_retcode,
        qois,
        convergence,
        counters = counters(sol),
    )
end

function accuracy_part(core, dir)
    core in CORES || error("CORE must be old or new, got $(core)")
    started = time()
    load = Sys.loadavg()[1]
    cm = cases_on(core)
    results = Any[]
    for id in cm.ids
        t = @elapsed x = case_accuracy(cm, id)
        push!(results, x)
        println(
            "[$(id) $(core)] $(x.retcode), $(count(q -> q.pass, x.qois))/$(length(x.qois)) quantities, $(round(t, digits = 1)) s",
        )
        flush(stdout)
    end
    mkpath(dir)
    serialize(
        joinpath(dir, "cases_accuracy_$(core).jls"),
        (;
            core,
            results,
            stamp = COMMIT_STAMP,
            load = (load, Sys.loadavg()[1]),
            elapsed = time() - started,
        ),
    )
    println("wrote the $(core) accuracy part to $(dir); $(round(Int, time() - started)) s")
end

# ----------------------------------------------------------------------------
# The timing part
# ----------------------------------------------------------------------------

function case_time(cm, id)
    case = cm.cases[id]
    run() = begin
        problem = Base.invokelatest(case.build, case.ci.n_cells)
        stats = @timed quiet(
            () -> latest(cm.m, :_case_solve, case, problem; case.solve_kwargs...),
        )
        stats.time
    end
    warm_up = run()
    n = warm_up > SLOW_WARMUP ? 3 : REPEATS
    times = Float64[]
    for _ = 1:n
        GC.gc()
        push!(times, run())
    end
    return (; id, median = median(times), low = minimum(times), high = maximum(times), n)
end

# The CPUs this process may run on, as `taskset` set them.
function affinity()
    try
        return strip(last(split(read(`taskset -cp $(getpid())`, String), ":")))
    catch
        return "unknown"
    end
end

function timing_part(dir)
    started = time()
    load = Sys.loadavg()
    results = Dict{String,Any}()
    for core in CORES
        cm = cases_on(core)
        results[core] = Any[]
        for id in filter(timed_case, cm.ids)
            x = case_time(cm, id)
            push!(results[core], x)
            println("[$(id) $(core)] $(round(x.median, sigdigits = 3)) s")
            flush(stdout)
        end
    end
    serialize(
        joinpath(dir, "cases_timing.jls"),
        (;
            results,
            stamp = COMMIT_STAMP,
            load = (load, Sys.loadavg()),
            elapsed = time() - started,
            cpus = affinity(),
        ),
    )
    println("wrote the timing part to $(dir); $(round(Int, time() - started)) s")
end

# ----------------------------------------------------------------------------
# The merge
# ----------------------------------------------------------------------------

function check_stamp(stamp, what)
    stamp == COMMIT_STAMP ||
        error("$(what) was measured at $(stamp), not at $(COMMIT_STAMP)")
    return nothing
end

_g(x; d = 3) =
    x === nothing ? "" :
    x isa AbstractString ? x :
    !isfinite(x) ? (isnan(x) ? "" : string(x)) :
    x == 0 ? "0" : string(round(x; sigdigits = d))

# The ratio of the finite-volume core's error to the old core's on one
# quantity, where both errors lie above round-off (`ROUND_OFF` of the exact
# value, or absolute where the case gives none); `NaN` otherwise. Read on the
# exact-solution cases only: the burnout cases' quantities include step and
# failure budgets, and their accuracy is `cutover_comparison.jl`'s.
const ROUND_OFF = 1e-12
function error_ratio(o, n)
    (o === nothing || n === nothing) && return NaN
    floor = ROUND_OFF * (isfinite(o.exact) ? max(1.0, abs(o.exact)) : 1.0)
    (isfinite(o.error) && isfinite(n.error) && o.error > floor && n.error > floor) ||
        return NaN
    return n.error / o.error
end
exact_solution_case(id) = !startswith(id, "V-BO-")

geomean(v) = isempty(v) ? NaN : exp(sum(log, v) / length(v))
steps(c) = c === nothing ? "" : "$(c.steps) ($(c.rejected))"

function merge_parts(dir; out = OUT)
    acc = Dict{String,Any}()
    loads = Dict{String,Any}()
    for core in CORES
        p = deserialize(joinpath(dir, "cases_accuracy_$(core).jls"))
        check_stamp(p.stamp, "the $(core) accuracy part")
        acc[core] = Dict(x.id => x for x in p.results)
        loads[core] = p
    end
    timing = deserialize(joinpath(dir, "cases_timing.jls"))
    check_stamp(timing.stamp, "the timing part")
    times = Dict(core => Dict(x.id => x for x in timing.results[core]) for core in CORES)
    ids = sort(collect(intersect(keys(acc["old"]), keys(acc["new"]))))

    mkpath(out)
    open(joinpath(out, "case_comparison.csv"), "w") do io
        println(
            io,
            "case,quantity,old_value,new_value,exact,old_error,new_error,tol,old_pass,new_pass",
        )
        for id in ids
            o, n = acc["old"][id], acc["new"][id]
            names = unique(vcat([q.name for q in o.qois], [q.name for q in n.qois]))
            for name in names
                qo = findfirst(q -> q.name == name, o.qois)
                qn = findfirst(q -> q.name == name, n.qois)
                a = qo === nothing ? nothing : o.qois[qo]
                b = qn === nothing ? nothing : n.qois[qn]
                exact = a !== nothing ? a.exact : b.exact
                tol = a !== nothing ? a.tol : b.tol
                f(x, k) = x === nothing ? "" : string(getfield(x, k))
                println(
                    io,
                    join(
                        [
                            id,
                            "\"$(replace(name, "\"" => "'"))\"",
                            f(a, :value),
                            f(b, :value),
                            exact,
                            f(a, :error),
                            f(b, :error),
                            tol,
                            f(a, :pass),
                            f(b, :pass),
                        ],
                        ",",
                    ),
                )
            end
        end
    end

    open(joinpath(out, "case_comparison.md"), "w") do io
        println(io, "# The verification cases on both cores\n")
        println(
            io,
            "Produced by `case_comparison.jl` at $(COMMIT_STAMP). Julia $(VERSION), one thread per process. ",
            "Accuracy: the old core at a one-minute load of $(_g(loads["old"].load[1])), the finite-volume core at $(_g(loads["new"].load[1])). ",
            "Timing: one process on CPUs `$(timing.cpus)`, one-minute load $(_g(timing.load[1][1])) at the start and $(_g(timing.load[2][1])) at the end, ",
            "$(round(Int, timing.elapsed)) s. The script's header defines every measure.\n",
        )

        # The summary per case.
        println(io, "## Per case\n")
        println(
            io,
            "Quantities passing on each core; on the exact-solution cases, the geometric mean of the finite-volume core's error ",
            "over the old core's on the quantities both report with both errors above round-off ($(ROUND_OFF), relative to the ",
            "exact value where the case gives one; below one: the finite-volume core closer to the exact solution), and their ",
            "count (the burnout cases' accuracy is `cutover_comparison.md`'s); the convergence ",
            "check's observed orders; each core's accepted (rejected) steps and Jacobians (blank where a case's solve returns ",
            "no integrator statistics); the warm medians [s] and their ratio.\n",
        )
        println(
            io,
            "| case | cells | passing old, new | error ratio (n) | orders old | orders new | steps old | steps new | Jacobians old, new | old [s] | new [s] | ratio |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
        ratios_all, times_all = Float64[], Float64[]
        for id in ids
            o, n = acc["old"][id], acc["new"][id]
            r = Float64[]
            for qn in n.qois
                exact_solution_case(id) || break
                k = findfirst(q -> q.name == qn.name, o.qois)
                k === nothing && continue
                x = error_ratio(o.qois[k], qn)
                isfinite(x) && x > 0 && push!(r, x)
            end
            append!(ratios_all, r)
            ord(x) = x.convergence === nothing ? "" : join(_g.(x.convergence.orders), ", ")
            to, tn = get(times["old"], id, nothing), get(times["new"], id, nothing)
            tr = (to === nothing || tn === nothing) ? NaN : tn.median / to.median
            isfinite(tr) && push!(times_all, tr)
            println(
                io,
                "| $(id) | $(o.cells) | $(count(q -> q.pass, o.qois))/$(length(o.qois)), $(count(q -> q.pass, n.qois))/$(length(n.qois)) | ",
                "$(_g(geomean(r))) ($(length(r))) | $(ord(o)) | $(ord(n)) | ",
                "$(steps(o.counters)) | $(steps(n.counters)) | ",
                "$(o.counters === nothing ? "" : o.counters.jacobians), $(n.counters === nothing ? "" : n.counters.jacobians) | ",
                "$(to === nothing ? "" : _g(to.median)) | $(tn === nothing ? "" : _g(tn.median)) | $(_g(tr)) |",
            )
        end
        println(
            io,
            "\nOver the exact-solution cases: the geometric mean of the error ratios $(_g(geomean(ratios_all))) over $(length(ratios_all)) quantities ",
            "(median $(_g(isempty(ratios_all) ? NaN : median(ratios_all)))), and of the time ratios $(_g(geomean(times_all))) over ",
            "$(length(times_all)) cases (from $(_g(isempty(times_all) ? NaN : minimum(times_all))) to $(_g(isempty(times_all) ? NaN : maximum(times_all)))).\n",
        )

        # The convergence checks.
        println(io, "## The convergence checks\n")
        println(
            io,
            "Each case's own error norm at the check's resolutions on each core, and the ratio of the finite-volume core's to the old core's.\n",
        )
        println(io, "| case | cells | old errors | new errors | ratios | expected order |")
        println(io, "|---|---|---|---|---|---|")
        for id in ids
            o, n = acc["old"][id], acc["new"][id]
            (o.convergence === nothing || n.convergence === nothing) && continue
            println(
                io,
                "| $(id) | $(join(o.convergence.n_cells, ", ")) | $(join(_g.(o.convergence.errors), ", ")) | ",
                "$(join(_g.(n.convergence.errors), ", ")) | $(join(_g.(n.convergence.errors ./ o.convergence.errors), ", ")) | ",
                "$(_g(o.convergence.expected_order)) |",
            )
        end

        # Every quantity.
        println(io, "\n## Every quantity\n")
        println(
            io,
            "As each case reports it: the value on each core, the exact or reference value (blank where the case has none), each core's error and the case's tolerance; failing values in bold.\n",
        )
        println(io, "| case | quantity | old | new | exact | old error | new error | tol |")
        println(io, "|---|---|---|---|---|---|---|---|")
        for id in ids
            o, n = acc["old"][id], acc["new"][id]
            for name in unique(vcat([q.name for q in o.qois], [q.name for q in n.qois]))
                k = findfirst(q -> q.name == name, o.qois)
                j = findfirst(q -> q.name == name, n.qois)
                a = k === nothing ? nothing : o.qois[k]
                b = j === nothing ? nothing : n.qois[j]
                e(x) = x === nothing ? "" : x.pass ? _g(x.error) : "**$(_g(x.error))**"
                v(x) = x === nothing ? "" : _g(x.value; d = 6)
                ex = a !== nothing ? a.exact : b.exact
                tol = a !== nothing ? a.tol : b.tol
                println(
                    io,
                    "| $(id) | $(name) | $(v(a)) | $(v(b)) | $(_g(ex; d = 6)) | $(e(a)) | $(e(b)) | $(_g(tol)) |",
                )
            end
        end
    end
    println("wrote $(out)/case_comparison.{csv,md}")
end

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

function comparison_main(args)
    mode = isempty(args) ? "" : args[1]
    if mode == "accuracy" && length(args) == 3
        accuracy_part(args[2], args[3])
    elseif mode == "timing" && length(args) == 2
        timing_part(args[2])
    elseif mode == "merge" && length(args) in (2, 3)
        merge_parts(args[2]; out = length(args) == 3 ? args[3] : OUT)
    else
        error("usage: case_comparison.jl accuracy CORE DIR | timing DIR | merge DIR [OUT]")
    end
    return nothing
end

abspath(PROGRAM_FILE) == @__FILE__() && comparison_main(ARGS)
