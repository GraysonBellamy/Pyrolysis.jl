# Gate item 3 of the Phase 1 plan (§7.3): the manufactured-solution orders on
# both cores.
#
# Spatial: each MMS case's own convergence check (its error norm on its cell
# ladder at its convergence tolerances), on the old core and on the new one in
# both energy forms, with the observed orders between consecutive meshes.
# V-MMS-02 is a free-stream exactness case with no order; its drifts are its
# quantities of interest, listed with the others.
#
# Temporal: V-MMS-01 at its report sweep (80 cells, fixed steps 0.5 to
# 0.0625 s) on each core's default integrator, the L2 and L∞ of the final
# temperature against the same mesh at an eighth of the smallest step (the
# spatial error cancels), with the observed orders.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/mms_orders.jl
# Writes test/newcore/gates/measurements/mms_orders.md.

using Pyrolysis
using Dates

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const MMS_CASES = ("V-MMS-01", "V-MMS-02", "V-MMS-03", "V-MMS-04")

# The case definitions on one core and energy form: casedefs.jl reads both
# from the environment when it is included, so each is its own module.
function load_cases(core::String, energy::String)
    ENV["PYROLYSIS_CORE"] = core
    ENV["PYROLYSIS_ENERGY_STATE"] = energy
    m = Module(Symbol("cases_", core, "_", energy))
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(m, joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))
    cases = Base.invokelatest(m.load_verification_cases)
    return m, Dict(c.id => c for c in cases if c.id in MMS_CASES)
end

_orders(errors, ns) =
    [log(errors[i-1] / errors[i]) / log(ns[i] / ns[i-1]) for i = 2:length(ns)]

function spatial(m, case)
    conv = case.ci.convergence
    conv === nothing && return nothing
    build = conv.build === nothing ? case.build : conv.build
    errors = Float64[]
    for n in conv.n_cells
        problem = Base.invokelatest(build, n)
        s = Base.invokelatest(m._case_solve, case, problem; conv.solve_kwargs...)
        push!(errors, Base.invokelatest(conv.error_fn, s))
    end
    return (; n_cells = conv.n_cells, errors, orders = _orders(errors, conv.n_cells))
end

function qois(m, case)
    problem = Base.invokelatest(case.build, case.ci.n_cells)
    sol = Base.invokelatest(m._case_solve, case, problem; case.solve_kwargs...)
    return sol.retcode, Base.invokelatest(case.qois, sol)
end

function temporal(m, case)
    ts = case.report.temporal_sweep
    kw = haskey(case.report, :solve_kwargs) ? case.report.solve_kwargs : case.solve_kwargs
    solve_dt(dt) = Base.invokelatest(
        m._case_solve,
        case,
        Base.invokelatest(case.build, ts.n_cells);
        kw...,
        dt = dt,
    )
    ref = solve_dt(minimum(ts.dts) / 8)
    L2, Linf = Float64[], Float64[]
    for dt in ts.dts
        e = solve_dt(dt).T[:, end] .- ref.T[:, end]
        push!(L2, sqrt(sum(abs2, e) / length(e)))
        push!(Linf, maximum(abs, e))
    end
    dts = collect(Float64, ts.dts)
    order(err) = [log(err[i-1] / err[i]) / log(dts[i-1] / dts[i]) for i = 2:length(dts)]
    return (;
        n_cells = ts.n_cells,
        dts,
        L2,
        Linf,
        L2_orders = order(L2),
        Linf_orders = order(Linf),
    )
end

load_start = Sys.loadavg()
started = now()
runs = [("old", "temperature"), ("new", "temperature"), ("new", "enthalpy")]
spatial_rows, qoi_rows, temporal_rows = [], [], []
for (core, energy) in runs
    label = core == "old" ? "old core" : "new core, $(energy) form"
    m, cases = load_cases(core, energy)
    for id in MMS_CASES
        case = cases[id]
        r = spatial(m, case)
        r === nothing || push!(spatial_rows, (; id, label, r...))
        r === nothing || println(id, " ", label, " ", r)
        retcode, qs = qois(m, case)
        for q in qs
            push!(qoi_rows, (; id, label, retcode, q.name, q.error, q.tol))
        end
        flush(stdout)
    end
    r = temporal(m, cases["V-MMS-01"])
    push!(temporal_rows, (; label, r...))
    println("V-MMS-01 temporal ", label, " ", r)
    flush(stdout)
end
load_end = Sys.loadavg()

_g(x; d = 3) = x isa AbstractFloat ? string(round(x; sigdigits = d)) : string(x)
_gs(v; d = 3) = join(_g.(v; d), ", ")
open(joinpath(OUT, "mms_orders.md"), "w") do io
    println(io, "# Gate item 3: the manufactured-solution orders\n")
    println(
        io,
        "Produced by `test/newcore/gates/mms_orders.jl`. Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. Each core on its default integrator.\n",
    )
    println(io, "## Spatial orders (each case's convergence check)\n")
    println(io, "| case | run | cells | errors | observed orders |")
    println(io, "|---|---|---|---|---|")
    for r in spatial_rows
        println(
            io,
            "| $(r.id) | $(r.label) | $(join(r.n_cells, ", ")) | $(_gs(r.errors)) | $(_gs(r.orders)) |",
        )
    end
    println(io, "\n## V-MMS-01 temporal sweep\n")
    println(io, "| run | cells | dt [s] | L2 [K] | L∞ [K] | L2 orders | L∞ orders |")
    println(io, "|---|---|---|---|---|---|---|")
    for r in temporal_rows
        println(
            io,
            "| $(r.label) | $(r.n_cells) | $(_gs(r.dts)) | $(_gs(r.L2)) | $(_gs(r.Linf)) | ",
            "$(_gs(r.L2_orders)) | $(_gs(r.Linf_orders)) |",
        )
    end
    println(io, "\n## Quantities of interest at the CI resolution\n")
    println(io, "| case | run | retcode | quantity | error | tolerance | passes |")
    println(io, "|---|---|---|---|---|---|---|")
    for r in qoi_rows
        println(
            io,
            "| $(r.id) | $(r.label) | $(r.retcode) | $(r.name) | $(_g(r.error)) | $(_g(r.tol)) | $(r.error < r.tol) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(now() - started, Second))",
)
