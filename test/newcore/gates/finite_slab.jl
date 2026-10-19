# Leftover 3 of the Phase 1 plan (§5.5): the finite-slab gate.
#
# V-HC-03 (constant flux into an adiabatic slab) and V-HC-02c (convective
# heating, Bi = 10) at 320 cells, at reltol/abstol 1e-6/1e-8 and 1e-8/1e-10:
# the three integrator candidates of the finite-volume core, each under its
# own policy (`kencarp4_policy`, `rodas5p_policy` with either restart,
# `fbdf_policy`), against the old core's default (`KenCarp4` with the
# smoothed estimate and the package's step-update callback, which raises the
# derivative-discontinuity flag), and its raw-estimate variant, in one session.
# Reported per run: steps, rejections, right-hand sides, Jacobians, the warm
# median of three solves, and the RMS (L2) and maximum (L∞) temperature error
# against the exact series at every output time.
#
# Pass criterion of §5.5, read at the spatial floor: for the default
# candidate and policy (`rodas5p_policy`), the largest L2 over the output
# times is no larger than the old core's default at the same tolerance or
# than the 320-cell spatial floor, whichever is larger. The floor is the
# largest L2 of the default at reltol 1e-11, abstol 1e-13, where the time error
# is gone; a run is at the floor within `FLOOR_MARGIN` of it (its time error
# at most that fraction of the spatial error). The old default's L2 can lie
# below the floor, its time error cancelling part of the spatial error. The
# script reports, for every candidate, whether it meets the criterion.
#
# The Phase 0 driver (scripts/mesh_phase0/baselines/integrator_flag.jl) writes
# into the Phase 0 results when it runs, so the old core is driven here
# through the harness's own probes (`timed_runs`) with the package's callback,
# its `current` variant, and the exact profiles are those of that driver.
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/finite_slab.jl
# Writes test/newcore/gates/measurements/finite_slab.{csv,md}.

# The harness's probes solve through the suite's core, which they read from
# PYROLYSIS_CORE when they are included: the old one here.
ENV["PYROLYSIS_CORE"] = "old"
const _HARNESS =
    normpath(joinpath(@__DIR__, "..", "..", "..", "scripts", "mesh_phase0", "baselines"))
include(joinpath(_HARNESS, "common.jl"))

using Dates
using Pyrolysis.Integration: kencarp4_policy, rodas5p_policy, fbdf_policy
import Pyrolysis.Outputs

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const N_CELLS = 320
const TOLERANCES = ((1e-6, 1e-8), (1e-8, 1e-10))   # (reltol, abstol)
const MAXITERS = 20_000
const REPEATS = 3
const FLOOR_TOLERANCE = (reltol = 1e-11, abstol = 1e-13)
const FLOOR_MARGIN = 0.01

# The exact profiles T(z, t) of the two cases (constants as in the case files).
function exact_profile(case_id::String)
    if case_id == "V-HC-03"
        ρ, c, k, L, T0, q = 1000.0, 1000.0, 0.5, 0.02, 300.0, 5000.0
        α = k / (ρ * c)
        return (z, t) -> T0 + (q * L / k) * slab_constflux_phi(z / L, α * t / L^2)
    elseif case_id == "V-HC-02c"
        ρ, c, k, L, T_i, T_inf, Bi = 1000.0, 1000.0, 0.5, 0.1, 300.0, 500.0, 10.0
        α = k / (ρ * c)
        λs = biot_eigenvalues(Bi, 60)
        return (z, t) ->
            T_inf + (T_i - T_inf) * slab_convective_theta(z / L, α * t / L^2, λs)
    end
    error("no exact profile for $case_id")
end

_positions(sol) =
    sol isa Outputs.CoreSolution ? Outputs.cell_positions(sol) : get_cell_positions(sol)

# RMS and maximum profile error at every output time after t = 0.
function profile_errors(sol, T_exact)
    z = _positions(sol)
    L2, Linf = Float64[], Float64[]
    for j = 2:length(sol.t)
        e = [sol.T[i, j] - T_exact(z[i], sol.t[j]) for i in eachindex(z)]
        push!(L2, sqrt(sum(abs2, e) / length(e)))
        push!(Linf, maximum(abs, e))
    end
    return L2, Linf
end

# A new-core solve: one warm-up, then the median of `REPEATS` timed solves.
function new_core_runs(case, policy, reltol, abstol)
    solve_once() = Pyrolysis.Integration.solve(
        case.build(N_CELLS);
        reltol,
        abstol,
        maxiters = MAXITERS,
        integrator_policy = policy,
    )
    sol = solve_once()
    times = Float64[]
    for _ = 1:REPEATS
        GC.gc()
        push!(times, @elapsed sol = solve_once())
    end
    return sol, median(times)
end

load_start = Sys.loadavg()
started = now()
cases =
    Dict(c.id => c for c in load_verification_cases() if c.id in ("V-HC-03", "V-HC-02c"))
rows = []
# The spatial floor of each case: the default at the floor tolerances.
floors = Dict(
    case_id => maximum(
        first(
            profile_errors(
                Pyrolysis.Integration.solve(
                    cases[case_id].build(N_CELLS);
                    FLOOR_TOLERANCE...,
                    maxiters = MAXITERS,
                ),
                exact_profile(case_id),
            ),
        ),
    ) for case_id in ("V-HC-03", "V-HC-02c")
)
for case_id in ("V-HC-03", "V-HC-02c"), (reltol, abstol) in TOLERANCES
    case = cases[case_id]
    T_exact = exact_profile(case_id)
    entries = Any[]
    for (name, make) in (
        (
            "old KenCarp4 (default)",
            () -> KenCarp4(autodiff = AutoFiniteDiff(), linsolve = KLUFactorization()),
        ),
        (
            "old KenCarp4 raw",
            () -> KenCarp4(
                smooth_est = false,
                autodiff = AutoFiniteDiff(),
                linsolve = KLUFactorization(),
            ),
        ),
    )
        r = timed_runs(
            case,
            N_CELLS;
            integrator = make(),
            reltol,
            abstol,
            maxiters = MAXITERS,
        )
        st = solve_stats(r.sol)
        L2, Linf =
            r.sol === nothing ? (Float64[NaN], Float64[NaN]) :
            profile_errors(r.sol, T_exact)
        push!(
            entries,
            (;
                name,
                retcode = r.retcode,
                steps = st.steps,
                rejects = st.rejects,
                nf = st.nf,
                njacs = st.njacs,
                time = r.warm_median_s,
                L2max = maximum(L2),
                Linfmax = maximum(Linf),
                L2final = L2[end],
            ),
        )
    end
    for (name, policy) in (
        ("new KenCarp4 raw", kencarp4_policy(KLUFactorization())),
        ("new Rodas5P (default)", rodas5p_policy(KLUFactorization())),
        (
            "new Rodas5P forced restart",
            rodas5p_policy(KLUFactorization(); restart = :forced_dt),
        ),
        ("new FBDF", fbdf_policy(KLUFactorization())),
    )
        sol, time = new_core_runs(case, policy, reltol, abstol)
        st = sol.statistics
        L2, Linf = profile_errors(sol, T_exact)
        push!(
            entries,
            (;
                name,
                retcode = sol.retcode,
                steps = st.accepted,
                rejects = st.rejected,
                nf = st.rhs_evaluations,
                njacs = st.jacobians,
                time,
                L2max = maximum(L2),
                Linfmax = maximum(Linf),
                L2final = L2[end],
            ),
        )
    end
    bound = max(first(entries).L2max, (1 + FLOOR_MARGIN) * floors[case_id])
    for e in entries
        row = (;
            case = case_id,
            reltol,
            abstol,
            e...,
            floor = floors[case_id],
            passes = e.L2max <= bound,
        )
        push!(rows, row)
        println(join(("$k=$v" for (k, v) in pairs(row)), " "))
        flush(stdout)
    end
end
load_end = Sys.loadavg()

open(joinpath(OUT, "finite_slab.csv"), "w") do io
    csv_row(io, string.(keys(first(rows)))...)
    for row in rows
        csv_row(io, values(row)...)
    end
end

_g(x; d = 3) = x isa AbstractFloat ? string(round(x; sigdigits = d)) : string(x)
open(joinpath(OUT, "finite_slab.md"), "w") do io
    println(io, "# Leftover 3: the finite-slab gate\n")
    println(
        io,
        "Produced by `test/newcore/gates/finite_slab.jl`. Julia $(VERSION), ",
        "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
        "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
        "at the end. $(N_CELLS) cells; warm medians of $(REPEATS) solves; L2 and L∞ the ",
        "temperature error against the exact series [K], the largest over the output ",
        "times. Spatial floor (the default at reltol $(FLOOR_TOLERANCE.reltol), abstol ",
        "$(FLOOR_TOLERANCE.abstol)): V-HC-03 $(_g(floors["V-HC-03"]; d = 6)) K, V-HC-02c ",
        "$(_g(floors["V-HC-02c"]; d = 6)) K. `passes`: L2 no larger than the old core's ",
        "default at the same tolerance or than the floor plus $(Int(100FLOOR_MARGIN)) %, ",
        "whichever is larger.\n",
    )
    println(
        io,
        "| case | reltol | run | retcode | steps | rejected | RHS | Jacobians | warm [s] | L2 max | L∞ max | passes |",
    )
    println(io, "|" * repeat("---|", 12))
    for r in rows
        println(
            io,
            "| $(r.case) | $(r.reltol) | $(r.name) | $(r.retcode) | $(r.steps) | $(r.rejects) | ",
            "$(r.nf) | $(r.njacs) | $(_g(r.time)) | $(_g(r.L2max; d = 6)) | $(_g(r.Linfmax; d = 6)) | $(r.passes) |",
        )
    end
end
println(
    "wrote $(OUT); load $(load_start) -> $(load_end); $(round(now() - started, Second))",
)
