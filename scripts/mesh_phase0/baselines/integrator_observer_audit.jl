# Integrator ranking on the performance-audit decks with an observer
# step-update callback.
#
# The audit's integrator ranking (results/newton_baseline.md, from the
# gitignored docs/performance_audit/solves.csv) ran every integrator with the
# package's step-update callback, which flags a derivative discontinuity after
# every accepted step; a multistep method (FBDF) then restarts at order one on
# every step, so its rows measure a crippled configuration. This driver
# replaces the step-update callback in-process (same technique as
# integrator_flag.jl) and re-runs the audit decks with:
#
#   KenCarp4      flag set (`current`), default smoothed error estimate
#   KenCarp4_raw  flag set, smooth_est = false
#   FBDF          flag set (`current`, as the audit ran it)
#   FBDF          observer: derivative_discontinuity!(integrator, false)
#   Rodas5P       observer
#
# Decks: `problem_for(kind, n)` copied from scripts/performance_audit/audit.jl
# (eulerian, ale, ale_darcy at 25 and 100 cells, 80 s, 3 outputs, the deck's
# own solve kwargs); Jacobian strategy and linear solver are the package
# defaults (Eulerian: Structured(DirectSolve) + KLU; ALE: Structured(ApproxSolve)
# + KLU). Tolerances abstol 1e-8, reltol 1e-6 (the `solve` defaults, as in the
# audit). Accuracy: maximum absolute temperature difference over cells at the
# final time against a reference solve of the same deck with KenCarp4 (flag
# set) at abstol 1e-10, reltol 1e-8, the audit's convention. Warm time: median
# of 3 solves after one warm-up. Every solve runs under `maxiters` and a
# wall-clock guard.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase0/baselines/integrator_observer_audit.jl
#
# Outputs (scripts/mesh_phase0/results/): integrator_observer_audit.csv (one
# row per deck × configuration, plus one reference row per deck, written
# incrementally) and integrator_observer_audit.md.

include(joinpath(@__DIR__, "common.jl"))

using Printf

const STEP_VARIANT = Ref(:current)
const MAXITERS = 200_000
const T_END = 80.0
const TOL = (abstol = 1e-8, reltol = 1e-6)
const REF_TOL = (abstol = 1e-10, reltol = 1e-8)

# In-process replacement of the step-update callback (same signature as the
# package method; `solve` builds it at every call, so the variant is read at
# solve time). The only addition to the package code is the wall-clock guard.
@eval Pyrolysis.Solver function create_step_update_callback(ws, material)
    mesh = ws.mesh
    sync_mesh! = FunctionWrapper{Nothing,Tuple{_StateVector}}(
        u -> begin
            layout = ws.layout
            mesh.cross_section_area =
                cross_section_area_from_state(mesh, material, u, layout)
            update_mesh_from_state!(mesh, u, layout)
        end,
    )
    condition = (u, t, integrator) -> true
    variant = Main.STEP_VARIANT[]
    deadline = Main.DEADLINE
    timed_out = Main.TIMED_OUT
    last_t = Main.LAST_T
    function affect!(integrator)
        sync_mesh!(integrator.u)
        last_t[] = integrator.t
        if time() > deadline[]
            timed_out[] = true
            terminate!(integrator)
        end
        if variant === :current
            derivative_discontinuity!(integrator, true)
        elseif variant === :observer
            derivative_discontinuity!(integrator, false)
        end
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

# Copy of `problem_for` (scripts/performance_audit/audit.jl). `solid_λ`
# overrides the solid's gas-transfer coefficient (the PMMA deck has none, the
# wood deck 1e-5); `lateral` gives the material a lateral shrinkage law (the
# WithChi mode); `flux` is the incident radiative flux (W/m²).
function problem_for(
    kind,
    n;
    t_end = T_END,
    outputs = 3,
    solid_λ = nothing,
    lateral = nothing,
    flux = 50e3,
)
    wood = kind in (:ale_darcy, :eulerian_darcy)
    mat = if wood
        Material(
            name = :audit_wood,
            components = (
                SolidComponent(
                    :wood,
                    ρ = 500.0,
                    ρ_intrinsic = 1500.0,
                    c = 1500.0,
                    k = 0.15,
                    ε = 0.9,
                    α = 5.0,
                    κ = 1e-14,
                    λ = something(solid_λ, 1e-5),
                ),
                GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5),
            ),
            reactions = (Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = 300e3),),
            lateral_shrinkage_law = lateral,
        )
    else
        Material(
            name = :audit_pmma,
            components = (
                SolidComponent(
                    :pmma,
                    ρ = 1190.0,
                    c = 1420.0,
                    k = 0.21,
                    ε = 0.86,
                    λ = something(solid_λ, 0.0),
                ),
                GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5),
            ),
            reactions = (Reaction(:pyrolysis, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3),),
        )
    end
    thermal = (
        top = RadiativeFluxBC(flux = flux) +
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
    )
    mass =
        wood ?
        (
            top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
            bottom = ImpermeableBC(),
        ) : (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC())
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, n, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(; thermal, mass),
        tspan = (0.0, t_end),
        dt_initial = 1e-3,
        output_times = collect(range(0.0, t_end, length = outputs)),
    )
    return problem,
    (;
        use_ale = kind in (:ale, :ale_darcy),
        radiation_model = wood ? BEER_LAMBERT : SURFACE_ABSORPTION,
        progress = false,
    )
end

const DECKS = (
    (:eulerian, 25),
    (:eulerian, 100),
    (:ale, 25),
    (:ale, 100),
    (:ale_darcy, 25),
    (:ale_darcy, 100),
)
# (label, step-update variant, `integrator` argument of `solve`): a constructor
# or constructor-like function, so `_resolve_integrator` attaches the package's
# default linear solver.
kencarp4_raw(; kw...) = KenCarp4(; smooth_est = false, kw...)
const CONFIGS = (
    ("KenCarp4", :current, KenCarp4),
    ("KenCarp4_raw", :current, kencarp4_raw),
    ("FBDF", :current, FBDF),
    ("FBDF", :observer, FBDF),
    ("Rodas5P", :observer, Rodas5P),
)

"""
    guarded_deck_solve(kind, n; integrator, abstol, reltol) -> NamedTuple

One solve of the deck under the wall-clock deadline; never throws.
"""
function guarded_deck_solve(
    kind,
    n;
    integrator,
    abstol,
    reltol,
    timeout_s = SOLVE_TIMEOUT_S,
)
    prob, kw = problem_for(kind, n)
    DEADLINE[] = time() + timeout_s
    TIMED_OUT[] = false
    LAST_T[] = NaN
    GC.gc()
    timed = @timed try
        Logging.with_logger(QuietLogger()) do
            Base.invokelatest(
                solve,
                deepcopy(prob);
                kw...,
                integrator,
                abstol,
                reltol,
                maxiters = MAXITERS,
            )
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
            t_reached = LAST_T[],
        )
    end
    return (;
        sol,
        retcode = sol.retcode,
        timed_out = TIMED_OUT[],
        message = "",
        wall_s = timed.time,
        t_reached = sol.retcode == :Success ? sol.t[end] : LAST_T[],
    )
end

finished(r) =
    r.sol !== nothing && r.retcode == :Success && !r.timed_out && r.t_reached == T_END

# One warm-up solve followed by `repeats` timed solves (skipped when the
# warm-up did not finish). Returns the last solution and the timings.
function timed_deck_runs(kind, n; repeats = 3, kwargs...)
    warm = guarded_deck_solve(kind, n; kwargs...)
    if !finished(warm)
        return (;
            sol = warm.sol,
            retcode = warm.retcode,
            timed_out = warm.timed_out,
            message = warm.message,
            t_reached = warm.t_reached,
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
        r = guarded_deck_solve(kind, n; kwargs...)
        if !finished(r)
            return (;
                sol = r.sol,
                retcode = r.retcode,
                timed_out = r.timed_out,
                message = "repeat: " * r.message,
                t_reached = r.t_reached,
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
        t_reached = last.t_reached,
        first_s = warm.wall_s,
        warm_median_s = median(times),
        warm_min_s = minimum(times),
        warm_max_s = maximum(times),
        repeats_consistent = consistent,
    )
end

_maxabs(a, b) =
    maximum((abs(x - y) for (x, y) in zip(a, b) if isfinite(x) && isfinite(y)); init = 0.0)

const HEADER = (
    "deck",
    "n_cells",
    "integrator",
    "callback_variant",
    "reltol",
    "abstol",
    "retcode",
    "timed_out",
    "message",
    "t_reached",
    "steps",
    "rejects",
    "nf",
    "njacs",
    "nw",
    "nsolve",
    "nonliniter",
    "nonlinfail",
    "first_s",
    "warm_median_s",
    "warm_min_s",
    "warm_max_s",
    "repeats_consistent",
    "ref_retcode",
    "ref_steps",
    "ref_wall_s",
    "err_T_final_K",
    "err_T_all_outputs_K",
    "err_z_final_m",
    "ref_cross_check_K",
)
const ROWS = Dict{String,Any}[]

function record!(io, vals)
    csv_row(io, vals...)
    push!(ROWS, Dict(zip(HEADER, vals)))
end

function run_matrix(io)
    for (kind, n) in DECKS
        # Reference: KenCarp4 with the flag set at the tight tolerance.
        STEP_VARIANT[] = :current
        say("REFERENCE ", (kind, n))
        ref = guarded_deck_solve(kind, n; integrator = KenCarp4, REF_TOL...)
        rst = solve_stats(ref.sol)
        ref_ok = finished(ref)
        say("   -> ", (ref.retcode, rst.steps, rst.nonlinfail, ref.wall_s))
        # Second tight solve with a different method family (FBDF as an
        # observer) at the same tolerances: the final-time difference between
        # the two bounds the reference's own error.
        STEP_VARIANT[] = :observer
        say("REFERENCE_CHECK ", (kind, n))
        ref2 = guarded_deck_solve(kind, n; integrator = FBDF, REF_TOL...)
        rst2 = solve_stats(ref2.sol)
        xcheck =
            ref_ok && finished(ref2) ?
            _maxabs(
                view(ref.sol.T, :, size(ref.sol.T, 2)),
                view(ref2.sol.T, :, size(ref2.sol.T, 2)),
            ) : NaN
        say(
            "   -> ",
            (ref2.retcode, rst2.steps, rst2.nonlinfail, ref2.wall_s),
            " cross-check K=",
            xcheck,
        )
        for (lab, var, rr, rs, e) in (
            ("KenCarp4", "reference", ref, rst, 0.0),
            ("FBDF", "reference_observer", ref2, rst2, xcheck),
        )
            record!(
                io,
                (
                    string(kind),
                    n,
                    lab,
                    var,
                    REF_TOL.reltol,
                    REF_TOL.abstol,
                    string(rr.retcode),
                    rr.timed_out,
                    rr.message,
                    rr.t_reached,
                    rs.steps,
                    rs.rejects,
                    rs.nf,
                    rs.njacs,
                    rs.nw,
                    rs.nsolve,
                    rs.nonliniter,
                    rs.nonlinfail,
                    rr.wall_s,
                    NaN,
                    NaN,
                    NaN,
                    false,
                    string(ref.retcode),
                    rst.steps,
                    ref.wall_s,
                    e,
                    e,
                    NaN,
                    xcheck,
                ),
            )
        end
        for (label, variant, alg) in CONFIGS
            STEP_VARIANT[] = variant
            say("RUN ", (kind, n, label, variant))
            r = timed_deck_runs(kind, n; integrator = alg, TOL...)
            st = solve_stats(r.sol)
            ok = finished(r) && ref_ok
            eT =
                ok ?
                _maxabs(
                    view(r.sol.T, :, size(r.sol.T, 2)),
                    view(ref.sol.T, :, size(ref.sol.T, 2)),
                ) : NaN
            eTall = ok ? _maxabs(r.sol.T, ref.sol.T) : NaN
            ez =
                ok && !isempty(ref.sol.z) ?
                _maxabs(
                    view(r.sol.z, :, size(r.sol.z, 2)),
                    view(ref.sol.z, :, size(ref.sol.z, 2)),
                ) : NaN
            record!(
                io,
                (
                    string(kind),
                    n,
                    label,
                    string(variant),
                    TOL.reltol,
                    TOL.abstol,
                    string(r.retcode),
                    r.timed_out,
                    r.message,
                    r.t_reached,
                    st.steps,
                    st.rejects,
                    st.nf,
                    st.njacs,
                    st.nw,
                    st.nsolve,
                    st.nonliniter,
                    st.nonlinfail,
                    r.first_s,
                    r.warm_median_s,
                    r.warm_min_s,
                    r.warm_max_s,
                    r.repeats_consistent,
                    string(ref.retcode),
                    rst.steps,
                    ref.wall_s,
                    eT,
                    eTall,
                    ez,
                    xcheck,
                ),
            )
            say(
                "   -> ",
                (r.retcode, st.steps, st.rejects, st.njacs, st.nonlinfail, r.warm_median_s),
                " errT_final=",
                eT,
            )
        end
    end
end

t_total = @elapsed open(joinpath(RESULTS_DIR, "integrator_observer_audit.csv"), "w") do io
    csv_row(io, HEADER...)
    Base.invokelatest(run_matrix, io)
end
say("TOTAL wall time (s): ", t_total)

# ----------------------------------------------------------------------------
# Markdown: per-deck tables, rankings by time and by error, FBDF's position
# ----------------------------------------------------------------------------
fmtv(x) =
    x isa AbstractString ? x :
    x isa Bool ? string(x) :
    x isa Integer ? string(x) :
    !isfinite(x) ? string(x) :
    (abs(x) >= 1e4 || (0 < abs(x) < 1e-3) ? @sprintf("%.3e", x) : @sprintf("%.4g", x))
config_label(r) =
    r["integrator"] * (r["callback_variant"] == "observer" ? " (observer)" : "")
const CONFIG_LABELS = [l * (v == :observer ? " (observer)" : "") for (l, v, _) in CONFIGS]
deck_rows(kind, n) = [
    r for r in ROWS if r["deck"] == string(kind) &&
        r["n_cells"] == n &&
        !startswith(r["callback_variant"], "reference")
]
ref2_row(kind, n) = ROWS[findfirst(
    r ->
        r["deck"] == string(kind) &&
        r["n_cells"] == n &&
        r["callback_variant"] == "reference_observer",
    ROWS,
)]
ref_row(kind, n) = ROWS[findfirst(
    r ->
        r["deck"] == string(kind) &&
        r["n_cells"] == n &&
        r["callback_variant"] == "reference",
    ROWS,
)]
row_ok(r) = r["retcode"] == "Success" && !r["timed_out"] && isfinite(r["warm_median_s"])

# Ranks (1 = best) of the finished configurations of one deck by `key`;
# unfinished ones get the last rank.
function ranks(rows, key)
    fin = filter(row_ok, rows)
    order = sortperm([r[key] for r in fin])
    rk = Dict{String,Int}(config_label(fin[order[i]]) => i for i in eachindex(order))
    for r in rows
        haskey(rk, config_label(r)) || (rk[config_label(r)] = length(rows))
    end
    return rk
end

# FBDF rows of the audit's ranking (results/newton_baseline.md, table (b)),
# for the comparison; absent when the file is not there.
function audit_fbdf_rows()
    path = joinpath(RESULTS_DIR, "newton_baseline.md")
    isfile(path) || return Dict{Tuple{String,Int},NamedTuple}()
    out = Dict{Tuple{String,Int},NamedTuple}()
    for line in eachline(path)
        startswith(line, "| ") || continue
        f = strip.(split(line, '|'))
        length(f) >= 14 || continue
        f[4] == "FBDF" || continue
        n = tryparse(Int, f[3])
        n === nothing && continue
        out[(f[2], n)] = (;
            steps = tryparse(Int, f[5]),
            median_s = tryparse(Float64, f[11]),
            max_T_err = tryparse(Float64, f[12]),
        )
    end
    return out
end

open(joinpath(RESULTS_DIR, "integrator_observer_audit.md"), "w") do md
    println(md, "# FBDF observer audit on the performance-audit decks\n")
    println(
        md,
        "Generated by `scripts/mesh_phase0/baselines/integrator_observer_audit.jl`; data in `integrator_observer_audit.csv`; log alongside. ",
        "Run as `julia --project=test --startup-file=no --threads=1 <script>` from the repository root.\n",
    )
    println(
        md,
        "Why: the audit's integrator ranking (`newton_baseline.md`, table (b), from the gitignored `docs/performance_audit/solves.csv`) ran every ",
        "integrator with the package's step-update callback, which flags a derivative discontinuity after every accepted step; a multistep ",
        "method restarts at order one on every such flag, so the audit's FBDF rows were suspected of measuring a crippled configuration. Here the step-update ",
        "callback is replaced in-process (same technique as `integrator_flag.jl`): `current` = mirror + `derivative_discontinuity!(integrator, true)` ",
        "(the package), `observer` = mirror + `derivative_discontinuity!(integrator, false)`. Decks: `problem_for(kind, n)` copied from ",
        "`scripts/performance_audit/audit.jl` (eulerian and ale: PMMA, surface absorption; ale_darcy: wood with Darcy flow and Beer–Lambert) at ",
        "25 and 100 cells, 80 s, 3 outputs, the deck's own solve kwargs; Jacobian strategy and linear solver are the package defaults ",
        "(Eulerian: `Structured(strategy = DirectSolve())` + KLU; ALE: `Structured(strategy = ApproxSolve())` + KLU), every integrator passed as a ",
        "constructor so `_resolve_integrator` attaches that linear solver and `autodiff = AutoFiniteDiff()` (unused: the package supplies `jac!`). ",
        "`KenCarp4_raw` = `KenCarp4(smooth_est = false)`. Tolerances abstol $(TOL.abstol), reltol $(TOL.reltol) (the `solve` defaults, as in the audit). ",
        "Reference per deck: KenCarp4 with the flag set at abstol $(REF_TOL.abstol), reltol $(REF_TOL.reltol) (the audit's convention; one solve); a second tight solve with FBDF as an observer at the same tolerances bounds the reference's own error (`cross-check`, the final-time difference between the two; candidate errors below it are not resolved). ",
        "`err T final` = maximum absolute temperature difference over cells at t = $(Int(T_END)) s against the reference; `err T all` = the same over the three outputs; ",
        "`err z final` = maximum node-position difference at the final time (ALE only). Warm time = median of 3 solves after one warm-up; ",
        "`maxiters = $(MAXITERS)`, wall-clock guard $(SOLVE_TIMEOUT_S) s per solve. Wall times are indicative (shared machine, another Julia process running concurrently).\n",
    )
    for (kind, n) in DECKS
        rr = ref_row(kind, n)
        r2 = ref2_row(kind, n)
        println(
            md,
            "### $kind, $n cells (reference KenCarp4: $(rr["retcode"]), $(rr["steps"]) steps, $(rr["nonlinfail"]) Newton failures, $(fmtv(rr["ref_wall_s"])) s; ",
            "cross-check FBDF observer at the same tolerances: $(r2["retcode"]), $(r2["steps"]) steps, $(fmtv(r2["first_s"])) s; ",
            "the two differ by $(fmtv(rr["ref_cross_check_K"])) K at the final time)\n",
        )
        println(
            md,
            "| integrator | callback | retcode (t reached) | steps | rej | nf | njacs | nonliniter | its/step | nlfail | first s | warm s | consistent | err T final K | err T all K | err z final m |",
        )
        println(md, "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
        for r in deck_rows(kind, n)
            println(
                md,
                "| ",
                join(
                    (
                        r["integrator"],
                        r["callback_variant"],
                        r["retcode"] *
                        (r["timed_out"] ? " (timed out)" : "") *
                        (r["retcode"] == "Success" ? "" : " (t=$(fmtv(r["t_reached"])))") *
                        (isempty(r["message"]) ? "" : " " * first(r["message"], 60)),
                        r["steps"],
                        r["rejects"],
                        r["nf"],
                        r["njacs"],
                        r["nonliniter"],
                        r["steps"] > 0 ? @sprintf("%.1f", r["nonliniter"] / r["steps"]) :
                        "–",
                        r["nonlinfail"],
                        fmtv(r["first_s"]),
                        fmtv(r["warm_median_s"]),
                        r["repeats_consistent"],
                        fmtv(r["err_T_final_K"]),
                        fmtv(r["err_T_all_outputs_K"]),
                        fmtv(r["err_z_final_m"]),
                    ),
                    " | ",
                ),
                " |",
            )
        end
        println(md)
    end

    # Rankings.
    time_rk = Dict{String,Vector{Int}}(l => Int[] for l in CONFIG_LABELS)
    err_rk = Dict{String,Vector{Int}}(l => Int[] for l in CONFIG_LABELS)
    for (kind, n) in DECKS
        rows = deck_rows(kind, n)
        tr = ranks(rows, "warm_median_s")
        er = ranks(rows, "err_T_final_K")
        for l in CONFIG_LABELS
            push!(time_rk[l], tr[l])
            push!(err_rk[l], er[l])
        end
    end
    println(md, "## Ranking over the six decks\n")
    println(
        md,
        "Rank 1 = fastest (by warm median time) or most accurate (by `err T final`); a configuration that did not finish gets the last rank on that deck. Ranks per deck in the order ",
        join(("$kind/$n" for (kind, n) in DECKS), ", "),
        ". Reference cross-check (K at the final time, same deck order): ",
        join((fmtv(ref_row(kind, n)["ref_cross_check_K"]) for (kind, n) in DECKS), ", "),
        "; error ranks between configurations whose errors differ by less than the deck's cross-check are not resolved.\n",
    )
    println(
        md,
        "| configuration | time: mean rank | first on | ranks | error: mean rank | first on | ranks |",
    )
    println(md, "|---|---|---|---|---|---|---|")
    for l in sort(CONFIG_LABELS, by = l -> mean(time_rk[l]))
        println(
            md,
            @sprintf(
                "| %s | %.2f | %d/%d | %s | %.2f | %d/%d | %s |",
                l,
                mean(time_rk[l]),
                count(==(1), time_rk[l]),
                length(DECKS),
                join(time_rk[l], " "),
                mean(err_rk[l]),
                count(==(1), err_rk[l]),
                length(DECKS),
                join(err_rk[l], " ")
            )
        )
    end

    # FBDF's position.
    println(md, "\n## FBDF: flag set versus observer\n")
    audit = audit_fbdf_rows()
    println(
        md,
        "| deck | cells | FBDF current: steps / warm s / err T final K | FBDF observer: steps / warm s / err T final K | steps× / time× / error× (observer / current) | audit's FBDF row: steps / median s / max T err K |",
    )
    println(md, "|---|---|---|---|---|---|")
    for (kind, n) in DECKS
        rows = deck_rows(kind, n)
        a = rows[findfirst(
            r -> r["integrator"] == "FBDF" && r["callback_variant"] == "current",
            rows,
        )]
        b = rows[findfirst(
            r -> r["integrator"] == "FBDF" && r["callback_variant"] == "observer",
            rows,
        )]
        au = get(audit, (string(kind), n), nothing)
        println(
            md,
            @sprintf(
                "| %s | %d | %s / %s / %s | %s / %s / %s | %s / %s / %s | %s |",
                kind,
                n,
                a["retcode"] == "Success" ? string(a["steps"]) :
                "$(a["retcode"]) after $(a["steps"]) steps at t = $(fmtv(a["t_reached"])) s",
                fmtv(a["warm_median_s"]),
                fmtv(a["err_T_final_K"]),
                b["steps"],
                fmtv(b["warm_median_s"]),
                fmtv(b["err_T_final_K"]),
                a["steps"] > 0 ? @sprintf("%.2f", b["steps"] / a["steps"]) : "–",
                fmtv(b["warm_median_s"] / a["warm_median_s"]),
                fmtv(b["err_T_final_K"] / a["err_T_final_K"]),
                au === nothing ? "–" :
                @sprintf("%d / %.3f / %s", au.steps, au.median_s, fmtv(au.max_T_err))
            )
        )
    end
    fc = "FBDF"
    fo = "FBDF (observer)"
    kc = "KenCarp4"
    # Observer step counts against the audit's FBDF rows.
    agree = String[]
    for (kind, n) in DECKS
        au = get(audit, (string(kind), n), nothing)
        au === nothing && continue
        b = deck_rows(kind, n)[findfirst(
            r -> r["integrator"] == "FBDF" && r["callback_variant"] == "observer",
            deck_rows(kind, n),
        )]
        push!(
            agree,
            @sprintf(
                "%s/%d %d vs %d (%+.1f%%)",
                kind,
                n,
                b["steps"],
                au.steps,
                100 * (b["steps"] / au.steps - 1)
            )
        )
    end
    isempty(agree) || println(
        md,
        "\nFBDF as an observer reproduces the audit's FBDF accepted-step counts (observer vs audit): ",
        join(agree, ", "),
        ". With the flag set FBDF does not finish any deck here (`MaxIters` at $(MAXITERS) steps, one Jacobian, one Newton iteration per step, reaching only t = 3–7 s of 80 s), so the audit's FBDF rows cannot have run under an effective per-step restart: they correspond to the observer configuration. ",
        "`newton_baseline.md` reports FBDF at mean rank 1.17 by time (first on 5 of 6 decks) and 2.17 by error.",
    )
    println(
        md,
        @sprintf(
            "\nBy time FBDF with the flag set ranks %.2f on average (first on %d of %d decks) and as an observer %.2f (first on %d); by final-time error %.2f and %.2f. KenCarp4 (flag set, the package default) ranks %.2f by time and %.2f by error. ",
            mean(time_rk[fc]),
            count(==(1), time_rk[fc]),
            length(DECKS),
            mean(time_rk[fo]),
            count(==(1), time_rk[fo]),
            mean(err_rk[fc]),
            mean(err_rk[fo]),
            mean(time_rk[kc]),
            mean(err_rk[kc])
        ),
        mean(time_rk[fo]) < mean(time_rk[fc]) ?
        "FBDF's position by time improves once it runs as an observer" :
        "FBDF's position by time does not improve as an observer",
        mean(err_rk[fo]) < mean(err_rk[fc]) ? "; its position by error improves as well." :
        "; its position by error does not improve.",
    )
end
say("wrote ", joinpath(RESULTS_DIR, "integrator_observer_audit.md"))
