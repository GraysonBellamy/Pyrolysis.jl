# Park every qualifying cell in one depletion-callback invocation.
#
# The package's depletion callback parks at most one cell per invocation; a
# cascade in which every cell crosses the threshold at once (V-BO-04) pays one
# park event per accepted step. This driver replaces the callback in-process
# with a variant whose affect! keeps calling `deplete!` while the callback's
# own event predicate (`has_event`: a depletion candidate exists, or on a fixed
# mesh an interior void) still holds after a park, then flags the
# discontinuity and proposes DT_AFTER_PARK once. Re-testing the predicate
# before every further `deplete!` keeps the package's termination semantics:
# `deplete!` checks `n_live <= min_cells` before it scans for a candidate, so
# calling it once more after a cascade that leaves exactly `min_cells` live
# cells would terminate the solve at the cascade instead of at the next
# threshold crossing. The burnout cases are measured under three variants:
#
#   current           the package method
#   current_override  the replacement with the loop disabled (must reproduce
#                     `current` step for step)
#   one_callback      the replacement with the loop enabled
#
# Cases: V-BO-04 (40 cells, min_cells = 2, the build default), V-BO-04 with
# min_cells = 20 (the first cascade leaves exactly min_cells live cells; the
# solve must terminate at the next crossing, retcode Terminated), V-BO-02 and
# V-BO-01. A final check re-solves V-BO-01 at three output spacings to show
# how much of the mass-ledger closure is output-grid quadrature.
#
# Outputs (scripts/mesh_phase0/results/): one_callback_parking.csv (one row
# per case × variant), one_callback_parking_qois.csv (one row per QoI),
# one_callback_parking_ledger.csv (output-spacing check) and
# one_callback_parking.md.

include(joinpath(@__DIR__, "common.jl"))

include(joinpath(ROOT, "test", "verification", "cases", "V-BO-01.jl"))
include(joinpath(ROOT, "test", "verification", "cases", "V-BO-02.jl"))
include(joinpath(ROOT, "test", "verification", "cases", "V-BO-04.jl"))

using Printf

const MAXITERS = 200_000
const PARK_LOOP = Ref(false)
const PARK_INVOCATIONS = Ref(0)     # affect! calls that parked at least one cell
const LAST_PARK_T = Ref(NaN)        # integrator time of the last park (any variant)
const PREV_T = Ref(NaN)             # accepted-step time before LAST_T (last step size)

# Step-update callback: the package behaviour plus the wall-clock guard.
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
    deadline = Main.DEADLINE
    timed_out = Main.TIMED_OUT
    last_t = Main.LAST_T
    prev_t = Main.PREV_T
    function affect!(integrator)
        sync_mesh!(integrator.u)
        prev_t[] = last_t[]
        last_t[] = integrator.t
        if time() > deadline[]
            timed_out[] = true
            terminate!(integrator)
        end
        derivative_discontinuity!(integrator, true)
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

# Installed after the `current` rows have been measured with the package method.
function install_depletion_override()
    @eval Pyrolysis.Solver function create_depletion_callback(
        ws,
        mesh::Mesh1D,
        material,
        cb_state::DepletionCallbackState;
        threshold::Float64 = 0.05,
        min_cells::Int = 2,
    )
        layout = cb_state.layout
        fixed_mesh = !(mesh_mode(layout.mode) isa ALE)

        has_event = FunctionWrapper{Bool,Tuple{_StateVector}}(
            u ->
                depletion_candidate(u, layout, mesh, material, cb_state, threshold) != 0 ||
                fixed_mesh && interior_void_cell(u, layout, material, cb_state) != 0,
        )

        deplete! = FunctionWrapper{DepletionEvent,Tuple{_StateVector,Float64}}(
            (u, t) -> begin
                if fixed_mesh
                    i_void = interior_void_cell(u, layout, material, cb_state)
                    if i_void != 0
                        # string(): a `$` inside this @eval body would be interpolated at eval time.
                        @warn string(
                            "Cell ",
                            i_void,
                            ", below the top live cell ",
                            cb_state.n_live,
                            ", is a void (matrix volume fraction < ",
                            MIXING_BLEND_LO,
                            ") at t=",
                            round(t, digits = 4),
                            " s: the material has gasified ",
                            "from within. The model carries no heat across a void, so the ",
                            "fixed-mesh solve stops here. A regressing material belongs on ",
                            "the ALE mesh (use_ale = true).",
                        )
                        return STOP_SOLVE
                    end
                end
                if cb_state.n_live <= min_cells
                    @warn string("At minimum cells (", min_cells, "), terminating simulation")
                    return STOP_SOLVE
                end
                i_park = depletion_candidate(
                    u,
                    layout,
                    mesh,
                    material,
                    cb_state,
                    threshold,
                )
                i_park == 0 && return NO_PARK
                n_live = park!(u, layout, mesh, material, cb_state, i_park, t)
                ws.top_face_ids[1] = mesh.boundary_faces[1][1]
                Main.LAST_PARK_T[] = t
                @info string(
                    "Parked cell ",
                    i_park,
                    " at t=",
                    round(t, digits = 4),
                    " s; ",
                    n_live,
                    " live cells",
                )
                return CELL_PARKED
            end,
        )

        condition = (u, t, integrator) -> has_event(u)

        loop = Main.PARK_LOOP[]
        invocations = Main.PARK_INVOCATIONS
        affect! = (integrator) -> begin
            u = integrator.u::_StateVector
            event = deplete!(u, integrator.t)
            parked = event == CELL_PARKED
            if loop
                # The new top cell may qualify as well: keep parking while
                # the callback's own event predicate still holds. The
                # predicate is re-tested before every further `deplete!`
                # because `deplete!` terminates at `n_live <= min_cells`
                # before it looks for a candidate; the fixed-mesh void
                # test is part of the predicate and `deplete!` turns it
                # into STOP_SOLVE as in the package method.
                while event == CELL_PARKED && has_event(u)
                    event = deplete!(u, integrator.t)
                end
            end
            parked && (invocations[] += 1)
            if event == STOP_SOLVE
                terminate!(integrator)
            elseif parked
                derivative_discontinuity!(integrator, true)
                set_proposed_dt!(integrator, DT_AFTER_PARK)
            end
            return nothing
        end

        return DiscreteCallback(condition, affect!; save_positions = (false, false))
    end
end

# A copy of `case` that solves with `min_cells = m` and expects the callback
# to terminate the solve.
function with_min_cells(case, m; id)
    return VerificationCase(
        id = id,
        title = case.title,
        tags = case.tags,
        references = case.references,
        build = case.build,
        solve_kwargs = merge(case.solve_kwargs, (; min_cells = m)),
        solve_fn = case.solve_fn,
        expected_retcode = :Terminated,
        qois = case.qois,
        ci = case.ci,
        report = case.report,
    )
end

const VBO04 = _verification_case_vbo04()
const CASES = (
    ("V-BO-04", VBO04),
    ("V-BO-04/min20", with_min_cells(VBO04, 20; id = "V-BO-04/min20")),
    ("V-BO-02", _verification_case_vbo02()),
    ("V-BO-01", _verification_case_vbo01()),
)
min_cells_of(case) = get(case.solve_kwargs, :min_cells, 2)

const HEADER = (
    "case",
    "n_cells",
    "min_cells",
    "variant",
    "retcode",
    "timed_out",
    "message",
    "t_last_accepted",
    "last_step_s",
    "last_park_t",
    "park_invocations",
    "cells_parked",
    "live_cells_end",
    "steps",
    "rejects",
    "nf",
    "njacs",
    "nonliniter",
    "nonlinfail",
    "first_s",
    "warm_median_s",
    "warm_min_s",
    "warm_max_s",
    "repeats_consistent",
    "mass_ledger",
    "qoi_count",
    "qoi_failures",
    "qoi_max_error_over_tol",
)
const QOI_HEADER = (
    "case",
    "n_cells",
    "min_cells",
    "variant",
    "qoi",
    "value",
    "exact",
    "error",
    "tol",
    "pass",
)
const LEDGER_HEADER = (
    "case",
    "n_cells",
    "variant",
    "output_spacing_s",
    "n_outputs",
    "mass_ledger",
    "mass_drop",
    "final_mass",
    "steps",
    "cells_parked",
)

const RUN_ROWS = Dict{String,Any}[]      # one per case × variant
const QOI_ROWS = Dict{String,Any}[]      # one per QoI
const LEDGER_ROWS = Dict{String,Any}[]   # output-spacing check

function run_variant(io, qio, variant::String)
    for (case_id, case) in CASES
        n = case.ci.n_cells
        mc = min_cells_of(case)
        say("RUN ", (case_id, n, mc, variant))
        r = timed_runs(case, n; maxiters = MAXITERS)
        st = solve_stats(r.sol)
        t_last = LAST_T[]
        last_step = LAST_T[] - PREV_T[]
        if r.sol === nothing
            vals = (
                case_id,
                n,
                mc,
                variant,
                string(r.retcode),
                r.timed_out,
                r.message,
                t_last,
                last_step,
                LAST_PARK_T[],
                -1,
                -1,
                -1,
                st.steps,
                st.rejects,
                st.nf,
                st.njacs,
                st.nonliniter,
                st.nonlinfail,
                r.first_s,
                r.warm_median_s,
                r.warm_min_s,
                r.warm_max_s,
                r.repeats_consistent,
                NaN,
                0,
                -1,
                NaN,
            )
            csv_row(io, vals...)
            push!(RUN_ROWS, Dict(zip(HEADER, vals)))
            continue
        end
        sol = r.sol
        cells_parked = length(sol.extras.park_times)
        last_park = isempty(sol.extras.park_times) ? NaN : last(sol.extras.park_times)
        # The package method parks at most one cell per invocation, so its
        # invocation count equals its park count.
        invocations = variant == "current" ? cells_parked : PARK_INVOCATIONS[]
        ledger = Base.invokelatest(_vbo_measures, sol).ledger
        qois = Base.invokelatest(case.qois, sol)
        for q in qois
            qv = (
                case_id,
                n,
                mc,
                variant,
                q.name,
                q.value,
                q.exact,
                q.error,
                q.tol,
                q.error < q.tol,
            )
            csv_row(qio, qv...)
            push!(QOI_ROWS, Dict(zip(QOI_HEADER, qv)))
        end
        vals = (
            case_id,
            n,
            mc,
            variant,
            string(r.retcode),
            r.timed_out,
            r.message,
            t_last,
            last_step,
            last_park,
            invocations,
            cells_parked,
            sol.extras.live_cells[end],
            st.steps,
            st.rejects,
            st.nf,
            st.njacs,
            st.nonliniter,
            st.nonlinfail,
            r.first_s,
            r.warm_median_s,
            r.warm_min_s,
            r.warm_max_s,
            r.repeats_consistent,
            ledger,
            length(qois),
            count(q -> !(q.error < q.tol), qois),
            maximum(q.error / q.tol for q in qois),
        )
        csv_row(io, vals...)
        push!(RUN_ROWS, Dict(zip(HEADER, vals)))
        say(
            "   -> ",
            (
                r.retcode,
                t_last,
                last_park,
                invocations,
                cells_parked,
                sol.extras.live_cells[end],
                st.steps,
                st.rejects,
                st.nonlinfail,
                r.warm_median_s,
                ledger,
            ),
        )
    end
end

# V-BO-01 at several output spacings: the ledger is a trapezoid of the
# sampled MLR on the output grid, so its closure depends on that grid.
function ledger_check(lio, variant::String)
    case = CASES[end][2]
    n = case.ci.n_cells
    for spacing in (0.5, 0.125, 0.03125)
        say("LEDGER ", (variant, spacing))
        problem = case.build(n)
        empty!(problem.output_times)
        append!(problem.output_times, 0.0:spacing:500.0)
        DEADLINE[] = time() + SOLVE_TIMEOUT_S
        TIMED_OUT[] = false
        PRE_SOLVE_HOOK[]()
        sol = Logging.with_logger(QuietLogger()) do
            Base.invokelatest(
                _case_solve,
                case,
                problem;
                case.solve_kwargs...,
                maxiters = MAXITERS,
            )
        end
        m = sol.extras.total_mass
        vals = (
            "V-BO-01",
            n,
            variant,
            spacing,
            length(problem.output_times),
            Base.invokelatest(_vbo_measures, sol).ledger,
            m[1] - m[end],
            m[end],
            sol.raw_solution.stats.naccept,
            length(sol.extras.park_times),
        )
        csv_row(lio, vals...)
        push!(LEDGER_ROWS, Dict(zip(LEDGER_HEADER, vals)))
        say("   -> ", vals[6:end])
    end
end

# `timed_runs` makes four solves per configuration; the counters are reset
# before each of them so the row reports the last solve only.
PRE_SOLVE_HOOK[] = () -> (PARK_INVOCATIONS[] = 0; LAST_PARK_T[] = NaN; PREV_T[] = NaN)

t_total = @elapsed begin
    io = open(joinpath(RESULTS_DIR, "one_callback_parking.csv"), "w")
    qio = open(joinpath(RESULTS_DIR, "one_callback_parking_qois.csv"), "w")
    lio = open(joinpath(RESULTS_DIR, "one_callback_parking_ledger.csv"), "w")
    csv_row(io, HEADER...)
    csv_row(qio, QOI_HEADER...)
    csv_row(lio, LEDGER_HEADER...)
    Base.invokelatest(run_variant, io, qio, "current")
    install_depletion_override()
    PARK_LOOP[] = false
    Base.invokelatest(run_variant, io, qio, "current_override")
    Base.invokelatest(ledger_check, lio, "current_override")
    PARK_LOOP[] = true
    Base.invokelatest(run_variant, io, qio, "one_callback")
    Base.invokelatest(ledger_check, lio, "one_callback")
    close(io)
    close(qio)
    close(lio)
end
say("TOTAL wall time (s): ", t_total)

# ----------------------------------------------------------------------------
# Markdown
# ----------------------------------------------------------------------------
fmtv(x) =
    x isa AbstractString ? x :
    x isa Bool ? string(x) :
    x isa Integer ? string(x) :
    !isfinite(x) ? string(x) :
    (abs(x) >= 1e4 || (0 < abs(x) < 1e-3) ? @sprintf("%.3e", x) : @sprintf("%.6g", x))
row_of(case_id, variant) =
    RUN_ROWS[findfirst(r -> r["case"] == case_id && r["variant"] == variant, RUN_ROWS)]
qoi_of(case_id, variant, name) = QOI_ROWS[findfirst(
    r -> r["case"] == case_id && r["variant"] == variant && r["qoi"] == name,
    QOI_ROWS,
)]
const VARIANT_NAMES = ("current", "current_override", "one_callback")

open(joinpath(RESULTS_DIR, "one_callback_parking.md"), "w") do md
    println(md, "# One-callback parking\n")
    println(
        md,
        "Generated by `scripts/mesh_phase0/baselines/one_callback_parking.jl` (shared helpers in `common.jl`); data in ",
        "`one_callback_parking.csv`, `one_callback_parking_qois.csv` and `one_callback_parking_ledger.csv`; log alongside. ",
        "Run as `julia --project=test --startup-file=no --threads=1 <script>` from the repository root.\n",
    )
    println(
        md,
        "Variants of the depletion callback (`create_depletion_callback`, replaced in-process): `current` = package method ",
        "(at most one cell per invocation); `current_override` = the replacement with its loop disabled (must reproduce `current`); ",
        "`one_callback` = the replacement looping `deplete!` while it returns `CELL_PARKED` **and** the callback's own event predicate ",
        "`has_event(u)` (a depletion candidate exists, or on a fixed mesh an interior void) still holds, then one ",
        "`derivative_discontinuity!(integrator, true)` and one `set_proposed_dt!(integrator, DT_AFTER_PARK)`; `STOP_SOLVE` still terminates. ",
        "The predicate is re-tested before every further `deplete!` because `deplete!` returns `STOP_SOLVE` at `n_live <= min_cells` ",
        "before it scans for a candidate: a loop that called it once more after a cascade leaving exactly `min_cells` live cells would ",
        "terminate at the cascade instead of at the next threshold crossing. ",
        "Each case at its CI cell count (V-BO-04: 40, the build default) through its own `solve_fn`/`solve_kwargs`; ",
        "`V-BO-04/min20` is V-BO-04 with `min_cells = 20` (expected retcode `Terminated`: the first cascade leaves 20 live cells and the ",
        "callback must stop the solve at the next crossing). `maxiters = $(MAXITERS)`; wall-clock guard $(SOLVE_TIMEOUT_S) s. ",
        "`t_last` = integrator time of the last accepted step (the termination time for `Terminated` runs); `last park` = time of the last park. ",
        "Park invocations of `current` equal its park count by construction; the other two are counted in the override. ",
        "Warm time = median of 3 solves after one warm-up; indicative (shared machine, another Julia process running concurrently). ",
        "Mass ledger = |Δm − ∫ṁ″dt| / Δm with the integral a trapezoid of the sampled MLR on the output times (`_vbo_measures`), ",
        "the same quantity as the \"mass ledger closure\" QoI of V-BO-01/02; it is a sampled diagnostic on the output grid, not a ",
        "conservation statement (see the output-spacing check below).\n",
    )
    println(md, "## Runs\n")
    println(
        md,
        "| case | cells | min_cells | variant | retcode | t_last s | last step s | last park s | park invocations | cells parked | live cells end | steps | rej | nf | njacs | nlfail | warm s | mass ledger | QoIs | QoI fails | max error/tol |",
    )
    println(
        md,
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    )
    for r in RUN_ROWS
        println(
            md,
            "| ",
            join(
                (
                    r["case"],
                    r["n_cells"],
                    r["min_cells"],
                    r["variant"],
                    r["retcode"] *
                    (r["timed_out"] ? " (timed out)" : "") *
                    (isempty(r["message"]) ? "" : " " * first(r["message"], 60)),
                    @sprintf("%.6f", r["t_last_accepted"]),
                    fmtv(r["last_step_s"]),
                    @sprintf("%.6f", r["last_park_t"]),
                    r["park_invocations"],
                    r["cells_parked"],
                    r["live_cells_end"],
                    r["steps"],
                    r["rejects"],
                    r["nf"],
                    r["njacs"],
                    r["nonlinfail"],
                    fmtv(r["warm_median_s"]),
                    fmtv(r["mass_ledger"]),
                    r["qoi_count"],
                    r["qoi_failures"],
                    fmtv(r["qoi_max_error_over_tol"]),
                ),
                " | ",
            ),
            " |",
        )
    end

    println(
        md,
        "\nQoIs per case (value and error per variant; `pass` = current / current_override / one_callback):\n",
    )
    for (case_id, _) in CASES
        names = unique(r["qoi"] for r in QOI_ROWS if r["case"] == case_id)
        isempty(names) && continue
        println(md, "**$case_id**\n")
        println(
            md,
            "| QoI | tol | exact | ",
            join(("$v: value" for v in VARIANT_NAMES), " | "),
            " | ",
            join(("$v: error" for v in VARIANT_NAMES), " | "),
            " | pass |",
        )
        println(md, "|---|---|---|---|---|---|---|---|---|---|")
        for name in names
            qs = [qoi_of(case_id, v, name) for v in VARIANT_NAMES]
            println(
                md,
                "| $name | ",
                fmtv(qs[1]["tol"]),
                " | ",
                fmtv(qs[1]["exact"]),
                " | ",
                join((@sprintf("%.8g", q["value"]) for q in qs), " | "),
                " | ",
                join((fmtv(q["error"]) for q in qs), " | "),
                " | ",
                join((string(q["pass"]) for q in qs), "/"),
                " |",
            )
        end
        println(md)
    end

    # Termination boundary case.
    println(md, "## Termination with min_cells = 20 (V-BO-04/min20)\n")
    S04 = _vbo04_setup(VBO04.ci.n_cells)
    t_cross = [log(2^k * 20) * S04.ρ_s / S04.r for k = 0:2]
    println(
        md,
        "Every cell of V-BO-04 shrinks as exp(−r t/ρ_s) and the depletion reference is the initial cell thickness, so the k-th ",
        "cascade crosses the 0.05 threshold at t_k = ln(2^k · 20) ρ_s/r = ",
        join((@sprintf("%.3f", t) for t in t_cross), ", "),
        " s. The first cascade leaves exactly 20 live cells; at the second crossing `deplete!` must find `n_live <= min_cells` and stop the solve. ",
        "A `DiscreteCallback` detects a crossing on the first accepted step after it, so the recorded termination time is the crossing plus up to one step.\n",
    )
    println(
        md,
        "| variant | retcode | t_last s | t_last − t₂ s | last step s | last park s | park invocations | cells parked | live cells end | steps |",
    )
    println(md, "|---|---|---|---|---|---|---|---|---|---|")
    for v in VARIANT_NAMES
        r = row_of("V-BO-04/min20", v)
        println(
            md,
            @sprintf(
                "| %s | %s | %.6f | %.4f | %s | %.6f | %d | %d | %d | %d |",
                v,
                r["retcode"],
                r["t_last_accepted"],
                r["t_last_accepted"] - t_cross[2],
                fmtv(r["last_step_s"]),
                r["last_park_t"],
                r["park_invocations"],
                r["cells_parked"],
                r["live_cells_end"],
                r["steps"]
            )
        )
    end
    let a = row_of("V-BO-04/min20", "current"), b = row_of("V-BO-04/min20", "one_callback")
        dt = b["t_last_accepted"] - a["t_last_accepted"]
        same_event =
            a["live_cells_end"] == b["live_cells_end"] &&
            a["cells_parked"] == b["cells_parked"] &&
            a["t_last_accepted"] > t_cross[2] &&
            b["t_last_accepted"] > t_cross[2] &&
            a["t_last_accepted"] - t_cross[2] <= a["last_step_s"] &&
            b["t_last_accepted"] - t_cross[2] <= b["last_step_s"]
        println(
            md,
            "\nBoth variants terminate at the second crossing (t₂ = ",
            @sprintf("%.3f", t_cross[2]),
            " s) with ",
            a["live_cells_end"],
            " / ",
            b["live_cells_end"],
            " live cells after ",
            a["cells_parked"],
            " / ",
            b["cells_parked"],
            " parks: `current` on its first accepted step after t₂ at ",
            @sprintf("%.6f", a["t_last_accepted"]),
            " s (",
            @sprintf("%.4f", a["t_last_accepted"] - t_cross[2]),
            " s after the crossing, last step ",
            fmtv(a["last_step_s"]),
            " s), `one_callback` at ",
            @sprintf("%.6f", b["t_last_accepted"]),
            " s (",
            @sprintf("%.4f", b["t_last_accepted"] - t_cross[2]),
            " s after, last step ",
            fmtv(b["last_step_s"]),
            " s). The ",
            fmtv(abs(dt)),
            " s between the two recorded times is ",
            same_event ?
            "within one accepted step of each variant, i.e. the discrete callback's detection delay on different step sequences (the cascade takes 20 steps of DT_AFTER_PARK in `current` and one in `one_callback`), not a different termination event." :
            "NOT explained by the detection delay alone; see the table.",
            " A loop without the predicate re-test (the submitted version) terminated at the first cascade instead, at 299.981468 s (reviewer's probe); the corrected loop does not.",
        )
    end

    # Ledger versus output spacing.
    println(md, "\n## Mass ledger versus output spacing (V-BO-01)\n")
    println(
        md,
        "The same solves with the saved output times at 0.5 s (the case's grid), 0.125 s and 0.03125 s. Steps and parks do not depend on the output grid; the ledger does.\n",
    )
    println(
        md,
        "| variant | output spacing s | outputs | mass ledger | mass drop kg/m² | final mass kg/m² | steps | cells parked |",
    )
    println(md, "|---|---|---|---|---|---|---|---|")
    for r in LEDGER_ROWS
        println(
            md,
            @sprintf(
                "| %s | %s | %d | %s | %.10f | %.10f | %d | %d |",
                r["variant"],
                r["output_spacing_s"],
                r["n_outputs"],
                fmtv(r["mass_ledger"]),
                r["mass_drop"],
                r["final_mass"],
                r["steps"],
                r["cells_parked"]
            )
        )
    end
    let lc = filter(r -> r["variant"] == "current_override", LEDGER_ROWS),
        lo = filter(r -> r["variant"] == "one_callback", LEDGER_ROWS)

        if !isempty(lc) && !isempty(lo)
            println(
                md,
                "\nHalving the spacing twice moves the ledger from ",
                fmtv(lc[1]["mass_ledger"]),
                " / ",
                fmtv(lo[1]["mass_ledger"]),
                " (current / one_callback at 0.5 s) to ",
                fmtv(lc[end]["mass_ledger"]),
                " / ",
                fmtv(lo[end]["mass_ledger"]),
                " (at 0.03125 s), while the final mass of the two variants agrees to ",
                fmtv(
                    maximum(
                        abs(lc[k]["final_mass"] - lo[k]["final_mass"]) for
                        k in eachindex(lc)
                    ),
                ),
                " kg/m² and the mass drop to ",
                fmtv(
                    maximum(
                        abs(lc[k]["mass_drop"] - lo[k]["mass_drop"]) for k in eachindex(lc)
                    ),
                ),
                " kg/m²: the ledger difference between the variants is output-grid quadrature of the sampled MLR, not a transfer loss.",
            )
        end
    end

    # What changed.
    println(md, "\n## What one-callback parking changed\n")
    let same = all(
            all(
                row_of(c, "current")[k] == row_of(c, "current_override")[k] for k in (
                    "steps",
                    "rejects",
                    "nf",
                    "njacs",
                    "nonlinfail",
                    "cells_parked",
                    "live_cells_end",
                )
            ) for (c, _) in CASES
        )
        println(
            md,
            "- `current_override` ",
            same ? "reproduces" : "does NOT reproduce",
            " `current` step for step on all $(length(CASES)) cases (steps, rejects, nf, njacs, Newton failures, parks, live cells), so the replacement is a faithful copy of the package method with the loop as its only addition.",
        )
    end
    for (case_id, case) in CASES
        a = row_of(case_id, "current")
        b = row_of(case_id, "one_callback")
        qa = [q for q in QOI_ROWS if q["case"] == case_id && q["variant"] == "current"]
        qb = [q for q in QOI_ROWS if q["case"] == case_id && q["variant"] == "one_callback"]
        moved = String[]
        for (x, y) in zip(qa, qb)
            x["error"] == y["error"] && continue
            push!(
                moved,
                @sprintf(
                    "%s: value %s → %s, error %s → %s (tol %s, %s)",
                    first(x["qoi"], 60),
                    @sprintf("%.8g", x["value"]),
                    @sprintf("%.8g", y["value"]),
                    fmtv(x["error"]),
                    fmtv(y["error"]),
                    fmtv(x["tol"]),
                    y["pass"] ? "within tolerance" : "FAILS"
                )
            )
        end
        println(
            md,
            @sprintf(
                "- %s: park invocations %d → %d, cells parked %d → %d, accepted steps %d → %d, f-evals %d → %d, Newton failures %d → %d, warm time %s → %s s, retcode %s → %s (t_last %.6f → %.6f s), QoI failures %d → %d. %s",
                case_id,
                a["park_invocations"],
                b["park_invocations"],
                a["cells_parked"],
                b["cells_parked"],
                a["steps"],
                b["steps"],
                a["nf"],
                b["nf"],
                a["nonlinfail"],
                b["nonlinfail"],
                fmtv(a["warm_median_s"]),
                fmtv(b["warm_median_s"]),
                a["retcode"],
                b["retcode"],
                a["t_last_accepted"],
                b["t_last_accepted"],
                a["qoi_failures"],
                b["qoi_failures"],
                isempty(moved) ? "No QoI error changed." :
                "QoIs whose error changed: " * join(moved, "; ") * "."
            ),
        )
    end
    println(
        md,
        "- V-BO-04/min20: the \"park cascade: 40 → 5 live cells\" QoI compares against 5 live cells and cannot pass with `min_cells = 20` under any variant; it is listed for completeness, the termination time is the quantity of interest there.",
    )
end
say("wrote ", joinpath(RESULTS_DIR, "one_callback_parking.md"))
