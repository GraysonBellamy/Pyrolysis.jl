# The accuracy the requested tolerance buys on the audit decks: each deck on
# a tolerance ladder against its own tight reference, in both energy forms,
# with the cost of every rung.
#
# The decks (`problem_for` below) are a PMMA column on a fixed mesh
# (`eulerian`) and on a moving one (`ale`) and a permeable wood column with
# Darcy–Fick transport on a moving mesh (`ale_darcy`), 80 s, three outputs,
# at 25 and 100 cells. The two PMMA decks give the solid a gas-transfer
# coefficient of 1e-6 m²/s: with none, gas transport vanishes in the solid,
# the gas forms a sharp front, and no tolerance converges the solution (at
# 100 cells reltol 1e-10 and 1e-11 differ by 0.02 K, above the error at the
# default tolerances), so an error against a reference measures nothing
# there. The original PMMA decks are timed at the default tolerances beside
# them.
#
# Each deck runs the default integrator on a tolerance ladder (reltol from
# 1e-4 to 1e-8, abstol a hundredth of it) in both energy forms. The achieved
# error is the largest temperature difference over the saved times and the
# cells against the same deck and form at reltol `REFERENCE` (1e-9), so it
# is the time error alone. The reference spread is the difference between
# that reference and one at a tenth of its tolerances; an error near it is
# not resolved. Per deck and form the summary gives the spread, the error,
# steps and warm time at the default tolerances (reltol 1e-6), and the slope
# of log error against log reltol over the rungs whose error is at least ten
# times the spread (one when the error is proportional to the tolerance).
#
# Warm medians of `REPEATS` solves after one warm-up, one thread.
#
# Run from the repository root when the machine is least loaded:
#   julia --project=test --startup-file=no --threads=1 test/gates/audit_decks.jl
# Writes test/gates/measurements/audit_decks.{csv,md}, audit_decks_summary.csv
# and audit_decks_degenerate.csv.

using Pyrolysis
using Statistics: median

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const REPEATS = 3
# (kind, cells, the solid's gas-transfer coefficient: `nothing` keeps the deck's)
const DECKS = (
    (:eulerian, 25, 1e-6),
    (:eulerian, 100, 1e-6),
    (:ale, 25, 1e-6),
    (:ale, 100, 1e-6),
    (:ale_darcy, 25, nothing),
    (:ale_darcy, 100, nothing),
)
const DEGENERATE_DECKS = ((:eulerian, 25), (:eulerian, 100), (:ale, 25), (:ale, 100))
const FORMS = (:temperature, :enthalpy)
const LADDER = (1e-4, 3e-5, 1e-5, 3e-6, 1e-6, 3e-7, 1e-7, 3e-8, 1e-8)
const DEFAULT_RELTOL = 1e-6
const REFERENCE = 1e-9
const RESOLVED = 10.0
_tol(reltol) = (reltol = reltol, abstol = reltol / 100)

# A deck: a 10 mm column under 50 kW/m² with re-radiation and convection at
# the top and an adiabatic back, `n` cells, `t_end` [s] with `outputs` saved
# times. PMMA (`eulerian`, `ale`) absorbs at its surface; the wood
# (`ale_darcy`) absorbs in depth and carries Darcy–Fick gas transport with a
# pressure condition at the top. Returns the problem and its solve keywords.
function problem_for(kind, n; t_end = 80.0, outputs = 3, solid_λ = nothing, flux = 50e3)
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
        mesh = create_uniform_mesh(0.01, n),
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
    )
end

# CSV output: a field holding a comma, a double quote or a line break is
# quoted (quotes doubled); everything else is written bare.
function csv_field(x::AbstractString)
    needs = any(c -> c == ',' || c == '"' || c == '\n' || c == '\r', x)
    return needs ? '"' * replace(x, '"' => "\"\"") * '"' : String(x)
end
csv_field(x) = csv_field(string(x))
csv_row(io, xs...) = println(io, join((csv_field(x) for x in xs), ','))

# The median of `REPEATS` timed solves, after a warm-up the caller ran.
function timed(f)
    times = Float64[]
    for _ = 1:REPEATS
        GC.gc()
        push!(times, @elapsed f())
    end
    return median(times)
end

time_error(sol, ref) = maximum(abs, sol.T .- ref.T)

# The ladder of one deck in one energy form, every rung timed.
function form_rows(deck, form, solve_at)
    ref = solve_at(_tol(REFERENCE))
    finer = solve_at(_tol(REFERENCE / 10))
    spread = time_error(ref, finer)
    rows = []
    for reltol in LADDER
        sol = solve_at(_tol(reltol))
        err = sol.retcode === :Success ? time_error(sol, ref) : NaN
        push!(
            rows,
            (;
                deck,
                form,
                reltol,
                retcode = sol.retcode,
                error = err,
                spread,
                steps = sol.statistics.accepted,
                time = timed(() -> solve_at(_tol(reltol))),
            ),
        )
        println(join(("$k=$v" for (k, v) in pairs(rows[end])), " "))
        flush(stdout)
    end
    return rows
end

_deck(kind, n, solid_λ) =
    solid_λ === nothing ? "$(kind) $(n)" : "$(kind) $(n), solid λ $(solid_λ)"

function deck_rows(kind, n, solid_λ)
    problem, kw = problem_for(kind, n; solid_λ)
    deck = _deck(kind, n, solid_λ)
    at(energy_state) = tol -> solve(problem; energy_state, kw..., tol...)
    return vcat((form_rows(deck, form, at(form)) for form in FORMS)...)
end

# The original PMMA decks at the default tolerances: warm times only.
function degenerate_rows(kind, n)
    problem, kw = problem_for(kind, n)
    rows = []
    for form in FORMS
        f = () -> solve(problem; energy_state = form, kw...)
        sol = f()
        push!(
            rows,
            (;
                deck = "$(kind) $(n)",
                form,
                retcode = sol.retcode,
                steps = sol.statistics.accepted,
                time = timed(f),
            ),
        )
        println(join(("$k=$v" for (k, v) in pairs(rows[end])), " "))
        flush(stdout)
    end
    return rows
end

# The least-squares slope of log error against log reltol over the rungs
# whose error is at least `RESOLVED` times the spread (`NaN` with fewer than
# two).
function tolerance_slope(rows)
    use = filter(r -> isfinite(r.error) && r.error >= RESOLVED * r.spread, rows)
    length(use) >= 2 || return NaN
    x = [log(r.reltol) for r in use]
    y = [log(r.error) for r in use]
    x̄, ȳ = sum(x) / length(x), sum(y) / length(y)
    return sum((x .- x̄) .* (y .- ȳ)) / sum(abs2, x .- x̄)
end

include(joinpath(@__DIR__, "stamp.jl"))
if abspath(PROGRAM_FILE) == @__FILE__()
    load_start = Sys.loadavg()
    started = time()
    rows = []
    for (kind, n, solid_λ) in DECKS
        append!(rows, deck_rows(kind, n, solid_λ))
    end
    degenerate = []
    for (kind, n) in DEGENERATE_DECKS
        append!(degenerate, degenerate_rows(kind, n))
    end
    load_end = Sys.loadavg()

    summary = []
    for (kind, n, solid_λ) in DECKS, form in FORMS
        deck = _deck(kind, n, solid_λ)
        ladder = filter(r -> r.deck == deck && r.form == form, rows)
        default = only(filter(r -> r.reltol == DEFAULT_RELTOL, ladder))
        push!(
            summary,
            (;
                deck,
                form,
                spread = default.spread,
                error_at_default = default.error,
                steps_at_default = default.steps,
                time_at_default = default.time,
                slope = tolerance_slope(ladder),
            ),
        )
        println(join(("$k=$v" for (k, v) in pairs(summary[end])), " "))
    end

    function csv_out(path, rs)
        open(path, "w") do io
            csv_row(io, string.(keys(first(rs)))...)
            for r in rs
                csv_row(io, values(r)...)
            end
        end
    end
    csv_out(joinpath(OUT, "audit_decks.csv"), rows)
    csv_out(joinpath(OUT, "audit_decks_summary.csv"), summary)
    csv_out(joinpath(OUT, "audit_decks_degenerate.csv"), degenerate)

    _g(x; d = 3) =
        x isa AbstractFloat ? (isnan(x) ? "" : string(round(x; sigdigits = d))) : string(x)
    open(joinpath(OUT, "audit_decks.md"), "w") do io
        println(io, "# The audit decks: accuracy per tolerance\n")
        println(
            io,
            "Produced by `test/gates/audit_decks.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
            "$(Threads.nthreads()) thread, $(Sys.CPU_THREADS) logical CPUs; load average ",
            "$(round.(load_start; digits = 1)) at the start and $(round.(load_end; digits = 1)) ",
            "at the end. The default integrator; error: the largest temperature difference [K] ",
            "over the saved times and cells against the same deck and form at reltol ",
            "$(REFERENCE) (abstol a hundredth of reltol); spread: that reference against one ",
            "at a tenth of its tolerances. Slope: of log error against log reltol over the rungs ",
            "whose error is at least $(Int(RESOLVED)) times the spread. Warm medians of ",
            "$(REPEATS) solves [s].\n",
        )
        println(
            io,
            "| deck | form | spread [K] | error at reltol $(DEFAULT_RELTOL) [K] | steps | warm [s] | slope |",
        )
        println(io, "|" * repeat("---|", 7))
        for s in summary
            println(
                io,
                "| $(s.deck) | $(s.form) | $(_g(s.spread)) | $(_g(s.error_at_default)) | ",
                "$(s.steps_at_default) | $(_g(s.time_at_default)) | $(_g(s.slope)) |",
            )
        end
        println(
            io,
            "\n## The original PMMA decks (solid λ = 0) at the default tolerances\n\n",
            "No tolerance converges these decks (see the script's header); warm times only.\n",
        )
        println(io, "| deck | form | retcode | steps | warm [s] |")
        println(io, "|---|---|---|---|---|")
        for r in degenerate
            println(
                io,
                "| $(r.deck) | $(r.form) | $(r.retcode) | $(r.steps) | $(_g(r.time)) |",
            )
        end
        println(io, "\n## Ladders\n")
        println(io, "| deck | form | reltol | retcode | error [K] | steps | warm [s] |")
        println(io, "|---|---|---|---|---|---|---|")
        for r in rows
            println(
                io,
                "| $(r.deck) | $(r.form) | $(r.reltol) | $(r.retcode) | $(_g(r.error)) | ",
                "$(r.steps) | $(_g(r.time)) |",
            )
        end
    end
    println(
        "wrote $(OUT); load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
    )

end
