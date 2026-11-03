# One complete parameter-to-objective derivative qualified by central
# differences, the supported route to sensitivities. The objective is the
# misfit of V-BO-01's mass-loss rate against the Gpyro V2 reference, the mean
# square of the saved gas-generation rate, interpolated linearly, minus
# Gpyro's at Gpyro's output times within the run [(g m⁻² s⁻¹)²]: the
# objective of `burnout_gradients.jl` and of an inverse fit to that curve.
#
# V-BO-01 at its CI resolution (60 PMMA cells), the temperature form,
# through the case's solve with its keywords, abstol scaled with reltol.
# Each of the PMMA absorption coefficient α, the incident flux q and the
# pre-exponential factor A is perturbed relatively, p = p₀ (1 + ε), ε ∈ {0,
# ±1, ±2, ±4, ±8}·1e-3, at reltol 1e-6, 1e-7 and 1e-8: 81 runs. The central
# differences D(h) = (J(h) − J(−h))/(2h) at h ∈ {1, 2, 4, 8}·1e-3 estimate
# dJ/d ln p; the flux's is also given per W/m².
#
# The verdict per parameter is `fd_convergence.jl`'s: at reltol 1e-8
# a window of three consecutive steps whose two Richardson extrapolants lie
# within 5 % of each other, the converged value the extrapolant of its smaller
# pair, and at each of its steps the reltol 1e-7 and 1e-8 differences within
# 5 %; the usable steps at each tolerance are those whose difference lies
# within 5 % of the converged value. The literal reading (the three largest
# steps' differences within 5 % of one another) is reported beside, and the
# runs' invocation counts and steps.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/gates/fd_gradient_qualification.jl
# Writes test/gates/measurements/fd_gradient_qualification.{csv,md}.

ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "stamp.jl"))
include(joinpath(@__DIR__, "fd_convergence.jl"))

using Logging
using Statistics: mean

const OUT = joinpath(@__DIR__, "measurements")
const CASE =
    only(filter(c -> c.id == "V-BO-01", Base.invokelatest(load_verification_cases)))
const N_CELLS = CASE.ci.n_cells
const EPS = [0.0, 1e-3, -1e-3, 2e-3, -2e-3, 4e-3, -4e-3, 8e-3, -8e-3]
const STEPS = (1e-3, 2e-3, 4e-3, 8e-3)
const RELTOLS = (1e-6, 1e-7, 1e-8)
const PARAMS = (:pmma_alpha, :flux, :A)
const P0 = Dict(:pmma_alpha => 2870.0 / VBO_RHO_PMMA, :flux => 50e3, :A => 1.35e11)
const LABEL = Dict(
    :pmma_alpha => "PMMA absorption coefficient α",
    :flux => "incident flux q",
    :A => "pre-exponential factor A",
)
const CASE_RELTOL = CASE.solve_kwargs.reltol

# The Gpyro V2 mass-loss rate [g m⁻² s⁻¹].
const GPYRO = let
    _, cols = Base.invokelatest(
        read_reference_csv,
        "gpyro/vbo01_gpyro_v2_pmma.csv";
        header_rows = 1,
    )
    (; t = cols[1], mlr = cols[2])
end

function lininterp(x, xs, ys)
    j = searchsortedlast(xs, x)
    j == 0 && return ys[1]
    j >= length(xs) && return ys[end]
    return ys[j] + (ys[j+1] - ys[j]) * (x - xs[j]) / (xs[j+1] - xs[j])
end

# The objective of a solve.
function misfit(sol)
    mlr = sol.extras.gas_generation_rate .* 1000
    sel = [k for k in eachindex(GPYRO.t) if GPYRO.t[k] <= sol.t[end]]
    return mean(abs2, lininterp(GPYRO.t[k], sol.t, mlr) - GPYRO.mlr[k] for k in sel)
end

case_keywords(reltol) = merge(
    CASE.solve_kwargs,
    (; reltol, abstol = CASE.solve_kwargs.abstol * reltol / CASE_RELTOL),
)

function perturbed(param, reltol, ε)
    p = NamedTuple{(param,)}((P0[param] * (1 + ε),))
    problem = Base.invokelatest(_vbo_problem, N_CELLS; p...)
    started = time()
    sol = Logging.with_logger(Logging.NullLogger()) do
        Base.invokelatest(CASE.solve_fn, problem; case_keywords(reltol)...)
    end
    wall = time() - started
    ε_m = VBO_STEP_STATS["V-BO-01"].threshold
    inv = Base.invokelatest(_vbo_invocations, sol, ε_m)
    return (;
        param,
        reltol,
        ε,
        retcode = sol.retcode,
        J = sol.retcode == :Success ? misfit(sol) : NaN,
        invocations = length(inv),
        t_invocation = isempty(inv) ? NaN : first(inv).t,
        steps = sol.statistics.accepted,
        rejected = sol.statistics.rejected,
        wall,
    )
end

# The check of one parameter at one tolerance.
check(rows, param, reltol) = convergence(
    Dict(r.ε => r.J for r in rows if r.param == param && r.reltol == reltol),
    STEPS,
)

function qualification_main()
    load_start, started = Sys.loadavg(), time()
    perturbed(:pmma_alpha, CASE_RELTOL, 0.0)
    rows = Any[]
    for param in PARAMS, reltol in RELTOLS, ε in EPS
        r = perturbed(param, reltol, ε)
        push!(rows, r)
        println(
            "$(param) reltol=$(reltol) ε=$(ε): $(r.retcode) J=$(r.J) invocations=$(r.invocations) ",
            "t_inv=$(r.t_invocation) steps=$(r.steps) wall=$(round(r.wall; digits = 1)) s",
        )
        flush(stdout)
    end
    results = Dict(
        param => begin
            checks = [check(rows, param, reltol) for reltol in RELTOLS]
            (; checks, q = qualification(checks))
        end for param in PARAMS
    )
    write_qualification(rows, results, load_start, Sys.loadavg(), time() - started)
    println(
        "wrote $(OUT)/fd_gradient_qualification.{csv,md}; ",
        join(
            ("$(p) $(results[p].q.pass ? "qualified" : "not qualified")" for p in PARAMS),
            ", ",
        ),
        "; $(round(Int, time() - started)) s",
    )
    return (; rows, results)
end

_g(x; d = 4) =
    x isa AbstractFloat ? (isfinite(x) ? string(round(x; sigdigits = d)) : string(x)) :
    string(x)
_yes(b) = b ? "pass" : "**fail**"
_step(h) = isnan(h) ? "none" : "$(round(Int, 1000h))e-3"
_steps(hs) = isempty(hs) ? "none" : join(_step.(hs), ", ")
const UNIT = Dict(:pmma_alpha => "m²/kg", :flux => "W/m²", :A => "s⁻¹")

function write_qualification(rows, results, load_start, load_end, elapsed)
    mkpath(OUT)
    open(joinpath(OUT, "fd_gradient_qualification.csv"), "w") do io
        println(
            io,
            "param,reltol,eps,retcode,J,invocations,t_invocation,steps,rejected,wall",
        )
        for r in rows
            println(
                io,
                join(
                    (
                        r.param,
                        r.reltol,
                        r.ε,
                        r.retcode,
                        r.J,
                        r.invocations,
                        r.t_invocation,
                        r.steps,
                        r.rejected,
                        r.wall,
                    ),
                    ",",
                ),
            )
        end
    end
    open(joinpath(OUT, "fd_gradient_qualification.md"), "w") do io
        println(io, "# The finite-difference gradient qualification\n")
        println(
            io,
            "Produced by `test/gates/fd_gradient_qualification.jl` at $(COMMIT_STAMP). Julia ",
            "$(VERSION), $(Threads.nthreads()) thread; load average $(round.(load_start; digits = 1)) at ",
            "the start and $(round.(load_end; digits = 1)) at the end; $(round(Int, elapsed)) s. V-BO-01 at ",
            "$(N_CELLS) PMMA cells, the temperature form, the case's solve and keywords (abstol scaled ",
            "with reltol). J: the mean square misfit of the mass-loss rate against Gpyro V2 ",
            "[(g m⁻² s⁻¹)²]. D(h) = (J(h) − J(−h))/(2h) ≈ dJ/d ln p. The script's header and ",
            "`fd_convergence.jl` define the rule.\n",
        )
        println(io, "## The verdict\n")
        println(
            io,
            "Per parameter, at reltol $(RELTOLS[end]): the converged window (three consecutive steps ",
            "whose two Richardson extrapolants lie within $(Int(100AGREE)) % of each other), the ",
            "converged dJ/d ln p (the extrapolant of the window's smaller pair), its error (the ",
            "difference of the two extrapolants) and per unit of the parameter; the reltol ",
            "$(RELTOLS[end-1]) and $(RELTOLS[end]) differences within $(Int(100AGREE)) % at each ",
            "step of the window; the usable steps at each tolerance (the difference within ",
            "$(Int(100AGREE)) % of the converged value); the verdict; and the literal reading (the three ",
            "largest steps within $(Int(100AGREE)) % of one another, and its value).\n",
        )
        println(
            io,
            "| parameter | J(0) | window | dJ/d ln p | error | per unit | tolerances agree | ",
            join(("usable at reltol $(r)" for r in RELTOLS), " | "),
            " | qualified | literal reading |",
        )
        println(
            io,
            "|---|---|---|---|---|---|---|",
            repeat("---|", length(RELTOLS)),
            "---|---|",
        )
        for p in PARAMS
            q = results[p].q
            J0 = only(
                r.J for r in rows if r.param == p && r.reltol == RELTOLS[end] && r.ε == 0
            )
            literal =
                "$(_yes(q.literal.pass))$(q.literal.agrees ? "" : " (the three largest disagree)")" *
                (isnan(q.literal.value) ? "" : ", $(_g(q.literal.value; d = 5))")
            println(
                io,
                "| $(LABEL[p]) | $(_g(J0; d = 6)) | $(_steps(q.window)) | $(_g(q.value; d = 5)) | ",
                "$(_g(q.error; d = 2)) | $(_g(q.value / P0[p])) per $(UNIT[p]) | $(_yes(q.across)) | ",
                "$(join(_steps.(q.usable), " | ")) | $(_yes(q.pass)) | $(literal) |",
            )
        end
        println(io, "\n## The differences\n")
        println(
            io,
            "D(h) at each step and tolerance; the Richardson extrapolants of successive pairs; the ",
            "tolerance's own window and converged value; the spread of the three largest differences ",
            "((max − min)/min |D|, the literal reading's); the runs' invocations and their accepted ",
            "steps at ε = 0.\n",
        )
        println(
            io,
            "| parameter | reltol | ",
            join(("D($(_step(h)))" for h in STEPS), " | "),
            " | extrapolants | window | converged | spread of the three largest | invocations | steps at ε = 0 |",
        )
        println(io, "|---|---|", repeat("---|", length(STEPS)), "---|---|---|---|---|---|")
        for p in PARAMS, (k, reltol) in enumerate(RELTOLS)
            c = results[p].checks[k]
            rs = [r for r in rows if r.param == p && r.reltol == reltol]
            inv = unique(r.invocations for r in rs)
            steps = only(r.steps for r in rs if r.ε == 0)
            println(
                io,
                "| $(LABEL[p]) | $(reltol) | $(join(_g.(c.D; d = 5), " | ")) | ",
                "$(join(_g.(c.extrapolants; d = 5), ", ")) | $(_steps(c.window)) | ",
                "$(_g(c.value; d = 5)) | $(_g(spread(c.D[end-2:end]); d = 2)) | $(join(inv, ", ")) | ",
                "$(steps) |",
            )
        end
        println(io, "\n## The runs\n")
        println(
            io,
            "| parameter | reltol | ε | retcode | J | invocations (first at) | steps (rejected) | wall [s] |",
        )
        println(io, "|---|---|---|---|---|---|---|---|")
        for r in rows
            println(
                io,
                "| $(LABEL[r.param]) | $(r.reltol) | $(r.ε) | $(r.retcode) | $(_g(r.J; d = 10)) | ",
                "$(r.invocations) ($(_g(r.t_invocation; d = 7))) | $(r.steps) ($(r.rejected)) | ",
                "$(_g(r.wall; d = 3)) |",
            )
        end
    end
end

abspath(PROGRAM_FILE) == @__FILE__() && qualification_main()
