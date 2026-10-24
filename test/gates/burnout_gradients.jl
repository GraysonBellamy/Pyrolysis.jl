# Whether finite-difference gradients of V-BO-01's event-timed quantities
# are smooth through burnout.
#
# Four configurations: the PMMA absorption coefficient at 60 and 120 PMMA
# cells, the incident flux and the pre-exponential factor at 60, each at the
# case's tolerance and again at reltol 1e-7 (abstol scaled with it). Each
# parameter is perturbed by ε ∈ {0, ±1, ±2, ±4, ±8}·1e-3 (relative), and the
# central differences D(h) = (Q(h) − Q(−h))/(2h) = dQ/d ln p are taken at the
# four steps h ∈ {1, 2, 4, 8}·1e-3. Each run solves the case through its own
# solve at the package defaults.
#
# Per quantity: the spread of the four differences, (max − min)/|mean|
# (smooth below 0.1, moderate below 1, erratic above); σ, the RMS residual of
# a least-squares quadratic in ε over the nine runs, and h10 = 10σ/|slope|,
# the relative step at which the noise makes a central difference about
# 10 % wrong. A configuration passes when the interpolated peak, its time,
# the interpolated burnout time and the invocation time are smooth and the
# invocation count is the same in every run. A safety-net fallback is an
# invocation whose smallest merged fraction lies below ε_m (1 − 1e-4): a
# located event takes its crossing cell at the threshold, the end-of-step
# event below it. The whole-curve misfits against the Gpyro V2 curves (what
# an inverse fit sees) are reported beside.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/gates/burnout_gradients.jl
# Writes test/gates/measurements/burnout_gradients.{csv,md}.

ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
include(joinpath(@__DIR__, "..", "verification", "casedefs.jl"))
include(joinpath(@__DIR__, "stamp.jl"))

using Logging
using Printf
using Statistics: mean

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const CASE =
    only(filter(c -> c.id == "V-BO-01", Base.invokelatest(load_verification_cases)))
const EPS = [0.0, 1e-3, -1e-3, 2e-3, -2e-3, 4e-3, -4e-3, 8e-3, -8e-3]
const H = [1e-3, 2e-3, 4e-3, 8e-3]
const CASE_RELTOL = CASE.solve_kwargs.reltol
const TIGHT_RELTOL = 1e-7
# A configuration: the parameter, the PMMA cells and the relative tolerance.
const CONFIGS = [
    (param, n, reltol) for reltol in (CASE_RELTOL, TIGHT_RELTOL) for
    (param, n) in ((:pmma_alpha, 60), (:flux, 60), (:A, 60), (:pmma_alpha, 120))
]
# The case's keywords at `reltol`, abstol scaled with it.
case_keywords(reltol) = merge(
    CASE.solve_kwargs,
    (; reltol, abstol = CASE.solve_kwargs.abstol * reltol / CASE_RELTOL),
)
const P0 = Dict(:pmma_alpha => 2870.0 / VBO_RHO_PMMA, :flux => 50e3, :A => 1.35e11)
const LABEL = Dict(
    :pmma_alpha => "PMMA absorption coefficient",
    :flux => "incident flux",
    :A => "pre-exponential factor",
)
const SMOOTH = 0.1
const JUDGED = (:peak_interp, :t_peak_interp, :t_bo_interp, :t_event_first)
const QOIS =
    (:peak_interp, :t_peak_interp, :t_bo_interp, :t_event_first, :misfit_mlr, :misfit_Tback)
const COUNTS = (:invocations, :fallbacks, :steps, :rejected)

# The Gpyro V2 curves: MLR [g m⁻² s⁻¹] and back temperature [K].
const GPYRO = let
    _, cols = Base.invokelatest(
        read_reference_csv,
        "gpyro/vbo01_gpyro_v2_pmma.csv";
        header_rows = 1,
    )
    (; t = cols[1], mlr = cols[2], T_back = cols[4] .+ 273.15)
end

# ----------------------------------------------------------------------------
# The quantities of a solve
# ----------------------------------------------------------------------------

# The peak by a parabola through the largest saved rate and its neighbours,
# and the burnout time by linear interpolation of the live thickness across
# the substrate's thickness plus 0.02 mm.
function interpolated(sol)
    mlr = sol.extras.gas_generation_rate .* 1000
    t = sol.t
    i = argmax(mlr)
    t_pk, pk = t[i], mlr[i]
    if 1 < i < length(mlr)
        y0, y1, y2 = mlr[i-1], mlr[i], mlr[i+1]
        d = y0 - 2y1 + y2
        δ = d == 0 ? 0.0 : 0.5 * (y0 - y2) / d
        t_pk = t[i] + δ * (t[i+1] - t[i])
        pk = y1 - 0.25 * (y0 - y2) * δ
    end
    th = sol.extras.thickness_history
    L_bo = VBO_L_KW + 0.02e-3
    j = findfirst(<(L_bo), th)
    t_bo =
        j === nothing || j == 1 ? NaN :
        t[j-1] + (th[j-1] - L_bo) / (th[j-1] - th[j]) * (t[j] - t[j-1])
    return (; peak_interp = pk, t_peak_interp = t_pk, t_bo_interp = t_bo)
end

function lininterp(x, xs, ys)
    j = searchsortedlast(xs, x)
    j == 0 && return ys[1]
    j >= length(xs) && return ys[end]
    return ys[j] + (ys[j+1] - ys[j]) * (x - xs[j]) / (xs[j+1] - xs[j])
end

# The mean square misfits of the MLR and the back temperature against Gpyro
# V2 over the solve.
function misfits(sol)
    mlr = sol.extras.gas_generation_rate .* 1000
    T_back = sol.extras.back_temperature
    sel = [k for k in eachindex(GPYRO.t) if GPYRO.t[k] <= sol.t[end]]
    return (;
        misfit_mlr = mean(
            abs2,
            lininterp(GPYRO.t[k], sol.t, mlr) - GPYRO.mlr[k] for k in sel
        ),
        misfit_Tback = mean(
            abs2,
            lininterp(GPYRO.t[k], sol.t, T_back) - GPYRO.T_back[k] for k in sel
        ),
    )
end

# The first invocation's time, the invocation count and the safety-net
# fallbacks.
function events(sol)
    ε_m = VBO_STEP_STATS["V-BO-01"].threshold
    inv = Base.invokelatest(_vbo_invocations, sol, ε_m)
    fallbacks = count(
        k ->
            minimum(e.fraction for e in sol.events if e.invocation == k) < ε_m * (1 - 1e-4),
        unique(e.invocation for e in sol.events),
    )
    t_first = isempty(inv) ? NaN : first(inv).t
    return (; t_event_first = t_first, invocations = length(inv), fallbacks)
end

function perturbed_solve(param, n, reltol, ε)
    p = NamedTuple{(param,)}((P0[param] * (1 + ε),))
    problem = Base.invokelatest(_vbo_problem, n; p...)
    kw = case_keywords(reltol)
    sol = Logging.with_logger(Logging.NullLogger()) do
        Base.invokelatest(CASE.solve_fn, problem; kw...)
    end
    st = sol.raw_solution.stats
    base = (;
        param,
        n,
        reltol,
        ε,
        retcode = sol.retcode,
        steps = st.naccept,
        rejected = st.nreject,
    )
    sol.retcode == :Success || return merge(
        base,
        (; (q => NaN for q in QOIS)..., invocations = -1, fallbacks = -1),
    )
    return merge(base, interpolated(sol), misfits(sol), events(sol))
end

# ----------------------------------------------------------------------------
# The differences
# ----------------------------------------------------------------------------

function spread(ds)
    all(isfinite, ds) || return NaN
    m = mean(ds)
    s = maximum(ds) - minimum(ds)
    return m == 0 ? (s == 0 ? 0.0 : Inf) : s / abs(m)
end

classify(s) = !isfinite(s) ? "n/a" : s < SMOOTH ? "smooth" : s < 1 ? "moderate" : "erratic"

# The slope and the residual RMS of a least-squares quadratic in ε.
function quadratic(x, y)
    X = hcat(ones(length(x)), x, x .^ 2)
    β = X \ y
    return β[2], sqrt(sum(abs2, y .- X * β) / (length(x) - 3))
end

function differences(rows, q)
    value(ε) = Float64(getproperty(only(r for r in rows if r.ε == ε), q))
    ds = [(value(h) - value(-h)) / (2h) for h in H]
    slope, σ = quadratic(EPS, value.(EPS))
    return (; q0 = value(0.0), ds, spread = spread(ds), slope, σ, h10 = 10σ / abs(slope))
end

# ----------------------------------------------------------------------------
# The runs and the tables
# ----------------------------------------------------------------------------

load_start = Sys.loadavg()
started = time()
perturbed_solve(:pmma_alpha, 60, CASE_RELTOL, 0.0)
rows = Any[]
for (param, n, reltol) in CONFIGS, ε in EPS
    r = perturbed_solve(param, n, reltol, ε)
    push!(rows, r)
    @printf(
        "%s n=%d reltol=%.0e ε=%+.0e: %s steps=%d invocations=%d fallbacks=%d t_bo=%.6f t_event=%.6f\n",
        param,
        n,
        reltol,
        ε,
        r.retcode,
        r.steps,
        r.invocations,
        r.fallbacks,
        r.t_bo_interp,
        r.t_event_first
    )
    flush(stdout)
end
load_end = Sys.loadavg()

open(joinpath(OUT, "burnout_gradients.csv"), "w") do io
    cols = (:param, :n, :reltol, :ε, :retcode, QOIS..., COUNTS...)
    println(io, join(string.(cols), ","))
    for r in rows
        println(io, join((string(getproperty(r, c)) for c in cols), ","))
    end
end

_g(x; d = 3) =
    x isa AbstractFloat ? (isfinite(x) ? string(round(x; sigdigits = d)) : string(x)) :
    string(x)

verdicts = Dict{Any,Bool}()
open(joinpath(OUT, "burnout_gradients.md"), "w") do io
    println(io, "# Smooth finite-difference gradients through burnout\n")
    println(
        io,
        "Produced by `test/gates/burnout_gradients.jl` at $(COMMIT_STAMP). Julia $(VERSION), ",
        "$(Threads.nthreads()) thread; load average $(round.(load_start; digits = 1)) at the start and ",
        "$(round.(load_end; digits = 1)) at the end. V-BO-01 (temperature form) ",
        "through the case's solve at the package defaults; ",
        "each configuration at the case's tolerance (reltol $(CASE_RELTOL)), then at ",
        "reltol $(TIGHT_RELTOL) (abstol scaled with it). ",
        "Each parameter perturbed by ε ∈ {0, ±1, ±2, ±4, ±8}·1e-3; D(h) = (Q(h) − Q(−h))/(2h) at ",
        "h ∈ {1, 2, 4, 8}·1e-3. Spread: (max − min)/|mean| of the four differences (smooth below ",
        "$(SMOOTH), moderate below 1); σ: the RMS residual of a quadratic in ε over the nine runs; ",
        "h10 = 10σ/|slope|. A configuration passes when the interpolated peak, its time, the interpolated ",
        "burnout time and the invocation time are smooth and the ",
        "invocation count is the same in every run. Misfits: mean squares against the Gpyro V2 MLR ",
        "[(g m⁻² s⁻¹)²] and back temperature [K²], reported only.\n",
    )
    for (param, n, reltol) in CONFIGS
        tolerance = reltol == CASE_RELTOL ? "the case's tolerance" : "reltol $(reltol)"
        println(io, "## $(LABEL[param]), $(n) PMMA cells, $(tolerance)\n")
        rs = filter(r -> r.param == param && r.n == n && r.reltol == reltol, rows)
        counts = unique(r.invocations for r in rs)
        fallbacks = sum(r.fallbacks for r in rs)
        steps = only(r.steps for r in rs if r.ε == 0.0)
        println(
            io,
            "Every run $(join(unique(string.(r.retcode for r in rs)), ", ")); ",
            "invocations $(counts), safety-net fallbacks $(fallbacks) of $(length(rs))",
            "; $(steps) steps at ε = 0 ($(minimum(r.steps for r in rs)) to $(maximum(r.steps for r in rs))).\n",
        )
        println(
            io,
            "| quantity | Q(0) | D(1e-3) | D(2e-3) | D(4e-3) | D(8e-3) | spread | class | σ | h10 |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|")
        smooth = true
        for q in QOIS
            d = differences(rs, q)
            q in JUDGED && (smooth &= d.spread < SMOOTH)
            println(
                io,
                "| $(q == :t_event_first ? "invocation time" : q) | ",
                "$(_g(d.q0; d = 6)) | $(join(_g.(d.ds), " | ")) | $(_g(d.spread)) | ",
                "$(classify(d.spread)) | $(_g(d.σ)) | $(_g(d.h10)) |",
            )
        end
        println(io)
        verdicts[(param, n, reltol)] = smooth && length(counts) == 1
        println(
            io,
            "This configuration: **$(verdicts[(param, n, reltol)] ? "pass" : "fail")**.\n",
        )
    end
end
println(
    "wrote $(OUT); the gradients at the case's tolerance ",
    all(v for (k, v) in verdicts if k[3] == CASE_RELTOL) ? "pass" : "fail",
    ", at reltol $(TIGHT_RELTOL) ",
    all(v for (k, v) in verdicts if k[3] == TIGHT_RELTOL) ? "pass" : "fail",
    "; load $(load_start) -> $(load_end); $(round(Int, time() - started)) s",
)
