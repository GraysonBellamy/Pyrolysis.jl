# One case: solve, closures, and a look against the data.
#   julia --project=scripts/verification scripts/fredlund_wet_wood/smoke.jl CASE [VARIANT]
# CASE is spruce14, pine24 or pine14_1985; VARIANT is law[/tag] (variants.jl);
# then key=value pairs for `simulate` (e.g. k_scale=1.4).
include(joinpath(@__DIR__, "FredlundWetWood.jl"))
using .FredlundWetWood
using .FredlundWetWood: PHASE_CHANGE_LAWS
using Printf
using TOML
include(joinpath(@__DIR__, "variants.jl"))

case = get(ARGS, 1, "spruce14")
variant = get(ARGS, 2, "sorption")
extra = Dict{Symbol,Any}()
for a in ARGS[3:end]
    k, v = split(a, '=')
    extra[Symbol(k)] = k == "front" ? Symbol(v) : parse(Float64, v)
end
t0 = time()
s = simulate(case; case_kwargs(variant, case)..., extra...)
@printf(
    "%s %s: %s in %.1f s (wall %.1f s); mass closure %.1e, energy %.1e\n",
    case,
    variant,
    s.retcode,
    s.solve_time,
    time() - t0,
    s.audit.mass.relative,
    s.audit.energy.relative
)
meas = temperatures(case)
println(
    " t/min | " *
    join((@sprintf("T%-2d mod meas", d) for d in (20, 30, 40, 50, 70)), " | ") *
    " | P10 P30 kPa | X30  X50",
)
for tm = 5:5:round(Int, s.t[end]/60)
    ts = 60.0tm
    m(d) =
        haskey(meas, d) ?
        (
            tt = meas[d];
            ts < tt[1][1] - 60 || ts > tt[1][end] + 60 ? NaN :
            interp(tt[1], tt[2], ts) - 273.15
        ) : NaN
    @printf("%5d  | ", tm)
    for d in (20, 30, 40, 50, 70)
        @printf("%5.0f %5.0f | ", interp(s.t, s.cols[Symbol("T$d")], ts) - 273.15, m(d))
    end
    @printf(
        "%5.0f %5.0f | %.3f %.3f\n",
        interp(s.t, s.cols[:P10], ts) / 1e3,
        interp(s.t, s.cols[:P30], ts) / 1e3,
        interp(s.t, s.cols[:X30], ts),
        interp(s.t, s.cols[:X50], ts)
    )
end
@printf("p_max %.0f kPa gauge, X_max %.3f\n", maximum(s.p_max) / 1e3, maximum(s.X_max))
