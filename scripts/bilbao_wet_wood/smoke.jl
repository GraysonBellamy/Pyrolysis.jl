# One run: solve, closures, and a look against the data.
#   julia --project=scripts/verification scripts/bilbao_wet_wood/smoke.jl run [key=value ...]
# keys: law, kinetics, shrinkage, flux_factor, k_scale, κ_scale, h_front, cells, ...
include(joinpath(@__DIR__, "BilbaoWetWood.jl"))
using .BilbaoWetWood
using Printf

run = get(ARGS, 1, "dry_E")
kw = Dict{Symbol,Any}()
for a in ARGS[2:end]
    k, v = split(a, '=')
    kw[Symbol(k)] =
        k in ("law", "kinetics", "front") ? Symbol(v) :
        (k == "cells" ? parse(Int, v) : parse(Float64, v))
end
r = RUNS[run]
t0 = time()
s = simulate(r; kw...)
@printf(
    "%s %s: %s in %.1f s (wall %.1f s); mass closure %.1e, energy %.1e\n",
    run,
    join(ARGS[2:end], " "),
    s.retcode,
    s.solve_time,
    time() - t0,
    s.audit.mass.relative,
    s.audit.energy.relative
)
meas = temperatures(run)
tX, X = conversion(run)
println(
    "  t   | T1 mod  meas | T9 mod  meas | T19 mod  meas | X_T mod  meas | Ts    L/mm  pmax/kPa",
)
for tm in (2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 24, 28, 32, 36, 40, 44, 50, 55, 60, 66)
    ts = 60.0tm
    ts > s.t[end] && break
    g(d) =
        haskey(meas, d) ?
        (i = findfirst(==(ts), meas[d][1]); i === nothing ? NaN : meas[d][2][i]) : NaN
    gx = (i = findfirst(==(ts), tX); i === nothing ? NaN : X[i])
    @printf(
        "%4d  | %5.0f %5.0f | %5.0f %5.0f | %5.0f  %5.0f  | %.3f  %.3f  | %5.0f %5.1f  %6.1f\n",
        tm,
        interp(s.t, s.T1, ts),
        g(1),
        interp(s.t, s.T9, ts),
        g(9),
        interp(s.t, s.T19, ts),
        g(19),
        interp(s.t, s.X_T, ts),
        gx,
        interp(s.t, s.T_surface, ts),
        1000interp(s.t, s.thickness, ts),
        (interp(s.t, s.p_max, ts) - 101325) / 1000
    )
end
rms = rms_by_depth(run, s)
@printf(
    "rms K: 1 mm %.1f, 9 mm %.1f, 19 mm %.1f\n",
    get(rms, 1, NaN),
    get(rms, 9, NaN),
    get(rms, 19, NaN)
)
if !isempty(tX)
    xm = [interp(s.t, s.X_T, t) for t in tX]
    @printf("rms X_T %.3f\n", sqrt(sum((xm .- X) .^ 2) / length(X)))
end
