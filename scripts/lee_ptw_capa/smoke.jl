# One run: solve, closures, and a look against the data.
#   julia --project=scripts/verification scripts/lee_ptw_capa/smoke.jl RUN [VARIANT]
# RUN is "63" or "29"; VARIANT is law[/tag] (variants.jl), default `sorption`.
include(joinpath(@__DIR__, "LeePtwCapa.jl"))
using .LeePtwCapa
using Printf
include(joinpath(@__DIR__, "variants.jl"))

run = get(ARGS, 1, "63")
variant = get(ARGS, 2, "sorption")
r = RUNS[run]
t0 = time()
s = simulate(r; variant_kwargs(variant)...)
@printf(
    "%s %s: %s in %.1f s (wall %.1f s); mass closure %.1e, energy %.1e\n",
    run,
    variant,
    s.retcode,
    s.solve_time,
    time() - t0,
    s.audit.mass.relative,
    s.audit.energy.relative
)
Tm = bottom_temperature(run)
Mm = mass_loss_rate(run)
println("   t  | T_bottom mod  meas  paper | MLR mod   meas    paper  | X_back  pmax/kPa")
for ts in (
    run == "63" ? (10, 20, 40, 60, 80, 100, 120, 150, 200, 250, 300, 400, 450) :
    (20, 50, 100, 130, 160, 200, 250, 300, 400, 600, 800, 1000, 1100)
)
    m(x) = interp(Tm.t, x, ts)
    @printf(
        "%5d | %6.1f  %6.1f %6.1f | %.5f %.5f %.5f | %.4f  %6.1f\n",
        ts,
        interp(s.t, s.T_bottom, ts),
        m(Tm.T),
        m(Tm.model),
        interp(s.t, s.MLR, ts),
        interp(Mm.t, Mm.MLR, ts),
        interp(Mm.t, Mm.model, ts),
        interp(s.t, s.X_back, ts),
        (interp(s.t, s.p_max, ts) - 101325) / 1000
    )
end
@printf("mass lost %.3f kg/m² of %.3f\n", s.mass[1] - s.mass[end], s.mass[1])
