# Smoke run: one species at each flux. Solves the species' mass tests (25, 35
# and 45 kW/m²) and temperature tests (35 and 45 kW/m²) on the chamber deck,
# checks that each solve succeeds and that its mass and energy ledgers close
# (ConservationRecorder), and prints the mass loss and two thermocouples
# against the data as a sanity check, not a comparison; then a few
# sensitivities on the 35 kW/m² tests (see `sensitivities`). Writes
# results/smoke/<run>.csv and results/figures/smoke_<species>.png.
#
#   julia --project=scripts/verification scripts/lemta_argon/smoke.jl [species]
include(joinpath(@__DIR__, "LemtaArgon.jl"))
using .LemtaArgon
using CairoMakie
using Printf

species = isempty(ARGS) ? "poplar" : ARGS[1]
check_data()
mr = filter(r -> r.species == species, mass_runs())
tr = filter(r -> r.species == species, temperature_runs())
outdir = joinpath(RESULTS, "smoke")
mkpath(outdir)
mkpath(joinpath(RESULTS, "figures"))

ledger(s) = @sprintf(
    "%s in %.0f s; mass defect %.2g (relative %.2g), energy defect %.3g J/m² (relative %.2g)",
    s.retcode,
    s.solve_time,
    s.audit.mass.defect,
    s.audit.mass.relative,
    s.audit.energy.defect,
    s.audit.energy.relative
)

at(t, y, t0) = interp(t, y, t0)

results = Dict{String,Any}()
println("== $species: mass tests")
for r in sort(mr, by = r -> r.flux)
    s = simulate(r)
    results[r.name] = s
    write_series(joinpath(outdir, r.name * ".csv"), s)
    println(r.name, " (", r.flux, " kW/m²): ", ledger(s))
    for t0 in (300.0, r.t[end])
        @printf(
            "   t = %4.0f s: m/m₀ model %.3f, data %.3f; MLRPUA model %.2f, data %.2f g/(m² s) (data ±15 s mean)\n",
            t0,
            at(s.t, s.m, t0),
            at(r.t, r.m, t0),
            at(s.t, s.mlr, t0),
            sum(r.mlr[abs.(r.t .- t0) .<= 15]) / count(abs.(r.t .- t0) .<= 15)
        )
    end
    flush(stdout)
end

println("== $species: temperature tests")
probes = (0.010, 0.020)
for r in sort(tr, by = r -> r.flux)
    s = simulate(r)
    results[r.name] = s
    write_series(joinpath(outdir, r.name * ".csv"), s)
    println(r.name, " (", r.flux, " kW/m²): ", ledger(s))
    for d in probes
        j = findfirst(x -> isapprox(x, d), DEPTHS)
        js = findfirst(x -> isapprox(x, d), s.depths)
        y = usable(r, j)
        for t0 in (300.0, 600.0)
            t0 > s.t[end] && continue
            @printf(
                "   %2.0f mm at %3.0f s: model %6.1f °C, data %6.1f °C\n",
                1000d,
                t0,
                at(s.t, view(s.T, :, js), t0) - 273.15,
                live(r, j) ? at(r.t, y, t0) - 273.15 : NaN
            )
        end
    end
    flush(stdout)
end

# Sensitivities on the 35 kW/m² tests, against the deck (the face under the
# cone and window, h calibrated on the silicate, 20 % shrinkage): no
# shrinkage, and the species' own shrinkage in the cone tests; the back held
# at the measured 40 mm record; volatiles without sensible heat; h at natural
# convection; and the gray front (front = :silcal), calibrated at the
# silicate absorptivity of front.toml and at the bounds 0.80 and 0.92, and
# with the thesis's h.
function sensitivities(species, mr, tr, results)
    println("== $species: sensitivities at 35 kW/m²")
    r = only(filter(r -> r.flux == 35, tr))
    m = only(filter(r -> r.flux == 35, mr))
    probe(s, d, t0) =
        at(s.t, view(s.T, :, findfirst(x -> isapprox(x, d), s.depths)), t0) - 273.15
    data(d, t0) = at(r.t, usable(r, findfirst(x -> isapprox(x, d), DEPTHS)), t0) - 273.15
    cal = LemtaArgon.TOML.parsefile(LemtaArgon.FRONT)
    cases = (
        ("deck (cone and window, shrinkage 0.20)", (;)),
        ("no shrinkage", (shrinkage = nothing,)),
        ("shrinkage of the cone tests", (shrinkage = :cone,)),
        ("back at the 40 mm record", (back = :tc40,)),
        ("volatiles without sensible heat", (volatiles_c = 0.0,)),
        (
            "cone and window, natural convection",
            (h_front = chamber()["window"]["h_window"],),
        ),
        ("gray front", (front = :silcal,)),
        (
            "gray front, silicate absorptivity 0.80",
            (front = :silcal, h_front = cal["h_at_absorptivity_80"]["35"]),
        ),
        (
            "gray front, silicate absorptivity 0.92",
            (front = :silcal, h_front = cal["h_at_absorptivity_92"]["35"]),
        ),
        ("gray front, thesis h", (front = :thesis,)),
    )
    @printf("%-40s %6s %6s %6s %8s %8s\n", "", "2 mm", "10 mm", "20 mm", "m/m₀", "MLRPUA")
    @printf("%-40s %6s %6s %6s %8s %8s\n", "", "600 s", "600 s", "600 s", "698 s", "300 s")
    for (tag, kw) in cases
        st = isempty(kw) ? results[r.name] : simulate(r; kw...)
        sm =
            isempty(kw) ? results[m.name] :
            (haskey(kw, :back) ? nothing : simulate(m; kw...))
        @printf(
            "%-40s %6.0f %6.0f %6.0f %8s %8s   %s\n",
            tag,
            probe(st, 0.002, 600.0),
            probe(st, 0.010, 600.0),
            probe(st, 0.020, 600.0),
            sm === nothing ? "–" : @sprintf("%.3f", sm.m[end]),
            sm === nothing ? "–" : @sprintf("%.2f", at(sm.t, sm.mlr, 300.0)),
            ledger(st)
        )
        flush(stdout)
    end
    @printf(
        "%-40s %6.0f %6.0f %6.0f %8.3f %8.2f\n",
        "data",
        data(0.002, 600.0),
        data(0.010, 600.0),
        data(0.020, 600.0),
        m.m[end],
        sum(m.mlr[abs.(m.t .- 300) .<= 15]) / count(abs.(m.t .- 300) .<= 15)
    )
end
sensitivities(species, mr, tr, results)

fig = Figure(size = (1200, 800))
ax1 = Axis(fig[1, 1]; title = "$species: m/m₀", xlabel = "time (s)")
ax2 = Axis(fig[1, 2]; title = "MLRPUA (g/(m² s))", xlabel = "time (s)")
for (i, r) in enumerate(sort(mr, by = r -> r.flux))
    s = results[r.name]
    c = Makie.wong_colors()[i]
    lines!(ax1, r.t, r.m; color = c, label = "$(round(Int, r.flux)) data")
    lines!(ax1, s.t, s.m; color = c, linestyle = :dash, label = "model")
    lines!(ax2, r.t, r.mlr; color = (c, 0.5))
    lines!(ax2, s.t, s.mlr; color = c, linestyle = :dash)
end
axislegend(ax1; position = :lb, labelsize = 9)
for (i, r) in enumerate(sort(tr, by = r -> r.flux))
    s = results[r.name]
    ax = Axis(fig[2, i]; title = "$(r.name)", xlabel = "time (s)", ylabel = "T (°C)")
    for (k, d) in enumerate((0.002, 0.010, 0.020, 0.040))
        j = findfirst(x -> isapprox(x, d), DEPTHS)
        js = findfirst(x -> isapprox(x, d), s.depths)
        c = Makie.wong_colors()[k]
        live(r, j) && lines!(
            ax,
            r.t,
            usable(r, j) .- 273.15;
            color = c,
            label = "$(round(Int, 1000d)) mm",
        )
        js === nothing ||
            lines!(ax, s.t, s.T[:, js] .- 273.15; color = c, linestyle = :dash)
    end
    axislegend(ax; position = :lt, labelsize = 9)
end
save(joinpath(RESULTS, "figures", "smoke_$species.png"), fig)
println("wrote smoke_$species.png")
