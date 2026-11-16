# The front face of the decks, calibrated on the chamber's calcium-silicate runs
# (decks/silcal.csv, Th Fig. 5.5) as the thesis did, but in 1D. Fitted, for each
# of two face models: one h per flux and one factor on the hot-wire
# conductivity of the calcium silicate, shared by the three fluxes (the thesis
# also fitted the silicate's k, 11.6 % from the hot wire), to the thermocouples
# at 2–40 mm over 10–600 s.
# - gray (`front = :silcal`): the face absorbs the nominal flux at a gray
#   absorptivity and re-radiates to T_a, so h closes the balance with that
#   sink, the window's return included. The silicate's absorptivity is not
#   published and moves h one for one; h is refitted at 0.80 and 0.92.
# - window (`front = :window`): the cone, the sapphire window and the face's
#   exchange with it from first principles (chamber_radiation.jl), the
#   silicate's two-band absorptivity; h is the convection left to fit. With h
#   at natural convection (6 W/(m² K)) nothing is fitted but k: that run is
#   reported too.
# Writes front.toml (read by `problem(...; front)`) and
# results/figures/front_calibration.png.
#
#   julia --project=scripts/verification scripts/lemta_argon/calibrate_front.jl
include(joinpath(@__DIR__, "LemtaArgon.jl"))
using .LemtaArgon
using CairoMakie
using Printf
using Pyrolysis

const L = LemtaArgon
const DATA = silcal_data()
const CH = chamber()
const GEOMETRY = L.chamber_geometry(CH)
const FLUXES = (25, 35, 45)
const SCORED = 0.002 .<= DATA[25].depths .<= 0.04   # the rear face is the boundary
const TIMES = 10.0:10.0:600.0
const E = CH["silcal"]["emissivity"]
ε_silcal(T) = (Tc = min(T, 1300.0); E[1] + E[2] * Tc + E[3] * Tc^2)

# temperature at (t, depth) of a solved silicate block
function evaluator(sol, nodes)
    string(sol.retcode) == "Success" || error("silicate solve: $(sol.retcode)")
    centres = 0.5 .* (nodes[1:(end-1)] .+ nodes[2:end])
    return (t, z) ->
        interp(centres, view(sol.T, :, findfirst(>=(t - 1e-9), sol.t)), nodes[end] - z)
end

function gray(flux; h, k_factor, absorptivity = CH["silcal"]["absorptivity"])
    prob, nodes = silcal_problem(flux; h, k_factor, absorptivity, ch = CH, data = DATA)
    return evaluator(solve(prob; reltol = 1e-6, abstol = 1e-6, dtmax = 5.0), nodes)
end

function window(flux; h, k_factor, passes = 2)
    q = 1000.0 * flux
    T_c = L.cone_temperature(q, GEOMETRY)
    α = L.absorptivity_behind_window(L.silcal_spectrum(CH), T_c, GEOMETRY.T_room)
    d = DATA[flux]
    T_a = t -> interp(d.t, d.T[:, end], t)
    build(s) = silcal_problem(
        flux;
        h,
        k_factor,
        absorptivity = α,
        ch = CH,
        data = DATA,
        series = s,
    )
    series = L.window_series(build, q, T_a, ε_silcal, 600.0, GEOMETRY; passes)
    prob, nodes = build(series)
    return evaluator(solve(prob; reltol = 1e-6, abstol = 1e-6, dtmax = 5.0), nodes),
    series,
    α
end

# rms [K] over the scored thermocouples and times, and the mean bias by depth
function score(flux, m)
    d = DATA[flux]
    err = Float64[]
    bias = Dict{Int,Float64}()
    for (j, z) in enumerate(d.depths)
        SCORED[j] || continue
        e = [m(t, z) - interp(d.t, view(d.T, :, j), t) for t in TIMES]
        append!(err, e)
        bias[round(Int, 1000z)] = sum(e) / length(e)
    end
    return sqrt(sum(abs2, err) / length(err)), bias
end

function golden(f, a, b; n = 14)
    g = (sqrt(5) - 1) / 2
    c, d = b - g * (b - a), a + g * (b - a)
    fc, fd = f(c), f(d)
    for _ = 1:n
        if fc < fd
            b, d, fd = d, c, fc
            c = b - g * (b - a)
            fc = f(c)
        else
            a, c, fc = c, d, fd
            d = a + g * (b - a)
            fd = f(d)
        end
    end
    return fc < fd ? c : d
end

# h per flux at a conductivity factor, for a face model `m(flux; h, k_factor)`
# returning an evaluator
function fit_h(m, k_factor)
    out = Dict{Int,Tuple{Float64,Float64,Dict{Int,Float64}}}()
    for flux in FLUXES
        h = golden(h -> score(flux, m(flux; h, k_factor))[1], 0.0, 40.0)
        r, b = score(flux, m(flux; h, k_factor))
        out[flux] = (h, r, b)
    end
    return out
end
total(fit) = sqrt(sum(v[2]^2 for v in values(fit)) / length(fit))
hs(fit) = join([@sprintf("%.1f", fit[f][1]) for f in FLUXES], " / ")

function best_fit(m, factors, tag)
    println("$tag: conductivity factor on the silicate's hot-wire k")
    best = nothing
    for kf in factors
        fit = fit_h(m, kf)
        @printf("  k × %.3f: rms %.2f K; h %s\n", kf, total(fit), hs(fit))
        flush(stdout)
        (best === nothing || total(fit) < total(best[2])) && (best = (kf, fit))
    end
    kf, fit = best
    println("  best: k × $kf")
    for f in FLUXES
        h, r, b = fit[f]
        @printf(
            "  %d kW/m²: h %.2f W/(m² K), rms %.2f K; bias %s K\n",
            f,
            h,
            r,
            join(
                ["$(d) mm $(@sprintf("%+.1f", b[d]))" for d in sort(collect(keys(b)))],
                ", ",
            )
        )
    end
    return kf, fit
end

kf_gray, fit_gray = best_fit(gray, 0.70:0.025:0.95, "gray face")
α_gray = CH["silcal"]["absorptivity"]
alt = Dict(
    a => fit_h((f; h, k_factor) -> gray(f; h, k_factor, absorptivity = a), kf_gray) for
    a in (0.80, 0.92)
)
for a in (0.80, 0.92)
    @printf("  absorptivity %.2f: h %s, rms %.2f K\n", a, hs(alt[a]), total(alt[a]))
end

window_eval(f; h, k_factor) = window(f; h, k_factor)[1]
kf_win, fit_win = best_fit(window_eval, 0.70:0.05:0.90, "cone and window")
natural = Dict{Int,Tuple{Float64,Float64,Dict{Int,Float64}}}()
for f in FLUXES
    r, b = score(f, window_eval(f; h = GEOMETRY.h_window, k_factor = kf_win))
    natural[f] = (GEOMETRY.h_window, r, b)
end
@printf(
    "  natural convection (h %.1f) at k × %.3f: rms %s K\n",
    GEOMETRY.h_window,
    kf_win,
    join([@sprintf("%.1f", natural[f][2]) for f in FLUXES], " / ")
)
win_detail = Dict(f => window(f; h = fit_win[f][1], k_factor = kf_win) for f in FLUXES)
for f in FLUXES
    _, s, α = win_detail[f]
    @printf(
        "  %d kW/m²: coil %.0f °C, silicate absorptivity %.3f, window %.0f / %.0f °C at 300 / 600 s, long-wave irradiation %.1f kW/m² at 600 s\n",
        f,
        s.T_c - 273.15,
        α,
        interp(s.t, s.T_w, 300.0) - 273.15,
        s.T_w[end] - 273.15,
        s.G[end] / 1e3
    )
end

open(FRONT, "w") do io
    println(
        io,
        """
# Front-face convection [W/(m² K)] of the decks, by nominal flux [kW/m²],
# written by calibrate_front.jl: the chamber's calcium-silicate runs
# (decks/silcal.csv, Th Fig. 5.5) fitted in 1D, one h per flux and one factor
# on the silicate's hot-wire conductivity shared by the three fluxes, over the
# thermocouples at 2–40 mm, 10–600 s.
#
# Gray face (front = :silcal): the nominal flux absorbed at the silicate's
# absorptivity below, re-radiation to T_a; h moves one for one with that
# absorptivity (h_at_absorptivity_*).""",
    )
    @printf(io, "absorptivity_silcal = %.2f\n", α_gray)
    @printf(io, "k_factor_silcal = %.3f\n", kf_gray)
    @printf(io, "rms_K = %.2f\n", total(fit_gray))
    println(io, "\n[h_front]")
    for f in FLUXES
        @printf(io, "%d = %.2f\n", f, fit_gray[f][1])
    end
    println(io, "\n[rms_by_flux_K]")
    for f in FLUXES
        @printf(io, "%d = %.2f\n", f, fit_gray[f][2])
    end
    for a in (0.80, 0.92)
        @printf(io, "\n[h_at_absorptivity_%d]\n", round(Int, 100a))
        for f in FLUXES
            @printf(io, "%d = %.2f\n", f, alt[a][f][1])
        end
    end
    println(
        io,
        """

# Under the cone and window (front = :window, chamber_radiation.jl): h is the
# convection left to fit once the window's radiation is computed.""",
    )
    println(io, "[window]")
    @printf(io, "k_factor_silcal = %.3f\n", kf_win)
    @printf(io, "rms_K = %.2f\n", total(fit_win))
    @printf(
        io,
        "rms_natural_convection_K = %.2f   # h = %.1f, nothing fitted but k\n",
        total(natural),
        GEOMETRY.h_window
    )
    println(io, "\n[window.h_front]")
    for f in FLUXES
        @printf(io, "%d = %.2f\n", f, fit_win[f][1])
    end
    println(io, "\n[window.rms_by_flux_K]")
    for f in FLUXES
        @printf(io, "%d = %.2f\n", f, fit_win[f][2])
    end
end
println("wrote ", FRONT)

fig = Figure(size = (1500, 850))
for (row, (tag, m, kf, fit)) in enumerate((
    ("gray face", gray, kf_gray, fit_gray),
    ("cone and window", window_eval, kf_win, fit_win),
))
    for (i, f) in enumerate(FLUXES)
        ev = m(f; h = fit[f][1], k_factor = kf)
        d = DATA[f]
        ax = Axis(
            fig[row, i];
            title = "$tag, $f kW/m² (h $(round(fit[f][1], digits = 1)))",
            xlabel = "time (s)",
            ylabel = "T (°C)",
        )
        cmap = cgrad(:viridis, length(d.depths), categorical = true)
        for (j, z) in enumerate(d.depths)
            lines!(
                ax,
                d.t,
                d.T[:, j] .- 273.15;
                color = cmap[j],
                label = "$(round(Int, 1000z)) mm",
            )
            SCORED[j] && lines!(
                ax,
                d.t,
                [ev(t, z) - 273.15 for t in d.t];
                color = cmap[j],
                linestyle = :dash,
            )
        end
        row == 1 && i == 1 && axislegend(ax; position = :lt, nbanks = 2, labelsize = 8)
    end
end
mkpath(joinpath(LemtaArgon.RESULTS, "figures"))
save(joinpath(LemtaArgon.RESULTS, "figures", "front_calibration.png"), fig)
println("wrote front_calibration.png")
