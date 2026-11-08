# Cost at matched accuracy and at practitioner settings, from a run_timing.py
# campaign.
#
#   julia --project=scripts/macfp_pmma scripts/macfp_pmma/timing_report.jl OUTDIR
#
# Every code is scored against its own finest run of the case: Pyrolysis.jl
# against its 80 cells/mm, reltol 1e-8 reference (time_pyrolysis.jl); FDS
# against its half-cell run; Gpyro and ThermaKin against their half-cell and
# half-step run (`cost_rows` in MacfpPMMA.jl).
include(joinpath(@__DIR__, "MacfpPMMA.jl"))
using .MacfpPMMA, CairoMakie, Printf, Statistics

const RES = RESULTS
const COLOUR = Dict(
    "Pyrolysis.jl" => colorant"#2a78d6",
    "fds" => colorant"#eb6834",
    "gpyro" => colorant"#1baf7a",
    "thermakin" => colorant"#e87ba4",
)
const NAME = Dict(
    "Pyrolysis.jl" => "Pyrolysis.jl",
    "fds" => "FDS 6.11.1",
    "gpyro" => "Gpyro V2",
    "thermakin" => "ThermaKin2Ds",
)

function main(out)
    rows = cost_rows(out)   # (code, set, case, setting, time, e)
    # Cheapest within tolerance, and the practitioner setting, per code and case.
    open(joinpath(out, "cost_summary.md"), "w") do io
        println(
            io,
            "| set | case | code | cheapest within tolerance (s) | its setting | practitioner setting (s) | practitioner, span-scaled to burnout (s) |",
        )
        println(io, "|---|---|---|---|---|---|---|")
        for key in unique((r[2], r[3]) for r in rows),
            code in ("Pyrolysis.jl", "fds", "gpyro", "thermakin")

            rs = filter(r -> (r[2], r[3]) == key && r[1] == code, rows)
            isempty(rs) && continue
            ok = filter(r -> r[6] <= 1, rs)
            best = isempty(ok) ? nothing : ok[argmin([r[5] for r in ok])]
            prac = filter(r -> r[4] in ("20/mm ×1", "base"), rs)
            ps = isempty(prac) ? NaN : span_scaled(code, key[1], key[2], "base", prac[1][5])
            @printf(
                io,
                "| %s | %s | %s | %s | %s | %s | %s |\n",
                key...,
                NAME[code],
                best === nothing ? "not reached" : @sprintf("%.2f", best[5]),
                best === nothing ? "—" : best[4],
                isempty(prac) ? "—" : @sprintf("%.2f", prac[1][5]),
                isnan(ps) ? "—" : @sprintf("%.2f", ps)
            )
        end
    end
    # Time against combined error, one panel per external-subset case. Marker
    # shape: Pyrolysis.jl cells per mm (filled at 10 × the default tolerances,
    # open at the default); external codes, base (circle) or refined (square).
    keys_ = unique((r[2], r[3]) for r in rows if r[1] != "Pyrolysis.jl")
    n = length(keys_)
    n == 0 && return
    ncol = 4
    CODES = ("Pyrolysis.jl", "fds", "gpyro", "thermakin")
    SHAPE = Dict("10" => :circle, "20" => :rect, "40" => :diamond)
    marker(r) =
        r[1] == "Pyrolysis.jl" ? SHAPE[first(split(r[4], "/"))] :
        r[4] == "base" ? :circle : :rect
    filled(r) = r[1] != "Pyrolysis.jl" || endswith(r[4], "×10")
    xt = ([1e-3, 1e-2, 1e-1, 1, 10], ["0.001", "0.01", "0.1", "1", "10"])
    yt = ([0.1, 0.3, 1, 3, 10, 30, 100], ["0.1", "0.3", "1", "3", "10", "30", "100"])
    fig = Figure(
        size = (1900, 400 * cld(n, ncol) + 110),
        fontsize = 14,
        backgroundcolor = :white,
    )
    for (k, key) in enumerate(keys_)
        ax = Axis(
            fig[cld(k, ncol)+1, mod1(k, ncol)];
            title = "$(key[1]), $(replace(key[2], "gasification_" => ""))",
            xscale = log10,
            yscale = log10,
            xticks = xt,
            yticks = yt,
            limits = ((5e-4, 30), (0.1, 100)),
            xlabel = "error / tolerance",
            ylabel = "solve time (s)",
            xgridcolor = colorant"#e4e3df",
            ygridcolor = colorant"#e4e3df",
        )
        vlines!(ax, [1.0]; color = :gray50, linestyle = :dash)
        for code in CODES,
            r in filter(r -> (r[2], r[3]) == key && r[1] == code && isfinite(r[6]), rows)

            c = COLOUR[code]
            scatter!(
                ax,
                [max(r[6], 1e-3)],
                [r[5]];
                marker = marker(r),
                markersize = 12,
                color = filled(r) ? c : :white,
                strokecolor = c,
                strokewidth = 2,
            )
        end
    end
    Legend(
        fig[1, 1:ncol],
        [
            [
                MarkerElement(color = COLOUR[c], marker = :circle, markersize = 12) for
                c in CODES
            ];
            [
                MarkerElement(color = :gray30, marker = SHAPE[s], markersize = 12) for
                s in ("10", "20", "40")
            ];
            [
                MarkerElement(
                    color = :white,
                    strokecolor = :gray30,
                    strokewidth = 2,
                    marker = :circle,
                    markersize = 12,
                ),
            ]
        ],
        [
            [NAME[c] for c in CODES];
            [
                "10 cells/mm (external: base)",
                "20 cells/mm (external: refined)",
                "40 cells/mm",
            ];
            ["Pyrolysis.jl at default tolerances (filled: × 10)"]
        ];
        orientation = :horizontal,
        framevisible = false,
    )
    save(joinpath(RES, "figures", "cost_vs_error.png"), fig)
    println("wrote $(joinpath(out, "cost_summary.md")) and figures/cost_vs_error.png")
end

main(ARGS[1])
