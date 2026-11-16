# Check the dataset against its pinned checksums, load every file with the
# defect list applied, and draw it: per species, the temperature at every
# depth (cut parts dashed and grey) and the mass at the three fluxes
# (results/figures/data_<species>.png); the runs in results/data_summary.csv.
#
#   julia --project=scripts/verification scripts/lemta_argon/check_data.jl
include(joinpath(@__DIR__, "LemtaArgon.jl"))
using .LemtaArgon
using CairoMakie
using Printf

check_data()
ms = mass_runs()
ts = temperature_runs()
figdir = joinpath(LemtaArgon.RESULTS, "figures")
mkpath(figdir)

open(joinpath(LemtaArgon.RESULTS, "data_summary.csv"), "w") do io
    println(
        io,
        "file,species,flux_kW_m2,rho_dry,duration_s,t_end_s,m_end_over_m0,live_TCs,T0_2mm_C,T0_40mm_C,notes",
    )
    for r in ms
        @printf(
            io,
            "%s,%s,%g,%g,%g,%g,%.4f,,,,\n",
            r.name,
            r.species,
            r.flux,
            r.rho,
            r.t[end],
            r.t[end],
            r.m[end]
        )
    end
    for r in ts
        @printf(
            io,
            "%s,%s,%g,%g,%g,%g,,%d,%.1f,%.1f,\"%s\"\n",
            r.name,
            r.species,
            r.flux,
            r.rho,
            r.t[end],
            r.t_end,
            count(j -> live(r, j), eachindex(DEPTHS)),
            r.T0[1] - 273.15,
            r.T0[end] - 273.15,
            join(r.notes, "; ")
        )
    end
end

for species in sort(unique(r.species for r in ms))
    tr = filter(r -> r.species == species, ts)
    mr = filter(r -> r.species == species, ms)
    fig = Figure(size = (1200, 330 * (length(tr) + 1)))
    for (i, r) in enumerate(tr)
        ax = Axis(fig[i, 1:2]; title = r.name, xlabel = "time (s)", ylabel = "T (°C)")
        cmap = cgrad(:viridis, length(DEPTHS), categorical = true)
        for j in eachindex(DEPTHS)
            live(r, j) || continue
            lines!(ax, r.t, r.T[:, j] .- 273.15; color = (:gray, 0.6), linestyle = :dash)
            lines!(
                ax,
                r.t,
                usable(r, j) .- 273.15;
                color = cmap[j],
                label = "$(round(Int, 1000DEPTHS[j])) mm",
            )
        end
        r.t_end < r.t[end] && vlines!(ax, [r.t_end]; color = :red, linestyle = :dot)
        axislegend(ax; position = :lt, nbanks = 4, labelsize = 9)
        ax2 = Axis(fig[i, 3]; title = "first 120 s", xlabel = "time (s)")
        for j in eachindex(DEPTHS)
            live(r, j) || continue
            k = r.t .<= 120
            lines!(ax2, r.t[k], r.T[k, j] .- 273.15; color = cmap[j])
        end
    end
    n = length(tr) + 1
    ax = Axis(fig[n, 1]; title = "$species: m/m₀", xlabel = "time (s)")
    ax2 = Axis(fig[n, 2]; title = "MLRPUA (g/(m² s))", xlabel = "time (s)")
    for r in mr
        lab = "$(round(Int, r.flux)) kW/m², ρ $(r.rho)"
        lines!(ax, r.t, r.m; label = lab)
        lines!(ax2, r.t, r.mlr; label = lab)
    end
    axislegend(ax; position = :lb, labelsize = 9)
    save(joinpath(figdir, "data_$species.png"), fig)
    println("wrote data_$species.png")
end
