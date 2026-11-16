# Step 1: load every run, apply the defect list, and plot the temperature at
# every depth and the mass against time, one figure per species
# (results/figures/data_<species>.png), with a table of the runs
# (results/data_summary.csv).
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/check_data.jl
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using CairoMakie
using Printf

check_data()
runs = load_runs(excluded = true)
figdir = joinpath(RESULTS, "figures")
mkpath(figdir)

open(joinpath(RESULTS, "data_summary.csv"), "w") do io
    println(
        io,
        "run,species,mc,rho_dry,thermocouples,duration_h,t_off_s,T_amb_C,m0_g,water_g,loss_2p5h_g,loss_end_g,defect",
    )
    for r in runs
        t25 = min(9000.0, r.t[end])
        @printf(
            io,
            "\"%s\",%s,%.4f,%.0f,%d,%.2f,%.0f,%.1f,%.1f,%.1f,%.2f,%.2f,%s\n",
            r.name,
            r.species,
            r.mc,
            r.rho,
            length(r.depths),
            r.t[end] / 3600,
            r.t_off,
            r.T_amb - 273.15,
            isempty(r.mass) ? NaN : r.mass[1],
            initial_water_g(r),
            mass_loss(r, t25),
            mass_loss(r, r.t[end]),
            r.defect
        )
    end
end

for species in unique(r.species for r in runs)
    rs = filter(r -> r.species == species, runs)
    n = length(rs)
    fig = Figure(size = (1100, 260 * (n + 1)))
    for (i, r) in enumerate(rs)
        title = "$(r.name)" * (r.defect == :none ? "" : "  [$(r.defect)]")
        ax = Axis(fig[i, 1:2]; title, xlabel = "time (h)", ylabel = "T (°C)")
        cmap = cgrad(:viridis, length(r.depths), categorical = true)
        for (j, d) in enumerate(r.depths)
            lines!(
                ax,
                r.t ./ 3600,
                r.T[:, j] .- 273.15;
                color = cmap[j],
                label = "$(round(Int, 1000d)) mm",
            )
        end
        isfinite(r.t_off) && vlines!(ax, [r.t_off / 3600]; color = :gray, linestyle = :dash)
        i == 1 && axislegend(ax; position = :rb, nbanks = 3, labelsize = 9)
    end
    ax = Axis(
        fig[n+1, 1:2];
        title = "$species: mass lost",
        xlabel = "time (h)",
        ylabel = "m₀ − m (g)",
    )
    for r in rs
        isempty(r.mass) && continue
        lines!(ax, r.t ./ 3600, r.mass[1] .- r.mass; label = r.name)
    end
    axislegend(ax; position = :rb, labelsize = 9)
    save(joinpath(figdir, "data_$species.png"), fig)
    println("wrote data_$species.png")
end
