# Step 1: load every run, apply the defect list, and plot the temperature at
# every depth and the mass against time, one figure per species
# (results/figures/data_<species>.png), with a table of the runs
# (results/data_summary.csv). With --along the along-grain runs
# (data_<species>_along.png, data_summary_along.csv), and the withdrawal read
# from the mass checked against the thermocouples' where a run has both.
#
#   julia --project=scripts/verification scripts/lemta_moist_wood/check_data.jl [--along]
include(joinpath(@__DIR__, "LemtaMoistWood.jl"))
using .LemtaMoistWood
using CairoMakie
using Printf

grain = "--along" in ARGS ? :along : :across
tag = grain === :along ? "_along" : ""
check_data(; grain)
runs = load_runs(; grain, excluded = true)
figdir = joinpath(RESULTS, "figures")
mkpath(figdir)

open(joinpath(RESULTS, "data_summary$tag.csv"), "w") do io
    println(
        io,
        "run,species,grain,mc,rho_dry,rho_from_mass,thermocouples,duration_h,t_off_s,T_amb_C,m0_g,water_g,loss_2p5h_g,loss_end_g,water_lost_2p5h,defect",
    )
    for r in runs
        t25 = min(9000.0, r.t[end])
        @printf(
            io,
            "\"%s\",%s,%s,%.4f,%.0f,%.0f,%d,%.2f,%.0f,%.1f,%.2f,%.2f,%.2f,%.2f,%.3f,%s\n",
            r.name,
            r.species,
            r.grain,
            r.mc,
            r.rho,
            isempty(r.mass) ? NaN : mass_density(r),
            length(r.depths),
            r.t[end] / 3600,
            r.t_off,
            r.T_amb - 273.15,
            isempty(r.mass) ? NaN : r.mass[1],
            measured_water_g(r),
            mass_loss(r, t25),
            isempty(r.mass) ? NaN : mass_loss(r, r.t_mass[end]),
            water_fraction_lost(r, t25),
            r.defect
        )
    end
end

if grain === :along
    println("withdrawal: thermocouples / mass rate (h)")
    for r in runs
        isempty(r.mass) && continue
        m = LemtaMoistWood.cone_off_mass(r.t_mass, r.mass)
        @printf(
            "  %-18s %s / %.3f\n",
            r.name,
            isempty(r.depths) ? "  –  " : @sprintf("%.3f", r.t_off / 3600),
            m / 3600
        )
    end
end

for species in unique(r.species for r in runs)
    rs = filter(r -> r.species == species, runs)
    shown = filter(r -> !isempty(r.depths), rs)
    n = length(shown)
    fig = Figure(size = (1100, 260 * (n + 1)))
    for (i, r) in enumerate(shown)
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
        lines!(ax, r.t_mass ./ 3600, r.mass[1] .- r.mass; label = r.name)
    end
    axislegend(ax; position = :rb, labelsize = 9)
    save(joinpath(figdir, "data_$species$tag.png"), fig)
    println("wrote data_$species$tag.png")
end
