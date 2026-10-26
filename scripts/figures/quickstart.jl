# Execute the canonical tutorial itself, so its deck cannot silently drift from
# the plotted example. Only the optional Plots import and progress UI are omitted.
function quickstart_data()
    guide = read(joinpath(@__DIR__, "..", "..", "docs", "User_Guide.md"), String)
    section =
        split(split(guide, "## 2.10 Complete script"; limit = 2)[2], "## 2.11"; limit = 2)[1]
    code = split(split(section, "```julia\n"; limit = 2)[2], "```"; limit = 2)[1]
    code = first(split(code, "# ---- Step 6:"; limit = 2))
    code = replace(code, "using Plots\n" => "", "progress = true" => "progress = false")
    tutorial = Module(:GuideFigureTutorial)
    sol = Base.include_string(tutorial, code * "\nsolution\n", "User_Guide_2_10.jl")
    @assert sol.retcode == :Success
    @assert sol.t[end] == 1600.0 && length(sol.t) == 601
    @assert all(isfinite, sol.surface_T) && all(isfinite, sol.mass_loss_rate)
    @assert maximum(sol.mass_loss_rate) > 0
    mkpath(joinpath(OUT, "data"))
    open(joinpath(OUT, "data", "quickstart_timeseries.csv"), "w") do io
        println(
            io,
            "time_s,surface_K,back_K,mlr_kg_m2_s,gas_generation_kg_m2_s,thickness_m",
        )
        for k in eachindex(sol.t)
            println(
                io,
                join(
                    (
                        sol.t[k],
                        sol.surface_T[k],
                        sol.extras.back_temperature[k],
                        sol.mass_loss_rate[k],
                        sol.extras.gas_generation_rate[k],
                        sol.extras.thickness_history[k],
                    ),
                    ',',
                ),
            )
        end
    end
    selected = [0.0, 200.0, 400.0, 800.0, 1200.0, 1600.0]
    indices = [argmin(abs.(sol.t .- t)) for t in selected]
    open(joinpath(OUT, "data", "quickstart_profiles.csv"), "w") do io
        println(io, "time_s,cell_id,current_z_m,temperature_K")
        for k in indices, i in axes(sol.T, 1)
            isfinite(sol.T[i, k]) || continue
            println(io, join((sol.t[k], i, sol.cell_centers[i, k], sol.T[i, k]), ','))
        end
    end
    open(joinpath(OUT, "data", "quickstart_provenance.toml"), "w") do io
        println(io, "julia_version = \"", VERSION, "\"")
        println(io, "tutorial = \"docs/User_Guide.md section 2.10\"")
        println(io, "executed_deck_sha256 = \"", bytes2hex(SHA.sha256(code)), "\"")
        root = normpath(joinpath(@__DIR__, "..", ".."))
        revision = readchomp(Cmd(["git", "-C", root, "rev-parse", "HEAD"]))
        println(io, "source_revision = \"", revision, "\"")
        println(io, "retcode = \"", sol.retcode, "\"")
        println(io, "outcome = \"", sol.outcome, "\"")
        println(io, "cells = ", size(sol.T, 1))
        println(io, "saved_times = ", length(sol.t))
        println(io, "final_time_s = ", sol.t[end])
        println(io, "peak_mlr_kg_m2_s = ", maximum(sol.mass_loss_rate))
        println(io, "peak_surface_K = ", maximum(sol.surface_T))
    end
    return sol, indices
end

# Redraw from committed data while adjusting presentation, without another solve.
function saved_quickstart_data()
    readcsv(path) = [parse.(Float64, split(row, ',')) for row in readlines(path)[2:end]]
    rows = readcsv(joinpath(OUT, "data", "quickstart_timeseries.csv"))
    data = reduce(hcat, rows)
    profiles = readcsv(joinpath(OUT, "data", "quickstart_profiles.csv"))
    n = Int(maximum(p[2] for p in profiles))
    T = fill(NaN, n, size(data, 2))
    centres = copy(T)
    indices = Int[]
    for p in profiles
        k = findfirst(==(p[1]), data[1, :])
        i = Int(p[2])
        T[i, k], centres[i, k] = p[4], p[3]
        k in indices || push!(indices, k)
    end
    sol = (;
        t = data[1, :],
        surface_T = data[2, :],
        mass_loss_rate = data[4, :],
        T,
        cell_centers = centres,
        extras = (;
            back_temperature = data[3, :],
            gas_generation_rate = data[5, :],
            thickness_history = data[6, :],
        ),
    )
    return sol, indices
end

function quickstart_figures(sol, indices)
    ticks = [(t, string(t)) for t = 0:400:1600]
    tticks = [(t, string(t)) for t = 300:100:700]
    colors = [MUTED, BLUE, TEAL, ORANGE, PURPLE, "#634641"]
    svg(
        "08_quickstart",
        "Figure 8. PMMA Cone-Calorimeter Simulation",
        "Tutorial section 2.10, 72 cells, 24.5325 kilowatts per square metre incident radiation, solved to 1600 seconds. Surface and back temperatures, actual boundary mass-loss rate and gas-generation proxy, and temperature profiles at current cell positions.";
        height = 840,
    ) do io
        label(io, 110, 90, "(a) Boundary Temperatures"; weight = "bold")
        X, Y = chart_axes(
            io,
            110,
            125,
            440,
            225;
            xlim = (0, 1600),
            ylim = (280, 720),
            xticks = ticks,
            yticks = tticks,
            xlabel = "Time t [s]",
            ylabel = "Temperature T [K]",
        )
        polyline(io, X.(sol.t), Y.(sol.surface_T); color = ORANGE)
        polyline(io, X.(sol.t), Y.(sol.extras.back_temperature); color = BLUE, dash = "9 5")
        key(io, 350, 300, "Surface", ORANGE)
        key(io, 350, 330, "Back Face", BLUE; dash = "9 5")
        label(io, 700, 90, "(b) Mass-Loss Rate"; weight = "bold")
        X, Y = chart_axes(
            io,
            700,
            125,
            440,
            225;
            xlim = (0, 1600),
            ylim = (0, 6),
            xticks = ticks,
            yticks = [(i, string(i)) for i = 0:2:6],
            xlabel = "Time t [s]",
            ylabel = "Mass Rate [g/(m²·s)]",
        )
        polyline(
            io,
            X.(sol.t),
            Y.(1000 .* sol.extras.gas_generation_rate);
            color = MUTED,
            width = 4,
            dash = "8 5",
        )
        polyline(io, X.(sol.t), Y.(1000 .* sol.mass_loss_rate); color = TEAL, width = 2.5)
        key(io, 850, 300, "Boundary MLR", TEAL)
        key(io, 850, 330, "Gas-Generation Proxy", MUTED; dash = "8 5")
        label(io, 110, 465, "(c) Temperature Profiles"; weight = "bold")
        X, Y = chart_axes(
            io,
            110,
            500,
            650,
            265;
            xlim = (0, 5.8),
            ylim = (280, 720),
            xticks = [(i, string(i)) for i = 0:1:5],
            yticks = tticks,
            xlabel = "Current Position z [mm]",
            ylabel = "Temperature T [K]",
        )
        for (k, color) in zip(indices, colors)
            polyline(io, X.(1000 .* sol.cell_centers[:, k]), Y.(sol.T[:, k]); color)
            # Surface is a face temperature at the current outer boundary.
            circle(
                io,
                X(1000sol.extras.thickness_history[k]),
                Y(sol.surface_T[k]),
                4;
                fill = color,
            )
        end
        for (j, k) in enumerate(indices)
            key(io, 830, 523 + 33(j-1), "$(Int(round(sol.t[k]))) s", colors[j])
        end
        circle(io, 836, 745, 4; fill = MUTED)
        label(io, 855, 751, "Surface Face"; size = 19)
    end
    svg(
        "08_quickstart_compact",
        "PMMA Cone-Calorimeter Simulation",
        "Compact view of the 1600-second User Guide tutorial: surface and back temperatures, and actual boundary mass-loss rate. This is the full tutorial deck, not the shorter home-page example.";
        height = 405,
    ) do io
        X, Y = chart_axes(
            io,
            105,
            105,
            440,
            220;
            xlim = (0, 1600),
            ylim = (280, 720),
            xticks = ticks,
            yticks = tticks,
            xlabel = "Time t [s]",
            ylabel = "Temperature T [K]",
        )
        polyline(io, X.(sol.t), Y.(sol.surface_T); color = ORANGE)
        polyline(io, X.(sol.t), Y.(sol.extras.back_temperature); color = BLUE, dash = "9 5")
        key(io, 350, 275, "Surface", ORANGE)
        key(io, 350, 305, "Back Face", BLUE; dash = "9 5")
        X, Y = chart_axes(
            io,
            705,
            105,
            440,
            220;
            xlim = (0, 1600),
            ylim = (0, 6),
            xticks = ticks,
            yticks = [(i, string(i)) for i = 0:2:6],
            xlabel = "Time t [s]",
            ylabel = "MLR [g/(m²·s)]",
        )
        polyline(io, X.(sol.t), Y.(1000 .* sol.mass_loss_rate); color = TEAL)
        key(io, 820, 300, "Boundary Mass-Loss Rate", TEAL)
    end
end
