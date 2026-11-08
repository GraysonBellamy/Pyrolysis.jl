# Pyrolysis.jl on the MaCFP-2 PMMA targets: every deck × gasification case
# (and its variants) and every deck × TGA rate.
#
#   julia --project=scripts/macfp_pmma scripts/macfp_pmma/run_pyrolysis.jl [SET ...] [--tga-only] [--cells=20]
#
# Writes results/pyrolysis/<set>/<variant>/{gasification_<q>kW_<L>mm,tga_<β>K}.csv
# and a run log per set (results/pyrolysis/<set>/runs.csv).
include(joinpath(@__DIR__, "MacfpPMMA.jl"))
using .MacfpPMMA, Pyrolysis, Printf

const OUT = joinpath(@__DIR__, "results", "pyrolysis")

# FDS's two-flux sweep crosses a cell of width Δz with transmission
# 1/(1 + 2κΔz); the diffuse two-stream model crosses it with e^(−1.5κΔz). The
# two agree cell by cell, in transmission and in emission, at an absorption
# coefficient of ln(1 + 2κΔz)/(1.5Δz), taken here at FDS's 50 µm cells.
const FDS_CELL = 5e-5

function fds_closure(deck; dz = FDS_CELL)
    d = deepcopy(deck)
    for c in d["components"]
        haskey(c, "kappa") && (c["kappa"] = log1p(2 * Float64(c["kappa"]) * dz) / (1.5dz))
    end
    return d
end

# The one item of a deck that FDS cannot express, aligned with FDS on top of its
# cell attenuation; `nothing` where FDS expresses the whole deck.
function fds_matched(deck)
    id = deck["id"]
    d = fds_closure(deck)
    if id in ("UMD", "UMET_TK")
        # FDS removes the pyrolyzate where it forms: no gas enthalpy in the solid.
        d["gas"]["cp"] = 0.0
    elseif id == "UMET_GP"
        # Reaction 2 (B → gas, n = 1.8) is first order in FDS: on a shrinking cell
        # of pure B the partial density stays at ρ_B, so its rate per unit mass is
        # A (ρ_B/M_ref)^(n−1), with M_ref = 0.962 ρ₀ the inventory of B that
        # reaction 1 forms.
        r = d["reactions"][2]
        ρB = Float64(only(filter(c -> c["id"] == r["reactant"], d["components"]))["rho"])
        Mref =
            Float64(d["reactions"][1]["residue_yield"]) *
            MacfpPMMA._initial_condensed_density(d)
        r["A"] = Float64(r["A"]) * (ρB / Mref)^(Float64(r["n"]) - 1)
        r["n"] = 1.0
    else
        return nothing
    end
    return d
end

as_surface(deck) = merge(deck, Dict("absorption" => "surface"))

# The deck as written ("base", diffuse heater) and its variants. Each variant
# moves one choice, or reproduces what a group evidently ran.
function variants(deck)
    base = (deck = deck, h = 0.0, T_amb = MacfpPMMA.T_AMBIENT, distribution = Diffuse())
    v = Pair{String,Any}["base"=>base]
    id = deck["id"]
    if deck["absorption"] == "in-depth"
        # The heater as one beam along the normal, absorbed at κ rather than at
        # the diffuse 1.5κ.
        push!(v, "collimated" => merge(base, (distribution = Collimated(),)))
        # FDS's attenuation on its own cells, every component's κ scaled.
        push!(v, "fds_closure" => merge(base, (deck = fds_closure(deck),)))
        # FDS's attenuation and the item FDS cannot express aligned with FDS.
        m = fds_matched(deck)
        m === nothing || push!(v, "fds_matched" => merge(base, (deck = m,)))
        # No radiation in depth: the heater absorbed at the face, which alone
        # re-radiates.
        push!(v, "surface" => merge(base, (deck = as_surface(deck),)))
    end
    if id in ("BUW-FZJ_A", "BUW-FZJ_B")
        # Their runs end in a steady state that needs a constant h = 5.33 W/(m² K)
        # at every flux: the submissions included convection to 293 K.
        push!(v, "h5.33" => merge(base, (h = 5.33,)))
    end
    if id in ("DBI_1", "DBI_calibrated")
        # Their FDS RAMPs evidently read the K-valued cp and k tables as °C.
        d = deepcopy(deck)
        for c in d["components"], key in ("k", "cp")
            c[key]["form"] == "table" && (c[key]["T"] = c[key]["T"] .+ 273.15)
        end
        push!(v, "tablesC" => merge(base, (deck = d,)))
    end
    if id == "Aalto_I"
        # Their submissions start (and re-radiate) at 303.15 K, not 293 K.
        push!(v, "T303" => merge(base, (T_amb = 303.15,)))
    end
    return v
end

function stats_row(sol)
    s = sol.statistics
    return join((string(getfield(s, f)) for f in fieldnames(typeof(s))), ",")
end

# Deactivation events: cells retired, the first and last event, and the largest
# number of cells retired at one time with that time.
function events_row(sol)
    ev = sol.extras.deactivation_times
    isempty(ev) && return "0,NaN,NaN,0,NaN"
    u = unique(ev)
    counts = [count(==(t), ev) for t in u]
    k = argmax(counts)
    return @sprintf(
        "%d,%.3f,%.3f,%d,%.3f",
        length(ev),
        first(ev),
        last(ev),
        counts[k],
        u[k]
    )
end

# The radiation budget of a radiative-transfer run [W/m²].
function extras_series(sol)
    x = sol.extras
    return (
        time_s = sol.t,
        radiation_absorbed = x.radiation_absorbed,
        radiation_emitted = x.radiation_emitted,
        radiation_back = x.radiation_back,
        surface_emissivity = x.surface_emissivity,
        surface_heat_flux_net = x.surface_heat_flux_net,
    )
end

function run_set(deck; tga_only = false, cells_per_mm = 20)
    id = deck["id"]
    tag = cells_per_mm == 20 ? "" : "_c$(cells_per_mm)"
    logpath = joinpath(OUT, id, "runs$tag.csv")
    mkpath(dirname(logpath))
    io = open(logpath, "w")
    fields = join(string.(fieldnames(Pyrolysis.SolverStatistics)), ",")
    println(
        io,
        "set,variant,case,cells_per_mm,retcode,outcome,t_last,solve_time,",
        "cells_retired,first_event,last_event,largest_batch,largest_batch_time,$fields",
    )
    function record(variant, case, sol)
        println(
            io,
            @sprintf(
                "%s,%s,%s,%d,%s,%s,%.3f,%.3f,%s,%s",
                id,
                variant,
                case,
                cells_per_mm,
                sol.retcode,
                sol.outcome,
                sol.t[end],
                sol.solve_time,
                events_row(sol),
                stats_row(sol)
            )
        )
        flush(io)
    end
    for β in (tag == "" ? TGA_RATES : ())
        sol = solve(tga_problem(deck, β); deactivation = true)
        write_series(joinpath(OUT, id, "base", "tga_$(β)K.csv"), tga_series(sol, β))
        record("base", "tga_$(β)K", sol)
        @printf("%-10s tga %3d K/min  %-10s %6.2f s\n", id, β, sol.retcode, sol.solve_time)
        flush(stdout)
    end
    if !(tga_only || get(deck, "kinetics_only", false))
        for (vname, v) in variants(deck), (q, L, t_end) in GASIFICATION_CASES
            case = "gasification_$(q)kW_$(L)mm"
            try
                p = gasification_problem(
                    v.deck,
                    q,
                    L;
                    t_end,
                    cells_per_mm,
                    h = v.h,
                    T_amb = v.T_amb,
                    distribution = v.distribution,
                )
                sol = solve(
                    p;
                    SOLVE_KWARGS...,
                    radiation_model = MacfpPMMA.radiation_model(v.deck),
                )
                write_series(
                    joinpath(OUT, id, vname * tag, "$case.csv"),
                    gasification_series(sol),
                )
                MacfpPMMA.radiation_model(v.deck) == RADIATIVE_TRANSFER && write_series(
                    joinpath(OUT, id, vname * tag, "$(case)_extras.csv"),
                    extras_series(sol),
                )
                record(vname, case, sol)
                @printf(
                    "%-10s %-11s %-22s %-10s %-13s t=%7.1f %7.2f s\n",
                    id,
                    vname,
                    case,
                    sol.retcode,
                    sol.outcome,
                    sol.t[end],
                    sol.solve_time
                )
                flush(stdout)
            catch err
                println(
                    io,
                    "$id,$vname,$case,$cells_per_mm,ERROR,$(replace(sprint(showerror, err), r"[\n,]" => " ")[1:min(end, 300)])",
                )
                flush(io)
                @warn "failed" id vname case exception = (err, catch_backtrace())
            end
        end
    end
    close(io)
end

function main(args)
    check_matl_db()
    tga_only = "--tga-only" in args
    ca = findfirst(a -> startswith(a, "--cells="), args)
    cells = ca === nothing ? 20 : parse(Int, split(args[ca], "=")[2])
    wanted = filter(a -> !startswith(a, "--"), args)
    decks = load_decks()
    isempty(wanted) || filter!(d -> d["id"] in wanted, decks)
    for d in decks
        run_set(d; tga_only, cells_per_mm = cells)
    end
end

main(ARGS)
