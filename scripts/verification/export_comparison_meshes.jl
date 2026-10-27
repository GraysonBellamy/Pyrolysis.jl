# Export actual report meshes without solving or loading the plotting stack.
# Used by audit_comparison_basis.py; n is a cell count, never a node count.
using TOML
using SHA
include(joinpath(@__DIR__, "..", "..", "test", "verification", "casedefs.jl"))

function export_comparison_meshes(path)
    rows = Dict{String,Any}[]
    reference = joinpath(@__DIR__, "..", "..", "test", "verification", "reference")
    active = Set([
        "V-HC-01",
        "V-HC-03",
        "V-RAD-01",
        "V-RX-01",
        "V-CPL-01",
        "V-GT-07a",
        "V-GT-07b",
        "V-GT-08a",
        "V-GT-08b",
        "V-RX-05c",
    ])
    for solver in ("fds", "gpyro", "thermakin")
        for entry in
            TOML.parsefile(joinpath(reference, "matched", "$(solver).toml"))["series"]
            push!(active, entry["case"])
        end
    end
    for case in load_verification_cases()
        case.id in active || continue
        isempty(pairs(case.report)) && continue
        ns = collect(Int, get(case.report, :n_cells_ladder, [case.ci.n_cells]))
        # Include all ladder rungs: a matching rung may already exist even
        # though the plotted curve is taken from the last rung.
        for n in ns
            problem = Base.invokelatest(case.build, n)
            faces = first.(problem.mesh.positions)
            widths = diff(faces)
            push!(
                rows,
                Dict(
                    "case" => case.id,
                    "requested_cells" => n,
                    "cells" => length(widths),
                    "plotted" => n == last(ns),
                    "faces_m" => faces,
                    "spacing_min_m" => minimum(widths),
                    "spacing_max_m" => maximum(widths),
                    "t_start_s" => first(problem.tspan),
                    "t_end_s" => last(problem.tspan),
                    "saved_times_s" => problem.output_times,
                ),
            )
        end
    end
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        root = normpath(joinpath(@__DIR__, "..", ".."))
        sources = vcat(
            [joinpath(root, "test", "verification", "casedefs.jl")],
            readdir(joinpath(root, "test", "verification", "cases"); join = true),
        )
        hashes = Dict(
            relpath(f, root) => bytes2hex(sha256(read(f))) for
            f in sources if endswith(f, ".jl")
        )
        TOML.print(io, Dict("meshes" => rows, "source_hashes" => hashes); sorted = true)
    end
end

length(ARGS) == 1 || error("usage: export_comparison_meshes.jl OUTPUT.toml")
export_comparison_meshes(only(ARGS))
