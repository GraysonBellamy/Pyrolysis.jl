# Gate item 5 of the Phase 2 plan (§7.3): the decks that must not change under
# the correction, hashed on the finite-volume core. Each deck is solved as the
# suite solves it (its builder at the CI resolution, its solve keywords) and
# the digest of `(T, m, z, t)` of its solution recorded; the Douglas fir deck
# of examples/07_douglas_fir_cone.jl over its first 60 s joins them.
#
# With `--record` the digests become the reference
# (measurements/burnout_identity_reference.csv, taken before any Phase 2 code
# landed); without it every digest is compared with the reference and the
# table says which decks are bit-identical.
#
# Run from the repository root:
#   julia --project=test --startup-file=no --threads=1 test/newcore/gates/burnout_identity.jl [--record]
# Writes test/newcore/gates/measurements/burnout_identity.md.

ENV["PYROLYSIS_CORE"] = "new"
ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
using Pyrolysis
using Printf
include(joinpath(@__DIR__, "stamp.jl"))

const OUT = joinpath(@__DIR__, "measurements")
mkpath(OUT)
const REFERENCE = joinpath(OUT, "burnout_identity_reference.csv")
const RECORD = "--record" in ARGS
const DECKS =
    ("V-CPL-02", "V-CPL-02b", "V-CPL-03", "V-CPL-06a", "V-HC-05", "V-RAD-01", "V-RX-02d")
const DF_SECONDS = 60.0

function load_cases()
    m = Module(:cases_identity)
    Core.eval(m, :(include(p) = Base.include($m, p)))
    Base.include(m, joinpath(@__DIR__, "..", "..", "verification", "casedefs.jl"))
    return m, Base.invokelatest(m.load_verification_cases)
end

# The old core's depletion keywords have no counterpart on the finite-volume
# core; `core_solve` strips them for the routed cases, and this does the same
# for a case the suite does not route.
function _core_keywords(kw)
    drop = (:depletion_threshold, :min_cells, :jacobian, :linear_solver)
    kept = [k => v for (k, v) in pairs(kw) if !(k in drop)]
    nt = (; kept...)
    if haskey(nt, :handle_depletion)
        d = nt.handle_depletion
        nt = (;
            (k => v for (k, v) in pairs(nt) if k != :handle_depletion)...,
            deactivation = d,
        )
    end
    return nt
end

digest(sol) = string(hash((sol.T, sol.m, sol.z, sol.t)); base = 16)

# The Douglas fir deck: the material, mesh and boundary conditions of the
# example, without its printing, plotting and solve.
function douglas_fir_problem()
    lines = readlines(
        joinpath(@__DIR__, "..", "..", "..", "examples", "07_douglas_fir_cone.jl"),
    )
    a = findfirst(startswith("const T_HFM_MAX"), lines)
    b = findfirst(startswith("solution = solve("), lines)
    out = String[]
    skip = false
    for line in lines[a:(b-1)]
        if skip
            startswith(line, ")") && (skip = false)
            continue
        end
        if startswith(line, "println(")
            skip = !endswith(rstrip(line), ")")
            continue
        end
        startswith(line, "flush(") && continue
        push!(out, line)
    end
    m = Module(:douglas_fir_deck)
    Core.eval(m, :(using Pyrolysis))
    Core.eval(m, Meta.parseall(join(out, "\n")))
    # Built by a second `eval`, which runs in the latest world, so the deck's
    # fresh bindings are visible.
    return Core.eval(
        m,
        :(PyrolysisProblem(
            mesh = mesh,
            material = douglas_fir,
            bc_set = bc_set,
            T_initial = T_initial_func,
            ξ_initial = ξ_initial_funcs,
            tspan = (0.0, $DF_SECONDS),
            dt_initial = dt_initial,
            output_times = collect(0.0:1.0:($DF_SECONDS)),
        )),
    )
end

function main()
    m, cases = load_cases()
    rows = []
    for id in DECKS
        case = only(filter(c -> c.id == id, cases))
        problem = Base.invokelatest(case.build, case.ci.n_cells)
        sol = if Base.invokelatest(m.runs_on_core, case)
            Base.invokelatest(m._case_solve, case, problem; case.solve_kwargs...)
        else
            Base.invokelatest(
                Pyrolysis.Integration.solve,
                problem;
                _core_keywords(case.solve_kwargs)...,
            )
        end
        push!(
            rows,
            (;
                deck = id,
                cells = length(problem.mesh.cells),
                digest = digest(sol),
                steps = sol.statistics.accepted,
                retcode = sol.retcode,
            ),
        )
        println(rows[end])
        flush(stdout)
    end
    df = douglas_fir_problem()
    sol = Base.invokelatest(
        Pyrolysis.Integration.solve,
        df;
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        abstol = 1e-9,
        reltol = 1e-7,
        dtmin = 1e-14,
    )
    push!(
        rows,
        (;
            deck = "Douglas fir, 60 s",
            cells = length(df.mesh.cells),
            digest = digest(sol),
            steps = sol.statistics.accepted,
            retcode = sol.retcode,
        ),
    )
    println(rows[end])

    reference = Dict{String,String}()
    if isfile(REFERENCE)
        for line in readlines(REFERENCE)[2:end]
            deck, dig = rsplit(line, ","; limit = 2)
            reference[strip(deck, '"')] = dig
        end
    end
    if RECORD
        open(REFERENCE, "w") do io
            println(io, "deck,digest")
            for r in rows
                println(io, "\"", r.deck, "\",", r.digest)
            end
        end
    end
    open(joinpath(OUT, "burnout_identity.md"), "w") do io
        println(io, "# Burnout identity: the decks the correction must not change\n")
        println(
            io,
            "Generated by `test/newcore/gates/burnout_identity.jl` at $(COMMIT_STAMP), finite-volume core, temperature form, each deck at its CI resolution with its solve keywords; the digest is `hash((T, m, z, t))` of the solution. The reference digests are `burnout_identity_reference.csv`$(RECORD ? ", recorded by this run" : "").\n",
        )
        println(io, "| deck | cells | steps | retcode | digest | equal to the reference |")
        println(io, "|---|---|---|---|---|---|")
        for r in rows
            ref = get(reference, r.deck, "")
            eq =
                RECORD ? "recorded" :
                (isempty(ref) ? "no reference" : (ref == r.digest ? "yes" : "**no**"))
            println(
                io,
                "| ",
                r.deck,
                " | ",
                r.cells,
                " | ",
                r.steps,
                " | ",
                r.retcode,
                " | `",
                r.digest,
                "` | ",
                eq,
                " |",
            )
        end
    end
    return rows
end

main()
