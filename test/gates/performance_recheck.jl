# Warm solve times and allocations for three representative configurations.
# Run all configurations in one process, one warm-up and five measured solves
# each, with a full collection before each measured solve. Run on a quiet
# machine with one Julia thread:
#   julia --project=test --threads=1 test/gates/performance_recheck.jl

using Pyrolysis
using Statistics: median
using Printf

module IdentityDecks
include("burnout_identity.jl")
end

module AuditDecks
include("audit_decks.jl")
end

const OUT = joinpath(@__DIR__, "measurements")
include("stamp.jl")

function measured(f)
    f()
    times = Float64[]
    bytes = Int[]
    steps = Int[]
    for _ = 1:5
        GC.gc()
        result = @timed f()
        succeeded(result.value) || error("solve ended with $(result.value.retcode)")
        push!(times, result.time)
        push!(bytes, result.bytes)
        push!(steps, result.value.statistics.accepted)
    end
    return (; times, bytes, steps)
end

function main()
    Threads.nthreads() == 1 || error("use one Julia thread")
    occursin("uncommitted", COMMIT_STAMP) && error("measure a clean commit")
    load_before = Sys.loadavg()
    load_before[1] < 2 || error("one-minute load must be below 2")
    _, cases = IdentityDecks.load_cases()
    case = only(c for c in cases if c.id == "V-BO-01")
    burnout = Base.invokelatest(case.build, 60)
    audit, audit_kw = AuditDecks.problem_for(:ale, 100; solid_λ = 1e-6, outputs = 81)
    fir = IdentityDecks.douglas_fir_problem()
    runs = (
        "V-BO-01, 60 PMMA cells" => (
            () -> Base.invokelatest(
                solve,
                burnout;
                case.solve_kwargs...,
                energy_state = :temperature,
            )
        ),
        "ale 100, solid gas-transfer coefficient 1e-6" =>
            (() -> solve(audit; audit_kw..., reltol = 1e-6, abstol = 1e-8)),
        "Douglas fir, 60 s" =>
            (() -> IdentityDecks.douglas_fir_solve(fir, :min_thickness => 0.5e-3)),
    )
    rows = [(; name, measured(f)...) for (name, f) in runs]
    load_after = Sys.loadavg()
    mkpath(OUT)
    open(joinpath(OUT, "performance_recheck.md"), "w") do io
        println(io, "# Warm solve performance\n")
        println(
            io,
            "Measured at $(COMMIT_STAMP), Julia $(VERSION), $(Sys.CPU_NAME), one thread.",
        )
        println(
            io,
            "Load before: $(load_before); after: $(load_after). Five solves after one warm-up; GC before each.\n",
        )
        println(
            io,
            "| configuration | median [s] | range [s] | allocated [MiB], range | accepted steps |",
        )
        println(io, "|---|---|---|---|---|")
        for r in rows
            @printf(
                io,
                "| %s | %.6g | %.6g–%.6g | %.6g–%.6g | %s |\n",
                r.name,
                median(r.times),
                minimum(r.times),
                maximum(r.times),
                minimum(r.bytes) / 2^20,
                maximum(r.bytes) / 2^20,
                join(unique(r.steps), ", ")
            )
        end
    end
    open(joinpath(OUT, "performance_recheck.csv"), "w") do io
        println(io, "configuration,repeat,seconds,bytes,accepted_steps")
        for (j, r) in enumerate(rows), i in eachindex(r.times)
            println(io, join((j, i, r.times[i], r.bytes[i], r.steps[i]), ','))
        end
    end
    foreach(println, rows)
    load_after[1] < 2 || error("load rose above 2; repeat on a quiet machine")
end

abspath(PROGRAM_FILE) == @__FILE__() && main()
