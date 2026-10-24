# The record a deletion commit is compared against: the solutions of the
# identity decks and the containment suite's Jacobians, as digests, so that
# a commit that deletes and renames only can show it changes no number.
#
# The decks are `burnout_identity.jl`'s (eight cases at their CI resolutions
# with their solve keywords, and the Douglas fir deck over its first 60 s),
# each in both energy forms. Per solve: the digest of the saved temperatures,
# masses, positions and times; the digest of the event records (every field
# of every event, vectors by value); the outcome and the return code; the
# accepted and rejected steps. The containment suite
# (`jacobian/test_containment.jl`) runs in a module whose `dense_jacobian`
# and `jacobian!` record the digest of every Jacobian they return, in call
# order, while its tests run as usual; its test counts are recorded beside.
#
# The digests are `hash` values, so the comparison holds in one Julia binary
# and one resolved manifest.
#
# With `--record` the record is written to
# measurements/deletion_identity_reference.csv; without it every digest is
# compared with the record and measurements/deletion_identity.md says which
# are equal.
#
# Run from the repository root:
#   julia --project=ENV --startup-file=no --threads=1 deletion_identity.jl [--record]

include(joinpath(@__DIR__, "burnout_identity.jl"))

using Test

const IDENTITY_REFERENCE = joinpath(OUT, "deletion_identity_reference.csv")
const CONTAINMENT = joinpath(@__DIR__, "..", "jacobian", "test_containment.jl")

events_digest(sol) = string(
    hash([
        Tuple(getfield(e, f) for f in fieldnames(Pyrolysis.Deactivation.DeactivationEvent)) for e in sol.events
    ]);
    base = 16,
)

# The cases loaded in `form`: the suite routes them to the finite-volume core
# in the energy state it reads when `casedefs.jl` loads.
function cases_in(form)
    ENV["PYROLYSIS_ENERGY_STATE"] = string(form)
    return load_cases()
end

function deck_records(form)
    m, cases = cases_in(form)
    rows = Any[]
    record(deck, sol) = begin
        r = (;
            deck,
            form,
            fields = digest(sol),
            events = events_digest(sol),
            outcome = string(sol.outcome),
            retcode = string(sol.retcode),
            accepted = sol.statistics.accepted,
            rejected = sol.statistics.rejected,
        )
        println(r)
        flush(stdout)
        push!(rows, r)
    end
    for id in DECKS
        case = only(filter(c -> c.id == id, cases))
        problem = Base.invokelatest(case.build, case.ci.n_cells)
        sol =
            Base.invokelatest(m.runs_on_core, case) ?
            Base.invokelatest(m._case_solve, case, problem; case.solve_kwargs...) :
            Base.invokelatest(
                Pyrolysis.Integration.solve,
                problem;
                energy_state = form,
                _core_keywords(case.solve_kwargs)...,
            )
        record(id, sol)
    end
    df = douglas_fir_problem()
    sol = Base.invokelatest(
        Pyrolysis.Integration.solve,
        df;
        energy_state = form,
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
        abstol = 1e-9,
        reltol = 1e-7,
        dtmin = 1e-14,
    )
    record("Douglas fir, 60 s", sol)
    return rows
end

# The containment suite with its Jacobians recorded: `(kind, digest)` in
# call order, and the suite's pass and fail counts.
function containment_record()
    m = Module(:containment_record)
    Core.eval(
        m,
        quote
            using Pyrolysis
            const JACOBIANS = Tuple{Symbol,String}[]
            function dense_jacobian(args...; kwargs...)
                J = Pyrolysis.ColoredJacobian.dense_jacobian(args...; kwargs...)
                push!(JACOBIANS, (:dense, string(hash(J); base = 16)))
                return J
            end
            function jacobian!(args...; kwargs...)
                J = Pyrolysis.ColoredJacobian.jacobian!(args...; kwargs...)
                push!(
                    JACOBIANS,
                    (:colored, string(hash((J.colptr, J.rowval, J.nzval)); base = 16)),
                )
                return J
            end
            include(p) = Base.include($m, p)
        end,
    )
    ts = @testset "containment" begin
        Base.include(m, CONTAINMENT)
    end
    counts = Test.get_test_counts(ts)
    jacobians = copy(Base.invokelatest(getglobal, m, :JACOBIANS))
    return (;
        jacobians,
        passes = counts.passes + counts.cumulative_passes,
        fails = counts.fails + counts.cumulative_fails,
        errors = counts.errors + counts.cumulative_errors,
    )
end

const DECK_FIELDS = (:fields, :events, :outcome, :retcode, :accepted, :rejected)

function read_identity(path)
    ref = Dict{String,String}()
    isfile(path) || return ref
    for line in readlines(path)[2:end]
        key, value = rsplit(line, ","; limit = 2)
        ref[strip(key, '"')] = value
    end
    return ref
end

function identity_main(args)
    record_mode = "--record" in args
    entries = Pair{String,String}[]
    decks = vcat(deck_records(:temperature), deck_records(:enthalpy))
    for r in decks, f in DECK_FIELDS
        push!(entries, "$(r.deck) | $(r.form) | $(f)" => string(getfield(r, f)))
    end
    c = containment_record()
    push!(entries, "containment | count" => string(length(c.jacobians)))
    for (k, (kind, d)) in enumerate(c.jacobians)
        push!(entries, "containment | $(k) | $(kind)" => d)
    end
    reference = read_identity(IDENTITY_REFERENCE)
    if record_mode
        open(IDENTITY_REFERENCE, "w") do io
            println(io, "entry,value")
            for (k, v) in entries
                println(io, "\"", k, "\",", v)
            end
        end
    end
    equal = [get(reference, k, nothing) == v for (k, v) in entries]
    open(joinpath(OUT, "deletion_identity.md"), "w") do io
        println(io, "# The deletion identity record\n")
        println(
            io,
            "Produced by `deletion_identity.jl` at $(COMMIT_STAMP), Julia $(VERSION). ",
            record_mode ?
            "This run recorded the reference `deletion_identity_reference.csv`. " :
            "Each entry against the reference `deletion_identity_reference.csv`. ",
            "The containment suite: $(c.passes) passed, $(c.fails) failed, $(c.errors) errored; ",
            "$(length(c.jacobians)) Jacobians recorded",
            record_mode ? ".\n" :
            ", $(count(equal[(end - length(c.jacobians)):end])) of $(length(c.jacobians) + 1) containment entries equal.\n",
        )
        println(
            io,
            "| deck | form | fields | events | outcome | retcode | accepted (rejected) | equal |",
        )
        println(io, "|---|---|---|---|---|---|---|---|")
        for r in decks
            same =
                record_mode ? "recorded" :
                all(
                    get(reference, "$(r.deck) | $(r.form) | $(f)", nothing) ==
                    string(getfield(r, f)) for f in DECK_FIELDS
                ) ? "yes" : "**no**"
            println(
                io,
                "| $(r.deck) | $(r.form) | `$(r.fields)` | `$(r.events)` | $(r.outcome) | $(r.retcode) | ",
                "$(r.accepted) ($(r.rejected)) | $(same) |",
            )
        end
        if !record_mode
            differing = [k for ((k, _), e) in zip(entries, equal) if !e]
            println(
                io,
                "\nEntries that differ from the reference: ",
                isempty(differing) ? "none." : join(differing, "; "),
            )
        end
    end
    println("wrote $(OUT)/deletion_identity.md", record_mode ? " and the reference" : "")
    return nothing
end

abspath(PROGRAM_FILE) == @__FILE__() && identity_main(ARGS)
