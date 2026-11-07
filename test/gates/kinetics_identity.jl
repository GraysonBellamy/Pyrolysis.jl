# The record that proves a change to the reaction kinetics leaves first-order
# material bit-identical: one energy form per process, each deck solved once,
# and the process fails on a digest that differs from its reference.
#
# The decks: a table of 3072 rates and heats through `evaluate_reactions!` at
# order one; the verification decks below at their CI resolutions with their
# solve keywords, built and selected as `burnout_identity.jl` builds them; and
# the three examples at their native resolutions over truncated spans, with
# each example's complete solve keywords. Per solve the digest of the saved
# temperatures, masses, positions and times, the return code and the accepted
# and rejected steps.
#
# The reference, measurements/kinetics_identity_reference.csv, holds one row
# per energy form and deck. Without `--record` a digest that differs from its
# row, a missing reference or row and an unexpected return code each exit
# nonzero; with `--record` the rows of the decks solved are replaced and an
# unexpected return code still exits nonzero. `--only=ID,ID` solves the named
# decks alone (`--list` prints the names). Including this file defines the
# decks and solves nothing.
#
# Run from the repository root, once per form:
#   PYROLYSIS_ENERGY_STATE=temperature julia --project=test --startup-file=no --threads=1 \
#       test/gates/kinetics_identity.jl [--record] [--only=ID,ID]

const FORM = let form = get(ENV, "PYROLYSIS_ENERGY_STATE", "temperature")
    form in ("temperature", "enthalpy") ||
        error("PYROLYSIS_ENERGY_STATE must be temperature or enthalpy, got $form")
    form
end
ENV["PYROLYSIS_SMOKE"] = "0"

# Sets `PYROLYSIS_ENERGY_STATE` as it loads and solves nothing when included.
module IdentityDecks
include(joinpath(@__DIR__, "burnout_identity.jl"))
end

# `casedefs.jl` captures the form when it loads, in `load_cases()`.
ENV["PYROLYSIS_ENERGY_STATE"] = FORM

using Pyrolysis
using Pyrolysis.Physics: evaluate_reactions!
using Printf

const REFERENCE = joinpath(@__DIR__, "measurements", "kinetics_identity_reference.csv")
const HEADER = "form,deck,digest,retcode,accepted,rejected"

const RATES = "rates-order-one"

# The verification decks. V-RX-02d holds an order other than one; V-RAD-06
# and V-RAD-07 solve under radiative transfer.
const CASES = (
    "V-RX-01",
    "V-RX-02d",
    "V-RX-03",
    "V-RX-04a",
    "V-RX-04b",
    "V-RX-05a",
    "V-RX-05b",
    "V-RX-05c",
    "V-CPL-01",
    "V-CPL-01b",
    "V-CPL-02",
    "V-CPL-02b",
    "V-CPL-03",
    "V-CPL-05",
    "V-CPL-06a",
    "V-HC-05",
    "V-BO-04",
    "V-RAD-06",
    "V-RAD-07",
)

# The examples: the deck between `marker` and the solve, without its printing,
# and the keywords of the example's solve.
const CONE_KEYWORDS = (;
    abstol = 1e-9,
    reltol = 1e-7,
    dtmin = 1e-14,
    use_ale = true,
    radiation_model = SURFACE_ABSORPTION,
    min_thickness = 0.5e-3,
)
const EXAMPLES = (
    (;
        id = "examples/06-300s",
        file = "06_cone_calorimeter.jl",
        marker = "pmma_cp = CallableProperty(",
        seconds = 300.0,
        keywords = (;
            abstol = 1e-11,
            reltol = 1e-9,
            pressure_tolerance = 1e-4,
            dtmin = 1e-14,
            use_ale = true,
            min_thickness = 0.01e-3,
            deactivation = true,
        ),
    ),
    (;
        id = "examples/07-60s",
        file = "07_douglas_fir_cone.jl",
        marker = "const T_HFM_MAX",
        seconds = 60.0,
        keywords = CONE_KEYWORDS,
    ),
    (;
        id = "examples/08-60s",
        file = "08_pressure_treated_wood_cone.jl",
        marker = "ptw_cp = CallableProperty(",
        seconds = 60.0,
        keywords = CONE_KEYWORDS,
    ),
)

# The digest of 384 states of a four-reaction material at order one, a rate and
# a heat per reaction: limiter thresholds 1e-3, 0.1 and 10 kg/m³, on and off, 16
# concentrations from −1 to 10³ kg/m³ at four temperatures, ungated, gated from
# below, from above and from both sides.
function rate_digest()
    solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    reactions = (
        Reaction(:open, 1 => 2; A = 1e10, E = 1.5e5, h = 5e5),
        Reaction(:low, 1 => 2; A = 1e8, E = 1.2e5, h = -2e5, T_min = 400.0),
        Reaction(:high, 1 => 2; A = 1e6, E = 9e4, h = 1e5, T_max = 800.0),
        Reaction(:both, 1 => 2; A = 1e4, E = 6e4, h = 3e5, T_min = 400.0, T_max = 800.0),
    )
    material = Material(name = :rates, components = (solid, gas), reactions = reactions)
    ξs = (-1.0, -1e-8, -0.0, 0.0, 1e-20, 1e-12, 1e-8, 1e-6)
    ξs = (ξs..., 1e-4, 1e-3, 1e-2, 0.1, 1.0, 10.0, 100.0, 1000.0)
    rates, heats = zeros(4), zeros(4)
    table = Float64[]
    for threshold in (1e-3, 0.1, 10.0), enabled in (true, false)
        limiter = DepletionLimiter(; threshold, enabled)
        for T in (300.0, 450.0, 700.0, 1200.0), ξ in ξs
            evaluate_reactions!(rates, heats, (ξ, 0.0), T, material, limiter)
            append!(table, rates)
            append!(table, heats)
        end
    end
    length(table) == 3072 || error("the rate table has $(length(table)) entries")
    return string(hash(table); base = 16)
end

# The example `e` at its native resolution over its truncated span.
function example_problem(e)
    lines = readlines(joinpath(@__DIR__, "..", "..", "examples", e.file))
    a = findfirst(startswith(e.marker), lines)
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
    m = Module(Symbol(replace(e.file, ".jl" => "")))
    Core.eval(m, :(using Pyrolysis))
    Core.eval(m, :(const SMOKE = false))
    Core.eval(m, Meta.parseall(join(out, "\n")))
    # A second `eval` runs in the latest world, where the deck's bindings exist.
    return Core.eval(
        m,
        :(PyrolysisProblem(
            mesh = problem.mesh,
            material = problem.material,
            bc_set = problem.bc_set,
            T_initial = problem.T_initial,
            ξ_initial = problem.ξ_initial,
            tspan = (0.0, $(e.seconds)),
            dt_initial = problem.dt_initial,
            output_times = collect(0.0:1.0:($(e.seconds))),
        )),
    )
end

# The solve of the example `e`: returns the solution.
function solve_example(e)
    problem = example_problem(e)
    return Base.invokelatest(
        Pyrolysis.solve,
        problem;
        e.keywords...,
        energy_state = Symbol(FORM),
    )
end

function read_reference()
    ref = Dict{Tuple{String,String},String}()
    isfile(REFERENCE) || return ref
    for line in readlines(REFERENCE)[2:end]
        form, deck, rest = split(line, ","; limit = 3)
        ref[(String(form), String(strip(deck, '"')))] = rest
    end
    return ref
end

function write_reference(ref)
    open(REFERENCE, "w") do io
        println(io, HEADER)
        for (form, deck) in sort!(collect(keys(ref)))
            println(io, form, ",\"", deck, "\",", ref[(form, deck)])
        end
    end
end

row(sol) = join(
    (
        IdentityDecks.digest(sol),
        sol.retcode,
        sol.statistics.accepted,
        sol.statistics.rejected,
    ),
    ",",
)

function main(args)
    record = "--record" in args
    only_arg = findfirst(startswith("--only="), args)
    names = (RATES, CASES..., (e.id for e in EXAMPLES)...)
    if "--list" in args
        foreach(println, names)
        return 0
    end
    selected =
        only_arg === nothing ? collect(names) :
        String.(split(args[only_arg][(length("--only=")+1):end], ","))
    unknown = setdiff(selected, names)
    isempty(unknown) || error("unknown decks: $(join(unknown, ", "))")

    reference = read_reference()
    if !record && !isfile(REFERENCE)
        println("FAIL: no reference at $REFERENCE")
        return 1
    end

    m, cases = IdentityDecks.load_cases()
    known_cases = (IdentityDecks.CHAR_LAYER, (c.id for c in cases)...)
    missing = setdiff(intersect(selected, CASES), known_cases)
    isempty(missing) || error("no case builds $(join(missing, ", "))")
    failures = String[]
    for id in selected
        expected, measured, seconds = :Success, "", 0.0
        wrong = false
        if id == RATES
            seconds = @elapsed measured = "$(rate_digest()),Success,0,0"
        else
            example = findfirst(e -> e.id == id, EXAMPLES)
            local sol
            if example !== nothing
                seconds = @elapsed sol = solve_example(EXAMPLES[example])
            else
                id == IdentityDecks.CHAR_LAYER ||
                    (expected = only(filter(c -> c.id == id, cases)).expected_retcode)
                _, solve, _ = IdentityDecks.deck(m, cases, id)
                seconds = @elapsed sol = solve()
            end
            measured = row(sol)
            if sol.retcode != expected
                push!(failures, "$id: return code $(sol.retcode), expected $expected")
                wrong = true
            end
        end
        known = get(reference, (FORM, id), nothing)
        status = if record
            "recorded"
        elseif known === nothing
            push!(failures, "$id: no reference row for the $FORM form")
            "NO REFERENCE"
        elseif first(split(known, ",")) != first(split(measured, ","))
            push!(failures, "$id: digest $measured, reference $known")
            "DIFFERS"
        else
            "equal"
        end
        # Re-read: the other form's process records into the same file.
        if record && !wrong
            reference = read_reference()
            reference[(FORM, id)] = measured
            write_reference(reference)
        end
        @printf("%-22s %-9s %8.1f s  %-12s %s\n", id, FORM, seconds, status, measured)
        flush(stdout)
    end
    for f in failures
        println("FAIL: ", f)
    end
    println(isempty(failures) ? "PASS" : "FAILED ($(length(failures)))")
    return isempty(failures) ? 0 : 1
end

abspath(PROGRAM_FILE) == @__FILE__() && exit(main(ARGS))
