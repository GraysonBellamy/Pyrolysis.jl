# Burnout through the event path against the exact solution: a column at a rate
# that reads no temperature, reactant A at order n on `basis = 1` with its
# initial conversion graded from 0.99 at the back to 0.05 at the exposed face, a
# filler B stripped at first order, no correction, deactivation at ε. Every
# material point has a closed form, so the column's condensed mass and its gas
# generation are known. The ladder runs n ∈ {0.25, 0.5, 0.75, 1.2} on 20, 40 and
# 80 cells at ε ∈ {1e-3, 1e-4}, with the n = 0 rows beside it, whose staircase is
# the discretisation of a zero-order front and is recorded only.
#
# Per run: the largest condensed-mass error over the initial mass, the largest
# gas-generation error after the filler is gone over the exact peak of that
# window, the integrated |Δrate| over the integrated exact rate, the events and
# the outcome. A run passes when it ends with `min_thickness` or `burned_out` at
# its full coverage of the burn and, where the table `EXPECTED` holds a row, within
# three times that row. The ladder stops at the first failing row's end with a
# nonzero exit; `--only-order=n` runs one order.
#
# Run from the repository root, once per form:
#   PYROLYSIS_ENERGY_STATE=temperature julia --project=test --startup-file=no --threads=1 \
#       test/gates/burnout_oracle.jl
# Writes test/gates/measurements/burnout_oracle_<form>.md.

using Pyrolysis
using Pyrolysis.Closures: NoCorrection
using Printf
include(joinpath(@__DIR__, "stamp.jl"))

const FORM = Symbol(get(ENV, "PYROLYSIS_ENERGY_STATE", "temperature"))
const OUT = joinpath(@__DIR__, "measurements")
const L0, T0, ρ0, FLOW, FHIGH = 10e-3, 600.0, 500.0, 0.05, 0.99
const KB = 2.0

# (rate constant, horizon) of each order.
const DECK = Dict(
    0.0 => (0.02, 52.0),
    0.25 => (0.02, 75.0),
    0.5 => (0.02, 110.0),
    0.75 => (0.02, 220.0),
    1.2 => (0.2, 300.0),
)

# (n, N, ε) => (mass error, rate error, L1 error) measured on the event path.
const EXPECTED = Dict(
    (0.5, 40, 1e-3) => (3.0e-5, 1.4e-3, 2.2e-4),
    (0.5, 40, 1e-4) => (2.0e-5, 3.7e-4, 9.0e-5),
    (0.5, 20, 1e-3) => (7.9e-5, 2.0e-3, 5.0e-4),
    (0.5, 80, 1e-3) => (2.7e-5, 7.1e-4, 1.2e-4),
    (0.25, 40, 1e-3) => (3.1e-5, 4.0e-3, 1.3e-3),
    (0.75, 40, 1e-3) => (4.8e-5, 2.0e-4, 6.3e-5),
    (1.2, 40, 1e-3) => (7.3e-5, 1.5e-4, 7.9e-5),
)

function oracle_problem(n, N, k, tend)
    A = SolidComponent(:A, ρ = ρ0, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    B = SolidComponent(:B, ρ = ρ0, φ = 0.3, c = 1000.0, k = 0.2, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    material = Material(
        name = :oracle,
        components = (A, B, gas),
        reactions = (
            Reaction(:burn, 1 => 3; A = k, E = 0.0, h = 0.0, n = n, basis = 1.0),
            Reaction(:strip, 2 => 3; A = KB, E = 0.0, h = 0.0),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1e-3, conversion_threshold = 1e-6),
    )
    f(z) = FHIGH - (FHIGH - FLOW) * z / L0
    return PyrolysisProblem(
        mesh = uniform_mesh(L0, N),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = T0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> ρ0 * f(z), z -> ρ0 * (1 - f(z)), z -> 0.0],
        tspan = (0.0, tend),
        output_times = collect(0.0:(tend/1000):tend),
    )
end

# The exact condensed mass and bare rate per unit initial condensed mass, the
# midpoint rule over the column's initial composition.
function exact(n, k, t; points = 40_000)
    M, r = 0.0, 0.0
    for q = 1:points
        f = FLOW + (FHIGH - FLOW) * (q - 0.5) / points
        p = f^(1 - n) - (1 - n) * k * t
        x = p > 0 ? p^(1 / (1 - n)) : 0.0
        b = (1 - f) * exp(-KB * t)
        M += x + b
        r += (x > 0 ? k * x^n : 0.0) + KB * b
    end
    return M / points, r / points
end

function run_row(n, N, ε)
    k, tend = DECK[n]
    p = oracle_problem(n, N, k, tend)
    wall = time()
    sol = solve(
        p;
        use_ale = true,
        deactivation = true,
        threshold = ε,
        reltol = 1e-6,
        abstol = 1e-8,
        correction = NoCorrection(),
        energy_state = FORM,
    )
    wall = time() - wall
    condensed = sol.m[1, :, :] .+ sol.m[2, :, :]
    mass = vec(sum(condensed, dims = 1)) ./ sum(condensed[:, 1])
    rate = sol.extras.gas_generation_rate ./ (ρ0 * L0)
    pairs = [exact(n, k, t) for t in sol.t]
    M, r = first.(pairs), last.(pairs)
    after = sol.t .>= 6.0
    peak = maximum(r[after])
    dt = diff(sol.t)
    finite = isfinite.(rate[1:(end-1)])
    return (;
        n,
        N,
        ε,
        retcode = sol.retcode,
        outcome = sol.outcome,
        t_end = last(sol.t),
        accepted = sol.statistics.accepted,
        rejected = sol.statistics.rejected,
        events = length(sol.events),
        mass = maximum(abs.(mass .- M)),
        rate = maximum(abs.(rate[after] .- r[after])) / peak,
        l1 = sum(abs.(rate[1:(end-1)][finite] .- r[1:(end-1)][finite]) .* dt[finite]) /
             sum(r[1:(end-1)] .* dt),
        wall,
    )
end

function verdict(row)
    row.outcome in (:min_thickness, :burned_out) ||
        return "failed: outcome $(row.outcome), $(row.retcode) at $(row.t_end) s"
    expected = get(EXPECTED, (row.n, row.N, row.ε), nothing)
    expected === nothing && return "ended"
    row.mass <= 3expected[1] || return "failed: mass error above 3 × $(expected[1])"
    row.rate <= 3expected[2] || return "failed: rate error above 3 × $(expected[2])"
    row.l1 <= 3expected[3] || return "failed: L1 error above 3 × $(expected[3])"
    return "within"
end

function main()
    orders = [0.25, 0.5, 0.75, 1.2]
    for a in ARGS
        startswith(a, "--only-order=") && (orders = [parse(Float64, a[14:end])])
    end
    lines = String[]
    failed = false
    for n in orders, N in (20, 40, 80), ε in (1e-3, 1e-4)
        row = run_row(n, N, ε)
        v = verdict(row)
        failed |= startswith(v, "failed")
        line = @sprintf(
            "| %.2f | %d | %.0e | %s | %.2f | %d/%d | %d | %.2e | %.2e | %.2e | %.0f | %s |",
            n,
            N,
            ε,
            row.outcome,
            row.t_end,
            row.accepted,
            row.rejected,
            row.events,
            row.mass,
            row.rate,
            row.l1,
            row.wall,
            v
        )
        println(line)
        flush(stdout)
        push!(lines, line)
    end
    for N in (20, 40), ε in (1e-3, 1e-4)
        row = run_row(0.0, N, ε)
        line = @sprintf(
            "| 0 | %d | %.0e | %s | %.2f | %d/%d | %d | %.2e | %.2e | %.2e | %.0f | recorded |",
            N,
            ε,
            row.outcome,
            row.t_end,
            row.accepted,
            row.rejected,
            row.events,
            row.mass,
            row.rate,
            row.l1,
            row.wall
        )
        println(line)
        flush(stdout)
        push!(lines, line)
    end
    open(joinpath(OUT, "burnout_oracle_$(FORM).md"), "w") do io
        println(io, "# Burnout against the exact solution, $(FORM) form\n")
        println(io, COMMIT_STAMP)
        println(
            io,
            "\n| n | N | ε | outcome | end [s] | steps | events | mass | rate | L1 | wall [s] | verdict |",
        )
        println(io, "|---|---|---|---|---|---|---|---|---|---|---|---|")
        foreach(l -> println(io, l), lines)
    end
    failed && exit(1)
end

main()
