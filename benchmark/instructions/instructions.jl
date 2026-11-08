# Retired user-space instructions of the hot operations and of short solves, on
# a machine that may be busy: unlike time, the count does not depend on the load.
# See README.md.
#
#   taskset -c 0-11 julia --project=test --startup-file=no --threads=1 \
#       benchmark/instructions/instructions.jl [--record F] [--compare F]

# The counters' package does not resolve with the package's dependencies, so it
# keeps an environment of its own, stacked on the test environment.
push!(LOAD_PATH, @__DIR__)

using LinearAlgebra
using LinuxPerf
using Printf
using Pyrolysis
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian: jac_prototype, jacobian!

include(joinpath(@__DIR__, "..", "decks.jl"))

# Hot operations run this many times per measurement.
const REPEATS = 200

_counts(s) = [Float64(c.value) for g in s.threads[1].groups for c in g]

# Instructions, cycles and bytes of `f()` after a warm-up.
function measure(f)
    f()
    GC.gc()
    bytes = Ref(0)
    s = LinuxPerf.@pstats "(instructions:u,cpu-cycles:u)" (bytes[] = @allocated f())
    instructions, cycles = _counts(s)
    instructions > 0 || error(
        "the counter read zero: the process is on a core without one; pin it with taskset",
    )
    return (; instructions, cycles, bytes = bytes[])
end

function rows()
    out = Tuple{String,NamedTuple}[]
    burnout = (; deactivation = true, use_ale = true)
    decks = (
        (
            name = "polymer_burnout",
            problem = () -> _column(_pmma(); t_final = 80.0, L = 1e-3, n = 30),
            options = burnout,
        ),
        (
            name = "wood_moving_radiative_transfer_darcy",
            problem = () -> _column(_wood(); darcy = true, t_final = 5.0),
            options = _CONFIGURATIONS.wood_moving_radiative_transfer_darcy.options,
        ),
        (
            name = "polymer_moving_lumped",
            problem = () -> _column(_pmma(); t_final = 5.0),
            options = (use_ale = true,),
        ),
        (
            name = "pmma_fixed_surface",
            problem = () -> _column(_pmma(); t_final = 5.0),
            options = (;),
        ),
    )
    for d in decks, energy_state in (:temperature, :enthalpy)
        problem = d.problem()
        cp = setup(problem; energy_state, d.options...)
        u = copy(cp.u0)
        du = similar(u)
        J = jac_prototype(cp.jacobian)
        id = "$(d.name)/$(energy_state)"
        push!(
            out,
            (
                id * "/residual",
                measure(
                    () -> foreach(_ -> residual!(du, u, cp.def, 0.0, cp.ws), 1:REPEATS),
                ),
            ),
        )
        push!(
            out,
            (
                id * "/jacobian",
                measure(() -> foreach(_ -> jacobian!(J, cp.jacobian, u, 0.0), 1:REPEATS)),
            ),
        )
        push!(
            out,
            (id * "/solve", measure(() -> solve(problem; energy_state, d.options...))),
        )
        println(stderr, "measured ", id)
    end
    return out
end

function read_record(path)
    record = Dict{String,Float64}()
    for line in Iterators.drop(eachline(path), 1)
        name, instructions = split(line, ",")[1:2]
        record[name] = parse(Float64, instructions)
    end
    return record
end

function main(args)
    Threads.nthreads() == 1 || error("run with --threads=1")
    BLAS.set_num_threads(1)
    option(name) = (k = findfirst(==(name), args); k === nothing ? nothing : args[k+1])
    measured = rows()
    compare = option("--compare")
    record = compare === nothing ? Dict{String,Float64}() : read_record(compare)
    @printf(
        "%-58s %14s %12s %12s %8s\n",
        "operation",
        "instructions",
        "cycles",
        "bytes",
        "ratio"
    )
    for (name, m) in measured
        ratio = haskey(record, name) ? @sprintf("%.4f", m.instructions / record[name]) : ""
        @printf(
            "%-58s %14.0f %12.0f %12d %8s\n",
            name,
            m.instructions,
            m.cycles,
            m.bytes,
            ratio
        )
    end
    path = option("--record")
    if path !== nothing
        open(path, "w") do io
            println(io, "operation,instructions,cycles,bytes")
            for (name, m) in measured
                println(io, join((name, m.instructions, m.cycles, m.bytes), ","))
            end
        end
    end
    return nothing
end

main(ARGS)
