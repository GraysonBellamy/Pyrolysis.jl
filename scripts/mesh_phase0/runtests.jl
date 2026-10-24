# Initial executable contracts for docs/Mesh_Framework_Design_2026-10-17.md.
# Run from any directory:
#   julia --project=/path/to/Pyrolysis.jl --startup-file=no /path/to/this/runtests.jl
#   julia --project=. scripts/mesh_phase0/runtests.jl --case=closure
#
# This is a prototype harness, NOT the complete Phase 0 gate. It does not yet
# audit production decks or color counts, validate production property laws,
# or integrate the future production core. Timings are prototype measurements.
# Known contract failures are @test_broken, explained in the output, and make
# the command exit nonzero. They must not be counted as a passing design gate.

const CASES = [
    "kinetics",
    "caloric",
    "closure",
    "transfer",
    "fields",
    "energy",
    "integration",
    "nonlinear_caloric",
    "caloric_comparison",
    "solvers",
    "restart_callbacks",
    "restarts",
]
if ARGS == ["--help"] || ARGS == ["--list"]
    println("Initial mesh Phase 0 contracts; not a production or complete Phase 0 gate.")
    println("Usage: julia --project=. scripts/mesh_phase0/runtests.jl [--case=NAME]")
    println("Cases: ", join(CASES, ", "))
    println("Exit status is nonzero for a failing test OR an unresolved contract.")
    exit(0)
end
selected = if isempty(ARGS)
    CASES
elseif length(ARGS) == 1 && startswith(only(ARGS), "--case=")
    name = split(only(ARGS), '='; limit = 2)[2]
    name in CASES || error("unknown case: $name (use --list)")
    [name]
else
    error("unexpected arguments (use --help)")
end

using Test, LinearAlgebra, ForwardDiff, Printf
using OrdinaryDiffEq, OrdinaryDiffEqSDIRK
import SciMLBase
include("prototypes.jl")
using .MeshPhase0

const unresolved = String[]
function metric(name, value)
    println("  ", name, " = ", value)
    flush(stdout)
end

# Pin the small prototype linear solves for reproducible CPU timings. This
# affects only this harness process, not the production solver configuration.
BLAS.set_num_threads(1)
println("Mesh Phase 0 initial contracts")
println("Julia ", VERSION, "; OrdinaryDiffEqSDIRK ", pkgversion(OrdinaryDiffEqSDIRK))
println("Julia threads ", Threads.nthreads(), "; BLAS threads ", BLAS.get_num_threads())
flush(stdout)
@testset "mesh Phase 0 prototype contracts" begin
    for name in selected
        include("test_$(name).jl")
    end
end

if !isempty(unresolved)
    println("\nUNRESOLVED CONTRACTS (these are not passing gate results):")
    foreach(s -> println("  - ", s), unresolved)
    exit(1)
end
println("\nSelected prototype contracts passed. Full Phase 0 remains incomplete.")
