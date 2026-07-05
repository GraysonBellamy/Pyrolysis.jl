"""
CI invariant — method-count regression budget.

Trait combinatorial explosion is the highest-impact pitfall here: each new
orthogonal type parameter on `residual!` multiplies the method table. All
axes are bundled inside `SimulationMode{Mesh, Rad, Trans, Chi}`. Method
count must not grow beyond what the bundling allows.
"""

using Test
using Pyrolysis
using Pyrolysis: residual!

# A single public residual entry point bundles simulation axes
# inside `SimulationMode`, keeping the method table small.
const RESIDUAL_METHOD_BUDGET = 16

@testset "method-count budget" begin
    n_inplace = length(methods(residual!))
    @info "residual! method count" n_inplace
    @test n_inplace <= RESIDUAL_METHOD_BUDGET
end
