"""
CI invariant — zero-allocation RHS.

Baseline + the three additional `SimulationMode` configurations
(`{ALE, NoRadiation, Fickian, NoChi}`,
`{Eulerian, SurfaceAbsorption, Fickian, NoChi}`,
`{ALE, BeerLambert, DarcyFick, NoChi}`). Each must allocate zero bytes
after warm-up — proves the trait-dispatched no-ops compile away and
the ALE/radiation/Darcy paths are alloc-free in steady state.
"""

using Test
using Pyrolysis

include(joinpath(@__DIR__, "perf_helpers.jl"))

const RHS_ALLOC_BUDGET = 0

# Run residual! once warm and check the allocation count of the *next* call.
function _check_zero_alloc(builder, label)
    @testset "$label" begin
        def, ws, u, du = builder()
        residual!(du, u, def, 0.0, ws)
        residual!(du, u, def, 0.0, ws)
        bytes = @allocated residual!(du, u, def, 0.0, ws)
        @test bytes <= RHS_ALLOC_BUDGET
    end
end

@testset "zero-allocation RHS" begin
    _check_zero_alloc(
        build_baseline_eulerian,
        "baseline {Eulerian, NoRadiation, Fickian, NoChi}",
    )
    _check_zero_alloc(build_ale_no_radiation, "{ALE, NoRadiation, Fickian, NoChi}")
    _check_zero_alloc(
        build_eulerian_surface_absorption,
        "{Eulerian, SurfaceAbsorption, Fickian, NoChi}",
    )
    _check_zero_alloc(build_ale_beer_lambert_darcy, "{ALE, BeerLambert, DarcyFick, NoChi}")
end
