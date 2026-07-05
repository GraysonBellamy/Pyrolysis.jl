"""
Performance-invariant suite.

Run via:
```
julia --project=. test/perf/runtests.jl
```

The CI workflow (`.github/workflows/perf.yml`) invokes this on every PR.
Tests follow a layered budget approach.
"""

using Test

# Build helpers once so each child test can reuse the canonical
# `SimulationMode{...}` configuration builders without re-defining them.
include(joinpath(@__DIR__, "perf_helpers.jl"))

function _load_jet_for_type_stability()
    try
        @eval using JET
        return true
    catch _
        @info "JET.jl not available; type-stability tests will be skipped"
        return false
    end
end

@testset "performance invariants" begin
    include("test_zero_allocation_rhs.jl")
    include("test_jacobian_allocations.jl")
    include("test_method_count.jl")
    if _load_jet_for_type_stability()
        include("test_type_stability.jl")
    else
        @testset "type stability via JET (skipped)" begin
            @test true
        end
    end
    # AD-vs-FD reference-correctness gate (the contract verify_jacobian relies on).
    include("test_ad_jacobian_correctness.jl")
    # Backend × config differential-correctness matrix.
    include("test_verify_jacobian.jl")
end
