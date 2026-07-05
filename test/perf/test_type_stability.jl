"""
CI invariant — type stability via JET.jl.

Loaded only by `runtests.jl` when JET is available; the conditional
include avoids the parse-time `@test_call` binding when JET is missing.

- `residual!` is the single core RHS entry point.
- The test scopes JET to `Pyrolysis.Problem`, where the unified state layout
  and workspace contract live.
"""

using Test
using JET
using Pyrolysis

@testset "type stability via JET" begin
    def, ws, u, du = build_baseline_eulerian()
    # Warm up.
    residual!(du, u, def, 0.0, ws)
    @test_call target_modules=(Pyrolysis.Residual,) residual!(du, u, def, 0.0, ws)
end
