"""
Operator-list topology test.

Every stub operator has been replaced by a real implementation. The
gates are:

1. **Topo-sort:** the full operator list produces a valid
   `ResidualContext` whose dependency order places every provider before
   its consumer (`BCThermalOperator` before `BCConvectionOperator` for
   `:T_s`).

2. **Real-operator no-op modes** — operators with conditional behavior
   (today: `RadiationOperator` under `NO_RADIATION`) return cleanly.
"""

using Test
using Pyrolysis
using Pyrolysis.Jacobian:
    build_residual_context,
    ReactionSourceOperator,
    ConductionOperator,
    GasTransportOperator,
    RadiationOperator,
    DarcyTransportOperator,
    BCThermalOperator,
    BCMassOperator,
    BCConvectionOperator,
    ALEMeshVelocityOperator

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

@testset "Operator topology" begin

    @testset "full operator list: topo-sort places provider before consumer" begin
        ops = [
            ReactionSourceOperator(),
            ConductionOperator(),
            GasTransportOperator(),
            RadiationOperator(),
            DarcyTransportOperator(),
            BCThermalOperator(),         # writes :T_s
            BCMassOperator(),
            BCConvectionOperator(),      # reads :T_s
            ALEMeshVelocityOperator(),   # writes :w_node, :w_cell, :theta
        ]
        ctx = build_residual_context(ops)

        # All operators present
        @test length(ctx.order) == 9
        @test findfirst(o -> o isa ALEMeshVelocityOperator, ctx.order) !== nothing

        # Provider-before-consumer for :T_s (BCConvection's blended boundary
        # gradient runs through the surface Newton).
        thermal_idx = findfirst(o -> o isa BCThermalOperator, ctx.order)
        conv_idx = findfirst(o -> o isa BCConvectionOperator, ctx.order)
        @test thermal_idx !== nothing && conv_idx !== nothing
        @test thermal_idx < conv_idx

        # :T_s provider lookup
        t_s_provider = Pyrolysis.Jacobian.provider_for(ctx, :T_s)
        @test t_s_provider isa BCThermalOperator
    end

    @testset "RadiationOperator NO_RADIATION is a no-op" begin
        op = RadiationOperator()
        def, ws, _, _ = build_baseline_eulerian()
        material = def.material
        dummy = zeros(ws.layout.n_cells)

        # NO_RADIATION: no-op, returns nothing.
        @test Pyrolysis.Jacobian.apply!(
            dummy,
            op,
            ws.state_cache,
            ws.mesh,
            ws.prop_cache,
            material,
            0.0,
            Pyrolysis.Materials.NO_RADIATION,
            true,
        ) === nothing
    end
end
