# Geometry derivatives via AD.
#
# Acceptance criterion: `apply_geometry_jacobian!` produces the geometry block
# to `rtol=1e-10` relative to a full-state ForwardDiff reference
# (`DenseAD_ForwardDiff`).

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis:
    JacobianSpec, DenseAD_ForwardDiff, compute_jacobian!, build_ad_workspace, residual!
using Pyrolysis.Jacobian: GeometryADCache, apply_geometry_jacobian!
const _state_length = Pyrolysis.Problem.state_length

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# Reference: DenseAD_ForwardDiff via the spec API. Returns the
# z-column subset of the full Jacobian.
function _reference_geometry_columns(def, ws, u, z_off, n_nodes)
    n_total = length(u)
    spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    J_full = zeros(Float64, n_total, n_total)
    compute_jacobian!(J_full, spec, u, def, 0.0, ws)
    z_range = (z_off+1):(z_off+n_nodes)
    return J_full[:, z_range]
end

@testset "apply_geometry_jacobian! vs ForwardDiff baseline" begin
    @testset "ALE, NoRadiation, Fickian" begin
        def, ws, u, _ = build_ale_no_radiation()
        cache = GeometryADCache(def, ws)
        n_total = _state_length(ws.layout)
        n_nodes = ws.layout.n_nodes
        z_off = ws.layout.n_cells * (1 + 2)  # NC = 2
        z_range = (z_off+1):(z_off+n_nodes)

        J_ref_geom = _reference_geometry_columns(def, ws, u, z_off, n_nodes)

        J = zeros(Float64, n_total, n_total)
        apply_geometry_jacobian!(J, cache, u, def, 0.0)

        J_test_geom = J[:, z_range]

        @test size(J_test_geom) == size(J_ref_geom)
        # Tight tolerance — both are ForwardDiff over the same pure residual,
        # differing only in whether non-z slots are perturbed.
        @test isapprox(J_test_geom, J_ref_geom; rtol = 1e-10, atol = 1e-12)
    end

    @testset "ALE, BeerLambert, DarcyFick (wood-like) — geometry pass alone" begin
        def, ws, u, _ = build_ale_beer_lambert_darcy()
        cache = GeometryADCache(def, ws)
        n_total = _state_length(ws.layout)
        n_nodes = ws.layout.n_nodes
        z_off = ws.layout.n_cells * (1 + 2)
        z_range = (z_off+1):(z_off+n_nodes)

        J = zeros(Float64, n_total, n_total)
        apply_geometry_jacobian!(J, cache, u, def, 0.0)
        # Sanity: the geometry block is non-trivial.
        @test count(!iszero, J[:, z_range]) >= n_nodes
    end

    @testset "Eulerian configurations have no z block" begin
        def, ws, u, _ = build_baseline_eulerian()
        @test ws.layout.n_cells * (1 + 2) == _state_length(ws.layout)
    end
end
