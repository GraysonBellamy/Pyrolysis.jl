"""
`StructuredCoupling` unit tests.

Checks `densify!` of the three embedded couplings, and `compose_into!` of a
`PrefixSumCoupling`, against dense references:

- [`EmbeddedNodePrefixSumCoupling`](@ref) — the ALE cumulative mesh-velocity
  coupling at the node rows.
- [`EmbeddedScaledSuffixSumCoupling`](@ref) — the Beer–Lambert suffix sum
  with a per-row scale.
- [`EmbeddedRank1Coupling`](@ref) — `u · vᵀ` (the WithChi θ_A channels).

The assembled couplings are also gated end-to-end against DenseAD by the
operator and backend tests; this file pins the algebra at the block level.
"""

using Test
using SparseArrays
using Pyrolysis.Jacobian:
    PrefixSumCoupling,
    EmbeddedNodePrefixSumCoupling,
    EmbeddedScaledSuffixSumCoupling,
    EmbeddedRank1Coupling,
    StructuredJacobian,
    compose_into!,
    densify!

@testset "StructuredCoupling blocks" begin
    v = [0.5, -1.0, 2.0, 0.25]
    n_global = 12
    row_off, col_off = 3, 1

    @testset "EmbeddedNodePrefixSumCoupling" begin
        n_nodes = length(v) + 1
        M = zeros(n_global, n_global)
        densify!(M, EmbeddedNodePrefixSumCoupling(v, row_off, col_off, n_nodes, n_global))
        ref = zeros(n_global, n_global)
        for k = 1:n_nodes, j = 1:length(v)
            j <= k - 1 && (ref[row_off+k, col_off+j] += v[j])
        end
        @test M == ref
        @test all(iszero, M[row_off+1, :])   # node 1 is the fixed bottom node
    end

    @testset "EmbeddedScaledSuffixSumCoupling" begin
        row_scale = [1.0, -2.0, 3.0, 0.5]
        N = length(v)
        M = zeros(n_global, n_global)
        densify!(
            M,
            EmbeddedScaledSuffixSumCoupling(v, row_scale, row_off, col_off, N, n_global),
        )
        ref = zeros(n_global, n_global)
        for i = 1:N, j = 1:N
            j > i && (ref[row_off+i, col_off+j] += row_scale[i] * v[j])
        end
        @test M == ref
        @test_throws DimensionMismatch EmbeddedScaledSuffixSumCoupling(
            v,
            row_scale[1:3],
            row_off,
            col_off,
            N,
            n_global,
        )
    end

    @testset "EmbeddedRank1Coupling" begin
        u = [1.0, -2.0, 0.0]
        w = [3.0, 0.5]
        M = zeros(n_global, n_global)
        densify!(M, EmbeddedRank1Coupling(u, w, row_off, col_off, n_global))
        ref = zeros(n_global, n_global)
        for i = 1:length(u), j = 1:length(w)
            ref[row_off+i, col_off+j] = u[i] * w[j]
        end
        @test M == ref
    end

    @testset "compose_into! embeds a PrefixSumCoupling at the node rows" begin
        jac = StructuredJacobian(spzeros(Float64, n_global, n_global))
        n_nodes = length(v) + 1
        compose_into!(jac, PrefixSumCoupling(v), row_off, col_off, n_nodes, n_global)
        @test length(jac.structured) == 1
        c = only(jac.structured)
        @test c isa EmbeddedNodePrefixSumCoupling
        @test c.v === v   # aliases the factors, no copy
        M = zeros(n_global, n_global)
        densify!(M, c)
        ref = zeros(n_global, n_global)
        densify!(ref, EmbeddedNodePrefixSumCoupling(v, row_off, col_off, n_nodes, n_global))
        @test M == ref
    end
end
