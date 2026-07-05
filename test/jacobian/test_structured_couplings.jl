"""
`StructuredCoupling` primitives unit tests.

Verifies `apply!` / `apply_inverse!` / `densify!` against dense
references for the two concrete couplings:

- [`PrefixSumCoupling`](@ref) — `L · diag(v)`, the prefix-sum operator
  that captures the ALE cumulative mesh-velocity coupling.
- [`Rank1Coupling`](@ref) — `u · vᵀ`, used for column-global modes
  (today: θ_A back-reaction in WithChi).

Both are wired into the real ALE assembly; this test pins the
math at the primitive level so wiring bugs don't get blamed on
algebra bugs (and vice versa).
"""

using Test
using LinearAlgebra
using Pyrolysis.Jacobian: PrefixSumCoupling, SuffixSumCoupling, Rank1Coupling
using Pyrolysis.Jacobian: apply!, densify!
import Pyrolysis.Jacobian: apply_inverse!

@testset "StructuredCoupling primitives" begin

    @testset "PrefixSumCoupling" begin
        N = 8
        v = collect(range(0.5, 4.0, length = N))
        c = PrefixSumCoupling(v)

        # Dense reference: L * diag(v) where L is exclusive lower-triangular.
        # (L · diag(v))[i, j] = v[j] if j ≤ i-1 else 0.
        M_ref = zeros(N, N)
        for j = 1:N, i = 1:N
            if j <= i - 1
                M_ref[i, j] = v[j]
            end
        end

        # `densify!` matches the dense reference.
        M = zeros(N, N)
        densify!(M, c)
        @test M ≈ M_ref

        # `apply!` matches the dense reference's matvec.
        x = randn(N)
        y = zeros(N)
        apply!(y, c, x)
        @test y ≈ M_ref * x

        # `apply_inverse!` round-trips for x in the column space — i.e.
        # any vector `M_ref * z` should map back to a vector y satisfying
        # `M_ref * y[2:end] = M_ref * z`. The first row is identically zero
        # (consistent system requires x[1] == 0); we ignore y[1].
        z = randn(N)
        x2 = M_ref * z
        y2 = zeros(N)
        apply_inverse!(y2, c, x2)
        @test M_ref * y2 ≈ x2 atol = 1e-10
    end

    @testset "SuffixSumCoupling" begin
        N = 8
        v = collect(range(0.5, 4.0, length = N))
        c = SuffixSumCoupling(v)

        # Dense reference: U * diag(v) where U is exclusive upper-triangular.
        # (U · diag(v))[i, j] = v[j] if j ≥ i+1 else 0.
        M_ref = zeros(N, N)
        for j = 1:N, i = 1:N
            if j >= i + 1
                M_ref[i, j] = v[j]
            end
        end

        # `densify!` matches the dense reference.
        M = zeros(N, N)
        densify!(M, c)
        @test M ≈ M_ref

        # `apply!` matches the dense reference's matvec.
        x = randn(N)
        y = zeros(N)
        apply!(y, c, x)
        @test y ≈ M_ref * x

        # `apply!` accumulates additively (matches docstring contract).
        y .= 1.0
        apply!(y, c, x)
        @test y ≈ ones(N) .+ M_ref * x

        # `apply_inverse!` round-trips for x in the column space — i.e.
        # any vector `M_ref * z` should map back to a vector y satisfying
        # `M_ref * y = M_ref * z`. The last row is identically zero
        # (consistent system requires x[N] == 0); the first index is in
        # the right-nullspace and is set to zero by convention.
        z = randn(N)
        x2 = M_ref * z
        y2 = zeros(N)
        apply_inverse!(y2, c, x2)
        @test M_ref * y2 ≈ x2 atol = 1e-10
    end

    @testset "Rank1Coupling" begin
        N = 6
        u = randn(N)
        v = randn(N)
        c = Rank1Coupling(u, v)

        M_ref = u * v'

        # `densify!`.
        M = zeros(N, N)
        densify!(M, c)
        @test M ≈ M_ref

        # `apply!` is one inner product + scaled add.
        x = randn(N)
        y = zeros(N)
        apply!(y, c, x)
        @test y ≈ M_ref * x

        # `apply!` accumulates additively (matches the docstring contract).
        y .= 1.0
        apply!(y, c, x)
        @test y ≈ ones(N) .+ M_ref * x
    end
end
