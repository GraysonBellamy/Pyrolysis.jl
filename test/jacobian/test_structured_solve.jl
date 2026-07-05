"""
`StructuredSolve` strategy unit tests.

Verifies that the Richardson-with-`A⁻¹`-preconditioner iteration in
[`StructuredSolve`](@ref) converges to the same answer as
[`DirectSolve`](@ref) (which densifies `A + Σⱼ Cⱼ` and factors with
`\\`). Both strategies should solve `(A + Σⱼ Cⱼ) y = b` exactly modulo
floating-point round-off.

This is a **strategy** unit test — it builds synthetic
[`StructuredJacobian`](@ref) instances directly, no operator wiring.
The operator-wired integration test lives alongside the ALE
mesh-velocity operator.
"""

using Test
using LinearAlgebra
using SparseArrays
using Pyrolysis
using Pyrolysis.Jacobian:
    StructuredJacobian,
    PrefixSumCoupling,
    SuffixSumCoupling,
    Rank1Coupling,
    DirectSolve,
    StructuredSolve,
    ApproxSolve,
    MatrixFreeSolve
import Pyrolysis.Jacobian: ldiv!

@testset "StructuredSolve" begin

    @testset "Empty structured list — degenerates to A⁻¹ b" begin
        # No structured corrections: StructuredSolve must equal A \ b.
        N = 10
        A =
            sparse(Diagonal(2.0 .+ 0.1 .* (1:N))) +
            sparse(spdiagm(-1 => -ones(N - 1), 1 => -ones(N - 1)))
        jac = StructuredJacobian(:test, A)
        b = randn(N)

        y_struct = zeros(N)
        ldiv!(y_struct, StructuredSolve(), jac, b)
        y_direct = zeros(N)
        ldiv!(y_direct, DirectSolve(), jac, b)

        @test y_struct ≈ y_direct rtol = 1e-12
    end

    @testset "A + PrefixSumCoupling — Richardson converges" begin
        # Diagonal-dominant A so (A⁻¹·C) has spectral radius < 1 and the
        # Richardson iteration converges.
        N = 12
        A = sparse(Diagonal(10.0 .+ 0.5 .* (1:N)))
        v = 0.05 .* collect(range(0.5, 1.5, length = N))
        c = PrefixSumCoupling(v)

        jac = StructuredJacobian(:test, copy(A), [c])

        # Generate a known solution and the corresponding rhs.
        x_true = randn(N)
        b = zeros(N)
        mul!(b, A, x_true)
        # additive apply: b += C · x_true
        Pyrolysis.Jacobian.apply!(b, c, x_true)

        # DirectSolve: the truth (assembles dense A+C and factors).
        y_direct = zeros(N)
        ldiv!(y_direct, DirectSolve(), jac, b)
        @test y_direct ≈ x_true rtol = 1e-10

        # StructuredSolve: same answer.
        y_struct = zeros(N)
        ldiv!(y_struct, StructuredSolve(), jac, b)
        @test y_struct ≈ y_direct rtol = 1e-10
    end

    @testset "A + SuffixSumCoupling" begin
        N = 12
        A = sparse(Diagonal(10.0 .+ 0.5 .* (1:N)))
        v = 0.05 .* collect(range(0.5, 1.5, length = N))
        c = SuffixSumCoupling(v)

        jac = StructuredJacobian(:test, copy(A), [c])

        x_true = randn(N)
        b = zeros(N)
        mul!(b, A, x_true)
        Pyrolysis.Jacobian.apply!(b, c, x_true)

        y_direct = zeros(N)
        ldiv!(y_direct, DirectSolve(), jac, b)
        y_struct = zeros(N)
        ldiv!(y_struct, StructuredSolve(), jac, b)

        @test y_struct ≈ y_direct rtol = 1e-10
        @test y_struct ≈ x_true rtol = 1e-10
    end

    @testset "A + PrefixSum + Rank1 (WithChi-style)" begin
        N = 16
        A = sparse(Diagonal(10.0 .+ 0.5 .* (1:N)))
        v_pref = 0.03 .* collect(range(0.5, 1.5, length = N))
        c_pref = PrefixSumCoupling(v_pref)
        u_r = 0.01 .* randn(N)
        v_r = randn(N)
        c_r1 = Rank1Coupling(u_r, v_r)

        jac = StructuredJacobian(:test, copy(A), [c_pref, c_r1])

        x_true = randn(N)
        b = zeros(N)
        mul!(b, A, x_true)
        Pyrolysis.Jacobian.apply!(b, c_pref, x_true)
        Pyrolysis.Jacobian.apply!(b, c_r1, x_true)

        y_direct = zeros(N)
        ldiv!(y_direct, DirectSolve(), jac, b)
        y_struct = zeros(N)
        ldiv!(y_struct, StructuredSolve(), jac, b)

        @test y_struct ≈ y_direct rtol = 1e-10
        @test y_struct ≈ x_true rtol = 1e-10
    end

    @testset "MatrixFreeSolve matches DirectSolve" begin
        N = 12
        A = sparse(Diagonal(8.0 .+ 0.25 .* (1:N)))
        c_pref = PrefixSumCoupling(0.02 .* collect(range(0.5, 1.5, length = N)))
        c_suff = SuffixSumCoupling(0.015 .* collect(range(1.5, 0.5, length = N)))
        c_r1 = Rank1Coupling(0.005 .* randn(N), randn(N))
        jac = StructuredJacobian(:test, copy(A), [c_pref, c_suff, c_r1])
        b = randn(N)

        y_direct = zeros(N)
        ldiv!(y_direct, DirectSolve(), jac, b)
        y_matrix_free = zeros(N)
        ldiv!(y_matrix_free, MatrixFreeSolve(), jac, b; tol = 1e-12)

        @test y_matrix_free ≈ y_direct rtol = 1e-10
    end

    @testset "StructuredSolve reports non-convergence" begin
        N = 12
        A = sparse(Diagonal(10.0 .+ 0.5 .* (1:N)))
        v = 0.05 .* collect(range(0.5, 1.5, length = N))
        c = PrefixSumCoupling(v)
        jac = StructuredJacobian(:test, copy(A), [c])
        b = ones(N)
        y = zeros(N)

        @test_throws ErrorException ldiv!(
            y,
            StructuredSolve(),
            jac,
            b;
            tol = 1e-16,
            max_iter = 1,
        )
    end

    @testset "ApproxSolve drops the structured terms" begin
        # When ApproxSolve drops the coupling, the result should match
        # A⁻¹ b — NOT (A + C)⁻¹ b. This pins the documented behaviour
        # (old reduced-ALE parity) so it doesn't drift.
        N = 10
        A = sparse(Diagonal(10.0 .+ 0.5 .* (1:N)))
        v = 0.05 .* collect(range(0.5, 1.5, length = N))
        c = PrefixSumCoupling(v)

        jac = StructuredJacobian(:test, copy(A), [c])
        b = randn(N)

        y_approx = zeros(N)
        ldiv!(y_approx, ApproxSolve(), jac, b)
        @test y_approx ≈ A \ b rtol = 1e-12

        # And ApproxSolve should NOT match DirectSolve when the coupling
        # is non-zero (sanity check that the test actually exercises the
        # coupling).
        y_direct = zeros(N)
        ldiv!(y_direct, DirectSolve(), jac, b)
        @test !isapprox(y_approx, y_direct; rtol = 1e-3)
    end
end
