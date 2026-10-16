# Jacobian interface skeleton.
#
# Verifies that the src/Jacobian/ tree exposes its public surface
# correctly: types are visible, JacobianSpec construction is type-stable,
# and the dispatcher errors cleanly for backends that are not implemented.

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis: Jacobian
using Pyrolysis.Jacobian:
    JacobianBackend,
    JacobianSpec,
    compute_jacobian!,
    _compute_jacobian!,
    _build_jacobian_spec,
    FastSparsityMap

@testset "Jacobian interface skeleton" begin
    @testset "public surface" begin
        @test isabstracttype(JacobianBackend)
        @test isconcretetype(FastSparsityMap)
        @test JacobianSpec isa UnionAll
        @test compute_jacobian! isa Function
        # `_compute_jacobian!` and `_build_jacobian_spec` are the
        # extension hooks. They must exist and be Functions so backend
        # files can specialize them.
        @test _compute_jacobian! isa Function
        @test _build_jacobian_spec isa Function
    end

    # A tiny problem for the (def, ws) constructor.
    mat = Material(
        name = :interface_test,
        components = (
            SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 0.2),
            GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5),
        ),
        reactions = (Reaction(:decomp, 1 => 2, A = 1e10, E = 150e3, h = -1e5),),
    )
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 4, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
        ),
        tspan = (0.0, 1.0),
    )
    def, ws, _ = initialize_state(problem)

    @testset "JacobianSpec construction is type-stable" begin
        # A locally-defined dummy backend exercises the interface without a
        # real backend: a method on `_build_jacobian_spec` builds a trivial
        # spec, and `JacobianSpec(::DummyBackend, def, ws)` must be
        # `@inferred`-stable.
        struct DummyBackend <: JacobianBackend end
        struct DummyCache end
        function Pyrolysis.Jacobian._build_jacobian_spec(
            ::DummyBackend,
            ::Pyrolysis.Problem.ProblemDef,
            ::Pyrolysis.Problem.Workspace;
            kwargs...,
        )
            return JacobianSpec(DummyBackend(), nothing, nothing, DummyCache())
        end

        spec = @inferred JacobianSpec(DummyBackend(), def, ws)
        @test spec.backend isa DummyBackend
        @test spec.cache isa DummyCache
        @test spec.sparsity_phys === nothing
        @test spec.sparsity_geom === nothing
    end

    @testset "compute_jacobian! errors cleanly on unimplemented backend" begin
        # A backend that has not provided `_compute_jacobian!` must hit the
        # fallback that explains why. No `MethodError` — a clean message
        # pointing the user at the available backends.
        struct UnshippedBackend <: JacobianBackend end
        spec = JacobianSpec(UnshippedBackend(), nothing, nothing, nothing)
        J = zeros(2, 2)
        err = try
            compute_jacobian!(J, spec, Float64[], nothing, 0.0, nothing)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("UnshippedBackend", err.msg)
        @test occursin("no `_compute_jacobian!`", err.msg)
    end

    @testset "JacobianSpec(backend, def, ws) errors cleanly when backend ships no _build_jacobian_spec" begin
        struct UnconstructibleBackend <: JacobianBackend end
        err = try
            JacobianSpec(UnconstructibleBackend(), def, ws)
            nothing
        catch e
            e
        end
        @test err isa ErrorException
        @test occursin("UnconstructibleBackend", err.msg)
        @test occursin("no `_build_jacobian_spec`", err.msg)
    end

    @testset "FastSparsityMap O(1) sparse-write helper" begin
        # Round-trip: build a FastSparsityMap from a sparse matrix, scatter
        # values through fast_sparse_add!, and check they land in nzval.
        J = sparse([1, 2, 3], [1, 2, 3], [1.0, 2.0, 3.0], 3, 3)
        fsm = FastSparsityMap(J)
        Pyrolysis.Jacobian.fast_sparse_add!(J, fsm, 2, 2, 10.0)
        @test J[2, 2] == 12.0
        # Out-of-pattern writes are silently dropped (existing contract).
        Pyrolysis.Jacobian.fast_sparse_add!(J, fsm, 1, 2, 99.0)
        @test J[1, 2] == 0.0
    end
end
