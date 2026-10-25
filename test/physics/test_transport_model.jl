using Pyrolysis.Integration: validate_transport_support
using Test, Pyrolysis, ForwardDiff
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Properties: transport_matrix_state, transport_coefficients, transport_domain
using Pyrolysis.Physics: pore_state

@testset "matrix occupancy and pore state" begin
    a = SolidComponent(
        :a;
        ρ = 500.0,
        ρ_intrinsic = 1000.0,
        c = 1000.0,
        k = 0.2,
        λ = 2e-5,
        κ = 1e-14,
    )
    b = SolidComponent(
        :b;
        ρ = 200.0,
        ρ_intrinsic = 1000.0,
        c = 1000.0,
        k = 0.2,
        γ = 0.5,
        λ = 6e-5,
        κ = 3e-14,
    )
    water = LiquidComponent(:water; ρ = 1000.0, c = 4000.0, k = 0.5, γ = 0.0)
    gas = GaseousComponent(:gas; M = 0.028, c = 1000.0, k = 0.03, γ = 1.0)
    mat = Material(
        name = :test,
        components = (a, b, water, gas),
        mixing = MaterialMixing(λ = PARALLEL, κ = PARALLEL),
    )
    for ξg in (0.0, 1e-10, 1.0, 1e5)
        ξ = (250.0, 100.0, 10.0, ξg)
        q, w = transport_matrix_state(ξ, 500.0, mat)
        @test q ≈ 0.75
        @test w == (2/3, 1/3, 0.0, 0.0)
        co = transport_coefficients(ξ, 500.0, mat, 0.64)
        @test co.λ ≈ (2*2e-5+6e-5)/3
        @test co.κ ≈ (2*1e-14+3e-14)/3
        ps = pore_state(ξ, 500.0, mat)
        @test ps.φ ≈ 0.64
        @test ps.p[4] ≈ ξg*R_GAS*500/(0.028*0.64)
        @test ps.P == sum(ps.p)
    end
    @test_throws ArgumentError Material(
        name = :missing,
        components = (SolidComponent(:s; ρ = 1.0, c = 1.0, k = 1.0), gas),
    )
    @test_throws ArgumentError LiquidComponent(
        :ignored;
        ρ = 1.0,
        c = 1.0,
        k = 1.0,
        γ = 0.0,
        λ = 1.0,
    )
    @test_throws MethodError GaseousComponent(
        :ignored;
        M = 0.028,
        c = 1.0,
        k = 1.0,
        λ = 1.0,
    )
    @test_throws ArgumentError DustyGas(
        tortuosity = 0.5,
        D_ref = 1e-5,
        T_ref = 300.0,
        P_ref = 1e5,
    )
end

@testset "matrix endpoint accounts for remaining skeleton" begin
    water = LiquidComponent(:water; ρ = 1000.0, c = 4000.0, k = 0.5, γ = 0.0)
    gas = GaseousComponent(:gas; M = 0.028, c = 1000.0, k = 0.03)
    endpoint = LumpedEndpoint(λ = 1e-4, κ = 0.0, max_skeletal_occupancy = 0.5)
    mat = Material(
        name = :occupied_endpoint,
        components = (water, gas),
        gas_transport = LumpedTransfer(endpoint = endpoint),
    )
    ps = pore_state((200.0, 1.0), 500.0, mat)
    @test ps.φ == 0.8
    co = transport_coefficients((200.0, 1.0), 500.0, mat, ps.φ)
    @test co.λ == 1e-4
    @test co.κ == 0
    @test !transport_domain((600.0, 1.0), 500.0, mat).admissible
end

@testset "reachable exact-series support" begin
    a=SolidComponent(:open; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, λ = 1e-5)
    b=SolidComponent(:barrier; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, λ = 0.0)
    precursor=SolidComponent(:precursor; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, λ = 1e-5)
    gas=GaseousComponent(:gas; M = 0.03, c = 1000.0, k = 0.03)
    rxn=Reaction(:make_barrier, 3=>2; A = 1.0, E = 0.0, h = 0.0)
    mat=Material(
        name = :support,
        components = (a, b, precursor, gas),
        reactions = (rxn,),
        mixing = MaterialMixing(λ = MixingSpec(SERIES)),
    )
    # No precursor exists or can be produced: no barrier can appear.
    @test isnothing(validate_transport_support(mat, reshape([500.0, 0.0, 0.0, 0.0], 4, 1)))
    @test_throws ArgumentError validate_transport_support(
        mat,
        reshape([400.0, 0.0, 100.0, 0.0], 4, 1),
    )
    closed=Material(
        name = :closed,
        components = (b, gas),
        reactions = (Reaction(:consume, 1=>2; A = 1.0, E = 0.0, h = 0.0),),
        mixing = MaterialMixing(λ = MixingSpec(SERIES)),
    )
    @test isnothing(validate_transport_support(closed, reshape([500.0, 0.0], 2, 1)))
    parallel=Material(
        name = :parallel,
        components = mat.components,
        reactions = mat.reactions,
        mixing = MaterialMixing(λ = MixingSpec(WEIGHTED, 1.0)),
    )
    @test isnothing(
        validate_transport_support(parallel, reshape([400.0, 0.0, 100.0, 0.0], 4, 1)),
    )
end
