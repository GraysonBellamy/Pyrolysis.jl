using Test, Pyrolysis, ForwardDiff
using Pyrolysis.Materials: R_GAS
using Pyrolysis.Properties: transport_coefficients

function capillary_problem(ε; n = 24)
    mat=Material(
        name = :capillary_limit,
        gas_viscosity = 2e-5,
        gas_transport = DustyGas(
            tortuosity = 2.0,
            D_ref = 2e-5,
            T_ref = 300.0,
            P_ref = 1e5,
            capillary_regularization = ε,
        ),
        components = (
            SolidComponent(:matrix; ρ = 500.0, φ = 0.5, c = 1000.0, k = 0.2, κ = 1e-14),
            GaseousComponent(:gas; M = 0.03, c = 0.0, k = 0.03),
        ),
    )
    PyrolysisProblem(
        mesh = create_uniform_mesh(0.002, n),
        material = mat,
        T_initial = z->500.0,
        ξ_initial = [z->500.0, z->0.5*(2e5-1e5*z/0.002)*0.03/(R_GAS*500)],
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                bottom = AmbientGasBC(
                    GasReservoir(P = 2e5, T = 500.0, mole_fractions = (gas = 1.0,)),
                ),
                top = AmbientGasBC(
                    GasReservoir(P = 1e5, T = 500.0, mole_fractions = (gas = 1.0,)),
                ),
            ),
        ),
        tspan = (0.0, 2.0),
        output_times = [0.0, 1.0, 2.0],
    )
end
function capillary_flux_reference(ε)
    setprecision(256) do
        κ, φ, τ, T, M, μ, L=BigFloat(1e-14),
        BigFloat(0.5),
        BigFloat(2),
        BigFloat(500),
        BigFloat(0.03),
        BigFloat(2e-5),
        BigFloat(0.002)
        R=BigFloat(R_GAS)
        PL, PR=BigFloat(2e5), BigFloat(1e5)
        dp=sqrt(32τ/φ)*κ/sqrt(κ+BigFloat(ε))
        K=φ*dp/(3τ)*sqrt(8R*T/(big(π)*M))
        Float64(M/(R*T*L)*(K*(PL-PR)+κ*(PL^2-PR^2)/(2μ)))
    end
end
@testset "capillary regularization scale in a production solve" begin
    scales=(0.0, 1e-15, 1e-16, 1e-17)
    fluxes=Float64[]
    for ε in scales
        p=capillary_problem(ε)
        sol=solve(p; reltol = 1e-9, abstol = 1e-11, pressure_tolerance = 0.0001)
        @test completed(sol)
        N=get_boundary_exchange(sol, :top).fluxes[1, end]
        push!(fluxes, N)
        @test abs(N/capillary_flux_reference(ε)-1)<0.02 # first-order Darcy donor discretization
        @test maximum(abs.(sol.T .- 500))<1e-9
    end
    # The zero-preserving closure approaches its unregularized positive-κ
    # limit linearly in ε/κ. Compare at identical mesh and temporal accuracy.
    changes=abs.(fluxes[2:4] ./ fluxes[1] .- 1)
    @test all(changes[2:3] .< changes[1:2]/8)
    @test changes[end]<3e-4
    @test changes[1]>0.001
    # A path through zero has finite one-sided first derivatives at fixed ε.
    ε=1e-16
    φ=0.5
    τ=2.0
    diameter=κ->sqrt(32τ/φ)*κ/sqrt(κ+ε)
    @test diameter(0.0)==0
    @test ForwardDiff.derivative(diameter, 0.0)≈sqrt(32τ/(φ*ε))
end
