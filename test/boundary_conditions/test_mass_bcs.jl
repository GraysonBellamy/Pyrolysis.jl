using Test, Pyrolysis, StaticArrays, ForwardDiff
using Pyrolysis.BoundaryConditions: validate_gas_boundary, _reservoir_values, mass_film_flux

const material=Material(
    name = :boundary_validation,
    components = (
        SolidComponent(
            :matrix;
            ρ = 500.0,
            ρ_intrinsic = 1000.0,
            c = 1000.0,
            k = 0.2,
            λ = 1e-5,
        ),
        GaseousComponent(:a; M = 0.028, c = 1000.0, k = 0.03),
        GaseousComponent(:b; M = 0.044, c = 1000.0, k = 0.03),
    ),
)
@testset "complete reservoir and named flux data" begin
    @test_throws ArgumentError GasReservoir(P = 1e5, T = 300.0)
    @test_throws ArgumentError GasReservoir(P = -1.0, T = 300.0)
    @test_throws ArgumentError GasReservoir(P = 0.0, T = 0.0)
    @test_throws ArgumentError GasReservoir(P = 1e5, T = 300.0, mole_fractions = (a = 0.9,))
    @test_throws ArgumentError GasReservoir(
        P = 1e5,
        T = 300.0,
        mole_fractions = (a = 1.1, b = -0.1),
    )
    vacuum=GasReservoir(P = 0.0, T = 300.0)
    @test _reservoir_values(vacuum, material, 0.0).x==(0.0, 0.0, 0.0)
    timed=GasReservoir(P = t->1e5*(1+t), T = t->300+t, mole_fractions = t->(a = t, b = 1-t))
    r=_reservoir_values(timed, material, 0.3)
    @test r.P==1.3e5 && r.T==300.3 && r.x==(0.0, 0.3, 0.7)
    @test_throws ArgumentError _reservoir_values(timed, material, 2.0)
    @test_throws ArgumentError validate_gas_boundary(
        AmbientGasBC(GasReservoir(P = 1e5, T = 300.0, mole_fractions = (other = 1.0,))),
        material,
        0.0,
    )
    @test_throws ArgumentError MassFluxBC(fluxes = (a = -1.0,))
    @test_throws ArgumentError MassFluxBC(fluxes = (a = t->t,))
    @test_throws ArgumentError MassFluxBC(
        fluxes = (a = -1.0,),
        inflow_temperature = (a = NaN,),
    )
    @test_throws ArgumentError validate_gas_boundary(
        MassFluxBC(fluxes = (matrix = 1.0,)),
        material,
        0.0,
    )
    @test isnothing(
        validate_gas_boundary(
            MassFluxBC(fluxes = (a = -0.01, b = 0.02), inflow_temperature = (a = 500.0,)),
            material,
            0.0,
        ),
    )
    @test_throws ArgumentError ConvectivePressureBC(vacuum; resistance = -1.0, flow = :mass)
    @test_throws ArgumentError ConvectivePressureBC(
        vacuum;
        resistance = 1.0,
        flow = :ambiguous,
    )
    @test_throws ArgumentError NeumannPressureBC(
        gradient = NaN,
        T = 300.0,
        mole_fractions = (a = 1.0,),
    )
end
@testset "independently integrated signed mass film" begin
    ys=SVector(0.2, 0.3, 0.5)
    ya=SVector(0.6, 0.4, 0.0)
    for m in (-1000.0, -1.0, -1e-10, 0.0, 1e-10, 1.0, 1000.0), g in (1e-8, 0.1, 1.0, 1e8)
        n=mass_film_flux(ys, ya, m, g)
        expected=setprecision(256) do
            # Reference from the original data before dependent arithmetic.
            M=BigFloat(m)
            G=BigFloat(g)
            B=M/G
            rate=M==0 ? G : M/expm1(B)
            Float64.(BigFloat.(ys)*M+rate*(BigFloat.(ys)-BigFloat.(ya)))
        end
        @test n≈expected rtol=4e-14 atol=2e-14*max(abs(m), g)
        @test sum(n)≈m atol=1e-14*max(abs(m), g)
    end
    @test mass_film_flux(ys, ya, 2.0, 0.0)==2ys
    @test mass_film_flux(ys, ya, -2.0, 0.0) == -2ya
    @test mass_film_flux(ys, ya, 0.0, 0.1)≈0.1*(ys-ya)
    @test all(isfinite, ForwardDiff.jacobian(x->mass_film_flux(ys, ya, x[1], 0.1), [0.0]))
end
