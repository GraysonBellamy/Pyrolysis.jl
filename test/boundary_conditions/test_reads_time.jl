using Test
using Pyrolysis
using Pyrolysis.BoundaryConditions

struct _OwnCondition <: BoundaryConditions.AbstractThermalBC end

@testset "which boundary conditions read the time" begin
    of_time = t -> 1.0 + t
    of_temperature = T -> 0.9
    air = (N2 = 1.0,)
    reservoir(; P = 1e5, T = 300.0, x = air) =
        GasReservoir(P = P, T = T, mole_fractions = x)

    for bc in (
        AdiabaticBC(),
        ImpermeableBC(),
        HeatFluxBC(1e4),
        ConvectiveBC(h = 10.0, T_inf = 300.0),
        RadiativeBC(ε = 0.9, T_inf = 300.0),
        RadiativeBC(ε = of_temperature, T_inf = 300.0),
        DirichletThermalBC(400.0),
        RadiativeFluxBC(1e4),
        RadiativeFluxBC(1e4, absorptivity = 0.9),
        HeatFluxBC(1e4) + RadiativeBC(ε = of_temperature, T_inf = 300.0),
        AmbientGasBC(reservoir()),
        AmbientGasBC(reservoir(), film = MassTransferFilm(of_temperature)),
        ConvectivePressureBC(reservoir(), resistance = 1.0, flow = :molar),
        NeumannPressureBC(gradient = 0.0, T = 300.0, mole_fractions = air),
        MassFluxBC(fluxes = (N2 = 0.1,)),
        BoundaryConditionSet(
            thermal = (top = RadiativeBC(ε = of_temperature, T_inf = 300.0),),
        ),
    )
        @test !reads_time(bc)
    end

    for bc in (
        HeatFluxBC(of_time),
        ConvectiveBC(h = of_time, T_inf = 300.0),
        ConvectiveBC(h = 10.0, T_inf = of_time),
        RadiativeBC(ε = 0.9, T_inf = of_time),
        RadiativeBC(ε = 0.9, T_inf = 300.0, F = of_time),
        DirichletThermalBC(of_time),
        RadiativeFluxBC(of_time),
        RadiativeFluxBC(1e4, absorptivity = of_time),
        HeatFluxBC(1e4) + ConvectiveBC(h = 10.0, T_inf = of_time),
        AmbientGasBC(reservoir(P = of_time)),
        AmbientGasBC(reservoir(T = of_time)),
        AmbientGasBC(reservoir(x = t -> air)),
        ConvectivePressureBC(reservoir(T = of_time), resistance = 1.0, flow = :molar),
        ConvectivePressureBC(reservoir(), resistance = of_time, flow = :molar),
        NeumannPressureBC(gradient = of_time, T = 300.0, mole_fractions = air),
        NeumannPressureBC(gradient = 0.0, T = of_time, mole_fractions = air),
        NeumannPressureBC(gradient = 0.0, T = 300.0, mole_fractions = t -> air),
        MassFluxBC(fluxes = (N2 = of_time,), inflow_temperature = (N2 = 300.0,)),
        MassFluxBC(fluxes = (N2 = -0.1,), inflow_temperature = (N2 = of_time,)),
        BoundaryConditionSet(thermal = (top = HeatFluxBC(of_time),)),
        BoundaryConditionSet(mass = (top = AmbientGasBC(reservoir(P = of_time)),)),
    )
        @test reads_time(bc)
    end

    # A condition with no method of its own is taken to read it.
    @test reads_time(_OwnCondition())
end
