"""
A solve that stops while the domain test keeps rejecting a gas-bearing cell
says so: the outcome is recorded and `solve` warns with the cell and the
quantity that left the domain, even when the step reaches `dtmin` before the
test gives up.

The deck converts a porous solid on a fixed mesh into a product whose skeleton
does not fit the cell: the top cell's pores close (`:pore_closure`).
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: R_GAS

function _closing(; A = 1e5)
    solid = SolidComponent(
        :solid;
        φ = 0.3,
        λ = 1e-5,
        κ = 1e-12,
        ρ = 1000.0,
        c = 1000.0,
        k = 0.2,
        ε = 0.9,
    )
    # Dense enough that the converted cell's skeleton outgrows the cell.
    foam = SolidComponent(
        :foam;
        φ = 0.0,
        λ = 1e-5,
        κ = 1e-12,
        ρ = 500.0,
        c = 1000.0,
        k = 0.2,
        ε = 0.9,
    )
    n2 = GaseousComponent(:N2; M = 0.028014, c = 1100.0, k = 0.04)
    material = Material(
        name = :closing,
        components = (solid, foam, n2),
        reactions = (Reaction(:grow, 1 => 2; A, E = 5e4, h = 0.0),),
    )
    purge = GasReservoir(P = 101325.0, T = 360.0, mole_fractions = (N2 = 1.0,))
    return PyrolysisProblem(
        mesh = uniform_mesh(0.005, 8),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(5e4) + ConvectiveBC(h = 11.2, T_inf = 360.0),
                bottom = AdiabaticBC(),
            ),
            mass = (top = AmbientGasBC(purge), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 360.0,
        ξ_initial = [
            z -> 1000.0,
            z -> 0.0,
            z -> 0.3 * 101325.0 * 0.028014 / (R_GAS * 360.0),
        ],
        tspan = (0.0, 20.0),
    )
end

@testset "a domain outcome names its cell" begin
    sol = @test_logs (
        :warn,
        r"solve\(\): stopped at t = .*pore_closure: cell \d+: raw porosity",
    ) match_mode = :any Pyrolysis.Integration.solve(_closing(); saveat = [0.0, 20.0])
    @test sol.retcode === :DtLessThanMin
    @test sol.outcome === :pore_closure
    @test !succeeded(sol)
    @test 0.0 < sol.t[end] < 20.0
end
