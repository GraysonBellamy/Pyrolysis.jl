"""
A solve that stops because the domain test keeps rejecting a cell under gas
transport says so: the outcome is recorded and `solve` warns with the cell and
the quantity that left the domain, even when the step reaches `dtmin` before
the test gives up.

The deck is a slab whose solid converts entirely to gas, under lumped transport
without an endpoint: its top cell burns through the matrix, which ends the
supported domain (`:matrix_exhausted`). A residue keeps the same deck inside it.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: R_GAS

function _charfree(; char = 0.0)
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
    vol = GaseousComponent(:volatiles; M = 0.070, c = 2000.0, k = 0.03)
    n2 = GaseousComponent(:N2; M = 0.028014, c = 1100.0, k = 0.04)
    pore_n2 = z -> 0.3 * 101325.0 * 0.028014 / (R_GAS * 360.0)
    components, reaction, ξ = if char > 0
        residue = SolidComponent(
            :char;
            φ = 0.3,
            λ = 1e-5,
            κ = 1e-12,
            ρ = 200.0,
            c = 1000.0,
            k = 0.1,
        )
        (
            (solid, residue, vol, n2),
            Reaction(
                :decomp,
                1 => (2, 3);
                A = 1e5,
                E = 5e4,
                h = 1e5,
                yields = (char, 1 - char),
            ),
            [z -> 1000.0, z -> 0.0, z -> 0.0, pore_n2],
        )
    else
        (
            (solid, vol, n2),
            Reaction(:decomp, 1 => 2; A = 1e5, E = 5e4, h = 1e5),
            [z -> 1000.0, z -> 0.0, pore_n2],
        )
    end
    material = Material(name = :charfree, components = components, reactions = (reaction,))
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
        ξ_initial = ξ,
        tspan = (0.0, 20.0),
    )
end

@testset "a domain outcome names its cell" begin
    sol = @test_logs (
        :warn,
        r"solve\(\): stopped at t = .*matrix_exhausted: cell \d+: matrix_exhausted at matrix occupancy",
    ) match_mode = :any Pyrolysis.Integration.solve(_charfree(); saveat = [0.0, 20.0])
    @test sol.retcode === :DtLessThanMin
    @test sol.outcome === :matrix_exhausted
    @test !succeeded(sol)
    @test 10.0 < sol.t[end] < 20.0

    # With a residue the matrix stays and the same deck completes, without the
    # warning.
    sol_char = @test_logs min_level = Base.CoreLogging.Warn Pyrolysis.Integration.solve(
        _charfree(char = 0.2);
        saveat = [0.0, 20.0],
    )
    @test sol_char.retcode === :Success
    @test sol_char.outcome === :none
end
