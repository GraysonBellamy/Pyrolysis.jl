"""
Repeated solves of one problem in one process are bit-identical, for each
default sparse factorization. The integrator allocates its stage matrix
without initializing it, so a linear solver that chose its factorization
pattern from the starting matrix would pick up whatever the memory held:
round-off, and with it the step sequence, would change from run to run. The
first solve in a process usually receives zero-filled memory and later ones
reused memory, so three solves are compared.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT

function _determinism_problem(n_cells::Int)
    mat = Material(
        name = :determinism_fixture,
        components = (
            SolidComponent(:s, ρ = 600.0, c = 1500.0, k = 0.15, α = 50.0, λ = 1e-6),
            GaseousComponent(:g, M = 0.05, c = 1200.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (Reaction(:r, 1 => 2, A = 1e10, E = 150e3, h = 300e3),),
    )
    mesh = create_uniform_mesh(0.005, n_cells, Val(2))
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) + ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc,
        tspan = (0.0, 60.0),
        output_times = collect(0.0:5.0:60.0),
    )
end

function _fingerprint(sol)
    st = sol.raw_solution.stats
    return (sol.retcode, st.naccept, st.nf, sol.t, sol.T, sol.ξ, sol.surface_T)
end

@testset "repeated solves are bit-identical" begin
    for (label, kw) in (
        ("fixed mesh, Beer–Lambert (KLU)", (use_ale = false,)),
        ("ALE, DirectSolve (UMFPACK)", (use_ale = true,)),
        (
            "ALE with depletion, ApproxSolve (KLU)",
            (use_ale = true, handle_depletion = true),
        ),
    )
        @testset "$label" begin
            fps = [
                _fingerprint(
                    solve(
                        _determinism_problem(30);
                        radiation_model = BEER_LAMBERT,
                        progress = false,
                        kw...,
                    ),
                ) for _ = 1:3
            ]
            @test fps[1][1] == :Success
            @test isequal(fps[1], fps[2])
            @test isequal(fps[1], fps[3])
        end
    end
end
