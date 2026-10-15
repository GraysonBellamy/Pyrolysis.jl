"""
The progress bar observes the integration without steering it: a solve with
`progress = true` takes exactly the steps of the same solve with
`progress = false`, on the fixed and on the ALE mesh.
"""

using Test
using Pyrolysis

function _progress_problem()
    mat = Material(
        name = :progress_fixture,
        components = (
            SolidComponent(:s, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86, λ = 1e-6),
            GaseousComponent(:g, M = 0.1, c = 1500.0, k = 0.02, λ = 1e-5),
        ),
        reactions = (Reaction(:r, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3),),
    )
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) + ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.005, 20, Val(2)),
        material = mat,
        bc_set = bc,
        tspan = (0.0, 60.0),
        output_times = collect(0.0:5.0:60.0),
    )
end

function _progress_fingerprint(sol)
    st = sol.raw_solution.stats
    return (sol.retcode, st.naccept, st.nreject, st.nf, sol.t, sol.T, sol.ξ, sol.z)
end

@testset "progress bar leaves the steps unchanged" begin
    for use_ale in (false, true)
        @testset "use_ale = $use_ale" begin
            kw = (use_ale = use_ale, radiation_model = SURFACE_ABSORPTION)
            quiet = solve(_progress_problem(); kw..., progress = false)
            shown = redirect_stderr(devnull) do
                solve(_progress_problem(); kw..., progress = true)
            end
            @test quiet.retcode == :Success
            @test isequal(_progress_fingerprint(shown), _progress_fingerprint(quiet))
        end
    end
end
