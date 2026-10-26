using Test
using Pyrolysis

# A solid→gas material, with a lateral shrinkage law on request, on a uniform
# 1 cm column under a 10 kW/m² flux with lumped wall venting at the top.
function _cross_section_problem(; lateral::Bool, tspan = (0.0, 0.02))
    solid = SolidComponent(
        :solid;
        φ = 0.25,
        λ = 0.0,
        ρ = 1000.0,
        c = 1000.0,
        k = 0.2,
        ε = 0.9,
        γ = 1.0,
    )
    gas = GaseousComponent(:gas; M = 0.050, c = 1000.0, k = 0.02, γ = 1.0)
    material = Material(
        name = :cross_section_test,
        components = (solid, gas),
        reactions = (Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5),),
        lateral_shrinkage_law = lateral ? (χ -> 1.0 - 0.2 * χ) : nothing,
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 5),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1000.0, z -> 0.0],
        tspan = tspan,
    )
end

@testset "Variable Cross Section" begin
    @testset "a lateral law adds the history pair, the area and the lateral rate" begin
        problem = _cross_section_problem(lateral = true)
        cp = Pyrolysis.Integration.setup(problem; use_ale = true)
        layout = cp.def.layout
        n = Pyrolysis.States.n_cells(layout)
        @test n == 5
        @test Pyrolysis.States.has_chi(layout)
        # Per cell T, two masses and the history pair; per node position and
        # velocity; the area and the lateral rate.
        @test length(cp.u0) == 5 * (1 + 2 + 2) + 6 * 2 + 2
        area = Pyrolysis.States.area_index(layout)
        A0 = cp.u0[area]

        # Each cell has turned half its initial dry solid into dry gas: χ̄ = 0.5,
        # so the law gives A = (1 − 0.2·0.5) A₀ = 0.9 A₀.
        u = copy(cp.u0)
        for i = 1:n
            M = u[Pyrolysis.States.history_index(layout, i, 2)]
            @test M > 0
            u[Pyrolysis.States.history_index(layout, i, 1)] = 0.5 * M
        end
        du = similar(u)
        Pyrolysis.Assembly.residual!(du, u, cp.def, 0.0, cp.ws)
        @test abs(du[area]) > 0.05 * A0
        u[area] = 0.9 * A0
        Pyrolysis.Assembly.residual!(du, u, cp.def, 0.0, cp.ws)
        @test abs(du[area]) <= 1e-14 * A0

        # The residual is finite, the cross-section rows included.
        @test all(isfinite, du)
        @test isfinite(du[Pyrolysis.States.lateral_rate_index(layout)])
    end

    @testset "without a lateral law the moving layout carries no cross-section rows" begin
        cp = Pyrolysis.Integration.setup(
            _cross_section_problem(lateral = false);
            use_ale = true,
        )
        @test !Pyrolysis.States.has_chi(cp.def.layout)
        @test length(cp.u0) == 5 * (1 + 2) + 6 * 2
    end

    # Only a moving mesh evolves the cross section, so a lateral shrinkage
    # law on the fixed mesh is refused wherever the problem is set up.
    @testset "a lateral law needs use_ale = true" begin
        problem = _cross_section_problem(lateral = true)
        message(f) =
            try
                f()
                ""
            catch err
                err isa ArgumentError ? err.msg : "unexpected $(typeof(err))"
            end
        for f in (
            () -> Pyrolysis.Integration.setup(problem; use_ale = false),
            () -> Pyrolysis.Integration.solve(problem),
            () -> Pyrolysis.Integration.solve(problem; validate = false),
        )
            msg = message(f)
            @test occursin("lateral_shrinkage_law", msg)
            @test occursin("use_ale = true", msg)
        end
        @test Pyrolysis.States.has_chi(
            Pyrolysis.Integration.setup(problem; use_ale = true).def.layout,
        )
    end

    @testset "a variable-cross-section solve records a finite area history" begin
        problem = _cross_section_problem(lateral = true)
        sol = Pyrolysis.Integration.solve(problem; use_ale = true, validate = false)
        @test sol.retcode == :Success
        @test all(isfinite, sol.extras.cross_section_area_history)
    end
end
