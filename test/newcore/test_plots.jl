using Test
using Pyrolysis
using Pyrolysis.Integration
using RecipesBase

# The plotting backend declares the attributes it supports; without one, every
# attribute is.
RecipesBase.is_key_supported(::Symbol) = true

@testset "plots recipes on the finite-volume core's solution" begin
    @test Base.get_extension(Pyrolysis, :PyrolysisPlotsExt) !== nothing
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 6),
        material = Material(;
            name = :inert,
            components = (SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2),),
            reactions = (),
        ),
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(1e4), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1000.0],
        tspan = (0.0, 5.0),
        output_times = collect(0.0:1.0:5.0),
    )
    sol = Integration.solve(problem)
    recipe(attributes) = only(RecipesBase.apply_recipe(attributes, sol))
    for (kind, series) in ((:surface_T, sol.surface_T), (:mlr, sol.mass_loss_rate))
        attributes = Dict{Symbol,Any}(:kind => kind)
        @test recipe(attributes).args == (sol.t, series)
        @test attributes[:seriestype] === :line
    end
    @test recipe(Dict{Symbol,Any}()).args == (sol.t, sol.surface_T)
    attributes = Dict{Symbol,Any}(:kind => :heatmap)
    @test recipe(attributes).args == (sol.t, axes(sol.T, 1), sol.T)
    @test attributes[:seriestype] === :heatmap
    @test_throws ArgumentError recipe(Dict{Symbol,Any}(:kind => :profile))
end
