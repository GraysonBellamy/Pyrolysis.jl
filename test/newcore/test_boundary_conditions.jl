using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT, SURFACE_ABSORPTION, NO_RADIATION
using Pyrolysis.States
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.Integration
using Pyrolysis.Integration: setup

# A solid with an emissivity and an absorption coefficient, and a gas with a
# gas-transfer coefficient.
function _bc_material(; α::Real = 0.0)
    solid = SolidComponent(:solid, ρ = 600.0, c = 1500.0, k = 0.15, ε = 0.9, α = α)
    gas = GaseousComponent(:gas, M = 0.05, c = 1200.0, k = 0.03, λ = 1e-5)
    return Material(name = :bc, components = (solid, gas), reactions = ())
end

_bc_problem(bc; n = 6, tspan = (0.0, 1.0), α = 0.0, kwargs...) = PyrolysisProblem(;
    mesh = create_uniform_mesh(0.01, n, Val(2)),
    material = _bc_material(; α),
    bc_set = bc,
    tspan = tspan,
    kwargs...,
)

@testset "the residual carries no surface-balance state between evaluations" begin
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) +
                  RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    for energy in (:temperature, :enthalpy)
        cp = setup(_bc_problem(bc); energy_state = energy)
        fresh = setup(_bc_problem(bc); energy_state = energy)
        layout = cp.def.layout
        u = copy(cp.u0)
        for i = 1:6
            u[mass_index(layout, i, 2)] = 1e-4 * i
        end
        du_fresh = zeros(length(u))
        residual!(du_fresh, u, fresh.def, 0.0, fresh.ws)
        @test all(isfinite, du_fresh)
        du = zeros(length(u))
        # A NaN inventory is refused by the inversion (`invert_caloric`), so
        # the inventory form is driven with a wild finite one only.
        wilds = energy === :temperature ? (-9.1e40, NaN) : (-9.1e40,)
        for wild in wilds
            u_wild = copy(u)
            u_wild[energy_index(layout, 6)] = wild
            residual!(similar(u), u_wild, cp.def, 0.0, cp.ws)
            residual!(du, u, cp.def, 0.0, cp.ws)
            @test all(isfinite, du)
            @test isapprox(du, du_fresh; rtol = 1e-12)
        end
    end
end

@testset "an ambient gas enters a gas-free column through a film condition" begin
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (
            top = ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.5)),
            bottom = ImpermeableBC(),
        ),
    )
    problem = _bc_problem(bc; T_initial = z -> 400.0, ξ_initial = [z -> 400.0, z -> 0.0])
    cp = setup(problem)
    layout = cp.def.layout
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    @test du[mass_index(layout, 6, 2)] > 0
    @test all(i -> du[mass_index(layout, i, 1)] == 0, 1:6)
    ds, c = buffers(cp.ws, cp.def, Float64)
    @test all(iszero, c.F_mass[1, :])
end

@testset "mass leaving through both faces closes against the inventories" begin
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.05), bottom = DiffusiveMassBC(0.05)),
    )
    problem = _bc_problem(
        bc;
        tspan = (0.0, 40.0),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 400.0, z -> 1.0],
        output_times = collect(0.0:0.01:40.0),
    )
    sol = Integration.solve(problem; reltol = 1e-8, abstol = 1e-10)
    @test sol.retcode === :Success
    m = vec(sum(sol.m; dims = (1, 2)))
    lost = m[1] - m[end]
    @test lost > 0
    # Both faces vent: the bottom flux is outward too and, in the symmetric
    # start, as large as the top one.
    out = sol.boundary_mass_flux[:top] .+ sol.boundary_mass_flux[:bottom]
    @test sol.boundary_mass_flux[:bottom][2] ≈ sol.boundary_mass_flux[:top][2] rtol = 1e-6
    t = sol.t
    ∫out = sum(0.5 * (out[i] + out[i+1]) * (t[i+1] - t[i]) for i = 1:(length(t)-1))
    @test ∫out ≈ lost rtol = 1e-3
end

@testset "a prescribed face temperature is the reported surface temperature" begin
    bc = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(600.0), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
    )
    sol = Integration.solve(
        _bc_problem(bc; tspan = (0.0, 5.0), output_times = [0.0, 2.5, 5.0]),
    )
    @test sol.retcode === :Success
    @test all(==(600.0), sol.surface_T)
end

@testset "in-depth radiation enters through the top faces only" begin
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = RadiativeFluxBC(flux = 10e3)),
    )
    problem = _bc_problem(bc; α = 200.0)
    @test_throws ArgumentError setup(problem; radiation_model = BEER_LAMBERT)
    sol = Integration.solve(problem; radiation_model = SURFACE_ABSORPTION)
    @test sol.retcode === :Success
    @test_throws ArgumentError setup(
        problem;
        radiation_model = Pyrolysis.Materials.P1_QUASI_STEADY,
    )
end

@testset "surface absorption is the surface balance of the incident flux" begin
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 30e3) + ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
    )
    a = setup(_bc_problem(bc); radiation_model = NO_RADIATION)
    b = setup(_bc_problem(bc); radiation_model = SURFACE_ABSORPTION)
    du_a, du_b = zeros(length(a.u0)), zeros(length(b.u0))
    residual!(du_a, a.u0, a.def, 0.0, a.ws)
    residual!(du_b, b.u0, b.def, 0.0, b.ws)
    @test du_a == du_b
    # Under BEER_LAMBERT the incident flux leaves the surface balance and is
    # absorbed in depth instead: the same total energy rate when nothing is
    # transmitted and nothing else depends on the surface temperature.
    alone = BoundaryConditionSet(
        thermal = (top = RadiativeFluxBC(flux = 30e3), bottom = AdiabaticBC()),
    )
    c = setup(_bc_problem(alone; α = 1e3); radiation_model = BEER_LAMBERT)
    s = setup(_bc_problem(alone; α = 1e3); radiation_model = SURFACE_ABSORPTION)
    rate(cp) = begin
        du = zeros(length(cp.u0))
        residual!(du, cp.u0, cp.def, 0.0, cp.ws)
        ds, _ = buffers(cp.ws, cp.def, Float64)
        sum(ds.cp_sum[i] * du[energy_index(cp.def.layout, i)] for i = 1:6)
    end
    @test rate(c) ≈ rate(s) rtol = 1e-10
end

@testset "a time-varying incident flux is sampled at every evaluation" begin
    bc_of(flux) = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = flux) + ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
    )
    q0 = 30e3
    prob(flux) = PyrolysisProblem(
        mesh = create_uniform_mesh(0.005, 8, Val(2)),
        material = _bc_material(α = 500.0),
        bc_set = bc_of(flux),
        tspan = (0.0, 1.0),
        output_times = [0.0, 0.5, 1.0],
    )
    # sin(2π·0) = sin(2π·1) = 0: sampling only the window's ends would freeze
    # the deposit at q₀ and make the two runs identical.
    var =
        Integration.solve(prob(t -> q0 * (1 + sin(2π * t))); radiation_model = BEER_LAMBERT)
    con = Integration.solve(prob(q0); radiation_model = BEER_LAMBERT)
    @test var.retcode === :Success && con.retcode === :Success
    @test var.surface_T[2] > con.surface_T[2] + 1.0
    @test all(isfinite, var.extras.volumetric_absorption)
    @test var.extras.volumetric_absorption[2] > con.extras.volumetric_absorption[2]
end

@testset "blowing lowers the convective loss of an outgassing surface" begin
    function run(model)
        solid = SolidComponent(:solid, ρ = 1000.0, c = 1500.0, k = 0.2, λ = 1e-5)
        gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 1e-5)
        mat = Material(
            name = :outgas,
            components = (solid, gas),
            reactions = (Reaction(:g, 1 => 2; A = 1e6, E = 1e5, h = 0.0),),
        )
        bc = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(50e3) +
                      ConvectiveBC(h = 15.0, T_inf = 300.0, blowing = model),
                bottom = AdiabaticBC(),
            ),
            mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
        )
        problem = PyrolysisProblem(
            mesh = create_uniform_mesh(0.005, 8, Val(2)),
            material = mat,
            bc_set = bc,
            T_initial = z -> 600.0,
            ξ_initial = [z -> 1000.0, z -> 0.0],
            tspan = (0.0, 5.0),
        )
        return Integration.solve(problem; saveat = [0.0, 2.5, 5.0])
    end
    off = run(NoBlowing())
    on = run(CouetteFilm())
    @test off.retcode === :Success && on.retcode === :Success
    @test on.surface_T[end] > off.surface_T[end]
    @test on.extras.surface_mass_flux[end] > off.extras.surface_mass_flux[end]
    @test all(==(1.0), off.extras.surface_blowing_ratio)
    @test on.extras.surface_blowing_ratio[end] < 1.0
end
