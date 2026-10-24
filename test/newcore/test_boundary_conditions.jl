using Test
using Logging
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT, SURFACE_ABSORPTION, NO_RADIATION
using Pyrolysis.Materials: STEFAN_BOLTZMANN, T_PROPERTY_MAX
using Pyrolysis.BoundaryConditions: BoundaryState, solve_boundary_temperature!
using Pyrolysis.Meshes: boundary_faces
using Pyrolysis.States
using Pyrolysis.Kernels: ThermalBoundaryKernel, find_kernel
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
    mesh = create_uniform_mesh(0.01, n),
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
    # On a fixed and on a moving mesh.
    for energy in (:temperature, :enthalpy), use_ale in (false, true)
        cp = setup(_bc_problem(bc); energy_state = energy, use_ale)
        fresh = setup(_bc_problem(bc); energy_state = energy, use_ale)
        layout = cp.def.layout
        u = copy(cp.u0)
        for i = 1:6
            u[mass_index(layout, i, 2)] = 1e-4 * i
        end
        du_fresh = zeros(length(u))
        residual!(du_fresh, u, fresh.def, 0.0, fresh.ws)
        @test all(isfinite, du_fresh)
        du = zeros(length(u))
        for wild in (-9.1e40, NaN)
            u_wild = copy(u)
            u_wild[energy_index(layout, 6)] = wild
            du_wild = similar(u)
            residual!(du_wild, u_wild, cp.def, 0.0, cp.ws)
            isnan(wild) && @test any(isnan, du_wild)
            residual!(du, u, cp.def, 0.0, cp.ws)
            @test all(isfinite, du)
            @test isapprox(du, du_fresh; rtol = 1e-12)
        end
    end
end

@testset "the surface balance reads the trial clamp: no clamp, no warning, finite derivatives" begin
    material = _bc_material()
    k, d = 0.15, 1e-3
    # The face temperature from the cell temperature `T_cell`, with or
    # without the trial clamp.
    function face(bc, T_cell; trial, J = (0.0, 0.0))
        state = BoundaryState(300.0, 0.0, 0.0, 300.0, k, d, false, 0)
        return solve_boundary_temperature!(
            state,
            T_cell,
            k,
            d,
            bc,
            material,
            J,
            0.0,
            0.9;
            trial,
        )
    end
    # An incident flux the hottest admissible surface cannot reradiate: the
    # balance has no root in the property bracket.
    hot = RadiativeFluxBC(flux = 1e7) + RadiativeBC(T_inf = 300.0)
    q_max = 0.9 * 1e7 - 0.9 * STEFAN_BOLTZMANN * (T_PROPERTY_MAX^4 - 300.0^4)
    r = @test_logs min_level = Logging.Warn face(hot, 2500.0; trial = trial_temperature)
    @test r.T_s ≈ 2500 + q_max * d / k rtol = 1e-10
    @test r.q_conducted ≈ q_max rtol = 1e-10
    # Beyond the bracket the condition is constant, so the face follows the
    # cell one for one.
    @test ForwardDiff.derivative(T -> face(hot, T; trial = trial_temperature).T_s, 2500.0) ≈
          1 rtol = 1e-12
    # The hard clamp pins the face at the bracket's end with zero derivative,
    # and warns.
    @test_logs (:warn, r"pinned at the physical clamp") match_mode = :any begin
        @test face(hot, 2500.0; trial = nothing).T_s == T_PROPERTY_MAX
        @test ForwardDiff.derivative(T -> face(hot, T; trial = nothing).T_s, 2500.0) == 0
    end
    # A NaN cell temperature: a NaN face, no warning.
    r = @test_logs min_level = Logging.Warn face(hot, NaN; trial = trial_temperature)
    @test isnan(r.T_s)
    # On a physical state the two agree bit for bit, derivatives included.
    mild =
        RadiativeFluxBC(flux = 5e4) +
        RadiativeBC(T_inf = 300.0) +
        ConvectiveBC(h = 10.0, T_inf = 300.0)
    for T_cell in (350.0, 600.0, 900.0), J in ((0.0, 0.0), (0.0, 2e-3))
        a = face(mild, T_cell; trial = trial_temperature, J)
        b = face(mild, T_cell; trial = nothing, J)
        @test (a.T_s, a.q_applied, a.q_conducted) === (b.T_s, b.q_applied, b.q_conducted)
        @test ForwardDiff.derivative(
            T -> face(mild, T; trial = trial_temperature, J).T_s,
            T_cell,
        ) === ForwardDiff.derivative(T -> face(mild, T; trial = nothing, J).T_s, T_cell)
    end
    # Through the residual: a trial state far outside the bracket neither
    # clamps nor warns.
    bc = BoundaryConditionSet(
        thermal = (top = hot, bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    cp = setup(_bc_problem(bc))
    for wild in (1e5, -1e5, NaN)
        u = copy(cp.u0)
        u[energy_index(cp.def.layout, 6)] = wild
        du = similar(u)
        @test_logs min_level = Logging.Warn residual!(du, u, cp.def, 0.0, cp.ws)
        @test isnan(wild) || all(isfinite, du)
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
    # The film law's value along the outward normal, −h_m (ξ_inf − ξ): the
    # ambient gas enters at h_m ξ_inf, also at a slightly negative
    # concentration, and the surface cell gains exactly that.
    top = only(boundary_faces(cp.def.setup.topology, :top))
    A = ds.geometry.face_area[top]
    @test c.F_mass[2, top] ≈ -0.05 * 0.5 rtol = 1e-12
    @test du[mass_index(layout, 6, 2)] ≈ A * 0.05 * 0.5 rtol = 1e-12
    u = copy(cp.u0)
    for i = 1:6
        u[mass_index(layout, i, 2)] = -1e-12 * ds.V[i]
    end
    residual!(du, u, cp.def, 0.0, cp.ws)
    @test c.F_mass[2, top] ≈ -0.05 * (0.5 + 1e-12) rtol = 1e-12
    # A uniform gas above the ambient: no interior flux, and the surface cell
    # loses h_m (ξ − ξ_inf).
    for i = 1:6
        u[mass_index(layout, i, 2)] = 1.0 * ds.V[i]
    end
    residual!(du, u, cp.def, 0.0, cp.ws)
    @test c.F_mass[2, top] ≈ 0.05 * (1.0 - 0.5) rtol = 1e-12
    @test du[mass_index(layout, 6, 2)] ≈ -A * 0.05 * 0.5 rtol = 1e-12
    @test all(i -> du[mass_index(layout, i, 2)] == 0, 1:5)
    # The solution reports the inflow of the gas-free start, on a fixed and on
    # a moving mesh.
    for use_ale in (false, true)
        sol = Integration.solve(problem; use_ale)
        @test sol.retcode === :Success
        @test sol.boundary_mass_flux[:top][1] ≈ -0.05 * 0.5 rtol = 1e-12
    end
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
    # The mass-loss rate counts both faces (unit cross section).
    @test sol.mass_loss_rate ≈ out rtol = 1e-12
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
    # The face's net flux is the conducted one: finite, and into the column
    # from its 300 K start.
    @test all(isfinite, sol.extras.surface_heat_flux_net)
    @test sol.extras.surface_heat_flux_net[2] > 0
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
        mesh = create_uniform_mesh(0.005, 8),
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
            mesh = create_uniform_mesh(0.005, 8),
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

# A heated convective face over a gas-free column, with a prescribed gas flux
# through the same face out of the column, into it, or none.
@testset "blowing reads the outward gas flux of its face, at :top and at :bottom" begin
    exposed =
        HeatFluxBC(3e4) + ConvectiveBC(h = 20.0, T_inf = 300.0, blowing = CouetteFilm())
    conditions =
        (outflow = MassFluxBC(-0.01), none = ImpermeableBC(), inflow = MassFluxBC(0.01))
    for (tag, other) in ((:top, :bottom), (:bottom, :top))
        face = map(conditions) do mass
            bc = BoundaryConditionSet(
                thermal = NamedTuple{(tag, other)}((exposed, AdiabaticBC())),
                mass = NamedTuple{(tag, other)}((mass, ImpermeableBC())),
            )
            problem = _bc_problem(
                bc;
                T_initial = z -> 700.0,
                ξ_initial = [z -> 600.0, z -> 0.0],
            )
            cp = setup(problem)
            residual!(zeros(length(cp.u0)), cp.u0, cp.def, 0.0, cp.ws)
            ds, c = buffers(cp.ws, cp.def, Float64)
            f = only(boundary_faces(cp.def.setup.topology, tag))
            thermal = find_kernel(ThermalBoundaryKernel, cp.def.kernels)
            (J = c.F_mass[2, f], χ = thermal.states[f].blowing_ratio, T_s = c.T_face[f])
        end
        # The flux leaves along the face's outward normal at either face.
        @test face.outflow.J ≈ 0.01 rtol = 1e-14
        @test face.none.J == 0
        @test face.inflow.J ≈ -0.01 rtol = 1e-14
        # Outflow suppresses the convective coefficient, inflow enhances it, so
        # the heated face runs hotter under outflow and cooler under inflow.
        @test face.none.χ == 1
        @test face.outflow.χ < 1 < face.inflow.χ
        @test face.outflow.T_s > face.none.T_s > face.inflow.T_s
    end
end
