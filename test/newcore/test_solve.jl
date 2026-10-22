using Test
using LinearAlgebra
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter
using Pyrolysis.States
using Pyrolysis.Meshes: n_cells
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration
using OrdinaryDiffEq: ODEFunction, ODEProblem
import OrdinaryDiffEq
using OrdinaryDiffEqSDIRK: KenCarp4

# Two closed cells exchanging heat by conduction, with a temperature-dependent
# heat capacity: summed enthalpy is a linear invariant.
function _two_cells()
    c = LinearProperty(400.0, 2.0)
    solid = SolidComponent(:s, ρ = 1000.0, c = c, k = 50.0)
    mat = Material(name = :pair, components = (solid,), reactions = ())
    mesh = create_uniform_mesh(0.02, 2)
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> z < 0.01 ? 1000.0 : 300.0,
        ξ_initial = [z -> 1000.0],
        tspan = (0.0, 40.0),
    )
end

function _enthalpy(cp, u)
    layout, s = cp.def.layout, cp.def.setup
    ds = cp.ws.derived
    derive!(ds, u, layout, s)
    return sum(
        ds.m[1, i] * sensible_enthalpy(s.laws.h[1], ds.T[i]) for i = 1:n_cells(layout)
    )
end

function _fixed_step(cp, dt)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    f = ODEFunction(
        (du, u, p, t) -> residual!(du, u, def, t, ws);
        jac = (J, u, p, t) -> jacobian!(J, cj, u, t),
        jac_prototype = jac_prototype(cj),
    )
    prob = ODEProblem(f, cp.u0, (0.0, 40.0))
    return OrdinaryDiffEq.solve(
        prob,
        KenCarp4(autodiff = OrdinaryDiffEq.AutoFiniteDiff());
        adaptive = false,
        dt = dt,
        reltol = 1e-12,
        abstol = 1e-14,
    )
end

@testset "energy forms: invariant to round-off, temperature form at fourth order" begin
    errors = Float64[]
    for dt in (0.5, 0.25, 0.125)
        cpH = Integration.setup(_two_cells(); energy_state = :enthalpy)
        H0 = _enthalpy(cpH, cpH.u0)
        solH = _fixed_step(cpH, dt)
        @test solH.retcode == OrdinaryDiffEq.ReturnCode.Success
        @test maximum(abs(_enthalpy(cpH, u) - H0) for u in solH.u) <= 1e-12 * abs(H0)
        cpT = Integration.setup(_two_cells(); energy_state = :temperature)
        solT = _fixed_step(cpT, dt)
        @test solT.retcode == OrdinaryDiffEq.ReturnCode.Success
        push!(errors, maximum(abs(_enthalpy(cpT, u) - H0) for u in solT.u) / abs(H0))
    end
    @test errors[1] > 1e-10
    orders = log2.(errors[1:2] ./ errors[2:3])
    @test all(p -> 3.5 < p < 4.6, orders)
end

# One reacting cell between Dirichlet faces, as the kinetics cases use it.
function _reacting_cell()
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    b = SolidComponent(:b, ρ = 1000.0, c = 1000.0, k = 0.5)
    mat = Material(
        name = :cell,
        components = (a, b),
        reactions = (Reaction(:r, 1 => 2; A = 0.5, E = 0.0, h = 0.0),),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    mesh = create_uniform_mesh(1e-3, 1)
    bcs = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(400.0), bottom = DirichletThermalBC(400.0)),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bcs,
        T_initial = z -> 400.0,
        ξ_initial = [z -> 1000.0, z -> 0.0],
        tspan = (0.0, 4.0),
        output_times = [0.0, 1.0, 2.0, 4.0],
    )
end

@testset "solve: kinetics through the adapter, both energy forms" begin
    for energy in (:temperature, :enthalpy)
        problem = _reacting_cell()
        positions = copy(problem.mesh.positions)
        sol = Integration.solve(
            problem;
            energy_state = energy,
            abstol = 1e-12,
            reltol = 1e-10,
        )
        @test sol.retcode === :Success
        @test sol.t == [0.0, 1.0, 2.0, 4.0]
        @test maximum(abs(sol.ξ[1, 1, k] - 1000.0 * exp(-0.5 * sol.t[k])) for k = 1:4) <
              1e-6
        @test maximum(abs(sum(sol.ξ[:, 1, k]) - 1000.0) for k = 1:4) < 1e-10
        @test maximum(abs.(sol.T .- 400.0)) < 1e-9
        @test sol.statistics.accepted > 0
        # The caller's problem is not touched.
        @test problem.mesh.positions == positions
    end
    @test_throws ArgumentError Integration.solve(_reacting_cell(); bogus = 1)
    @test_throws ArgumentError Integration.solve(
        _reacting_cell();
        use_ale = true,
        w_override = (x, t) -> 0.0,
    )
    @test_throws ArgumentError Integration.solve(
        _reacting_cell();
        use_ale = true,
        min_thickness = 0.0,
    )
    @test_throws ArgumentError Integration.solve(_reacting_cell(); energy_state = :entropy)
end

@testset "outputs: face temperatures from the boundary balances" begin
    solid = SolidComponent(:s, ρ = 800.0, c = 1200.0, k = 0.4)
    mat = Material(name = :heated, components = (solid,), reactions = ())
    L, n, q = 0.01, 10, 2e4
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(L, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(q), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 800.0],
        tspan = (0.0, 20.0),
        output_times = [0.0, 10.0, 20.0],
    )
    sol = Integration.solve(problem; abstol = 1e-10, reltol = 1e-8)
    @test sol.retcode === :Success
    # Prescribed flux: T_s = T_cell + q d / k with d the half cell; adiabatic:
    # the face temperature is the cell temperature.
    @test sol.surface_T ≈ sol.T[end, :] .+ q * (L / n / 2) / 0.4 rtol = 1e-12
    @test sol.boundary_T[:bottom] == sol.T[1, :]
    @test sol.surface_T == sol.boundary_T[:top]
    @test get_cell_positions(sol) ≈ collect(range(L / 2n, L - L / 2n, length = n))
    # The prescribed flux is what the surface applies and conducts.
    @test sol.extras.surface_heat_flux_net == fill(q, 3)
    @test sol.extras.surface_heat_flux_conduction ≈ fill(q, 3) rtol = 1e-12
    @test sol.extras.back_temperature == sol.boundary_T[:bottom]
    @test sol.extras.back_heat_flux == zeros(3)
    @test sol.mass_loss_rate == zeros(3)
    @test sol.extras.live_cells == fill(n, 3)
    @test isempty(sol.extras.park_times)
end

@testset "outputs: the surface-flux decomposition and the mass-loss rate" begin
    solid = SolidComponent(:s, ρ = 800.0, c = 1200.0, k = 0.4, ε = 0.9, α = 2.0, λ = 1e-6)
    gas = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(:gasify, 1 => 2; A = 1e8, E = 1.2e5, h = 1e5)
    mat = Material(name = :charring, components = (solid, gas), reactions = (rxn,))
    h_back, T_back_inf = 5.0, 300.0
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 12),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 4e4) +
                      RadiativeBC(T_inf = 300.0) +
                      ConvectiveBC(h = 10.0, T_inf = 300.0),
                bottom = ConvectiveBC(h = h_back, T_inf = T_back_inf),
            ),
            mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 500.0,
        ξ_initial = [z -> 800.0, z -> 0.0],
        tspan = (0.0, 60.0),
        output_times = collect(0.0:10.0:60.0),
    )
    for radiation_model in (SURFACE_ABSORPTION, BEER_LAMBERT)
        sol = Integration.solve(problem; radiation_model, abstol = 1e-8, reltol = 1e-6)
        @test sol.retcode === :Success
        x = sol.extras
        # The decomposition sums to the applied flux, which the surface balance
        # conducts into the cell.
        @test x.surface_heat_flux_net ≈
              x.surface_heat_flux_absorbed .- x.surface_heat_flux_reradiation .-
              x.surface_heat_flux_convection rtol = 1e-12
        @test x.surface_heat_flux_net ≈ x.surface_heat_flux_conduction rtol = 1e-8
        @test all(>(0), x.surface_heat_flux_reradiation[2:end])
        @test all(>(0), x.surface_heat_flux_convection[2:end])
        if radiation_model === SURFACE_ABSORPTION
            @test x.surface_heat_flux_absorbed ≈ 4e4 .* x.surface_emissivity rtol = 1e-12
            @test all(iszero, x.volumetric_absorption)
        else
            # The incident flux enters in depth, not through the face.
            @test all(iszero, x.surface_heat_flux_absorbed)
            @test all(>(0), x.volumetric_absorption)
        end
        @test x.back_temperature == sol.boundary_T[:bottom]
        @test x.back_heat_flux ≈ h_back .* (T_back_inf .- x.back_temperature) rtol = 1e-12
        # Gas leaves through the top face only.
        @test sol.mass_loss_rate ≈
              sol.boundary_mass_flux[:top] .+ sol.boundary_mass_flux[:bottom] rtol = 1e-12
        @test x.surface_mass_flux == sol.boundary_mass_flux[:top]
        @test sol.mass_loss_rate[end] > 0
        @test x.MLR_gauge == sol.mass_loss_rate
        @test all(0 .< x.surface_emissivity .<= 1)
        @test all(>(0), x.surface_conductivity)
        @test x.live_cells == fill(12, length(sol.t))
        @test isempty(x.park_times)
    end
end
