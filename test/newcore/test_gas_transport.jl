using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter
using Pyrolysis.Discretization: boundary_darcy_velocity
using Pyrolysis.Properties: harmonic_mean_weighted
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: GasTransportKernel, outward_darcy_velocity, MOBILITY_GUARD
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration: setup

# A decomposing slab with a product gas and a background gas: temperature-
# dependent properties, gas-transfer coefficients on every component and, for
# Darcy–Fick transport, an impermeable virgin solid and a permeable char, so
# the char-free cells near the back carry zero mobility and the face
# porosities differ across the char front.
function _gas_slab(n; darcy::Bool, thermal, mass)
    virgin = SolidComponent(
        :virgin,
        ρ = 600.0,
        c = (1200.0, 2.0),
        k = (0.15, 1e-4),
        λ = 2e-6,
        ρ_intrinsic = 1400.0,
    )
    char = SolidComponent(
        :char,
        ρ = 150.0,
        c = 900.0,
        k = 0.08,
        λ = 1e-5,
        κ = darcy ? 1e-11 : 0.0,
        ρ_intrinsic = 1800.0,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = (1000.0, 0.4), k = 0.03, λ = 3e-5)
    air = GaseousComponent(:air, M = 0.029, c = 1006.0, k = 0.026, λ = 2e-5)
    mat = Material(
        name = :gas_slab,
        components = (virgin, char, gas, air),
        reactions = (
            Reaction(
                :pyro,
                1 => (2, 3);
                A = 5e8,
                E = 1.3e5,
                h = 5e5,
                yields = (0.25, 0.75),
            ),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
    )
    L = 0.01
    return PyrolysisProblem(
        mesh = create_stretched_mesh(L, n, Val(4), 1.15),
        material = mat,
        bc_set = BoundaryConditionSet(thermal = thermal, mass = mass),
        T_initial = z -> 500.0 + 3e4 * z,
        ξ_initial = [
            z -> 600.0 * (1 - 0.6 * max(0.0, z / L - 0.4) / 0.6),
            z -> 150.0 * max(0.0, z / L - 0.4) / 0.6,
            z -> 0.3 + 40.0 * z,
            z -> 0.5 - 20.0 * z,
        ],
        tspan = (0.0, 10.0),
    )
end

const _convective_top = ConvectiveBC(h = 15.0, T_inf = 1000.0, blowing = CouetteFilm())
const _gas_configurations = (
    fick_diffusive = (
        darcy = false,
        thermal = (
            top = _convective_top + RadiativeBC(T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.02), bottom = DiffusiveMassBC(Inf)),
    ),
    fick_film = (
        darcy = false,
        thermal = (top = HeatFluxBC(2e4), bottom = ConvectiveBC(h = 5.0, T_inf = 350.0)),
        mass = (top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0), bottom = MassFluxBC(1e-4)),
    ),
    darcy_venting = (
        darcy = true,
        thermal = (
            top = _convective_top + RadiativeBC(T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DirichletPressureBC(101325.0), bottom = ImpermeableBC()),
    ),
    darcy_combined = (
        darcy = true,
        thermal = (top = _convective_top, bottom = DirichletThermalBC(450.0)),
        mass = (
            top = DiffusiveMassBC(0.02) + ConvectivePressureBC(1e-6, 101325.0),
            bottom = NeumannPressureBC(-1e5),
        ),
    ),
)

# States on and off the physical manifold.
function _states(cp)
    u0 = cp.u0
    layout = cp.def.layout
    hot = copy(u0)
    for i = 1:n_cells(layout)
        hot[energy_index(layout, i)] =
            energy_state(layout) isa TemperatureState ? 700.0 + 20i :
            u0[energy_index(layout, i)] * 3
    end
    perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0)))
    zero_state = zeros(length(u0))
    negative = -abs.(u0)
    return (; u0, hot, perturbed, zero_state, negative)
end

@testset "gas transport: $(energy) form, $(n) cells, $(name)" for energy in
                                                                  (:temperature, :enthalpy),
    n in (6, 8),
    name in keys(_gas_configurations)

    cfg = _gas_configurations[name]
    cp = setup(
        _gas_slab(n; darcy = cfg.darcy, thermal = cfg.thermal, mass = cfg.mass);
        energy_state = energy,
    )
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    @test def.setup.gas isa (cfg.darcy ? DarcyFickTransport : FickTransport)
    @test count(k -> k isa GasTransportKernel, def.kernels) == 1
    layout = def.layout
    pattern = cj.pattern
    @test cj.n_colors == 3 * (n_components(layout) + 1)

    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    @test all(isfinite, du)
    residual!(du, cp.u0, def, 0.0, ws)
    @test (@allocated residual!(du, cp.u0, def, 0.0, ws)) == 0

    J = jac_prototype(cj)
    for (state, u) in pairs(_states(cp))
        @testset "$state" begin
            residual!(du, u, def, 1.0, ws)
            @test all(isfinite, du)
            Jd = dense_jacobian(def, ws, u, 1.0)
            @test all(isfinite, Jd)
            @test isempty(missing_entries(Jd, pattern))
            jacobian!(J, cj, u, 1.0)
            cmp = compare_entries(Matrix(J), Jd; rtol = 1e-12)
            @test cmp.passed
        end
    end
end

@testset "gas transport: the mass fluxes reach the gas rows only, and the vent carries gas out" begin
    cfg = _gas_configurations.darcy_venting
    cp = setup(_gas_slab(8; darcy = true, thermal = cfg.thermal, mass = cfg.mass))
    def, ws = cp.def, cp.ws
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    ds, c = buffers(ws, def, Float64)
    @test all(iszero, c.F_mass[1:2, :])
    top = only(boundary_faces(def.setup.topology, :top))
    bottom = only(boundary_faces(def.setup.topology, :bottom))
    # Overpressure against the ambient at the top: gas leaves; sealed back.
    @test ds.P[end] > 101325.0
    @test all(>(0), c.F_mass[3:4, top])
    @test all(iszero, c.F_mass[:, bottom])
    # The char-free back cells are impermeable: no Darcy flow between them.
    @test ds.mob[1] == 0 && ds.mob[2] == 0
end

# Uniform temperature, a nonuniform gas, Darcy flow and no reactions or
# conduction gradient: every face carries the enthalpy of its own mass flux
# at the common temperature, so the temperature row is exactly at rest and
# the inventory row is exactly the enthalpy of the mass change.
@testset "gas transport: isothermal transport leaves the temperature at rest ($energy form)" for energy in
                                                                                                 (
    :temperature,
    :enthalpy,
)

    virgin = SolidComponent(
        :virgin,
        ρ = 600.0,
        c = 1200.0,
        k = 0.15,
        λ = 2e-6,
        κ = 1e-12,
        ρ_intrinsic = 1500.0,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = (1000.0, 0.4), k = 0.03, λ = 3e-5)
    air = GaseousComponent(:air, M = 0.029, c = 1006.0, k = 0.026, λ = 2e-5)
    mat = Material(name = :iso, components = (virgin, gas, air), reactions = ())
    problem = PyrolysisProblem(
        mesh = create_stretched_mesh(0.01, 9, Val(3), 1.2),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = DirichletPressureBC(101325.0), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 650.0,
        ξ_initial = [z -> 600.0, z -> 0.2 + 80.0 * z, z -> 0.6 - 30.0 * z],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem; energy_state = energy)
    def, ws = cp.def, cp.ws
    layout = def.layout
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    ds, c = buffers(ws, def, Float64)
    @test maximum(abs, c.F_mass) > 0
    for i = 1:n_cells(layout)
        rate = du[energy_index(layout, i)]
        if energy === :temperature
            @test abs(rate) <= 1e-12 * 650.0
        else
            dH = sum(du[mass_index(layout, i, j)] * ds.h[j, i] for j = 1:3)
            scale = sum(abs(du[mass_index(layout, i, j)] * ds.h[j, i]) for j = 1:3)
            @test abs(rate - dH) <= 1e-12 * scale
        end
    end
end

@testset "outward Darcy velocity equals the tagged kernel along the outward normal" begin
    mob, d = 2e-7, 3e-4
    for bc in (
            DirichletPressureBC(101325.0),
            NeumannPressureBC(-2e5),
            ConvectivePressureBC(1e-6, 1e5),
            DiffusiveMassBC(0.01),
            DiffusiveMassBC(0.01) + DirichletPressureBC(9e4),
        ),
        P in (8e4, 101325.0, 1.4e5)

        κ, μ = mob * 2e-5, 2e-5
        v_top = boundary_darcy_velocity(bc, P, κ, μ, d, Val(:top))
        v_bottom = boundary_darcy_velocity(bc, P, κ, μ, d, Val(:bottom))
        @test outward_darcy_velocity(bc, P, mob, d) ≈ v_top rtol = 1e-14 atol = 1e-300
        @test outward_darcy_velocity(bc, P, mob, d) ≈ -v_bottom rtol = 1e-14 atol = 1e-300
    end
    # An impermeable cell cannot vent through a film condition; the gate is
    # smooth in the mobility and exactly one at a physical mobility.
    film = ConvectivePressureBC(1e-6, 1e5)
    @test outward_darcy_velocity(film, 2e5, 0.0, d) == 0
    @test outward_darcy_velocity(film, 2e5, 2e-7, d) == 1e-6 * (2e5 - 1e5)
    @test outward_darcy_velocity(film, 2e5, MOBILITY_GUARD, d) ≈ 0.05 rtol = 1e-14
    @test ForwardDiff.derivative(m -> outward_darcy_velocity(film, 2e5, m, d), 0.0) ≈
          0.1 / MOBILITY_GUARD rtol = 1e-14
end

@testset "gas transport: an interior face carries the old core's Fick flux, and gas leaving through :bottom leaves the bottom cell" begin
    thermal = _gas_configurations.fick_diffusive.thermal
    open = setup(
        _gas_slab(
            2;
            darcy = false,
            thermal,
            mass = (top = DiffusiveMassBC(0.02), bottom = DiffusiveMassBC(Inf)),
        ),
    )
    closed = setup(
        _gas_slab(
            2;
            darcy = false,
            thermal,
            mass = (top = DiffusiveMassBC(0.02), bottom = ImpermeableBC()),
        ),
    )
    layout = open.def.layout
    du_open, du_closed = zeros(length(open.u0)), zeros(length(closed.u0))
    residual!(du_open, open.u0, open.def, 0.0, open.ws)
    residual!(du_closed, closed.u0, closed.def, 0.0, closed.ws)
    ds, c = open.ws.derived, open.ws.contributions
    g, topo = ds.geometry, open.def.setup.topology
    d1 = center_to_face_distance(g, topo, 2, 1)
    d2 = center_to_face_distance(g, topo, 2, 2)
    λ_f = harmonic_mean_weighted(ds.λ[1], ds.λ[2], d1, d2)
    for j in (3, 4)
        F = Pyrolysis.Physics.gas_flux(
            λ_f,
            ds.ξ[j, 1],
            ds.ξ[j, 2],
            ds.T_trial[1],
            ds.T_trial[2],
            d1 + d2,
        )
        @test c.F_mass[j, 2] ≈ F rtol = 1e-14
        # Gas leaves through :bottom along its outward normal, and the bottom
        # cell's gas row loses exactly that.
        @test c.F_mass[j, 1] > 0
        k = mass_index(layout, 1, j)
        @test du_open[k] - du_closed[k] ≈ -g.face_area[1] * c.F_mass[j, 1] rtol = 1e-10
    end
end

# A permeable solid with two gases of one molar mass in alternating
# proportions: the pressure is uniform to the bit, so every interior face
# carries `v = 0` between unequal upwind states. With `outlet` the top vents
# against exactly the cells' pressure, so its outward velocity is zero too.
function _pressure_tie(n; outlet::Bool, energy)
    solid = SolidComponent(
        :solid;
        ρ = 600.0,
        ρ_intrinsic = 1200.0,
        c = 1000.0,
        k = 0.1,
        κ = 1e-12,
        λ = 1e-5,
    )
    a = GaseousComponent(:a; M = 0.032, c = 1000.0, k = 0.03, λ = 1e-5)
    b = GaseousComponent(:b; M = 0.032, c = 1000.0, k = 0.03, λ = 1e-5)
    mat = Material(name = :pressure_tie, components = (solid, a, b), reactions = ())
    y(z) = isodd(clamp(floor(Int, z * n) + 1, 1, n)) ? 0.25 : 0.75
    problem(top) = PyrolysisProblem(
        mesh = create_uniform_mesh(1.0, n, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = top, bottom = ImpermeableBC()),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [z -> 600.0, y, z -> 1 - y(z)],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem(ImpermeableBC()); energy_state = energy)
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    P = cp.ws.derived.P
    outlet || return cp, P
    return setup(problem(DirichletPressureBC(first(P))); energy_state = energy), P
end

@testset "gas transport: at zero Darcy velocity every seed takes the v ≥ 0 side ($energy form, $n cells, outlet $outlet)" for energy in
                                                                                                                              (
        :temperature,
        :enthalpy,
    ),
    n in (4, 8),
    outlet in (false, true)

    cp, P = _pressure_tie(n; outlet, energy)
    @test maximum(P) == minimum(P)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    @test def.setup.gas isa DarcyFickTransport
    Jd = dense_jacobian(def, ws, cp.u0, 0.0)
    J = jac_prototype(cj)
    jacobian!(J, cj, cp.u0, 0.0)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
end

@testset "gas transport: the upwind side at v = 0 does not depend on the seed's sign" begin
    for s in (1.0, -1.0)
        F = Pyrolysis.Kernels._advective_flux(ForwardDiff.Dual(0.0, s), 0.5, 0.25, 0.75)
        @test ForwardDiff.value(F) == 0
        @test ForwardDiff.partials(F)[1] == s * 0.25 / 0.5
    end
end
