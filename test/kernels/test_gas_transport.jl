using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, R_GAS
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: GasTransportKernel, outward_darcy_velocity, MOBILITY_GUARD
using Pyrolysis.Kernels: ThermalBoundaryKernel, find_kernel
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
        mesh = create_stretched_mesh(L, n, 1.15),
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

# The other blowing models on a top face hotter than its free stream, with gas
# leaving through it: a film condition (Fick), or a pressure film gentle
# enough that the outflow leaves the blowing factor well above zero
# (Darcy–Fick). The pooling model then solves for its pooling factor.
const _pooling = PoolingFilm(NITROGEN; diffusion_volume = (gas = 60.0,))
_blowing_configuration(model; darcy::Bool) = (
    darcy = darcy,
    thermal = (
        top = ConvectiveBC(h = 15.0, T_inf = 300.0, blowing = model) +
              RadiativeBC(T_inf = 300.0),
        bottom = AdiabaticBC(),
    ),
    mass = (
        top = darcy ? ConvectivePressureBC(1e-8, 101325.0) : DiffusiveMassBC(0.02),
        bottom = ImpermeableBC(),
    ),
)
const _blowing_configurations = (
    fick_couette_film = _blowing_configuration(CouetteFilm(0.31); darcy = false),
    fick_property_ratio = _blowing_configuration(
        PropertyRatioFilm(NITROGEN);
        darcy = false,
    ),
    fick_pooling = _blowing_configuration(_pooling; darcy = false),
    darcy_couette_film = _blowing_configuration(CouetteFilm(0.31); darcy = true),
    darcy_property_ratio = _blowing_configuration(
        PropertyRatioFilm(NITROGEN);
        darcy = true,
    ),
    darcy_pooling = _blowing_configuration(_pooling; darcy = true),
)

const _transport_configurations = (
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
const _gas_configurations = merge(_transport_configurations, _blowing_configurations)

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

@testset "gas transport: $(name) blows through a face hotter than its free stream" for name in
                                                                                       keys(
    _blowing_configurations,
)
    cfg = _blowing_configurations[name]
    cp = setup(_gas_slab(6; darcy = cfg.darcy, thermal = cfg.thermal, mass = cfg.mass))
    residual!(zeros(length(cp.u0)), cp.u0, cp.def, 0.0, cp.ws)
    ds, c = buffers(cp.ws, cp.def, Float64)
    top = only(boundary_faces(cp.def.setup.topology, :top))
    state = find_kernel(ThermalBoundaryKernel, cp.def.kernels).states[top]
    # Gas leaves through the face and the blowing factor stays off its limits,
    # so the model's film law shapes the surface balance.
    @test all(>(0), c.F_mass[3:4, top])
    @test c.T_face[top] > 300.0
    @test 0 < state.blowing_ratio < 1
    # The pooling factor comes from the root solve, which only runs with
    # outflow through a face hotter than the free stream.
    if name in (:fick_pooling, :darcy_pooling)
        @test state.pooling_factor < 1
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
        mesh = create_stretched_mesh(0.01, 9, 1.2),
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

@testset "outward Darcy velocity of each mass condition in closed form" begin
    mob, d = 2e-7, 3e-4
    for P in (8e4, 101325.0, 1.4e5)
        # A prescribed pressure: the one-sided gradient, outward under an
        # overpressure, and nothing out of an impermeable cell.
        @test outward_darcy_velocity(DirichletPressureBC(101325.0), P, mob, d) ==
              mob * (P - 101325.0) / d
        @test outward_darcy_velocity(DirichletPressureBC(101325.0), P, 0.0, d) == 0
        # A prescribed gradient along the outward normal: a pressure falling
        # outward drives outflow, and a zero gradient drives none.
        @test outward_darcy_velocity(NeumannPressureBC(-2e5), P, mob, d) == mob * 2e5
        @test outward_darcy_velocity(NeumannPressureBC(2e5), P, mob, d) == -mob * 2e5
        @test outward_darcy_velocity(NeumannPressureBC(0.0), P, mob, d) == 0
        # A film condition drives none; a combined condition sums its members.
        @test outward_darcy_velocity(DiffusiveMassBC(0.01), P, mob, d) == 0
        combined = DiffusiveMassBC(0.01) + DirichletPressureBC(9e4)
        @test outward_darcy_velocity(combined, P, mob, d) ==
              outward_darcy_velocity(DirichletPressureBC(9e4), P, mob, d)
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

# A uniform, hot char column of porosity one half holding gas well above the
# ambient pressure, its mass conditions `top` and `bottom`: every interior
# flux vanishes, so each boundary cell's gas row is its face's outflow alone.
function _vent_column(κ; top, bottom)
    mat = Material(
        name = :vent_column,
        components = (
            SolidComponent(:char, ρ = 500.0, c = 1000.0, k = 0.1, κ = κ),
            GaseousComponent(:gas, M = 0.028, c = 1000.0, k = 0.03, λ = 1e-6),
        ),
        reactions = (),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 6),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
            mass = (top = top, bottom = bottom),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [z -> 250.0, z -> 2.0],
        tspan = (0.0, 1.0),
    )
end

@testset "gas transport: a pressure condition vents its cell at v ξ/φ beside the film flux, faster through a more permeable solid" begin
    film = DiffusiveMassBC(0.05)
    vent = film + DirichletPressureBC(101325.0)
    function surface(κ, bc)
        cp = setup(_vent_column(κ; top = bc, bottom = ImpermeableBC()))
        du = zeros(length(cp.u0))
        residual!(du, cp.u0, cp.def, 0.0, cp.ws)
        ds, c = buffers(cp.ws, cp.def, Float64)
        topo = cp.def.setup.topology
        f = only(boundary_faces(topo, :top))
        n = n_cells(topo)
        return (; cp, du, ds, c, f, n, d = center_to_face_distance(ds.geometry, topo, f, n))
    end
    (; cp, du, ds, c, f, n, d) = surface(1e-12, vent)
    k = mass_index(cp.def.layout, n, 2)
    # The skeleton fills half the cell; the pressure and the mobility are the
    # ideal gas in the pores and κ/μ.
    @test ds.φ[n] ≈ 0.5 rtol = 1e-12
    P = 2.0 * R_GAS * 600.0 / (0.028 * 0.5)
    mob = 1e-12 / cp.def.setup.material.gas_viscosity(600.0)
    @test ds.P[n] ≈ P rtol = 1e-12
    @test ds.mob[n] ≈ mob rtol = 1e-12
    # The overpressure vents along the outward normal at the one-sided
    # velocity, carrying the intrinsic density ξ/φ = 2ξ, and the film flux
    # λ h ξ/(λ + h d) adds to it.
    v = mob * (P - 101325.0) / d
    @test v > 0
    ξ, λ = ds.ξ[2, n], ds.λ[n]
    J_film = λ * 0.05 * ξ / (λ + 0.05 * d)
    @test c.F_mass[2, f] ≈ 2v * ξ + J_film rtol = 1e-12
    @test du[k] ≈ -ds.geometry.face_area[f] * c.F_mass[2, f] rtol = 1e-10
    # Without the pressure condition only the film law acts.
    @test surface(1e-12, film).c.F_mass[2, f] ≈ J_film rtol = 1e-12
    # The outflow follows the permeability, which the film law alone does not.
    rate = Dict(κ => surface(κ, vent).du[k] for κ in (1e-14, 1e-12, 1e-11))
    @test rate[1e-12] < 0
    @test abs(rate[1e-12]) > 10 * abs(rate[1e-14])
    @test abs(rate[1e-11]) > 5 * abs(rate[1e-12])
end

@testset "gas transport: a pressure condition on :bottom vents the bottom cell along the outward normal" begin
    sealed = setup(_vent_column(1e-12; top = ImpermeableBC(), bottom = ImpermeableBC()))
    du_sealed = zeros(length(sealed.u0))
    residual!(du_sealed, sealed.u0, sealed.def, 0.0, sealed.ws)
    k = mass_index(sealed.def.layout, 1, 2)
    # An overpressure against a prescribed pressure, and a pressure falling
    # outward: both drive gas out through the bottom face.
    for bc in (DirichletPressureBC(101325.0), NeumannPressureBC(-1e5))
        cp = setup(_vent_column(1e-12; top = ImpermeableBC(), bottom = bc))
        du = zeros(length(cp.u0))
        residual!(du, cp.u0, cp.def, 0.0, cp.ws)
        ds, c = buffers(cp.ws, cp.def, Float64)
        f = only(boundary_faces(cp.def.setup.topology, :bottom))
        d = center_to_face_distance(ds.geometry, cp.def.setup.topology, f, 1)
        v = outward_darcy_velocity(bc, ds.P[1], ds.mob[1], d)
        @test v > 0
        @test c.F_mass[2, f] ≈ v * ds.ξ[2, 1] / ds.φ[1] rtol = 1e-12
        @test du[k] - du_sealed[k] ≈ -ds.geometry.face_area[f] * c.F_mass[2, f] rtol = 1e-10
    end
end

@testset "gas transport: an interior face carries the Fick flux of the harmonic-mean transfer coefficient, and gas leaving through :bottom leaves the bottom cell" begin
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
    # The transfer coefficients in series over the two half-cells.
    λ_f = (d1 + d2) / (d1 / ds.λ[1] + d2 / ds.λ[2])
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
        mesh = create_uniform_mesh(1.0, n),
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

# The interior faces of a set-up column after one residual, each gas `j` with
# the parts of its Darcy–Fick flux in closed form: the Fick flux at the trial
# temperatures with the transfer coefficients in series, and the upwind
# advection of the intrinsic density ξ/φ_f at the velocity the mobilities in
# series drive.
function _interior_fluxes(cp, gases)
    residual!(zeros(length(cp.u0)), cp.u0, cp.def, 0.0, cp.ws)
    ds, c = buffers(cp.ws, cp.def, Float64)
    topo, g = cp.def.setup.topology, ds.geometry
    series(a, b, d_a, d_b) =
        iszero(a) || iszero(b) ? zero(a) : (d_a + d_b) / (d_a / a + d_b / b)
    faces = []
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        nb == 0 && continue
        d_o = center_to_face_distance(g, topo, f, o)
        d_n = center_to_face_distance(g, topo, f, nb)
        d = d_o + d_n
        λ_f = series(ds.λ[o], ds.λ[nb], d_o, d_n)
        v = -series(ds.mob[o], ds.mob[nb], d_o, d_n) * (ds.P[nb] - ds.P[o]) / d
        φ_f = (ds.φ[o] + ds.φ[nb]) / 2
        for j in gases
            ξ_o, ξ_n = ds.ξ[j, o], ds.ξ[j, nb]
            fick =
                Pyrolysis.Physics.gas_flux(λ_f, ξ_o, ξ_n, ds.T_trial[o], ds.T_trial[nb], d)
            advection = v * (v >= 0 ? ξ_o : ξ_n) / φ_f
            push!(faces, (; f, j, o, nb, d, φ_f, v, fick, advection, F = c.F_mass[j, f]))
        end
    end
    return faces, ds, c
end

@testset "gas transport: an interior Darcy–Fick face carries the Fick flux plus the upwind advection of ξ/φ, and no advection across the char front" begin
    cfg = _gas_configurations.darcy_venting
    cp = setup(_gas_slab(8; darcy = true, thermal = cfg.thermal, mass = cfg.mass))
    faces, ds, c = _interior_fluxes(cp, (3, 4))
    for r in faces
        @test abs(r.F - (r.fick + r.advection)) <= 1e-12 * (abs(r.fick) + abs(r.advection))
    end
    @test all(iszero, c.F_mass[1:2, :])
    # Between char-bearing cells both parts act.
    flowing = filter(r -> ds.mob[r.o] > 0 && ds.mob[r.nb] > 0, faces)
    @test !isempty(flowing)
    @test all(r -> r.v != 0 && r.fick != 0 && r.advection != 0, flowing)
    # The impermeable virgin cell beside the first char-bearing one stops the
    # Darcy flow through their face: the Fick flux alone, where any mean of
    # the two mobilities that lets the char side through would carry far more.
    front = filter(r -> (ds.mob[r.o] == 0) != (ds.mob[r.nb] == 0), faces)
    @test length(front) == 2
    for r in front
        @test r.v == 0
        @test r.F ≈ r.fick rtol = 1e-12
        leaky = -(ds.mob[r.o] + ds.mob[r.nb]) / 2 * (ds.P[r.nb] - ds.P[r.o]) / r.d
        ξ_up = leaky >= 0 ? ds.ξ[r.j, r.o] : ds.ξ[r.j, r.nb]
        @test abs(leaky * ξ_up / r.φ_f) > 10 * abs(r.fick)
    end
    # At a pressure uniform to the bit no face carries a velocity: the Fick
    # flux alone between unequal gas concentrations.
    tie, _ = _pressure_tie(8; outlet = false, energy = :temperature)
    faces, _, _ = _interior_fluxes(tie, (2, 3))
    @test all(r -> r.v == 0 && r.fick != 0, faces)
    @test all(r -> isapprox(r.F, r.fick; rtol = 1e-14), faces)
end

# A sealed, isothermal Darcy–Fick column with a pressure peak in the middle:
# the gas flows out of the peak both ways and stays in the column.
function _sealed_column()
    L = 0.01
    solid = SolidComponent(
        :solid,
        ρ = 1000.0,
        c = 1500.0,
        k = 0.3,
        κ = 1e-11,
        ρ_intrinsic = 1500.0,
    )
    gas = GaseousComponent(:gas, M = 0.044, c = 1000.0, k = 0.025, λ = 1e-5)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, 20),
        material = Material(
            name = :sealed_column,
            components = (solid, gas),
            reactions = (),
        ),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 800.0, z -> 0.1 + 2.0 * exp(-((z / L - 0.5) / 0.1)^2)],
        tspan = (0.0, 1.0),
        output_times = [0.0, 0.5, 1.0],
    )
end

@testset "gas transport: a sealed Darcy–Fick column keeps its gas and its solid ($energy form)" for energy in
                                                                                                    (
    :temperature,
    :enthalpy,
)
    problem = _sealed_column()
    cp = setup(problem; energy_state = energy)
    layout, topo = cp.def.layout, cp.def.setup.topology
    @test cp.def.setup.gas isa DarcyFickTransport
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, cp.def, 0.0, cp.ws)
    ds, c = buffers(cp.ws, cp.def, Float64)
    # A permeable solid under a steep pressure gradient: the gas flows.
    @test all(>(0), ds.mob)
    @test maximum(ds.P) > 10 * minimum(ds.P)
    @test maximum(abs, c.F_mass) > 0
    # The sealed faces carry nothing, the solid does not move, and what the
    # gas rows lose and gain sums to zero.
    for tag in (:top, :bottom)
        @test all(iszero, c.F_mass[:, only(boundary_faces(topo, tag))])
    end
    @test all(i -> du[mass_index(layout, i, 1)] == 0, 1:n_cells(layout))
    rates = [du[mass_index(layout, i, 2)] for i = 1:n_cells(layout)]
    @test abs(sum(rates)) <= 1e-12 * sum(abs, rates)
    # Over a solve the gas inventory stays where it started.
    sol = Pyrolysis.Integration.solve(
        problem;
        energy_state = energy,
        reltol = 1e-8,
        abstol = 1e-10,
    )
    @test sol.retcode === :Success
    gas = [sum(sol.m[2, :, k]) for k in eachindex(sol.t)]
    @test maximum(abs, gas .- gas[1]) <= 1e-12 * gas[1]
end

# Two sealed isothermal cells at 600 K with unequal gas concentrations, the gas
# carrying the transfer law `λ`.
function _sealed_pair(λ)
    solid = SolidComponent(:solid; ρ = 600.0, c = 1000.0, k = 0.1)
    gas = GaseousComponent(:gas; M = 0.032, c = 1000.0, k = 0.03, λ)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(1.0, 2),
        material = Material(name = :sealed_pair, components = (solid, gas), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 600.0,
        ξ_initial = [z -> 600.0, z -> z < 0.5 ? 1.0 : 2.0],
        tspan = (0.0, 1.0),
    )
end

@testset "gas transport: a transfer law zero at sampled temperatures still transports" begin
    rates = map((
        TableProperty((200.0, 300.0, 600.0, 1000.0, 3000.0), (0.0, 0.0, 1e-5, 0.0, 0.0)),
        1e-5,
    )) do λ
        cp = setup(_sealed_pair(λ))
        @test cp.def.setup.gas isa FickTransport
        du = zeros(length(cp.u0))
        residual!(du, cp.u0, cp.def, 0.0, cp.ws)
        [du[mass_index(cp.def.layout, i, 2)] for i = 1:2]
    end
    @test rates[1] ≈ rates[2] rtol = 1e-12
    @test rates[1][1] > 0
    @test setup(_sealed_pair(0.0)).def.setup.gas isa NoGasTransport
end

# The Carman–Kozeny permeability reads the state mapping's smoothly floored
# porosity, so a trial state whose skeleton overfills its cell keeps the
# positive mobility of the floor, joined C¹ to the physical states.
@testset "Carman–Kozeny permeability at the smooth porosity" begin
    ℓ = 1e-4
    solid =
        SolidComponent(:s, ρ = 500.0, c = 1500.0, k = 0.2, λ = 1e-5, ρ_intrinsic = 1000.0)
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.05, λ = 1e-5)
    mat = Material(;
        name = :ck,
        components = (solid, gas),
        reactions = (),
        mixing = MaterialMixing(κ = MixingSpec(CARMAN_KOZENY; length_scale = ℓ)),
    )
    T = 400.0
    for ξ in ((500.0, 0.1), (900.0, 0.1))
        φ = 1 - ξ[1] / 1000.0
        @test Pyrolysis.Properties.effective_permeability(ξ, T, mat, φ) ≈
              min(ℓ^2 * φ^3 / (180 * (1 - φ)^2), ℓ^2) rtol = 1e-14
    end
    # Overfilled: the floored porosity keeps a positive permeability.
    κ(s) = Pyrolysis.Properties.effective_permeability(
        (s, oftype(s, 0.1)),
        T,
        mat,
        Pyrolysis.States.smooth_floor(1 - s / 1000.0, 1e-12),
    )
    @test κ(1100.0) > 0
    @test κ(1100.0) == κ(1000.0 * (1 + 1e-9))
    @test ForwardDiff.derivative(κ, 1100.0) == 0
    # Inside the floor's band the permeability still falls with the solid.
    band = 1000.0 * (1 - 1.5e-12)
    @test ForwardDiff.derivative(κ, band) < 0
    @test isfinite(ForwardDiff.derivative(κ, band))
end
