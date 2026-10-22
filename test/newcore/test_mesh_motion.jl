using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, BEER_LAMBERT
using Pyrolysis.Meshes
using Pyrolysis.Meshes: MeshGeometry
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: Contributions
using Pyrolysis.Closures
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration
using Pyrolysis.Integration: setup

# A decomposing slab whose char is lighter than the virgin solid, so the
# reaction changes the volume; hot enough at the front that the strain rates
# are large. With `law` the material also shrinks laterally with conversion.
function _slab(
    n;
    law = nothing,
    top = HeatFluxBC(3e4),
    bottom = DirichletThermalBC(400.0),
    radiating = false,
    λ = 0.0,
    κ = 0.0,
    mass = nothing,
    coordinates = Planar(),
)
    virgin = SolidComponent(
        :virgin,
        ρ = 1100.0,
        c = (900.0, 1.5),
        k = (0.2, 2e-4),
        ε = (0.8, 5e-5),
        absorptivity = 0.7,
        α = radiating ? 1.5 : 0.0,
        λ = λ / 10,
        ρ_intrinsic = 1400.0,
    )
    char = SolidComponent(
        :char,
        ρ = 300.0,
        c = PolynomialProperty((700.0, 0.8, -1e-4)),
        k = 0.1,
        ε = 0.95,
        α = radiating ? 20.0 : 0.0,
        λ = λ,
        κ = κ,
        ρ_intrinsic = 1800.0,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = λ)
    mat = Material(
        name = :slab,
        components = (virgin, char, gas),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
        lateral_shrinkage_law = law,
    )
    bcs =
        mass === nothing ? BoundaryConditionSet(thermal = (top = top, bottom = bottom)) :
        BoundaryConditionSet(thermal = (top = top, bottom = bottom), mass = mass)
    return PyrolysisProblem(
        mesh = create_stretched_mesh(0.01, n, 1.15; coordinates),
        material = mat,
        bc_set = bcs,
        T_initial = z -> 500.0 + 4e4 * z,
        ξ_initial = [z -> 1100.0 * (1 - 30z), z -> 300.0 * 30z, z -> 0.5],
        tspan = (0.0, 10.0),
    )
end

const _configurations = (
    planar = (problem = _slab, options = (;)),
    surface_balance = (
        problem = n -> _slab(
            n;
            top = HeatFluxBC(3e4) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0) +
                  RadiativeBC(T_inf = 300.0),
            bottom = ConvectiveBC(h = 5.0, T_inf = 350.0),
        ),
        options = (;),
    ),
    cylinder = (problem = n -> _slab(n; coordinates = Cylindrical()), options = (;)),
    sphere = (problem = n -> _slab(n; coordinates = Spherical()), options = (;)),
    chi = (problem = n -> _slab(n; law = χ -> 1 - 0.4χ), options = (;)),
    beer_lambert = (
        problem = n -> _slab(
            n;
            radiating = true,
            top = RadiativeFluxBC(flux = t -> 3e4 * (1 + 0.5 * sin(t))) +
                  RadiativeBC(T_inf = 300.0),
        ),
        options = (radiation_model = BEER_LAMBERT,),
    ),
    fick = (
        problem = n -> _slab(
            n;
            λ = 3e-5,
            top = ConvectiveBC(h = 15.0, T_inf = 1000.0, blowing = CouetteFilm()),
            mass = (top = DiffusiveMassBC(0.02), bottom = ImpermeableBC()),
        ),
        options = (;),
    ),
    darcy = (
        problem = n -> _slab(
            n;
            λ = 3e-5,
            κ = 1e-11,
            law = χ -> 1 - 0.3χ,
            top = ConvectiveBC(h = 15.0, T_inf = 1000.0),
            mass = (top = DirichletPressureBC(101325.0), bottom = ImpermeableBC()),
        ),
        options = (;),
    ),
)

# A conversion history that is not zero, so the lateral law and its
# derivatives bite: every cell has released some of its dry solid as gas.
function _with_history(u, layout)
    v = copy(u)
    has_chi(layout) || return v
    for i = 1:n_cells(layout)
        v[history_index(layout, i, 1)] = (0.2 + 0.05i) * v[history_index(layout, i, 2)]
    end
    return v
end

# States on and off the physical manifold.
function _states(cp)
    layout = cp.def.layout
    u0 = _with_history(cp.u0, layout)
    hot = copy(u0)
    for i = 1:n_cells(layout)
        hot[energy_index(layout, i)] =
            energy_state(layout) isa TemperatureState ? 700.0 + 20i :
            u0[energy_index(layout, i)] * 3
    end
    perturbed = u0 .* (1 .+ 0.2 .* sin.(1:length(u0)))
    return (; u0, hot, perturbed, zero_state = zeros(length(u0)), negative = -abs.(u0))
end

@testset "moving mesh: $(name), $(energy) form, $(n) cells" for name in
                                                                keys(_configurations),
    energy in (:temperature, :enthalpy),
    n in (6, 8)

    config = _configurations[name]
    cp = setup(config.problem(n); energy_state = energy, use_ale = true, config.options...)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    layout = def.layout
    @test is_moving(layout)
    @test same_structure(jac_prototype(cj), cj.pattern)
    @test all(cj.pattern[r, r] for r in axes(cj.pattern, 1))

    # The start is consistent: every algebraic row is satisfied.
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, def, 0.0, ws)
    @test all(isfinite, du)
    w = [cp.u0[velocity_index(layout, nd)] for nd = 1:n_nodes(layout)]
    @test maximum(w) > 0
    @test maximum(abs(du[velocity_index(layout, nd)]) for nd = 1:n_nodes(layout)) <=
          1e-12 * maximum(abs, w)
    if has_chi(layout)
        @test abs(du[area_index(layout)]) <= 1e-14
        @test abs(du[lateral_rate_index(layout)]) <= 1e-14
    end
    @test (@allocated residual!(du, cp.u0, def, 0.0, ws)) == 0

    J = jac_prototype(cj)
    for (state, u) in pairs(_states(cp))
        @testset "$state" begin
            if has_chi(layout) && state in (:u0, :hot, :perturbed)
                ΣG = sum(u[history_index(layout, i, 1)] for i = 1:n)
                @test ΣG > 0
            end
            residual!(du, u, def, 1.0, ws)
            @test all(isfinite, du)
            Jd = dense_jacobian(def, ws, u, 1.0)
            @test all(isfinite, Jd)
            @test isempty(missing_entries(Jd, cj.pattern))
            jacobian!(J, cj, u, 1.0)
            @test same_structure(J, cj.pattern)
            @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
            if state === :hot
                # The velocity rows read the cells' state through the strain
                # rate, and under a variable cross section the lateral rate.
                v = velocity_index(layout, n_nodes(layout))
                @test Jd[v, energy_index(layout, n)] != 0
                @test Jd[v, mass_index(layout, n, 1)] != 0
                has_chi(layout) && @test Jd[v, lateral_rate_index(layout)] != 0
                has_chi(layout) &&
                    @test Jd[lateral_rate_index(layout), energy_index(layout, 1)] != 0
            end
        end
    end
end

# A slab that converts at the same rate everywhere: isothermal, uniform, the
# activation energy zero. Its strain rate is uniform.
function _uniform_slab(n; law = nothing, ρ_b = 400.0, coordinates = Planar())
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    b = SolidComponent(:b, ρ = ρ_b, c = 1000.0, k = 0.5)
    g = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 0.0)
    mat = Material(
        name = :uniform,
        components = (a, b, g),
        reactions = (
            Reaction(:r, 1 => (2, 3); A = 0.5, E = 0.0, h = 0.0, yields = (0.5, 0.5)),
        ),
        depletion_limiter = DepletionLimiter(enabled = false),
        lateral_shrinkage_law = law,
    )
    return PyrolysisProblem(
        mesh = create_stretched_mesh(0.02, n, 1.1; coordinates),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(400.0), bottom = DirichletThermalBC(400.0)),
        ),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 1000.0, z -> 0.0, z -> 0.0],
        tspan = (0.0, 4.0),
    )
end

@testset "balance rows: a uniform strain rate gives θ x, θ r/2, θ r/3" begin
    for (coords, d) in ((Planar(), 1), (Cylindrical(), 2), (Spherical(), 3))
        cp = setup(_uniform_slab(7; coordinates = coords); use_ale = true)
        layout = cp.def.layout
        θ = cp.ws.contributions.θ
        @test all(≈(θ[1]; rtol = 1e-13), θ)
        @test θ[1] > 0
        for nd = 1:n_nodes(layout)
            x = cp.u0[position_index(layout, nd)]
            @test cp.u0[velocity_index(layout, nd)] ≈ θ[1] * x / d rtol = 1e-12 atol = 1e-18
        end
    end
    # A lateral strain rate takes its share: the axial velocity is (θ − θ_A) x.
    cp = setup(_uniform_slab(7; law = χ -> 1 - 0.5χ); use_ale = true)
    layout = cp.def.layout
    θ, θ_A = cp.ws.contributions.θ[1], cp.u0[lateral_rate_index(layout)]
    @test θ_A < 0
    # dA/dt = A θ_A: the conversion rate is the dry gas yield of the reaction rate.
    @test θ_A ≈ -0.5 * (0.5 * 0.5) rtol = 1e-12
    for nd = 1:n_nodes(layout)
        x = cp.u0[position_index(layout, nd)]
        @test cp.u0[velocity_index(layout, nd)] ≈ (θ - θ_A) * x rtol = 1e-12 atol = 1e-18
    end
end

@testset "closure: rank, consistent sweep and refusals" begin
    mesh = line_mesh(collect(range(0.0, 1.0, length = 6)))
    topo = mesh.topology
    g = MeshGeometry{Float64}(topo)
    geometry!(g, [p[1] for p in mesh.positions], 1.0, topo)
    cl = BalanceRowClosure(topo)
    @test rank_check(cl, topo, g)
    @test cl.row_node == collect(2:6)
    @test cl.order == collect(1:5)
    @test rank_check(
        BalanceRowClosure(topo; correction = PrescribedCorrection((x, t) -> x)),
        topo,
        g,
    )
    @test_throws ArgumentError BalanceRowClosure(
        line_mesh([0.0, 1.0]; fixed = true).topology,
    )
    # The sphere's center node has zero face area; it is the pinned node, so
    # every balance row keeps a nonzero coefficient.
    sphere = line_mesh(collect(range(0.0, 1.0, length = 6)); coordinates = Spherical())
    gs = MeshGeometry{Float64}(sphere.topology)
    geometry!(gs, collect(range(0.0, 1.0, length = 6)), 4π, sphere.topology)
    @test rank_check(BalanceRowClosure(sphere.topology), sphere.topology, gs)

    reacting = _uniform_slab(4)
    @test_throws ArgumentError setup(reacting; use_ale = true, w_override = (x, t) -> 0.0)
    @test_throws ArgumentError setup(reacting; w_override = (x, t) -> 0.0)
    @test_throws ArgumentError setup(_uniform_slab(4; law = χ -> 1 - 0.5χ))
    @test_throws ArgumentError setup(
        _uniform_slab(4; law = χ -> 1 - 0.5χ, coordinates = Cylindrical());
        use_ale = true,
    )
    @test_throws ArgumentError setup(
        _uniform_slab(4; law = χ -> 0.9 - 0.5χ);
        use_ale = true,
    )
    @test_throws ArgumentError setup(_uniform_slab(4; law = χ -> 1 - 1.5χ); use_ale = true)
end

@testset "lateral rows: the lateral rate is the area row's logarithmic rate" begin
    law = χ -> 1 - 0.4χ + 0.1χ^2
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    g = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 0.0)
    mat = Material(
        name = :m,
        components = (a, g),
        reactions = (),
        lateral_shrinkage_law = law,
    )
    lat = LateralShrinkage(mat, 2.0, 3.0)
    @test lat.dry_gas == (false, true)
    # Along a history path, the area row's derivative is A θ_A.
    path(t) = (0.6 + 0.3t, 3.0 + 0.1t)
    t = 0.4
    ΣG, ΣM = path(t)
    ΣĠ, ΣṀ = 0.3, 0.1
    A, θ_A = consistent_lateral(lat, ΣG, ΣM, ΣĠ, ΣṀ)
    @test A ≈ 2.0 * law(ΣG / ΣM)
    area(t) = consistent_lateral(lat, path(t)..., ΣĠ, ΣṀ)[1]
    @test ForwardDiff.derivative(area, t) ≈ A * θ_A rtol = 1e-12
    r = lateral_rows(lat, Pyrolysis.Closures.SVector(ΣG, ΣM, ΣĠ, ΣṀ, θ_A, A))
    @test maximum(abs, r) <= 1e-14
    # Trial states: the conversion is clamped into the law's band and the
    # dry-solid sum floored, so the rows stay finite with finite derivatives.
    for s in (
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        (-5.0, -1.0, 2.0, -3.0, 1.0, -1.0),
        (9.0, 1.0, 0.0, 0.0, 0.0, 1.0),
    )
        x = Pyrolysis.Closures.SVector(s...)
        @test all(isfinite, lateral_rows(lat, x))
        @test all(isfinite, ForwardDiff.jacobian(y -> lateral_rows(lat, y), x))
    end
    @test_throws ArgumentError LateralShrinkage(mat, 2.0, 0.0)
    no_law = Material(name = :m, components = (a, g), reactions = ())
    @test_throws ArgumentError LateralShrinkage(no_law, 2.0, 3.0)
end

@testset "solve: the swelling column's thickness follows its conversion" begin
    # Uniform first-order conversion of a into b (lighter) and g: the volume
    # change is linear in the converted mass, so L(t)/L₀ = 1 + (V_f/V₀ − 1) X(t)
    # with X = 1 − e^{−kt} and V_f/V₀ = 0.5 ρ_a/ρ_b.
    for energy in (:temperature, :enthalpy)
        problem = _uniform_slab(6)
        sol = Integration.solve(
            problem;
            energy_state = energy,
            use_ale = true,
            abstol = 1e-10,
            reltol = 1e-9,
            saveat = [0.0, 1.0, 2.0, 4.0],
        )
        @test sol.retcode === :Success
        L0 = sol.extras.thickness_history[1]
        @test L0 ≈ 0.02
        for (k, t) in enumerate(sol.t)
            X = 1 - exp(-0.5t)
            @test sol.extras.thickness_history[k] / L0 ≈ 1 + (0.5 * 1000 / 400 - 1) * X rtol =
                1e-7
        end
        @test sol.extras.total_mass ≈ fill(sol.extras.total_mass[1], 4) rtol = 1e-12
        @test sol.extras.cross_section_area_history == fill(1.0, 4)
        @test size(sol.z) == (7, 4)
        @test sol.z[1, :] == zeros(4)
    end
    # A shrinking column ends its solve when it thins below min_thickness.
    sol = Integration.solve(
        _uniform_slab(6; ρ_b = 4000.0);
        use_ale = true,
        min_thickness = 0.01,
        abstol = 1e-10,
        reltol = 1e-9,
        saveat = Float64[],
    )
    @test sol.retcode === :Terminated
    @test sol.extras.thickness_history[end] ≈ 0.01 rtol = 1e-6
    @test sol.t[end] ≈ -log(1 - 0.5 / (1 - 0.5 * 1000 / 4000)) / 0.5 rtol = 1e-6
end
