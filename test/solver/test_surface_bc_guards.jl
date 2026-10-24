"""
Regression tests for BC / surface-Newton and radiation guards. Each testset
exercises the failing regime it guards against.
"""

using Test
using Pyrolysis
using Pyrolysis.Problem: validate_problem
using Pyrolysis.Residual: residual!
using Pyrolysis.Materials: SURFACE_ABSORPTION, BEER_LAMBERT, NO_RADIATION
using Pyrolysis.BoundaryConditions: solve_boundary_temperature
using Pyrolysis.BoundaryConditions: ambient_concentration

function _minors_material(; α::Real = 0.0)
    solid = SolidComponent(:solid, ρ = 600.0, c = 1500.0, k = 0.15, ε = 0.9, α = α)
    gas = GaseousComponent(:gas, M = 0.05, c = 1200.0, k = 0.03, λ = 1e-5)
    return Material(name = :minors, components = (solid, gas), reactions = ())
end

# ============================================================================
# Surface Newton non-convergence must warn (never silently fall through to
# the final flux evaluation at the unconverged iterate)
# ============================================================================
@testset "Newton non-convergence warns and is recorded" begin
    mat = _minors_material()
    bc = RadiativeBC(T_inf = 300.0) + ConvectiveBC(h = 10.0, T_inf = 300.0)
    J0 = (0.0, 0.0)

    # max_iter = 1 cannot converge the cubic radiation balance from a cold
    # guess: the solve must warn and record it in the `converged` flag.
    state =
        @test_logs (:warn, r"did not converge") match_mode = :any solve_boundary_temperature(
            800.0,
            0.15,
            1e-3,
            bc,
            mat,
            J0,
            0.0,
            0.9;
            max_iter = 1,
        )
    @test state.converged == false
    @test state.iterations == 1

    # The same call with the default budget converges silently.
    state_ok = solve_boundary_temperature(800.0, 0.15, 1e-3, bc, mat, J0, 0.0, 0.9)
    @test state_ok.converged
end

# ============================================================================
# The residual carries no surface-Newton state between evaluations: after a
# wild trial state (the kind a diverging outer Newton produces), the residual
# of a good state equals its residual on a fresh workspace
# ============================================================================
@testset "residual is independent of earlier evaluations" begin
    mat = _minors_material()
    bc = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = 50e3) +
                  RadiativeBC(T_inf = 300.0) +
                  ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    n = 6
    fresh_setup(use_ale) = initialize_state(
        PyrolysisProblem(
            mesh = create_uniform_mesh(0.01, n, Val(2)),
            material = mat,
            bc_set = bc,
            tspan = (0.0, 1.0),
        );
        use_ale = use_ale,
    )
    for use_ale in (false, true)
        def, ws, u = fresh_setup(use_ale)
        for i = 1:n
            u[i] = 300.0 + 60.0 * i
            u[2n+i] = 0.1 * i
        end
        du_fresh = similar(u)
        def_f, ws_f, _ = fresh_setup(use_ale)
        residual!(du_fresh, u, def_f, 0.0, ws_f)
        @test all(isfinite, du_fresh)

        du = similar(u)
        for T_wild in (-9.1e40, NaN)
            u_wild = copy(u)
            u_wild[n] = T_wild
            residual!(similar(u), u_wild, def, 0.0, ws)
            residual!(du, u, def, 0.0, ws)
            @test isapprox(du, du_fresh; rtol = 1e-12)
        end
    end
end

# ============================================================================
# A final iterate pinned at the [200, 3000] K clamp warns (the balance has
# no root in the bracket; AD sensitivities of T_s are zero there)
# ============================================================================
@testset "clamped final iterate warns" begin
    mat = _minors_material()
    # A huge prescribed flux against weak conduction pushes the root far above
    # 3000 K: every iterate is clamped, including the last.
    bc = HeatFluxBC(1e9)
    J0 = (0.0, 0.0)

    state =
        @test_logs (:warn, r"pinned at the physical clamp") match_mode = :any solve_boundary_temperature(
            300.0,
            0.15,
            1e-3,
            bc,
            mat,
            J0,
            0.0,
            0.9,
        )
    @test state.T_surface == 3000.0
end

# ============================================================================
# Mass BCs write only gas rows of the boundary-face flux matrix — the solid
# rows must never carry the phantom film-law flux
# ============================================================================
@testset "no solid boundary flux in J" begin
    mat = _minors_material()
    mesh = create_uniform_mesh(0.01, 6, Val(2))
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 1.0))
    def, ws, u = initialize_state(problem; use_ale = false)
    n = 6
    for i = 1:n
        u[i] = 600.0
        u[n+i] = 400.0   # solid — the film law must not write h_eff·ξ_solid ≈ O(10) into J
        u[2n+i] = 2.0    # gas
    end
    residual!(similar(u), u, def, 0.0, ws)

    J = ws.rhs_cache.J_gas
    @test all(iszero, @view J[1, :])      # solid row untouched everywhere
    @test any(!iszero, @view J[2, :])     # gas row carries the film-law venting
end

# ============================================================================
# ConvectiveMassBC ξ_inf is species-resolvable
# ============================================================================
@testset "per-species ambient concentrations" begin
    # Scalar broadcasts to every species, tuple resolves per species,
    # function resolves in time.
    bc_scalar = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.3)
    @test ambient_concentration(bc_scalar, 0.0, 1) == 0.3
    @test ambient_concentration(bc_scalar, 0.0, 2) == 0.3

    bc_tuple = ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.27))
    @test ambient_concentration(bc_tuple, 0.0, 1) == 0.0
    @test ambient_concentration(bc_tuple, 0.0, 2) == 0.27

    bc_fun = ConvectiveMassBC(h_m = 0.05, ξ_inf = t -> 0.1 * (1 + t))
    @test ambient_concentration(bc_fun, 1.0, 2) ≈ 0.2

    # Negative ambient concentrations are rejected at construction.
    @test_throws ArgumentError ConvectiveMassBC(h_m = 0.05, ξ_inf = (-0.1, 0.0))

    # Tuple length must match NC — validated at problem validation.
    mat = _minors_material()
    mesh = create_uniform_mesh(0.01, 4, Val(2))
    bad_bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (
            top = ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.1, 0.2)),
            bottom = ImpermeableBC(),
        ),
    )
    bad = PyrolysisProblem(mesh = mesh, material = mat, bc_set = bad_bc, tspan = (0.0, 1.0))
    @test_throws Exception validate_problem(bad)

    # Physical check: ambient gas (species 2 only) flows INTO a column that
    # starts gas-free — a regime a scalar ξ_inf cannot express (a scalar
    # would inject "ambient solid" through the same film law).
    mesh2 = create_uniform_mesh(0.01, 6, Val(2))
    bc2 = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (
            top = ConvectiveMassBC(h_m = 0.05, ξ_inf = (0.0, 0.5)),
            bottom = ImpermeableBC(),
        ),
    )
    problem2 =
        PyrolysisProblem(mesh = mesh2, material = mat, bc_set = bc2, tspan = (0.0, 1.0))
    def2, ws2, u2 = initialize_state(problem2; use_ale = false)
    n = 6
    for i = 1:n
        u2[i] = 400.0
        u2[n+i] = 400.0
        u2[2n+i] = 0.0   # gas-free: ambient ingress must be positive
    end
    du2 = similar(u2)
    residual!(du2, u2, def2, 0.0, ws2)
    @test du2[2n+n] > 0          # gas enters the surface cell
    @test du2[n+n] == 0          # solid sees no "ambient solid" injection
end

# ============================================================================
# mass_loss_rate counts the bottom film flux — ∫MLR dt must close against
# the domain mass for a permeable back face
# ============================================================================
@testset "MLR includes bottom venting" begin
    mat = _minors_material()
    mesh = create_uniform_mesh(0.01, 6, Val(2))
    # Vent the initial gas inventory through BOTH faces; an MLR that counted
    # only the top face would under-report the mass loss ~2×.
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.05), bottom = DiffusiveMassBC(0.05)),
    )
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc,
        tspan = (0.0, 40.0),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 400.0, z -> 1.0],
        output_times = collect(0.0:0.2:40.0),
    )
    sol = solve(problem)
    @test sol.retcode == :Success

    t = sol.t
    mlr = sol.mass_loss_rate
    # Bottom face actually vents: at t=0 the symmetric setup loses through
    # both faces, so MLR must exceed the top-face film law alone.
    m_lost = sol.extras.total_mass[1] - sol.extras.total_mass[end]
    @test m_lost > 0
    ∫mlr = sum(0.5 * (mlr[i] + mlr[i+1]) * (t[i+1] - t[i]) for i = 1:(length(t)-1))
    # Top face only would give ∫MLR dt ≈ m_lost/2 (symmetric venting) — well
    # below the 0.8·m_lost floor. The upper bound allows the trapezoidal
    # quadrature's overshoot on the convex early decay.
    @test 0.8 * m_lost < ∫mlr < 1.3 * m_lost
end

# ============================================================================
# Standalone DirichletThermalBC updates the boundary state and reports the
# conducted flux (no NaN poisoning, no 300 K placeholder)
# ============================================================================
@testset "Dirichlet boundary diagnostics" begin
    mat = _minors_material()
    mesh = create_uniform_mesh(0.01, 6, Val(2))
    bc = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(600.0), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc,
        tspan = (0.0, 5.0),
        output_times = [0.0, 2.5, 5.0],
    )
    sol = solve(problem)
    @test sol.retcode == :Success

    # surface_T and surface_heat_flux_* report the solved Dirichlet face (no
    # constructor placeholder, no NaN).
    # (The solved T_s sits ~q/1e8 ≈ 5e-4 K below the target: the Newton path
    # enforces Dirichlet through the h = 1e8 penalty.)
    @test all(T -> isapprox(T, 600.0; atol = 0.01), sol.surface_T)
    @test all(isfinite, sol.extras.surface_heat_flux_net)
    # The column starts at 300 K below a 600 K face: heat flows in.
    @test sol.extras.surface_heat_flux_net[2] > 0
end

# ============================================================================
# Bottom RadiativeFluxBC under BEER_LAMBERT errors loudly
# ============================================================================
@testset "BEER_LAMBERT rejects bottom RadiativeFluxBC" begin
    mat = _minors_material(α = 200.0)
    mesh = create_uniform_mesh(0.01, 6, Val(2))
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(0.0), bottom = RadiativeFluxBC(flux = 10e3)),
        mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
    )
    problem = PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 1.0))
    @test_throws ArgumentError solve(problem; radiation_model = BEER_LAMBERT)

    # SURFACE_ABSORPTION (the suggested alternative) still accepts it.
    sol = solve(problem; radiation_model = SURFACE_ABSORPTION)
    @test sol.retcode == :Success
end

# ============================================================================
# A time-varying radiation flux that happens to agree at t = 0 and t = 1 is
# NOT frozen at I(0) — q₀(1 + sin 2πt) is exactly such a flux
# ============================================================================
@testset "time-varying flux is re-sampled every RHS" begin
    mat = _minors_material(α = 500.0)
    bc_of(flux) = BoundaryConditionSet(
        thermal = (
            top = RadiativeFluxBC(flux = flux) + ConvectiveBC(h = 10.0, T_inf = 300.0),
            bottom = AdiabaticBC(),
        ),
        mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
    )
    q0 = 30e3
    prob_of(flux) = PyrolysisProblem(
        mesh = create_uniform_mesh(0.005, 8, Val(2)),
        material = mat,
        bc_set = bc_of(flux),
        tspan = (0.0, 1.0),
        output_times = [0.0, 0.5, 1.0],
    )

    # sin(2π·0) = sin(2π·1) = 0: a time-dependence check that samples only
    # t = 0 and t = 1 sees equal values and would freeze the volumetric
    # deposit at I(0) = q₀, making the varying run identical to the
    # constant-q₀ run.
    sol_var = solve(prob_of(t -> q0 * (1 + sin(2π * t))); radiation_model = BEER_LAMBERT)
    sol_const = solve(prob_of(q0); radiation_model = BEER_LAMBERT)
    @test sol_var.retcode == :Success && sol_const.retcode == :Success

    # By t = 0.5 the varying flux has deposited ∫₀^½ q₀ sin(2πt) dt = q₀/π
    # extra energy: the surface must be measurably hotter than the constant
    # run.
    i_mid = findfirst(≈(0.5), sol_var.t)
    @test sol_var.surface_T[i_mid] > sol_const.surface_T[i_mid] + 1.0
end
