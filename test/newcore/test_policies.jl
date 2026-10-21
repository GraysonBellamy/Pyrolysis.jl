using Test
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter, BEER_LAMBERT, SURFACE_ABSORPTION
using Pyrolysis.States
using Pyrolysis.Integration
using Pyrolysis.Integration: setup, kencarp4_policy, rodas5p_policy, fbdf_policy
using Pyrolysis.Integration: DomainCheck, NewtonCap, cap_step!, FAILURE_MEMORY_GROWTH
using Pyrolysis.Integration: FAILURE_MEMORY_STEPS
using Pyrolysis.Assembly: residual!
import OrdinaryDiffEq
using OrdinaryDiffEq: ODEProblem, DiscreteCallback, CallbackSet, init, step!
using OrdinaryDiffEq: derivative_discontinuity!
using Pyrolysis.Solver: KenCarp4
using OrdinaryDiffEq.SciMLBase: get_proposed_dt
using LinearSolve: LinearSolve, KLUFactorization, UMFPACKFactorization
using SparseArrays

# A decomposing slab whose char is lighter than the virgin solid (the mesh
# moves), porous with intrinsic densities, permeable when `κ > 0`, and
# absorbing in depth when `radiating`.
function _slab(n; κ = 0.0, radiating = false)
    virgin = SolidComponent(
        :virgin,
        ρ = 1100.0,
        c = (900.0, 1.5),
        k = (0.2, 2e-4),
        ε = 0.8,
        α = radiating ? 1.5 : 0.0,
        λ = 3e-6,
        ρ_intrinsic = 1400.0,
    )
    char = SolidComponent(
        :char,
        ρ = 300.0,
        c = 1100.0,
        k = 0.1,
        ε = 0.95,
        α = radiating ? 20.0 : 0.0,
        λ = 3e-5,
        κ = κ,
        ρ_intrinsic = 1800.0,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03, λ = 3e-5)
    mat = Material(
        name = :slab,
        components = (virgin, char, gas),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
    )
    top = radiating ? RadiativeFluxBC(flux = 3e4) : HeatFluxBC(3e4)
    return PyrolysisProblem(
        mesh = create_stretched_mesh(0.01, n, Val(3), 1.15),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = top, bottom = DirichletThermalBC(500.0)),
            mass = (top = DirichletPressureBC(101325.0), bottom = ImpermeableBC()),
        ),
        T_initial = z -> 500.0 + 4e4 * z,
        ξ_initial = [z -> 1100.0 * (1 - 30z), z -> 300.0 * 30z, z -> 0.5],
        tspan = (0.0, 10.0),
    )
end

@testset "a policy names a known restart and flag rule" begin
    alg = OrdinaryDiffEq.Rodas5P()
    @test IntegratorPolicy(alg, :forced_dt, :keep, 1e-6).restart === :forced_dt
    @test IntegratorPolicy(alg, :heuristic, :clear, 1e-6).restart === :heuristic
    @test rodas5p_policy(KLUFactorization()).restart === :heuristic
    @test default_policy(KLUFactorization()).restart === :heuristic
    @test_throws ArgumentError IntegratorPolicy(alg, :immediate, :keep, 1e-6)
    @test_throws ArgumentError IntegratorPolicy(alg, :solver_chosen, :raise, 1e-6)
    @test_throws ArgumentError rodas5p_policy(KLUFactorization(); restart = :immediate)
end

@testset "domain check: admissibility, porosity and constraint residuals" begin
    cp = setup(
        _slab(6; κ = 1e-12, radiating = true);
        use_ale = true,
        radiation_model = BEER_LAMBERT,
    )
    layout = cp.def.layout
    check = DomainCheck(cp; max_consecutive = 3)
    m = check.monitor
    @test !check(cp.u0, nothing, 0.0)
    @test m.last <= 1e-3
    # Crossed nodes.
    u = copy(cp.u0)
    u[position_index(layout, 3)] = u[position_index(layout, 2)]
    @test check(u, nothing, 0.0)
    @test m.rejections == 1
    # A porosity below zero: more virgin solid than its intrinsic density fits.
    u = copy(cp.u0)
    V = cp.def.setup.reference.volume[2]
    u[mass_index(layout, 2, 1)] = 1.2 * 1400.0 * V
    @test check(u, nothing, 0.0)
    # A velocity off its balance row, and an intensity off its attenuation row.
    for row in (velocity_index(layout, 4), intensity_index(layout, 2))
        u = copy(cp.u0)
        u[row] += 1e3 * cp.abstol[row] / cp.reltol
        @test check(u, nothing, 0.0)
        @test m.last > 1
        @test !check(cp.u0, nothing, 0.0)
    end
    # After `max_consecutive` constraint rejections the state is accepted and
    # the monitor marked failed, for the observer to end the solve.
    u = copy(cp.u0)
    row = velocity_index(layout, 6)
    u[row] += 1e3 * cp.abstol[row] / cp.reltol
    @test all(_ -> check(u, nothing, 0.0), 1:3)
    @test !check(u, nothing, 0.0)
    @test m.failed
    @test occursin("consecutive", m.message)
end

# A 1 mm inert slab losing 500 kW/m² through its top: its sensible energy is
# gone at 0.6 s, and no physical state lies beyond.
function _cooling_slab()
    solid = SolidComponent(:solid; ρ = 1000.0, c = 1000.0, k = 1.0)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.001, 1, Val(1)),
        material = Material(name = :cooling, components = (solid,), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(-5e5), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [z -> 1000.0],
        tspan = (0.0, 1.0),
        dt_initial = 0.01,
        output_times = [0.0, 0.5, 1.0],
    )
end

@testset "a cell at or below absolute zero is inadmissible: $(energy) form" for energy in (
    :temperature,
    :enthalpy,
)
    cp = setup(_cooling_slab(); energy_state = energy)
    def, layout = cp.def, cp.def.layout
    m = (cp.u0[mass_index(layout, 1, 1)],)
    function at(T)
        u = copy(cp.u0)
        u[energy_index(layout, 1)] =
            energy_state_value(energy_state(layout), def.setup, T, m, 1)
        return u
    end
    @test below_absolute_zero(at(1.0), layout, def.setup) == 0
    @test below_absolute_zero(at(0.0), layout, def.setup) == 1
    @test inadmissible(at(-10.0), layout, def.setup, cp.abstol)
    # The raw-state test agrees with the inversion, which the trial
    # temperature keeps evaluable.
    du = zeros(length(cp.u0))
    residual!(du, at(-10.0), def, 0.0, cp.ws)
    @test cp.ws.derived.T[1] ≈ -10.0 rtol = 1e-10
    @test all(isfinite, du)
    check = DomainCheck(cp)
    @test check(at(-10.0), nothing, 0.0)
    @test check.monitor.cold == 1
    @test !check(cp.u0, nothing, 0.0)
    @test check.monitor.cold == 0
    @test_throws ArgumentError integrate(cp; u0 = at(-10.0))

    # Adaptive steps stop where the energy runs out, and say why.
    sol = @test_logs (:warn, r"absolute zero") match_mode = :any integrate(cp)
    @test sol.retcode !== :Success
    @test 0.55 < sol.t[end] <= 0.6 + 1e-9
    @test minimum(sol.T) > 0
    # A fixed step, which the domain test never sees, ends the solve at the
    # first accepted state below zero.
    sol = @test_logs (:warn, r"absolute zero") match_mode = :any integrate(
        cp;
        adaptive = false,
        dt = 0.07,
    )
    @test sol.retcode === :Failure
    @test sol.t[end] ≈ 0.63
end

# Stiff scalar relaxation towards 1: after the transient the controller grows
# the step by the largest factor it allows. The first right-hand side past
# `t_fail` returns NaN, which fails that stage's Newton once.
function _relaxation_problem(; t_fail = Inf)
    armed = Ref(isfinite(t_fail))
    f(u, p, t) =
        if armed[] && t > t_fail
            armed[] = false
            oftype(u, NaN)
        else
            -50.0 * (u - 1.0) * (1.0 + u^2)
        end
    return ODEProblem(f, 0.0, (0.0, 20.0))
end

_cap_callback(cap) = DiscreteCallback(
    (u, t, integ) -> true,
    integ -> (cap_step!(integ, cap); derivative_discontinuity!(integ, true));
    save_positions = (false, false),
)

@testset "Newton-failure step cap" begin
    g, M = FAILURE_MEMORY_GROWTH, FAILURE_MEMORY_STEPS
    integ = init(
        _relaxation_problem(t_fail = 0.2),
        KenCarp4();
        callback = _cap_callback(NewtonCap()),
        abstol = 1e-8,
        reltol = 1e-6,
    )
    rows = NamedTuple{(:dt, :proposed, :nlfail),Tuple{Float64,Float64,Int}}[]
    while integ.t < 20.0 && length(rows) < 400
        step!(integ)
        push!(
            rows,
            (
                dt = integ.t - integ.tprev,
                proposed = get_proposed_dt(integ),
                nlfail = integ.stats.nnonlinconvfail,
            ),
        )
    end
    i_fail = findfirst(r -> r.nlfail > 0, rows)
    @test i_fail !== nothing
    cap = g * rows[i_fail].dt
    for i = i_fail:(i_fail+M-1)
        @test rows[i].proposed <= cap * (1 + 1e-12)
    end
    @test any(i -> rows[i].proposed == cap, i_fail:(i_fail+M-1))
    @test rows[i_fail+M].proposed > cap

    # Without a failure the trajectory is the free one.
    kw = (; abstol = 1e-8, reltol = 1e-6)
    noop = DiscreteCallback(
        (u, t, integ) -> true,
        integ -> derivative_discontinuity!(integ, true);
        save_positions = (false, false),
    )
    free = OrdinaryDiffEq.solve(_relaxation_problem(), KenCarp4(); callback = noop, kw...)
    capped = OrdinaryDiffEq.solve(
        _relaxation_problem(),
        KenCarp4();
        callback = _cap_callback(NewtonCap()),
        kw...,
    )
    @test free.stats.nnonlinconvfail == 0
    @test capped.t == free.t
    @test capped.u == free.u
end

# A slab that converts at the same rate everywhere: isothermal, uniform, the
# activation energy zero, so its thickness follows the conversion exactly.
function _uniform_slab(n)
    a = SolidComponent(:a, ρ = 1000.0, c = 1000.0, k = 0.5)
    b = SolidComponent(:b, ρ = 400.0, c = 1000.0, k = 0.5)
    g = GaseousComponent(:g, M = 0.03, c = 1000.0, k = 0.03, λ = 0.0)
    mat = Material(
        name = :uniform,
        components = (a, b, g),
        reactions = (
            Reaction(:r, 1 => (2, 3); A = 0.5, E = 0.0, h = 0.0, yields = (0.5, 0.5)),
        ),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    return PyrolysisProblem(
        mesh = create_stretched_mesh(0.02, n, Val(3), 1.1),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(400.0), bottom = DirichletThermalBC(400.0)),
        ),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 1000.0, z -> 0.0, z -> 0.0],
        tspan = (0.0, 4.0),
    )
end

@testset "integrator policies: $name, $(energy) form" for (name, policy) in (
        ("KenCarp4", kencarp4_policy(KLUFactorization())),
        ("Rodas5P", rodas5p_policy(KLUFactorization())),
        (
            "Rodas5P, forced restart",
            rodas5p_policy(KLUFactorization(); restart = :forced_dt),
        ),
        (
            "Rodas5P, solver-chosen restart",
            rodas5p_policy(KLUFactorization(); restart = :solver_chosen),
        ),
        ("FBDF", fbdf_policy(KLUFactorization())),
    ),
    energy in (:temperature, :enthalpy)

    @test policy.restart in (:solver_chosen, :forced_dt, :heuristic)
    sol = Integration.solve(
        _uniform_slab(6);
        energy_state = energy,
        use_ale = true,
        abstol = 1e-10,
        reltol = 1e-9,
        saveat = [0.0, 1.0, 2.0, 4.0],
        integrator_policy = policy,
    )
    @test sol.retcode === :Success
    L0 = sol.extras.thickness_history[1]
    for (k, t) in enumerate(sol.t)
        X = 1 - exp(-0.5t)
        @test sol.extras.thickness_history[k] / L0 ≈ 1 + (0.5 * 1000 / 400 - 1) * X rtol =
            1e-6
    end
    @test sol.statistics.largest_constraint_residual <= 1
end

# Two decks that differ in component count, property types, reactions,
# boundary conditions, radiation model and gas transport.
function _plain_deck()
    mat = Material(
        name = :plain_deck,
        components = (
            SolidComponent(:resin, ρ = 1300.0, c = 1700.0, k = 0.20, ε = 0.9, λ = 2e-5),
            GaseousComponent(:vapor, M = 0.080, c = 1400.0, k = 0.03, λ = 1e-5),
        ),
        reactions = (Reaction(:cure, 1 => 2, A = 1e12, E = 1.6e5, h = 300e3),),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(5e-3, 8, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(40e3) + RadiativeBC(ε = 0.9, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        ),
        tspan = (0.0, 5.0),
        output_times = [0.0, 5.0],
    )
end

function _rich_deck()
    virgin = SolidComponent(
        :virgin,
        ρ = 500.0,
        c = CallableProperty(T -> 1500.0 + 0.5 * (T - 300.0)),
        k = CallableProperty(T -> 0.15 * (T / 300.0)^0.5),
        ε = 0.9,
        α = 5.0,
        κ = 1e-14,
        ρ_intrinsic = 1500.0,
        λ = 1e-5,
    )
    char = SolidComponent(
        :char,
        ρ = 150.0,
        c = 1100.0,
        k = 0.08,
        ε = 0.95,
        α = 5.0,
        κ = 1e-12,
        ρ_intrinsic = 1500.0,
        λ = 1e-5,
    )
    gas = GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
    mat = Material(
        name = :rich_deck,
        components = (virgin, char, gas),
        reactions = (
            Reaction(
                :pyrolysis,
                1 => (2, 3);
                A = 1e10,
                E = 150e3,
                h = 300e3,
                yields = (0.25, 0.75),
            ),
        ),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, 8, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 50e3) +
                      RadiativeBC(T_inf = 300.0) +
                      ConvectiveBC(h = 10.0, T_inf = 300.0),
                bottom = RadiativeBC(ε = 0.9, T_inf = 300.0),
            ),
            mass = (
                top = DiffusiveMassBC(0.05) + DirichletPressureBC(101325.0),
                bottom = ImpermeableBC(),
            ),
        ),
        tspan = (0.0, 5.0),
        output_times = [0.0, 5.0],
    )
end

# The integrator's type at the first step; every solve gets a probe of the
# same type.
function _solve_types(problem; kw...)
    seen = DataType[]
    probe = DiscreteCallback(
        (u, t, integrator) -> isempty(seen),
        integrator -> push!(seen, typeof(integrator));
        save_positions = (false, false),
    )
    sol = Integration.solve(problem; callback = probe, kw...)
    @test sol.retcode === :Success
    return only(seen), typeof(sol.raw_solution)
end

# On a fixed mesh in-depth radiation adds algebraic rows, and so a mass
# matrix: a different solver configuration, compared with itself.
@testset "the integrator's type does not depend on the problem: $name" for (
    name,
    plain,
    rich,
) in (
    ("fixed mesh", (;), (radiation_model = SURFACE_ABSORPTION,)),
    (
        "fixed mesh with in-depth radiation",
        (radiation_model = BEER_LAMBERT,),
        (radiation_model = BEER_LAMBERT,),
    ),
    ("moving mesh", (use_ale = true,), (use_ale = true, radiation_model = BEER_LAMBERT)),
    (
        "moving mesh with deactivation",
        (use_ale = true, deactivation = true),
        (use_ale = true, radiation_model = BEER_LAMBERT, deactivation = true),
    ),
)
    integrator_plain, solution_plain = _solve_types(_plain_deck(); plain...)
    integrator_rich, solution_rich = _solve_types(_rich_deck(); rich...)
    @test integrator_plain === integrator_rich
    @test solution_plain === solution_rich
end

@testset "a read-only user callback leaves every policy's trajectory unchanged" begin
    problem = _uniform_slab(12)
    probe = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> nothing;
        save_positions = (false, false),
    )
    for policy in (
        rodas5p_policy(KLUFactorization()),
        kencarp4_policy(KLUFactorization()),
        fbdf_policy(KLUFactorization()),
    )
        kwargs = (
            use_ale = true,
            integrator_policy = policy,
            saveat = Float64[],
            abstol = 1e-8,
            reltol = 1e-6,
        )
        plain = Integration.solve(problem; kwargs...)
        @test plain.retcode === :Success
        for callback in (probe, CallbackSet(probe))
            probed = Integration.solve(problem; kwargs..., callback)
            @test probed.raw_solution.t == plain.raw_solution.t
            @test probed.raw_solution.u == plain.raw_solution.u
            @test probed.statistics.accepted == plain.statistics.accepted
            @test probed.statistics.rhs_evaluations == plain.statistics.rhs_evaluations
        end
    end
end

@testset "repeated solves are bit-identical under each linear solver" begin
    problem = _uniform_slab(12)
    for linear_solver in (StructuredKLU(), KLUFactorization(), UMFPACKFactorization())
        sols = [
            Integration.solve(
                problem;
                use_ale = true,
                linear_solver,
                saveat = Float64[],
                abstol = 1e-8,
                reltol = 1e-6,
            ) for _ = 1:2
        ]
        @test sols[1].retcode === :Success
        @test sols[2].raw_solution.t == sols[1].raw_solution.t
        @test sols[2].raw_solution.u == sols[1].raw_solution.u
        @test sols[2].statistics.accepted == sols[1].statistics.accepted
    end
end

# An ablating solid heated at the top, its heat capacity rising with the
# temperature so that the caloric inversion is nonlinear: its cells thin and
# deactivate.
function _ablating(n)
    ρ, c, k = 1000.0, 1000.0, 1.0
    R, E, Tp, T0 = 8.314462618, 1.0e5, 700.0, 300.0
    v = 2.0e-3
    A = v^2 * E * (Tp - T0) * exp(E / (R * Tp)) / (k / (ρ * c) * R * Tp^2)
    mat = Material(
        name = :ablator,
        components = (
            SolidComponent(
                :solid,
                ρ = ρ,
                c = PolynomialProperty((700.0, 0.8, -1e-4)),
                k = k,
                λ = 1e-6,
            ),
            GaseousComponent(:gas, M = 0.05, c = c, k = 0.02, λ = 1e-4),
        ),
        reactions = (Reaction(:ablate, 1 => 2; A = A, E = E, h = 1.0e6),),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(3e-3, n, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(2.8e6), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0 + 1e5 * z,
        ξ_initial = [z -> ρ, z -> 0.0],
        tspan = (0.0, 0.8),
    )
end

@testset "every integration starts the caloric inversion from the set-up temperatures" begin
    problem = _ablating(30)
    kw = (energy_state = :enthalpy, use_ale = true, deactivation = true)
    cp = setup(problem; kw...)
    @test cp.ws.T_warm == cp.T0
    once = integrate(cp; saveat = Float64[])
    @test !isempty(once.events)
    # An evaluation outside the integrator, at enthalpies 1 % above the initial
    # ones, moves the inversion's starting point.
    layout = cp.def.layout
    v = copy(cp.u0)
    for i in eachindex(cp.T0)
        v[energy_index(layout, i)] *= 1.01
    end
    residual!(similar(v), v, cp.def, 0.0, cp.ws)
    @test cp.ws.T_warm != cp.T0
    again = integrate(cp; saveat = Float64[])
    fresh = integrate(setup(problem; kw...); saveat = Float64[])
    for sol in (again, fresh)
        @test sol.raw_solution.t == once.raw_solution.t
        @test sol.raw_solution.u == once.raw_solution.u
    end
end

@testset "the structured KLU factors the union of W's nonzeros" begin
    # A prototype with a full diagonal and two off-diagonal entries per column.
    n = 6
    rows = Int[]
    cols = Int[]
    for j = 1:n, i in (j - 1, j, j + 1)
        1 <= i <= n && (push!(rows, i); push!(cols, j))
    end
    prototype = sparse(rows, cols, ones(length(rows)), n, n)
    tridiagonal(lower, upper) = sparse(
        vcat(1:n, 2:n, 1:(n-1)),
        vcat(1:n, 1:(n-1), 2:n),
        vcat(fill(4.0, n), lower, upper),
        n,
        n,
    )
    b = collect(1.0:n)
    cache = LinearSolve.init(LinearSolve.LinearProblem(copy(prototype), b), StructuredKLU())
    c = cache.cacheval
    @test c isa Integration.StructuredKLUCache
    @test nnz(c.union) == n
    # The first W grows the union by its nonzeros; a stored zero stays out.
    W1 = tridiagonal(fill(-1.0, n - 1), vcat(0.0, fill(-1.0, n - 2)))
    cache.A = W1
    @test LinearSolve.solve!(cache).u ≈ Matrix(W1) \ b
    @test c.growths == 1
    @test nnz(c.union) == nnz(dropzeros(W1))
    # A W within the union refactors without growing it.
    W2 = tridiagonal(fill(-0.5, n - 1), vcat(0.0, fill(-2.0, n - 2)))
    cache.A = W2
    @test LinearSolve.solve!(cache).u ≈ Matrix(W2) \ b
    @test c.growths == 1
    # A nonzero outside the union grows it once more, to the prototype here.
    W3 = tridiagonal(fill(-1.0, n - 1), fill(-1.0, n - 1))
    cache.A = W3
    @test LinearSolve.solve!(cache).u ≈ Matrix(W3) \ b
    @test c.growths == 2
    @test nnz(c.union) == nnz(prototype)
    # An entry outside the prototype is an error.
    cache.A = sparse([1, n], [n, 1], [1.0, 1.0], n, n) + W3
    @test_throws ErrorException LinearSolve.solve!(cache)
end

# The KLU factorization in Rosenbrock's linear cache, a new object exactly
# when KLU analyses afresh.
_klu(cacheval) = cacheval
_klu(cacheval::Integration.StructuredKLUCache) = cacheval.inner.cacheval

@testset "one KLU analysis per growth of the union, use_ale = $ale" for ale in (true, false)
    problem = _slab(12)
    @test default_linear_solver() isa StructuredKLU
    function factored(linear_solver)
        factorizations = UInt[]
        last = Ref{Any}(nothing)
        probe = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator -> begin
                local cacheval = integrator.cache.linsolve.cacheval
                push!(factorizations, objectid(_klu(cacheval)))
                last[] = cacheval
            end;
            save_positions = (false, false),
        )
        sol = Integration.solve(
            problem;
            use_ale = ale,
            linear_solver,
            callback = probe,
            saveat = Float64[],
            abstol = 1e-8,
            reltol = 1e-6,
        )
        return sol, length(unique(factorizations)), last[]
    end
    structured, n_structured, c = factored(nothing)
    plain, _, _ = factored(KLUFactorization())
    @test structured.retcode === :Success
    @test c isa Integration.StructuredKLUCache
    @test 1 <= n_structured <= c.growths
    @test count(c.in_union) == nnz(c.union) < nnz(c.prototype)
    @test all(j -> j in rowvals(c.union)[nzrange(c.union, j)], 1:size(c.union, 2))
    @test structured.statistics.accepted == plain.statistics.accepted
    @test structured.statistics.rejected == plain.statistics.rejected
    @test maximum(abs, structured.T .- plain.T) <= 1e-9
end
