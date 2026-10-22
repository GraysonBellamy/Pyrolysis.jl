using Test
using Pyrolysis
using Pyrolysis.Outputs
using Pyrolysis.Integration
using Pyrolysis.Integration: rodas5p_policy, fbdf_policy, kencarp4_policy, default_policy
using Pyrolysis.Integration: default_linear_solver
using OrdinaryDiffEq: DiscreteCallback, terminate!

# An ablating solid, as in the deactivation tests: it decomposes to a gas that
# takes no volume, so on a moving mesh a cell's thickness falls with its
# condensed mass and the cells at the heated top deactivate one by one.
function _ablator(n; L = 3e-3, t_end = 0.8, from_bottom = false, output_times = nothing)
    ρ, c, k = 1000.0, 1000.0, 1.0
    R, E, Tp, T0 = 8.314462618, 1.0e5, 700.0, 300.0
    v = 2.0e-3
    A = v^2 * E * (Tp - T0) * exp(E / (R * Tp)) / (k / (ρ * c) * R * Tp^2)
    solid = SolidComponent(:solid, ρ = ρ, c = c, k = k, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.05, c = c, k = 0.02, λ = 1e-4)
    mat = Material(
        name = :ablator,
        components = (solid, gas),
        reactions = (Reaction(:ablate, 1 => 2; A = A, E = E, h = 1.0e6),),
    )
    heated = HeatFluxBC(2.8e6)
    thermal =
        from_bottom ? (top = AdiabaticBC(), bottom = heated) :
        (top = heated, bottom = AdiabaticBC())
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = thermal,
            mass = (
                top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = from_bottom ? (z -> T0 + 1e5 * (L - z)) : (z -> T0 + 1e5 * z),
        ξ_initial = [z -> ρ, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = something(output_times, collect(0.0:0.05:t_end)),
    )
end

const _TOL = (abstol = 1e-7, reltol = 1e-5)

@testset "accessors on a moving deactivating deck" begin
    problem = _ablator(12)
    sol = Integration.solve(problem; use_ale = true, deactivation = true, _TOL...)
    @test sol.retcode === :Success
    @test length(sol.events) >= 3
    nt = length(sol.t)
    mid = nt ÷ 2
    @test get_times(sol) === sol.t
    @test get_temperatures(sol) === sol.T
    @test get_concentrations(sol) === sol.ξ
    @test get_gas_generation_rate(sol) === sol.extras.gas_generation_rate
    @test n_timesteps(sol) == nt

    # Cell positions at a saved time, the initial ones at the first, NaN for
    # the deactivated cells, which are those beyond the outermost live node.
    z0 = [p[1] for p in problem.mesh.positions]
    @test get_cell_positions(sol; at = 1) == [(z0[i] + z0[i+1]) / 2 for i = 1:12]
    for k in (1, mid, nt)
        x = get_cell_positions(sol; at = k)
        live = sol.extras.live_cells[k]
        @test all(isfinite, x[1:live]) && all(isnan, x[(live+1):end])
        z = sol.z[:, k]
        @test x[1:live] == [(z[i] + z[i+1]) / 2 for i = 1:live]
    end
    @test isequal(get_cell_positions(sol), get_cell_positions(sol; at = nt))
    @test sol.extras.live_cells[end] < sol.extras.live_cells[1]

    # Profiles: every field at a saved time, deactivated cells NaN.
    live = sol.extras.live_cells[nt]
    @test isequal(get_profile(sol, :T, nt), sol.T[:, nt])
    @test isequal(get_profile(sol, :ξ, nt), sol.ξ[:, :, nt])
    @test isequal(get_profile(sol, :ξ2, nt), sol.ξ[2, :, nt])
    m = get_profile(sol, :m, nt)
    @test m[:, 1:live] == sol.m[:, 1:live, nt]
    @test all(isnan, m[:, (live+1):end])
    @test isequal(get_profile(sol, :m1, nt), m[1, :])
    z = get_profile(sol, :z, nt)
    @test z[1:(live+1)] == sol.z[1:(live+1), nt] && all(isnan, z[(live+2):end])
    @test_throws ArgumentError get_profile(sol, :P, nt)
    fixed = Integration.solve(_ablator(6; t_end = 0.1); _TOL...)
    @test get_profile(fixed, :z, 1) == [p[1] for p in fixed.problem.mesh.positions]

    # The fields at a time: at a saved time the saved ones; between two saved
    # states no event separates, the reconstruction of the integrator's
    # interpolant, whose masses are the states; refused across an event.
    for k in (1, mid, nt)
        f = sol(sol.t[k])
        @test f.t == sol.t[k]
        @test isequal(f.T, sol.T[:, k])
        @test isequal(f.ξ, sol.ξ[:, :, k])
        @test f.m == sol.m[:, :, k]
        @test isequal(f.z, sol.z[:, k])
        @test isequal(f.cell_centers, sol.cell_centers[:, k])
    end
    separated(a, b) = any(e -> a < e.t <= b, sol.events)
    quiet = findfirst(
        k -> !separated(sol.t[k], sol.t[k+1]) && sol.extras.live_cells[k] < 12,
        1:(nt-1),
    )
    @test quiet !== nothing
    a, b = sol.t[quiet], sol.t[quiet+1]
    t = (2a + b) / 3
    f = interpolate_solution(sol, t)
    θ = (t - a) / (b - a)
    live = sol.extras.live_cells[quiet]
    @test f.m[:, 1:live] ≈
          (1 - θ) .* sol.m[:, 1:live, quiet] .+ θ .* sol.m[:, 1:live, quiet+1] rtol = 1e-12
    @test all(isnan, f.T[(live+1):end]) && all(isfinite, f.T[1:live])
    @test f.z[1:(live+1)] ≈
          (1 - θ) .* sol.z[1:(live+1), quiet] .+ θ .* sol.z[1:(live+1), quiet+1] rtol =
        1e-12
    @test all(isnan, f.z[(live+2):end])
    crossed = findfirst(k -> separated(sol.t[k], sol.t[k+1]), 1:(nt-1))
    @test crossed !== nothing
    t = (sol.t[crossed] + sol.t[crossed+1]) / 2
    err = try
        sol(t)
        ""
    catch e
        sprint(showerror, e)
    end
    @test occursin("deactivation of cell", err)
    @test_throws ArgumentError sol(sol.t[end] + 1.0)
    # The solution's mask is the last event's again after a query.
    @test isequal(sol(sol.t[end]).T, sol.T[:, end])

    text = sprint(show, MIME"text/plain"(), sol)
    @test occursin("Outcome: running", text)
    @test occursin("Invocations: $(sol.statistics.invocations)", text)
    @test occursin("outcome :running", sprint(show, sol))
end

@testset "succeeded and completed by return code and outcome" begin
    options = (; deactivation = true, _TOL...)
    ok = Integration.solve(_ablator(6; t_end = 0.1); options...)
    @test ok.retcode === :Success && succeeded(ok) && completed(ok)

    burned = Integration.solve(
        _ablator(4; L = 1e-3, t_end = 20.0, output_times = [0.0, 20.0]);
        use_ale = true,
        options...,
    )
    @test burned.retcode === :Terminated && burned.outcome === :burned_out
    @test succeeded(burned) && !completed(burned)

    thin = Integration.solve(
        _ablator(6; t_end = 5.0);
        use_ale = true,
        min_thickness = 2e-3,
        options...,
    )
    @test thin.retcode === :Terminated && thin.outcome === :min_thickness
    @test succeeded(thin) && !completed(thin)
    @test thin.extras.thickness_history[end] ≈ 2e-3 rtol = 1e-3
    # Cells went before the stop, so the thickness it ends on is the live
    # span, up to the outermost live node.
    live = thin.extras.live_cells[end]
    @test !isempty(thin.events) && live < 6
    @test thin.extras.thickness_history[end] == thin.z[live+1, end] - thin.z[1, end]

    void = Integration.solve(
        _ablator(6; from_bottom = true, t_end = 5.0);
        use_ale = true,
        options...,
    )
    @test void.retcode === :Terminated && void.outcome === :interior
    @test !succeeded(void) && !completed(void)

    stop = DiscreteCallback((u, t, integrator) -> t > 0.05, terminate!)
    user = Integration.solve(_ablator(6; t_end = 0.1); callback = stop, _TOL...)
    @test user.retcode === :Terminated && user.outcome === :running
    @test !succeeded(user) && !completed(user)

    failed = Integration.solve(_ablator(6; t_end = 0.1); maxiters = 3, _TOL...)
    @test failed.retcode !== :Success && failed.retcode !== :Terminated
    @test !succeeded(failed) && !completed(failed)
end

_digest(sol) = hash((sol.T, sol.m, sol.z, sol.t, [e.t for e in sol.events], sol.statistics))

@testset "the progress bar is inert: $(name)" for (name, policy) in (
    ("default", default_policy(default_linear_solver())),
    ("Rodas5P", rodas5p_policy(default_linear_solver())),
    ("FBDF", fbdf_policy(default_linear_solver())),
    ("KenCarp4", kencarp4_policy(default_linear_solver())),
)
    for (problem, options) in (
        (_ablator(8; t_end = 0.1), (;)),
        (_ablator(8; t_end = 0.4), (use_ale = true, deactivation = true)),
    )
        quiet = Integration.solve(problem; integrator_policy = policy, options..., _TOL...)
        drawn = redirect_stderr(devnull) do
            Integration.solve(
                problem;
                integrator_policy = policy,
                progress = true,
                options...,
                _TOL...,
            )
        end
        @test _digest(drawn) == _digest(quiet)
        @test drawn.retcode === quiet.retcode
    end
end
