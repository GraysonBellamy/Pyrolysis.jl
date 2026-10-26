using Test
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT, SURFACE_ABSORPTION
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Closures
using Pyrolysis.Assembly: residual!
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration
using Pyrolysis.Integration: setup, integrate
using Pyrolysis.Outputs: ConservationRecorder, conservation_audit

# ============================================================================
# Decks
# ============================================================================

# V-BO-01's configuration restricted so that every rate scales with `s`: the
# incident flux, the convective coefficient, the conductivities, the
# gas-transfer coefficients and the pre-exponential factor times `s`, no
# re-radiation (the emitted flux does not scale), an adiabatic and
# impermeable back, and the times divided by `s`. With `s` a power of two
# every scaled input is exact.
function _restricted(s; n_top = 8, n_bottom = 4, t_end = 200.0)
    L_top, L_bottom = 2e-3, 3e-3
    top = SolidComponent(
        :top,
        φ = 0.3,
        ρ = 1200.0,
        c = 1500.0,
        k = s * 0.2,
        ε = 0.95,
        α = 2.0,
        λ = s * 1e-6,
    )
    sub = SolidComponent(
        :sub,
        φ = 0.8,
        ρ = 250.0,
        c = 1100.0,
        k = s * 0.06,
        ε = 1.0,
        α = 30.0,
        λ = s * 1e-6,
    )
    gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = s * 0.02)
    mat = Material(
        name = :restricted,
        components = (top, sub, gas),
        reactions = (Reaction(:gasify, 1 => 3; A = s * 1e11, E = 1.6e5, h = 8e5),),
    )
    z = vcat(
        range(0.0, L_bottom, length = n_bottom + 1),
        range(L_bottom, L_bottom + L_top, length = n_top + 1)[2:end],
    )
    return PyrolysisProblem(
        mesh = create_mesh_from_nodes(z),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = s * 5e4, absorptivity = 0.95) +
                      ConvectiveBC(h = s * 10.0, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [
            z -> z > L_bottom ? 1200.0 : 0.0,
            z -> z > L_bottom ? 0.0 : 250.0,
            z -> 0.0,
        ],
        tspan = (0.0, t_end / s),
        dt_initial = 1e-3 / s,
        output_times = collect(range(0.0, t_end / s, length = 101)),
    )
end

# Two layers that gasify completely, a fast one over a slow one, on an inert
# substrate of `n_sub` cells.
function _two_vanishing(; n_fast = 6, n_slow = 6, n_sub = 4, t_end = 400.0)
    fast = SolidComponent(
        :fast,
        φ = 0.3,
        ρ = 1200.0,
        c = 1500.0,
        k = 0.2,
        ε = 0.95,
        α = 2.0,
        λ = 1e-6,
    )
    slow = SolidComponent(
        :slow,
        φ = 0.3,
        ρ = 1000.0,
        c = 1500.0,
        k = 0.2,
        ε = 0.95,
        α = 2.0,
        λ = 1e-6,
    )
    sub = SolidComponent(
        :sub,
        φ = 0.8,
        ρ = 250.0,
        c = 1100.0,
        k = 0.06,
        ε = 1.0,
        α = 30.0,
        λ = 1e-6,
    )
    gas = GaseousComponent(:gas, M = 0.1, c = 1500.0, k = 0.02)
    mat = Material(
        name = :two_vanishing,
        components = (fast, slow, sub, gas),
        reactions = (
            Reaction(:fast, 1 => 4; A = 1e11, E = 1.6e5, h = 8e5),
            Reaction(:slow, 2 => 4; A = 2e10, E = 1.6e5, h = 8e5),
        ),
    )
    z1, z2, z3 = 2e-3, 3e-3, 4e-3
    z = vcat(
        range(0.0, z1, length = n_sub + 1),
        range(z1, z2, length = n_slow + 1)[2:end],
        range(z2, z3, length = n_fast + 1)[2:end],
    )
    return PyrolysisProblem(
        mesh = create_mesh_from_nodes(z),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = 5e4, absorptivity = 0.95) +
                      RadiativeBC(T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 300.0,
        ξ_initial = [
            z -> z > z2 ? 1200.0 : 0.0,
            z -> z1 < z <= z2 ? 1000.0 : 0.0,
            z -> z <= z1 ? 250.0 : 0.0,
            z -> 0.0,
        ],
        tspan = (0.0, t_end),
        output_times = collect(range(0.0, t_end, length = 101)),
    )
end

const _burnout = (use_ale = true, deactivation = true)

_set_up(problem; τ = ActivationTime(), radiation_model = BEER_LAMBERT, kw...) = setup(
    problem;
    _burnout...,
    radiation_model,
    correction = EquidistributionCorrection(; τ),
    kw...,
)

# A solve of `problem` and its bound correction afterwards.
function _solved(problem; kw...)
    cp = _set_up(problem; kw...)
    return integrate(cp; min_thickness = 1e-4), cp.def.closure.correction
end

# The largest differences between two solutions of the restricted deck in
# dimensionless time, `a` at scale `sa` and `b` at `sb`: the saved
# temperatures, node positions and event times, the dead cells' masks equal.
function _differences(a, sa, b, sb)
    @test length(a.t) == length(b.t)
    @test a.t .* sa ≈ b.t .* sb
    @test isnan.(a.T) == isnan.(b.T)
    @test length(a.events) == length(b.events)
    return (
        T = maximum(abs, filter(isfinite, a.T .- b.T)),
        z = maximum(abs, filter(isfinite, a.z .- b.z)),
        events = maximum(
            abs,
            [e.t * sa for e in a.events] .- [e.t * sb for e in b.events];
            init = 0.0,
        ),
    )
end

# The state of `cp` with the cell `i` at `ratio` of its initial thickness,
# the nodes above it shifted with it, and the node velocities `w`.
function _state(cp, i, ratio, w)
    layout = cp.def.layout
    h0 = cp.def.setup.reference.thickness
    u = copy(cp.u0)
    shift = (ratio - 1) * h0[i]
    for nd = (i+1):n_nodes(layout)
        u[position_index(layout, nd)] += shift
    end
    for nd = 1:n_nodes(layout)
        u[velocity_index(layout, nd)] = w[nd]
    end
    return u
end

# ============================================================================
# The rescaled deck
# ============================================================================

@testset "the rescaled deck equals the unscaled one in dimensionless time" begin
    # The transport equation scales with s because κ=0 and both Λ and every
    # thermal/kinetic rate scale with s. Resolve vacuum-start pressure and
    # caloric transients before comparing this time-change invariant.
    controls = (reltol = 1e-8, abstol = 1e-10, pressure_tolerance = 0.001)
    scaled_solve(p; kw...) = _solved(p; controls..., kw...)
    # The solve's own error at its tolerance bounds the rescaled pairs.
    self = _differences(
        first(scaled_solve(_restricted(1.0); τ = 720.0)),
        1.0,
        first(scaled_solve(_restricted(1.0); τ = 720.0, reltol = 1e-9, abstol = 1e-11)),
        1.0,
    )
    ref, _ = scaled_solve(_restricted(1.0); τ = 720.0)
    ref_b, corr_b = scaled_solve(_restricted(1.0))
    @test only(corr_b.stage) == 2
    for s in (2.0, 0.5)
        # The core: the layer time scaled by hand.
        sol, _ = scaled_solve(_restricted(s); τ = 720.0 / s)
        d = _differences(ref, 1.0, sol, s)
        @test d.T <= 5 * self.T
        @test d.T < 0.01 # K; independent physical tolerance
        @test d.z <= max(5 * self.z, 1e-12) # m; roundoff floor
        @test d.events <= max(5 * self.events, 1e-6) # s
        # Candidate B scales by itself.
        sol, corr = scaled_solve(_restricted(s))
        d = _differences(ref_b, 1.0, sol, s)
        @test d.T <= 5 * self.T
        @test d.T < 0.01 # K; independent physical tolerance
        @test d.z <= max(5 * self.z, 1e-12) # m; roundoff floor
        @test d.events <= max(5 * self.events, 1e-6) # s
        @test only(corr.time) * s ≈ only(corr_b.time) rtol = 1e-6
        @test only(corr.set_at) * s ≈ only(corr_b.set_at) rtol = 1e-6
        # One fixed time does not.
        sol, _ = scaled_solve(_restricted(s); τ = 720.0)
        d = _differences(ref, 1.0, sol, s)
        @test d.T > 20 * self.T
        @test d.events > 20 * self.events
    end
end

# ============================================================================
# The estimate
# ============================================================================

@testset "binding: the moving layers, the fallback before the estimate" begin
    cp = _set_up(_restricted(1.0))
    corr = cp.def.closure.correction
    @test corr.layers == [(5, 13)]
    @test corr.stage == [0]
    @test corr.time == [EQUIDISTRIBUTION_TIME]
    @test all(isnan, corr.set_at)
    @test corr.gain[6:12] == fill(64 / EQUIDISTRIBUTION_TIME, 7)
    @test all(iszero, corr.gain[[1:5; 13]])
    @test estimates_layer_times(corr)
    fixed = _set_up(_restricted(1.0); τ = 100.0).def.closure.correction
    @test !estimates_layer_times(fixed)
    @test fixed.time == [100.0] && fixed.stage == [2] && fixed.set_at == [0.0]
    @test fixed.gain[6:12] == fill(64 / 100.0, 7)
    @test_throws ArgumentError ActivationTime(λ = 0.0)
    @test_throws ArgumentError ActivationTime(fallback = -1.0)
    @test_throws ArgumentError ActivationTime(floor = -1.0)
end

@testset "τ = λ h/|ḣ| at the first activation; each guard is dimensionless" begin
    cp = _set_up(_restricted(1.0))
    corr = cp.def.closure.correction
    layout, ref = cp.def.layout, cp.def.setup.reference
    nn = n_nodes(layout)
    h0 = ref.thickness
    margin(u) = layer_time_margin(corr, layout, ref, u)
    estimate!(u, t; horizon = 200.0) =
        estimate_layer_times!(corr, layout, ref, u, t, horizon; slack = 1e-6)
    # Floor: the layer's 2 mm over the 200 s horizon.
    floor = 2e-3 / 200
    w(ḣ) = [nd == nn ? ḣ : 0.0 for nd = 1:nn]

    @test margin(cp.u0) ≈ 1 / 0.3 - 1
    above = _state(cp, 12, 0.31, w(-3floor))
    @test margin(above) > 0
    estimate!(above, 10.0)
    @test corr.stage == [0] && corr.time == [EQUIDISTRIBUTION_TIME]

    u = _state(cp, 12, 0.29, w(-3floor))
    @test margin(u) ≈ 0.29 / 0.3 - 1
    estimate!(u, 10.0)
    h = u[position_index(layout, nn)] - u[position_index(layout, 5)]
    @test corr.stage == [2]
    @test only(corr.time) == ACTIVATION_TIME_FACTOR * h / (3floor)
    @test corr.set_at == [10.0]
    @test corr.gain[6:12] == fill(64 / only(corr.time), 7)
    @test all(iszero, corr.gain[[1:5; 13]])
    # Estimated once: later crossings change nothing.
    @test margin(_state(cp, 12, 0.05, w(-floor))) == Inf
    estimate!(_state(cp, 12, 0.05, w(-10floor)), 20.0)
    @test corr.set_at == [10.0]

    # The same state with every rate doubled over half the horizon: half the time.
    τ1 = only(corr.time)
    reset_layer_times!(corr)
    @test corr.stage == [0] && all(isnan, corr.set_at)
    @test corr.gain[6:12] == fill(64 / EQUIDISTRIBUTION_TIME, 7)
    estimate!(_state(cp, 12, 0.29, w(-6floor)), 5.0; horizon = 100.0)
    @test only(corr.time) == τ1 / 2

    # A rate at or below the floor waits for the second threshold, where the
    # estimate is taken once more; the floor scales with the horizon.
    for (s, ḣ, horizon) in ((1.0, -0.5floor, 200.0), (2.0, -floor, 100.0))
        reset_layer_times!(corr)
        estimate!(_state(cp, 12, 0.29, w(ḣ)), 10.0 / s; horizon)
        @test corr.stage == [1]
        @test corr.time == [EQUIDISTRIBUTION_TIME] && all(isnan, corr.set_at)
        u = _state(cp, 12, 0.2, w(ḣ))
        @test margin(u) ≈ 0.2 / 0.1 - 1
        estimate!(u, 11.0 / s; horizon)
        @test corr.stage == [1]
        u = _state(cp, 12, 0.09, w(4ḣ))
        @test margin(u) ≈ 0.09 / 0.1 - 1
        estimate!(u, 12.0 / s; horizon)
        h = u[position_index(layout, nn)] - u[position_index(layout, 5)]
        @test corr.stage == [2]
        @test only(corr.time) ≈ ACTIVATION_TIME_FACTOR * h / (4 * abs(ḣ))
        @test corr.set_at == [12.0 / s]
    end
    # Still floored at the second threshold: the fallback, once.
    reset_layer_times!(corr)
    estimate!(_state(cp, 12, 0.29, w(-0.5floor)), 10.0)
    estimate!(_state(cp, 12, 0.09, w(-0.5floor)), 12.0)
    @test corr.stage == [2] && corr.time == [EQUIDISTRIBUTION_TIME]
    @test corr.set_at == [12.0]
    # A thickening layer keeps the fallback.
    reset_layer_times!(corr)
    estimate!(_state(cp, 12, 0.29, w(floor)), 10.0)
    @test corr.stage == [2] && corr.time == [EQUIDISTRIBUTION_TIME]
    @test corr.set_at == [10.0]
    # A dead cell's thickness is not read: the margin is over the live cells.
    reset_layer_times!(corr)
    u = _state(cp, 12, 0.05, w(-floor))
    layout.mask[energy_index(layout, 12)] = false
    try
        @test margin(u) ≈ 1 / 0.3 - 1
    finally
        fill!(layout.mask, true)
    end
end

# ============================================================================
# Solves
# ============================================================================

@testset "two layers: each takes its time at its own activation" begin
    for radiation_model in (BEER_LAMBERT, SURFACE_ABSORPTION)
        sol, corr = _solved(_two_vanishing(); radiation_model)
        @test sol.retcode == :Success
        @test corr.layers == [(5, 11), (11, 17)]
        @test corr.stage == [2, 2]
        @test all(>(0), corr.set_at)
        @test corr.set_at[1] != corr.set_at[2]
        @test corr.time[1] != corr.time[2]
        @test all(<(EQUIDISTRIBUTION_TIME / 10), corr.time)
        for (k, (lo, hi)) in enumerate(corr.layers)
            @test corr.gain[(lo+1):(hi-1)] == fill((hi - lo)^2 / corr.time[k], hi - lo - 1)
        end
        # The substrate cannot vanish and keeps its mesh with the solid.
        @test all(iszero, corr.gain[1:5])
        @test all(iszero, sol.correction_velocity[1:5, :])
        # One invocation per layer.
        @test sol.statistics.invocations == 2
    end
end

@testset "containment with two layers of different τ: $(energy) form" for energy in (
    :temperature,
    :enthalpy,
)
    cp = _set_up(
        _two_vanishing();
        energy_state = energy,
        correction = EquidistributionCorrection(
            τ = ActivationTime(),
            activation = (1.5, 0.5),
            monitor = (1.0, 1.0),
        ),
    )
    corr = cp.def.closure.correction
    layout = cp.def.layout
    for (k, τ) in enumerate((7.0, 300.0))
        lo, hi = corr.layers[k]
        corr.gain[(lo+1):(hi-1)] .= (hi - lo)^2 / τ
    end
    u0 = copy(cp.u0)
    for nd = 2:(n_nodes(layout)-1)
        u0[position_index(layout, nd)] *= 1 + 0.05 * sin(nd)
    end
    L = length(u0)
    for u in (u0, u0 .* (1 .+ 0.2 .* sin.(1:L)), -abs.(u0))
        du = zeros(L)
        residual!(du, u, cp.def, 1.0, cp.ws)
        @test all(isfinite, du)
        Jd = dense_jacobian(cp.def, cp.ws, u, 1.0)
        @test isempty(missing_entries(Jd, cp.jacobian.pattern))
        J = jacobian!(jac_prototype(cp.jacobian), cp.jacobian, u, 1.0)
        @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
    end
    reset_layer_times!(corr)
end

@testset "the event: its time against the tolerance, admissibility, repeated solves" begin
    problem = _restricted(1.0)
    runs = [_solved(problem; reltol) for reltol in (1e-6, 1e-8)]
    (a, ca), (b, cb) = runs
    @test a.retcode == b.retcode == :Success
    @test abs(only(ca.set_at) - only(cb.set_at)) <= 1e-5 * only(cb.set_at)
    @test only(ca.time) ≈ only(cb.time) rtol = 1e-5
    # The state at the event is the trajectory's: the solve runs on through it
    # and deactivates the layer in one invocation.
    @test a.statistics.invocations == 1
    @test only(ca.set_at) < first(a.events).t
    # Every call starts before the first activation.
    cp = _set_up(problem)
    first_sol = integrate(cp; min_thickness = 1e-4, saveat = Float64[])
    τ_first = copy(cp.def.closure.correction.time)
    again = integrate(cp; min_thickness = 1e-4, saveat = Float64[])
    @test cp.def.closure.correction.time == τ_first
    @test again.raw_solution.t == first_sol.raw_solution.t
    @test again.raw_solution.u == first_sol.raw_solution.u
end

@testset "the conservation recorder integrates up to the event" begin
    rec = ConservationRecorder()
    cp = _set_up(_restricted(1.0; t_end = 70.0))
    sol = integrate(cp; min_thickness = 1e-4, callback = rec)
    corr = cp.def.closure.correction
    @test only(corr.stage) == 2 && only(corr.set_at) < 70.0
    audit = conservation_audit(sol, rec)
    @test abs(audit.mass.relative) <= 1e-6
    @test abs(audit.energy.relative) <= 1e-6
end
