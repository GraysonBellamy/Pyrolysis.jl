using Test
using Pyrolysis
using Pyrolysis.Materials: DepletionLimiter
using Pyrolysis.Outputs
using Pyrolysis.Outputs: _extension
using Pyrolysis.Integration
using Pyrolysis.Integration: setup, integrate
using Pyrolysis.Assembly: residual!, buffers
using OrdinaryDiffEq: DiscreteCallback

# The char-layer column of V-CPL-06 with a reaction heat (the energy-ledger
# test's lumped wall deck): a reactive band low in a hot, gas-permeable matrix
# gasifies, and the gas diffuses up through the hot layer to leave through a
# film at the top, carrying its enthalpy out.
function _column(n)
    L = 0.01
    matrix = SolidComponent(
        :matrix,
        ρ = 2000.0,
        c = (800.0, 0.4),
        k = 1.0,
        λ = 1e-4,
        ρ_intrinsic = 2600.0,
    )
    fuel = SolidComponent(
        :fuel,
        ρ = 2000.0,
        c = 1000.0,
        k = 1.0,
        λ = 1e-4,
        ρ_intrinsic = 2600.0,
    )
    gas = GaseousComponent(:gas, M = 0.05, c = (900.0, 0.3), k = 0.02)
    rxn = Reaction(:gasify, 2 => 3; A = 9.3e3, E = 8.0e4, h = 3e5)
    mat = Material(name = :column, components = (matrix, fuel, gas), reactions = (rxn,))
    band(z) = 100.0 * exp(-((z / L - 0.3) / 0.08)^2)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = ConvectiveBC(h = 1000.0, T_inf = 900.0),
                bottom = ConvectiveBC(h = 1000.0, T_inf = 300.0),
            ),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 550.0 + 3e4 * z,
        ξ_initial = [z -> 400.0, band, z -> 0.05],
        tspan = (0.0, 20.0),
        output_times = collect(0.0:2.0:20.0),
    )
end

# The ablating solid of the deactivation tests: its cells deactivate at the
# heated top, each event venting the dead cell's gas.
function _ablator(n; L = 3e-3, t_end = 0.8)
    ρ, c, k = 1000.0, 1000.0, 1.0
    R, E, Tp, T0 = 8.314462618, 1.0e5, 700.0, 300.0
    v = 2.0e-3
    A = v^2 * E * (Tp - T0) * exp(E / (R * Tp)) / (k / (ρ * c) * R * Tp^2)
    solid = SolidComponent(:solid, φ = 0.05, ρ = ρ, c = c, k = k, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.05, c = c, k = 0.02)
    mat = Material(
        name = :ablator,
        components = (solid, gas),
        reactions = (Reaction(:ablate, 1 => 2; A = A, E = E, h = 1.0e6),),
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(2.8e6), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0 + 1e5 * z,
        ξ_initial = [z -> ρ, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = collect(0.0:0.05:t_end),
    )
end

const _TOL = (abstol = 1e-8, reltol = 1e-6)

_digest(sol) = hash((sol.T, sol.m, sol.z, sol.t, [e.t for e in sol.events], sol.statistics))

# The energy ledger gate's per-step closure: the trapezoid of the ledger rates
# over the accepted steps, from a read-only observer in the solve's own
# workspace, against the live inventories at the first and last steps.
function _per_step_trapezoid(problem; energy_state)
    cp = setup(problem; energy_state, _TOL...)
    def, ws = cp.def, cp.ws
    ext = _extension(def)
    du = zeros(length(cp.u0))
    acc, prev, last_t = zeros(2), zeros(2), Ref(0.0)
    residual!(du, cp.u0, def, 0.0, ws)
    ds, c = buffers(ws, def, Float64)
    prev .= ledger_rates(def, c, ds)
    first_node = collect(live_inventories(def, ext, ds, cp.u0))
    last_node = copy(first_node)
    step = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> begin
            t = integrator.t
            residual!(du, integrator.u, def, t, ws)
            ds, c = buffers(ws, def, Float64)
            cur = collect(ledger_rates(def, c, ds))
            acc .+= 0.5 .* (prev .+ cur) .* (t - last_t[])
            prev .= cur
            last_t[] = t
            last_node .= live_inventories(def, ext, ds, integrator.u)
            return nothing
        end;
        save_positions = (false, false),
    )
    integrate(cp; callback = step)
    ΔE, ΔM = last_node .- first_node
    return (; energy = ΔE - acc[1], mass = ΔM + acc[2], inflow = acc[1], ΔE)
end

@testset "the recorder on the gasifying column: $(energy) form" for energy in (
    :temperature,
    :enthalpy,
)
    problem = _column(20)
    rec = ConservationRecorder()
    sol = Integration.solve(problem; energy_state = energy, callback = rec, _TOL...)
    @test sol.retcode === :Success
    # Read-only: the trajectory of a solve without it.
    @test _digest(sol) ==
          _digest(Integration.solve(problem; energy_state = energy, _TOL...))
    @test rec.steps == sol.statistics.accepted
    audit = conservation_audit(sol, rec)
    @test audit.method === :recorder
    @test audit.mass.lost == 0 && audit.energy.lost == 0
    # The gas leaving at the top carries enthalpy out: the terms are not small.
    @test audit.mass.inflow < 0
    trapezoid = _per_step_trapezoid(problem; energy_state = energy)
    @info "column ($energy): recorder vs per-step trapezoid" audit.energy audit.mass trapezoid
    @test audit.energy.change ≈ trapezoid.ΔE rtol = 1e-12
    @test abs(audit.energy.inflow - trapezoid.inflow) <= 2 * abs(trapezoid.energy)
    @test abs(audit.energy.relative) <= 1e-6
    # The mass changes by a thousandth of the inventory: the defect is at the
    # integrator's tolerance on the inventory.
    @test abs(audit.mass.defect) <= 1e-8 * audit.mass.initial

    # Without a recorder, the trapezoid over the saved states of rates the
    # test evaluates itself.
    cp = setup(problem; energy_state = energy, _TOL...)
    def, ws = cp.def, cp.ws
    du = zeros(length(cp.u0))
    rates = map(zip(sol.t, sol.raw_solution.u)) do (t, u)
        residual!(du, u, def, t, ws)
        collect(
            ledger_rates(def, buffers(ws, def, Float64)[2], buffers(ws, def, Float64)[1]),
        )
    end
    inflow = sum(
        (rates[k-1][1] + rates[k][1]) * (sol.t[k] - sol.t[k-1]) / 2 for k = 2:length(sol.t)
    )
    outflow = sum(
        (rates[k-1][2] + rates[k][2]) * (sol.t[k] - sol.t[k-1]) / 2 for k = 2:length(sol.t)
    )
    saved = conservation_audit(sol)
    @test saved.method === :trapezoid
    @test saved.energy.inflow ≈ inflow rtol = 1e-12
    @test saved.mass.inflow ≈ -outflow rtol = 1e-12
    @test saved.energy.change ≈ audit.energy.change rtol = 1e-12
    # The saving interval's error is the trapezoid's, far above the recorder's.
    @test abs(saved.energy.relative) > abs(audit.energy.relative)
end

# On a fixed mesh an event vents the dead cell's gas and its enthalpy; on a
# moving one the receiver keeps them.
@testset "the recorder across deactivation events: $(moving ? "moving" : "fixed") mesh, $(energy) form" for moving in
                                                                                                            (
        true,
        false,
    ),
    energy in (:temperature, :enthalpy)

    problem = _ablator(12)
    options = (; energy_state = energy, use_ale = moving, deactivation = true)
    audits = map(((1e-6, 1e-8), (1e-8, 1e-10))) do (reltol, abstol)
        rec = ConservationRecorder()
        sol = Integration.solve(
            problem;
            callback = (rec,),
            options...,
            reltol,
            abstol,
            pressure_tolerance = 1e6*reltol,
        )
        @test sol.retcode === :Success
        @test length(sol.events) >= 3
        return conservation_audit(sol, rec), sol, rec
    end
    (coarse, _, _), (fine, sol, rec) = audits
    @test rec.steps == sol.statistics.accepted
    @test _digest(sol) == _digest(
        Integration.solve(
            problem;
            options...,
            reltol = 1e-8,
            abstol = 1e-10,
            pressure_tolerance = 0.01,
        ),
    )
    @info "ablator ($(moving ? "moving" : "fixed"), $energy): recorder across events" coarse.mass coarse.energy fine.mass fine.energy
    # The events' terms are the events' records.
    @test fine.mass.lost == sum(sum(e.lost_mass) for e in sol.events)
    @test fine.energy.lost == sum(e.lost_enthalpy + e.regularization for e in sol.events)
    # The defect is the integrator's error: it falls with the tolerance.
    @test abs(fine.mass.defect) < abs(coarse.mass.defect) / 10
    @test abs(fine.energy.defect) < abs(coarse.energy.defect) / 10
    @test abs(fine.mass.relative) <= 1e-7
    @test abs(fine.energy.relative) <= 1e-6
    if moving
        @test fine.mass.lost == 0
    else
        # The vented gas and its enthalpy close the ledger; without them it
        # would be off by them.
        @test fine.mass.lost > 0
        @test abs(fine.mass.defect) < 1e-3 * fine.mass.lost
        @test abs(fine.energy.defect) < 1e-2 * abs(fine.energy.lost)
    end
    # The recorder ends where the solve ends, after the last event.
    @test fine.mass.final ≈ sol.extras.live_mass[end] rtol = 1e-12
end

# A decomposing slab whose char is permeable and whose cross section shrinks
# with the conversion: the reaction raises the pore pressure, and the gas
# vents through a pressure condition at the heated top while the face area
# falls.
function _shrinking_vent(n)
    virgin = SolidComponent(
        :virgin,
        ρ = 1100.0,
        c = (900.0, 1.5),
        k = (0.2, 2e-4),
        λ = 3e-6,
        ρ_intrinsic = 1400.0,
    )
    char = SolidComponent(
        :char,
        ρ = 300.0,
        c = 700.0,
        k = 0.1,
        λ = 3e-5,
        κ = 1e-11,
        ρ_intrinsic = 1800.0,
    )
    gas = GaseousComponent(:gas, M = 0.03, c = 1100.0, k = 0.03)
    mat = Material(
        name = :shrinking_vent,
        mixing = MaterialMixing(κ = PARALLEL),
        components = (virgin, char, gas),
        reactions = (
            Reaction(:pyro, 1 => (2, 3); A = 5e9, E = 1.4e5, h = 8e5, yields = (0.3, 0.7)),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1.0),
        lateral_shrinkage_law = χ -> 1 - 0.3χ,
    )
    return PyrolysisProblem(
        mesh = create_uniform_mesh(0.01, n),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(3e4) + ConvectiveBC(h = 15.0, T_inf = 300.0),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(
                    GasReservoir(P = 101325.0, T = 300.0, mole_fractions = (gas = 1.0,)),
                ),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 500.0 + 3e4 * z,
        ξ_initial = [z -> 1100.0 * (1 - 30z), z -> 300.0 * 30z, z -> 0.3],
        tspan = (0.0, 10.0),
        output_times = collect(0.0:1.0:10.0),
    )
end

@testset "the recorder on a venting column whose cross section shrinks: $(energy) form" for energy in
                                                                                            (
    :temperature,
    :enthalpy,
)
    rec = ConservationRecorder()
    sol = Integration.solve(
        _shrinking_vent(12);
        energy_state = energy,
        use_ale = true,
        callback = rec,
        reltol = 1e-8,
        abstol = 1e-10,
    )
    @test sol.retcode === :Success
    @test rec.steps == sol.statistics.accepted
    A = sol.extras.cross_section_area_history
    @test A[end] < 0.99 * A[1]
    audit = conservation_audit(sol, rec)
    @test audit.mass.inflow < 0
    @test abs(audit.mass.relative) <= 1e-7
    @test abs(audit.energy.relative) <= 1e-6
    # The outflow is integrated at the face's current area: taken at the
    # initial one it would be off by a share of the area's fall, far above
    # the defect.
    @test abs(audit.mass.inflow) * (1 - A[end] / A[1]) > 1e3 * abs(audit.mass.defect)
end
