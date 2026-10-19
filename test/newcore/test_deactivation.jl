using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: ReactivitySwitch, reactivity, SWITCH_WIDTH
using Pyrolysis.Assembly
using Pyrolysis.Assembly: residual!, buffers
using Pyrolysis.Deactivation
using Pyrolysis.Deactivation: ConservativeTransfer
using Pyrolysis.ColoredJacobian
using Pyrolysis.Integration
using Pyrolysis.Integration: setup, consistent_state!

# An ablating solid: it decomposes to a gas that takes no volume, so a cell's
# thickness falls with its condensed mass, as in V-CPL-02. `radiating` makes
# the solid absorb in depth.
function _ablator(n; L = 3e-3, q = 2.8e6, t_end = 0.8, radiating = false)
    ρ, c, k = 1000.0, 1000.0, 1.0
    R, E, Tp, T0 = 8.314462618, 1.0e5, 700.0, 300.0
    v = 2.0e-3
    A = v^2 * E * (Tp - T0) * exp(E / (R * Tp)) / (k / (ρ * c) * R * Tp^2)
    solid =
        SolidComponent(:solid, ρ = ρ, c = c, k = k, λ = 1e-6, α = radiating ? 50.0 : 0.0)
    gas = GaseousComponent(:gas, M = 0.05, c = c, k = 0.02, λ = 1e-4)
    mat = Material(
        name = :ablator,
        components = (solid, gas),
        reactions = (Reaction(:ablate, 1 => 2; A = A, E = E, h = 1.0e6),),
    )
    top = radiating ? RadiativeFluxBC(flux = q) : HeatFluxBC(q)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n, Val(2)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = top, bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T0 + 1e5 * z,
        ξ_initial = [z -> ρ, z -> 0.0],
        tspan = (0.0, t_end),
        output_times = collect(0.0:0.05:t_end),
    )
end

_total_masses(u, layout) = [
    sum(u[mass_index(layout, i, j)] for i = 1:n_cells(layout)) for
    j = 1:n_components(layout)
]

# Σ m h(T) over the cells, in either energy form.
function _total_enthalpy(cp, u)
    layout, s = cp.def.layout, cp.def.setup
    energy_state(layout) isa EnthalpyState &&
        return sum(u[energy_index(layout, i)] for i = 1:n_cells(layout))
    ds = cp.ws.derived
    derive!(ds, u, layout, s)
    return sum(
        ds.m[j, i] * sensible_enthalpy(s.laws.h[j], ds.T[i]) for
        i = 1:n_cells(layout), j = 1:n_components(layout)
    )
end

# The state with cell `i`'s solid at the fraction `f` of its initial mass and
# its products in the gas.
function _deplete!(u, cp, i, f)
    layout = cp.def.layout
    m0 = cp.deactivation.reference[i]
    u[mass_index(layout, i, 2)] += (1 - f) * m0
    u[mass_index(layout, i, 1)] = f * m0
    return u
end

@testset "a transfer conserves the inventories and the receiver absorbs the volume: $(energy) form, $(moving ? "moving" : "fixed") mesh" for energy in
                                                                                                                                             (
        :temperature,
        :enthalpy,
    ),
    moving in (true, false)

    cp = setup(_ablator(6); energy_state = energy, use_ale = moving, deactivation = true)
    def, layout = cp.def, cp.def.layout
    cells = cp.deactivation
    u = _deplete!(copy(cp.u0), cp, 6, 0.5 * cells.threshold)
    consistent_state!(u, def, cp.ws, 0.0)
    masses, H = _total_masses(u, layout), _total_enthalpy(cp, u)
    x_top = moving ? u[position_index(layout, 7)] : NaN
    @test any_below(cells, layout, u)
    @test deactivate!(u, cells, def, 0.3) === :transferred
    @test top_live_cell(layout) == 5
    @test _total_masses(u, layout) ≈ masses rtol = 1e-14
    @test _total_enthalpy(cp, u) ≈ H rtol = 1e-12
    @test all(u[mass_index(layout, 6, j)] == 0 for j = 1:2)
    @test !any(layout.mask[cell_rows(layout, 6)])
    @test all(layout.mask[cell_rows(layout, 5)])
    e = only(cells.events)
    @test (e.t, e.cell, e.receiver) == (0.3, 6, 5)
    @test e.created_volume == 0
    if moving
        # The receiver's far face is the dead cell's far face; both of the dead
        # cell's nodes are there.
        @test u[position_index(layout, 6)] ≈ x_top rtol = 1e-14
        @test u[position_index(layout, 7)] == u[position_index(layout, 6)]
        @test !any(layout.mask[node_rows(layout, 7)])
    end
    consistent_state!(u, def, cp.ws, 0.3)
    du = zeros(length(u))
    residual!(du, u, def, 0.3, cp.ws)
    @test all(isfinite, du)
    # The thermal condition of :top now acts on the receiver's far face.
    @test live_face(def.setup.topology, cp.ws.derived.active, 7) == 6
    Jd = dense_jacobian(def, cp.ws, u, 0.3)
    @test isempty(missing_entries(Jd, cp.jacobian.pattern))
    J = jacobian!(jac_prototype(cp.jacobian), cp.jacobian, u, 0.3)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
end

@testset "every cell below the threshold goes in one pass, from the exposed side" begin
    cp = setup(_ablator(8); use_ale = true, deactivation = true)
    def, layout, cells = cp.def, cp.def.layout, cp.deactivation
    u = copy(cp.u0)
    for i = 6:8
        _deplete!(u, cp, i, 0.2 * cells.threshold)
    end
    masses = _total_masses(u, layout)
    @test deactivate!(u, cells, def, 1.0) === :transferred
    @test [e.cell for e in cells.events] == [8, 7, 6]
    @test [e.receiver for e in cells.events] == [7, 6, 5]
    @test top_live_cell(layout) == 5
    @test _total_masses(u, layout) ≈ masses rtol = 1e-14
    # The receiver keeps its own reference: with 0.6 % of three cells' solid
    # added it stays live.
    @test !below_threshold(cells, layout, u, 5)

    # A cell below the threshold inward of a live one is not handled.
    reset_cells!(cells, layout)
    u = _deplete!(copy(cp.u0), cp, 4, 0.2 * cells.threshold)
    @test deactivate!(u, cells, def, 1.0) === :interior
    @test cells.outcome[] === :interior
    @test isempty(cells.events)

    # The last cell goes with no receiver: its inventories are lost.
    one_cell = setup(_ablator(1); use_ale = true, deactivation = true)
    u = _deplete!(copy(one_cell.u0), one_cell, 1, 0.5 * one_cell.deactivation.threshold)
    m = [u[mass_index(one_cell.def.layout, 1, j)] for j = 1:2]
    @test deactivate!(u, one_cell.deactivation, one_cell.def, 2.0) === :burned_out
    e = only(one_cell.deactivation.events)
    @test e.receiver == 0
    @test e.lost_mass ≈ m
    @test e.lost_enthalpy > 0
end

# Two unit cells of a planar column (unit area) and the pore-preserving merge
# of design §3.5: component A's intrinsic density falls with temperature
# (ρ_A = 300/T), B's is one, and the heat capacities are equal (large enough
# that the cells sit far above the heat-capacity floor).
function _merge_pair(state_a, state_b)
    c = 1e6
    a = SolidComponent(:A, ρ = 1.0, c = c, k = 1.0, ρ_intrinsic = T -> 300 / T)
    b = SolidComponent(:B, ρ = 1.0, c = c, k = 1.0, ρ_intrinsic = 1.0, κ = 1e-12)
    g = GaseousComponent(:G, M = 0.028, c = c, k = 0.03, λ = 1e-5)
    mat = Material(name = :pair, components = (a, b, g), reactions = ())
    problem = PyrolysisProblem(
        mesh = create_uniform_mesh(2.0, 2, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
        ),
        T_initial = z -> 400.0,
        ξ_initial = [z -> 0.5, z -> 0.4, z -> 1e-4],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem; use_ale = true, deactivation = true)
    layout = cp.def.layout
    u = copy(cp.u0)
    for (i, s) in ((1, state_a), (2, state_b))
        for j = 1:3
            u[mass_index(layout, i, j)] = s.m[j]
        end
        u[energy_index(layout, i)] = s.T
    end
    Pyrolysis.Deactivation._transfer!(
        u,
        cp.deactivation,
        ConservativeTransfer(),
        cp.def,
        2,
        1,
        0.0,
    )
    ds = cp.ws.derived
    derive!(ds, u, layout, cp.def.setup)
    V = u[position_index(layout, 2)] - u[position_index(layout, 1)]
    return (;
        m = [u[mass_index(layout, 1, j)] for j = 1:3],
        T = ds.T[1],
        V,
        P = ds.P[1],
        φ = ds.φ[1],
        created = only(cp.deactivation.events).created_volume,
    )
end

@testset "a merge never destroys pore space" begin
    @testset "the pore-free counterexample" begin
        r = _merge_pair((m = [1.0, 0.0, 0.0], T = 300.0), (m = [0.0, 1.0, 0.0], T = 500.0))
        @test r.T ≈ 400 atol = 1e-10
        # Summed volumes 2, merged skeletal volume 4/3 + 1: the created volume is
        # the skeletal expansion at the merged temperature.
        @test r.V ≈ 7 / 3 atol = 1e-12
        @test r.created ≈ 1 / 3 atol = 1e-12
        @test r.m == [1.0, 1.0, 0.0]
    end
    @testset "gas-containing merges keep their pore volume at a finite pressure" begin
        for φ in (0.01, 0.05, 0.1), gas in (1e-6, 1e-3)
            r = _merge_pair(
                (m = [1 - φ, 0.0, gas], T = 300.0),
                (m = [0.0, 1 - φ, gas], T = 500.0),
            )
            @test r.T ≈ 400 atol = 1e-10
            @test r.V ≈ 7 / 3 - φ / 3 atol = 1e-12
            @test r.created ≈ 1 / 3 - φ / 3 atol = 1e-12
            @test r.φ ≈ 2φ / (7 / 3 - φ / 3) rtol = 1e-10
            @test r.P ≈ (2gas / 0.028) * 8.314462618 * 400 / (2φ) rtol = 1e-8
        end
    end
    @testset "a pore-free donor and a porous receiver" begin
        r = _merge_pair(
            (m = [1.0, 0.0, 0.0], T = 300.0),
            (m = [0.0, 0.95, 1e-3], T = 500.0),
        )
        @test r.T ≈ 300 + 0.951 * 200 / 1.951 atol = 1e-10
        @test r.V - (1 / (300 / r.T) + 0.95) ≈ 0.05 atol = 1e-12
        @test isfinite(r.P) && r.P > 0
    end
end

@testset "masked rows are identity rows and the live block is the compacted problem's: $(name), $(energy) form" for name in
                                                                                                                    (
        :fixed,
        :moving,
        :beer_lambert,
    ),
    energy in (:temperature, :enthalpy)

    options =
        name === :fixed ? (;) :
        name === :moving ? (use_ale = true,) :
        (use_ale = true, radiation_model = BEER_LAMBERT)
    radiating = name === :beer_lambert
    full = setup(
        _ablator(6; radiating);
        energy_state = energy,
        deactivation = true,
        options...,
    )
    compact = setup(
        _ablator(5; L = 2.5e-3, radiating);
        energy_state = energy,
        deactivation = true,
        options...,
    )
    layout = full.def.layout
    Lc = length(compact.u0)
    u = zeros(length(full.u0))
    u[1:Lc] .= compact.u0
    if is_moving(layout)
        u[position_index(layout, 7)] = u[position_index(layout, 6)]
    end
    dead = vcat(collect(cell_rows(layout, 6)), collect(node_rows(layout, 7)))
    layout.mask[dead] .= false
    for (state_name, v, vc) in (
        ("set up", u, compact.u0),
        ("perturbed", u .* (1 .+ 0.1 .* sin.(1:length(u))), nothing),
    )
        @testset "$state_name" begin
            vc = vc === nothing ? v[1:Lc] : vc
            Jf = dense_jacobian(full.def, full.ws, v, 0.2)
            Jc = dense_jacobian(compact.def, compact.ws, vc, 0.2)
            @test compare_entries(Jf[1:Lc, 1:Lc], Jc; rtol = 1e-11).passed
            @test all(iszero, Jf[1:Lc, (Lc+1):end])
            for r in dead
                row = Jf[r, :]
                if layout.differential[r]
                    @test all(iszero, row)
                else
                    @test row[r] == 1 && count(!iszero, row) == 1
                end
            end
            J = jacobian!(jac_prototype(full.jacobian), full.jacobian, v, 0.2)
            @test compare_entries(Matrix(J), Jf; rtol = 1e-12).passed
        end
    end
    fill!(layout.mask, true)
end

@testset "the reactivity switch" begin
    cp = setup(_ablator(4); use_ale = true, deactivation = true)
    cells = cp.deactivation
    sw = ReactivitySwitch(cells.threshold, cells.reference, cells.solid)
    ds = buffers(cp.ws, cp.def, Float64)[1]
    derive!(ds, cp.u0, cp.def.layout, cp.def.setup)
    m0 = cells.reference[2]
    s(f) = (ds.m[1, 2] = f * m0; reactivity(sw, ds, 2))
    edge = SWITCH_WIDTH * cells.threshold
    @test s(1.0) === 1.0
    @test s(edge) === 1.0
    @test s(cells.threshold) ≈ 20 / 27 rtol = 1e-12
    @test s(0.0) == 0
    @test s(-1.0) == 0
    # C¹: zero slope at both ends of the band, finite everywhere.
    D = ForwardDiff.Dual{Nothing}
    slope(f) = (
        ds2 = buffers(cp.ws, cp.def, D{Float64,1})[1];
        ds2.m[1, 2] = D(f * m0, m0);
        ds2.m[2, 2] = D(0.0, 0.0);
        ForwardDiff.partials(reactivity(sw, ds2, 2), 1)
    )
    @test abs(slope(edge * (1 - 1e-9))) < 1e-6
    @test abs(slope(1e-12)) < 1e-6
    @test all(isfinite, slope.(range(-0.1, 2edge; length = 50)))
    # Without deactivation the kernel has no switch; above the band the rates
    # are those of the plain kernel.
    plain = setup(_ablator(4); use_ale = true)
    du_s, du_p = zeros(length(cp.u0)), zeros(length(plain.u0))
    residual!(du_s, cp.u0, cp.def, 0.0, cp.ws)
    residual!(du_p, plain.u0, plain.def, 0.0, plain.ws)
    @test du_s == du_p
end

@testset "a solve through deactivations: $(energy) form, $(moving ? "moving" : "fixed") mesh" for energy in
                                                                                                  (
        :temperature,
        :enthalpy,
    ),
    moving in (true, false)

    sol = Integration.solve(
        _ablator(30);
        energy_state = energy,
        use_ale = moving,
        deactivation = true,
        abstol = 1e-7,
        reltol = 1e-5,
    )
    @test sol.retcode === :Success
    @test sol.outcome === :running
    cells = [e.cell for e in sol.events]
    @test length(cells) >= 5
    @test cells == 30:-1:(31-length(cells))
    @test all(e.receiver == e.cell - 1 for e in sol.events)
    @test sol.statistics.largest_constraint_residual <= 1
    last_live = 30 - length(cells)
    @test all(isnan, sol.T[(last_live+1):end, end])
    @test all(isfinite, sol.T[1:last_live, end])
    if moving
        z = sol.z[:, end]
        @test all(isnan, z[(last_live+2):end])
        @test sol.extras.thickness_history[end] ≈ z[last_live+1] - z[1]
        @test all(diff(sol.extras.thickness_history) .<= 1e-12 * 3e-3)
    end
    # The total mass only falls, by the gas leaving through the top face.
    m = sol.extras.total_mass
    @test all(diff(m) .<= 1e-12 * m[1])
    @test all(>(0), sol.extras.gas_generation_rate[2:end])
    # One park time per deactivated cell, and the live cells at every saved
    # time (a state saved at an event's time is the one after the event).
    @test sol.extras.park_times == [e.t for e in sol.events]
    @test sol.extras.live_cells == [30 - count(e -> e.t <= t, sol.events) for t in sol.t]
    @test sol.extras.live_cells[end] == last_live
end

@testset "a terminated solve saves its last state, and a state saved at an event's time is the one after it" begin
    for moving in (true, false)
        problem = _ablator(4; L = 1e-3, t_end = 20.0)
        sol = Integration.solve(
            problem;
            use_ale = moving,
            deactivation = true,
            saveat = [0.0, 20.0],
            abstol = 1e-7,
            reltol = 1e-5,
        )
        @test sol.retcode === :Terminated
        @test sol.outcome === :burned_out
        @test length(sol.events) == 4
        @test sol.t[end] == sol.extras.park_times[end]
        @test length(sol.t) == 2
        @test all(isnan, sol.T[:, end])
        @test sol.extras.live_cells[end] == 0
        # A state saved exactly at an event's time carries the post-event mask.
        t_event = sol.events[1].t
        sol2 = Integration.solve(
            problem;
            use_ale = moving,
            deactivation = true,
            saveat = [0.0, t_event, 20.0],
            abstol = 1e-7,
            reltol = 1e-5,
        )
        k = findfirst(==(t_event), sol2.t)
        @test k !== nothing
        @test sol2.extras.live_cells[k] == 3
        @test isnan(sol2.T[4, k])
        @test all(isfinite, sol2.T[1:3, k])
    end
end
