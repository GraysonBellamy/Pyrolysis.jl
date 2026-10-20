using Test
using ForwardDiff
using SparseArrays
using Pyrolysis
using Pyrolysis.Materials: BEER_LAMBERT, T_REF
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
using Pyrolysis.Integration: setup, consistent_state!, integrate
using Pyrolysis.Integration: fbdf_policy, rodas5p_policy, kencarp4_policy
using OrdinaryDiffEq: DiscreteCallback
using LinearSolve: KLUFactorization, UMFPACKFactorization

# An ablating solid: it decomposes to a gas that takes no volume, so a cell's
# thickness falls with its condensed mass, as in V-CPL-02. `radiating` makes
# the solid absorb in depth.
function _ablator(
    n;
    L = 3e-3,
    q = 2.8e6,
    t_end = 0.8,
    radiating = false,
    from_bottom = false,
)
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
    heated = radiating ? RadiativeFluxBC(flux = q) : HeatFluxBC(q)
    thermal =
        from_bottom ? (top = AdiabaticBC(), bottom = heated) :
        (top = heated, bottom = AdiabaticBC())
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, n, Val(2)),
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
    gas = u[mass_index(layout, 6, 2)]
    @test deactivate!(u, cells, def, 0.3) === :transferred
    @test top_live_cell(layout) == 5
    e = only(cells.events)
    # A fixed mesh vents the dead cell's gas to the event.
    @test e.lost_mass == (moving ? [0.0, 0.0] : [0.0, gas])
    @test _total_masses(u, layout) .+ e.lost_mass ≈ masses rtol = 1e-14
    @test _total_enthalpy(cp, u) + e.lost_enthalpy ≈ H rtol = 1e-12
    @test (e.lost_enthalpy > 0) == !moving
    @test abs(e.regularization) <= 1e-12 * abs(H)
    @test all(u[mass_index(layout, 6, j)] == 0 for j = 1:2)
    @test !any(layout.mask[cell_rows(layout, 6)])
    @test all(layout.mask[cell_rows(layout, 5)])
    @test (e.t, e.invocation, e.cell, e.receiver) == (0.3, 1, 6, 5)
    @test e.fraction ≈ 0.5 * cells.threshold rtol = 1e-14
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
    if moving
        ds, c = cp.ws.derived, cp.ws.contributions
        scale = sum(abs, c.θ .* ds.V_eval) + sum(abs, ds.geometry.face_area .* ds.w)
        @test all(abs(du[velocity_index(layout, nd)]) <= 1e-10 * scale for nd = 1:7)
    end
    # The thermal condition of :top now acts on the receiver's far face.
    @test live_face(def.setup.topology, cp.ws.derived.active, 7) == 6
    Jd = dense_jacobian(def, cp.ws, u, 0.3)
    @test isempty(missing_entries(Jd, cp.jacobian.pattern))
    J = jacobian!(jac_prototype(cp.jacobian), cp.jacobian, u, 0.3)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
end

# The extended energy of the live cells, in either energy form.
function _live_energy(cp, u)
    layout, s, ext = cp.def.layout, cp.def.setup, cp.deactivation.caloric
    ds = cp.ws.derived
    derive!(ds, u, layout, s)
    return sum(1:top_live_cell(layout)) do i
        m = ntuple(j -> u[mass_index(layout, i, j)], n_components(layout))
        energy_state(layout) isa EnthalpyState ? u[energy_index(layout, i)] :
        extended_energy(ext, s.laws, ds.T[i], m, i)
    end
end

# The live extended energy plus what the events wrote as lost or regularized.
_ledger(cp, u, events) =
    _live_energy(cp, u) +
    sum(e.lost_enthalpy + e.regularization for e in events; init = 0.0)

@testset "a run below the threshold merges into the cell beneath it in one invocation: $(energy) form" for energy in
                                                                                                           (
    :temperature,
    :enthalpy,
)
    cp = setup(_ablator(8); energy_state = energy, use_ale = true, deactivation = true)
    def, layout, cells = cp.def, cp.def.layout, cp.deactivation
    u = copy(cp.u0)
    for i = 6:8
        _deplete!(u, cp, i, 0.2 * cells.threshold)
    end
    consistent_state!(u, def, cp.ws, 0.0)
    masses, E = _total_masses(u, layout), _live_energy(cp, u)
    x_top = u[position_index(layout, 9)]
    @test deactivate!(u, cells, def, 1.0) === :transferred
    @test [e.cell for e in cells.events] == [8, 7, 6]
    @test all(e.receiver == 5 && e.invocation == 1 && e.t == 1.0 for e in cells.events)
    @test all(e.fraction ≈ 0.2 * cells.threshold for e in cells.events)
    @test top_live_cell(layout) == 5
    @test _total_masses(u, layout) ≈ masses rtol = 1e-14
    @test _ledger(cp, u, cells.events) ≈ E rtol = 1e-12
    # The merge's created volume and regularization are written to its first
    # event; the receiver's far face is the run's far face, where every node
    # of the run now sits.
    @test all(e.created_volume == 0 && e.regularization == 0 for e in cells.events[2:end])
    V0 = def.setup.topology.cross_section * 3e-3
    @test abs(cells.events[1].created_volume) <= 1e-12 * V0
    @test u[position_index(layout, 6)] ≈ x_top rtol = 1e-12
    @test all(u[position_index(layout, nd)] == u[position_index(layout, 6)] for nd = 7:9)
    @test all(u[velocity_index(layout, nd)] == 0 for nd = 7:9)
    @test !any(
        layout.mask[vcat(collect(cell_rows(layout, 6)), collect(node_rows(layout, 7)))],
    )
    # The receiver keeps its own reference: with 0.6 % of three cells' solid
    # added it stays live, and the state is admissible and consistent.
    @test !below_threshold(cells, layout, u, 5)
    consistent_state!(u, def, cp.ws, 1.0)
    @test !inadmissible(u, layout, cp.abstol)
    du = zeros(length(u))
    residual!(du, u, def, 1.0, cp.ws)
    @test all(isfinite, du)
    Jd = dense_jacobian(def, cp.ws, u, 1.0)
    @test isempty(missing_entries(Jd, cp.jacobian.pattern))

    # A second pass is the solve's next invocation.
    _deplete!(u, cp, 5, 0.5 * cells.threshold)
    @test deactivate!(u, cells, def, 2.0) === :transferred
    @test last(cells.events).invocation == 2
    @test (last(cells.events).cell, last(cells.events).receiver) == (5, 4)

    # A cell below the threshold inward of a live one is not handled.
    reset_cells!(cells, layout)
    u = _deplete!(copy(cp.u0), cp, 4, 0.2 * cells.threshold)
    @test deactivate!(u, cells, def, 1.0) === :interior
    @test cells.outcome[] === :interior
    @test isempty(cells.events)
end

@testset "a column below the band leaves its innermost cell as the survivor: $(energy) form" for energy in
                                                                                                 (
    :temperature,
    :enthalpy,
)
    cp = setup(_ablator(8); energy_state = energy, use_ale = true, deactivation = true)
    def, layout, cells = cp.def, cp.def.layout, cp.deactivation
    u = copy(cp.u0)
    for i = 1:8
        _deplete!(u, cp, i, (0.4 + 0.1 * isodd(i)) * cells.threshold)
    end
    consistent_state!(u, def, cp.ws, 0.0)
    masses, E = _total_masses(u, layout), _live_energy(cp, u)
    x_top = u[position_index(layout, 9)]
    @test deactivate!(u, cells, def, 1.0) === :transferred
    @test [e.cell for e in cells.events] == 8:-1:2
    @test all(e.receiver == 1 && e.invocation == 1 for e in cells.events)
    @test top_live_cell(layout) == 1
    # The survivor holds the column's inventories: 3.6 of its own threshold.
    @test condensed_mass(cells, layout, u, 1) ≈ 3.6 * cells.threshold * cells.reference[1] rtol =
        1e-12
    @test [u[mass_index(layout, 1, j)] for j = 1:2] ≈ masses rtol = 1e-14
    @test _ledger(cp, u, cells.events) ≈ E rtol = 1e-12
    @test u[position_index(layout, 2)] ≈ x_top rtol = 1e-12

    # A survivor still below the threshold is the last cell, and goes.
    reset_cells!(cells, layout)
    u = copy(cp.u0)
    for i = 1:8
        _deplete!(u, cp, i, 0.1 * cells.threshold)
    end
    masses = _total_masses(u, layout)
    @test deactivate!(u, cells, def, 1.0) === :burned_out
    @test cells.outcome[] === :burned_out
    @test [e.cell => e.receiver for e in cells.events] == [((8:-1:2) .=> 1)..., 1 => 0]
    @test all(e.invocation == 1 for e in cells.events)
    @test last(cells.events).lost_mass ≈ masses rtol = 1e-14
    @test top_live_cell(layout) == 0

    # The last cell goes with no receiver: its inventories are lost.
    one_cell = setup(_ablator(1); use_ale = true, deactivation = true)
    u = _deplete!(copy(one_cell.u0), one_cell, 1, 0.5 * one_cell.deactivation.threshold)
    m = [u[mass_index(one_cell.def.layout, 1, j)] for j = 1:2]
    @test deactivate!(u, one_cell.deactivation, one_cell.def, 2.0) === :burned_out
    e = only(one_cell.deactivation.events)
    @test (e.receiver, e.invocation) == (0, 1)
    @test e.lost_mass ≈ m
    @test e.lost_enthalpy > 0
end

@testset "the near-threshold band" begin
    # Absolute condensed fractions, outermost first: below ε_m, inside the
    # band twice, above it.
    fractions(ε_m) = (8 => 0.9ε_m, 7 => 0.015, 6 => 0.05, 5 => 0.2)
    for ε_m in (1e-2, 3e-3)
        cp = setup(
            _ablator(8);
            use_ale = true,
            deactivation = true,
            threshold = ε_m,
            band = 0.1,
        )
        u = copy(cp.u0)
        for (i, f) in fractions(ε_m)
            _deplete!(u, cp, i, f)
        end
        @test deactivate!(u, cp.deactivation, cp.def, 1.0) === :transferred
        events = cp.deactivation.events
        @test [e.cell for e in events] == [8, 7, 6]
        @test all(e.receiver == 5 for e in events)
        @test [e.fraction for e in events] ≈ [0.9ε_m, 0.015, 0.05] rtol = 1e-12
    end
    # At its default the band is the threshold: the run is the cells below it.
    cp = setup(_ablator(8); use_ale = true, deactivation = true)
    @test cp.deactivation.band == cp.deactivation.threshold
    u = copy(cp.u0)
    for (i, f) in fractions(1e-2)
        _deplete!(u, cp, i, f)
    end
    @test deactivate!(u, cp.deactivation, cp.def, 1.0) === :transferred
    @test [e.cell => e.receiver for e in cp.deactivation.events] == [8 => 7]
    # A run in the band without a cell below the threshold stays.
    cp = setup(_ablator(8); use_ale = true, deactivation = true, band = 0.1)
    u = copy(cp.u0)
    for i = 6:8
        _deplete!(u, cp, i, 0.05)
    end
    @test !any_below(cp.deactivation, cp.def.layout, u)
    @test deactivate!(u, cp.deactivation, cp.def, 1.0) === :transferred
    @test isempty(cp.deactivation.events)
    # A fixed mesh keeps the outermost-cell rule: the receiver rises above the
    # threshold and the pass stops.
    cp = setup(_ablator(8); deactivation = true, band = 0.1)
    u = copy(cp.u0)
    for (i, f) in fractions(1e-2)
        _deplete!(u, cp, i, f)
    end
    @test deactivate!(u, cp.deactivation, cp.def, 1.0) === :transferred
    @test [e.cell => e.receiver for e in cp.deactivation.events] == [8 => 7]
    @test top_live_cell(cp.def.layout) == 7
    # The band lies in [threshold, 1).
    @test_throws ArgumentError setup(_ablator(4); deactivation = true, band = 5e-3)
    @test_throws ArgumentError setup(_ablator(4); deactivation = true, band = 1.0)
    @test setup(_ablator(4); band = 5e-3).deactivation === nothing
end

@testset "a transfer writes the regularization the receiver does not hold: $(energy) form, $(moving ? "moving" : "fixed") mesh" for energy in
                                                                                                                                    (
        :temperature,
        :enthalpy,
    ),
    moving in (true, false)

    cp = setup(_ablator(6); energy_state = energy, use_ale = moving, deactivation = true)
    def, layout, cells = cp.def, cp.def.layout, cp.deactivation
    # The dead cell's solid at a tenth of the threshold and no gas: its
    # guaranteed capacity is the floor, half of where the switch begins.
    u = _deplete!(copy(cp.u0), cp, 6, 0.1 * cells.threshold)
    u[mass_index(layout, 6, 2)] = 0
    consistent_state!(u, def, cp.ws, 0.0)
    m6 = ntuple(j -> u[mass_index(layout, 6, j)], 2)
    κ = regularization(cells.caloric, m6, 6)
    @test κ > 0
    E = _live_energy(cp, u)
    T6 = cp.ws.derived.T[6]
    @test deactivate!(u, cells, def, 0.3) === :transferred
    e = only(cells.events)
    # The receiver is far above the floor before and after: the event holds
    # the dead cell's regularization.
    @test e.regularization ≈ κ * (T6 - T_REF) rtol = 1e-9
    @test _live_energy(cp, u) + e.lost_enthalpy + e.regularization ≈ E rtol = 1e-12
end

# Two unit cells of a planar column (unit area) and the pore-preserving merge
# of design §3.5, or on a fixed mesh the transfer: component A's intrinsic
# density falls with temperature (ρ_A = 300/T), B's is one, and the heat
# capacities are equal (large enough that the cells sit far above the
# heat-capacity floor).
function _merge_pair(state_a, state_b; use_ale = true)
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
    cp = setup(problem; use_ale, deactivation = true)
    layout = cp.def.layout
    u = copy(cp.u0)
    for (i, s) in ((1, state_a), (2, state_b))
        for j = 1:3
            u[mass_index(layout, i, j)] = s.m[j]
        end
        u[energy_index(layout, i)] = s.T
    end
    Pyrolysis.Deactivation._merge!(
        u,
        cp.deactivation,
        ConservativeTransfer(),
        cp.def,
        1,
        2:2,
        0.0,
        1,
    )
    ds = cp.ws.derived
    derive!(ds, u, layout, cp.def.setup)
    V = use_ale ? u[position_index(layout, 2)] - u[position_index(layout, 1)] : 1.0
    e = only(cp.deactivation.events)
    return (;
        m = [u[mass_index(layout, 1, j)] for j = 1:3],
        T = ds.T[1],
        V,
        P = ds.P[1],
        φ = ds.φ[1],
        created = e.created_volume,
        lost = e.lost_mass,
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
    @testset "on a fixed mesh the receiver keeps its volume and none of the dead cell's gas" begin
        receiver = (m = [0.3, 0.3, 1e-4], T = 400.0)
        r0 = _merge_pair(receiver, (m = [0.0, 0.0, 0.0], T = 400.0); use_ale = false)
        r = _merge_pair(receiver, (m = [0.003, 0.003, 2e-4], T = 500.0); use_ale = false)
        @test r.lost == [0.0, 0.0, 2e-4]
        @test r.m == [0.3 + 0.003, 0.3 + 0.003, 1e-4]
        @test r.created == 0
        # Equal heat capacities: the merged temperature is the mass-weighted
        # mean of what the receiver takes.
        @test r.T ≈ (0.6001 * 400 + 0.006 * 500) / 0.6061 atol = 1e-10
        φ = 1 - (r.m[1] / (300 / r.T) + r.m[2])
        @test r.φ ≈ φ rtol = 1e-10
        @test r.P ≈ (1e-4 / 0.028) * 8.314462618 * r.T / φ rtol = 1e-10
        # Without the vent the receiver would hold three times its gas; with
        # it the pressure moves only with the condensed volume and the
        # temperature.
        @test r.P < 1.05 * r0.P
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
    # At the threshold the cell reacts at 20/27 of the plain kernel's rates.
    u = copy(cp.u0)
    u[mass_index(cp.def.layout, 2, 1)] = cells.threshold * m0
    residual!(du_s, u, cp.def, 0.0, cp.ws)
    residual!(du_p, u, plain.def, 0.0, plain.ws)
    S_s, S_p = cp.ws.contributions.S_mass, plain.ws.contributions.S_mass
    @test S_p[1, 2] < 0
    @test S_s[1, 2] / S_p[1, 2] ≈ 20 / 27 rtol = 1e-12
    @test S_s[2, 2] / S_p[2, 2] ≈ 20 / 27 rtol = 1e-12
    @test S_s[1, 1] == S_p[1, 1]
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
    @test all(e.receiver < e.cell for e in sol.events)
    @test [e.invocation for e in sol.events] == cumsum([1; diff([e.t for e in sol.events]) .> 0])
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

@testset "the consistent sweep after a transfer leaves the velocity and intensity rows at round-off" begin
    q = 2.8e6
    cp = setup(
        _ablator(6; q, radiating = true);
        use_ale = true,
        radiation_model = BEER_LAMBERT,
        deactivation = true,
    )
    def, layout, cells = cp.def, cp.def.layout, cp.deactivation
    u = _deplete!(copy(cp.u0), cp, 6, 0.5 * cells.threshold)
    consistent_state!(u, def, cp.ws, 0.0)
    @test deactivate!(u, cells, def, 0.2) === :transferred
    consistent_state!(u, def, cp.ws, 0.2)
    du = zeros(length(u))
    residual!(du, u, def, 0.2, cp.ws)
    ds, c = cp.ws.derived, cp.ws.contributions
    scale = sum(abs, c.θ .* ds.V_eval) + sum(abs, ds.geometry.face_area .* ds.w)
    @test all(abs(du[velocity_index(layout, nd)]) <= 1e-10 * scale for nd = 1:6)
    @test all(abs(du[intensity_index(layout, f)]) <= 1e-12 * q for f = 1:6)
    @test !layout.mask[intensity_index(layout, 7)]
end

@testset "the consistent sweep reports a velocity row it cannot satisfy" begin
    cp = setup(_ablator(4); use_ale = true)
    u = copy(cp.u0)
    @test Integration._consistent_motion!(u, cp.def, cp.ws, 0.0) === nothing
    u[mass_index(cp.def.layout, 2, 1)] = NaN
    @test Integration._consistent_motion!(u, cp.def, cp.ws, 0.0) isa Tuple
    @test_throws ErrorException consistent_state!(u, cp.def, cp.ws, 0.0)
end

@testset "the domain test rejects proposals through a solve and gives up on a run of constraint rejections" begin
    options = (use_ale = true, radiation_model = BEER_LAMBERT, deactivation = true)
    # The solver-chosen restart keeps the step through an event, so the domain
    # test sees the proposals that follow one; the heuristic restart avoids them.
    policy = rodas5p_policy(KLUFactorization(); restart = :solver_chosen)
    sol = integrate(
        setup(_ablator(8; radiating = true); options...);
        saveat = Float64[],
        integrator_policy = policy,
    )
    @test sol.retcode === :Success
    @test sol.statistics.domain_rejections > 0
    @test sol.statistics.largest_constraint_residual <= 1
    # At tight tolerances a constraint residual stays above its tolerance on
    # two consecutive proposals; with a bound of one the domain test gives up.
    cp = setup(_ablator(8; radiating = true); options..., reltol = 1e-10, abstol = 1e-12)
    sol =
        @test_logs (:warn, r"constraint residual stayed above") match_mode = :any integrate(
            cp;
            saveat = Float64[],
            max_constraint_rejections = 1,
            integrator_policy = policy,
        )
    @test sol.retcode === :Failure
    @test sol.t[end] < 0.8
    @test sol.statistics.domain_rejections > 0
end

@testset "the restart after an event: $(name)" for (name, policy) in (
    ("FBDF", fbdf_policy(KLUFactorization())),
    ("Rodas5P, forced restart", rodas5p_policy(KLUFactorization(); restart = :forced_dt)),
    ("Rodas5P, heuristic restart", rodas5p_policy(KLUFactorization())),
    ("KenCarp4", kencarp4_policy(KLUFactorization())),
)
    steps = Tuple{Float64,Float64}[]
    probe = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> push!(steps, (integrator.tprev, integrator.t));
        save_positions = (false, false),
    )
    sol = Integration.solve(
        _ablator(30);
        use_ale = true,
        deactivation = true,
        abstol = 1e-7,
        reltol = 1e-5,
        integrator_policy = policy,
        callback = probe,
    )
    @test sol.retcode === :Success
    @test length(sol.events) >= 5
    @test [e.cell for e in sol.events] == 30:-1:(31-length(sol.events))
    if policy.restart === :forced_dt
        # The step after each event is the restart step.
        after = [steps[findfirst(s -> s[1] == e.t, steps)] for e in sol.events]
        @test all(s -> s[2] - s[1] ≈ DT_AFTER_EVENT, after)
    end
end

# Two isothermal cells: the outer one holds a hundredth of the inner one's
# mass, of a solid that decomposes 2 % faster, so it crosses its threshold
# just before the inner one and its merge barely delays the inner crossing.
function _staggered_pair()
    ρ, c, R, E = 1000.0, 1000.0, 8.314462618, 1.0e5
    A = 40 * exp(E / (R * 700))
    fast = SolidComponent(:fast, ρ = ρ, c = c, k = 1.0, λ = 1e-6)
    slow = SolidComponent(:slow, ρ = ρ, c = c, k = 1.0, λ = 1e-6)
    gas = GaseousComponent(:gas, M = 0.05, c = c, k = 0.02, λ = 1e-4)
    mat = Material(
        name = :pair,
        components = (fast, slow, gas),
        reactions = (
            Reaction(:fast, 1 => 3; A = 1.02 * A, E = E, h = 0.0),
            Reaction(:slow, 2 => 3; A = A, E = E, h = 0.0),
        ),
        depletion_limiter = DepletionLimiter(threshold = 1e-3),
    )
    L = 2e-3
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L, 2, Val(3)),
        material = mat,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> 700.0,
        ξ_initial = [z -> z > L / 2 ? ρ / 100 : 0.0, z -> z > L / 2 ? 0.0 : ρ, z -> 0.0],
        tspan = (0.0, 0.5),
        output_times = collect(0.0:0.05:0.5),
    )
end

@testset "a located event's time converges with the tolerance: $(moving ? "moving" : "fixed") mesh" for moving in
                                                                                                        (
    true,
    false,
)
    problem = _ablator(2; L = 2e-4, t_end = 0.3)
    events(located, reltol) = Integration.solve(
        problem;
        use_ale = moving,
        deactivation = true,
        abstol = 1e-2 * reltol,
        reltol,
        located_events = located,
    ).events
    reference = first(events(true, 1e-9)).t
    shift(located, reltol) = abs(first(events(located, reltol)).t - reference)
    located = [shift(true, reltol) for reltol in (1e-5, 1e-6, 1e-7)]
    @test all(diff(located) .< 0)
    @test located[1] < shift(false, 1e-5) / 10
    # The located cell goes at its threshold, the end-of-step event below it.
    @test all(e -> isapprox(e.fraction, 1e-2; rtol = 1e-12), events(true, 1e-5))
    @test all(e -> e.fraction < 0.999e-2, events(false, 1e-5))
end

@testset "location and the end-of-step test deactivate the same cells: $(moving ? "moving" : "fixed") mesh" for moving in
                                                                                                                (
    true,
    false,
)
    options = (use_ale = moving, deactivation = true, abstol = 1e-7, reltol = 1e-5)
    discrete = Integration.solve(_ablator(30); options..., located_events = false)
    located = Integration.solve(_ablator(30); options..., located_events = true)
    @test located.retcode === discrete.retcode === :Success
    @test length(located.events) >= 15
    @test [e.cell for e in located.events] == [e.cell for e in discrete.events]
    @test [e.invocation for e in located.events] == [e.invocation for e in discrete.events]
    # Every invocation came from a located root: its crossing cell at the threshold.
    for k in unique(e.invocation for e in located.events)
        fractions = [e.fraction for e in located.events if e.invocation == k]
        @test minimum(fractions) ≈ 1e-2 rtol = 1e-12
    end
    @test all(e -> e.t > 0, located.events)
end

@testset "the end-of-step safety net takes a crossing location cannot see: $(moving ? "moving" : "fixed") mesh" for moving in
                                                                                                                    (
    true,
    false,
)
    cp = setup(
        _ablator(4);
        use_ale = moving,
        deactivation = true,
        abstol = 1e-7,
        reltol = 1e-5,
    )
    cells, layout = cp.deactivation, cp.def.layout
    # The outer cell starts exactly at its threshold, where a located root
    # would leave it: its margin is zero and never crosses zero downward.
    u0 = _deplete!(copy(cp.u0), cp, 4, cells.threshold)
    consistent_state!(u0, cp.def, cp.ws, 0.0)
    @test threshold_margin(cells, layout, u0) == 0.0
    @test !any_below(cells, layout, u0)
    @test below_threshold(cells, layout, u0, 4; slack = 1e-5)
    @test !below_threshold(cells, layout, u0, 3; slack = 1e-5)
    ends = Float64[]
    probe = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> push!(ends, integrator.t);
        save_positions = (false, false),
    )
    sol = integrate(cp; u0, located_events = true, callback = probe)
    @test sol.retcode === :Success
    @test first(sol.events).cell == 4
    @test first(sol.events).t == first(ends)
    @test first(sol.events).fraction < cells.threshold
end

@testset "two crossings in one step are located at their own roots" begin
    ends = Tuple{Float64,Float64}[]
    probe = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> push!(ends, (integrator.tprev, integrator.t));
        save_positions = (false, false),
    )
    options = (deactivation = true, abstol = 1e-7, reltol = 1e-5)
    discrete = Integration.solve(
        _staggered_pair();
        options...,
        located_events = false,
        callback = probe,
    )
    located = Integration.solve(_staggered_pair(); options..., located_events = true)
    # The end-of-step test takes both cells at the end of one step.
    @test discrete.outcome === :burned_out
    @test [e.cell for e in discrete.events] == [2, 1]
    @test [e.invocation for e in discrete.events] == [1, 1]
    step = ends[findfirst(s -> s[2] == first(discrete.events).t, ends)]
    # Location takes each at its own root inside that step.
    @test located.outcome === :burned_out
    @test [e.cell for e in located.events] == [2, 1]
    @test [e.invocation for e in located.events] == [1, 2]
    t1, t2 = (e.t for e in located.events)
    @test step[1] < t1 < t2 < step[2]
    @test all(e -> isapprox(e.fraction, 1e-2; rtol = 1e-12), located.events)
end

@testset "a located event leaves an admissible state with the velocity and intensity rows at round-off" begin
    q = 2.8e6
    cp = setup(
        _ablator(8; q, radiating = true);
        use_ale = true,
        radiation_model = BEER_LAMBERT,
        deactivation = true,
        abstol = 1e-7,
        reltol = 1e-5,
    )
    def, layout, cells = cp.def, cp.def.layout, cp.deactivation
    du = zeros(length(cp.u0))
    seen = Ref(0)
    results = NamedTuple[]
    # Runs at the end of the step a located event cut short, on the state the
    # event and its sweep left.
    probe = DiscreteCallback(
        (u, t, integrator) -> length(cells.events) > seen[],
        integrator -> begin
            seen[] = length(cells.events)
            u, t = integrator.u, integrator.t
            residual!(du, u, def, t, cp.ws)
            ds, c = cp.ws.derived, cp.ws.contributions
            scale =
                sum(abs, c.θ .* ds.V_eval) + sum(abs, ds.geometry.face_area .* ds.w)
            velocity = [
                abs(du[velocity_index(layout, nd)]) / scale for
                nd = 1:n_nodes(layout) if layout.mask[velocity_index(layout, nd)]
            ]
            intensity = [
                abs(du[intensity_index(layout, f)]) / q for
                f = 1:(n_cells(layout)+1) if layout.mask[intensity_index(layout, f)]
            ]
            push!(
                results,
                (;
                    t,
                    admissible = !inadmissible(u, layout, cp.abstol),
                    velocity = maximum(velocity),
                    intensity = maximum(intensity),
                ),
            )
        end;
        save_positions = (false, false),
    )
    sol = integrate(cp; located_events = true, callback = probe)
    @test sol.retcode === :Success
    @test length(results) == length(unique(e.invocation for e in sol.events)) >= 3
    @test [r.t for r in results] == unique(e.t for e in sol.events)
    @test all(r -> r.admissible, results)
    @test all(r -> r.velocity <= 1e-10, results)
    @test all(r -> r.intensity <= 1e-12, results)
end

@testset "a void under a live cell ends the solve" begin
    for moving in (true, false)
        sol = Integration.solve(
            _ablator(6; from_bottom = true, t_end = 5.0);
            use_ale = moving,
            deactivation = true,
            abstol = 1e-7,
            reltol = 1e-5,
        )
        @test sol.retcode === :Terminated
        @test sol.outcome === :interior
        @test isempty(sol.events)
        @test sol.t[end] < 5.0
    end
end

@testset "deactivation under UMFPACK" begin
    options = (use_ale = true, deactivation = true, abstol = 1e-7, reltol = 1e-5)
    klu = Integration.solve(_ablator(30); options...)
    umfpack =
        Integration.solve(_ablator(30); options..., linear_solver = UMFPACKFactorization())
    @test umfpack.retcode === :Success
    @test [e.cell for e in umfpack.events] == [e.cell for e in klu.events]
    @test [e.t for e in umfpack.events] ≈ [e.t for e in klu.events] rtol = 1e-6
end
