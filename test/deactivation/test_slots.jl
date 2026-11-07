using Test
using Pyrolysis
using Pyrolysis.Materials: RADIATIVE_TRANSFER
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.Kernels:
    ReactionKernel,
    ConductionKernel,
    GasTransportKernel,
    ThermalBoundaryKernel,
    RadiativeTransferKernel,
    RelativeFluxKernel,
    find_kernel
using Pyrolysis.Closures: NoCorrection, EquidistributionCorrection
using Pyrolysis.Deactivation
using Pyrolysis.Integration
using Pyrolysis.Integration: setup, integrate
using OrdinaryDiffEq: DiscreteCallback

# A column of two materials on a graded mesh, every cell with its own
# thickness, temperature and composition, so that cells differ in every table:
# `h` the thicknesses, `T` the temperatures, `layer` 1 for the reacting solid
# and 2 for the inert one, `ρ` the density of each cell's solid. The value of a
# field in a cell follows from the cell the position falls in, so removing a
# cell's entries gives the deck without it.
function _column(h, T, layer, ρ; radiating = false)
    a = SolidComponent(:a; φ = 0.3, ρ = 800.0, c = 1000.0, k = 0.2, λ = 1e-6, α = 50.0)
    b = SolidComponent(:b; φ = 0.6, ρ = 900.0, c = 1300.0, k = 0.1, λ = 1e-6, α = 50.0)
    gas = GaseousComponent(:gas, M = 0.03, c = 1000.0, k = 0.02)
    material = Material(
        name = :column,
        components = (a, b, gas),
        reactions = (Reaction(:burn, 1 => 3; A = 1e5, E = 1e5, h = 0.0),),
    )
    z = cumsum(vcat(0.0, h))
    cell(x) = clamp(searchsortedlast(z, x), 1, length(h))
    return PyrolysisProblem(
        mesh = line_mesh(z),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = radiating ? RadiativeFluxBC(1e4) : HeatFluxBC(1e4),
                bottom = AdiabaticBC(),
            ),
            mass = (
                top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = x -> T[cell(x)],
        ξ_initial = [
            x -> layer[cell(x)] == 1 ? ρ[cell(x)] : 0.0,
            x -> layer[cell(x)] == 2 ? ρ[cell(x)] : 0.0,
            x -> 0.0,
        ],
        tspan = (0.0, 1.0),
    )
end

const _H = [2.0, 1.5, 1.1, 0.9, 0.7, 0.6, 0.5, 0.45] .* 1e-4
const _T = [640.0, 630.0, 625.0, 610.0, 600.0, 590.0, 585.0, 570.0]
const _L = [2, 2, 2, 1, 1, 1, 1, 1]
const _R = [500.0, 480.0, 450.0, 400.0, 390.0, 370.0, 360.0, 340.0]

_set_up(p, energy; model = SURFACE_ABSORPTION) = setup(
    p;
    energy_state = energy,
    use_ale = true,
    deactivation = true,
    correction = NoCorrection(),
    radiation_model = model,
)

# The registered tables by name: (kind, live array).
function _named_tables(cp)
    r = cp.def.setup.reference
    k = cp.def.kernels
    out = Dict{String,Tuple{Symbol,Any}}(
        "volume" => (:cell, r.volume),
        "solid_mass" => (:cell, r.solid_mass),
        "mass_scale" => (:cell, r.mass_scale),
        "volume_floor" => (:cell, r.volume_floor),
        "capacity_floor" => (:cell, r.capacity_floor),
        "thickness" => (:cell, r.thickness),
        "deactivation reference" => (:cell, cp.deactivation.reference),
        "conduction floor" => (:face, find_kernel(ConductionKernel, k).d_floor),
        "boundary floor" => (:face, find_kernel(ThermalBoundaryKernel, k).d_floor),
        "lineage" => (:matrix, find_kernel(ReactionKernel, k).reference),
    )
    gas = find_kernel(GasTransportKernel, k)
    gas === nothing || (out["gas floor"] = (:face, gas.d_floor))
    transfer = find_kernel(RadiativeTransferKernel, k)
    transfer === nothing || (out["transfer floor"] = (:matrix, transfer.d_floor))
    return out
end

# What the slot `s` of a table of `kind` holds with `origin` from the initial `a`.
_at(::Val{:cell}, a, origin, s) = a[origin[s]]
_at(::Val{:matrix}, a, origin, s) = a[:, origin[s]]
_at(::Val{:face}, a, origin, s) = a[origin[s]+1]

_kill!(cells, layout, slot) = (
    layout.mask[cell_rows(layout, slot)] .= false;
    layout.mask[node_rows(layout, slot + 1)] .= false
)

@testset "the slots start with every cell in its own and the writer changes nothing there: $(energy) form" for energy in
                                                                                                               (
    :temperature,
    :enthalpy,
)
    cp = _set_up(_column(_H, _T, _L, _R), energy)
    t = cp.deactivation.tables
    @test t.origin == 1:8
    before = Dict(n => copy(a) for (n, (_, a)) in _named_tables(cp))
    tol = copy(cp.abstol)
    permute_slots!(t, cp.def.layout)
    @test all(before[n] == a for (n, (_, a)) in _named_tables(cp))
    @test cp.abstol == tol
end

@testset "a hand-built permutation writes every table from its initial copy: $(energy) form" for energy in
                                                                                                 (
    :temperature,
    :enthalpy,
)
    cp = _set_up(_column(_H, _T, _L, _R), energy)
    layout, cells = cp.def.layout, cp.deactivation
    t = cells.tables
    initial = Dict(n => copy(a) for (n, (_, a)) in _named_tables(cp))
    tol0 = copy(cp.abstol)
    origin = [1, 2, 4, 5, 6, 7, 8, 3]
    _kill!(cells, layout, 8)
    arrange_slots!(t, layout, origin, Tuple{Int,Int,Float64}[])
    @test t.origin == origin
    for (n, (kind, a)) in _named_tables(cp)
        for s = 1:7
            @test (kind === :matrix ? a[:, s] : kind === :face ? a[s+1] : a[s]) ==
                  _at(Val(kind), initial[n], origin, s)
        end
    end
    # The rows move with their cells; a node between new neighbours has the
    # tolerances of the thinner of them.
    for s = 1:7
        @test cp.abstol[cell_rows(layout, s)] == tol0[cell_rows(layout, origin[s])]
    end
    h = cp.def.setup.reference.thickness
    gap = position_index(layout, 3)
    @test cp.abstol[gap] ≈ cp.reltol * 1e-2 * min(h[2], h[3]) rtol = 1e-14
    @test cp.abstol[position_index(layout, 5)] == tol0[position_index(layout, origin[5])]
end

@testset "the reset restores the identity and every table" begin
    cp = _set_up(_column(_H, _T, _L, _R), :temperature)
    layout, cells = cp.def.layout, cp.deactivation
    before = Dict(n => copy(a) for (n, (_, a)) in _named_tables(cp))
    tol = copy(cp.abstol)
    kernel = find_kernel(ReactionKernel, cp.def.kernels)
    _kill!(cells, layout, 8)
    update = [(1, 2, 3.0)]
    arrange_slots!(cells.tables, layout, [2, 1, 4, 5, 6, 7, 8, 3], update)
    @test kernel.reference[1, 1] == 3.0
    reset_cells!(cells, layout)
    @test cells.tables.origin == 1:8
    @test isempty(cells.tables.held)
    @test all(layout.mask)
    @test all(before[n] == a for (n, (_, a)) in _named_tables(cp))
    @test cp.abstol == tol
end

@testset "a replay reproduces the permutation and the updates of its time" begin
    cp = _set_up(_column(_H, _T, _L, _R), :enthalpy)
    layout, cells = cp.def.layout, cp.deactivation
    kernel = find_kernel(ReactionKernel, cp.def.kernels)
    event(t, cell, slot, origin, updates) = DeactivationEvent(
        t,
        1,
        cell,
        slot,
        0.001,
        vcat(collect(cell_rows(layout, slot)), collect(node_rows(layout, slot + 1))),
        0.0,
        [0.0, 0.0, 0.0],
        0.0,
        0.0,
        origin,
        updates,
        false,
    )
    push!(cells.events, event(1.0, 3, 8, [1, 2, 4, 5, 6, 7, 8, 3], [(1, 4, 2.0)]))
    push!(cells.events, event(2.0, 6, 7, [1, 2, 4, 5, 7, 8, 6, 3], [(1, 5, 7.0)]))
    replay_events!(cells, layout, 0.5)
    @test cells.tables.origin == 1:8
    @test all(layout.mask)
    @test all(iszero, kernel.reference)
    replay_events!(cells, layout, 1.5)
    @test cells.tables.origin == [1, 2, 4, 5, 6, 7, 8, 3]
    @test top_live_cell(layout) == 7
    @test kernel.reference[1, 3] == 2.0 && kernel.reference[1, 4] == 0.0
    @test cp.def.setup.reference.volume[3] == first(cells.tables.cell)[1][4]
    replay_events!(cells, layout, 2.0)
    @test top_live_cell(layout) == 6
    @test kernel.reference[1, 3] == 2.0 && kernel.reference[1, 4] == 7.0
    replay_events!(cells, layout, 1.5)
    @test kernel.reference[1, 4] == 0.0
    replay_events!(cells, layout, Inf)
    @test cells.tables.origin == [1, 2, 4, 5, 7, 8, 6, 3]
    @test kernel.reference[1, 4] == 7.0
end

@testset "the tolerances and the caloric floors are one array each" begin
    cp = _set_up(_column(_H, _T, _L, _R), :enthalpy)
    cells = cp.deactivation
    r = cp.def.setup.reference
    @test DomainCheck(cp).abstol === cp.abstol
    @test only(cells.tables.rows).live === cp.abstol
    @test cells.caloric.δ === r.mass_scale
    @test cells.caloric.floor === r.capacity_floor
    @test cp.def.setup.caloric.δ === r.mass_scale
    @test cp.def.setup.caloric.floor === r.capacity_floor
    seen = Ref{Any}(nothing)
    callback = DiscreteCallback(
        (u, t, integrator) -> true,
        integrator -> (seen[] = integrator.opts.abstol; nothing);
        save_positions = (false, false),
    )
    integrate(cp; callback, dt = 1e-6, maxiters = 3, saveat = [1.0])
    @test seen[] === cp.abstol
end

# The table of a problem with the cell `k` taken out of the mesh and the same
# table after the cell left its slot: equal to round-off, slot by slot, so that
# a table a change adds and forgets to register fails here.
function _compare_with_fresh(cp, fresh, top)
    layout, flayout = cp.def.layout, fresh.def.layout
    a, b = _named_tables(cp), _named_tables(fresh)
    @test keys(a) == keys(b)
    for (n, (kind, table)) in a
        other = b[n][2]
        if kind === :cell
            @test table[1:top] ≈ other[1:top] rtol = 1e-12
        elseif kind === :matrix
            @test table[:, 1:top] ≈ other[:, 1:top] rtol = 1e-12
        else
            @test table[1:(top+1)] ≈ other[1:(top+1)] rtol = 1e-12
        end
    end
    for s = 1:top
        @test cp.abstol[cell_rows(layout, s)] ≈ fresh.abstol[cell_rows(flayout, s)] rtol =
            1e-12
    end
    for s = 1:(top+1)
        @test cp.abstol[node_rows(layout, s)] ≈ fresh.abstol[node_rows(flayout, s)] rtol =
            1e-12
    end
end

@testset "a removed cell's tables equal a fresh setup without it: cell $k, $(energy) form, $(model)" for k in
                                                                                                         (
        1,
        3,
        5,
        7,
    ),
    energy in (:temperature, :enthalpy),
    model in (SURFACE_ABSORPTION, RADIATIVE_TRANSFER)

    radiating = model !== SURFACE_ABSORPTION
    cp = _set_up(_column(_H, _T, _L, _R; radiating), energy; model)
    keep = [i for i = 1:8 if i != k]
    fresh =
        _set_up(_column(_H[keep], _T[keep], _L[keep], _R[keep]; radiating), energy; model)
    layout, cells = cp.def.layout, cp.deactivation
    _kill!(cells, layout, 8)
    arrange_slots!(cells.tables, layout, vcat(keep, k), Tuple{Int,Int,Float64}[])
    _compare_with_fresh(cp, fresh, 7)
end

@testset "a cell a pass takes from under live cells leaves the tables of a fresh setup: cell $k, $(energy) form, $(model)" for k in
                                                                                                                               (
        1,
        3,
        5,
        7,
    ),
    energy in (:temperature, :enthalpy),
    model in (SURFACE_ABSORPTION, RADIATIVE_TRANSFER)

    radiating = model !== SURFACE_ABSORPTION
    cp = _set_up(_column(_H, _T, _L, _R; radiating), energy; model)
    keep = [i for i = 1:8 if i != k]
    fresh =
        _set_up(_column(_H[keep], _T[keep], _L[keep], _R[keep]; radiating), energy; model)
    layout, cells = cp.def.layout, cp.deactivation
    u = copy(cp.u0)
    # Cell k empties under spent cells.
    for j = 1:n_components(layout)
        u[mass_index(layout, k, j)] *= 0.5 * cells.threshold
        for i = (k+1):8
            u[mass_index(layout, i, j)] *= 2 * cells.threshold
        end
    end
    @test deactivate!(u, cells, cp.def, 0.1) === :transferred
    @test cells.tables.origin == vcat(keep, k)
    @test top_live_cell(layout) == 7
    _compare_with_fresh(cp, fresh, 7)
end

@testset "the slots of a bound correction follow its cells" begin
    h = fill(1e-4, 10)
    layer = [2, 2, 2, 2, 1, 1, 1, 1, 1, 1]
    T = fill(600.0, 10)
    ρ = fill(400.0, 10)
    cp = setup(
        _column(h, T, layer, ρ);
        use_ale = true,
        deactivation = true,
        correction = EquidistributionCorrection(layers = :all),
    )
    layout, cells = cp.def.layout, cp.deactivation
    corr = cp.def.closure.correction
    @test corr.origin === cells.tables.origin
    gain = copy(corr.gain)
    @test iszero(gain[5]) && all(!iszero, gain[[2, 3, 4, 6, 7, 8, 9, 10]])
    kernel = find_kernel(RelativeFluxKernel, cp.def.kernels)
    η = copy(kernel.η)
    floors = copy(kernel.d_floor)
    # The cell at the interface and the one beneath it leave: the layer
    # interface moves down with the cells above it.
    _kill!(cells, layout, 10)
    _kill!(cells, layout, 9)
    origin = [1, 2, 3, 6, 7, 8, 9, 10, 4, 5]
    arrange_slots!(cells.tables, layout, origin, Tuple{Int,Int,Float64}[])
    @test corr.gain[4] == 0
    @test corr.gain[[2, 3]] == gain[[2, 3]]
    @test corr.gain[5:8] == gain[7:10]
    @test kernel.η[4] == η[6]
    @test kernel.d_floor[1, 5] == floors[1, 7] && kernel.d_floor[2, 4] == floors[2, 6]
    # Two layers of live cells, each in its own contiguous slots.
    @test Pyrolysis.Closures._live_cells(corr, layout, 1) == 1:3
    @test Pyrolysis.Closures._live_cells(corr, layout, 2) == 4:8
    reset_cells!(cells, layout)
    @test corr.gain == gain
    @test kernel.d_floor == floors
end
