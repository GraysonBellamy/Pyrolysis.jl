using Test
using LinearAlgebra
using SparseArrays
using StaticArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.States
using Pyrolysis.States: state_length
using Pyrolysis.Assembly: residual!
using Pyrolysis.Closures: PotentialClosure, rank_check
using Pyrolysis.ColoredJacobian
using Pyrolysis.Deactivation: element_death!
using Pyrolysis.Integration: setup

function _wood()
    virgin = SolidComponent(
        :wood,
        ρ_intrinsic = 1500.0,
        λ = 0.0,
        ρ = 500.0,
        c = 1500.0,
        k = 0.15,
        ε = 0.9,
    )
    char = SolidComponent(
        :char,
        ρ_intrinsic = 1800.0,
        λ = 0.0,
        ρ = 150.0,
        c = 1100.0,
        k = 0.08,
        ε = 0.95,
    )
    gas = GaseousComponent(:gas, M = 0.06, c = 1800.0, k = 0.03)
    rxn = Reaction(
        :pyrolysis,
        1 => (2, 3);
        A = 1e10,
        E = 1.5e5,
        h = -3e5,
        yields = (0.25, 0.75),
    )
    return Material(name = :wood, components = (virgin, char, gas), reactions = (rxn,))
end

_rot(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
const _A = _rot(π / 6) * (@SMatrix [4.0 0.0; 0.0 1.0]) * _rot(π / 6)'
const _slip = (:left => 1, :right => 1)
const _names = ((:left, :right), (:bottom, :top))
_inside(x) = SVector(x[1] + 0.001sin(π * x[1] / 0.03) * sin(2π * x[2] / 0.03), x[2])

# The surface balance on the exposed face: incident flux, re-radiation and
# convection.
_bc_set() = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 5e4) +
              RadiativeBC(T_inf = 300.0) +
              ConvectiveBC(h = 10.0, T_inf = 300.0),
        bottom = AdiabaticBC(),
        left = AdiabaticBC(),
        right = AdiabaticBC(),
    ),
)

_problem(mesh) = PyrolysisProblem(
    mesh = mesh,
    material = _wood(),
    bc_set = _bc_set(),
    T_initial = x -> 600.0 + 1e4 * x[2] + 3e3 * x[1],
    ξ_initial = [x -> 450.0 - 2e3 * x[1], x -> 15.0, x -> 0.5],
    tspan = (0.0, 1.0),
)

_options = (use_ale = true, slip = _slip, tensor = _A)

# The mesh of the live cells of `full` (`live` flags them), its nodes
# renumbered in order: a side another live cell shares stays interior, a side
# a dead cell shared is tagged `:top`, the exposed tag, and every other side
# keeps its tag. Returns the mesh and the maps from its cells and nodes to
# the full mesh's.
function _without(full::Mesh, live)
    topo = full.topology
    cells = findall(live)
    nodes = sort!(unique(nd for i in cells for nd in topo.cells[i].nodes))
    renumber = Dict(nd => k for (k, nd) in enumerate(nodes))
    tag(r, k, s) = begin
        I = cells[r]
        f = topo.cells[I].faces[2(k-1)+s+1]
        o, nb = topo.faces[f].cells
        other = o == I ? nb : o
        other != 0 && !live[other] ? :top : topo.faces[f].tag
    end
    mesh = cell_mesh(
        full.positions[nodes],
        [Tuple(renumber[nd] for nd in topo.cells[i].nodes) for i in cells];
        tag,
    )
    return mesh, cells, nodes
end

# A state of the full mesh restricted to the reduced mesh's unknowns.
function _restrict(u, full_layout, layout, cells, nodes)
    v = zeros(state_length(layout))
    D = dimension(layout)
    for (r, I) in enumerate(cells)
        for j = 1:n_components(layout)
            v[mass_index(layout, r, j)] = u[mass_index(full_layout, I, j)]
        end
        v[energy_index(layout, r)] = u[energy_index(full_layout, I)]
        v[multiplier_index(layout, r)] = u[multiplier_index(full_layout, I)]
    end
    for (k, nd) in enumerate(nodes), a = 1:D
        v[position_index(layout, k, a)] = u[position_index(full_layout, nd, a)]
    end
    return v
end

# The largest difference between the full residual's live rows and the
# reduced residual's, relative to the largest of each row kind.
function _live_rows_difference(cp, cpr, u, cells, nodes; t = 0.5)
    L, Lr = cp.def.layout, cpr.def.layout
    du = zeros(state_length(L))
    residual!(du, u, cp.def, t, cp.ws)
    v = _restrict(u, L, Lr, cells, nodes)
    dv = zeros(state_length(Lr))
    residual!(dv, v, cpr.def, t, cpr.ws)
    pairs = Dict{Symbol,Vector{Tuple{Float64,Float64}}}()
    for (r, I) in enumerate(cells)
        for j = 1:n_components(L)
            push!(
                get!(pairs, :mass, []),
                (du[mass_index(L, I, j)], dv[mass_index(Lr, r, j)]),
            )
        end
        push!(get!(pairs, :energy, []), (du[energy_index(L, I)], dv[energy_index(Lr, r)]))
        push!(
            get!(pairs, :balance, []),
            (du[multiplier_index(L, I)], dv[multiplier_index(Lr, r)]),
        )
    end
    for (k, nd) in enumerate(nodes), a = 1:2
        push!(
            get!(pairs, :position, []),
            (du[position_index(L, nd, a)], dv[position_index(Lr, k, a)]),
        )
    end
    return Dict(
        kind =>
            maximum(abs(a - b) for (a, b) in p) /
            max(maximum(abs(b) for (a, b) in p), 1e-300) for (kind, p) in pairs
    )
end

# The pattern contains the dense Jacobian's nonzeros and the colored Jacobian
# equals it per entry.
function _contained(cp, u; t = 0.5)
    def, ws, cj = cp.def, cp.ws, cp.jacobian
    du = zeros(length(u))
    residual!(du, u, def, t, ws)
    @test all(isfinite, du)
    Jd = dense_jacobian(def, ws, u, t)
    @test isempty(missing_entries(Jd, cj.pattern))
    J = jacobian!(jac_prototype(cj), cj, u, t)
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
    return Jd
end

@testset "element death on a 3×3 mesh, a corner cell then an edge cell sharing its nodes: $(energy) form" for energy in
                                                                                                              (
    :temperature,
    :enthalpy,
)
    mesh = box_mesh(0:0.01:0.03, 0:0.01:0.03; transform = _inside)
    cp = setup(_problem(mesh); energy_state = energy, _options...)
    layout = cp.def.layout
    topo = mesh.topology
    pattern = copy(cp.jacobian.pattern)
    live = trues(n_cells(topo))
    u = copy(cp.u0)
    for dead in (9, 8)
        element_death!(layout, topo, (dead,))
        live[dead] = false
        ds = cp.ws.derived
        derive!(ds, u, layout, cp.def.setup)
        @test ds.active == live
        # The node only the dead cells held leaves; the nodes they share with
        # live cells stay.
        held = Set(nd for i in findall(live) for nd in topo.cells[i].nodes)
        for nd = 1:n_nodes(topo)
            @test layout.mask[position_index(layout, nd, 1)] == (nd in held)
        end
        @test all(!layout.mask[r] for r in cell_rows(layout, dead))
        # Every live cell admissible, the closure determined on the new boundary.
        @test all(
            ds.V[i] > 0 && all(>(0), ds.m[:, i]) && ds.T[i] > 0 for i in findall(live)
        )
        @test rank_check(cp.def.closure, topo, ds.geometry; active = live)
        # The live rows are those of the mesh without the dead cells, its new
        # faces under the exposed condition.
        reduced, cells, nodes = _without(mesh, live)
        cpr = setup(_problem(reduced); energy_state = energy, _options...)
        @test cpr.def.closure isa PotentialClosure
        hot = copy(u)
        for i = 1:n_cells(topo)
            hot[energy_index(layout, i)] = energy_state_value(
                energy_state(layout),
                cp.def.setup,
                900.0 + 10i,
                Tuple(u[mass_index(layout, i, j)] for j = 1:3),
                i,
            )
        end
        for state in (u, hot)
            difference = _live_rows_difference(cp, cpr, state, cells, nodes)
            @test all(<(1e-12), values(difference))
        end
        # The pattern does not change, and still contains the Jacobian.
        @test cp.jacobian.pattern == pattern
        _contained(cp, u)
        _contained(cp, hot)
    end
end
