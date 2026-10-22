# ============================================================================
# Sparsity pattern from connectivity
# ============================================================================

"""
	jacobian_pattern(def::CoreDefinition) -> SparseMatrixCSC{Bool,Int}
	jacobian_pattern(layout, topology; closure = nothing, reach = 1)

The structural pattern of the residual's Jacobian, from the connectivity
alone: the cells within `reach` faces of a cell, their nodes, the nodes'
cells. The definition supplies the closure, with its correction, and the
reach of the cell rows: two faces where a kernel reads a neighbour's
neighbour (inventory fluxes relative to the mesh, whose face states are
reconstructed from two cells upwind, `RelativeFluxKernel`; conduction with
its cross-diffusion part, which reads the neighbours' least-squares
gradients), one otherwise.

- A cell row (a mass, the energy, the history pair) reaches the masses,
  energy and history of the cells within the reach and, on a moving mesh,
  the positions of their nodes: conduction reads the neighbouring cell
  centers, the concentrations the neighbours' volumes, the relative fluxes
  everything the corrections of the cell's nodes read; under a prescribed
  mesh velocity with relative fluxes, the velocities of its nodes. On a line
  it reads the intensities of its two faces and, under a variable cross
  section, the area.
- A face-intensity row (a line) reaches the intensities of nodes
  `f − 1 … f + 1`, the unknowns of the cells on either side of its face (the
  attenuation of the cell it leaves, the absorptivity of the exposed cell),
  on a moving mesh the positions of those nodes, and the area.
- A node-position row reaches its node's velocity, or, when the closure's
  unknowns are the cells' multipliers, the multipliers of the cells around
  the node and the positions of their nodes (`w = W⁻¹ Dᵀ φ`).
- A velocity row reaches its node's velocity and position. On a line a
  balance row also reaches the velocities and positions of its cell's other
  node, the masses and energy of its cell (through the strain rate), the
  area and the lateral rate, and under a [`LayerCorrection`](@ref) what the
  corrections of its cell's two nodes read: the positions of the nodes of
  the cells beside it and, when the monitor reads the masses, the masses of
  those cells. On a column of more than one dimension a lateral or pinned
  component reaches itself, a tied component its leader's, a balance row
  every velocity component and position of its cell's nodes and its cell's
  masses and energy; a prescribed velocity its node's position.
- A multiplier row reaches the multipliers of the cells that share a node
  with its cell, the positions of their nodes (through the velocities of its
  cell's nodes), and its cell's masses and energy.
- A deformation-gradient row also reaches the velocities of its cell's
  nodes, or the multipliers and node positions those velocities read.
- Under a variable cross section the area row reaches every cell's history
  pair and the area, and the lateral-rate row everything the history rows
  reach, the history pair, the area and the lateral rate: the two dense rows
  of [`dense_rows`](@ref).

Every row holds its diagonal. The pattern is fixed for a solve: dead rows
keep their entries and take zero values.
"""
function jacobian_pattern(def::CoreDefinition)
    relative = find_kernel(RelativeFluxKernel, def.kernels) !== nothing
    conduction = find_kernel(ConductionKernel, def.kernels)
    cross = conduction !== nothing && conduction.tensor !== nothing
    return jacobian_pattern(
        def.layout,
        def.setup.topology;
        closure = def.closure,
        reach = relative || cross ? 2 : 1,
        reconstruction = relative,
    )
end

function jacobian_pattern(
    layout::StateLayout,
    topo::MeshTopology;
    closure = nothing,
    reach::Int = 1,
    reconstruction::Bool = false,
)
    N = n_cells(layout)
    L = state_length(layout)
    D = dimension(layout)
    moving, bl, chi = is_moving(layout), has_intensities(layout), has_chi(layout)
    multipliers = moving && closure_unknowns(layout) === :multipliers
    correction = closure === nothing ? NoCorrection() : closure.correction
    adjacency = _adjacency(topo)
    rows = Int[]
    cols = Int[]
    add!(r, c) = (push!(rows, r); push!(cols, c))
    add_positions!(r, nd) =
        for k = 1:D
            add!(r, position_index(layout, nd, k))
        end
    for i = 1:N
        near = _within(adjacency, i, reach)
        for r in _physical_rows(layout, i)
            for k in near, c in _physical_rows(layout, k)
                add!(r, c)
            end
            if moving
                for nd in _nodes_of(topo, near)
                    add_positions!(r, nd)
                end
                if reconstruction && correction isa PrescribedCorrection
                    for nd in topo.cells[i].nodes, k = 1:D
                        add!(r, velocity_index(layout, nd, k))
                    end
                end
            end
            if bl
                add!(r, intensity_index(layout, i))
                add!(r, intensity_index(layout, i + 1))
            end
            chi && add!(r, area_index(layout))
        end
    end
    if bl
        for f = 1:(N+1)
            r = intensity_index(layout, f)
            for nd = max(1, f-1):min(N+1, f+1)
                add!(r, intensity_index(layout, nd))
                moving && add!(r, position_index(layout, nd))
            end
            for k = max(1, f-1):min(N, f), c in cell_rows(layout, k)
                add!(r, c)
            end
            chi && add!(r, area_index(layout))
        end
    end
    if has_deformation(layout)
        for i = 1:N, b = 1:D, a = 1:D
            r = deformation_index(layout, i, a, b)
            if multipliers
                around = sort!(
                    unique(
                        j for nd in topo.cells[i].nodes for j in adjacency.node_cells[nd]
                    ),
                )
                for j in around
                    add!(r, multiplier_index(layout, j))
                end
                for m in _nodes_of(topo, around)
                    add_positions!(r, m)
                end
            else
                for nd in topo.cells[i].nodes, k = 1:D
                    add!(r, velocity_index(layout, nd, k))
                end
            end
        end
    end
    if multipliers
        node_cells = adjacency.node_cells
        for nd = 1:n_nodes(layout)
            around = node_cells[nd]
            for k = 1:D
                r = position_index(layout, nd, k)
                for j in around
                    add!(r, multiplier_index(layout, j))
                end
                for m in _nodes_of(topo, around)
                    add_positions!(r, m)
                end
            end
        end
        for i = 1:N
            r = multiplier_index(layout, i)
            around = sort!(unique(j for nd in topo.cells[i].nodes for j in node_cells[nd]))
            for j in around
                add!(r, multiplier_index(layout, j))
            end
            for m in _nodes_of(topo, around)
                add_positions!(r, m)
            end
            for j = 1:n_components(layout)
                add!(r, mass_index(layout, i, j))
            end
            add!(r, energy_index(layout, i))
        end
    elseif moving
        _velocity_pattern!(add!, layout, topo, closure, correction, adjacency)
    end
    if chi
        history = Int[]
        for i = 1:N, k = 1:2
            push!(history, history_index(layout, i, k))
        end
        is_history = falses(L)
        is_history[history] .= true
        reached = Set(c for (r, c) in zip(rows, cols) if is_history[r])
        for c in history
            add!(area_index(layout), c)
        end
        add!(area_index(layout), area_index(layout))
        for c in union(reached, history)
            add!(lateral_rate_index(layout), c)
        end
        add!(lateral_rate_index(layout), area_index(layout))
    end
    for r = 1:L
        add!(r, r)
    end
    return sparse(rows, cols, trues(length(rows)), L, L, |)
end

# The cell-to-cell adjacency across faces and the cells around each node.
function _adjacency(topo::MeshTopology)
    N = n_cells(topo)
    neighbours = [Int[] for _ = 1:N]
    for face in topo.faces
        o, nb = face.cells
        nb == 0 && continue
        push!(neighbours[o], nb)
        push!(neighbours[nb], o)
    end
    node_cells = [Int[] for _ = 1:n_nodes(topo)]
    for i = 1:N, nd in topo.cells[i].nodes
        push!(node_cells[nd], i)
    end
    return (; neighbours, node_cells)
end

# The cells within `reach` faces of cell `i`, `i` included, sorted.
function _within(adjacency, i::Int, reach::Int)
    near = [i]
    front = [i]
    for _ = 1:reach
        front = [j for k in front for j in adjacency.neighbours[k] if !(j in near)]
        unique!(front)
        append!(near, front)
    end
    return sort!(near)
end

# The nodes of the cells `cells`, sorted.
_nodes_of(topo::MeshTopology, cells) =
    sort!(unique(nd for i in cells for nd in topo.cells[i].nodes))

# The rows of a cell's block the kernels write: the masses, the energy state
# and the history pair; the multiplier's row is the closure's.
function _physical_rows(layout::StateLayout, i::Int)
    rows = cell_rows(layout, i)
    multipliers = is_moving(layout) && closure_unknowns(layout) === :multipliers
    return multipliers ? rows[1:(end-1)] : rows
end

# The position and velocity rows of a mesh whose closure's unknowns are the
# node velocities.
function _velocity_pattern!(
    add!,
    layout,
    topo::MeshTopology{1},
    closure,
    correction,
    adjacency,
)
    N = n_cells(layout)
    balance = _balance_cells(closure, N)
    for nd = 1:(N+1)
        add!(position_index(layout, nd), velocity_index(layout, nd))
        r = velocity_index(layout, nd)
        add!(r, position_index(layout, nd))
        k = balance[nd]
        k == 0 && continue
        for m in topo.cells[k].nodes
            add!(r, velocity_index(layout, m))
            add!(r, position_index(layout, m))
        end
        for j = 1:n_components(layout)
            add!(r, mass_index(layout, k, j))
        end
        add!(r, energy_index(layout, k))
        _correction_reach!(add!, r, layout, topo, correction, adjacency, k)
        if has_chi(layout)
            add!(r, area_index(layout))
            add!(r, lateral_rate_index(layout))
        end
    end
    return nothing
end

function _velocity_pattern!(
    add!,
    layout,
    topo::MeshTopology{D},
    closure,
    correction,
    adjacency,
) where {D}
    for nd = 1:n_nodes(layout), k = 1:D
        add!(position_index(layout, nd, k), velocity_index(layout, nd, k))
        r = velocity_index(layout, nd, k)
        add!(r, velocity_index(layout, nd, k))
        if correction isa PrescribedCorrection
            for m = 1:D
                add!(r, position_index(layout, nd, m))
            end
        end
    end
    correction isa PrescribedCorrection && return nothing
    for (nd, leader) in closure.ties
        add!(velocity_index(layout, nd, D), velocity_index(layout, leader, D))
    end
    for k in closure.order
        r = velocity_index(layout, closure.row_node[k], D)
        for m in topo.cells[k].nodes, a = 1:D
            add!(r, velocity_index(layout, m, a))
            add!(r, position_index(layout, m, a))
        end
        for j = 1:n_components(layout)
            add!(r, mass_index(layout, k, j))
        end
        add!(r, energy_index(layout, k))
    end
    return nothing
end

# What the balance row `r` of cell `k` reads through the corrections of the
# cell's two nodes: the positions of the nodes of the cells beside it and,
# when the monitor reads the masses, those cells' masses.
_correction_reach!(add!, r, layout::StateLayout, topo, correction, adjacency, k) = nothing
function _correction_reach!(
    add!,
    r,
    layout::StateLayout,
    topo,
    corr::LayerCorrection,
    adjacency,
    k::Int,
)
    near = _within(adjacency, k, 1)
    for nd in _nodes_of(topo, near)
        add!(r, position_index(layout, nd))
    end
    iszero(corr.β) && return nothing
    for i in near, j = 1:n_components(layout)
        add!(r, mass_index(layout, i, j))
    end
    return nothing
end

# The cell whose balance each node's velocity row carries (zero for the pinned
# node, and for every node under a prescribed velocity).
_balance_cells(::Nothing, N) = zeros(Int, N + 1)
function _balance_cells(cl::BalanceRowClosure, N)
    out = zeros(Int, N + 1)
    cl.correction isa PrescribedCorrection && return out
    for k = 1:N
        out[cl.row_node[k]] = k
    end
    return out
end

"""
	dense_rows(layout) -> Vector{Int}

The rows that read the whole column — the area and lateral-rate rows of a
variable cross section — and are therefore kept out of the column coloring:
one dense row would make every pair of columns conflict.
"""
dense_rows(layout::StateLayout) =
    has_chi(layout) ? [area_index(layout), lateral_rate_index(layout)] : Int[]
