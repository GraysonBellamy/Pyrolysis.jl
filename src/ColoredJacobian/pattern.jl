# ============================================================================
# Sparsity pattern from connectivity
# ============================================================================

"""
	jacobian_pattern(def::CoreDefinition) -> SparseMatrixCSC{Bool,Int}
	jacobian_pattern(layout, topology; closure = nothing, reconstruction = false)

The structural pattern of the residual's Jacobian on a 1D mesh, from the
connectivity alone; the definition supplies the closure and whether inventory
fluxes relative to the mesh (`RelativeFluxKernel`) are present.

- A cell row reaches the unknowns of cells `i − 1 … i + 1` and, on a moving
  mesh, the positions of nodes `i − 1 … i + 2` (conduction reads the
  neighbouring cell centers, the concentrations the neighbours' volumes);
  with relative fluxes, whose face states are reconstructed from two cells
  upwind, cells `i − 2 … i + 2`, the positions of nodes `i − 2 … i + 3` and
  the velocities of its own two nodes. It reads the intensities of its two
  faces and, under a variable cross section, the area.
- A face-intensity row reaches the intensities of nodes `f − 1 … f + 1`, the
  unknowns of the cells on either side of its face (the attenuation of the
  cell it leaves, the absorptivity of the exposed cell), on a moving mesh the
  positions of those nodes, and the area.
- A node-position row reaches its node's velocity. A velocity row reaches
  its node's velocity and position; a balance row also the velocities and
  positions of its cell's other node, the masses and energy of its cell
  (through the strain rate), the area and the lateral rate.
- Under a variable cross section the area row reaches every cell's history
  pair and the area, and the lateral-rate row everything the history rows
  reach, the history pair, the area and the lateral rate: the two dense rows
  of [`dense_rows`](@ref).

Every row holds its diagonal. The pattern is fixed for a solve: dead rows
keep their entries and take zero values.
"""
function jacobian_pattern(def::CoreDefinition)
    reconstruction = find_kernel(RelativeFluxKernel, def.kernels) !== nothing
    return jacobian_pattern(
        def.layout,
        def.setup.topology;
        closure = def.closure,
        reconstruction,
    )
end

function jacobian_pattern(
    layout::StateLayout,
    topo::MeshTopology{1};
    closure = nothing,
    reconstruction::Bool = false,
)
    N = n_cells(layout)
    L = state_length(layout)
    moving, bl, chi = is_moving(layout), has_intensities(layout), has_chi(layout)
    reach = reconstruction ? 2 : 1
    rows = Int[]
    cols = Int[]
    add!(r, c) = (push!(rows, r); push!(cols, c))
    for i = 1:N
        cells = max(1, i-reach):min(N, i+reach)
        for r in cell_rows(layout, i)
            for k in cells, c in cell_rows(layout, k)
                add!(r, c)
            end
            if moving
                for nd = max(1, i-reach):min(N+1, i+1+reach)
                    add!(r, position_index(layout, nd))
                end
                if reconstruction
                    add!(r, velocity_index(layout, i))
                    add!(r, velocity_index(layout, i + 1))
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
    if moving
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
            if chi
                add!(r, area_index(layout))
                add!(r, lateral_rate_index(layout))
            end
        end
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

# The cell whose balance each node's velocity row carries (zero for the pinned
# node, and for every node under a prescribed velocity).
_balance_cells(::Nothing, N) = zeros(Int, N + 1)
function _balance_cells(cl::BalanceRowClosure, N)
    out = zeros(Int, N + 1)
    cl.override === nothing || return out
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
