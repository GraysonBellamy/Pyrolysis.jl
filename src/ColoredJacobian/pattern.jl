# ============================================================================
# Sparsity pattern from connectivity
# ============================================================================

"""
	jacobian_pattern(layout, topology) -> SparseMatrixCSC{Bool,Int}

The structural pattern of the residual's Jacobian on a 1D mesh, from the
connectivity alone. A cell row reaches the unknowns of cells `i − 1 … i + 1`
and, on a moving mesh, of nodes `i − 1 … i + 2` (conduction reads the
neighbouring cell centers, and the concentrations the neighbours' volumes);
on a fixed mesh it reaches only its own two nodes, whose face intensities
its absorbed radiation reads. A face-intensity row reaches the intensities
of nodes `f − 1 … f + 1` and the unknowns of the cells on either side of its
face (the attenuation of the cell it leaves, the absorptivity of the exposed
cell). A node-position row reaches its node's velocity; every row holds its
diagonal. The pattern is fixed for a solve: dead rows keep their entries and
take zero values.
"""
function jacobian_pattern(layout::StateLayout, topo::MeshTopology{1})
    N = n_cells(layout)
    L = state_length(layout)
    rows = Int[]
    cols = Int[]
    add!(r, c) = (push!(rows, r); push!(cols, c))
    for i = 1:N
        cells = max(1, i-1):min(N, i+1)
        nodes = is_moving(layout) ? (max(1, i-1):min(N+1, i+2)) : (i:(i+1))
        for r in cell_rows(layout, i)
            for k in cells, c in cell_rows(layout, k)
                add!(r, c)
            end
            for nd in nodes, c in node_rows(layout, nd)
                add!(r, c)
            end
        end
    end
    if has_intensities(layout)
        for f = 1:(N+1)
            r = intensity_index(layout, f)
            for nd = max(1, f-1):min(N+1, f+1)
                add!(r, intensity_index(layout, nd))
            end
            for k = max(1, f-1):min(N, f), c in cell_rows(layout, k)
                add!(r, c)
            end
        end
    end
    if is_moving(layout)
        for nd = 1:(N+1)
            add!(position_index(layout, nd), velocity_index(layout, nd))
        end
    end
    for r = 1:L
        add!(r, r)
    end
    return sparse(rows, cols, trues(length(rows)), L, L, |)
end
