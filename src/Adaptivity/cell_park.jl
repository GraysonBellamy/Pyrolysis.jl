"""
Parking a depleted cell.

The mesh arrays keep their length for the whole solve. Live cells occupy
`1:mesh.n_active_cells`, live nodes and faces `1:mesh.n_active_cells+1`; the
slots above are ghosts that no kernel reads. Parking a cell merges its mass
and condensed enthalpy into a neighbour and moves the freed slot into the
ghost range, so the live range stays contiguous and face
`mesh.n_active_cells + 1` is always the exposed (top) boundary face.

Under ALE the ghost nodes sit at the live surface position, so a ghost cell
has zero thickness and the live thickness is `z[n_live+1] - z[1]`. Under
Eulerian mesh motion nodes never move; only the top live cell can be parked,
and its slot keeps its geometry while the boundary condition moves down to
face `n_live + 1`.
"""

using ..Geometry: Face1D, Node1D

# Shift `v[from+1:to]` down one slot; `v[to]` is left holding the old `v[to]`
# (callers overwrite it with the ghost value).
@inline function _shift_down!(v::AbstractVector, from::Int, to::Int)
    @inbounds for i = from:(to-1)
        v[i] = v[i+1]
    end
    return v
end

# Face `f_top` becomes the exposed boundary face (cell `f_top - 1` below,
# nothing above); the previous top face `f_old` becomes an inert ghost face.
function _set_top_boundary_face!(mesh::Mesh1D, f_top::Int, f_old::Int)
    old = mesh.faces[f_old]
    mesh.faces[f_old] = Face1D(f_old, old.node_id, (f_old - 1, 0), false, 0)
    top = mesh.faces[f_top]
    mesh.faces[f_top] = Face1D(f_top, top.node_id, (f_top - 1, 0), true, 1)
    n_old = mesh.nodes[f_old]
    mesh.nodes[f_old] = Node1D(f_old, n_old.z, n_old.z_ref, false, 0)
    n_top = mesh.nodes[f_top]
    mesh.nodes[f_top] = Node1D(f_top, n_top.z, n_top.z_ref, true, 1)
    mesh.boundary_faces[1] = [f_top]
    return nothing
end

"""
	park_cell!(mesh::Mesh1D{NC}, i_park::Int, i_keep::Int;
	           material::Material{NC}, ale::Bool) -> Int

Merge live cell `i_park` into its live neighbour `i_keep` and move the freed
slot into the ghost range. Returns the new live-cell count.

Conservation: the component masses `V·ξⱼ` of both cells are summed into the
keeper (`ξ_keep = (V_p ξ_p + V_k ξ_k) / V_keep`, with `V_keep = V_p + V_k`
under ALE and the unchanged `V_k` on a fixed mesh), the temperature solves the
condensed-enthalpy balance (`merge_temperatures_conservatively`), and
`mesh.m_dry_solid_initial` is summed into the keeper. The ghost slot carries
the merged temperature, zero concentrations and — under ALE — zero volume.

Under ALE the cells above `i_park` shift down one index and the nodes above
the shared node shift down one index, so the live geometry is unchanged and
the ghost node lands on the live surface. On a fixed mesh (`ale = false`)
only the top live cell may be parked; nothing shifts.

The exposed boundary moves to face `n_live + 1` (`Face1D`/`Node1D` boundary
flags and `mesh.boundary_faces[1]`), `mesh.n_active_cells` and
`mesh.active_cell_mask` are updated, and the geometry cache is rebuilt.
"""
function park_cell!(
    mesh::Mesh1D{NC},
    i_park::Int,
    i_keep::Int;
    material::Material{NC},
    ale::Bool,
) where {NC}
    n_live = mesh.n_active_cells
    n_live >= 2 || error("park_cell!: cannot park with $(n_live) live cell(s)")
    (1 <= i_park <= n_live && 1 <= i_keep <= n_live) ||
        error("park_cell!: cells $(i_park), $(i_keep) must both be live (1:$(n_live))")
    abs(i_park - i_keep) == 1 ||
        error("park_cell!: cells $(i_park) and $(i_keep) are not adjacent")
    (ale || i_park == n_live) || error(
        "park_cell!: on a fixed mesh only the top live cell ($(n_live)) can be parked, got $(i_park)",
    )

    state_p = mesh.cell_states[i_park]
    state_k = mesh.cell_states[i_keep]
    V_p = mesh.cell_volumes[i_park]
    V_k = mesh.cell_volumes[i_keep]
    V_keep = ale ? V_p + V_k : V_k
    T_m = merge_temperatures_conservatively(V_p, state_p, V_k, state_k, material)
    ξ_m = ntuple(j -> (V_p * state_p.ξ[j] + V_k * state_k.ξ[j]) / V_keep, Val(NC))
    ξ_ghost = ntuple(_ -> 0.0, Val(NC))
    mesh.cell_states[i_keep] = CellState{NC}(
        T_m,
        ξ_m,
        V_keep,
        state_k.volume_ref + (ale ? state_p.volume_ref : 0.0),
    )
    mesh.cell_volumes[i_keep] = V_keep
    has_m_init = !isempty(mesh.m_dry_solid_initial)
    if has_m_init
        mesh.m_dry_solid_initial[i_keep] += mesh.m_dry_solid_initial[i_park]
        mesh.m_dry_solid_initial[i_park] = 0.0
    end

    n_top_old = n_live + 1          # top node / face before the park
    if ale
        q = max(i_park, i_keep)     # the node the two cells shared
        z_top = mesh.nodes[n_top_old].z
        z_ref_top = mesh.nodes[n_top_old].z_ref
        _shift_down!(mesh.cell_states, i_park, n_live)
        _shift_down!(mesh.cell_volumes, i_park, n_live)
        has_m_init && _shift_down!(mesh.m_dry_solid_initial, i_park, n_live)
        @inbounds for i = q:(n_top_old-1)
            nxt = mesh.nodes[i+1]
            cur = mesh.nodes[i]
            mesh.nodes[i] = Node1D(i, nxt.z, nxt.z_ref, cur.is_boundary, cur.boundary_id)
        end
        _shift_down!(mesh.node_velocities, q, n_top_old)
        ghost_node = mesh.nodes[n_top_old]
        mesh.nodes[n_top_old] = Node1D(
            n_top_old,
            z_top,
            z_ref_top,
            ghost_node.is_boundary,
            ghost_node.boundary_id,
        )
        mesh.node_velocities[n_top_old] = mesh.node_velocities[n_top_old-1]
    end

    ghost = n_live                  # the freed slot is the last live index
    V_ghost = ale ? 0.0 : V_p
    mesh.cell_states[ghost] = CellState{NC}(T_m, ξ_ghost, V_ghost, V_ghost)
    mesh.cell_volumes[ghost] = V_ghost
    has_m_init && (mesh.m_dry_solid_initial[ghost] = 0.0)
    mesh.active_cell_mask[ghost] = false
    mesh.n_active_cells = n_live - 1

    _set_top_boundary_face!(mesh, n_live, n_top_old)
    update_geometry!(mesh)
    assert_amr_invariants(mesh)
    return mesh.n_active_cells
end
