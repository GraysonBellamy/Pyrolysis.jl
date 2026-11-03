# ============================================================================
# The tables a line reads by slot
# ============================================================================

"""
	RowTable

A table of the state's rows, the absolute tolerances: its copy as set up, by
initial cell and node, and the live array, by slot. A node's position and
velocity tolerances follow the thinner of its two cells, so a node between
cells that were not neighbours takes them from `thickness`, the live
thicknesses: `scale` times the thinner, over `time_scale` for the velocity.
"""
struct RowTable
    initial::Vector{Float64}
    live::Vector{Float64}
    thickness::Vector{Float64}
    scale::Float64
    time_scale::Float64
end

"""
	SlotTables

Every table a line's solve reads by cell, node or row, in the two forms that
let cells change slots: its copy as set up, indexed by the cell of the initial
mesh, and the live array, indexed by slot. `origin[s]` is the initial cell in
slot `s`, a permutation of the initial cells that lists the dead ones last,
and [`permute_slots!`](@ref) is the one writer of the live arrays, `top` the
outermost live slot it last wrote them for. `cell`,
`node`, `rows` and `matrix` (by column) are read at their own slot; `face` and
`sides` hold per-face values of a face's owner, and `sides` of its neighbour
too, so they follow those cells; `lineage` holds the reaction references as set
up beside their initial table, the entry of `matrix` that a merge updates, and
`held` the updates that table carries; `gains` are the functions that rewrite
what a bound correction reads by slot. A node has no origin of its own: the
node below slot `s` is the node below cell `origin[s]`.
"""
struct SlotTables
    origin::Vector{Int}
    top::Base.RefValue{Int}
    cell::Vector{Pair{Vector{Float64},Vector{Float64}}}
    node::Vector{Pair{Vector{Float64},Vector{Float64}}}
    face::Vector{Pair{Vector{Float64},Vector{Float64}}}
    rows::Vector{RowTable}
    matrix::Vector{Pair{Matrix{Float64},Matrix{Float64}}}
    sides::Vector{Pair{Matrix{Float64},Matrix{Float64}}}
    lineage::Vector{Pair{Matrix{Float64},Matrix{Float64}}}
    held::Vector{Tuple{Int,Int,Float64}}
    gains::Vector{Function}
end

"""
	SlotTables(origin)

The tables of a line of `length(origin)` cells, with none registered yet.
"""
SlotTables(origin::Vector{Int}) = SlotTables(
    origin,
    Ref(length(origin)),
    Pair{Vector{Float64},Vector{Float64}}[],
    Pair{Vector{Float64},Vector{Float64}}[],
    Pair{Vector{Float64},Vector{Float64}}[],
    RowTable[],
    Pair{Matrix{Float64},Matrix{Float64}}[],
    Pair{Matrix{Float64},Matrix{Float64}}[],
    Pair{Matrix{Float64},Matrix{Float64}}[],
    Tuple{Int,Int,Float64}[],
    Function[],
)

"""
	register_tables!(t::SlotTables, kind, tables...) -> t

Register each live array of `tables` under `kind` (`:cell`, `:node`, `:face`,
`:matrix` or `:sides`), its current values becoming the copy as set up.
"""
function register_tables!(t::SlotTables, kind::Symbol, tables...)
    target =
        kind === :cell ? t.cell :
        kind === :node ? t.node :
        kind === :face ? t.face :
        kind === :matrix ? t.matrix :
        kind === :sides ? t.sides : throw(ArgumentError("no table kind :$kind"))
    for live in tables
        push!(target, copy(live) => live)
    end
    return t
end

"""
	register_lineage!(t::SlotTables, references) -> t

Register the reactions × cells table of lineage references: a matrix table
whose initial form a merge can update.
"""
function register_lineage!(t::SlotTables, references::Matrix{Float64})
    initial = copy(references)
    push!(t.matrix, initial => references)
    push!(t.lineage, copy(references) => initial)
    return t
end

"""
	register_rows!(t::SlotTables, tolerances, thickness; scale, time_scale) -> t

Register the absolute tolerances; `thickness` is the registered live table of
the cells' initial sizes, `scale` a node's position tolerance per metre of its
thinner cell and `time_scale` the time that takes it to its velocity tolerance.
"""
function register_rows!(
    t::SlotTables,
    tolerances::Vector{Float64},
    thickness::Vector{Float64};
    scale::Float64,
    time_scale::Float64,
)
    push!(t.rows, RowTable(copy(tolerances), tolerances, thickness, scale, time_scale))
    return t
end

# The live node table: the node below each live slot's cell, the far node of
# the top live cell, and the far nodes of the dead cells, in their order.
function _write_nodes!(live, initial, origin, top)
    for s = 1:top
        live[s] = initial[origin[s]]
    end
    for s = (top+1):length(live)
        live[s] = initial[s == 1 ? origin[1] : origin[s-1]+1]
    end
    return live
end

# The face above slot `s` is owned by the cell in `s`, and the bottom face by
# the cell in the first slot, which also owns the face above it.
function _write_faces!(live, initial, origin)
    live[1] = origin[1] == 1 ? initial[1] : initial[origin[1]+1]
    for s in eachindex(origin)
        live[s+1] = initial[origin[s]+1]
    end
    return live
end

function _write_sides!(live, initial, origin)
    live[1, 1] = origin[1] == 1 ? initial[1, 1] : initial[1, origin[1]+1]
    for s in eachindex(origin)
        live[1, s+1] = initial[1, origin[s]+1]
    end
    for s = 2:length(origin)
        live[2, s] = initial[2, origin[s]]
    end
    return live
end

function _write_rows!(r::RowTable, layout::StateLayout, origin, top)
    live, initial = r.live, r.initial
    for s in eachindex(origin)
        live[cell_rows(layout, s)] .= view(initial, cell_rows(layout, origin[s]))
    end
    for s = 1:(n_cells(layout)+1)
        nd = s <= top ? origin[s] : s == 1 ? origin[1] : origin[s-1] + 1
        live[node_rows(layout, s)] .= view(initial, node_rows(layout, nd))
    end
    is_moving(layout) || return r
    for s = 1:top
        (s == 1 ? 1 : origin[s-1] + 1) == origin[s] && continue
        h = s == 1 ? r.thickness[1] : min(r.thickness[s-1], r.thickness[s])
        live[position_index(layout, s)] = r.scale * h
        live[velocity_index(layout, s)] = r.scale * h / r.time_scale
    end
    return r
end

"""
	permute_slots!(t::SlotTables, layout) -> t

Write every live table of `t` from its initial form through `t.origin`, the
mask of `layout` naming the live slots: the one mechanism of an event, of
[`reset_cells!`](@ref) and of [`replay_events!`](@ref).
"""
function permute_slots!(t::SlotTables, layout::StateLayout)
    origin = t.origin
    top = top_live_cell(layout)
    t.top[] = top
    for (initial, live) in t.cell
        for s in eachindex(origin)
            live[s] = initial[origin[s]]
        end
    end
    for (initial, live) in t.node
        _write_nodes!(live, initial, origin, top)
    end
    for (initial, live) in t.face
        _write_faces!(live, initial, origin)
    end
    for (initial, live) in t.matrix
        for s in eachindex(origin)
            live[:, s] .= view(initial, :, origin[s])
        end
    end
    for (initial, live) in t.sides
        _write_sides!(live, initial, origin)
    end
    for r in t.rows
        _write_rows!(r, layout, origin, top)
    end
    for rewrite in t.gains
        rewrite()
    end
    return t
end

"""
	arrange_slots!(t::SlotTables, layout, origin, updates) -> t

The slots in the order `origin`, with the initial reference table carrying the
updates `updates` and no others, for the live slots of `layout`. Nothing is
written when `t` is already there: the outermost live slot enters the node and
row tables only once a cell has left its own slot.
"""
function arrange_slots!(
    t::SlotTables,
    layout::StateLayout,
    origin::AbstractVector{Int},
    updates::Vector{Tuple{Int,Int,Float64}},
)
    moved = any(s -> origin[s] != s, eachindex(origin))
    t.origin == origin &&
        t.held == updates &&
        (!moved || t.top[] == top_live_cell(layout)) &&
        return t
    copyto!(t.origin, origin)
    if t.held != updates
        for (set_up, initial) in t.lineage
            copyto!(initial, set_up)
            for (r, cell, value) in updates
                initial[r, cell] = value
            end
        end
        copy!(t.held, updates)
    end
    return permute_slots!(t, layout)
end
