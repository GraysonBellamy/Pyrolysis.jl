# ============================================================================
# Element death
# ============================================================================

"""
	element_death!(layout, topology, cells) -> layout

Element death on a conforming mesh of more than one dimension: the cells
`cells` keep their geometry and leave the system, their rows masked (identity
rows of the Newton matrix), as are the rows of every node no live cell keeps;
a node a dead cell shares with a live one stays live and moves with it. A
face between a live cell and a dead one becomes a boundary face of the live
cell, carrying the thermal boundary kernel's exposed condition, and the
closure balances the live cells alone. No node moves at the event and no
inventory is transferred: the dead cells' inventories leave the system. The
state length and the Jacobian's pattern never change.
"""
function element_death!(layout::StateLayout, topo::MeshTopology{D}, cells) where {D}
    D >= 2 || throw(ArgumentError("a line's cells leave through `deactivate!`"))
    for i in cells, r in cell_rows(layout, i)
        layout.mask[r] = false
    end
    kept = falses(n_nodes(topo))
    for i = 1:n_cells(topo)
        layout.mask[energy_index(layout, i)] || continue
        for nd in topo.cells[i].nodes
            kept[nd] = true
        end
    end
    for nd = 1:n_nodes(topo)
        kept[nd] && continue
        for r in node_rows(layout, nd)
            layout.mask[r] = false
        end
    end
    return layout
end
