"""
	Deactivation

The terminal event of a cell: a cell whose reactive inventory, the mass of
every reaction's condensed reactant, falls below a fraction `ε_m` of its
reference is deactivated, on a moving mesh wherever it sits, together with the
cells near the threshold around it, and on a fixed mesh from the exposed side.
Its inventories go to the live cell beside it (on a fixed mesh its gas is vented
to the event instead), the live cells above it move down into its slot with
every table a solve reads by slot, the outermost slots' rows are masked
(identity rows of the Newton matrix), and the boundary conditions move to the
face of the outermost live cell. The state length and the Jacobian's
pattern never change. [`element_death!`](@ref) is the building block on a
conforming mesh of more than one dimension, where a cell's live faces become
boundary faces; a solve does not reach it.
"""
module Deactivation

using ..Materials
using ..Properties: skeletal_volume_fraction
using ..Meshes
using ..States
using ..Assembly

include("slots.jl")
include("events.jl")
include("element_death.jl")

export SlotTables, RowTable, register_tables!, register_lineage!, register_rows!
export permute_slots!, arrange_slots!
export DeactivationEvent, CellDeactivation
export condensed_mass, reactive_mass
export below_threshold, top_live_cell, any_below, threshold_margin
export deactivate!, reset_cells!, VOID_VOLUME
export replay_events!
export element_death!

end # module Deactivation
