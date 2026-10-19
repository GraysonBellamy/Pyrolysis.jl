"""
	Deactivation

The terminal event of a cell of the finite-volume core: a cell whose condensed
mass falls below a fraction `ε_m` of its reference condensed mass is
deactivated, its inventories go to the live cell inward of it (on a fixed
mesh its gas is vented to the event instead), its rows are masked (identity
rows of the Newton matrix), and the boundary conditions move to the face of
the outermost live cell. The state length and the Jacobian's
pattern never change.
"""
module Deactivation

using ..Materials
using ..Properties: skeletal_volume_fraction
using ..Meshes
using ..States
using ..Assembly

include("events.jl")

export DeactivationPolicy, ConservativeTransfer, DeactivationEvent, CellDeactivation
export condensed_mass, below_threshold, top_live_cell, any_below, deactivate!, reset_cells!
export replay_events!

end # module Deactivation
