"""
h-coarsening (cell merging).

Reference: Design Document Section 4.4

h-coarsening merges sibling cells when error indicators are low.
Only cells that were created by a previous refinement can be coarsened.

## Conservation Property
When children a and b merge back into parent p:
  V_p · φ_p = V_a · φ_a + V_b · φ_b

This is exact: φ_p = (V_a · φ_a + V_b · φ_b) / V_p
"""

using ..Geometry: Mesh1D, Node1D, Cell1D, Face1D, CellState
using ..Geometry: cell_size, update_geometry!, n_cells
using ..Core: AbstractMaterial, Material, SOLID

# ============================================================================
# Energy-Conserving Temperature Merge
# ============================================================================

"""
	merge_temperatures_conservatively(V_a::Float64, state_a::CellState{NC},
									   V_b::Float64, state_b::CellState{NC},
									   material::Material{NC}) where NC -> Float64

Compute the merged temperature using energy conservation.

# Principle
Instead of simple volume-weighted averaging (which doesn't conserve energy when
cells have different heat capacities), we conserve total internal energy:

```
E_total = V_a × ρcₚ_a × T_a + V_b × ρcₚ_b × T_b
T_merged = E_total / (V_merged × ρcₚ_merged)
```

This is thermodynamically correct because it preserves the total thermal energy
even when cells have vastly different compositions (e.g., solid vs. gas).

# Arguments
- `V_a`: Volume of cell a [m³]
- `state_a`: Cell state for cell a (T, ξ)
- `V_b`: Volume of cell b [m³]
- `state_b`: Cell state for cell b (T, ξ)
- `material`: Material definition with component heat capacity functions

# Returns
Merged temperature [K] that conserves total internal energy

# Example Error from Non-Conservative Merge
Consider merging:
- Cell A: V=1, T=800K, ρcₚ=2000 kJ/(m³·K) (solid)
- Cell B: V=1, T=400K, ρcₚ=500 kJ/(m³·K) (mostly gas)

Volume-weighted (WRONG): T = (800 + 400) / 2 = 600 K
Energy-conserving (CORRECT): T = (1,600,000 + 200,000) / (1250 × 2) = 720 K

Error: 120 K (20% of temperature difference)
"""
function merge_temperatures_conservatively(V_a::Float64, state_a::CellState{NC},
	V_b::Float64, state_b::CellState{NC},
	material::Material{NC}) where NC
	# Compute volumetric heat capacities: ρcₚ = Σⱼ ξⱼ × cₚ,ⱼ(T)
	ρcₚ_a = 0.0
	ρcₚ_b = 0.0
	@inbounds for j in 1:NC
		c_j = material.components[j].heat_capacity
		ρcₚ_a += state_a.ξ[j] * c_j(state_a.T)
		ρcₚ_b += state_b.ξ[j] * c_j(state_b.T)
	end

	# Compute total internal energy
	E_a = V_a * ρcₚ_a * state_a.T
	E_b = V_b * ρcₚ_b * state_b.T
	E_total = E_a + E_b

	# Compute merged volume and concentrations
	V_total = V_a + V_b

	# Handle edge case where both cells have essentially zero heat capacity
	# (this shouldn't happen in practice, but be defensive)
	if ρcₚ_a < 1e-20 && ρcₚ_b < 1e-20
		# Fall back to simple volume-weighted average
		return (V_a * state_a.T + V_b * state_b.T) / V_total
	end

	# Estimate merged concentration for heat capacity calculation
	# Use simple volume-weighted average for concentrations (mass-conserving)
	ξ_merged = ntuple(NC) do j
		(V_a * state_a.ξ[j] + V_b * state_b.ξ[j]) / V_total
	end

	# Initial temperature estimate for evaluating merged heat capacity
	# Use energy-weighted estimate (better than volume-weighted for T)
	T_estimate = (E_a + E_b) / (V_a * ρcₚ_a + V_b * ρcₚ_b + 1e-20)

	# Compute merged volumetric heat capacity at estimated temperature
	ρcₚ_merged = 0.0
	@inbounds for j in 1:NC
		c_j = material.components[j].heat_capacity
		ρcₚ_merged += ξ_merged[j] * c_j(T_estimate)
	end

	# Avoid division by zero
	if ρcₚ_merged < 1e-20
		return T_estimate
	end

	# Compute final merged temperature from energy conservation
	T_merged = E_total / (V_total * ρcₚ_merged)

	return T_merged
end

"""
	merge_temperatures_conservatively_from_values(V_a::Float64, T_a::Float64, ξ_a::NTuple{NC,Float64},
												   V_b::Float64, T_b::Float64, ξ_b::NTuple{NC,Float64},
												   material::Material{NC}) where NC -> Float64

Alternate interface for energy-conserving temperature merge using raw values instead of CellState.
Useful when working directly with state vectors.
"""
function merge_temperatures_conservatively_from_values(
	V_a::Float64, T_a::Float64, ξ_a::NTuple{NC, Float64},
	V_b::Float64, T_b::Float64, ξ_b::NTuple{NC, Float64},
	material::Material{NC}) where NC

	# Create temporary cell states
	state_a = CellState{NC}(T_a, ξ_a, V_a, V_a)
	state_b = CellState{NC}(T_b, ξ_b, V_b, V_b)

	return merge_temperatures_conservatively(V_a, state_a, V_b, state_b, material)
end

# ============================================================================
# Single Cell Pair Coarsening
# ============================================================================

"""
	can_coarsen(mesh::Mesh1D, cell_id::Int) -> Bool

Check if a cell can be coarsened (merged with sibling).

A cell can be coarsened if:
1. It has a parent (is not a root cell)
2. It is a leaf (has no children)
3. Its sibling is also a leaf
"""
function can_coarsen(mesh::Mesh1D, cell_id::Int)
	cell = mesh.cells[cell_id]

	# Must have a parent
	if cell.parent == 0
		return false
	end

	# Must be a leaf
	if cell.children != (0, 0)
		return false
	end

	# Sibling must also be a leaf
	parent = mesh.cells[cell.parent]
	sibling_id = parent.children[1] == cell_id ? parent.children[2] : parent.children[1]
	sibling = mesh.cells[sibling_id]

	if sibling.children != (0, 0)
		return false
	end

	return true
end

"""
	get_sibling(mesh::Mesh1D, cell_id::Int) -> Int

Get the sibling cell ID (the other child of the same parent).
Returns 0 if no sibling exists.
"""
function get_sibling(mesh::Mesh1D, cell_id::Int)
	cell = mesh.cells[cell_id]

	if cell.parent == 0
		return 0
	end

	parent = mesh.cells[cell.parent]
	if parent.children[1] == cell_id
		return parent.children[2]
	else
		return parent.children[1]
	end
end

"""
	coarsen_cells!(mesh::Mesh1D{NC}, cell_a_id::Int, cell_b_id::Int; material=nothing) -> Int

Merge two sibling cells back into their parent.

# Arguments
- `mesh`: The mesh to modify
- `cell_a_id`: First child cell index
- `cell_b_id`: Second child cell index (sibling)
- `material=nothing`: Optional material for energy-conserving temperature merge.
  If provided, uses ρcₚ-weighted energy conservation. Otherwise, uses simple
  volume-weighted averaging.

# Returns
- Parent cell index

# Algorithm
1. Verify cells are siblings and both are leaves
2. Compute conservative average of child states
3. Update parent state
4. Restore parent face connectivity
5. Mark children as inactive
6. Mark parent as active

# Energy Conservation Note
When `material` is provided, temperature is computed via energy conservation:
  E_total = V_a × ρcₚ_a × T_a + V_b × ρcₚ_b × T_b
  T_parent = E_total / (V_parent × ρcₚ_parent)

This is critical when merging cells with different compositions (e.g., solid vs gas),
as simple volume-weighted averaging can introduce significant temperature errors.

# Note
This does NOT remove the child cells or intermediate node from the mesh arrays.
They are simply marked inactive. This avoids index invalidation.
"""
function coarsen_cells!(mesh::Mesh1D{NC}, cell_a_id::Int, cell_b_id::Int;
	material::Union{Nothing, Material{NC}} = nothing) where NC
	cell_a = mesh.cells[cell_a_id]
	cell_b = mesh.cells[cell_b_id]

	# Verify they are siblings
	if cell_a.parent != cell_b.parent || cell_a.parent == 0
		error("Cells $cell_a_id and $cell_b_id are not siblings")
	end

	# Verify both are leaves
	if cell_a.children != (0, 0) || cell_b.children != (0, 0)
		error("Cannot coarsen non-leaf cells")
	end

	parent_id = cell_a.parent
	parent = mesh.cells[parent_id]

	# Determine which is left (lower z) and right (higher z)
	if mesh.cell_centers[cell_a_id] < mesh.cell_centers[cell_b_id]
		left_child_id = cell_a_id
		right_child_id = cell_b_id
	else
		left_child_id = cell_b_id
		right_child_id = cell_a_id
	end

	left_child = mesh.cells[left_child_id]
	right_child = mesh.cells[right_child_id]

	# ========================================================================
	# 1. Get child states and volumes
	# ========================================================================
	state_a = mesh.cell_states[cell_a_id]
	state_b = mesh.cell_states[cell_b_id]

	V_a = mesh.cell_volumes[cell_a_id]
	V_b = mesh.cell_volumes[cell_b_id]
	V_parent = V_a + V_b

	# ========================================================================
	# 2. Compute conservative average
	# ========================================================================
	# Temperature: use energy-conserving merge if material is provided
	if material !== nothing
		T_parent = merge_temperatures_conservatively(V_a, state_a, V_b, state_b, material)
	else
		# Fallback: volume-weighted average (less accurate but works without material)
		T_parent = (V_a * state_a.T + V_b * state_b.T) / V_parent
	end

	# Concentrations: conservative (mass-preserving)
	# ξ_parent = (V_a * ξ_a + V_b * ξ_b) / V_parent
	ξ_parent = ntuple(NC) do j
		(V_a * state_a.ξ[j] + V_b * state_b.ξ[j]) / V_parent
	end

	# ========================================================================
	# 3. Update parent state
	# ========================================================================
	mesh.cell_states[parent_id] = CellState{NC}(
		T_parent,
		ξ_parent,
		V_parent,
		mesh.cell_states[parent_id].volume_ref,
	)
	mesh.cell_volumes[parent_id] = V_parent

	# ========================================================================
	# 4. Restore parent face connectivity
	# ========================================================================
	# The parent's faces should point to parent, not children

	# Left face of parent (was pointing to left child)
	left_face_id = parent.face_ids[1]
	left_face = mesh.faces[left_face_id]
	mesh.faces[left_face_id] = Face1D(
		left_face.id,
		left_face.node_id,
		(left_face.cell_ids[1], parent_id),  # Restore parent as right cell
		left_face.is_boundary,
		left_face.boundary_id,
	)

	# Right face of parent (was pointing to right child)
	right_face_id = parent.face_ids[2]
	right_face = mesh.faces[right_face_id]
	mesh.faces[right_face_id] = Face1D(
		right_face.id,
		right_face.node_id,
		(parent_id, right_face.cell_ids[2]),  # Restore parent as left cell
		right_face.is_boundary,
		right_face.boundary_id,
	)

	# ========================================================================
	# 5. Update parent cell to have no children
	# ========================================================================
	mesh.cells[parent_id] = Cell1D(
		parent.id,
		parent.node_ids,
		parent.face_ids,
		parent.level,
		parent.parent,
		(0, 0),  # No children anymore
	)

	# ========================================================================
	# 6. Update active cell mask
	# ========================================================================
	mesh.active_cell_mask[cell_a_id] = false
	mesh.active_cell_mask[cell_b_id] = false
	mesh.active_cell_mask[parent_id] = true

	mesh.n_active_cells -= 1  # Net loss of 1 (removed 2, added 1)

	# Note: The intermediate face and node are left in the arrays but are now orphaned.
	# A full mesh compaction would be needed to remove them, but that would invalidate
	# indices throughout the code. For now, we just leave them.

	return parent_id
end

# ============================================================================
# Batch Coarsening
# ============================================================================

"""
	mark_cells_for_coarsening(mesh::Mesh1D, η::Vector{Float64},
							  θ_coarsen::Float64, h_max::Float64)
							  -> Vector{Tuple{Int,Int}}

Mark sibling cell pairs for coarsening based on error indicators.

# Arguments
- `mesh`: The mesh
- `η`: Error indicator values
- `θ_coarsen`: Coarsening threshold (coarsen if η < θ_coarsen * max(η))
- `h_max`: Maximum allowed cell size

# Returns
- Vector of (cell_a, cell_b) tuples representing sibling pairs to coarsen
"""
function mark_cells_for_coarsening(mesh::Mesh1D, η::Vector{Float64},
	θ_coarsen::Float64, h_max::Float64)
	pairs = Tuple{Int, Int}[]
	processed = Set{Int}()

	η_max = maximum(η[i] for i in eachindex(η) if mesh.active_cell_mask[i]; init = 0.0)

	if η_max ≤ 0
		return pairs
	end

	threshold = θ_coarsen * η_max

	@inbounds for i in eachindex(η)
		if !mesh.active_cell_mask[i] || i in processed
			continue
		end

		if !can_coarsen(mesh, i)
			continue
		end

		cell = mesh.cells[i]
		sibling_id = get_sibling(mesh, i)

		if sibling_id == 0 || sibling_id in processed
			continue
		end

		# Check both cells have low error
		if η[i] < threshold && η[sibling_id] < threshold
			# Check resulting cell size
			h_combined = cell_size(mesh, i) + cell_size(mesh, sibling_id)
			if h_combined ≤ h_max
				push!(pairs, (i, sibling_id))
				push!(processed, i)
				push!(processed, sibling_id)
			end
		end
	end

	return pairs
end

"""
	perform_h_coarsening!(mesh::Mesh1D{NC}, pairs::Vector{Tuple{Int,Int}}) where NC

Perform h-coarsening on all marked sibling pairs.

# Arguments
- `mesh`: The mesh to modify
- `pairs`: Vector of (cell_a, cell_b) tuples to coarsen

# Notes
- After coarsening, call `update_geometry!(mesh)` to refresh caches
"""
function perform_h_coarsening!(mesh::Mesh1D{NC}, pairs::Vector{Tuple{Int, Int}}) where NC
	if isempty(pairs)
		return nothing
	end

	for (cell_a_id, cell_b_id) in pairs
		# Double-check both are still active (in case of complex scenarios)
		if mesh.active_cell_mask[cell_a_id] && mesh.active_cell_mask[cell_b_id]
			coarsen_cells!(mesh, cell_a_id, cell_b_id)
		end
	end

	# Update geometry after all coarsening
	update_geometry!(mesh)

	return nothing
end

# ============================================================================
# Conservation Verification  
# ============================================================================

"""
	verify_coarsening_conservation(V_a, state_a, V_b, state_b, V_p, state_p, NC) -> NamedTuple

Verify that coarsening conserved mass and energy.

# Arguments
- Pre-coarsening: V_a, state_a, V_b, state_b (children)
- Post-coarsening: V_p, state_p (parent)
- NC: Number of components

# Returns
- NamedTuple with fields: mass_error, energy_error
"""
function verify_coarsening_conservation(
	V_a::Float64, state_a::CellState{NC},
	V_b::Float64, state_b::CellState{NC},
	V_p::Float64, state_p::CellState{NC}) where NC

	# Mass conservation for each component
	max_mass_error = 0.0
	for j in 1:NC
		mass_children = V_a * state_a.ξ[j] + V_b * state_b.ξ[j]
		mass_parent = V_p * state_p.ξ[j]
		error = abs(mass_parent - mass_children) / max(mass_children, 1e-10)
		max_mass_error = max(max_mass_error, error)
	end

	# Energy conservation (using T*V as proxy)
	E_children = V_a * state_a.T + V_b * state_b.T
	E_parent = V_p * state_p.T
	energy_error = abs(E_parent - E_children) / max(E_children, 1e-10)

	return (mass_error = max_mass_error, energy_error = energy_error)
end

# ============================================================================
# Adjacent Cell Merging (Extended h-coarsening for surface depletion)
# ============================================================================

"""
	find_shared_face(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int) -> Int

Find the face shared between two adjacent cells.

# Returns
- Face ID of shared face, or 0 if cells are not adjacent
"""
function find_shared_face(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int)
	cell_a = mesh.cells[cell_a_id]
	cell_b = mesh.cells[cell_b_id]

	# Check each face of cell_a
	for face_id in cell_a.face_ids
		face = mesh.faces[face_id]
		# Check if this face connects to cell_b
		if face.cell_ids[1] == cell_b_id || face.cell_ids[2] == cell_b_id
			return face_id
		end
	end

	return 0  # Not adjacent
end

"""
	are_cells_adjacent(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int) -> Bool

Check if two cells share a face.
"""
function are_cells_adjacent(mesh::Mesh1D, cell_a_id::Int, cell_b_id::Int)
	return find_shared_face(mesh, cell_a_id, cell_b_id) != 0
end

"""
	merge_adjacent_cells!(mesh::Mesh1D{NC}, cell_to_remove_id::Int, 
						  cell_to_keep_id::Int; material=nothing) where NC -> Int

Merge two adjacent cells, removing one and expanding the other.

This is an extended h-coarsening that works on ANY adjacent pair,
not just sibling cells from refinement. The key difference from
`coarsen_cells!` is that we don't require a parent relationship.

# Arguments
- `mesh`: The mesh to modify
- `cell_to_remove_id`: Cell to remove (typically depleted surface cell)
- `cell_to_keep_id`: Cell to expand (will absorb the removed cell's volume/mass)
- `material=nothing`: Optional material for energy-conserving temperature merge.
  If provided, uses ρcₚ-weighted energy conservation. Otherwise, uses simple
  volume-weighted averaging.

# Returns
- ID of the merged cell (same as cell_to_keep_id)

# Algorithm
1. Verify cells are adjacent
2. Conservative merge: combine states with mass/energy conservation
3. Expand the kept cell's geometry to encompass both
4. Update face connectivity
5. Mark removed cell as inactive

# Physical Use Case
During pyrolysis, surface cells become depleted as solid decomposes.
When the surface cell is mostly gas, we merge it with the interior
neighbor to remove the singularity in dT/dt = Q/(ρcₚ).

# Conservation
- Mass: ξ_merged = (V_a * ξ_a + V_b * ξ_b) / V_merged
- Energy: When `material` is provided, uses ρcₚ-weighted energy conservation.
  Otherwise, uses simple volume-weighted temperature averaging.

# Note
Unlike `coarsen_cells!`, this function:
- Does NOT require a parent-child relationship
- Can work on any adjacent cell pair
- Is designed for removing depleted cells, not low-error coarsening
"""
function merge_adjacent_cells!(mesh::Mesh1D{NC}, cell_to_remove_id::Int,
	cell_to_keep_id::Int; material::Union{Nothing, Material{NC}} = nothing) where NC

	# ========================================================================
	# 1. Validation
	# ========================================================================
	if !mesh.active_cell_mask[cell_to_remove_id]
		error("Cell $cell_to_remove_id is not active")
	end
	if !mesh.active_cell_mask[cell_to_keep_id]
		error("Cell $cell_to_keep_id is not active")
	end

	shared_face_id = find_shared_face(mesh, cell_to_remove_id, cell_to_keep_id)
	if shared_face_id == 0
		error("Cells $cell_to_remove_id and $cell_to_keep_id are not adjacent")
	end

	cell_remove = mesh.cells[cell_to_remove_id]
	cell_keep = mesh.cells[cell_to_keep_id]
	shared_face = mesh.faces[shared_face_id]

	# ========================================================================
	# 2. Determine spatial arrangement
	# ========================================================================
	# In 1D: left = lower z, right = higher z
	z_remove = mesh.cell_centers[cell_to_remove_id]
	z_keep = mesh.cell_centers[cell_to_keep_id]

	if z_remove > z_keep
		# Removing the right cell (surface cell)
		# Keep cell will expand to the right
		is_removing_right = true
	else
		# Removing the left cell (unusual but handle it)
		# Keep cell will expand to the left
		is_removing_right = false
	end

	# ========================================================================
	# 3. Conservative state merging
	# ========================================================================
	state_remove = mesh.cell_states[cell_to_remove_id]
	state_keep = mesh.cell_states[cell_to_keep_id]

	V_remove = mesh.cell_volumes[cell_to_remove_id]
	V_keep = mesh.cell_volumes[cell_to_keep_id]
	V_merged = V_remove + V_keep

	# Temperature: use energy-conserving merge if material is provided
	if material !== nothing
		T_merged = merge_temperatures_conservatively(V_remove, state_remove, V_keep, state_keep, material)
	else
		# Fallback: volume-weighted average (less accurate but works without material)
		T_merged = (V_remove * state_remove.T + V_keep * state_keep.T) / V_merged
	end

	# Concentrations: conservative (mass-preserving)
	ξ_merged = ntuple(NC) do j
		(V_remove * state_remove.ξ[j] + V_keep * state_keep.ξ[j]) / V_merged
	end

	# Update the kept cell's state
	mesh.cell_states[cell_to_keep_id] = CellState{NC}(
		T_merged,
		ξ_merged,
		V_merged,
		mesh.cell_states[cell_to_keep_id].volume_ref + mesh.cell_states[cell_to_remove_id].volume_ref,
	)
	mesh.cell_volumes[cell_to_keep_id] = V_merged

	# ========================================================================
	# 4. Update node positions (expand kept cell)
	# ========================================================================
	if is_removing_right
		# The kept cell expands to take over the right node of removed cell
		new_right_node_id = cell_remove.node_ids[2]  # Right node of removed cell
		left_node_id = cell_keep.node_ids[1]  # Left node stays the same

		# Debug: check for node corruption
		if left_node_id == new_right_node_id
			@error "Node corruption in merge! is_removing_right=true" cell_to_remove_id cell_to_keep_id cell_remove.node_ids cell_keep.node_ids left_node_id new_right_node_id
		end

		# Update kept cell's node_ids
		mesh.cells[cell_to_keep_id] = Cell1D(
			cell_keep.id,
			(left_node_id, new_right_node_id),
			cell_keep.face_ids,  # Will update faces below
			cell_keep.level,
			cell_keep.parent,
			cell_keep.children,
		)
	else
		# The kept cell expands to take over the left node of removed cell
		new_left_node_id = cell_remove.node_ids[1]  # Left node of removed cell
		right_node_id = cell_keep.node_ids[2]  # Right node stays the same

		# Debug: check for node corruption
		if new_left_node_id == right_node_id
			@error "Node corruption in merge! is_removing_right=false" cell_to_remove_id cell_to_keep_id cell_remove.node_ids cell_keep.node_ids new_left_node_id right_node_id
		end

		# Update kept cell's node_ids
		mesh.cells[cell_to_keep_id] = Cell1D(
			cell_keep.id,
			(new_left_node_id, right_node_id),
			cell_keep.face_ids,  # Will update faces below
			cell_keep.level,
			cell_keep.parent,
			cell_keep.children,
		)
	end

	# ========================================================================
	# 5. Update face connectivity
	# ========================================================================
	if is_removing_right
		# Find the boundary face of the removed cell (its right face)
		# This becomes the new right face of the kept cell
		right_face_id_remove = cell_remove.face_ids[2]
		right_face = mesh.faces[right_face_id_remove]

		# Update the boundary face to point to the kept cell
		mesh.faces[right_face_id_remove] = Face1D(
			right_face.id,
			right_face.node_id,
			(cell_to_keep_id, right_face.cell_ids[2]),  # cell_ids[2] is boundary (0)
			right_face.is_boundary,
			right_face.boundary_id,
		)

		# Update kept cell's face_ids
		mesh.cells[cell_to_keep_id] = Cell1D(
			mesh.cells[cell_to_keep_id].id,
			mesh.cells[cell_to_keep_id].node_ids,
			(cell_keep.face_ids[1], right_face_id_remove),  # New right face
			mesh.cells[cell_to_keep_id].level,
			mesh.cells[cell_to_keep_id].parent,
			mesh.cells[cell_to_keep_id].children,
		)
	else
		# Find the boundary face of the removed cell (its left face)
		left_face_id_remove = cell_remove.face_ids[1]
		left_face = mesh.faces[left_face_id_remove]

		# Update the boundary face to point to the kept cell
		mesh.faces[left_face_id_remove] = Face1D(
			left_face.id,
			left_face.node_id,
			(left_face.cell_ids[1], cell_to_keep_id),  # cell_ids[1] is boundary (0)
			left_face.is_boundary,
			left_face.boundary_id,
		)

		# Update kept cell's face_ids
		mesh.cells[cell_to_keep_id] = Cell1D(
			mesh.cells[cell_to_keep_id].id,
			mesh.cells[cell_to_keep_id].node_ids,
			(left_face_id_remove, cell_keep.face_ids[2]),  # New left face
			mesh.cells[cell_to_keep_id].level,
			mesh.cells[cell_to_keep_id].parent,
			mesh.cells[cell_to_keep_id].children,
		)
	end

	# ========================================================================
	# 6. Update boundary face dictionaries
	# ========================================================================
	# If the removed cell had a boundary face, update the dictionary
	for (boundary_id, face_ids) in mesh.boundary_faces
		for (idx, face_id) in enumerate(face_ids)
			# The face connectivity was already updated above
			# No need to change the boundary_faces dict entries
		end
	end

	# ========================================================================
	# 7. Mark removed cell as inactive
	# ========================================================================
	mesh.active_cell_mask[cell_to_remove_id] = false
	mesh.n_active_cells -= 1

	# ========================================================================
	# 8. Update geometry cache
	# ========================================================================
	# Update cell center of kept cell
	node_left = mesh.nodes[mesh.cells[cell_to_keep_id].node_ids[1]]
	node_right = mesh.nodes[mesh.cells[cell_to_keep_id].node_ids[2]]
	mesh.cell_centers[cell_to_keep_id] = (node_left.z + node_right.z) / 2

	return cell_to_keep_id
end

"""
	merge_surface_cell!(mesh::Mesh1D{NC}; material=nothing) where NC -> Int

Convenience function to merge the depleted surface cell with its interior neighbor.

# Arguments
- `mesh`: The mesh to modify
- `material=nothing`: Optional material for energy-conserving temperature merge.
  If provided, uses ρcₚ-weighted energy conservation. Otherwise, uses simple
  volume-weighted averaging.

# Returns
- ID of the merged cell, or 0 if no merge was possible

# Usage
```julia
using Pyrolysis.Adaptivity

# When surface depletion is detected
if is_surface_depleted(mesh, material, 0.01)
	merged_id = merge_surface_cell!(mesh; material=material)
	if merged_id > 0
		# Update state vector etc.
	end
end
```
"""
function merge_surface_cell!(mesh::Mesh1D{NC}; material::Union{Nothing, Material{NC}} = nothing) where NC
	# Import from surface_depletion.jl (same module)
	surface_id = find_surface_cell(mesh)
	if surface_id == 0
		return 0
	end

	neighbor_id = find_interior_neighbor(mesh, surface_id)
	if neighbor_id == 0
		return 0  # Surface cell has no interior neighbor (only 1 cell left)
	end

	return merge_adjacent_cells!(mesh, surface_id, neighbor_id; material = material)
end

"""
	compact_mesh_after_surface_merge!(mesh::Mesh1D{NC}) where NC

Remove the inactive surface cell and orphaned surface node/face from mesh arrays.

This function should be called AFTER `merge_surface_cell!` to actually remove
the inactive cell from the arrays. This is necessary because:
1. `merge_surface_cell!` only marks the surface cell inactive
2. The ODE state vector actually removes elements  
3. The mesh arrays must match the state vector length

# Algorithm
After surface merge, the mesh has:
- Surface cell marked inactive (last active becomes new surface)
- Surface node orphaned (not connected to any active cell)
- Interior face between merged cells orphaned

This function:
1. Identifies orphaned elements (cell, node, face)
2. Removes them from arrays
3. Renumbers all indices consistently
4. Updates connectivity
5. Updates geometry caches

# Returns
Nothing. Modifies mesh in place.

# Note
This is specifically for surface cell removal where the surface cell is always
the highest-indexed active cell. For general compaction with arbitrary inactive
cells, a more complex algorithm would be needed.
"""
function compact_mesh_after_surface_merge!(mesh::Mesh1D{NC}) where NC
	n_cells_old = length(mesh.cells)
	n_nodes_old = length(mesh.nodes)
	n_faces_old = length(mesh.faces)

	if n_cells_old < 2
		error("Cannot compact mesh with fewer than 2 cells")
	end

	# The surface cell (last one) should be inactive
	if mesh.active_cell_mask[n_cells_old]
		error("Last cell is active - did you call merge_surface_cell! first?")
	end

	# ========================================================================
	# Step 0: Identify what needs to be removed
	# ========================================================================
	removed_cell_id = n_cells_old
	removed_cell = mesh.cells[removed_cell_id]

	# The new surface cell is the one that absorbed the removed cell
	new_surface_cell_id = n_cells_old - 1

	# After merge_adjacent_cells!, the kept cell's node_ids span both old cells.
	# The orphaned node is the one between the two merged cells - this is the
	# node that belonged to both cells before merging, which is:
	#   - Right node of kept cell (new_surface_cell) BEFORE merge = left node of removed cell BEFORE merge
	#
	# From merge_adjacent_cells!, when merging surface cell (cell_remove) with its neighbor (cell_keep):
	#   cell_keep's new node_ids = (cell_keep.node_ids[1], cell_remove.node_ids[2])
	# So the orphaned node is cell_keep.node_ids[2] from BEFORE the merge, which equals
	# cell_remove.node_ids[1] from BEFORE the merge.
	#
	# BUT merge_adjacent_cells! already updated cell_keep's nodes, so we need to figure
	# out the orphaned node differently.

	# The orphaned node is the one NOT used by any active cell
	used_nodes = Set{Int}()
	for i in 1:n_cells_old
		if mesh.active_cell_mask[i]
			push!(used_nodes, mesh.cells[i].node_ids[1])
			push!(used_nodes, mesh.cells[i].node_ids[2])
		end
	end

	orphaned_node_id = 0
	for i in 1:n_nodes_old
		if !(i in used_nodes)
			orphaned_node_id = i
			break
		end
	end

	if orphaned_node_id == 0
		# This shouldn't happen - there should be exactly one orphaned node
		@warn "No orphaned node found after surface merge - this is unexpected"
		orphaned_node_id = n_nodes_old  # Fallback: assume surface node
	end

	# The orphaned face is the interior face that connected the two merged cells
	# This is the face that has the removed cell in its cell_ids
	orphaned_face_id = 0
	for i in 1:n_faces_old
		face = mesh.faces[i]
		if removed_cell_id in face.cell_ids && !face.is_boundary
			orphaned_face_id = i
			break
		end
	end

	# The boundary face needs to be kept but its cell_id updated
	# It's the right face of the removed surface cell
	boundary_face_id = removed_cell.face_ids[2]

	# ========================================================================
	# Step 1: Build index remapping
	# ========================================================================
	# We'll use a swap-and-pop approach: swap the removed element with the last,
	# then pop. This requires tracking how indices change.

	# Cell remapping: cell removed_cell_id is removed
	# Since removed_cell_id == n_cells_old, it's already last, just pop
	cell_remap = collect(1:n_cells_old)
	cell_remap[removed_cell_id] = 0  # Indicates removed

	# Node remapping
	node_remap = collect(1:n_nodes_old)
	if orphaned_node_id == n_nodes_old
		# Orphaned is last, just pop
		node_remap[orphaned_node_id] = 0
	else
		# Swap orphaned with last, then pop
		node_remap[orphaned_node_id] = 0  # Removed
		node_remap[n_nodes_old] = orphaned_node_id  # Last moves to orphaned's slot
	end

	# Face remapping
	face_remap = collect(1:n_faces_old)
	if orphaned_face_id == 0
		# No orphaned face found - this can happen if merge was at boundary
		@warn "No orphaned face found - skipping face removal"
	elseif orphaned_face_id == n_faces_old
		face_remap[orphaned_face_id] = 0
	else
		face_remap[orphaned_face_id] = 0
		face_remap[n_faces_old] = orphaned_face_id
	end

	# ========================================================================
	# Step 2: Update the boundary face to point to new surface cell
	# ========================================================================
	if boundary_face_id > 0 && boundary_face_id <= n_faces_old
		boundary_face = mesh.faces[boundary_face_id]
		mesh.faces[boundary_face_id] = Face1D(
			boundary_face.id,
			boundary_face.node_id,
			(new_surface_cell_id, 0),  # (left_cell=new_surface, right_cell=boundary)
			boundary_face.is_boundary,
			boundary_face.boundary_id,
		)
	end

	# ========================================================================
	# Step 3: Remove the cell (last cell)
	# ========================================================================
	pop!(mesh.cells)
	pop!(mesh.cell_states)
	pop!(mesh.cell_centers)
	pop!(mesh.cell_volumes)
	pop!(mesh.active_cell_mask)

	# ========================================================================
	# Step 4: Remove the orphaned face
	# ========================================================================
	if orphaned_face_id > 0
		if orphaned_face_id == n_faces_old
			pop!(mesh.faces)
			pop!(mesh.face_positions)
			pop!(mesh.face_areas)
		else
			# Swap with last then pop
			mesh.faces[orphaned_face_id] = mesh.faces[n_faces_old]
			mesh.face_positions[orphaned_face_id] = mesh.face_positions[n_faces_old]
			mesh.face_areas[orphaned_face_id] = mesh.face_areas[n_faces_old]

			# Update the moved face's ID
			moved_face = mesh.faces[orphaned_face_id]
			mesh.faces[orphaned_face_id] = Face1D(
				orphaned_face_id,
				moved_face.node_id,
				moved_face.cell_ids,
				moved_face.is_boundary,
				moved_face.boundary_id,
			)

			pop!(mesh.faces)
			pop!(mesh.face_positions)
			pop!(mesh.face_areas)
		end

		n_faces_new = length(mesh.faces)

		# Update cell face references
		for (i, cell) in enumerate(mesh.cells)
			new_face_ids = Tuple(
				f == n_faces_old ? orphaned_face_id : (f > n_faces_new ? f : f)
				for f in cell.face_ids
			)
			# Simpler: just remap faces
			new_face_ids = Tuple(face_remap[f] for f in cell.face_ids)
			# Make sure no face is 0 (removed)
			if any(f == 0 for f in new_face_ids)
				# The cell references the orphaned face - this shouldn't happen
				# after merge_adjacent_cells! properly updated things
				@warn "Cell $i references orphaned face"
			end
			if new_face_ids != cell.face_ids
				mesh.cells[i] = Cell1D(cell.id, cell.node_ids, new_face_ids,
					cell.level, cell.parent, cell.children)
			end
		end

		# Update boundary_faces dict
		for (bid, fids) in mesh.boundary_faces
			mesh.boundary_faces[bid] = [face_remap[f] for f in fids if face_remap[f] != 0]
		end
	end

	# ========================================================================
	# Step 5: Remove the orphaned node
	# ========================================================================
	n_nodes_new = n_nodes_old - 1

	if orphaned_node_id == n_nodes_old
		# Orphaned node is last - just pop
		pop!(mesh.nodes)
		pop!(mesh.node_velocities)
	else
		# Swap with last then pop
		mesh.nodes[orphaned_node_id] = mesh.nodes[n_nodes_old]
		mesh.node_velocities[orphaned_node_id] = mesh.node_velocities[n_nodes_old]

		# Update the moved node's ID
		moved_node = mesh.nodes[orphaned_node_id]
		mesh.nodes[orphaned_node_id] = Node1D(
			orphaned_node_id,
			moved_node.z,
			moved_node.z_ref,
			moved_node.is_boundary,
			moved_node.boundary_id,
		)

		pop!(mesh.nodes)
		pop!(mesh.node_velocities)

		# Update cell node references (node n_nodes_old is now at orphaned_node_id)
		for (i, cell) in enumerate(mesh.cells)
			new_node_ids = Tuple(
				n == n_nodes_old ? orphaned_node_id : n
				for n in cell.node_ids
			)
			if new_node_ids != cell.node_ids
				mesh.cells[i] = Cell1D(cell.id, new_node_ids, cell.face_ids,
					cell.level, cell.parent, cell.children)
			end
		end

		# Update face node references
		for (i, face) in enumerate(mesh.faces)
			if face.node_id == n_nodes_old
				mesh.faces[i] = Face1D(face.id, orphaned_node_id, face.cell_ids,
					face.is_boundary, face.boundary_id)
			end
		end
	end

	# ========================================================================
	# Step 6: Update counts and verify
	# ========================================================================
	mesh.n_active_cells = sum(mesh.active_cell_mask)

	# ========================================================================
	# Step 7: Update geometry cache from current node positions
	# ========================================================================
	update_geometry!(mesh)

	# ========================================================================
	# Step 8: Verify consistency
	# ========================================================================
	n_cells_final = length(mesh.cells)
	n_nodes_final = length(mesh.nodes)
	n_faces_final = length(mesh.faces)

	@assert n_nodes_final == n_cells_final + 1 "Node count mismatch: expected $(n_cells_final + 1) nodes, got $n_nodes_final"
	@assert n_faces_final == n_cells_final + 1 "Face count mismatch: expected $(n_cells_final + 1) faces, got $n_faces_final"
	@assert all(mesh.active_cell_mask) "Inactive cells remain after compaction"

	return nothing
end

"""
	compact_mesh_after_merge!(mesh::Mesh1D{NC}, removed_cell_id::Int) where NC

Remove an inactive cell and its orphaned node/face from mesh arrays.
Works for any cell position (surface, interior, or bottom-adjacent).

This is the generalized version of `compact_mesh_after_surface_merge!`.

# Algorithm
After merge_adjacent_cells!, the mesh has:
- One cell marked inactive (the removed cell)
- One orphaned node (not connected to any active cell)
- One orphaned interior face (was between the merged cells)

This function:
1. Identifies orphaned elements
2. Removes elements using SHIFT (O(N) per element, maintains ordering)
3. Renumbers all indices consistently (simple decrement for elements after removed)
4. Updates connectivity and geometry

Note: We use shift instead of swap-and-pop to maintain the invariant that
cell/node/face indices correspond to spatial ordering (lower index = lower z).
This is important for 1D mesh operations that rely on left/right semantics.

# Arguments
- `mesh`: Mesh to compact (modified in place)
- `removed_cell_id`: ID of the cell that was merged away (should be inactive)

# Returns
Nothing. Modifies mesh in place.
"""
function compact_mesh_after_merge!(mesh::Mesh1D{NC}, removed_cell_id::Int) where NC
	n_cells_old = length(mesh.cells)
	n_nodes_old = length(mesh.nodes)
	n_faces_old = length(mesh.faces)

	if n_cells_old < 2
		error("Cannot compact mesh with fewer than 2 cells")
	end

	if removed_cell_id < 1 || removed_cell_id > n_cells_old
		error("Invalid removed_cell_id: $removed_cell_id (mesh has $n_cells_old cells)")
	end

	if mesh.active_cell_mask[removed_cell_id]
		error("Cell $removed_cell_id is still active - did you call merge_adjacent_cells! first?")
	end

	# ========================================================================
	# Step 1: Identify orphaned node (not used by any active cell)
	# ========================================================================
	used_nodes = Set{Int}()
	for i in 1:n_cells_old
		if mesh.active_cell_mask[i]
			push!(used_nodes, mesh.cells[i].node_ids[1])
			push!(used_nodes, mesh.cells[i].node_ids[2])
		end
	end

	orphaned_node_id = 0
	for i in 1:n_nodes_old
		if !(i in used_nodes)
			orphaned_node_id = i
			break
		end
	end

	if orphaned_node_id == 0
		@warn "No orphaned node found after merge - this may indicate mesh corruption" n_cells_old n_nodes_old
	end

	# ========================================================================
	# Step 2: Identify orphaned interior face (connected to removed cell, not boundary)
	# ========================================================================
	orphaned_face_id = 0
	for i in 1:n_faces_old
		face = mesh.faces[i]
		if removed_cell_id in face.cell_ids && !face.is_boundary
			orphaned_face_id = i
			break
		end
	end

	# ========================================================================
	# Step 3: Build index remappings for SHIFT-based removal
	# With shift: elements after removed index get decremented by 1
	# ========================================================================

	# Cell remapping: old_id -> new_id (0 = removed)
	cell_remap = collect(1:n_cells_old)
	cell_remap[removed_cell_id] = 0  # Removed cell
	for i in (removed_cell_id+1):n_cells_old
		cell_remap[i] = i - 1  # Shift down by 1
	end

	# Node remapping
	node_remap = collect(1:n_nodes_old)
	if orphaned_node_id > 0
		node_remap[orphaned_node_id] = 0  # Removed node
		for i in (orphaned_node_id+1):n_nodes_old
			node_remap[i] = i - 1  # Shift down by 1
		end
	end

	# Face remapping
	face_remap = collect(1:n_faces_old)
	if orphaned_face_id > 0
		face_remap[orphaned_face_id] = 0  # Removed face
		for i in (orphaned_face_id+1):n_faces_old
			face_remap[i] = i - 1  # Shift down by 1
		end
	end

	# ========================================================================
	# Step 4: Remove cell using deleteat! (shift)
	# ========================================================================
	deleteat!(mesh.cells, removed_cell_id)
	deleteat!(mesh.cell_states, removed_cell_id)
	deleteat!(mesh.cell_centers, removed_cell_id)
	deleteat!(mesh.cell_volumes, removed_cell_id)
	deleteat!(mesh.active_cell_mask, removed_cell_id)

	# ========================================================================
	# Step 5: Remove orphaned face using deleteat! (shift)
	# ========================================================================
	if orphaned_face_id > 0
		deleteat!(mesh.faces, orphaned_face_id)
		deleteat!(mesh.face_positions, orphaned_face_id)
		deleteat!(mesh.face_areas, orphaned_face_id)
	end

	# ========================================================================
	# Step 6: Remove orphaned node using deleteat! (shift)
	# ========================================================================
	if orphaned_node_id > 0
		deleteat!(mesh.nodes, orphaned_node_id)
		deleteat!(mesh.node_velocities, orphaned_node_id)
	end

	# ========================================================================
	# Step 7: Reindex all cell references (nodes, faces, id)
	# ========================================================================
	for i in eachindex(mesh.cells)
		cell = mesh.cells[i]
		# Apply remap: if remap[x] > 0, use remap[x]; if remap[x] == 0, ERROR
		new_node_ids = Tuple(begin
			mapped = node_remap[n]
			if mapped == 0
				error("Cell $i references removed node $n")
			end
			mapped
		end for n in cell.node_ids)

		new_face_ids = Tuple(begin
			mapped = face_remap[f]
			if mapped == 0
				error("Cell $i references removed face $f")
			end
			mapped
		end for f in cell.face_ids)

		new_parent = cell.parent > 0 ? cell_remap[cell.parent] : 0

		mesh.cells[i] = Cell1D(i, new_node_ids, new_face_ids,
			cell.level, new_parent, cell.children)
	end

	# ========================================================================
	# Step 8: Reindex all face references (nodes, cells, id)
	# ========================================================================
	for i in eachindex(mesh.faces)
		face = mesh.faces[i]
		new_node_id = node_remap[face.node_id]
		if new_node_id == 0
			error("Face $i references removed node $(face.node_id)")
		end

		# Handle cell_ids: 0 = boundary (keep as 0), otherwise apply remap
		new_cell_1 = face.cell_ids[1] == 0 ? 0 : cell_remap[face.cell_ids[1]]
		new_cell_2 = face.cell_ids[2] == 0 ? 0 : cell_remap[face.cell_ids[2]]

		mesh.faces[i] = Face1D(i, new_node_id, (new_cell_1, new_cell_2),
			face.is_boundary, face.boundary_id)
	end

	# ========================================================================
	# Step 9: Reindex node IDs
	# ========================================================================
	for i in eachindex(mesh.nodes)
		node = mesh.nodes[i]
		mesh.nodes[i] = Node1D(i, node.z, node.z_ref, node.is_boundary, node.boundary_id)
	end

	# ========================================================================
	# Step 10: Update boundary_faces dict
	# ========================================================================
	for (bid, fids) in mesh.boundary_faces
		mesh.boundary_faces[bid] = [
			face_remap[f]
			for f in fids if face_remap[f] != 0
		]
	end

	# ========================================================================
	# Step 11: Update counts and geometry
	# ========================================================================
	mesh.n_active_cells = sum(mesh.active_cell_mask)
	update_geometry!(mesh)

	# ========================================================================
	# Step 12: Verify consistency
	# ========================================================================
	n_cells_final = length(mesh.cells)
	n_nodes_final = length(mesh.nodes)
	n_faces_final = length(mesh.faces)

	@assert n_nodes_final == n_cells_final + 1 "Node count mismatch: expected $(n_cells_final + 1) nodes, got $n_nodes_final"
	@assert n_faces_final == n_cells_final + 1 "Face count mismatch: expected $(n_cells_final + 1) faces, got $n_faces_final"
	@assert all(mesh.active_cell_mask) "Inactive cells remain after compaction"

	return nothing
end
