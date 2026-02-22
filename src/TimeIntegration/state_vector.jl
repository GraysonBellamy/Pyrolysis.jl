"""
State vector management for ODE integration.

Reference: Design Document Section 7.1

The state vector packs temperature and mass concentrations from mesh cells
into a flat vector format compatible with DifferentialEquations.jl solvers.

Layout: u = [T₁, T₂, ..., Tₙ, ξ₁₁, ξ₁₂, ..., ξ₁ₙ, ξ₂₁, ..., ξ_NC,n]
where n = number of cells, NC = number of components
"""

# ============================================================================
# State Vector Type
# ============================================================================

"""
	StateVector{NC}

Wrapper providing structured access to the flat ODE state vector.

# Type Parameters
- `NC`: Number of components

# Fields
- `n_cells::Int`: Number of cells in mesh
- `T_offset::Int`: Offset to temperature data (always 0)
- `ξ_offset::Int`: Offset to concentration data (= n_cells)
"""
struct StateVector{NC}
	n_cells::Int
	T_offset::Int   # Always 0
	ξ_offset::Int   # = n_cells
end

"""
	StateVector(n_cells::Int, NC::Int)

Create a StateVector descriptor for a mesh with `n_cells` cells and `NC` components.
"""
function StateVector(n_cells::Int, NC::Int)
	StateVector{NC}(n_cells, 0, n_cells)
end

"""
	StateVector(mesh::Mesh1D{NC}, material::Material{NC}) where NC

Create a StateVector descriptor from mesh and material.
The `NC` parameter is inferred from the mesh type.
"""
function StateVector(mesh::Mesh1D{NC}, material::Material{NC}) where NC
	n_cells_val = length(mesh.cells)
	StateVector{NC}(n_cells_val, 0, n_cells_val)
end

"""
	state_vector_length(n_cells::Int, NC::Int)

Compute the total length of the state vector.

Returns `n_cells * (1 + NC)` for temperature + NC concentration fields.
"""
@inline state_vector_length(n_cells::Int, NC::Int) = n_cells * (1 + NC)

"""
	state_vector_length(sv::StateVector{NC}) where NC

Compute the total length from a StateVector descriptor.
"""
@inline state_vector_length(sv::StateVector{NC}) where NC = sv.n_cells * (1 + NC)

# ============================================================================
# Create State Vector
# ============================================================================

"""
	create_state_vector(mesh::Mesh1D{NC}, ::Val{NC}) where NC -> Vector{Float64}

Create and initialize a state vector from mesh cell states.

# Arguments
- `mesh`: Mesh with populated cell_states
- `::Val{NC}`: Component count as type parameter

# Returns
A flat vector of length `n_cells * (1 + NC)` containing packed state.

The `NC` type parameter enables the compiler to specialize for the specific
number of components.
"""
function create_state_vector(mesh::Mesh1D{NC}, ::Val{NC}) where NC
	n_cells_val = length(mesh.cells)
	u = Vector{Float64}(undef, state_vector_length(n_cells_val, NC))
	sv = StateVector(n_cells_val, NC)
	pack_state!(u, mesh, sv)
	return u
end

"""
	create_state_vector(mesh::Mesh1D{NC}, NC::Int) where NC -> Vector{Float64}

Convenience method that accepts NC as an integer.
Converts to Val{NC} for dispatch.
"""
function create_state_vector(mesh::Mesh1D{NC}, NC_int::Int) where NC
	@assert NC == NC_int "Component count mismatch: mesh has $NC components but $NC_int was specified"
	create_state_vector(mesh, Val(NC))
end

"""
	create_state_vector(mesh::Mesh1D{NC}, material::Material{NC}) where NC -> Vector{Float64}

Create a state vector from mesh and material.
The `NC` parameter is inferred from the mesh and material types.
"""
function create_state_vector(mesh::Mesh1D{NC}, material::Material{NC}) where NC
	create_state_vector(mesh, Val(NC))
end

# ============================================================================
# Pack State (Mesh → Vector)
# ============================================================================

"""
	pack_state!(u::AbstractVector{Float64}, mesh::Mesh1D{NC}, sv::StateVector{NC}) where NC

Pack mesh cell states into flat ODE vector.

Layout:
- `u[1:n]` = temperatures T₁, T₂, ..., Tₙ
- `u[n+1:2n]` = ξ₁ (component 1 concentrations)
- `u[2n+1:3n]` = ξ₂ (component 2 concentrations)
- ... and so on for all NC components

The `NC` type parameter enables the compiler to specialize for the specific
number of components, eliminating type instability and allocations.
"""
function pack_state!(u::AbstractVector{Float64}, mesh::Mesh1D{NC}, sv::StateVector{NC}) where NC
	n = sv.n_cells

	@inbounds for i in 1:n
		state = mesh.cell_states[i]

		# Pack temperature
		u[i] = state.T

		# Pack concentrations
		for j in 1:NC
			u[n+(j-1)*n+i] = state.ξ[j]
		end
	end

	return u
end

"""
	pack_state!(u::AbstractVector{Float64}, mesh::Mesh1D{NC}, material::Material{NC}) where NC

Pack state using material to determine NC.
The `NC` parameter is inferred from the mesh and material types.
"""
function pack_state!(u::AbstractVector{Float64}, mesh::Mesh1D{NC}, material::Material{NC}) where NC
	sv = StateVector(mesh, material)
	pack_state!(u, mesh, sv)
end

# ============================================================================
# Unpack State (Vector → Mesh)
# ============================================================================

"""
	unpack_state!(mesh::Mesh1D{NC}, u::AbstractVector, sv::StateVector{NC}) where NC

Unpack flat ODE vector into mesh cell states.

This updates `mesh.cell_states[i].T` and `mesh.cell_states[i].ξ` for all cells.

The element type of `u` can be Float64 or a ForwardDiff.Dual type for
automatic differentiation compatibility with implicit solvers.

Type-stable implementation using Val(NC) for compile-time dispatch.
The `NC` type parameter enables the compiler to specialize for the specific
number of components, eliminating type instability and allocations.
"""
function unpack_state!(mesh::Mesh1D{NC}, u::AbstractVector, sv::StateVector{NC}) where NC
	n = sv.n_cells

	@inbounds for i in 1:n
		state = mesh.cell_states[i]

		# Unpack temperature - convert to Float64 to handle Dual types
		state.T = Float64(real(u[i]))

		# Unpack concentrations - build tuple and assign, converting to Float64
		# Use Val(NC) for type stability
		ξ_new = ntuple(Val(NC)) do j
			Float64(real(u[n+(j-1)*n+i]))
		end
		state.ξ = ξ_new
	end

	return mesh
end

"""
	unpack_state!(mesh::Mesh1D{NC}, u::AbstractVector, material::Material{NC}) where NC

Unpack state using material to determine NC.
The `NC` parameter is inferred from the mesh and material types.
"""
function unpack_state!(mesh::Mesh1D{NC}, u::AbstractVector, material::Material{NC}) where NC
	sv = StateVector(mesh, material)
	unpack_state!(mesh, u, sv)
end

# ============================================================================
# Accessors for State Vector
# ============================================================================

"""
	get_temperature(u::AbstractVector{Float64}, sv::StateVector, i::Int)

Get temperature of cell i from state vector.
"""
@inline function get_temperature(u::AbstractVector{Float64}, sv::StateVector, i::Int)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError(u, i))
	@inbounds u[i]
end

"""
	get_concentration(u::AbstractVector{Float64}, sv::StateVector{NC}, i::Int, j::Int) where NC

Get concentration of component j in cell i from state vector.
"""
@inline function get_concentration(u::AbstractVector{Float64}, sv::StateVector{NC}, i::Int, j::Int) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	@inbounds u[n+(j-1)*n+i]
end

"""
	get_concentrations(u::AbstractVector{Float64}, sv::StateVector{NC}, i::Int) where NC

Get all concentrations for cell i as an NTuple.

Type-stable implementation using Val(NC).
"""
@inline function get_concentrations(u::AbstractVector{Float64}, sv::StateVector{NC}, i::Int) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	n = sv.n_cells
	@inbounds ntuple(Val(NC)) do j
		u[n+(j-1)*n+i]
	end
end

"""
	set_temperature!(u::AbstractVector{Float64}, sv::StateVector, i::Int, T::Float64)

Set temperature of cell i in state vector.
"""
@inline function set_temperature!(u::AbstractVector{Float64}, sv::StateVector, i::Int, T::Float64)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@inbounds u[i] = T
end

"""
	set_concentration!(u::AbstractVector{Float64}, sv::StateVector{NC}, i::Int, j::Int, ξ::Float64) where NC

Set concentration of component j in cell i.
"""
@inline function set_concentration!(u::AbstractVector{Float64}, sv::StateVector{NC}, i::Int, j::Int, ξ::Float64) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	@inbounds u[n+(j-1)*n+i] = ξ
end

# ============================================================================
# Views into State Vector
# ============================================================================

"""
	temperature_view(u::AbstractVector{Float64}, sv::StateVector)

Get a view of the temperature portion of the state vector.
"""
@inline function temperature_view(u::AbstractVector{Float64}, sv::StateVector)
	@view u[1:sv.n_cells]
end

"""
	concentration_view(u::AbstractVector{Float64}, sv::StateVector{NC}, j::Int) where NC

Get a view of the concentration of component j.
"""
@inline function concentration_view(u::AbstractVector{Float64}, sv::StateVector{NC}, j::Int) where NC
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	start_idx = n + (j-1)*n + 1
	end_idx = n + j*n
	@view u[start_idx:end_idx]
end

# ============================================================================
# Utility Functions
# ============================================================================

"""
	validate_state_vector(u::AbstractVector{Float64}, sv::StateVector{NC}) where NC

Check that state vector has correct length.
"""
function validate_state_vector(u::AbstractVector{Float64}, sv::StateVector{NC}) where NC
	expected = state_vector_length(sv)
	actual = length(u)
	if actual != expected
		throw(DimensionMismatch("State vector length $actual does not match expected $expected"))
	end
	return true
end

"""
	copy_state!(dest::AbstractVector{Float64}, src::AbstractVector{Float64}, sv::StateVector)

Copy state from src to dest.
"""
function copy_state!(dest::AbstractVector{Float64}, src::AbstractVector{Float64}, sv::StateVector)
	validate_state_vector(src, sv)
	validate_state_vector(dest, sv)
	copyto!(dest, src)
	return dest
end

# ============================================================================
# ALE State Vector (Extended with Node Positions)
# ============================================================================

"""
	StateVectorALE{NC}

Extended state vector descriptor including node positions for ALE formulation.

The ALE state vector includes mesh geometry (node positions) alongside the
physics state, allowing the ODE integrator to evolve mesh motion simultaneously
with the physics.

# Layout
```
u = [T₁, ..., Tₙ, ξ₁₁, ..., ξ_NC,n, z₁, ..., z_{n+1}]
```

Where:
- T: Cell temperatures (n values)
- ξ: Mass concentrations (NC × n values, column-major by component)
- z: Node positions (n+1 values for 1D mesh with n cells)

# Type Parameters
- `NC`: Number of components

# Fields
- `n_cells::Int`: Number of cells
- `n_components::Int`: Number of material components (NC)
- `n_nodes::Int`: Number of nodes (n_cells + 1 for 1D)
- `T_offset::Int`: Offset to temperature data (always 0)
- `ξ_offset::Int`: Offset to concentration data (= n_cells)
- `z_offset::Int`: Offset to node position data (= n_cells + NC*n_cells)

# Example
For n=4 cells, NC=2 components, the state vector has 17 elements:
- u[1:4]: Temperatures T₁...T₄
- u[5:8]: Component 1 concentrations ξ₁₁...ξ₁₄
- u[9:12]: Component 2 concentrations ξ₂₁...ξ₂₄
- u[13:17]: Node positions z₁...z₅
"""
struct StateVectorALE{NC}
	n_cells::Int
	n_components::Int  # Redundant with NC but useful for runtime access
	n_nodes::Int
	T_offset::Int
	ξ_offset::Int
	z_offset::Int
end

"""
	StateVectorALE(n_cells::Int, ::Val{NC}) where NC

Create a StateVectorALE descriptor for a mesh with `n_cells` cells and `NC` components.
"""
function StateVectorALE(n_cells::Int, ::Val{NC}) where NC
	n_nodes = n_cells + 1  # 1D mesh
	T_offset = 0
	ξ_offset = n_cells
	z_offset = n_cells + NC * n_cells
	StateVectorALE{NC}(n_cells, NC, n_nodes, T_offset, ξ_offset, z_offset)
end

"""
	StateVectorALE(mesh::Mesh1D{NC}) where NC

Create a StateVectorALE descriptor from a mesh.
"""
StateVectorALE(mesh::Mesh1D{NC}) where NC = StateVectorALE(length(mesh.cells), Val(NC))

"""
	StateVectorALE(mesh::Mesh1D{NC}, material::Material{NC}) where NC

Create a StateVectorALE descriptor from mesh and material.
"""
StateVectorALE(mesh::Mesh1D{NC}, ::Material{NC}) where NC = StateVectorALE(mesh)

"""
	state_vector_length_ale(n_cells::Int, NC::Int) -> Int

Compute the total length of the ALE state vector.

Returns `n_cells * (1 + NC) + (n_cells + 1)` for temperature + NC concentration fields + node positions.
"""
@inline function state_vector_length_ale(n_cells::Int, NC::Int)
	n_cells * (1 + NC) + (n_cells + 1)
end

"""
	state_vector_length_ale(sv::StateVectorALE{NC}) where NC -> Int

Compute the total length from a StateVectorALE descriptor.
"""
@inline function state_vector_length_ale(sv::StateVectorALE{NC}) where NC
	sv.n_cells * (1 + NC) + sv.n_nodes
end

# ============================================================================
# Create ALE State Vector
# ============================================================================

"""
	create_state_vector_ale(mesh::Mesh1D{NC}) where NC -> Vector{Float64}

Create and initialize an ALE state vector from mesh, including node positions.

# Arguments
- `mesh`: Mesh with populated cell_states and nodes

# Returns
A flat vector of length `n_cells * (1 + NC) + n_nodes` containing packed state
including node positions.
"""
function create_state_vector_ale(mesh::Mesh1D{NC}) where NC
	sv = StateVectorALE(mesh)
	u = Vector{Float64}(undef, state_vector_length_ale(sv))
	pack_state_ale!(u, mesh, sv)
	return u
end

"""
	create_state_vector_ale(mesh::Mesh1D{NC}, material::Material{NC}) where NC -> Vector{Float64}

Create an ALE state vector from mesh and material.
"""
function create_state_vector_ale(mesh::Mesh1D{NC}, ::Material{NC}) where NC
	create_state_vector_ale(mesh)
end

# ============================================================================
# Pack ALE State (Mesh → Vector)
# ============================================================================

"""
	pack_state_ale!(u::AbstractVector{Float64}, mesh::Mesh1D{NC}, sv::StateVectorALE{NC}) where NC

Pack mesh cell states and node positions into flat ODE vector.

Layout:
- `u[1:n]` = temperatures T₁, T₂, ..., Tₙ
- `u[n+1:2n]` = ξ₁ (component 1 concentrations)
- `u[2n+1:3n]` = ξ₂ (component 2 concentrations)
- ... and so on for all NC components
- `u[n*(1+NC)+1:end]` = node positions z₁, z₂, ..., z_{n+1}
"""
function pack_state_ale!(u::AbstractVector{Float64}, mesh::Mesh1D{NC},
	sv::StateVectorALE{NC}) where NC
	n = sv.n_cells

	# Pack temperatures
	@inbounds for i in 1:n
		u[i] = mesh.cell_states[i].T
	end

	# Pack concentrations (column-major by component)
	@inbounds for i in 1:n
		for j in 1:NC
			u[sv.ξ_offset+(j-1)*n+i] = mesh.cell_states[i].ξ[j]
		end
	end

	# Pack node positions
	@inbounds for i in 1:sv.n_nodes
		u[sv.z_offset+i] = mesh.nodes[i].z
	end

	return u
end

"""
	pack_state_ale!(u::AbstractVector{Float64}, mesh::Mesh1D{NC}, material::Material{NC}) where NC

Pack ALE state using material to determine NC.
"""
function pack_state_ale!(u::AbstractVector{Float64}, mesh::Mesh1D{NC},
	::Material{NC}) where NC
	sv = StateVectorALE(mesh)
	pack_state_ale!(u, mesh, sv)
end

# ============================================================================
# Unpack ALE State (Vector → Mesh)
# ============================================================================

"""
	_unpack_concentrations_ale(u, ξ_offset, n, i, ::Val{NC}) where NC

Generated function to unpack concentrations from state vector without closure allocations.
Returns an NTuple{NC, Float64} with all concentrations for cell i.

PERFORMANCE: This @generated function creates the tuple at compile time without
capturing variables in a closure, avoiding the ~1MB allocation per call that
occurs with ntuple(Val(NC)) do j ... end.
"""
@generated function _unpack_concentrations_ale(u::AbstractVector, ξ_offset::Int,
	n::Int, i::Int, ::Val{NC}) where NC
	# Build tuple construction at compile time
	# Note: j-1 is computed at compile time as a literal integer
	exprs = [:(Float64(real(u[ξ_offset+$(j-1)*n+i]))) for j in 1:NC]
	return :(($(exprs...),))
end

"""
	unpack_state_ale!(mesh::Mesh1D{NC}, u::AbstractVector, sv::StateVectorALE{NC}) where NC

Unpack flat ODE vector into mesh cell states and node positions.

This updates:
- `mesh.cell_states[i].T` and `mesh.cell_states[i].ξ` for all cells
- `mesh.nodes[i].z` for all nodes
- All cached geometry (volumes, centers, face positions) via `update_geometry_from_nodes!`

The element type of `u` can be Float64 or a ForwardDiff.Dual type for
automatic differentiation compatibility with implicit solvers.

Note: During ForwardDiff AD passes, u contains Dual numbers. We only update the mesh
when u contains real Float64 values, not during AD Jacobian computation.
"""
function unpack_state_ale!(mesh::Mesh1D{NC}, u::AbstractVector{T},
	sv::StateVectorALE{NC}) where {NC, T}

	# Skip mesh updates during AD (when T is a Dual type)
	# The mesh state is only needed for Float64 evaluations
	T <: ForwardDiff.Dual && return mesh

	n = sv.n_cells
	ξ_offset = sv.ξ_offset

	# Unpack temperatures and concentrations to cell_states
	@inbounds for i in 1:n
		state = mesh.cell_states[i]

		# Unpack temperature
		state.T = Float64(real(u[i]))

		# Unpack concentrations using generated function (allocation-free)
		state.ξ = _unpack_concentrations_ale(u, ξ_offset, n, i, Val(NC))
	end

	# Unpack node positions and update geometry
	@inbounds for i in 1:sv.n_nodes
		z_new = Float64(real(u[sv.z_offset+i]))
		# Update node position in-place (Node1D is now mutable - no allocation)
		mesh.nodes[i].z = z_new
	end

	# Recompute geometry from new node positions
	update_geometry_from_nodes!(mesh)

	return mesh
end

"""
	unpack_state_ale!(mesh::Mesh1D{NC}, u::AbstractVector, material::Material{NC}) where NC

Unpack ALE state using material to determine NC.
"""
function unpack_state_ale!(mesh::Mesh1D{NC}, u::AbstractVector,
	::Material{NC}) where NC
	sv = StateVectorALE(mesh)
	unpack_state_ale!(mesh, u, sv)
end

# ============================================================================
# ALE Accessors for State Vector
# ============================================================================

"""
	get_temperature_ale(u::AbstractVector, sv::StateVectorALE, i::Int)

Get temperature of cell i from ALE state vector.
"""
@inline function get_temperature_ale(u::AbstractVector, sv::StateVectorALE, i::Int)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError(u, i))
	@inbounds u[i]
end

"""
	get_concentration_ale(u::AbstractVector, sv::StateVectorALE{NC}, i::Int, j::Int) where NC

Get concentration of component j in cell i from ALE state vector.
"""
@inline function get_concentration_ale(u::AbstractVector, sv::StateVectorALE{NC},
	i::Int, j::Int) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	@inbounds u[sv.ξ_offset+(j-1)*n+i]
end

"""
	get_concentrations_ale(u::AbstractVector, sv::StateVectorALE{NC}, i::Int) where NC

Get all concentrations for cell i as an NTuple from ALE state vector.
"""
@inline function get_concentrations_ale(u::AbstractVector, sv::StateVectorALE{NC},
	i::Int) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	n = sv.n_cells
	@inbounds ntuple(Val(NC)) do j
		u[sv.ξ_offset+(j-1)*n+i]
	end
end

"""
	get_node_position(u::AbstractVector, sv::StateVectorALE, i::Int)

Get position of node i from ALE state vector.
"""
@inline function get_node_position(u::AbstractVector, sv::StateVectorALE, i::Int)
	@boundscheck 1 <= i <= sv.n_nodes || throw(BoundsError())
	@inbounds u[sv.z_offset+i]
end

"""
	get_node_positions(u::AbstractVector, sv::StateVectorALE)

Get a view of all node positions from ALE state vector.

Returns a view (not a copy) for efficiency.
"""
@inline function get_node_positions(u::AbstractVector, sv::StateVectorALE)
	@view u[(sv.z_offset+1):(sv.z_offset+sv.n_nodes)]
end

"""
	set_temperature_ale!(u::AbstractVector, sv::StateVectorALE, i::Int, T)

Set temperature of cell i in ALE state vector.
"""
@inline function set_temperature_ale!(u::AbstractVector, sv::StateVectorALE,
	i::Int, T)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@inbounds u[i] = T
end

"""
	set_concentration_ale!(u::AbstractVector, sv::StateVectorALE{NC}, i::Int, j::Int, ξ) where NC

Set concentration of component j in cell i in ALE state vector.
"""
@inline function set_concentration_ale!(u::AbstractVector, sv::StateVectorALE{NC},
	i::Int, j::Int, ξ) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	@inbounds u[sv.ξ_offset+(j-1)*n+i] = ξ
end

"""
	set_node_position!(u::AbstractVector, sv::StateVectorALE, i::Int, z)

Set position of node i in ALE state vector.
"""
@inline function set_node_position!(u::AbstractVector, sv::StateVectorALE,
	i::Int, z)
	@boundscheck 1 <= i <= sv.n_nodes || throw(BoundsError())
	@inbounds u[sv.z_offset+i] = z
end

"""
	set_node_velocity!(du::AbstractVector, sv::StateVectorALE, i::Int, w)

Set node velocity (dz/dt) for node i in the derivative vector.

This is used when assembling the ODE RHS to set dz/dt = w.
"""
@inline function set_node_velocity!(du::AbstractVector, sv::StateVectorALE,
	i::Int, w)
	@boundscheck 1 <= i <= sv.n_nodes || throw(BoundsError())
	@inbounds du[sv.z_offset+i] = w
end

# ============================================================================
# ALE Views into State Vector
# ============================================================================

"""
	temperature_view_ale(u::AbstractVector, sv::StateVectorALE)

Get a view of the temperature portion of the ALE state vector.
"""
@inline function temperature_view_ale(u::AbstractVector, sv::StateVectorALE)
	@view u[1:sv.n_cells]
end

"""
	concentration_view_ale(u::AbstractVector, sv::StateVectorALE{NC}, j::Int) where NC

Get a view of the concentration of component j from ALE state vector.
"""
@inline function concentration_view_ale(u::AbstractVector, sv::StateVectorALE{NC},
	j::Int) where NC
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	start_idx = sv.ξ_offset + (j-1)*n + 1
	end_idx = sv.ξ_offset + j*n
	@view u[start_idx:end_idx]
end

"""
	node_positions_view(u::AbstractVector, sv::StateVectorALE)

Get a view of the node positions portion of the ALE state vector.

Alias for `get_node_positions`.
"""
@inline node_positions_view(u::AbstractVector, sv::StateVectorALE) =
	get_node_positions(u, sv)

# ============================================================================
# ALE Utility Functions
# ============================================================================

"""
	validate_state_vector_ale(u::AbstractVector, sv::StateVectorALE{NC}) where NC

Check that ALE state vector has correct length.
"""
function validate_state_vector_ale(u::AbstractVector, sv::StateVectorALE{NC}) where NC
	expected = state_vector_length_ale(sv)
	actual = length(u)
	if actual != expected
		throw(DimensionMismatch("ALE state vector length $actual does not match expected $expected"))
	end
	return true
end

"""
	get_thickness_from_state(u::AbstractVector, sv::StateVectorALE) -> Float64

Compute current domain thickness from node positions in state vector.

Returns z_top - z_bottom (should be positive).
"""
@inline function get_thickness_from_state(u::AbstractVector, sv::StateVectorALE)
	z_bottom = get_node_position(u, sv, 1)
	z_top = get_node_position(u, sv, sv.n_nodes)
	return z_top - z_bottom
end

"""
	compute_volumes_from_state(u::AbstractVector, sv::StateVectorALE, A::Float64) -> Vector{Float64}

Compute cell volumes from node positions in state vector.

# Arguments
- `u`: ALE state vector
- `sv`: State vector descriptor
- `A`: Cross-sectional area

# Returns
Vector of cell volumes (length n_cells).
"""
function compute_volumes_from_state(u::AbstractVector, sv::StateVectorALE,
	A::Float64)
	n = sv.n_cells
	volumes = Vector{Float64}(undef, n)

	@inbounds for i in 1:n
		z_left = get_node_position(u, sv, i)
		z_right = get_node_position(u, sv, i + 1)
		volumes[i] = abs(z_right - z_left) * A
	end

	return volumes
end

"""
	compute_cell_centers_from_state(u::AbstractVector, sv::StateVectorALE) -> Vector{Float64}

Compute cell centers from node positions in state vector.
"""
function compute_cell_centers_from_state(u::AbstractVector, sv::StateVectorALE)
	n = sv.n_cells
	centers = Vector{Float64}(undef, n)

	@inbounds for i in 1:n
		z_left = get_node_position(u, sv, i)
		z_right = get_node_position(u, sv, i + 1)
		centers[i] = 0.5 * (z_left + z_right)
	end

	return centers
end

# ============================================================================
# State Vector Resizing for Cell Merging
# ============================================================================

"""
	remove_surface_cell_from_state_vector(u::AbstractVector{Float64}, 
										   sv::StateVectorALE{NC}) where NC 
		-> Tuple{Vector{Float64}, StateVectorALE{NC}}

Remove the surface cell (highest z) from the ALE state vector.

This is used after merging a depleted surface cell with its interior neighbor.
The mesh topology has already been updated, now we need to resize the state vector.

# Arguments
- `u`: Current ALE state vector
- `sv`: Current state vector descriptor

# Returns
- `(u_new, sv_new)`: Resized state vector and updated descriptor

# Layout Change
Before (n cells, n+1 nodes):
  u = [T₁, ..., Tₙ, ξ₁₁, ..., ξ_NC,n, z₁, ..., z_{n+1}]

After (n-1 cells, n nodes):
  u_new = [T₁, ..., T_{n-1}, ξ₁₁, ..., ξ_NC,{n-1}, z₁, ..., zₙ]

# Notes
- The surface cell is at index n (last cell)
- The surface node is at index n+1 (last node)
- This function assumes the mesh has already been updated via merge_adjacent_cells!
- The merged state values should already be in the interior cell's position

# Example
```julia
# After merging depleted surface cell
merged_id = merge_adjacent_cells!(mesh, surface_id, neighbor_id)

# Resize state vector
u_new, sv_new = remove_surface_cell_from_state_vector(u, sv)
```
"""
function remove_surface_cell_from_state_vector(u::AbstractVector{Float64},
	sv::StateVectorALE{NC}) where NC

	n_old = sv.n_cells
	n_new = n_old - 1

	if n_new < 1
		error("Cannot remove cell: only 1 cell remaining")
	end

	# Create new state vector with reduced size
	sv_new = StateVectorALE(n_new, Val(NC))
	u_new = Vector{Float64}(undef, state_vector_length_ale(sv_new))

	# Copy temperatures (skip last cell)
	@inbounds for i in 1:n_new
		u_new[i] = u[i]
	end

	# Copy concentrations (skip last cell for each component)
	@inbounds for j in 1:NC
		for i in 1:n_new
			# Source index in old vector
			src_idx = sv.ξ_offset + (j - 1) * n_old + i
			# Destination index in new vector
			dst_idx = sv_new.ξ_offset + (j - 1) * n_new + i
			u_new[dst_idx] = u[src_idx]
		end
	end

	# Copy node positions (skip last node)
	@inbounds for i in 1:(n_new+1)  # n_new+1 = n_old nodes (minus the surface)
		# Wait - we keep n_new nodes but actually we need n_new+1 nodes for n_new cells
		# So we copy all but the last node
		u_new[sv_new.z_offset+i] = u[sv.z_offset+i]
	end

	return (u_new, sv_new)
end

"""
	remove_cell_from_state_vector(u::AbstractVector{Float64}, 
								   sv::StateVectorALE{NC},
								   cell_id::Int,
								   node_id::Int) where NC
		-> Tuple{Vector{Float64}, StateVectorALE{NC}}

Remove a specific cell and node from the ALE state vector.

This is a more general version that can remove any cell (not just surface).

# Arguments
- `u`: Current ALE state vector
- `sv`: Current state vector descriptor
- `cell_id`: Index of cell to remove (in the state vector, not mesh!)
- `node_id`: Index of node to remove (in the state vector, not mesh!)

# Returns
- `(u_new, sv_new)`: Resized state vector and updated descriptor

# Warning
This function assumes cell_id and node_id refer to positions in the STATE VECTOR,
which may differ from mesh indices if cells have been previously removed.
For simple surface cell removal, use `remove_surface_cell_from_state_vector` instead.
"""
function remove_cell_from_state_vector(u::AbstractVector{Float64},
	sv::StateVectorALE{NC},
	cell_id::Int,
	node_id::Int) where NC

	n_old = sv.n_cells
	n_nodes_old = sv.n_nodes
	n_new = n_old - 1
	n_nodes_new = n_nodes_old - 1

	if n_new < 1
		error("Cannot remove cell: only 1 cell remaining")
	end

	# Create new state vector with reduced size
	sv_new = StateVectorALE(n_new, Val(NC))
	u_new = Vector{Float64}(undef, state_vector_length_ale(sv_new))

	# Copy temperatures, skipping cell_id
	dst_idx = 0
	@inbounds for i in 1:n_old
		if i != cell_id
			dst_idx += 1
			u_new[dst_idx] = u[i]
		end
	end

	# Copy concentrations, skipping cell_id for each component
	@inbounds for j in 1:NC
		dst_idx = 0
		for i in 1:n_old
			if i != cell_id
				dst_idx += 1
				# Source index in old vector
				src_idx = sv.ξ_offset + (j - 1) * n_old + i
				# Destination index in new vector
				dst_base = sv_new.ξ_offset + (j - 1) * n_new
				u_new[dst_base+dst_idx] = u[src_idx]
			end
		end
	end

	# Copy node positions, skipping node_id
	dst_idx = 0
	@inbounds for i in 1:n_nodes_old
		if i != node_id
			dst_idx += 1
			u_new[sv_new.z_offset+dst_idx] = u[sv.z_offset+i]
		end
	end

	return (u_new, sv_new)
end

"""
	create_merged_state_vector(u::AbstractVector{Float64},
								sv::StateVectorALE{NC},
								cell_to_remove::Int,
								cell_to_keep::Int,
								V_remove::Float64,
								V_keep::Float64;
								material=nothing) where NC
		-> Tuple{Vector{Float64}, StateVectorALE{NC}}

Create a new state vector with merged cell states.

This performs both the conservative merge AND the state vector resizing in one step.
Use this when you need to merge cell states before resizing.

# Arguments
- `u`: Current ALE state vector
- `sv`: Current state vector descriptor
- `cell_to_remove`: Index of cell to remove (state vector index)
- `cell_to_keep`: Index of cell to keep (will have merged state)
- `V_remove`: Volume of cell being removed
- `V_keep`: Volume of cell being kept
- `material=nothing`: Optional material for energy-conserving temperature merge.
  If provided, uses ρcₚ-weighted energy conservation. Otherwise, uses simple
  volume-weighted averaging.

# Returns
- `(u_new, sv_new)`: New state vector with merged cell and updated descriptor

# Notes
When `material` is provided, temperature is computed via energy conservation:
  E_total = V_remove × ρcₚ_remove × T_remove + V_keep × ρcₚ_keep × T_keep
  T_merged = E_total / (V_total × ρcₚ_merged)

This is critical when merging cells with different compositions (e.g., solid vs gas),
as simple volume-weighted averaging can introduce significant temperature errors.

Concentrations are always computed conservatively (mass-preserving):
  ξ_merged = (V_remove * ξ_remove + V_keep * ξ_keep) / V_total
"""
function create_merged_state_vector(u::AbstractVector{Float64},
	sv::StateVectorALE{NC},
	cell_to_remove::Int,
	cell_to_keep::Int,
	V_remove::Float64,
	V_keep::Float64;
	material = nothing) where NC

	n_old = sv.n_cells
	n_new = n_old - 1
	V_total = V_remove + V_keep

	if n_new < 1
		error("Cannot merge: only 1 cell remaining")
	end

	# First compute merged state values
	T_remove = get_temperature_ale(u, sv, cell_to_remove)
	T_keep = get_temperature_ale(u, sv, cell_to_keep)

	# Get concentrations
	ξ_remove = ntuple(Val(NC)) do j
		get_concentration_ale(u, sv, cell_to_remove, j)
	end
	ξ_keep = ntuple(Val(NC)) do j
		get_concentration_ale(u, sv, cell_to_keep, j)
	end

	# Compute merged concentrations (always mass-conserving)
	ξ_merged = ntuple(Val(NC)) do j
		(V_remove * ξ_remove[j] + V_keep * ξ_keep[j]) / V_total
	end

	# Compute merged temperature
	if material !== nothing
		# Energy-conserving merge: compute ρcₚ-weighted energy
		ρcₚ_remove = 0.0
		ρcₚ_keep = 0.0
		@inbounds for j in 1:NC
			c_j = material.components[j].heat_capacity
			ρcₚ_remove += ξ_remove[j] * c_j(T_remove)
			ρcₚ_keep += ξ_keep[j] * c_j(T_keep)
		end

		# Total internal energy
		E_remove = V_remove * ρcₚ_remove * T_remove
		E_keep = V_keep * ρcₚ_keep * T_keep
		E_total = E_remove + E_keep

		# Merged volumetric heat capacity at estimated temperature
		T_estimate = (E_remove + E_keep) / (V_remove * ρcₚ_remove + V_keep * ρcₚ_keep + 1e-20)
		ρcₚ_merged = 0.0
		@inbounds for j in 1:NC
			c_j = material.components[j].heat_capacity
			ρcₚ_merged += ξ_merged[j] * c_j(T_estimate)
		end

		# Energy-conserving temperature
		if ρcₚ_merged > 1e-20
			T_merged = E_total / (V_total * ρcₚ_merged)
		else
			T_merged = T_estimate
		end
	else
		# Fallback: simple volume-weighted average
		T_merged = (V_remove * T_remove + V_keep * T_keep) / V_total
	end

	# Update the kept cell in the old state vector
	set_temperature_ale!(u, sv, cell_to_keep, T_merged)
	for j in 1:NC
		set_concentration_ale!(u, sv, cell_to_keep, j, ξ_merged[j])
	end

	# Determine which node to remove
	# After merging cells (cell_to_remove, cell_to_keep), the kept cell absorbs
	# both cells' volumes. The node that becomes orphaned is the one BETWEEN them.
	# 
	# In 1D: cell i uses nodes (i, i+1)
	# For surface merge where cell_to_remove > cell_to_keep:
	#   - cell_to_remove (n) uses nodes (n, n+1) where n+1 is surface
	#   - cell_to_keep (n-1) uses nodes (n-1, n)
	#   - After merge: kept cell uses nodes (n-1, n+1), spanning both cells
	#   - Node n (the shared interior node) is orphaned
	#   - Node n+1 (surface) is kept!
	if cell_to_remove > cell_to_keep
		node_to_remove = cell_to_remove  # Interior node between cells (not surface!)
	else
		node_to_remove = cell_to_keep  # Interior node between cells
	end

	# Now remove the cell from the state vector
	return remove_cell_from_state_vector(u, sv, cell_to_remove, node_to_remove)
end

"""
	map_cell_indices_after_removal(n_cells_old::Int, removed_cell::Int) -> Dict{Int,Int}

Create a mapping from old cell indices to new cell indices after a cell removal.

# Arguments
- `n_cells_old`: Number of cells before removal
- `removed_cell`: Index of the cell that was removed

# Returns
Dictionary mapping old_index => new_index

# Example
```julia
# If cell 5 was removed from 10 cells:
mapping = map_cell_indices_after_removal(10, 5)
# mapping[1] = 1, mapping[2] = 2, ..., mapping[4] = 4
# mapping[5] doesn't exist (removed)
# mapping[6] = 5, mapping[7] = 6, ..., mapping[10] = 9
```
"""
function map_cell_indices_after_removal(n_cells_old::Int, removed_cell::Int)
	mapping = Dict{Int, Int}()

	new_idx = 0
	for old_idx in 1:n_cells_old
		if old_idx != removed_cell
			new_idx += 1
			mapping[old_idx] = new_idx
		end
	end

	return mapping
end

# ============================================================================
# P1 Radiation State Vector Support
# ============================================================================

"""
	StateVectorALERadiation{NC}

Extended ALE state vector descriptor that includes radiation intensity G.

For P1 radiation, the state vector layout is:
```
u = [T₁, ..., Tₙ, G₁, ..., Gₙ, ξ₁₁, ..., ξ_NC,n, z₁, ..., z_{n+1}]
```

Where:
- T: Cell temperatures (n values)
- G: Incident radiation intensity (n values) - NEW for P1
- ξ: Mass concentrations (NC × n values, column-major by component)
- z: Node positions (n+1 values for 1D mesh with n cells)

# Type Parameters
- `NC`: Number of components

# Fields
- `n_cells::Int`: Number of cells
- `n_components::Int`: Number of material components
- `n_nodes::Int`: Number of nodes (n_cells + 1 for 1D)
- `T_offset::Int`: Offset to temperature data (always 0)
- `G_offset::Int`: Offset to radiation intensity data (= n_cells)
- `ξ_offset::Int`: Offset to concentration data (= 2*n_cells)
- `z_offset::Int`: Offset to node position data (= 2*n_cells + NC*n_cells)

# Example
For n=4 cells, NC=2 components, the state vector has 21 elements:
- u[1:4]: Temperatures T₁...T₄
- u[5:8]: Radiation intensities G₁...G₄
- u[9:12]: Component 1 concentrations ξ₁₁...ξ₁₄
- u[13:16]: Component 2 concentrations ξ₂₁...ξ₂₄
- u[17:21]: Node positions z₁...z₅
"""
struct StateVectorALERadiation{NC}
	n_cells::Int
	n_components::Int
	n_nodes::Int
	T_offset::Int      # Always 0
	G_offset::Int      # = n_cells
	ξ_offset::Int      # = 2*n_cells
	z_offset::Int      # = 2*n_cells + NC*n_cells
end

"""
	StateVectorALERadiation(n_cells::Int, ::Val{NC}) where NC

Create a StateVectorALERadiation descriptor.
"""
function StateVectorALERadiation(n_cells::Int, ::Val{NC}) where NC
	n_nodes = n_cells + 1
	T_offset = 0
	G_offset = n_cells
	ξ_offset = 2 * n_cells
	z_offset = 2 * n_cells + NC * n_cells
	StateVectorALERadiation{NC}(n_cells, NC, n_nodes, T_offset, G_offset, ξ_offset, z_offset)
end

"""
	StateVectorALERadiation(mesh::Mesh1D{NC}) where NC

Create a StateVectorALERadiation descriptor from a mesh.
"""
StateVectorALERadiation(mesh::Mesh1D{NC}) where NC = StateVectorALERadiation(length(mesh.cells), Val(NC))

"""
	state_vector_length_ale_radiation(n_cells::Int, NC::Int) -> Int

Compute the total length of the ALE+radiation state vector.

Returns `2*n_cells + NC*n_cells + (n_cells+1)` for T + G + ξ + z.
"""
@inline function state_vector_length_ale_radiation(n_cells::Int, NC::Int)
	2 * n_cells + NC * n_cells + (n_cells + 1)
end

"""
	state_vector_length_ale_radiation(sv::StateVectorALERadiation{NC}) where NC -> Int

Compute the total length from a StateVectorALERadiation descriptor.
"""
@inline function state_vector_length_ale_radiation(sv::StateVectorALERadiation{NC}) where NC
	2 * sv.n_cells + NC * sv.n_cells + sv.n_nodes
end

# ============================================================================
# Create P1 Radiation State Vector
# ============================================================================

"""
	create_state_vector_ale_radiation(mesh::Mesh1D{NC}, rad_cache) where NC -> Vector{Float64}

Create and initialize an ALE state vector with P1 radiation.

# Arguments
- `mesh`: Mesh with populated cell_states and nodes
- `rad_cache`: RadiationCache with G field

# Returns
A flat vector containing packed state including G and node positions.
"""
function create_state_vector_ale_radiation(mesh::Mesh1D{NC}, rad_cache) where NC
	sv = StateVectorALERadiation(mesh)
	u = Vector{Float64}(undef, state_vector_length_ale_radiation(sv))
	pack_state_ale_radiation!(u, mesh, rad_cache, sv)
	return u
end

# ============================================================================
# Pack P1 Radiation State (Mesh → Vector)
# ============================================================================

"""
	pack_state_ale_radiation!(u::AbstractVector{Float64}, mesh::Mesh1D{NC}, 
							  rad_cache, sv::StateVectorALERadiation{NC}) where NC

Pack mesh cell states, radiation intensity, and node positions into ODE vector.

Layout:
- `u[1:n]` = temperatures T₁, ..., Tₙ
- `u[n+1:2n]` = radiation intensities G₁, ..., Gₙ
- `u[2n+1:2n+NC*n]` = concentrations (column-major by component)
- `u[2n+NC*n+1:end]` = node positions z₁, ..., z_{n+1}
"""
function pack_state_ale_radiation!(u::AbstractVector{Float64}, mesh::Mesh1D{NC},
	rad_cache, sv::StateVectorALERadiation{NC}) where NC
	n = sv.n_cells

	# Pack temperatures
	@inbounds for i in 1:n
		u[sv.T_offset+i] = mesh.cell_states[i].T
	end

	# Pack radiation intensities
	@inbounds for i in 1:n
		u[sv.G_offset+i] = rad_cache.G[i]
	end

	# Pack concentrations (column-major by component)
	@inbounds for j in 1:NC
		for i in 1:n
			u[sv.ξ_offset+(j-1)*n+i] = mesh.cell_states[i].ξ[j]
		end
	end

	# Pack node positions
	@inbounds for i in 1:sv.n_nodes
		u[sv.z_offset+i] = mesh.nodes[i].z
	end

	return u
end

# ============================================================================
# Unpack P1 Radiation State (Vector → Mesh)
# ============================================================================

"""
	_unpack_concentrations_ale_radiation(u, ξ_offset, n, i, ::Val{NC}) where NC

Generated function to unpack concentrations from radiation state vector.
"""
@generated function _unpack_concentrations_ale_radiation(u::AbstractVector, ξ_offset::Int,
	n::Int, i::Int, ::Val{NC}) where NC
	exprs = [:(Float64(real(u[ξ_offset+$(j-1)*n+i]))) for j in 1:NC]
	return :(($(exprs...),))
end

"""
	unpack_state_ale_radiation!(mesh::Mesh1D{NC}, rad_cache, u::AbstractVector,
								sv::StateVectorALERadiation{NC}) where NC

Unpack flat ODE vector into mesh cell states, radiation cache, and node positions.
"""
function unpack_state_ale_radiation!(mesh::Mesh1D{NC}, rad_cache, u::AbstractVector{T},
	sv::StateVectorALERadiation{NC}) where {NC, T}

	# Skip mesh updates during AD (when T is a Dual type)
	T <: ForwardDiff.Dual && return mesh

	n = sv.n_cells
	ξ_offset = sv.ξ_offset

	# Unpack temperatures and concentrations to cell_states
	@inbounds for i in 1:n
		state = mesh.cell_states[i]

		# Unpack temperature
		state.T = Float64(real(u[sv.T_offset+i]))

		# Unpack concentrations using generated function
		state.ξ = _unpack_concentrations_ale_radiation(u, ξ_offset, n, i, Val(NC))
	end

	# Unpack radiation intensities to cache
	@inbounds for i in 1:n
		rad_cache.G[i] = Float64(real(u[sv.G_offset+i]))
	end

	# Unpack node positions and update geometry
	@inbounds for i in 1:sv.n_nodes
		z_new = Float64(real(u[sv.z_offset+i]))
		node = mesh.nodes[i]
		mesh.nodes[i] = Node1D(node.id, z_new, node.z_ref, node.is_boundary, node.boundary_id)
	end

	# Recompute geometry from new node positions
	update_geometry_from_nodes!(mesh)

	return mesh
end

# ============================================================================
# P1 Radiation State Vector Accessors
# ============================================================================

"""
	get_temperature_ale_rad(u::AbstractVector, sv::StateVectorALERadiation, i::Int)

Get temperature of cell i from radiation ALE state vector.
"""
@inline function get_temperature_ale_rad(u::AbstractVector, sv::StateVectorALERadiation, i::Int)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError(u, i))
	@inbounds u[sv.T_offset+i]
end

"""
	get_radiation_intensity(u::AbstractVector, sv::StateVectorALERadiation, i::Int)

Get radiation intensity G of cell i from state vector.
"""
@inline function get_radiation_intensity(u::AbstractVector, sv::StateVectorALERadiation, i::Int)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@inbounds u[sv.G_offset+i]
end

"""
	get_concentration_ale_rad(u::AbstractVector, sv::StateVectorALERadiation{NC}, i::Int, j::Int) where NC

Get concentration of component j in cell i from radiation ALE state vector.
"""
@inline function get_concentration_ale_rad(u::AbstractVector, sv::StateVectorALERadiation{NC},
	i::Int, j::Int) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	@inbounds u[sv.ξ_offset+(j-1)*n+i]
end

"""
	get_node_position_rad(u::AbstractVector, sv::StateVectorALERadiation, i::Int)

Get position of node i from radiation ALE state vector.
"""
@inline function get_node_position_rad(u::AbstractVector, sv::StateVectorALERadiation, i::Int)
	@boundscheck 1 <= i <= sv.n_nodes || throw(BoundsError())
	@inbounds u[sv.z_offset+i]
end

"""
	set_temperature_ale_rad!(u::AbstractVector, sv::StateVectorALERadiation, i::Int, T)

Set temperature of cell i in radiation ALE state vector.
"""
@inline function set_temperature_ale_rad!(u::AbstractVector, sv::StateVectorALERadiation, i::Int, T)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@inbounds u[sv.T_offset+i] = T
end

"""
	set_radiation_intensity!(u::AbstractVector, sv::StateVectorALERadiation, i::Int, G)

Set radiation intensity G of cell i in state vector.
"""
@inline function set_radiation_intensity!(u::AbstractVector, sv::StateVectorALERadiation, i::Int, G)
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@inbounds u[sv.G_offset+i] = G
end

"""
	set_concentration_ale_rad!(u::AbstractVector, sv::StateVectorALERadiation{NC}, i::Int, j::Int, ξ) where NC

Set concentration of component j in cell i in radiation ALE state vector.
"""
@inline function set_concentration_ale_rad!(u::AbstractVector, sv::StateVectorALERadiation{NC},
	i::Int, j::Int, ξ) where NC
	@boundscheck 1 <= i <= sv.n_cells || throw(BoundsError())
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	@inbounds u[sv.ξ_offset+(j-1)*n+i] = ξ
end

"""
	set_node_position_rad!(u::AbstractVector, sv::StateVectorALERadiation, i::Int, z)

Set position of node i in radiation ALE state vector.
"""
@inline function set_node_position_rad!(u::AbstractVector, sv::StateVectorALERadiation, i::Int, z)
	@boundscheck 1 <= i <= sv.n_nodes || throw(BoundsError())
	@inbounds u[sv.z_offset+i] = z
end

# ============================================================================
# P1 Radiation Views
# ============================================================================

"""
	temperature_view_ale_rad(u::AbstractVector, sv::StateVectorALERadiation)

Get a view of the temperature portion of the radiation ALE state vector.
"""
@inline function temperature_view_ale_rad(u::AbstractVector, sv::StateVectorALERadiation)
	@view u[(sv.T_offset+1):(sv.T_offset+sv.n_cells)]
end

"""
	radiation_view(u::AbstractVector, sv::StateVectorALERadiation)

Get a view of the radiation intensity (G) portion of the state vector.
"""
@inline function radiation_view(u::AbstractVector, sv::StateVectorALERadiation)
	@view u[(sv.G_offset+1):(sv.G_offset+sv.n_cells)]
end

"""
	concentration_view_ale_rad(u::AbstractVector, sv::StateVectorALERadiation{NC}, j::Int) where NC

Get a view of the concentration of component j from radiation ALE state vector.
"""
@inline function concentration_view_ale_rad(u::AbstractVector, sv::StateVectorALERadiation{NC},
	j::Int) where NC
	@boundscheck 1 <= j <= NC || throw(BoundsError())
	n = sv.n_cells
	start_idx = sv.ξ_offset + (j - 1) * n + 1
	end_idx = sv.ξ_offset + j * n
	@view u[start_idx:end_idx]
end

"""
	node_positions_view_rad(u::AbstractVector, sv::StateVectorALERadiation)

Get a view of the node positions portion of the radiation ALE state vector.
"""
@inline function node_positions_view_rad(u::AbstractVector, sv::StateVectorALERadiation)
	@view u[(sv.z_offset+1):(sv.z_offset+sv.n_nodes)]
end

"""
	validate_state_vector_ale_radiation(u::AbstractVector, sv::StateVectorALERadiation{NC}) where NC

Check that radiation ALE state vector has correct length.
"""
function validate_state_vector_ale_radiation(u::AbstractVector, sv::StateVectorALERadiation{NC}) where NC
	expected = state_vector_length_ale_radiation(sv)
	actual = length(u)
	if actual != expected
		throw(DimensionMismatch("Radiation ALE state vector length $actual does not match expected $expected"))
	end
	return true
end
