"""
Mesh container types.

Reference: Design Document Section 6.3
"""

# ============================================================================# Boundary State (for surface temperature solve)
# ============================================================================

"""
	BoundaryState

Stores solved boundary temperature and flux balance information.
Used to properly resolve the surface temperature that satisfies energy balance:
	q_applied(T_s) = q_conducted(T_s, T_cell)

# Fields
- `T_surface::Float64`: Solved surface temperature [K]
- `q_applied::Float64`: Heat flux from boundary condition [W/m²]
- `q_conducted::Float64`: Heat flux conducted into material [W/m²]
- `T_cell::Float64`: Adjacent cell temperature [K]
- `k_eff::Float64`: Effective conductivity at boundary [W/(m·K)]
- `Δz_half::Float64`: Distance from boundary to cell center [m]
- `converged::Bool`: Did Newton iteration converge?
- `iterations::Int`: Number of iterations used
"""
mutable struct BoundaryState
	T_surface::Float64
	q_applied::Float64
	q_conducted::Float64
	T_cell::Float64
	k_eff::Float64
	Δz_half::Float64
	converged::Bool
	iterations::Int
end

# Default constructor
BoundaryState() = BoundaryState(300.0, 0.0, 0.0, 300.0, 0.2, 1e-5, false, 0)

# ============================================================================# 1D Mesh
# ============================================================================

"""
	Mesh1D <: AbstractMesh{1}

A 1D finite volume mesh.

# Fields
- `nodes::Vector{Node1D}`: Mesh nodes
- `cells::Vector{Cell1D}`: Mesh cells
- `faces::Vector{Face1D}`: Mesh faces
- `cell_states::Vector{CellState}`: Cell state variables
- `cell_centers::Vector{Float64}`: Cell center positions
- `cell_volumes::Vector{Float64}`: Cell volumes
- `face_positions::Vector{Float64}`: Face positions
- `face_areas::Vector{Float64}`: Face areas (cross-sectional)
- `boundary_faces::Dict{Int,Vector{Int}}`: Boundary ID to face indices
- `max_level::Int`: Maximum refinement level
- `n_active_cells::Int`: Number of active (leaf) cells
- `active_cell_mask::BitVector`: Mask for active cells
- `node_velocities::Vector{Float64}`: ALE mesh velocities
- `mesh_time::Float64`: Time of current mesh configuration
- `cross_section_area::Float64`: Cross-sectional area (for 1D slab)
- `boundary_states::Dict{Int,BoundaryState}`: Boundary ID to solved boundary state

# Type Parameter
- `NC`: Number of components (for type-stable cell state access)
"""
mutable struct Mesh1D{NC} <: AbstractMesh{1}
	# Topology
	nodes::Vector{Node1D}
	cells::Vector{Cell1D}
	faces::Vector{Face1D}

	# State
	cell_states::Vector{CellState{NC}}

	# Geometry cache (recomputed after mesh motion)
	cell_centers::Vector{Float64}
	cell_volumes::Vector{Float64}
	face_positions::Vector{Float64}
	face_areas::Vector{Float64}      # = cross-sectional area

	# Connectivity helpers
	boundary_faces::Dict{Int, Vector{Int}}  # boundary_id => face indices

	# Refinement tracking
	max_level::Int
	n_active_cells::Int              # excludes refined-away cells
	active_cell_mask::BitVector

	# ALE state
	node_velocities::Vector{Float64}
	mesh_time::Float64               # time of current mesh configuration

	# Domain properties
	cross_section_area::Float64      # for 1D: unit area typically 1.0

	# Boundary state (for surface temperature solve)
	boundary_states::Dict{Int, BoundaryState}  # boundary_id => state
end

# Default constructor (parameterized by NC)
function Mesh1D{NC}() where NC
	Mesh1D{NC}(
		Node1D[],
		Cell1D[],
		Face1D[],
		CellState{NC}[],
		Float64[],
		Float64[],
		Float64[],
		Float64[],
		Dict{Int, Vector{Int}}(),
		0,
		0,
		BitVector(),
		Float64[],
		0.0,
		1.0,
		Dict{Int, BoundaryState}(1 => BoundaryState(), 2 => BoundaryState()),
	)
end

# Convenience constructor from Val{NC}
Mesh1D(::Val{NC}) where NC = Mesh1D{NC}()

# ============================================================================
# 2D Mesh (Axisymmetric)
# ============================================================================

"""
	Mesh2D{NC} <: AbstractMesh{2}

A 2D axisymmetric finite volume mesh.

# Fields
Similar to Mesh1D but with 2D types and additional axisymmetric factors.

# Type Parameter
- `NC`: Number of components (for type-stable cell state access)
"""
mutable struct Mesh2D{NC} <: AbstractMesh{2}
	# Topology
	nodes::Vector{Node2D}
	cells::Vector{Cell2D}
	faces::Vector{Face2D}
	cell_states::Vector{CellState{NC}}

	# Geometry cache
	cell_centers::Vector{Tuple{Float64, Float64}}
	cell_volumes::Vector{Float64}
	face_centers::Vector{Tuple{Float64, Float64}}
	face_areas::Vector{Float64}
	face_normals::Vector{Tuple{Float64, Float64}}

	# Axisymmetric factors
	r_factors::Vector{Float64}       # 2πr for volume integration

	# Connectivity
	boundary_faces::Dict{Int, Vector{Int}}

	# Refinement
	max_level::Int
	n_active_cells::Int
	active_cell_mask::BitVector

	# ALE
	node_velocities::Vector{Tuple{Float64, Float64}}
	mesh_time::Float64

	# Boundary state (for surface temperature solve)
	boundary_states::Dict{Int, BoundaryState}  # boundary_id => state
end

# Default constructor (parameterized by NC)
function Mesh2D{NC}() where NC
	Mesh2D{NC}(
		Node2D[],
		Cell2D[],
		Face2D[],
		CellState{NC}[],
		Tuple{Float64, Float64}[],
		Float64[],
		Tuple{Float64, Float64}[],
		Float64[],
		Tuple{Float64, Float64}[],
		Float64[],
		Dict{Int, Vector{Int}}(),
		0,
		0,
		BitVector(),
		Tuple{Float64, Float64}[],
		0.0,
		Dict{Int, BoundaryState}(1 => BoundaryState(), 2 => BoundaryState()),
	)
end

# Convenience constructor from Val{NC}
Mesh2D(::Val{NC}) where NC = Mesh2D{NC}()

# ============================================================================
# Mesh Accessors
# ============================================================================

"""
	n_nodes(mesh::AbstractMesh)

Get the number of nodes in the mesh.
"""
@inline n_nodes(mesh::Mesh1D) = length(mesh.nodes)
@inline n_nodes(mesh::Mesh2D) = length(mesh.nodes)

"""
	n_cells(mesh::AbstractMesh)

Get the total number of cells in the mesh (including inactive).
"""
@inline n_cells(mesh::Mesh1D) = length(mesh.cells)
@inline n_cells(mesh::Mesh2D) = length(mesh.cells)

"""
	n_active_cells(mesh::AbstractMesh)

Get the number of active (leaf) cells in the mesh.
"""
@inline n_active_cells(mesh::Mesh1D) = mesh.n_active_cells
@inline n_active_cells(mesh::Mesh2D) = mesh.n_active_cells

"""
	n_faces(mesh::AbstractMesh)

Get the number of faces in the mesh.
"""
@inline n_faces(mesh::Mesh1D) = length(mesh.faces)
@inline n_faces(mesh::Mesh2D) = length(mesh.faces)

"""
	get_node(mesh::AbstractMesh, i::Int)

Get the i-th node.
"""
@inline get_node(mesh::Mesh1D, i::Int) = mesh.nodes[i]
@inline get_node(mesh::Mesh2D, i::Int) = mesh.nodes[i]

"""
	get_cell(mesh::AbstractMesh, i::Int)

Get the i-th cell.
"""
@inline get_cell(mesh::Mesh1D, i::Int) = mesh.cells[i]
@inline get_cell(mesh::Mesh2D, i::Int) = mesh.cells[i]

"""
	get_face(mesh::AbstractMesh, i::Int)

Get the i-th face.
"""
@inline get_face(mesh::Mesh1D, i::Int) = mesh.faces[i]
@inline get_face(mesh::Mesh2D, i::Int) = mesh.faces[i]

"""
	get_state(mesh::AbstractMesh, i::Int)

Get the state of the i-th cell.
"""
@inline get_state(mesh::Mesh1D, i::Int) = mesh.cell_states[i]
@inline get_state(mesh::Mesh2D, i::Int) = mesh.cell_states[i]

"""
	is_active(mesh::AbstractMesh, i::Int)

Check if the i-th cell is active.
"""
@inline is_active(mesh::Mesh1D, i::Int) = mesh.active_cell_mask[i]
@inline is_active(mesh::Mesh2D, i::Int) = mesh.active_cell_mask[i]

"""
	boundary_face_ids(mesh::AbstractMesh, boundary_id::Int)

Get face indices for a given boundary.
"""
function boundary_face_ids(mesh::Union{Mesh1D, Mesh2D}, boundary_id::Int)
	haskey(mesh.boundary_faces, boundary_id) ? mesh.boundary_faces[boundary_id] : Int[]
end

"""
	domain_length(mesh::Mesh1D)

Get the domain length (thickness) for a 1D mesh.
"""
function domain_length(mesh::Mesh1D)
	if isempty(mesh.nodes)
		return 0.0
	end
	z_min = minimum(n.z for n in mesh.nodes)
	z_max = maximum(n.z for n in mesh.nodes)
	return z_max - z_min
end
