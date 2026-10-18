"""
Mesh container types.
"""

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
- `cross_section_area::Float64`: Current cross-sectional area (time-varying once
  variable cross-section is enabled; equal to `cross_section_area_initial` until then)
- `cross_section_area_initial::Float64`: Calibration-reference cross-sectional area,
  fixed at problem setup. Used by per-original-area diagnostics (e.g. MLR_gauge)
  and as the `A_0` in the lateral shrinkage law `A(t)/A_0 = law(χ̄)`.
- `m_dry_solid_initial::Vector{Float64}`: Initial dry-solid mass per cell [kg],
  fixed at problem setup. Denominator in the per-cell pyrolysis-progress variable
  χ_i = m_released_dry,i / m_dry_solid_initial,i (variable cross-section v2).
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
    boundary_faces::Dict{Int,Vector{Int}}  # boundary_id => face indices

    # Refinement tracking
    max_level::Int
    n_active_cells::Int              # excludes refined-away cells
    active_cell_mask::BitVector

    # ALE state
    node_velocities::Vector{Float64}
    mesh_time::Float64               # time of current mesh configuration

    # Domain properties
    cross_section_area::Float64          # current (time-varying) cross-sectional area
    cross_section_area_initial::Float64  # calibration reference A_0 (set once at setup)

    # Pyrolysis-progress denominator (per cell, set once at setup; updated by adaptivity)
    m_dry_solid_initial::Vector{Float64}

    # Boundary state (for surface temperature solve)
    boundary_states::Dict{Int,BoundaryState}  # boundary_id => state
end

# Default constructor (parameterized by NC)
function Mesh1D{NC}() where {NC}
    Mesh1D{NC}(
        Node1D[],
        Cell1D[],
        Face1D[],
        CellState{NC}[],
        Float64[],
        Float64[],
        Float64[],
        Float64[],
        Dict{Int,Vector{Int}}(),
        0,
        0,
        BitVector(),
        Float64[],
        0.0,
        1.0,
        1.0,
        Float64[],
        Dict{Int,BoundaryState}(1 => BoundaryState(), 2 => BoundaryState()),
    )
end

# ============================================================================
# Mesh Accessors
# ============================================================================

"""
	n_nodes(mesh::Mesh1D)

Get the number of nodes in the mesh.
"""
@inline n_nodes(mesh::Mesh1D) = length(mesh.nodes)

"""
	n_cells(mesh::Mesh1D)

Get the total number of cells in the mesh (including inactive).
"""
@inline n_cells(mesh::Mesh1D) = length(mesh.cells)

"""
	n_faces(mesh::Mesh1D)

Get the number of faces in the mesh.
"""
@inline n_faces(mesh::Mesh1D) = length(mesh.faces)
