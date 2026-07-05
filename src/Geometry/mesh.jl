"""
Mesh container types.
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
- `q_transpired::Float64`: Surface transpiration enthalpy flux `Σ_g ṁ_g·Δh_g(T_s,T_∞)`
  at the converged `T_s` [W/m²]. Zero under the default (transpiration-off) convention.
- `T_cell::Float64`: Adjacent cell temperature [K]
- `k_eff::Float64`: Effective conductivity at boundary [W/(m·K)]
- `Δz_half::Float64`: Distance from boundary to cell center [m]
- `converged::Bool`: Did Newton iteration converge?
- `iterations::Int`: Number of iterations used
- `blowing_mcp::Float64`: Net blowing conductance `ṁ″·c_p,g = Σ_g ṁ_g·cp_g(T_s)`
  [W/(m²·K)] at the converged `T_s` (see `Discretization.blowing_factor`).
  Zero when the blowing correction is off. Diagnostic consumers recover the
  per-BC correction factor as `blowing_factor(blowing_mcp / h)`.
"""
mutable struct BoundaryState
    T_surface::Float64
    q_applied::Float64
    q_conducted::Float64
    q_transpired::Float64
    T_cell::Float64
    k_eff::Float64
    Δz_half::Float64
    converged::Bool
    iterations::Int
    blowing_mcp::Float64
end

# Backward-compatible 9-argument constructor (every pre-blowing call site).
function BoundaryState(
    T_surface::Real,
    q_applied::Real,
    q_conducted::Real,
    q_transpired::Real,
    T_cell::Real,
    k_eff::Real,
    Δz_half::Real,
    converged::Bool,
    iterations::Integer,
)
    return BoundaryState(
        T_surface,
        q_applied,
        q_conducted,
        q_transpired,
        T_cell,
        k_eff,
        Δz_half,
        converged,
        iterations,
        0.0,
    )
end

# Default constructor
BoundaryState() = BoundaryState(300.0, 0.0, 0.0, 0.0, 300.0, 0.2, 1e-5, false, 0)

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

# Convenience constructor from Val{NC}
Mesh1D(::Val{NC}) where {NC} = Mesh1D{NC}()

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
