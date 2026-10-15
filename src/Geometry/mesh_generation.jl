"""
Mesh generation utilities.
"""

# ============================================================================
# 1D Mesh Generation
# ============================================================================

"""
	create_uniform_mesh(L::Real, n_cells::Int, ::Val{NC}; 
						T_initial::Real=300.0, 
						ξ_initial=nothing,
						cross_section_area::Real=1.0) -> Mesh1D{NC}

Create a uniform 1D mesh.

# Arguments
- `L`: Domain length (thickness) [m]
- `n_cells`: Number of cells
- `Val{NC}`: Number of components as a type parameter (use `Val(NC)`)
- `T_initial=300.0`: Initial temperature [K]
- `ξ_initial=nothing`: Initial concentrations (NTuple or scalar)
- `cross_section_area=1.0`: Cross-sectional area [m²]

# Returns
A `Mesh1D{NC}` with uniform cell spacing.

# Convention
- z = 0 is the BOTTOM surface (substrate, boundary_id = 2)
- z = L is the TOP surface (exposed surface, boundary_id = 1)
- Cell 1 is at z ≈ 0 (bottom), Cell N is at z ≈ L (top)
- Heat enters from the top (z = L), material shrinks toward bottom

# Example
```julia
mesh = create_uniform_mesh(0.01, 100, Val(2))  # 1 cm, 100 cells, 2 components
```
"""
function create_uniform_mesh(
    L::Real,
    n_cells::Int,
    ::Val{NC};
    T_initial::Real = 300.0,
    ξ_initial = nothing,
    cross_section_area::Real = 1.0,
) where {NC}
    L = Float64(L)
    n_nodes = n_cells + 1
    Δz = L / n_cells
    A = Float64(cross_section_area)

    # Initialize mesh
    mesh = Mesh1D{NC}()
    mesh.cross_section_area = A
    mesh.cross_section_area_initial = A

    # Create nodes (z=0 at bottom, z=L at top)
    for i = 1:n_nodes
        z = (i - 1) * Δz
        if i == 1
            push!(mesh.nodes, Node1D(i, z, 2))  # bottom boundary (substrate)
        elseif i == n_nodes
            push!(mesh.nodes, Node1D(i, z, 1))  # top boundary (exposed surface)
        else
            push!(mesh.nodes, Node1D(i, z, 0))  # interior
        end
    end

    # Create faces (n_cells + 1 faces, including boundaries)
    for i = 1:n_nodes
        if i == 1
            # Bottom boundary face (only right cell exists, left_id=0)
            push!(mesh.faces, Face1D(i, i, (0, 1), 2))
        elseif i == n_nodes
            # Top boundary face (only left cell exists, right_id=0)
            push!(mesh.faces, Face1D(i, i, (n_cells, 0), 1))
        else
            # Interior face
            push!(mesh.faces, Face1D(i, i, (i-1, i)))
        end
    end

    # Create cells
    for i = 1:n_cells
        cell = Cell1D(i, (i, i+1), (i, i+1))
        push!(mesh.cells, cell)
    end

    # Initialize cell states
    if ξ_initial === nothing
        ξ_init = ntuple(_ -> 0.0, NC)
    elseif ξ_initial isa Real
        ξ_init = ntuple(_ -> Float64(ξ_initial), NC)
    else
        ξ_init = ntuple(j -> Float64(ξ_initial[j]), NC)
    end

    V_cell = Δz * A
    for i = 1:n_cells
        state = CellState{NC}(Float64(T_initial), ξ_init, V_cell)
        push!(mesh.cell_states, state)
    end

    # Initialize geometry cache
    resize!(mesh.cell_centers, n_cells)
    resize!(mesh.cell_volumes, n_cells)
    resize!(mesh.face_positions, n_nodes)
    resize!(mesh.face_areas, n_nodes)

    # Initialize ALE state
    resize!(mesh.node_velocities, n_nodes)
    fill!(mesh.node_velocities, 0.0)
    mesh.mesh_time = 0.0

    # Initialize refinement tracking
    mesh.max_level = 0
    mesh.n_active_cells = n_cells
    mesh.active_cell_mask = trues(n_cells)

    # Pyrolysis-progress denominator (populated by problem setup)
    resize!(mesh.m_dry_solid_initial, n_cells)
    fill!(mesh.m_dry_solid_initial, 0.0)

    # Set up boundary face dictionary
    mesh.boundary_faces[1] = [n_nodes]  # top (exposed surface at z=L)
    mesh.boundary_faces[2] = [1]        # bottom (substrate at z=0)

    # Compute geometry
    update_geometry!(mesh)

    return mesh
end

"""
	create_stretched_mesh(L::Real, n_cells::Int, ::Val{NC}, stretch::Real;
						  T_initial::Real=300.0,
						  ξ_initial=nothing,
						  cross_section_area::Real=1.0,
						  stretch_from::Symbol=:top) -> Mesh1D{NC}

Create a geometrically stretched 1D mesh with the **finest cells at the
boundary named by `stretch_from`** and cell sizes growing geometrically
(ratio `stretch`) away from it.

# Arguments
- `L`: Domain length [m]
- `n_cells`: Number of cells
- `Val{NC}`: Number of components as a type parameter (use `Val(NC)`)
- `stretch`: Geometric size ratio between adjacent cells (must be ≥ 1;
  1 = uniform). Each cell is `stretch` times thicker than its neighbor on
  the refined side.
- `T_initial=300.0`: Initial temperature [K]
- `ξ_initial=nothing`: Initial concentrations
- `cross_section_area=1.0`: Cross-sectional area [m²]
- `stretch_from=:top`: Which boundary gets the finest cells — `:top`
  (exposed surface at z = L, the default) or `:bottom` (substrate at z = 0)

# Returns
A `Mesh1D{NC}` with geometric stretching.

# Example
```julia
# Fine mesh near top surface (z=L), coarser toward bottom (z=0)
mesh = create_stretched_mesh(0.01, 50, Val(2), 1.1; stretch_from=:top)
```
"""
function create_stretched_mesh(
    L::Real,
    n_cells::Int,
    ::Val{NC},
    stretch::Real;
    T_initial::Real = 300.0,
    ξ_initial = nothing,
    cross_section_area::Real = 1.0,
    stretch_from::Symbol = :top,
) where {NC}
    L = Float64(L)
    stretch = Float64(stretch)
    n_nodes = n_cells + 1
    A = Float64(cross_section_area)

    stretch_from in (:top, :bottom) || throw(
        ArgumentError("stretch_from must be :top or :bottom, got $(repr(stretch_from))"),
    )
    stretch >= 1.0 - 1e-10 || throw(
        ArgumentError(
            "stretch must be ≥ 1 (got $stretch); use stretch_from to choose " *
            "which boundary gets the finest cells",
        ),
    )

    # Compute node positions with geometric stretching, built fine-end-first
    # from z = 0: the first cell is smallest and cells grow geometrically.
    if abs(stretch - 1.0) < 1e-10
        # Uniform case
        positions = collect(range(0.0, L, length = n_nodes))
    else
        # Geometric series: Δz_i = Δz_0 * stretch^i
        # Sum: L = Δz_0 * (stretch^n - 1) / (stretch - 1)
        Δz_0 = L * (stretch - 1.0) / (stretch^n_cells - 1.0)

        positions = zeros(n_nodes)
        positions[1] = 0.0
        for i = 2:n_nodes
            Δz_i = Δz_0 * stretch^(i-2)
            positions[i] = positions[i-1] + Δz_i
        end

        # Normalize to exactly L
        positions .*= L / positions[end]
    end

    # The series above puts the finest cells at z = 0 (the bottom). Reflect
    # for :top so the finest cells sit at the exposed surface z = L.
    if stretch_from == :top
        positions = L .- reverse(positions)
    end
    # Pin the endpoints exactly (renormalization/reflection can leave ≤1 ulp)
    positions[1] = 0.0
    positions[end] = L

    # Initialize mesh
    mesh = Mesh1D{NC}()
    mesh.cross_section_area = A
    mesh.cross_section_area_initial = A

    # Create nodes (z=0 at bottom, z=L at top)
    for i = 1:n_nodes
        z = positions[i]
        if i == 1
            push!(mesh.nodes, Node1D(i, z, 2))  # bottom boundary (substrate)
        elseif i == n_nodes
            push!(mesh.nodes, Node1D(i, z, 1))  # top boundary (exposed surface)
        else
            push!(mesh.nodes, Node1D(i, z, 0))  # interior
        end
    end

    # Create faces (z=0 at bottom, z=L at top)
    for i = 1:n_nodes
        if i == 1
            # Bottom boundary face (left_id=0)
            push!(mesh.faces, Face1D(i, i, (0, 1), 2))
        elseif i == n_nodes
            # Top boundary face (right_id=0)
            push!(mesh.faces, Face1D(i, i, (n_cells, 0), 1))
        else
            push!(mesh.faces, Face1D(i, i, (i-1, i)))
        end
    end

    # Create cells
    for i = 1:n_cells
        cell = Cell1D(i, (i, i+1), (i, i+1))
        push!(mesh.cells, cell)
    end

    # Initialize cell states
    if ξ_initial === nothing
        ξ_init = ntuple(_ -> 0.0, NC)
    elseif ξ_initial isa Real
        ξ_init = ntuple(_ -> Float64(ξ_initial), NC)
    else
        ξ_init = ntuple(j -> Float64(ξ_initial[j]), NC)
    end

    for i = 1:n_cells
        Δz = positions[i+1] - positions[i]
        V_cell = Δz * A
        state = CellState{NC}(Float64(T_initial), ξ_init, V_cell)
        push!(mesh.cell_states, state)
    end

    # Initialize geometry cache
    resize!(mesh.cell_centers, n_cells)
    resize!(mesh.cell_volumes, n_cells)
    resize!(mesh.face_positions, n_nodes)
    resize!(mesh.face_areas, n_nodes)

    # Initialize ALE state
    resize!(mesh.node_velocities, n_nodes)
    fill!(mesh.node_velocities, 0.0)
    mesh.mesh_time = 0.0

    # Initialize refinement tracking
    mesh.max_level = 0
    mesh.n_active_cells = n_cells
    mesh.active_cell_mask = trues(n_cells)

    # Pyrolysis-progress denominator (populated by problem setup)
    resize!(mesh.m_dry_solid_initial, n_cells)
    fill!(mesh.m_dry_solid_initial, 0.0)

    # Set up boundary faces
    mesh.boundary_faces[1] = [n_nodes]  # top (exposed surface at z=L)
    mesh.boundary_faces[2] = [1]        # bottom (substrate at z=0)

    # Compute geometry
    update_geometry!(mesh)

    return mesh
end

"""
    create_mesh_from_nodes(positions, Val(NC); T_initial=300., ξ_initial=nothing,
                           cross_section_area=1.)

Create a 1D mesh with prescribed, strictly increasing node positions from zero.
All current/reference cell volumes, geometry caches, and ALE bookkeeping are
initialized from these positions. Useful for preserving material interfaces
while refining each region independently.
"""
function create_mesh_from_nodes(
    node_positions::AbstractVector{<:Real},
    ::Val{NC};
    T_initial::Real = 300.0,
    ξ_initial = nothing,
    cross_section_area::Real = 1.0,
) where {NC}
    NC isa Integer && NC > 0 || throw(ArgumentError("NC must be a positive integer"))
    positions = Float64.(node_positions)
    length(positions) >= 2 ||
        throw(ArgumentError("at least two node positions are required"))
    all(isfinite, positions) || throw(ArgumentError("node positions must be finite"))
    first(positions) == 0.0 || throw(ArgumentError("node positions must start at zero"))
    all(>(0.0), diff(positions)) ||
        throw(ArgumentError("node positions must strictly increase"))
    isfinite(T_initial) && T_initial > 0 ||
        throw(ArgumentError("T_initial must be finite and positive"))
    A = Float64(cross_section_area)
    isfinite(A) && A > 0 ||
        throw(ArgumentError("cross_section_area must be finite and positive"))
    if ξ_initial !== nothing
        values = ξ_initial isa Real ? (ξ_initial,) : ξ_initial
        ξ_initial isa Real ||
            length(values) == NC ||
            throw(ArgumentError("ξ_initial must have NC entries"))
        all(isfinite, values) || throw(ArgumentError("ξ_initial must be finite"))
    end
    n_nodes = length(positions)
    n_cells = n_nodes - 1
    # Initialize mesh
    mesh = Mesh1D{NC}()
    mesh.cross_section_area = A
    mesh.cross_section_area_initial = A

    # Create nodes (z=0 at bottom, z=L at top)
    for i = 1:n_nodes
        z = positions[i]
        if i == 1
            push!(mesh.nodes, Node1D(i, z, 2))  # bottom boundary (substrate)
        elseif i == n_nodes
            push!(mesh.nodes, Node1D(i, z, 1))  # top boundary (exposed surface)
        else
            push!(mesh.nodes, Node1D(i, z, 0))  # interior
        end
    end

    # Create faces (z=0 at bottom, z=L at top)
    for i = 1:n_nodes
        if i == 1
            # Bottom boundary face (left_id=0)
            push!(mesh.faces, Face1D(i, i, (0, 1), 2))
        elseif i == n_nodes
            # Top boundary face (right_id=0)
            push!(mesh.faces, Face1D(i, i, (n_cells, 0), 1))
        else
            push!(mesh.faces, Face1D(i, i, (i-1, i)))
        end
    end

    # Create cells
    for i = 1:n_cells
        cell = Cell1D(i, (i, i+1), (i, i+1))
        push!(mesh.cells, cell)
    end

    # Initialize cell states
    if ξ_initial === nothing
        ξ_init = ntuple(_ -> 0.0, NC)
    elseif ξ_initial isa Real
        ξ_init = ntuple(_ -> Float64(ξ_initial), NC)
    else
        ξ_init = ntuple(j -> Float64(ξ_initial[j]), NC)
    end

    for i = 1:n_cells
        Δz = positions[i+1] - positions[i]
        V_cell = Δz * A
        state = CellState{NC}(Float64(T_initial), ξ_init, V_cell)
        push!(mesh.cell_states, state)
    end

    # Initialize geometry cache
    resize!(mesh.cell_centers, n_cells)
    resize!(mesh.cell_volumes, n_cells)
    resize!(mesh.face_positions, n_nodes)
    resize!(mesh.face_areas, n_nodes)

    # Initialize ALE state
    resize!(mesh.node_velocities, n_nodes)
    fill!(mesh.node_velocities, 0.0)
    mesh.mesh_time = 0.0

    # Initialize refinement tracking
    mesh.max_level = 0
    mesh.n_active_cells = n_cells
    mesh.active_cell_mask = trues(n_cells)

    # Pyrolysis-progress denominator (populated by problem setup)
    resize!(mesh.m_dry_solid_initial, n_cells)
    fill!(mesh.m_dry_solid_initial, 0.0)

    # Set up boundary faces
    mesh.boundary_faces[1] = [n_nodes]  # top (exposed surface at z=L)
    mesh.boundary_faces[2] = [1]        # bottom (substrate at z=0)

    # Compute geometry
    update_geometry!(mesh)

    return mesh
end
