"""
Mesh generation utilities.

Reference: Design Document Section 2.5, Stage 2
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
function create_uniform_mesh(L::Real, n_cells::Int, ::Val{NC};
	T_initial::Real = 300.0,
	ξ_initial = nothing,
	cross_section_area::Real = 1.0) where NC
	L = Float64(L)
	n_nodes = n_cells + 1
	Δz = L / n_cells
	A = Float64(cross_section_area)

	# Initialize mesh
	mesh = Mesh1D{NC}()
	mesh.cross_section_area = A

	# Create nodes (z=0 at bottom, z=L at top)
	for i in 1:n_nodes
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
	for i in 1:n_nodes
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
	for i in 1:n_cells
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
	for i in 1:n_cells
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

Create a geometrically stretched 1D mesh.

# Arguments
- `L`: Domain length [m]
- `n_cells`: Number of cells
- `Val{NC}`: Number of components as a type parameter (use `Val(NC)`)
- `stretch`: Geometric stretching factor (>1 = cells grow toward bottom)
- `T_initial=300.0`: Initial temperature [K]
- `ξ_initial=nothing`: Initial concentrations
- `cross_section_area=1.0`: Cross-sectional area [m²]
- `stretch_from=:top`: Direction of stretching (:top or :bottom)

# Returns
A `Mesh1D{NC}` with geometric stretching.

# Example
```julia
# Fine mesh near top surface (z=L), coarser toward bottom (z=0)
mesh = create_stretched_mesh(0.01, 50, Val(2), 1.1; stretch_from=:top)
```
"""
function create_stretched_mesh(L::Real, n_cells::Int, ::Val{NC}, stretch::Real;
	T_initial::Real = 300.0,
	ξ_initial = nothing,
	cross_section_area::Real = 1.0,
	stretch_from::Symbol = :top) where NC
	L = Float64(L)
	stretch = Float64(stretch)
	n_nodes = n_cells + 1
	A = Float64(cross_section_area)

	# Compute node positions with geometric stretching
	# For stretch > 1: first cell is smallest, cells grow geometrically
	if abs(stretch - 1.0) < 1e-10
		# Uniform case
		positions = collect(range(0.0, L, length = n_nodes))
	else
		# Geometric series: Δz_i = Δz_0 * stretch^i
		# Sum: L = Δz_0 * (stretch^n - 1) / (stretch - 1)
		Δz_0 = L * (stretch - 1.0) / (stretch^n_cells - 1.0)

		positions = zeros(n_nodes)
		positions[1] = 0.0
		for i in 2:n_nodes
			Δz_i = Δz_0 * stretch^(i-2)
			positions[i] = positions[i-1] + Δz_i
		end

		# Normalize to exactly L
		positions .*= L / positions[end]
	end

	# Reverse if stretching from bottom
	if stretch_from == :bottom
		positions = L .- reverse(positions)
	end

	# Initialize mesh
	mesh = Mesh1D{NC}()
	mesh.cross_section_area = A

	# Create nodes (z=0 at bottom, z=L at top)
	for i in 1:n_nodes
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
	for i in 1:n_nodes
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
	for i in 1:n_cells
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

	for i in 1:n_cells
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

	# Set up boundary faces
	mesh.boundary_faces[1] = [n_nodes]  # top (exposed surface at z=L)
	mesh.boundary_faces[2] = [1]        # bottom (substrate at z=0)

	# Compute geometry
	update_geometry!(mesh)

	return mesh
end
