"""
Tests for surface cell depletion detection and cell merging.

Tests cover:
1. Surface cell detection (find_surface_cell, find_bottom_cell)
2. Neighbor finding (find_interior_neighbor, find_exterior_neighbor)
3. Depletion detection (is_surface_depleted, is_cell_depleted)
4. Initial concentration computation
5. Adjacent cell merging (merge_adjacent_cells!)
6. State vector resizing for cell removal

References:
- docs/Surface_Cell_Depletion_Implementation.md
- src/Adaptivity/surface_depletion.jl
- src/Adaptivity/h_coarsening.jl
"""

using Test
using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Geometry
using Pyrolysis.Adaptivity
using Pyrolysis.TimeIntegration
using DiffEqCallbacks: DiscreteCallback

# ============================================================================
# Helper: Create test materials
# ============================================================================

"""Create a simple 2-component material (solid + gas) for testing."""
function create_test_material()
	# Solid component with constant properties
	solid = Component(
		:virgin,
		SOLID,
		0.0,  # molar mass (not used for solids)
		ConstantProperty(1200.0),  # density
		ConstantProperty(1500.0),  # heat capacity
		ConstantProperty(0.25),    # thermal conductivity
		ConstantProperty(0.0),     # gas transfer (not used)
		ConstantProperty(0.0),     # permeability (no Darcy)
		0.9,   # emissivity
		0.0,   # absorption coefficient
		0.0,    # swelling factor
	)

	# Gaseous component
	gas = Component(
		:pyrolyzate,
		GAS,
		0.044,  # molar mass
		ConstantProperty(100.0),   # gas density (will be overwritten)
		ConstantProperty(1000.0),  # heat capacity
		ConstantProperty(0.02),    # thermal conductivity
		ConstantProperty(1e-5),    # gas transfer coefficient
		nothing,                    # permeability (gases have no permeability)
		0.0,   # emissivity
		0.0,   # absorption coefficient
		1.0,    # swelling factor
	)

	# Simple decomposition reaction
	reaction = Reaction(
		:decomp,
		1, 0,      # reactant1, reactant2
		1.0, 0.0,  # order1, order2
		2, 0,      # product1, product2
		(1.0, 0.0, 1.0, 0.0),  # stoichiometry
		1e12,      # pre-exponential
		150e3,     # activation energy
		ConstantProperty(500e3),  # heat of reaction
		300.0,     # T_min
		2000.0,     # T_max
	)

	return Material(
		name = :test_material,
		components = (solid, gas),
		reactions = (reaction,),
		conductivity_mixing = WEIGHTED,
		β = 0.5,
		gas_transfer_mixing = WEIGHTED,
		β_λ = 0.5,
	)
end

# ============================================================================
# Helper: Create test mesh with varied composition
# ============================================================================

"""Create a mesh where surface cells are depleted of solid."""
function create_depleted_surface_mesh(; n_cells::Int = 10, L::Float64 = 0.01, NC::Int = 2)
	mesh = create_uniform_mesh(L, n_cells, Val(NC))

	# Set up composition: 
	# - Interior cells: mostly solid (virgin)
	# - Surface cell: mostly gas (depleted)
	ξ_solid_initial = 1200.0  # Full solid concentration

	for i in 1:n_cells
		z = mesh.cell_centers[i]
		T = 300.0 + 500.0 * z / L  # Temperature gradient

		if i == n_cells  # Surface cell (highest z)
			# Depleted: only 0.5% solid remaining
			ξ_solid = ξ_solid_initial * 0.005
			ξ_gas = ξ_solid_initial - ξ_solid
		elseif i == n_cells - 1  # Near-surface cell
			# Partially depleted: 50% solid remaining
			ξ_solid = ξ_solid_initial * 0.5
			ξ_gas = ξ_solid_initial - ξ_solid
		else  # Interior cells
			# Full solid
			ξ_solid = ξ_solid_initial
			ξ_gas = 0.0
		end

		ξ = (ξ_solid, ξ_gas)
		mesh.cell_states[i] = CellState{NC}(T, ξ, mesh.cell_volumes[i], mesh.cell_volumes[i])
	end

	return mesh
end

"""Create a simple uniform mesh for basic tests."""
function create_test_mesh(; n_cells::Int = 10, L::Float64 = 0.01, NC::Int = 2)
	mesh = create_uniform_mesh(L, n_cells, Val(NC))

	for i in 1:n_cells
		z = mesh.cell_centers[i]
		T = 300.0 + 500.0 * z / L
		ξ = (1200.0, 0.0)  # All solid, no gas
		mesh.cell_states[i] = CellState{NC}(T, ξ, mesh.cell_volumes[i], mesh.cell_volumes[i])
	end

	return mesh
end

# ============================================================================
# Test Surface Cell Detection
# ============================================================================

@testset "Surface Cell Detection" begin

	@testset "find_surface_cell" begin
		mesh = create_test_mesh(n_cells = 10)

		surface_id = find_surface_cell(mesh)

		# Surface cell should have the highest z center
		@test surface_id == 10  # Last cell
		@test mesh.cell_centers[surface_id] == maximum(mesh.cell_centers)
	end

	@testset "find_bottom_cell" begin
		mesh = create_test_mesh(n_cells = 10)

		bottom_id = find_bottom_cell(mesh)

		# Bottom cell should have the lowest z center
		@test bottom_id == 1  # First cell
		@test mesh.cell_centers[bottom_id] == minimum(mesh.cell_centers)
	end

	@testset "find_interior_neighbor" begin
		mesh = create_test_mesh(n_cells = 10)

		# Surface cell's interior neighbor is cell 9
		surface_id = find_surface_cell(mesh)
		neighbor_id = find_interior_neighbor(mesh, surface_id)
		@test neighbor_id == 9

		# Bottom cell has no interior neighbor
		bottom_id = find_bottom_cell(mesh)
		@test find_interior_neighbor(mesh, bottom_id) == 0
	end

	@testset "find_exterior_neighbor" begin
		mesh = create_test_mesh(n_cells = 10)

		# Cell 5's exterior neighbor is cell 6
		neighbor_id = find_exterior_neighbor(mesh, 5)
		@test neighbor_id == 6

		# Surface cell has no exterior neighbor
		surface_id = find_surface_cell(mesh)
		@test find_exterior_neighbor(mesh, surface_id) == 0
	end
end

# ============================================================================
# Test Adjacent Cell Merging
# ============================================================================

@testset "Adjacent Cell Merging" begin

	@testset "find_shared_face" begin
		mesh = create_test_mesh(n_cells = 10)

		# Adjacent cells share a face
		face_id = find_shared_face(mesh, 5, 6)
		@test face_id > 0

		# Non-adjacent cells don't share a face
		@test find_shared_face(mesh, 1, 5) == 0
	end

	@testset "are_cells_adjacent" begin
		mesh = create_test_mesh(n_cells = 10)

		@test are_cells_adjacent(mesh, 5, 6) == true
		@test are_cells_adjacent(mesh, 6, 5) == true  # Symmetric
		@test are_cells_adjacent(mesh, 1, 10) == false
	end

	@testset "merge_adjacent_cells! basic" begin
		mesh = create_test_mesh(n_cells = 10)

		# Record pre-merge state
		n_cells_before = mesh.n_active_cells
		V_9 = mesh.cell_volumes[9]
		V_10 = mesh.cell_volumes[10]
		T_9 = mesh.cell_states[9].T
		T_10 = mesh.cell_states[10].T
		ξ1_9 = mesh.cell_states[9].ξ[1]
		ξ1_10 = mesh.cell_states[10].ξ[1]

		# Merge surface cell (10) into neighbor (9)
		merged_id = merge_adjacent_cells!(mesh, 10, 9)

		# Check result
		@test merged_id == 9
		@test mesh.n_active_cells == n_cells_before - 1
		@test mesh.active_cell_mask[10] == false
		@test mesh.active_cell_mask[9] == true

		# Check conservation
		V_merged = mesh.cell_volumes[9]
		@test V_merged ≈ V_9 + V_10 atol=1e-10

		# Temperature should be volume-weighted average
		T_expected = (V_9 * T_9 + V_10 * T_10) / V_merged
		@test mesh.cell_states[9].T ≈ T_expected atol=1e-10

		# Concentration should be conservative
		ξ1_expected = (V_9 * ξ1_9 + V_10 * ξ1_10) / V_merged
		@test mesh.cell_states[9].ξ[1] ≈ ξ1_expected atol=1e-10
	end

	@testset "merge_adjacent_cells! connectivity" begin
		mesh = create_test_mesh(n_cells = 10)

		# Merge surface cell into neighbor
		merged_id = merge_adjacent_cells!(mesh, 10, 9)

		# The merged cell should now be at the boundary
		cell = mesh.cells[merged_id]
		right_face = mesh.faces[cell.face_ids[2]]

		@test right_face.is_boundary == true
		@test right_face.boundary_id == 1  # Top boundary
		@test right_face.cell_ids[1] == merged_id
		@test right_face.cell_ids[2] == 0  # No neighbor (boundary)
	end

	@testset "merge_surface_cell!" begin
		mesh = create_test_mesh(n_cells = 10)
		n_before = mesh.n_active_cells

		merged_id = merge_surface_cell!(mesh)

		@test merged_id == 9  # Surface was 10, merged into 9
		@test mesh.n_active_cells == n_before - 1
	end

	@testset "multiple merges" begin
		mesh = create_test_mesh(n_cells = 10)

		# Merge surface cell twice
		merge_surface_cell!(mesh)
		@test mesh.n_active_cells == 9

		merge_surface_cell!(mesh)
		@test mesh.n_active_cells == 8

		# Check new surface cell is correct
		surface_id = find_surface_cell(mesh)
		@test surface_id == 8
	end
end

# ============================================================================
# Test State Vector Resizing
# ============================================================================

@testset "State Vector Resizing" begin

	@testset "remove_surface_cell_from_state_vector" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		# Record original values
		T_orig = [get_temperature_ale(u, sv, i) for i in 1:sv.n_cells]
		z_orig = [get_node_position(u, sv, i) for i in 1:sv.n_nodes]

		# Remove surface cell
		u_new, sv_new = remove_surface_cell_from_state_vector(u, sv)

		# Check sizes
		@test sv_new.n_cells == sv.n_cells - 1
		@test sv_new.n_nodes == sv.n_nodes - 1
		@test length(u_new) == state_vector_length_ale(sv_new)

		# Check temperatures (should be same except last)
		for i in 1:sv_new.n_cells
			@test get_temperature_ale(u_new, sv_new, i) ≈ T_orig[i]
		end

		# Check node positions (should be same except last)
		for i in 1:sv_new.n_nodes
			@test get_node_position(u_new, sv_new, i) ≈ z_orig[i]
		end
	end

	@testset "remove_cell_from_state_vector" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		# Record original values
		T_orig = [get_temperature_ale(u, sv, i) for i in 1:sv.n_cells]

		# Remove cell 5 and node 6
		u_new, sv_new = remove_cell_from_state_vector(u, sv, 5, 6)

		# Check sizes
		@test sv_new.n_cells == sv.n_cells - 1
		@test length(u_new) == state_vector_length_ale(sv_new)

		# Check temperatures: should have all except cell 5
		@test get_temperature_ale(u_new, sv_new, 1) ≈ T_orig[1]
		@test get_temperature_ale(u_new, sv_new, 4) ≈ T_orig[4]
		@test get_temperature_ale(u_new, sv_new, 5) ≈ T_orig[6]  # Shifted
		@test get_temperature_ale(u_new, sv_new, 9) ≈ T_orig[10]
	end

	@testset "create_merged_state_vector" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		V_9 = mesh.cell_volumes[9]
		V_10 = mesh.cell_volumes[10]
		T_9 = get_temperature_ale(u, sv, 9)
		T_10 = get_temperature_ale(u, sv, 10)

		# Merge cells 10 into 9
		u_new, sv_new = create_merged_state_vector(u, sv, 10, 9, V_10, V_9)

		# Check size
		@test sv_new.n_cells == sv.n_cells - 1

		# Check merged temperature
		V_total = V_9 + V_10
		T_expected = (V_9 * T_9 + V_10 * T_10) / V_total
		@test get_temperature_ale(u_new, sv_new, 9) ≈ T_expected
	end

	@testset "map_cell_indices_after_removal" begin
		mapping = map_cell_indices_after_removal(10, 5)

		@test mapping[1] == 1
		@test mapping[4] == 4
		@test !haskey(mapping, 5)  # Removed
		@test mapping[6] == 5  # Shifted
		@test mapping[10] == 9
	end
end

# ============================================================================
# Test Conservation During Merge
# ============================================================================

@testset "Conservation During Merge" begin

	@testset "mass conservation" begin
		mesh = create_depleted_surface_mesh()

		# Compute total mass before merge
		mass_before = 0.0
		for i in 1:mesh.n_active_cells
			if mesh.active_cell_mask[i]
				V = mesh.cell_volumes[i]
				ξ_total = sum(mesh.cell_states[i].ξ)
				mass_before += V * ξ_total
			end
		end

		# Merge surface cell
		merge_surface_cell!(mesh)

		# Compute total mass after merge
		mass_after = 0.0
		for i in eachindex(mesh.cells)
			if mesh.active_cell_mask[i]
				V = mesh.cell_volumes[i]
				ξ_total = sum(mesh.cell_states[i].ξ)
				mass_after += V * ξ_total
			end
		end

		@test mass_after ≈ mass_before rtol=1e-10
	end

	@testset "energy conservation (proxy)" begin
		mesh = create_test_mesh(n_cells = 10)

		# Compute thermal energy proxy before merge: Σ V_i * T_i
		energy_before = sum(mesh.cell_volumes[i] * mesh.cell_states[i].T
							for i in 1:10 if mesh.active_cell_mask[i])

		# Merge surface cell
		merge_surface_cell!(mesh)

		# Compute thermal energy proxy after merge
		energy_after = sum(mesh.cell_volumes[i] * mesh.cell_states[i].T
						   for i in eachindex(mesh.cells) if mesh.active_cell_mask[i])

		@test energy_after ≈ energy_before rtol=1e-10
	end
end

# ============================================================================
# Test Error Handling
# ============================================================================

@testset "Error Handling" begin

	@testset "merge non-adjacent cells" begin
		mesh = create_test_mesh(n_cells = 10)
		@test_throws ErrorException merge_adjacent_cells!(mesh, 1, 5)
	end

	@testset "merge inactive cell" begin
		mesh = create_test_mesh(n_cells = 10)
		mesh.active_cell_mask[5] = false
		@test_throws ErrorException merge_adjacent_cells!(mesh, 5, 6)
	end

	@testset "remove from single-cell state vector" begin
		mesh = create_test_mesh(n_cells = 1)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)
		@test_throws ErrorException remove_surface_cell_from_state_vector(u, sv)
	end
end

# ============================================================================
# Integration with ODE Solver
# ============================================================================

using Pyrolysis.Solvers: create_depletion_callback, DepletionCallbackState
using Pyrolysis.Properties

@testset "Depletion Callback Structure" begin

	@testset "DepletionCallbackState creation" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		cb_state = DepletionCallbackState(sv, mesh)

		@test cb_state.n_merges == 0
		@test isempty(cb_state.merge_times)
		@test cb_state.sv === sv
		@test cb_state.initial_cell_thickness > 0.0  # Should capture initial thickness
	end

	@testset "callback creation" begin
		mesh = create_test_mesh(n_cells = 10)
		material = create_test_material()
		sv = StateVectorALE(mesh)
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50000.0),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		# Create ODE function
		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		# Create callback with thickness fraction threshold (5% of initial)
		cb = create_depletion_callback(ode_fn, mesh, material, sv;
			threshold = 0.05, min_cells = 2)

		@test cb isa DiscreteCallback
	end
end

# ============================================================================
# Tests for Thickness-Based Depletion Detection
# ============================================================================

using Pyrolysis.Solvers: is_cell_thin, get_cell_thickness, get_surface_cell_thickness,
	find_thin_cells

@testset "Thickness-Based Depletion" begin

	@testset "get_cell_thickness" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		# Check that we can get thickness for each cell
		for i in 1:10
			thickness = get_cell_thickness(u, sv, i)
			@test thickness > 0.0
			# For uniform mesh of L=0.01 with 10 cells, each cell is 0.001
			@test thickness ≈ 0.001 atol=1e-10
		end
	end

	@testset "get_surface_cell_thickness" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		thickness = get_surface_cell_thickness(u, sv)
		@test thickness ≈ 0.001 atol=1e-10  # 0.01 / 10 cells
	end

	@testset "is_cell_thin with uniform mesh" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		initial_thickness = 0.001  # Each cell is 0.01/10
		surface_cell = 10

		# Cell is not thin when at initial thickness
		@test !is_cell_thin(u, sv, surface_cell, initial_thickness; threshold = 0.05)
		@test !is_cell_thin(u, sv, surface_cell, initial_thickness; threshold = 0.5)

		# Cell would be thin if threshold was > 1.0 (impossible)
		@test is_cell_thin(u, sv, surface_cell, initial_thickness; threshold = 1.1)
	end

	@testset "is_cell_thin with shrunk mesh" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		initial_thickness = 0.001
		surface_cell = 10

		# Shrink the surface cell by moving the top node down
		# Original: node 11 at z=0.01, node 10 at z=0.009
		# Node positions start at index sv.z_offset + 1 (1-based)
		z_start = sv.z_offset

		# Move top node (node 11) to make surface cell 10% of original
		# New cell thickness should be 0.0001 (10% of 0.001)
		# u[z_start + 11] is node 11 (the top surface node)
		u[z_start+11] = 0.0091  # 0.0091 - 0.009 = 0.0001 thickness

		# Now cell should be thin at 50% threshold but not at 5%
		@test is_cell_thin(u, sv, surface_cell, initial_thickness; threshold = 0.5)
		@test is_cell_thin(u, sv, surface_cell, initial_thickness; threshold = 0.15)
		@test !is_cell_thin(u, sv, surface_cell, initial_thickness; threshold = 0.05)
	end

	@testset "find_thin_cells with no thin cells" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		initial_thickness = 0.001  # Each cell is 0.01/10

		# With uniform mesh at initial state, no cells should be thin
		thin_cells = find_thin_cells(u, sv, initial_thickness, 0.05)
		@test isempty(thin_cells)
	end

	@testset "find_thin_cells with one thin cell" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		initial_thickness = 0.001
		z_start = sv.z_offset

		# Shrink surface cell (cell 10) to 4% of original
		# Node 11 at z=0.01, node 10 at z=0.009 -> make cell 0.00004 thick
		u[z_start+11] = 0.00904  # 0.00904 - 0.009 = 0.00004

		thin_cells = find_thin_cells(u, sv, initial_thickness, 0.05)
		@test length(thin_cells) == 1
		@test thin_cells[1] == 10  # Surface cell
	end

	@testset "find_thin_cells with multiple thin cells sorted surface-first" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		initial_thickness = 0.001
		z_start = sv.z_offset

		# Shrink cells 3, 7, and 10 to be thin
		# Cell i spans from node i to node i+1
		# To shrink cell i, move node i+1 closer to node i
		# Note: changing node i+1 affects both cell i and cell i+1

		# Original node positions: node i at z = (i-1) * 0.001

		# Make cell 3 thin: move node 4 from 0.003 to 0.00204 (cell 3 = 4e-5)
		u[z_start+4] = 0.00204

		# Make cell 7 thin: move node 8 from 0.007 to 0.00604 (cell 7 = 4e-5)
		# But node 8 is now affected by changes above, so use absolute position
		u[z_start+8] = 0.00604

		# Make cell 10 thin: move node 11 from 0.01 to 0.00904 (cell 10 = 4e-5)
		u[z_start+11] = 0.00904

		thin_cells = find_thin_cells(u, sv, initial_thickness, 0.05)

		# Should find 3 thin cells, sorted by z descending (10, 7, 3)
		@test length(thin_cells) == 3
		@test thin_cells[1] == 10  # Highest z (surface)
		@test thin_cells[2] == 7   # Middle
		@test thin_cells[3] == 3   # Lowest z
	end

	@testset "find_thin_cells respects threshold" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		initial_thickness = 0.001
		z_start = sv.z_offset

		# Shrink surface cell to 10% of original
		u[z_start+11] = 0.0091  # 0.0091 - 0.009 = 0.0001 (10%)

		# At 5% threshold, cell is not thin (10% > 5%)
		thin_cells = find_thin_cells(u, sv, initial_thickness, 0.05)
		@test isempty(thin_cells)

		# At 15% threshold, cell is thin (10% < 15%)
		thin_cells = find_thin_cells(u, sv, initial_thickness, 0.15)
		@test length(thin_cells) == 1
		@test thin_cells[1] == 10
	end
end

@testset "resize_ode_function_ale!" begin
	mesh = create_test_mesh(n_cells = 10)
	material = create_test_material()
	bc_set = BoundaryConditionSet(
		thermal_top = HeatFluxBC(50000.0),
		thermal_bottom = AdiabaticBC(),
		mass_top = ImpermeableBC(),
		mass_bottom = ImpermeableBC(),
	)

	ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

	# Initial sizes
	@test length(ode_fn.dT_work) == 10
	@test size(ode_fn.dξ_work) == (2, 10)
	@test length(ode_fn.w_work) == 11  # n_nodes = n_cells + 1
	# NEW: Cell-centered advection arrays (not face-based flux arrays)
	@test length(ode_fn.ale_advection) == 10

	# Resize to smaller
	resize_ode_function_ale!(ode_fn, 8, 9)

	@test length(ode_fn.dT_work) == 8
	@test size(ode_fn.dξ_work) == (2, 8)
	@test length(ode_fn.w_work) == 9
	@test length(ode_fn.ale_advection) == 8
	@test size(ode_fn.ale_species_advection) == (2, 8)
	@test ode_fn.sv.n_cells == 8
end

@testset "resize_geometry_cache!" begin
	cache = GeometryCache(10, 11)

	@test length(cache.volumes) == 10
	@test length(cache.face_positions) == 11

	resize_geometry_cache!(cache, 8, 9)

	@test length(cache.volumes) == 8
	@test length(cache.centers) == 8
	@test length(cache.face_positions) == 9
	@test length(cache.face_distances) == 9
end

# ============================================================================
# Tests for Priority 1 Fixes: Node Removal Consistency & Energy Conservation
# ============================================================================

@testset "Node Removal Consistency" begin
	@testset "mesh and state vector dimensions match after merge" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		# Initial dimensions
		@test length(mesh.cells) == 10
		@test length(mesh.nodes) == 11
		@test sv.n_cells == 10
		@test sv.n_nodes == 11

		# Merge surface cell
		surface_id = find_surface_cell(mesh)
		neighbor_id = find_interior_neighbor(mesh, surface_id)
		V_surface = mesh.cell_volumes[surface_id]
		V_neighbor = mesh.cell_volumes[neighbor_id]

		# Merge in mesh
		merged_id = merge_adjacent_cells!(mesh, surface_id, neighbor_id)
		@test merged_id == neighbor_id

		# Compact mesh (using general function with surface_id as removed cell)
		compact_mesh_after_merge!(mesh, surface_id)

		# Create new state vector
		u_new, sv_new = create_merged_state_vector(u, sv, surface_id, neighbor_id,
			V_surface, V_neighbor)

		# Verify dimensions match
		@test length(mesh.cells) == sv_new.n_cells
		@test length(mesh.nodes) == sv_new.n_nodes
		@test length(u_new) == state_vector_length_ale(sv_new)

		# Verify expected sizes (9 cells, 10 nodes after merging 1)
		@test length(mesh.cells) == 9
		@test length(mesh.nodes) == 10
		@test sv_new.n_cells == 9
		@test sv_new.n_nodes == 10
	end

	@testset "multiple consecutive merges maintain consistency" begin
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		# Perform 3 consecutive merges
		for merge_num in 1:3
			n_before = length(mesh.cells)
			surface_id = n_before  # Current surface
			neighbor_id = n_before - 1

			V_surface = mesh.cell_volumes[surface_id]
			V_neighbor = mesh.cell_volumes[neighbor_id]

			# Merge in mesh
			merged_id = merge_adjacent_cells!(mesh, surface_id, neighbor_id)
			compact_mesh_after_merge!(mesh, surface_id)

			# Update state vector
			u, sv = create_merged_state_vector(u, sv, surface_id, neighbor_id,
				V_surface, V_neighbor)

			# Verify dimensions match after each merge
			@test length(mesh.cells) == sv.n_cells
			@test length(mesh.nodes) == sv.n_nodes
			@test length(u) == state_vector_length_ale(sv)
		end

		# Final size should be 7 cells
		@test length(mesh.cells) == 7
	end

	@testset "interior cell merge and compaction" begin
		# Test that we can merge and compact an interior cell (not surface)
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		# Merge cell 5 into cell 4 (interior cells)
		cell_to_remove = 5
		cell_to_keep = 4

		V_remove = mesh.cell_volumes[cell_to_remove]
		V_keep = mesh.cell_volumes[cell_to_keep]

		# Merge in mesh
		merged_id = merge_adjacent_cells!(mesh, cell_to_remove, cell_to_keep)
		@test merged_id == cell_to_keep

		# Compact mesh using general function
		compact_mesh_after_merge!(mesh, cell_to_remove)

		# Create new state vector
		u_new, sv_new = create_merged_state_vector(u, sv, cell_to_remove, cell_to_keep,
			V_remove, V_keep)

		# Verify dimensions match
		@test length(mesh.cells) == sv_new.n_cells
		@test length(mesh.nodes) == sv_new.n_nodes
		@test length(u_new) == state_vector_length_ale(sv_new)

		# Verify expected sizes (9 cells, 10 nodes after merging 1)
		@test length(mesh.cells) == 9
		@test length(mesh.nodes) == 10
		@test sv_new.n_cells == 9
		@test sv_new.n_nodes == 10

		# Verify mesh connectivity is still valid
		for i in 1:length(mesh.cells)
			cell = mesh.cells[i]
			@test cell.node_ids[1] >= 1
			@test cell.node_ids[2] <= length(mesh.nodes)
			@test cell.face_ids[1] >= 1
			@test cell.face_ids[2] <= length(mesh.faces)
		end
	end

	@testset "bottom-adjacent cell merge and compaction" begin
		# Test merging cell near bottom boundary
		mesh = create_test_mesh(n_cells = 10)
		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		# Merge cell 2 into cell 3 (cell 2 is next to bottom, merge toward interior)
		cell_to_remove = 2
		cell_to_keep = 3

		V_remove = mesh.cell_volumes[cell_to_remove]
		V_keep = mesh.cell_volumes[cell_to_keep]

		# Merge in mesh
		merged_id = merge_adjacent_cells!(mesh, cell_to_remove, cell_to_keep)
		@test merged_id == cell_to_keep

		# Compact mesh using general function
		compact_mesh_after_merge!(mesh, cell_to_remove)

		# Create new state vector
		u_new, sv_new = create_merged_state_vector(u, sv, cell_to_remove, cell_to_keep,
			V_remove, V_keep)

		# Verify dimensions
		@test length(mesh.cells) == 9
		@test length(mesh.nodes) == 10
		@test sv_new.n_cells == 9

		# Verify bottom cell still exists and has boundary connectivity
		bottom_cell = findfirst(i -> mesh.cells[i].node_ids[1] == 1, 1:length(mesh.cells))
		@test bottom_cell !== nothing

		# Check that boundary faces are still correct
		@test haskey(mesh.boundary_faces, 2)  # Bottom boundary
		@test length(mesh.boundary_faces[2]) >= 1
	end
end

@testset "Energy-Conserving Temperature Merge" begin
	material = create_test_material()

	@testset "identical cells: energy conservation equals volume average" begin
		# Two identical cells should give same result either way
		V_a, V_b = 1e-6, 1e-6  # Equal volumes
		T_a, T_b = 500.0, 500.0  # Same temperature
		ξ_a = (1200.0, 0.0)  # Same composition
		ξ_b = (1200.0, 0.0)

		state_a = CellState{2}(T_a, ξ_a, V_a, V_a)
		state_b = CellState{2}(T_b, ξ_b, V_b, V_b)

		T_merged = merge_temperatures_conservatively(V_a, state_a, V_b, state_b, material)
		T_simple = (V_a * T_a + V_b * T_b) / (V_a + V_b)

		@test T_merged ≈ T_simple rtol=1e-10
		@test T_merged ≈ 500.0 rtol=1e-10
	end

	@testset "different compositions: energy conservation vs volume average" begin
		# Cell A: solid-rich (high ρcₚ)
		# Cell B: gas-rich (low ρcₚ)
		V_a, V_b = 1e-6, 1e-6  # Equal volumes
		T_a = 800.0  # Hot cell
		T_b = 400.0  # Cold cell
		ξ_a = (1200.0, 0.0)    # All solid, ρcₚ ≈ 1200 * 1500 = 1.8e6
		ξ_b = (120.0, 1080.0)  # Mostly gas, ρcₚ ≈ 120 * 1500 + 1080 * 1000 = 1.26e6

		state_a = CellState{2}(T_a, ξ_a, V_a, V_a)
		state_b = CellState{2}(T_b, ξ_b, V_b, V_b)

		T_merged = merge_temperatures_conservatively(V_a, state_a, V_b, state_b, material)
		T_simple = (V_a * T_a + V_b * T_b) / (V_a + V_b)

		# Simple volume average: (800 + 400) / 2 = 600 K
		@test T_simple ≈ 600.0 rtol=1e-10

		# Energy-conserving should be different because cell A has higher ρcₚ
		# and higher temperature, so it contributes more energy
		# T_merged should be > 600 K (closer to the solid-rich cell's temp)
		@test T_merged > T_simple

		# Verify actual energy conservation
		c_solid = 1500.0  # From material definition
		c_gas = 1000.0
		ρcₚ_a = ξ_a[1] * c_solid + ξ_a[2] * c_gas
		ρcₚ_b = ξ_b[1] * c_solid + ξ_b[2] * c_gas
		E_before = V_a * ρcₚ_a * T_a + V_b * ρcₚ_b * T_b

		V_total = V_a + V_b
		ξ_merged = ((V_a * ξ_a[1] + V_b * ξ_b[1]) / V_total,
			(V_a * ξ_a[2] + V_b * ξ_b[2]) / V_total)
		ρcₚ_merged = ξ_merged[1] * c_solid + ξ_merged[2] * c_gas
		E_after = V_total * ρcₚ_merged * T_merged

		@test E_after ≈ E_before rtol=0.01  # Within 1% (due to T-dependent cₚ approximation)
	end

	@testset "state vector merge uses energy conservation when material provided" begin
		mesh = create_test_mesh(n_cells = 5)
		material = create_test_material()

		# Set up cells with different compositions
		mesh.cell_states[5] = CellState{2}(800.0, (1200.0, 0.0), mesh.cell_volumes[5], mesh.cell_volumes[5])
		mesh.cell_states[4] = CellState{2}(400.0, (120.0, 1080.0), mesh.cell_volumes[4], mesh.cell_volumes[4])

		sv = StateVectorALE(mesh)
		u = create_state_vector_ale(mesh)

		V_5 = mesh.cell_volumes[5]
		V_4 = mesh.cell_volumes[4]

		# Merge WITH material (energy-conserving) - use a copy of u
		u_copy1 = copy(u)
		u_energy, sv_energy = create_merged_state_vector(u_copy1, sv, 5, 4, V_5, V_4;
			material = material)

		# Merge WITHOUT material (volume-weighted) - use another copy
		u_copy2 = copy(u)
		u_simple, sv_simple = create_merged_state_vector(u_copy2, sv, 5, 4, V_5, V_4)

		# Get merged temperatures
		T_energy = get_temperature_ale(u_energy, sv_energy, 4)
		T_simple = get_temperature_ale(u_simple, sv_simple, 4)

		# Simple volume average should be (800 + 400) / 2 = 600
		@test T_simple ≈ 600.0 rtol=1e-10

		# Energy-conserving should differ from simple average
		@test T_energy != T_simple

		# Energy-conserving should be higher (solid-rich cell was hotter and had higher ρcₚ)
		@test T_energy > T_simple
	end
end
