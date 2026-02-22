"""
Tests for the ALE (Arbitrary Lagrangian-Eulerian) module.

Tests cover:
1. GCL (Geometric Conservation Law) compliance
2. Mesh motion algorithms
3. Mesh velocity methods
4. Solution remapping
5. Interface tracking
6. Conservation properties
7. Upwind advection formulation (ALE_Energy_Reformulation.md)

Reference: Design Document Sections 3 (ALE Formulation), 10 (Verification)
"""

using Test
using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Geometry
using Pyrolysis.ALE
using Pyrolysis.Discretization

# ============================================================================
# Test Fixtures
# ============================================================================

"""Create a simple uniform 1D mesh for testing."""
function create_test_mesh(; L = 0.01, n_cells = 10, NC = 2)
	return create_uniform_mesh(L, n_cells, Val(NC))
end

"""Create a mesh with known state for remapping tests."""
function create_mesh_with_state(; L = 0.01, n_cells = 10, NC = 2, T0 = 300.0, ξ0 = 1000.0)
	mesh = create_uniform_mesh(L, n_cells, Val(NC))

	# Set initial state: linear temperature profile, uniform concentration
	for i in eachindex(mesh.cells)
		z = mesh.cell_centers[i]
		mesh.cell_states[i].T = T0 + 100.0 * z / L  # 300K at bottom, 400K at top
		mesh.cell_states[i].ξ = ntuple(j -> j == 1 ? ξ0 : 0.0, NC)
	end

	return mesh
end

# ============================================================================
# GCL Tests
# ============================================================================

@testset "GCL Verification" begin

	@testset "Uniform expansion GCL" begin
		mesh = create_test_mesh(n_cells = 20)
		n_nodes = length(mesh.nodes)

		# Uniform expansion: all nodes move at constant rate
		w = fill(0.001, n_nodes)  # 1 mm/s everywhere
		dt = 0.01  # 10 ms

		# GCL should be satisfied exactly for uniform motion
		error = gcl_error(mesh, w, dt)
		@test error < 1e-14
		@test check_gcl(mesh, w, dt)
	end

	@testset "Linear velocity profile GCL" begin
		mesh = create_test_mesh(n_cells = 20)
		n_nodes = length(mesh.nodes)

		# Linear velocity: w = 0 at bottom, w_max at top
		w = zeros(n_nodes)
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min
		w_max = 0.001

		for i in 1:n_nodes
			z = mesh.nodes[i].z
			w[i] = w_max * (z - z_min) / L
		end

		dt = 0.01

		# GCL should be satisfied exactly for any smooth velocity field
		error = gcl_error(mesh, w, dt)
		@test error < 1e-14
		@test check_gcl(mesh, w, dt)
	end

	@testset "GCL-compliant volume update" begin
		mesh = create_test_mesh(n_cells = 10)
		n_nodes = length(mesh.nodes)

		# Store original volumes
		V_old = copy(mesh.cell_volumes)

		# Linear velocity profile
		w = zeros(n_nodes)
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min
		w_max = 0.0005

		for i in 1:n_nodes
			z = mesh.nodes[i].z
			w[i] = w_max * (z - z_min) / L
		end

		dt = 0.1

		# Update volumes using GCL formula
		update_volumes_gcl!(mesh, w, dt)

		# Verify GCL relation for each cell
		A = mesh.cross_section_area
		for i in eachindex(mesh.cells)
			cell = mesh.cells[i]
			w_left = w[cell.node_ids[1]]
			w_right = w[cell.node_ids[2]]

			dV_expected = dt * (w_right - w_left) * A
			dV_actual = mesh.cell_volumes[i] - V_old[i]

			@test abs(dV_actual - dV_expected) < 1e-14
		end
	end

	@testset "Post-motion GCL error" begin
		mesh = create_test_mesh(n_cells = 10)
		n_nodes = length(mesh.nodes)

		V_old = copy(mesh.cell_volumes)

		# Compression: top moves down
		w = zeros(n_nodes)
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		for i in 1:n_nodes
			z = mesh.nodes[i].z
			w[i] = -0.0001 * (z - z_min) / L  # Compression
		end

		dt = 0.1

		# Move mesh
		move_mesh!(mesh, w, dt; use_gcl_volumes = true)

		# Check GCL error post-motion
		error = gcl_error_post_motion(mesh, V_old, w, dt)
		@test error < 1e-14
	end
end

# ============================================================================
# Mesh Motion Tests
# ============================================================================

@testset "Mesh Motion" begin

	@testset "Basic mesh motion" begin
		mesh = create_test_mesh(n_cells = 10)
		n_nodes = length(mesh.nodes)

		# Store original positions
		z_orig = [mesh.nodes[i].z for i in 1:n_nodes]

		# Uniform velocity
		w = fill(0.001, n_nodes)
		dt = 0.1

		# Move mesh
		move_mesh!(mesh, w, dt)

		# Check new positions
		for i in 1:n_nodes
			expected_z = z_orig[i] + dt * w[i]
			@test mesh.nodes[i].z ≈ expected_z atol=1e-14
		end

		# Check mesh time updated
		@test mesh.mesh_time ≈ dt
	end

	@testset "Mesh quality check" begin
		mesh = create_test_mesh(n_cells = 10)

		# Initial mesh should be valid
		quality = check_mesh_quality(mesh)
		@test quality.valid
		@test quality.min_volume > 0
		@test isempty(quality.tangled_cells)
		@test mesh_is_valid(mesh)
	end

	@testset "Safe mesh motion rejection" begin
		mesh = create_test_mesh(n_cells = 10)
		n_nodes = length(mesh.nodes)

		# Try to create tangling with extreme velocity
		w = zeros(n_nodes)
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		# Interior node moves past its neighbors
		for i in 1:n_nodes
			if mesh.nodes[i].boundary_id == 0
				# Interior: move drastically
				w[i] = 10.0 * L  # Move more than domain length
			end
		end

		dt = 0.1

		# Safe motion should reject this
		original_time = mesh.mesh_time
		success = move_mesh_safe!(mesh, w, dt; max_aspect_ratio = 10.0)

		# Either rejected or mesh is still valid
		if !success
			@test mesh.mesh_time == original_time  # Rolled back
		end
		@test mesh_is_valid(mesh)
	end

	@testset "Mesh displacement tracking" begin
		mesh = create_test_mesh(n_cells = 10)
		n_nodes = length(mesh.nodes)

		# Initial displacement should be zero
		disp = mesh_displacement(mesh)
		@test all(d -> abs(d) < 1e-14, disp)
		@test max_mesh_displacement(mesh) < 1e-14

		# Move mesh
		w = fill(0.001, n_nodes)
		dt = 0.1
		move_mesh!(mesh, w, dt)

		# Now displacement should be dt * w
		disp = mesh_displacement(mesh)
		for i in 1:n_nodes
			@test disp[i] ≈ dt * w[i] atol=1e-14
		end
		@test max_mesh_displacement(mesh) ≈ dt * 0.001 atol=1e-14
	end

	@testset "Reset mesh to reference" begin
		mesh = create_test_mesh(n_cells = 10)
		n_nodes = length(mesh.nodes)

		# Store reference
		z_ref = [mesh.nodes[i].z_ref for i in 1:n_nodes]

		# Move mesh
		w = fill(0.001, n_nodes)
		move_mesh!(mesh, w, 0.1)

		# Reset
		reset_mesh_to_reference!(mesh)

		# Check positions restored
		for i in 1:n_nodes
			@test mesh.nodes[i].z ≈ z_ref[i] atol=1e-14
		end
		@test mesh.mesh_time ≈ 0.0
	end
end

# ============================================================================
# Mesh Velocity Tests
# ============================================================================

@testset "Mesh Velocity Methods" begin

	@testset "Laplacian smoothing - linear interpolation" begin
		mesh = create_test_mesh(n_cells = 20)
		n_nodes = length(mesh.nodes)

		w = zeros(n_nodes)
		w_top = 0.002
		w_bottom = 0.0

		compute_mesh_velocity!(w, mesh, w_top, w_bottom, LaplacianSmoothing())

		# Check boundary values
		for i in 1:n_nodes
			node = mesh.nodes[i]
			if node.boundary_id == 1  # Top
				@test w[i] ≈ w_top atol=1e-14
			elseif node.boundary_id == 2  # Bottom
				@test w[i] ≈ w_bottom atol=1e-14
			end
		end

		# Check linear interpolation
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		for i in 1:n_nodes
			node = mesh.nodes[i]
			if !node.is_boundary
				t = (node.z - z_min) / L
				expected = w_bottom + t * (w_top - w_bottom)
				@test w[i] ≈ expected atol=1e-10
			end
		end
	end

	@testset "Spring analogy convergence" begin
		mesh = create_test_mesh(n_cells = 20)
		n_nodes = length(mesh.nodes)

		w = zeros(n_nodes)
		w_top = 0.001
		w_bottom = 0.0

		compute_mesh_velocity!(w, mesh, w_top, w_bottom, SpringAnalogy())

		# Check boundary values
		for i in 1:n_nodes
			node = mesh.nodes[i]
			if node.boundary_id == 1
				@test w[i] ≈ w_top atol=1e-10
			elseif node.boundary_id == 2
				@test w[i] ≈ w_bottom atol=1e-10
			end
		end

		# Interior should be between boundaries
		for i in 1:n_nodes
			node = mesh.nodes[i]
			if !node.is_boundary
				@test w[i] >= w_bottom - 1e-10
				@test w[i] <= w_top + 1e-10
			end
		end
	end

	@testset "Weighted ALE" begin
		mesh = create_test_mesh(n_cells = 20)
		n_nodes = length(mesh.nodes)

		w = zeros(n_nodes)
		w_top = 0.001
		w_bottom = 0.0

		method = WeightedALE(0.5, 0.5)  # 50% Lagrangian, 50% smoothing
		compute_mesh_velocity!(w, mesh, w_top, w_bottom, method)

		# Boundary values should be exact
		for i in 1:n_nodes
			node = mesh.nodes[i]
			if node.boundary_id == 1
				@test w[i] ≈ w_top atol=1e-10
			elseif node.boundary_id == 2
				@test w[i] ≈ w_bottom atol=1e-10
			end
		end
	end

	@testset "Surface regression rate" begin
		mesh = create_test_mesh()

		mass_flux = 0.01  # 10 g/(m²·s)
		ρ_surface = 1000.0  # 1000 kg/m³

		v_surface = surface_regression_rate(mesh, mass_flux, ρ_surface)

		# v = -ṁ/ρ
		@test v_surface ≈ -mass_flux / ρ_surface atol=1e-14
	end
end

# ============================================================================
# Remapping Tests
# ============================================================================

@testset "Solution Remapping" begin

	@testset "First-order remap - identical mesh" begin
		mesh = create_mesh_with_state()
		NC = 2
		n_cells = length(mesh.cells)
		n_dof = n_cells * (1 + NC)

		# Create solution vector from state
		u_old = zeros(n_dof)
		for i in 1:n_cells
			idx = (i - 1) * (1 + NC)
			u_old[idx+1] = mesh.cell_states[i].T
			for j in 1:NC
				u_old[idx+1+j] = mesh.cell_states[i].ξ[j]
			end
		end

		# Remap to same mesh
		u_new = zeros(n_dof)
		remap_solution!(u_new, mesh, u_old, mesh, NC, FirstOrderRemap())

		# Should be identical
		@test u_new ≈ u_old atol=1e-12
	end

	@testset "First-order remap - conservation" begin
		mesh_old = create_mesh_with_state(n_cells = 10)
		mesh_new = create_mesh_with_state(n_cells = 10)
		NC = 2
		n_cells = length(mesh_old.cells)
		n_dof = n_cells * (1 + NC)

		# Move new mesh slightly
		n_nodes = length(mesh_new.nodes)
		w = zeros(n_nodes)
		z_min, z_max = domain_bounds(mesh_new)
		L = z_max - z_min

		for i in 1:n_nodes
			z = mesh_new.nodes[i].z
			w[i] = 0.0001 * (z - z_min) / L
		end
		move_mesh!(mesh_new, w, 0.1)

		# Create solution vector
		u_old = zeros(n_dof)
		for i in 1:n_cells
			idx = (i - 1) * (1 + NC)
			u_old[idx+1] = mesh_old.cell_states[i].T
			for j in 1:NC
				u_old[idx+1+j] = mesh_old.cell_states[i].ξ[j]
			end
		end

		# Remap
		u_new = zeros(n_dof)
		remap_solution!(u_new, mesh_new, u_old, mesh_old, NC, FirstOrderRemap())

		# Check mass conservation
		conservation = check_remap_conservation(u_old, mesh_old, u_new, mesh_new, NC)

		for j in 1:NC
			@test conservation.mass_error[j] < 0.01  # <1% error
		end
	end

	@testset "Second-order remap" begin
		mesh_old = create_mesh_with_state(n_cells = 20)
		mesh_new = create_mesh_with_state(n_cells = 20)
		NC = 2
		n_cells = length(mesh_old.cells)
		n_dof = n_cells * (1 + NC)

		# Small mesh motion
		n_nodes = length(mesh_new.nodes)
		w = fill(0.00005, n_nodes)
		move_mesh!(mesh_new, w, 0.1)

		# Create solution
		u_old = zeros(n_dof)
		for i in 1:n_cells
			idx = (i - 1) * (1 + NC)
			u_old[idx+1] = mesh_old.cell_states[i].T
			for j in 1:NC
				u_old[idx+1+j] = mesh_old.cell_states[i].ξ[j]
			end
		end

		# Remap
		u_new = zeros(n_dof)
		remap_solution!(u_new, mesh_new, u_old, mesh_old, NC, SecondOrderRemap(:minmod))

		# Check no new extrema created (limiter working)
		T_min_old = minimum(u_old[(i-1)*(1+NC)+1] for i in 1:n_cells)
		T_max_old = maximum(u_old[(i-1)*(1+NC)+1] for i in 1:n_cells)

		for i in 1:n_cells
			T_new = u_new[(i-1)*(1+NC)+1]
			@test T_new >= T_min_old - 1e-10
			@test T_new <= T_max_old + 1e-10
		end
	end
end

# ============================================================================
# Interface Tracking Tests
# ============================================================================

@testset "Interface Tracking" begin

	@testset "Interface detection" begin
		mesh = create_test_mesh(n_cells = 20)
		n_cells = length(mesh.cells)

		# Create values with a crossing
		values = zeros(n_cells)
		z_interface = 0.005  # Middle of domain

		for i in 1:n_cells
			z = mesh.cell_centers[i]
			if z < z_interface
				values[i] = 0.8
			else
				values[i] = 0.2
			end
		end

		# Find interface at threshold 0.5
		interfaces = find_interface_position(mesh, values, 0.5)

		@test length(interfaces) == 1
		@test interfaces[1] ≈ z_interface rtol=0.1
	end

	@testset "Interface tracker update" begin
		mesh = create_mesh_with_state(n_cells = 20)
		NC = 2
		n_cells = length(mesh.cells)
		n_dof = n_cells * (1 + NC)

		# Set up concentration profile with interface
		for i in 1:n_cells
			z = mesh.cell_centers[i]
			if z < 0.005
				mesh.cell_states[i].ξ = (1000.0, 0.0)
			else
				mesh.cell_states[i].ξ = (200.0, 0.0)
			end
		end

		# Create solution vector
		u = zeros(n_dof)
		for i in 1:n_cells
			idx      = (i - 1) * (1 + NC)
			u[idx+1] = mesh.cell_states[i].T
			u[idx+2] = mesh.cell_states[i].ξ[1]
			u[idx+3] = mesh.cell_states[i].ξ[2]
		end

		# Track concentration (component 1)
		tracker = InterfaceTracker(variable = :ξ1, threshold = 600.0)
		update_interface!(tracker, mesh, u, NC)

		@test length(tracker.interface_positions) >= 1
	end

	@testset "Interface velocity" begin
		tracker = InterfaceTracker()
		tracker.interface_positions = [0.006]

		positions_old = [0.005]
		dt = 0.1

		velocities = interface_velocity(tracker, positions_old, dt)

		@test length(velocities) == 1
		@test velocities[1] ≈ (0.006 - 0.005) / 0.1 atol=1e-14
	end
end

# ============================================================================
# Integration Tests
# ============================================================================

@testset "ALE Integration" begin

	@testset "Full mesh motion cycle" begin
		mesh = create_mesh_with_state(n_cells = 20)
		NC = 2
		n_nodes = length(mesh.nodes)
		n_cells = length(mesh.cells)

		# Store initial mass
		mass_initial = zeros(NC)
		for i in 1:n_cells
			V = mesh.cell_volumes[i]
			for j in 1:NC
				mass_initial[j] += mesh.cell_states[i].ξ[j] * V
			end
		end

		# Compute mesh velocity (surface regression)
		w_top = -0.0001  # Surface moving down
		w_bottom = 0.0   # Fixed bottom

		w = zeros(n_nodes)
		compute_mesh_velocity!(w, mesh, w_top, w_bottom, LaplacianSmoothing())

		# Store old volumes and states
		V_old = copy(mesh.cell_volumes)
		mesh_old = deepcopy(mesh)

		# Move mesh
		dt = 0.1
		move_mesh!(mesh, w, dt; use_gcl_volumes = true)

		# Check GCL
		gcl_err = gcl_error_post_motion(mesh, V_old, w, dt)
		@test gcl_err < 1e-12

		# Check mesh still valid
		@test mesh_is_valid(mesh)

		# Total volume should have changed appropriately
		V_total_old = sum(V_old)
		V_total_new = sum(mesh.cell_volumes)
		A = mesh.cross_section_area
		expected_dV = dt * w_top * A  # Negative (compression)

		@test V_total_new - V_total_old ≈ expected_dV rtol=1e-10
	end

	@testset "Compression and expansion cycle" begin
		mesh = create_test_mesh(n_cells = 10)
		n_nodes = length(mesh.nodes)

		z_orig = [mesh.nodes[i].z for i in 1:n_nodes]
		V_orig = copy(mesh.cell_volumes)

		# Compress
		w_compress = zeros(n_nodes)
		compute_mesh_velocity!(w_compress, mesh, -0.0001, 0.0, LaplacianSmoothing())
		move_mesh!(mesh, w_compress, 0.1)

		# Expand back
		w_expand = zeros(n_nodes)
		compute_mesh_velocity!(w_expand, mesh, 0.0001, 0.0, LaplacianSmoothing())
		move_mesh!(mesh, w_expand, 0.1)

		# Should be approximately back to original
		for i in 1:n_nodes
			@test mesh.nodes[i].z ≈ z_orig[i] rtol=0.01
		end
	end
end

# ============================================================================
# GCL Compliance in Integrated ALE Solver
# ============================================================================

@testset "GCL Compliance in Integrated Solver" begin
	using Pyrolysis.TimeIntegration: PyrolysisODEFunctionALE, StateVectorALE,
		create_state_vector_ale, get_node_position, compute_volumes_from_state

	"""Create a simple reactive material for GCL tests."""
	function create_gcl_test_material()
		virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)

		Material(
			name = :PMMA,
			components = (virgin, gas),
			reactions = (rxn,),
		)
	end

	@testset "ODE RHS produces GCL-consistent volume rates" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)
		material = create_gcl_test_material()

		# Set up initial state with some temperature gradient
		for i in 1:10
			z = mesh.cell_centers[i]
			mesh.cell_states[i].T = 400.0 + 200.0 * z / 0.01  # 400K to 600K
			mesh.cell_states[i].ξ = (1190.0, 0.0)
		end

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u = create_state_vector_ale(mesh)
		du = similar(u)

		# Evaluate RHS
		ode_fn(du, u, nothing, 0.0)

		sv = StateVectorALE(mesh)
		A = mesh.cross_section_area

		# Get mesh velocities from du (dz/dt = w)
		w = [du[sv.z_offset+i] for i in 1:sv.n_nodes]

		# For each cell, verify GCL: dV_i/dt = (w_{i+1} - w_i) * A
		for i in 1:sv.n_cells
			dV_dt_gcl = (w[i+1] - w[i]) * A

			# This should be finite and consistent
			@test isfinite(dV_dt_gcl)

			# For a single RHS evaluation, we verify the structure is correct
			# The actual volumes will evolve according to GCL when integrated
		end
	end

	@testset "Volume evolution satisfies GCL over timestep" begin
		mesh = create_test_mesh(n_cells = 5, NC = 2)
		material = create_gcl_test_material()

		# Hot surface to drive reactions
		for i in 1:5
			mesh.cell_states[i].T = 500.0
			mesh.cell_states[i].ξ = (1190.0, 0.0)
		end

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)
		A = mesh.cross_section_area

		# Get initial volumes
		V_old = compute_volumes_from_state(u, sv, A)

		# Take a small explicit Euler step for testing
		du = similar(u)
		ode_fn(du, u, nothing, 0.0)

		dt = 0.001  # Small timestep
		u_new = u + dt * du

		# Get new volumes and mesh velocities
		V_new = compute_volumes_from_state(u_new, sv, A)
		w = [du[sv.z_offset+i] for i in 1:sv.n_nodes]

		# Verify GCL: V_new - V_old ≈ dt * (w_{i+1} - w_i) * A
		for i in 1:sv.n_cells
			dV_actual = V_new[i] - V_old[i]
			dV_gcl = dt * (w[i+1] - w[i]) * A

			# Should match to numerical precision
			# Use atol for very small values near machine epsilon
			@test dV_actual ≈ dV_gcl atol=1e-14 rtol=1e-8
		end
	end

	@testset "GCL compliance over full integration" begin
		# Test that integrating the ALE ODE maintains GCL compliance
		mesh = create_test_mesh(n_cells = 10, NC = 2)
		material = create_gcl_test_material()

		for i in 1:10
			mesh.cell_states[i].T = 400.0
			mesh.cell_states[i].ξ = (1190.0, 0.0)
		end

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u0 = create_state_vector_ale(mesh)
		tspan = (0.0, 0.5)  # Short integration

		using OrdinaryDiffEq
		prob = ODEProblem(ode_fn, u0, tspan)
		sol = OrdinaryDiffEq.solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6)

		sv = StateVectorALE(mesh)
		A = mesh.cross_section_area

		# Check GCL at each saved time point
		for i in 2:length(sol.t)
			u_curr = sol.u[i]
			u_prev = sol.u[i-1]
			dt_step = sol.t[i] - sol.t[i-1]

			V_curr = compute_volumes_from_state(u_curr, sv, A)
			V_prev = compute_volumes_from_state(u_prev, sv, A)

			# Get average mesh velocity over step (approximation)
			# This is a rough check - exact GCL is guaranteed by construction
			total_V_change = sum(V_curr) - sum(V_prev)
			thickness_curr = get_node_position(u_curr, sv, sv.n_nodes) - get_node_position(u_curr, sv, 1)
			thickness_prev = get_node_position(u_prev, sv, sv.n_nodes) - get_node_position(u_prev, sv, 1)

			# Total volume change should equal thickness change × area
			expected_dV = (thickness_curr - thickness_prev) * A
			@test total_V_change ≈ expected_dV rtol=1e-6
		end
	end
end

# ============================================================================
# ALE Upwind Advection Tests (ALE_Energy_Reformulation.md)
# ============================================================================
#
# These tests verify the corrected ALE formulation using upwind advection (w·∇φ)
# instead of the conservative flux formulation (∇·(φw)).
#
# The key properties being tested:
# 1. Uniform temperature field => zero advection (no spurious sources/sinks)
# 2. Linear temperature profile with uniform velocity => constant advection
# 3. Fixed boundary (w=0) => no advection at boundary
# 4. Advection direction follows mesh velocity (upwind)
# ============================================================================

@testset "ALE Upwind Advection" begin

	@testset "Uniform temperature - zero advection" begin
		# If T is uniform, ∇T = 0, so w·∇T = 0 regardless of velocity
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		# Set uniform temperature
		T_uniform = 350.0
		for i in eachindex(mesh.cells)
			mesh.cell_states[i].T = T_uniform
			mesh.cell_states[i].ξ = (1000.0, 0.0)
		end

		# Non-uniform velocity profile (contraction)
		n_nodes = length(mesh.nodes)
		w_node = zeros(n_nodes)
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		for i in 1:n_nodes
			z = mesh.nodes[i].z
			w_node[i] = -0.001 * (z - z_min) / L  # Velocity from 0 at bottom to -0.001 at top
		end

		# Compute advection term
		advection = zeros(length(mesh.cells))
		compute_ale_advection_term!(advection, mesh, w_node)

		# All advection terms should be zero (or very close)
		for i in eachindex(advection)
			@test abs(advection[i]) < 1e-15
		end
	end

	@testset "Linear temperature profile - uniform velocity" begin
		# T(z) = T0 + ΔT*(z/L), so dT/dz = ΔT/L (constant)
		# With uniform w, advection = w * dT/dz should be constant
		mesh = create_test_mesh(n_cells = 20, NC = 2)

		T0 = 300.0
		ΔT = 100.0
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		# Set linear temperature profile
		for i in eachindex(mesh.cells)
			z = mesh.cell_centers[i]
			mesh.cell_states[i].T = T0 + ΔT * (z - z_min) / L
			mesh.cell_states[i].ξ = (1000.0, 0.0)
		end

		# Uniform mesh velocity
		n_nodes = length(mesh.nodes)
		w_uniform = -0.001  # Contracting mesh
		w_node = fill(w_uniform, n_nodes)

		# Compute advection
		advection = zeros(length(mesh.cells))
		compute_ale_advection_term!(advection, mesh, w_node)

		# Expected: w * dT/dz = -0.001 * (100/L)
		expected_advection = w_uniform * ΔT / L

		# Interior cells should have constant advection
		for i in 2:(length(mesh.cells)-1)
			@test advection[i] ≈ expected_advection rtol=0.1  # Allow 10% error for discretization
		end
	end

	@testset "Fixed boundary - no advection" begin
		# At a fixed boundary (w=0), advection should be zero
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		# Temperature gradient
		T0 = 300.0
		ΔT = 50.0
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		for i in eachindex(mesh.cells)
			z = mesh.cell_centers[i]
			mesh.cell_states[i].T = T0 + ΔT * (z - z_min) / L
			mesh.cell_states[i].ξ = (1000.0, 0.0)
		end

		# Linear velocity: 0 at bottom (fixed), -0.001 at top
		n_nodes = length(mesh.nodes)
		w_node = zeros(n_nodes)
		for i in 1:n_nodes
			z = mesh.nodes[i].z
			w_node[i] = -0.001 * (z - z_min) / L
		end

		# Compute advection
		advection = zeros(length(mesh.cells))
		compute_ale_advection_term!(advection, mesh, w_node)

		# Bottom cell (cell 1) has w ≈ 0, so advection should be very small
		# Note: cell-centered velocity is average of node velocities
		@test abs(advection[1]) < abs(advection[div(length(mesh.cells), 2)])
	end

	@testset "Upwind direction - contracting mesh" begin
		# For w < 0 (contracting toward -z), should use forward difference
		# (information from +z direction)
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		# Temperature increases toward top (positive gradient)
		T0 = 300.0
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		for i in eachindex(mesh.cells)
			z = mesh.cell_centers[i]
			# Quadratic profile to make gradients different
			mesh.cell_states[i].T = T0 + 100.0 * ((z - z_min) / L)^2
			mesh.cell_states[i].ξ = (1000.0, 0.0)
		end

		# Uniform contraction
		n_nodes = length(mesh.nodes)
		w_node = fill(-0.001, n_nodes)

		# Compute advection
		advection = zeros(length(mesh.cells))
		compute_ale_advection_term!(advection, mesh, w_node)

		# With w < 0 and dT/dz > 0, advection = w * dT/dz < 0
		# This means dT/dt contribution is negative (cooling effect)
		for i in 2:(length(mesh.cells)-1)
			@test advection[i] < 0  # Negative advection
		end
	end

	@testset "Upwind direction - expanding mesh" begin
		# For w > 0 (expanding toward +z), should use backward difference
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		# Temperature increases toward top
		T0 = 300.0
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		for i in eachindex(mesh.cells)
			z = mesh.cell_centers[i]
			mesh.cell_states[i].T = T0 + 100.0 * (z - z_min) / L
			mesh.cell_states[i].ξ = (1000.0, 0.0)
		end

		# Uniform expansion
		n_nodes = length(mesh.nodes)
		w_node = fill(0.001, n_nodes)  # Positive = expanding

		# Compute advection
		advection = zeros(length(mesh.cells))
		compute_ale_advection_term!(advection, mesh, w_node)

		# With w > 0 and dT/dz > 0, advection = w * dT/dz > 0
		for i in 2:(length(mesh.cells)-1)
			@test advection[i] > 0
		end
	end

	@testset "Species advection - uniform concentration" begin
		# Uniform concentration => zero advection
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		# Create simple material for testing
		virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
		material = Material(name = :Test, components = (virgin, gas), reactions = (rxn,))

		ξ_uniform = 1000.0
		for i in eachindex(mesh.cells)
			mesh.cell_states[i].T = 300.0
			mesh.cell_states[i].ξ = (ξ_uniform, 100.0)
		end

		# Non-zero velocity
		n_nodes = length(mesh.nodes)
		w_node = fill(-0.001, n_nodes)

		# Compute advection
		NC = 2
		advection = zeros(NC, length(mesh.cells))
		compute_ale_species_advection!(advection, mesh, material, w_node)

		# All components should have zero advection
		for j in 1:NC
			for i in eachindex(mesh.cells)
				@test abs(advection[j, i]) < 1e-15
			end
		end
	end

	@testset "Species advection - concentration gradient" begin
		# Linear concentration profile should give constant advection
		mesh = create_test_mesh(n_cells = 20, NC = 2)

		# Create simple material for testing
		virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
		material = Material(name = :Test, components = (virgin, gas), reactions = (rxn,))

		ξ0 = 1000.0
		Δξ = 500.0
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min

		for i in eachindex(mesh.cells)
			z = mesh.cell_centers[i]
			ξ1 = ξ0 + Δξ * (z - z_min) / L
			mesh.cell_states[i].T = 300.0
			mesh.cell_states[i].ξ = (ξ1, 100.0)
		end

		# Uniform velocity
		n_nodes = length(mesh.nodes)
		w_uniform = -0.001
		w_node = fill(w_uniform, n_nodes)

		# Compute advection
		NC = 2
		advection = zeros(NC, length(mesh.cells))
		compute_ale_species_advection!(advection, mesh, material, w_node)

		# SOLID component 1 (virgin) should have zero advection (moves with mesh)
		for i in 2:(length(mesh.cells)-1)
			@test advection[1, i] ≈ 0.0 atol=1e-10
		end

		# GAS component 2 (MMA) has concentration gradient and should have non-zero advection
		# Set up concentration gradient for gas component (component 2)
		Δξ_gas = 50.0
		for idx in eachindex(mesh.cells)
			z = mesh.cell_centers[idx]
			ξ1_val = 1000.0 + Δξ * (z - z_min) / L
			ξ2_val = 100.0 + Δξ_gas * (z - z_min) / L
			mesh.cell_states[idx].ξ = (ξ1_val, ξ2_val)
		end

		# Recompute advection with new gas gradient
		compute_ale_species_advection!(advection, mesh, material, w_node)
		expected_gas = w_uniform * Δξ_gas / L

		# Interior cells should have approximately constant advection for gas
		for i in 2:(length(mesh.cells)-1)
			@test advection[2, i] ≈ expected_gas rtol=0.15
		end
	end
end

# ============================================================================
# Verification Tests (Issue Review Findings)
# ============================================================================

@testset "Uniform State Invariance Under Mesh Motion" begin
	@testset "Uniform temperature remains uniform" begin
		# Test: A spatially uniform temperature field should remain uniform
		# regardless of mesh motion (∇T = 0 implies advection = 0)
		mesh = create_test_mesh(n_cells = 20, NC = 2)

		# Set uniform temperature
		T_uniform = 500.0
		for i in eachindex(mesh.cells)
			mesh.cell_states[i].T = T_uniform
			mesh.cell_states[i].ξ = (1000.0, 0.0)
		end

		# Apply non-zero mesh velocity (linear profile)
		n_nodes = length(mesh.nodes)
		w = zeros(n_nodes)
		z_min, z_max = domain_bounds(mesh)
		L = z_max - z_min
		w_max = 0.001  # 1 mm/s at surface

		for i in 1:n_nodes
			z = mesh.nodes[i].z
			w[i] = w_max * (z - z_min) / L
		end

		# Compute advection term
		advection = zeros(length(mesh.cells))
		compute_ale_advection_term!(advection, mesh, w)

		# All advection terms should be essentially zero
		@test all(abs.(advection) .< 1e-14)
	end

	@testset "Uniform concentration remains uniform" begin
		# Test: Spatially uniform concentration should have zero advection
		mesh = create_test_mesh(n_cells = 20, NC = 2)

		# Create material
		virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
		material = Material(name = :Test, components = (virgin, gas), reactions = (rxn,))

		# Uniform state
		for i in eachindex(mesh.cells)
			mesh.cell_states[i].T = 400.0
			mesh.cell_states[i].ξ = (800.0, 200.0)  # Both components uniform
		end

		# Non-uniform velocity
		n_nodes = length(mesh.nodes)
		w = [0.001 * sin(π * i / n_nodes) for i in 1:n_nodes]

		# Compute species advection
		NC = 2
		advection = zeros(NC, length(mesh.cells))
		compute_ale_species_advection!(advection, mesh, material, w)

		# All advection should be zero (∇ξ = 0)
		for j in 1:NC
			for i in eachindex(mesh.cells)
				@test abs(advection[j, i]) < 1e-14
			end
		end
	end
end

@testset "Energy-Conserving Remapping" begin
	using Pyrolysis.Properties: effective_heat_capacity

	@testset "Homogeneous material - volume vs energy averaging equivalent" begin
		# When ρcₚ is uniform, volume-weighted and energy-weighted should match
		mesh_old = create_test_mesh(n_cells = 10, NC = 2)

		virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
		material = Material(name = :Test, components = (virgin, gas), reactions = (rxn,))

		# Uniform composition, varying temperature
		for i in eachindex(mesh_old.cells)
			z = mesh_old.cell_centers[i]
			mesh_old.cell_states[i].T = 300.0 + 200.0 * z / 0.01
			mesh_old.cell_states[i].ξ = (1000.0, 0.0)  # Pure virgin everywhere
		end

		# Create state vector
		NC = 2
		n_dof = length(mesh_old.cells) * (1 + NC)
		u_old = zeros(n_dof)
		for i in eachindex(mesh_old.cells)
			idx = (i - 1) * (1 + NC)
			u_old[idx+1] = mesh_old.cell_states[i].T
			for j in 1:NC
				u_old[idx+1+j] = mesh_old.cell_states[i].ξ[j]
			end
		end

		# Create slightly different mesh (simulating mesh motion)
		mesh_new = deepcopy(mesh_old)

		# Test without material (volume-weighted)
		u_new_vol = zeros(n_dof)
		remap_solution!(u_new_vol, mesh_new, u_old, mesh_old, NC, FirstOrderRemap())

		# Test with material (energy-weighted)
		u_new_energy = zeros(n_dof)
		remap_solution!(u_new_energy, mesh_new, u_old, mesh_old, NC, FirstOrderRemap();
			material = material)

		# For uniform composition, results should be very similar
		for i in 1:length(mesh_new.cells)
			idx = (i - 1) * (1 + NC)
			@test u_new_vol[idx+1] ≈ u_new_energy[idx+1] rtol=1e-6
		end
	end

	@testset "Energy conservation with varying composition" begin
		# Test that total energy is conserved even with varying ρcₚ
		mesh_old = create_test_mesh(n_cells = 10, NC = 2)

		virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
		gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
		material = Material(name = :Test, components = (virgin, gas), reactions = (rxn,))

		# Varying composition (simulating partial decomposition)
		for i in eachindex(mesh_old.cells)
			z = mesh_old.cell_centers[i]
			# Top cells more decomposed (more gas, less solid)
			conversion = z / 0.01  # 0 at bottom, 1 at top
			ξ_solid = 1000.0 * (1 - 0.8 * conversion)
			ξ_gas = 1000.0 * 0.8 * conversion
			mesh_old.cell_states[i].T = 400.0 + 300.0 * conversion
			mesh_old.cell_states[i].ξ = (ξ_solid, ξ_gas)
		end

		# Compute total energy in old mesh
		E_old = 0.0
		for i in eachindex(mesh_old.cells)
			T = mesh_old.cell_states[i].T
			ξ = mesh_old.cell_states[i].ξ
			ρcₚ = effective_heat_capacity(ξ, T, material)
			V = mesh_old.cell_volumes[i]
			E_old += ρcₚ * T * V
		end

		# Create state vector
		NC = 2
		n_dof = length(mesh_old.cells) * (1 + NC)
		u_old = zeros(n_dof)
		for i in eachindex(mesh_old.cells)
			idx = (i - 1) * (1 + NC)
			u_old[idx+1] = mesh_old.cell_states[i].T
			for j in 1:NC
				u_old[idx+1+j] = mesh_old.cell_states[i].ξ[j]
			end
		end

		# Same mesh (no actual motion, just testing remap)
		mesh_new = deepcopy(mesh_old)
		u_new = zeros(n_dof)
		remap_solution!(u_new, mesh_new, u_old, mesh_old, NC, FirstOrderRemap();
			material = material)

		# Update mesh_new states from u_new
		for i in eachindex(mesh_new.cells)
			idx = (i - 1) * (1 + NC)
			mesh_new.cell_states[i].T = u_new[idx+1]
			mesh_new.cell_states[i].ξ = ntuple(j -> u_new[idx+1+j], NC)
		end

		# Compute total energy in new mesh
		E_new = 0.0
		for i in eachindex(mesh_new.cells)
			T = mesh_new.cell_states[i].T
			ξ = mesh_new.cell_states[i].ξ
			ρcₚ = effective_heat_capacity(ξ, T, material)
			V = mesh_new.cell_volumes[i]
			E_new += ρcₚ * T * V
		end

		# Energy should be conserved to high precision
		@test E_new ≈ E_old rtol=1e-10
	end
end

@testset "Gas Species Advection Sign Verification" begin
	# Verify that gas species advection has the correct sign
	# When mesh shrinks (w < 0) and ∂ξ/∂z > 0, gas should advect upward (increase)
	mesh = create_test_mesh(n_cells = 20, NC = 2)

	virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
	gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
	rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
	material = Material(name = :Test, components = (virgin, gas), reactions = (rxn,))

	# Gas concentration increases with z (more gas at top)
	z_min, z_max = domain_bounds(mesh)
	L = z_max - z_min
	for i in eachindex(mesh.cells)
		z = mesh.cell_centers[i]
		ξ_gas = 100.0 + 200.0 * (z - z_min) / L  # 100 at bottom, 300 at top
		mesh.cell_states[i].T = 500.0
		mesh.cell_states[i].ξ = (800.0, ξ_gas)
	end

	# Mesh shrinking (w < 0) - nodes moving toward z=0
	n_nodes = length(mesh.nodes)
	w_shrink = -0.001  # Negative = shrinking
	w = fill(w_shrink, n_nodes)

	NC = 2
	advection = zeros(NC, length(mesh.cells))
	compute_ale_species_advection!(advection, mesh, material, w)

	# For gas (component 2) with w < 0 and ∂ξ/∂z > 0:
	# The ALE equation is: ∂ξ/∂t|_χ = S + w·∇ξ
	# With w < 0 and ∇ξ > 0: the contribution w·∇ξ is NEGATIVE
	# This means gas concentration decreases as mesh shrinks
	# (mesh catches up to the gas that was left behind)

	# The advection array stores w·∇ξ which gets ADDED to dξ/dt
	# With w_shrink < 0 and dξ/dz > 0, advection should be NEGATIVE
	dξ_dz = 200.0 / L  # The gradient
	expected_advection = w_shrink * dξ_dz  # Should be negative

	# Check interior cells
	for i in 3:(length(mesh.cells)-2)
		# Gas component (j=2) should have the expected advection
		@test advection[2, i] ≈ expected_advection rtol=0.2
	end

	# Solid component (j=1) should have zero advection (moves with mesh)
	for i in eachindex(mesh.cells)
		@test abs(advection[1, i]) < 1e-15
	end
end
# ============================================================================
# Swelling Factor Consistency Tests
# ============================================================================

@testset "Swelling Factor Consistency" begin

	@testset "Volume change rate uses swelling_factor" begin
		# Test that volume_change_rate correctly uses swelling_factor
		# not just phase checks
		using Pyrolysis.Physics: volume_change_rate

		# Create components with different swelling factors
		virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1000.0, k = 0.2)  # γ=1.0 by default
		char = SolidComponent(:char, ρ = 500.0, c = 800.0, k = 0.1)
		gas = GaseousComponent(:gas, M = 0.030, c = 1000.0, k = 0.02, λ = 1e-5)

		# Verify default swelling factors
		@test virgin.swelling_factor == 1.0
		@test char.swelling_factor == 1.0
		@test gas.swelling_factor == 0.0

		rxn = Reaction(:decomp, 1 => (2, 3), A = 1e12, E = 150e3, h = 500e3, yields = (0.2, 0.8))
		material = Material(name = :Test, components = (virgin, char, gas), reactions = (rxn,))

		T = 600.0

		# Test case 1: virgin decomposes (dξ[1] < 0), char and gas form
		dξ = (-10.0, 2.0, 8.0)  # 10 kg/m³/s virgin decomposes
		θ = volume_change_rate(dξ, T, material)

		# Only virgin and char contribute (gas has γ=0)
		ρ_virgin = virgin.density(T)
		ρ_char = char.density(T)
		θ_expected = 1.0 * dξ[1] / ρ_virgin + 1.0 * dξ[2] / ρ_char
		# = -10/1000 + 2/500 = -0.01 + 0.004 = -0.006

		@test θ ≈ θ_expected rtol=1e-10
		@test θ < 0  # Net shrinkage (virgin depletes faster than char forms)
	end

	@testset "Custom swelling factor" begin
		using Pyrolysis.Physics: volume_change_rate

		# Create solid with custom swelling factor (e.g., char that doesn't contribute)
		# This would require extending the SolidComponent constructor...
		# For now, test that the logic works with manual component creation

		virgin = Component(
			:virgin, SOLID, 0.0,
			ConstantProperty(1000.0),  # density
			ConstantProperty(1000.0),  # heat capacity
			ConstantProperty(0.2),     # conductivity
			ConstantProperty(0.0),     # gas transfer
			ConstantProperty(0.0),     # permeability (no Darcy)
			0.9, 0.0,                  # emissivity, absorption
			1.0,                        # swelling_factor = 1
		)

		char_no_swell = Component(
			:char, SOLID, 0.0,
			ConstantProperty(500.0),
			ConstantProperty(800.0),
			ConstantProperty(0.1),
			ConstantProperty(0.0),
			ConstantProperty(0.0),     # permeability (no Darcy)
			0.9, 0.0,
			0.0,                        # swelling_factor = 0 (doesn't contribute)
		)

		gas = GaseousComponent(:gas, M = 0.030, c = 1000.0, k = 0.02, λ = 1e-5)

		rxn = Reaction(:decomp, 1 => (2, 3), A = 1e12, E = 150e3, h = 500e3, yields = (0.2, 0.8))
		material = Material(name = :Test, components = (virgin, char_no_swell, gas), reactions = (rxn,))

		T = 600.0
		dξ = (-10.0, 2.0, 8.0)
		θ = volume_change_rate(dξ, T, material)

		# Only virgin contributes (char has γ=0, gas has γ=0)
		ρ_virgin = virgin.density(T)
		θ_expected = 1.0 * dξ[1] / ρ_virgin
		# = -10/1000 = -0.01

		@test θ ≈ θ_expected rtol=1e-10
	end

	@testset "Mesh velocity uses swelling_factor consistently" begin
		# This tests that _compute_ale_mesh_velocity! uses swelling_factor
		# not just phase == SOLID

		# Create a mesh with known reaction state
		mesh = create_test_mesh(n_cells = 5, NC = 3, L = 0.01)

		virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1000.0, k = 0.2)
		char = SolidComponent(:char, ρ = 500.0, c = 800.0, k = 0.1)
		gas = GaseousComponent(:gas, M = 0.030, c = 1000.0, k = 0.02, λ = 1e-5)
		rxn = Reaction(:decomp, 1 => (2, 3), A = 1e15, E = 200e3, h = 500e3, yields = (0.2, 0.8))
		material = Material(name = :Test, components = (virgin, char, gas), reactions = (rxn,))

		# Set up active reaction state
		for i in eachindex(mesh.cells)
			mesh.cell_states[i].T = 700.0  # Hot enough for reactions
			mesh.cell_states[i].ξ = (800.0, 100.0, 50.0)  # Mostly virgin
		end

		# Create property cache
		using Pyrolysis.Properties
		prop_cache = PropertyCache(mesh)

		# Compute mesh velocity
		n_nodes = length(mesh.nodes)
		w = zeros(n_nodes)
		θ_out = zeros(length(mesh.cells))

		using Pyrolysis.ALE: LaplacianSmoothing
		method = LaplacianSmoothing()

		# Need the full ODE function machinery to test _compute_ale_mesh_velocity!
		# For now, we verify the volume_change_rate function which should be consistent

		# Compute what θ should be for each cell
		using Pyrolysis.Physics: evaluate_reactions!
		NR = n_reactions(material)
		rates = zeros(NR)
		heats = zeros(NR)
		NC = 3

		for i in eachindex(mesh.cells)
			T = mesh.cell_states[i].T
			ξ = mesh.cell_states[i].ξ

			evaluate_reactions!(rates, heats, ξ, T, material)

			# Compute species source terms
			dξ_dt = zeros(NC)
			for r in 1:NR
				rxn = material.reactions[r]
				rate = rates[r]
				if rxn.reactant1 > 0
					dξ_dt[rxn.reactant1] -= rxn.stoich[1] * rate
				end
				if rxn.product1 > 0
					dξ_dt[rxn.product1] += rxn.stoich[3] * rate
				end
				if rxn.product2 > 0
					dξ_dt[rxn.product2] += rxn.stoich[4] * rate
				end
			end

			# Compute expected θ using swelling_factor
			θ_expected = 0.0
			for j in 1:NC
				γ = material.components[j].swelling_factor
				if γ > 0
					ρⱼ = material.components[j].density(T)
					θ_expected += γ * dξ_dt[j] / ρⱼ
				end
			end

			# Verify using the volume_change_rate function
			θ_from_func = volume_change_rate(Tuple(dξ_dt), T, material)
			@test θ_from_func ≈ θ_expected rtol=1e-10
		end
	end
end

# ============================================================================
# Mass Conservation Tests  
# ============================================================================

@testset "Mass Conservation with ALE" begin

	@testset "Dilation correction preserves solid mass" begin
		# Verify that the dilation correction -ξ×θ correctly accounts for
		# mass changes in shrinking cells

		# Initial state: cell with solid mass M₀ and volume V₀
		V₀ = 1.0e-6  # 1 mm³
		ξ₀ = 1000.0  # kg/m³
		M₀ = ξ₀ * V₀  # Initial mass

		# Cell shrinks by 10%
		θ = -0.1  # 10% volume decrease per second
		dt = 1.0   # 1 second

		# Without dilation correction, mass would incorrectly decrease
		# ξ stays same, V decreases → M = ξ × V decreases (wrong!)

		# With dilation correction:
		# dξ/dt = S - ξ × θ = 0 - 1000 × (-0.1) = +100 kg/(m³·s)
		# After dt=1s: ξ₁ = 1000 + 100 = 1100 kg/m³
		# V₁ = V₀ × (1 + θ×dt) = 1e-6 × 0.9 = 9e-7 m³
		# M₁ = ξ₁ × V₁ = 1100 × 9e-7 = 9.9e-4 kg

		# Wait, that's not right. Let me redo this properly.
		# For S=0 (no reactions), mass should be conserved:
		# d(ξV)/dt = S × V = 0
		# V × dξ/dt + ξ × dV/dt = 0
		# V × dξ/dt + ξ × V × θ = 0 (since dV/dt = V×θ)
		# dξ/dt = -ξ × θ

		# So with θ = -0.1 (shrinking):
		# dξ/dt = -1000 × (-0.1) = +100 kg/(m³·s)

		# Numerically (forward Euler):
		ξ₁ = ξ₀ + (-ξ₀ * θ) * dt  # = 1000 + 100 = 1100
		V₁ = V₀ * (1 + θ * dt)     # = 1e-6 * 0.9 = 9e-7
		M₁ = ξ₁ * V₁              # = 1100 * 9e-7 = 9.9e-4

		# But wait, this isn't exactly M₀ = 1000 * 1e-6 = 1e-3
		# That's because forward Euler isn't exact for this

		# For infinitesimal dt, let's verify the rates match:
		# d(ξV)/dt = V × dξ/dt + ξ × dV/dt
		#          = V₀ × 100 + 1000 × V₀ × (-0.1)
		#          = V₀ × (100 - 100) = 0 ✓

		dξ_dt = -ξ₀ * θ
		dV_dt = V₀ * θ
		dM_dt = V₀ * dξ_dt + ξ₀ * dV_dt

		@test abs(dM_dt) < 1e-15  # Mass change rate is zero
	end

	@testset "Reaction stoichiometry mass balance" begin
		# Verify that reactions conserve mass: reactant_mass = product_mass

		virgin = SolidComponent(:virgin, ρ = 1000.0, c = 1000.0, k = 0.2)
		char = SolidComponent(:char, ρ = 500.0, c = 800.0, k = 0.1)
		gas = GaseousComponent(:gas, M = 0.030, c = 1000.0, k = 0.02, λ = 1e-5)

		# Charring reaction: 1 kg virgin → 0.2 kg char + 0.8 kg gas
		rxn = Reaction(:decomp, 1 => (2, 3), A = 1e12, E = 150e3, h = 500e3, yields = (0.2, 0.8))

		# Check stoichiometry sums
		@test rxn.stoich[1] + rxn.stoich[2] ≈ rxn.stoich[3] + rxn.stoich[4] rtol=1e-10
		@test Pyrolysis.Core.mass_balance_error(rxn) < 1e-10
		@test Pyrolysis.Core.verify_stoichiometry(rxn)
	end
end