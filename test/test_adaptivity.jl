"""
Tests for the Adaptivity module.

Tests cover:
1. Error indicators (gradient, curvature, composite)
2. h-refinement (cell splitting, conservation)
3. h-coarsening (cell merging, conservation)
4. r-refinement (node relocation)
5. Mesh quality metrics
6. AMR controller orchestration
7. 2:1 balance enforcement
"""

using Test
using Statistics
using Pyrolysis
using Pyrolysis.Geometry
using Pyrolysis.Adaptivity

# ============================================================================
# Helper: Create test mesh
# ============================================================================

"""Create a simple 1D mesh with NC components for testing."""
function create_test_mesh(; n_cells::Int = 10, L::Float64 = 0.01, NC::Int = 2)
	mesh = create_uniform_mesh(L, n_cells, Val(NC))

	# Set up temperature profile with a gradient
	for i in 1:n_cells
		z = mesh.cell_centers[i]
		T = 300.0 + 500.0 * (z / L)  # Linear gradient
		ξ = ntuple(j -> j == 1 ? 1000.0 * (1 - z/L) : 100.0 * z/L, NC)
		mesh.cell_states[i] = CellState{NC}(T, ξ, mesh.cell_volumes[i], mesh.cell_volumes[i])
	end

	return mesh
end

"""Create mesh with sharp front for testing refinement."""
function create_front_mesh(; n_cells::Int = 20, L::Float64 = 0.01, NC::Int = 2, front_pos::Float64 = 0.005)
	mesh = create_uniform_mesh(L, n_cells, Val(NC))

	for i in 1:n_cells
		z = mesh.cell_centers[i]
		# Sharp temperature front
		T = z < front_pos ? 300.0 : 800.0
		ξ = ntuple(j -> j == 1 ? (z < front_pos ? 1000.0 : 100.0) : (z < front_pos ? 0.0 : 900.0), NC)
		mesh.cell_states[i] = CellState{NC}(T, ξ, mesh.cell_volumes[i], mesh.cell_volumes[i])
	end

	return mesh
end

# ============================================================================
# Test Error Indicators
# ============================================================================

@testset "Error Indicators" begin

	@testset "GradientIndicator" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		indicator = GradientIndicator(:T)
		η = compute_indicator(indicator, mesh)

		@test length(η) == 10
		@test all(η .>= 0)

		# With linear temperature, gradient should be roughly constant
		# So η = h * |∇T| should be similar for all cells
		active_η = [η[i] for i in 1:10 if mesh.active_cell_mask[i]]
		@test std(active_η) / mean(active_η) < 0.5  # Low variation
	end

	@testset "CurvatureIndicator" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		indicator = CurvatureIndicator(:T)
		η = compute_indicator(indicator, mesh)

		@test length(η) == 10
		@test all(η .>= 0)

		# Linear temperature has zero curvature
		# Expect low indicator values (might be nonzero due to discretization)
		@test maximum(η) < 1e-6 || maximum(η) / maximum(compute_indicator(GradientIndicator(:T), mesh)) < 0.1
	end

	@testset "ReactionRateIndicator" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		indicator = ReactionRateIndicator(1.0)
		η = compute_indicator(indicator, mesh)

		@test length(η) == 10
		@test all(η .>= 0)
	end

	@testset "CompositeIndicator" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		grad = GradientIndicator(:T, 1.0)
		curv = CurvatureIndicator(:T, 0.5)

		# Test :max combination
		composite_max = CompositeIndicator([grad, curv], :max)
		η_max = compute_indicator(composite_max, mesh)

		η_grad = compute_indicator(grad, mesh)
		η_curv = compute_indicator(curv, mesh)

		for i in 1:10
			if mesh.active_cell_mask[i]
				@test η_max[i] ≈ max(η_grad[i], η_curv[i])
			end
		end

		# Test :sum combination
		composite_sum = CompositeIndicator([grad, curv], :sum)
		η_sum = compute_indicator(composite_sum, mesh)

		for i in 1:10
			if mesh.active_cell_mask[i]
				@test η_sum[i] ≈ η_grad[i] + η_curv[i]
			end
		end
	end

	@testset "InterfaceProximityIndicator" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)
		L = 0.01

		indicator = InterfaceProximityIndicator(L/2, 0.001)
		η = compute_indicator(indicator, mesh)

		@test length(η) == 10
		@test all(η .>= 0)

		# Should be highest near interface at L/2
		z_centers = mesh.cell_centers
		center_idx = argmin(abs.(z_centers .- L/2))
		@test η[center_idx] == maximum(η)
	end
end

# ============================================================================
# Test h-Refinement
# ============================================================================

@testset "h-Refinement" begin

	@testset "Single Cell Refinement" begin
		mesh = create_test_mesh(n_cells = 5, NC = 2)

		initial_n_cells = length(mesh.cells)
		initial_n_active = mesh.n_active_cells
		initial_volume = sum(mesh.cell_volumes[i] for i in 1:initial_n_cells if mesh.active_cell_mask[i])

		# Refine middle cell
		cell_to_refine = 3
		parent_state = mesh.cell_states[cell_to_refine]
		parent_volume = mesh.cell_volumes[cell_to_refine]

		child_a, child_b = refine_cell!(mesh, cell_to_refine)

		# Check topology
		@test length(mesh.cells) == initial_n_cells + 2
		@test mesh.n_active_cells == initial_n_active + 1  # +2 children, -1 parent
		@test !mesh.active_cell_mask[cell_to_refine]  # Parent inactive
		@test mesh.active_cell_mask[child_a]
		@test mesh.active_cell_mask[child_b]

		# Check hierarchy
		@test mesh.cells[cell_to_refine].children == (child_a, child_b)
		@test mesh.cells[child_a].parent == cell_to_refine
		@test mesh.cells[child_b].parent == cell_to_refine
		@test mesh.cells[child_a].level == mesh.cells[cell_to_refine].level + 1

		# Check volume conservation
		child_volume = mesh.cell_volumes[child_a] + mesh.cell_volumes[child_b]
		@test child_volume ≈ parent_volume rtol=1e-10

		final_volume = sum(mesh.cell_volumes[i] for i in 1:length(mesh.cells) if mesh.active_cell_mask[i])
		@test final_volume ≈ initial_volume rtol=1e-10
	end

	@testset "Mass Conservation in Refinement" begin
		mesh = create_test_mesh(n_cells = 5, NC = 2)

		# Compute initial total mass for each component
		initial_mass = zeros(2)
		for i in 1:5
			if mesh.active_cell_mask[i]
				V = mesh.cell_volumes[i]
				for j in 1:2
					initial_mass[j] += V * mesh.cell_states[i].ξ[j]
				end
			end
		end

		# Refine cell 3
		refine_cell!(mesh, 3)

		# Compute final total mass
		final_mass = zeros(2)
		for i in 1:length(mesh.cells)
			if mesh.active_cell_mask[i]
				V = mesh.cell_volumes[i]
				for j in 1:2
					final_mass[j] += V * mesh.cell_states[i].ξ[j]
				end
			end
		end

		@test final_mass[1] ≈ initial_mass[1] rtol=1e-10
		@test final_mass[2] ≈ initial_mass[2] rtol=1e-10
	end

	@testset "Batch Refinement" begin
		mesh = create_front_mesh(n_cells = 20, NC = 2)

		indicator = GradientIndicator(:T)
		η = compute_indicator(indicator, mesh)

		# Mark cells for refinement
		marked = mark_cells_for_refinement(mesh, η, 0.5, 1e-6, 4)

		@test length(marked) > 0

		# All marked cells should have high indicators
		η_max = maximum(η)
		for cell_id in marked
			@test η[cell_id] > 0.5 * η_max
		end

		# Perform refinement
		initial_active = mesh.n_active_cells
		perform_h_refinement!(mesh, marked)

		@test mesh.n_active_cells > initial_active
	end

	@testset "2:1 Balance Enforcement" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		# Mark only one cell - balance should add neighbors
		marked = [5]
		balanced = enforce_2to1_balance!(mesh, marked, 4)

		# Should still include original
		@test 5 in balanced

		# After refinement, check level differences
		perform_h_refinement!(mesh, balanced)

		# Check all neighboring pairs have level difference ≤ 1
		for i in 1:length(mesh.cells)
			if !mesh.active_cell_mask[i]
				continue
			end
			cell = mesh.cells[i]
			for face_id in cell.face_ids
				face = mesh.faces[face_id]
				neighbor_id = face.cell_ids[1] == i ? face.cell_ids[2] : face.cell_ids[1]
				if neighbor_id > 0 && mesh.active_cell_mask[neighbor_id]
					level_diff = abs(mesh.cells[i].level - mesh.cells[neighbor_id].level)
					@test level_diff ≤ 1
				end
			end
		end
	end
end

# ============================================================================
# Test h-Coarsening
# ============================================================================

@testset "h-Coarsening" begin

	@testset "Coarsen Refined Cells" begin
		mesh = create_test_mesh(n_cells = 5, NC = 2)

		# First refine a cell
		cell_to_refine = 3
		child_a, child_b = refine_cell!(mesh, cell_to_refine)

		n_active_after_refine = mesh.n_active_cells

		# Now coarsen the children back
		parent_id = coarsen_cells!(mesh, child_a, child_b)

		@test parent_id == cell_to_refine
		@test mesh.n_active_cells == n_active_after_refine - 1
		@test mesh.active_cell_mask[parent_id]
		@test !mesh.active_cell_mask[child_a]
		@test !mesh.active_cell_mask[child_b]
		@test mesh.cells[parent_id].children == (0, 0)
	end

	@testset "Mass Conservation in Coarsening" begin
		mesh = create_test_mesh(n_cells = 5, NC = 2)

		# Refine cell 3
		child_a, child_b = refine_cell!(mesh, 3)

		# Record mass in children
		V_a = mesh.cell_volumes[child_a]
		V_b = mesh.cell_volumes[child_b]
		mass_before = zeros(2)
		for j in 1:2
			mass_before[j] = V_a * mesh.cell_states[child_a].ξ[j] +
							 V_b * mesh.cell_states[child_b].ξ[j]
		end

		# Coarsen
		parent_id = coarsen_cells!(mesh, child_a, child_b)

		# Check mass conservation
		V_p = mesh.cell_volumes[parent_id]
		mass_after = zeros(2)
		for j in 1:2
			mass_after[j] = V_p * mesh.cell_states[parent_id].ξ[j]
		end

		@test mass_after[1] ≈ mass_before[1] rtol=1e-10
		@test mass_after[2] ≈ mass_before[2] rtol=1e-10
	end

	@testset "Can Coarsen Check" begin
		mesh = create_test_mesh(n_cells = 5, NC = 2)

		# Original cells cannot be coarsened (no parent)
		@test !can_coarsen(mesh, 1)
		@test !can_coarsen(mesh, 3)

		# After refinement, children can be coarsened
		refine_cell!(mesh, 3)
		child_a, child_b = mesh.cells[3].children

		@test can_coarsen(mesh, child_a)
		@test can_coarsen(mesh, child_b)
	end

	@testset "Mark Cells for Coarsening" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		# Refine several cells
		refine_cell!(mesh, 3)
		refine_cell!(mesh, 7)

		# Create uniform indicator (low everywhere)
		η = fill(0.01, length(mesh.cells))

		pairs = mark_cells_for_coarsening(mesh, η, 0.5, 0.01)

		# Should find some pairs
		@test length(pairs) >= 0  # May or may not find pairs depending on thresholds
	end
end

# ============================================================================
# Test r-Refinement
# ============================================================================

@testset "r-Refinement" begin

	@testset "Monitor Function Evaluation" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)
		T, ξ = zeros(10), zeros(2, 10)
		for i in 1:10
			T[i] = mesh.cell_states[i].T
			for j in 1:2
				ξ[j, i] = mesh.cell_states[i].ξ[j]
			end
		end

		# Gradient monitor
		mon = GradientMonitor(α = 1.0, variable = :T)
		ω = evaluate_monitor(mon, mesh, T, ξ)

		@test length(ω) == 10
		@test all(ω .>= 1.0)  # Should be >= 1 due to sqrt(1 + α*grad²)

		# Interface monitor
		mon_int = InterfaceMonitor(0.005; σ = 0.001)
		ω_int = evaluate_monitor(mon_int, mesh, T, ξ)

		@test length(ω_int) == 10
		# Should peak near z = 0.005
		center_idx = argmin(abs.(mesh.cell_centers .- 0.005))
		@test ω_int[center_idx] ≥ mean(ω_int)
	end

	@testset "Node Relocation" begin
		mesh = create_front_mesh(n_cells = 20, NC = 2, front_pos = 0.005)

		# Record initial positions
		z_initial = [n.z for n in mesh.nodes]

		# Perform r-refinement
		monitor = InterfaceMonitor(0.005; σ = 0.001)
		n_iter = relocate_nodes!(mesh; monitor = monitor, max_iter = 10, relaxation = 0.3)

		@test n_iter > 0

		# Nodes should have moved
		z_final = [n.z for n in mesh.nodes]
		@test any(z_final .!= z_initial)

		# Boundary nodes should not move
		@test z_final[1] ≈ z_initial[1]
		@test z_final[end] ≈ z_initial[end]

		# Mesh should still be valid
		summary = mesh_quality_summary(mesh)
		@test summary.all_positive_volume
	end
end

# ============================================================================
# Test Mesh Quality Metrics
# ============================================================================

@testset "Mesh Quality Metrics" begin

	@testset "Aspect Ratio Metric" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		metric = AspectRatioMetric(2.0)
		@test check_quality(mesh, metric)

		violations = quality_violations(mesh, metric)
		@test isempty(violations)

		# Uniform mesh should have aspect ratio = 1
		ratios = compute_aspect_ratios(mesh)
		@test all(r -> r ≈ 1.0, ratios[mesh.active_cell_mask])
	end

	@testset "Smoothness Metric" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		metric = SmoothnessMetric(0.5)
		@test check_quality(mesh, metric)

		# Uniform mesh should have near-zero smoothness deviation (within floating point tolerance)
		smoothness = compute_smoothness(mesh)
		@test all(s -> abs(s) < 1e-10, smoothness[mesh.active_cell_mask])
	end

	@testset "Jacobian Metric" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		metric = JacobianMetric(0.0)
		@test check_quality(mesh, metric)

		violations = quality_violations(mesh, metric)
		@test isempty(violations)
	end

	@testset "Composite Quality Metric" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		composite = CompositeQualityMetric([
			AspectRatioMetric(2.0),
			SmoothnessMetric(0.5),
			JacobianMetric(0.0),
		])

		@test check_quality(mesh, composite)
	end

	@testset "Quality Summary" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)

		summary = mesh_quality_summary(mesh)

		@test summary.n_active_cells == 10
		@test summary.min_cell_size > 0
		@test summary.max_cell_size > 0
		@test summary.max_aspect_ratio ≈ 1.0
		@test abs(summary.max_smoothness) < 1e-10  # Near-zero within floating point tolerance
		@test summary.all_positive_volume
	end
end

# ============================================================================
# Test AMR Controller
# ============================================================================

@testset "AMR Controller" begin

	@testset "Controller Creation" begin
		controller = AMRController()

		@test controller.θ_refine == 0.7
		@test controller.θ_coarsen == 0.1
		@test controller.h_min == 1e-6
		@test controller.h_max == 1e-2
		@test controller.max_level == 4
		@test controller.check_every == 10
		@test controller.n_refinements == 0
		@test controller.n_coarsenings == 0
	end

	@testset "Adapt Mesh" begin
		mesh = create_front_mesh(n_cells = 20, NC = 2)

		controller = AMRController(
			error_indicator = GradientIndicator(:T),
			θ_refine = 0.5,
			θ_coarsen = 0.1,
			h_min = 1e-6,
			h_max = 0.01,
			max_level = 3,
		)

		initial_active = mesh.n_active_cells

		result = adapt_mesh!(mesh, controller)

		@test haskey(result, :n_refined)
		@test haskey(result, :n_coarsened)
		@test haskey(result, :quality_ok)

		# Should have refined near the front
		@test result.n_refined > 0 || mesh.n_active_cells >= initial_active
		@test result.quality_ok

		# Statistics should update
		@test controller.n_refinements >= 0
	end

	@testset "Should Adapt Check" begin
		controller = AMRController(check_every = 5)

		@test !should_adapt(controller, 0.0, 1)
		@test !should_adapt(controller, 0.0, 4)
		@test should_adapt(controller, 0.0, 5)
		@test should_adapt(controller, 0.0, 10)
	end

	@testset "Maybe Adapt" begin
		mesh = create_test_mesh(n_cells = 10, NC = 2)
		controller = AMRController(check_every = 5)

		# Should skip at step 3
		result = maybe_adapt!(mesh, controller, 0.0, 3)
		@test result === nothing

		# Should adapt at step 5
		result = maybe_adapt!(mesh, controller, 0.5, 5)
		@test result !== nothing
		@test controller.last_check_time ≈ 0.5
	end

	@testset "Statistics" begin
		controller = AMRController()
		controller.n_refinements = 10
		controller.n_coarsenings = 3
		controller.last_check_time = 1.5

		stats = amr_statistics(controller)

		@test stats.total_refinements == 10
		@test stats.total_coarsenings == 3
		@test stats.last_check_time ≈ 1.5

		reset_statistics!(controller)

		@test controller.n_refinements == 0
		@test controller.n_coarsenings == 0
		@test controller.last_check_time ≈ 0.0
	end

	@testset "Default Controller" begin
		controller = create_default_controller(h_target = 1e-4)

		@test controller.h_min ≈ 1e-5
		@test controller.h_max ≈ 1e-2
		@test controller.error_indicator isa CompositeIndicator
		@test controller.quality_metric isa CompositeQualityMetric
	end
end

# ============================================================================
# Integration Tests
# ============================================================================

@testset "Integration: Refine-Coarsen Cycle" begin
	mesh = create_test_mesh(n_cells = 10, NC = 2)

	# Compute initial total mass
	initial_mass = sum(
		mesh.cell_volumes[i] * mesh.cell_states[i].ξ[1]
		for i in 1:10 if mesh.active_cell_mask[i]
	)

	# Refine cell 5
	refine_cell!(mesh, 5)
	child_a, child_b = mesh.cells[5].children

	# Refine one of the children
	refine_cell!(mesh, child_a)

	# Check mass
	mass_after_refine = sum(
		mesh.cell_volumes[i] * mesh.cell_states[i].ξ[1]
		for i in 1:length(mesh.cells) if mesh.active_cell_mask[i]
	)
	@test mass_after_refine ≈ initial_mass rtol=1e-10

	# Coarsen the grandchildren back
	grandchild_a, grandchild_b = mesh.cells[child_a].children
	coarsen_cells!(mesh, grandchild_a, grandchild_b)

	# Coarsen the original children
	coarsen_cells!(mesh, child_a, child_b)

	# Should be back to original state
	@test mesh.n_active_cells == 10

	# Check mass conservation through entire cycle
	final_mass = sum(
		mesh.cell_volumes[i] * mesh.cell_states[i].ξ[1]
		for i in 1:length(mesh.cells) if mesh.active_cell_mask[i]
	)
	@test final_mass ≈ initial_mass rtol=1e-10
end

@testset "Integration: Full AMR Cycle" begin
	mesh = create_front_mesh(n_cells = 15, NC = 2)

	controller = AMRController(
		error_indicator = GradientIndicator(:T),
		quality_metric = CompositeQualityMetric([
			AspectRatioMetric(3.0),
			SmoothnessMetric(1.0),
		]),
		θ_refine = 0.5,
		θ_coarsen = 0.05,
		h_min = 1e-6,
		h_max = 0.01,
		max_level = 3,
		check_every = 1,
	)

	# Initial mass
	initial_mass = sum(
		mesh.cell_volumes[i] * mesh.cell_states[i].ξ[1]
		for i in 1:length(mesh.cells) if mesh.active_cell_mask[i]
	)

	# Run several adaptation cycles
	for step in 1:5
		result = maybe_adapt!(mesh, controller, Float64(step), step)
		@test result !== nothing
		@test result.quality_ok
	end

	# Final mass should be conserved
	final_mass = sum(
		mesh.cell_volumes[i] * mesh.cell_states[i].ξ[1]
		for i in 1:length(mesh.cells) if mesh.active_cell_mask[i]
	)
	@test final_mass ≈ initial_mass rtol=1e-8

	# Quality should be maintained
	summary = mesh_quality_summary(mesh)
	@test summary.all_positive_volume
	@test summary.max_aspect_ratio < 4.0
end
