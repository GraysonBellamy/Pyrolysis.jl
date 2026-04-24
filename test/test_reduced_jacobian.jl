# ============================================================================
# Tests for Approximate Jacobian (reduced sparsity ALE solver)
# ============================================================================

using Test
using Pyrolysis
using Pyrolysis.Core
using Pyrolysis.Geometry
using Pyrolysis.TimeIntegration
using Pyrolysis.TimeIntegration: PyrolysisODEFunctionALE, create_state_vector_ale,
	JacobianCacheALE, compute_analytical_jacobian_ale_reduced!,
	jacobian_sparsity_ale_reduced, StateVectorALE
using SparseArrays

# Helper: create a simple ALE test problem (PMMA-like, 2 components)
function create_approx_jac_test_problem(; n_cells = 10, T_initial = 400.0)
	material = Material(
		name = :test_pmma,
		components = (
			SolidComponent(:virgin,
				ρ = ConstantProperty(1190.0),
				c = LinearProperty(1420.0, 2.5),
				k = ConstantProperty(0.21),
				ε = 0.86,
				α = 1.5,
				γ = 1.0),
			GaseousComponent(:volatiles,
				M = 0.100,
				c = ConstantProperty(1500.0),
				k = ConstantProperty(0.02),
				λ = ConstantProperty(2e-5),
				ε = 0.0,
				α = 0.0),
		),
		reactions = (
			Reaction(:decomposition, 1 => 2,
				A = 8.5e12,
				E = 188e3,
				h = 870e3,
				T_min = 300.0,
				T_max = 1500.0),
		),
	)

	mesh = create_uniform_mesh(0.005, n_cells, Val(2);
		T_initial = T_initial,
		ξ_initial = (1190.0, 0.0))

	bc_set = BoundaryConditionSet(
		thermal_top = HeatFluxBC(50e3),
		thermal_bottom = AdiabaticBC(),
		mass_top = DiffusiveMassBC(0.01),
		mass_bottom = ImpermeableBC(),
	)

	return material, mesh, bc_set
end

@testset "Approximate Jacobian (Reduced Sparsity ALE)" begin

	# ======================================================================
	# Test 1: Reduced sparsity pattern is sparser than full ALE
	# ======================================================================
	@testset "Reduced sparsity is sparser than full ALE" begin
		for n_cells in [10, 30, 60]
			NC = 2
			jac_reduced = jacobian_sparsity_ale_reduced(n_cells, NC)
			jac_ale = Pyrolysis.TimeIntegration.jacobian_sparsity_ale(n_cells, NC)

			@test nnz(jac_reduced) < nnz(jac_ale)
			@test size(jac_reduced) == size(jac_ale)

			if n_cells >= 30
				savings = 1.0 - nnz(jac_reduced) / nnz(jac_ale)
				@test savings > 0.3
			end
		end
	end

	# ======================================================================
	# Test 2: Reduced Jacobian matches FD on physics block
	# ======================================================================
	@testset "Reduced Jacobian vs finite differences (physics block)" begin
		material, mesh, bc_set = create_approx_jac_test_problem(n_cells = 8, T_initial = 400.0)
		NC = 2
		NR = 1
		n_cells = length(mesh.cells)
		n_faces = length(mesh.faces)
		n_nodes = length(mesh.nodes)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)
		u0 = create_state_vector_ale(mesh, material)
		n_dof = length(u0)

		# Compute reduced analytical Jacobian
		jac_cache = JacobianCacheALE(n_cells, n_faces, n_nodes, NC, NR)
		J_analytical = zeros(n_dof, n_dof)
		compute_analytical_jacobian_ale_reduced!(J_analytical, u0, nothing, 0.0, ode_fn, jac_cache)

		# Compute finite difference Jacobian of the FULL monolithic RHS
		# (the reduced Jacobian should match on physics-physics entries,
		#  but intentionally omit ALE coupling entries)
		J_fd = zeros(n_dof, n_dof)
		du_base = similar(u0)
		du_pert = similar(u0)
		u_pert = copy(u0)

		ode_fn(du_base, u0, nothing, 0.0)

		for j in 1:n_dof
			u_pert .= u0
			h = max(abs(u0[j]) * 1e-7, 1e-7)
			u_pert[j] += h
			ode_fn(du_pert, u_pert, nothing, 0.0)
			J_fd[:, j] .= (du_pert .- du_base) ./ h
		end

		# Physics block (T,ξ) × (T,ξ) should match well
		sv = ode_fn.sv
		n = n_cells
		physics_range = 1:(n * (1 + NC))

		for row in physics_range
			for col in physics_range
				a = J_analytical[row, col]
				fd = J_fd[row, col]
				if abs(fd) < 1e-4 && abs(a) < 1e-4
					continue
				end
				scale = max(abs(fd), abs(a))
				err = abs(a - fd) / scale
				@test err < 0.05
			end
		end

		# z-rows should be all zero (ALE coupling omitted)
		z_start = sv.z_offset + 1
		for row in z_start:(z_start + n_nodes - 1)
			@test all(J_analytical[row, :] .== 0.0)
		end
	end

	# ======================================================================
	# Test 3: Sparsity pattern includes z diagonal for W = I/(γdt) - J
	# ======================================================================
	@testset "Sparsity has z diagonal for W matrix" begin
		n_cells = 10
		NC = 2
		jac = jacobian_sparsity_ale_reduced(n_cells, NC)
		n_physics = n_cells * (1 + NC)
		n_nodes = n_cells + 1

		for k in 1:n_nodes
			z_idx = n_physics + k
			@test jac[z_idx, z_idx] != 0
		end
	end

	# ======================================================================
	# Test 4: use_reduced_jacobian validation errors
	# ======================================================================
	@testset "use_reduced_jacobian validation" begin
		material, mesh, bc_set = create_approx_jac_test_problem(n_cells = 5)

		problem = PyrolysisProblem(
			mesh = mesh,
			material = material,
			bc_set = bc_set,
			tspan = (0.0, 1.0),
			output_times = [0.5, 1.0],
		)

		@test_throws ErrorException solve(problem;
			use_reduced_jacobian = true, use_ale = false, use_analytical_jacobian = true)

		@test_throws ErrorException solve(problem;
			use_reduced_jacobian = true, use_ale = true, use_analytical_jacobian = false)
	end
end
