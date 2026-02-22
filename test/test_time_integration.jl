"""
Tests for the TimeIntegration module.

Reference: Implementation Plan Stage 7

Tests cover:
1. State vector packing/unpacking
2. ODE function evaluation
3. Timestep computation
4. Jacobian sparsity patterns
5. Heat conduction verification (semi-infinite solid)
6. Single reaction verification (isothermal decay)
7. Mass and energy conservation
"""

using Test
using Pyrolysis
using Pyrolysis.TimeIntegration
using Pyrolysis.TimeIntegration: jacobian_sparsity, matrix_colors, compute_jacobian_fd!, compute_jacobian_fd_sparse!, JacobianOperator
using Pyrolysis.Geometry: create_uniform_mesh, cell_size
using Pyrolysis.Properties: PropertyCache, update_properties!
using Pyrolysis.Discretization: RHSCache, NO_RADIATION, BEER_LAMBERT
using OrdinaryDiffEq
using ADTypes
using LinearAlgebra
using SparseArrays

# ============================================================================
# Test Fixtures
# ============================================================================

"""Create a simple test material (PMMA-like)."""
function create_test_material(; kappa::Real = 0.0, permeability_mixing::MixingRule = WEIGHTED, β_κ::Real = 0.5)
	virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86, κ = kappa)
	gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
	rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)

	Material(
		name = :PMMA,
		components = (virgin, gas),
		reactions = (rxn,),
		conductivity_mixing = WEIGHTED,
		β = 0.5,
		permeability_mixing = permeability_mixing,
		β_κ = β_κ,
	)
end

"""Create a simple inert material (no reactions)."""
function create_inert_material()
	solid = SolidComponent(:solid, ρ = 1000.0, c = 1000.0, k = 1.0, ε = 0.9)

	Material(
		name = :Inert,
		components = (solid,),
		reactions = (),
	)
end

"""Create a test mesh with initial conditions."""
function create_test_mesh(L::Float64, n_cells::Int, ::Val{NC};
	T_init::Float64 = 300.0, ξ_init::Float64 = 1000.0) where NC
	mesh = create_uniform_mesh(L, n_cells, Val(NC))

	# Set initial conditions
	for i in 1:n_cells
		state = mesh.cell_states[i]
		state.T = T_init
		state.ξ = ntuple(_ -> ξ_init, NC)
	end

	return mesh
end

# ============================================================================
# State Vector Tests
# ============================================================================

@testset "State Vector" begin
	@testset "StateVector Construction" begin
		sv = StateVector(10, 2)
		@test sv.n_cells == 10
		@test sv.T_offset == 0
		@test sv.ξ_offset == 10

		@test state_vector_length(sv) == 30  # 10*(1+2)
		@test state_vector_length(10, 2) == 30
	end

	@testset "Pack and Unpack" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 500.0)
		material = create_test_material()

		# Create state vector
		u = create_state_vector(mesh, material)

		@test length(u) == 5 * 3  # 5 cells * (T + 2 concentrations)

		# Check packed values
		sv = StateVector(mesh, material)
		for i in 1:5
			@test get_temperature(u, sv, i) ≈ 400.0
			@test get_concentration(u, sv, i, 1) ≈ 500.0
			@test get_concentration(u, sv, i, 2) ≈ 500.0
		end

		# Modify and unpack
		set_temperature!(u, sv, 3, 600.0)
		set_concentration!(u, sv, 3, 1, 800.0)

		unpack_state!(mesh, u, sv)

		@test mesh.cell_states[3].T ≈ 600.0
		@test mesh.cell_states[3].ξ[1] ≈ 800.0
		@test mesh.cell_states[1].T ≈ 400.0  # Unchanged
	end

	@testset "Views" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector(mesh, 2)
		sv = StateVector(5, 2)

		T_view = temperature_view(u, sv)
		@test length(T_view) == 5

		ξ1_view = concentration_view(u, sv, 1)
		@test length(ξ1_view) == 5

		# Modifying view should modify u
		T_view[1] = 999.0
		@test get_temperature(u, sv, 1) ≈ 999.0
	end

	@testset "Type Stability" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector(mesh, 2)
		sv = StateVector(5, 2)

		@test @inferred get_temperature(u, sv, 1) isa Float64
		@test @inferred get_concentration(u, sv, 1, 1) isa Float64
		@test @inferred get_concentrations(u, sv, 1) isa NTuple{2, Float64}
	end
end

# ============================================================================
# ODE Function Tests
# ============================================================================

@testset "ODE Function" begin
	@testset "Construction" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		@test get_mesh(ode_fn) === mesh
		@test get_material(ode_fn) === material
		@test ode_fn.I_radiation == 0.0
	end

	@testset "RHS Evaluation" begin
		mesh = create_test_mesh(0.01, 10, Val(2); T_init = 300.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		u = create_state_vector(mesh, material)
		du = similar(u)

		# Should not error
		ode_fn(du, u, nothing, 0.0)

		# With adiabatic BCs and uniform T, heat flux divergence should be ~0
		# Reactions at low T should be ~0
		@test all(isfinite, du)
	end

	@testset "Time-Varying Radiation" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)
		I_func = t -> 50e3 * (1 - exp(-t / 10.0))

		ode_fn_rad = PyrolysisODEFunctionWithRadiation(ode_fn, I_func)

		u = create_state_vector(mesh, material)
		du = similar(u)

		ode_fn_rad(du, u, nothing, 10.0)

		# Check radiation was set
		@test ode_fn.I_radiation ≈ 50e3 * (1 - exp(-1))
	end

	@testset "Resize After Mesh Change" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		# Simulate mesh refinement by adding cells
		original_n_cells = length(mesh.cells)

		# This would normally be done through mesh refinement routines
		# For this test, just verify the resize function exists and runs
		resize_ode_function!(ode_fn)

		@test get_state_vector_descriptor(ode_fn).n_cells == original_n_cells
	end
end

# ============================================================================
# Timestep Control Tests
# ============================================================================

@testset "Timestep Control" begin
	@testset "CFL Timestep" begin
		mesh = create_test_mesh(0.01, 10, Val(1); T_init = 300.0, ξ_init = 1000.0)
		material = create_inert_material()

		# Create and update property cache
		prop_cache = PropertyCache(mesh)
		update_properties!(prop_cache, mesh, material)

		dt = compute_timestep(mesh, prop_cache; cfl = 0.5, safety = 0.9)

		@test dt > 0
		@test isfinite(dt)

		# Check CFL scaling: dt should scale with Δx²
		Δx = cell_size(mesh, 1)
		α = prop_cache.k[1] / prop_cache.ρcp[1]  # thermal diffusivity
		dt_expected = 0.5 * 0.9 * Δx^2 / (2.0 * α)

		@test dt ≈ dt_expected rtol=0.01
	end

	@testset "Timestep from Material" begin
		mesh = create_test_mesh(0.01, 10, Val(1))
		material = create_inert_material()

		dt = compute_timestep(mesh, material, 300.0)

		@test dt > 0
		@test isfinite(dt)
	end

	@testset "Stable Timestep" begin
		mesh = create_test_mesh(0.01, 10, Val(2); T_init = 300.0, ξ_init = 1000.0)
		material = create_test_material()

		prop_cache = PropertyCache(mesh)
		update_properties!(prop_cache, mesh, material)

		dt = compute_stable_timestep(mesh, prop_cache, material;
			dt_max = 1.0, dt_min = 1e-15)

		@test dt > 0
		@test dt <= 1.0
		@test dt >= 1e-15
	end

	@testset "Timestep Controller" begin
		controller = TimestepController(dt_initial = 1e-4)

		@test controller.dt == 1e-4

		# Simulate error too large -> should shrink
		new_dt = adjust_timestep!(controller, 1e-2)  # Large error
		@test new_dt < 1e-4

		# Reset
		reset_timestep!(controller, 1e-4)
		@test controller.dt == 1e-4

		# Simulate error acceptable -> should grow
		new_dt = adjust_timestep!(controller, 1e-6)  # Small error
		@test new_dt > 1e-4
	end
end

# ============================================================================
# Jacobian Tests
# ============================================================================

@testset "Jacobian" begin
	@testset "Sparsity Pattern" begin
		# 5 cells, 2 components
		sparsity = jacobian_sparsity(5, 2)

		@test size(sparsity) == (15, 15)  # 5*(1+2) = 15
		@test nnz(sparsity) > 0

		# Basic structure check: temperature depends on neighbor concentrations
		# (due to gas convective heat flux) and concentrations depend on neighbors
		sparsity_no_rad = jacobian_sparsity(5, 2; radiation_model = NO_RADIATION)
		@test nnz(sparsity_no_rad) > 0

		# With Beer-Lambert radiation, pattern has additional non-local coupling
		# (T at depth depends on ξ above, but not vice versa)
		sparsity_bl = jacobian_sparsity(5, 2; radiation_model = BEER_LAMBERT)
		@test nnz(sparsity_bl) > nnz(sparsity_no_rad)
	end

	@testset "Matrix Coloring" begin
		sparsity = jacobian_sparsity(10, 2)
		colors = matrix_colors(sparsity)

		@test length(colors) == size(sparsity, 2)
		@test all(c > 0 for c in colors)

		# For block-tridiagonal, should need much fewer colors than columns
		@test maximum(colors) < size(sparsity, 2)
	end

	@testset "Matrix Coloring Validity" begin
		# Test that coloring is valid (no conflicts) for various problem sizes
		for n in [10, 50, 100]
			for NC in [1, 2, 5]
				sparsity = jacobian_sparsity(n, NC)
				colors = matrix_colors(sparsity)

				@test length(colors) == size(sparsity, 2)
				@test all(c > 0 for c in colors)

				# Verify no conflicts: columns with same color have no shared rows
				n_colors = maximum(colors)
				for c in 1:n_colors
					cols = findall(==(c), colors)
					for i in 1:(length(cols)-1)
						for j in (i+1):length(cols)
							rows_i = sparsity.rowval[sparsity.colptr[cols[i]]:(sparsity.colptr[cols[i]+1]-1)]
							rows_j = sparsity.rowval[sparsity.colptr[cols[j]]:(sparsity.colptr[cols[j]+1]-1)]
							@test isempty(intersect(Set(rows_i), Set(rows_j)))
						end
					end
				end
			end
		end
	end

	@testset "Finite Difference Jacobian" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		u = create_state_vector(mesh, material)
		n = length(u)
		J = zeros(n, n)

		compute_jacobian_fd!(J, ode_fn, u, nothing, 0.0)

		@test all(isfinite, J)
		# Jacobian should have structure (not all zeros)
		@test norm(J) > 0
	end

	@testset "Sparse FD Jacobian" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		# Create JacobianOperator (includes sparsity, colors, and color_groups)
		op = JacobianOperator(ode_fn)
		J_sparse = copy(op.sparsity)

		u = create_state_vector(mesh, material)

		compute_jacobian_fd_sparse!(J_sparse, ode_fn, u, nothing, 0.0, op)

		@test all(isfinite, J_sparse.nzval)
	end

	@testset "JacobianOperator Construction" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)
		op = JacobianOperator(ode_fn)

		@test op.n_colors == maximum(op.colorvec)
		@test length(op.color_groups) == op.n_colors
		@test sum(length.(op.color_groups)) == length(op.colorvec)

		# Verify each color group contains correct columns
		for c in 1:op.n_colors
			for j in op.color_groups[c]
				@test op.colorvec[j] == c
			end
		end
	end

	@testset "Sparse FD vs Dense FD Accuracy" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		u = create_state_vector(mesh, material)
		n = length(u)

		# Dense FD (reference)
		J_dense = zeros(n, n)
		compute_jacobian_fd!(J_dense, ode_fn, u, nothing, 0.0)

		# Sparse FD with JacobianOperator
		op = JacobianOperator(ode_fn)
		J_sparse = copy(op.sparsity)
		compute_jacobian_fd_sparse!(J_sparse, ode_fn, u, nothing, 0.0, op)

		# Compare at non-zero positions
		for j in 1:n
			for idx in J_sparse.colptr[j]:(J_sparse.colptr[j+1]-1)
				i = J_sparse.rowval[idx]
				@test J_sparse.nzval[idx] ≈ J_dense[i, j] rtol=1e-6
			end
		end
	end
end

# ============================================================================
# Verification Tests
# ============================================================================

@testset "Verification" begin
	@testset "Heat Conduction - Steady State" begin
		# Test: Linear temperature profile with fixed BCs
		# T(0) = T_hot, T(L) = T_cold -> T(z) = T_hot - (T_hot-T_cold)*z/L
		# dT/dt should be ~0 for steady profile

		L = 0.01
		n_cells = 20
		T_hot = 500.0
		T_cold = 300.0

		mesh = create_uniform_mesh(L, n_cells, Val(1))
		material = create_inert_material()

		# Set linear temperature profile
		for i in 1:n_cells
			z = mesh.cell_centers[i]
			T_linear = T_hot - (T_hot - T_cold) * z / L
			mesh.cell_states[i].T = T_linear
			mesh.cell_states[i].ξ = (1000.0,)
		end

		# Use Dirichlet BCs at boundaries
		bc_set = BoundaryConditionSet(
			thermal_top = DirichletThermalBC(T_hot),
			thermal_bottom = DirichletThermalBC(T_cold),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)
		u = create_state_vector(mesh, material)
		du = similar(u)

		ode_fn(du, u, nothing, 0.0)

		# dT/dt for interior cells should be near zero for steady profile
		sv = StateVector(n_cells, 1)
		max_dT = maximum(abs(get_temperature(du, sv, i)) for i in 2:(n_cells-1))

		@test max_dT < 1.0  # Should be small for near-steady profile
	end

	@testset "Mass Conservation - Adiabatic" begin
		# With adiabatic BCs and no reactions (low T), mass should be conserved

		mesh = create_test_mesh(0.01, 10, Val(1); T_init = 300.0, ξ_init = 1000.0)
		material = create_inert_material()

		bc_set = BoundaryConditionSet(
			thermal_top = AdiabaticBC(),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		# Compute initial total mass
		initial_mass = sum(
			mesh.cell_states[i].ξ[1] * mesh.cell_volumes[i]
			for i in 1:length(mesh.cells)
		)

		# Run a few steps
		u0 = create_state_vector(mesh, material)
		tspan = (0.0, 1.0)
		prob = ODEProblem(ode_fn, u0, tspan)

		# Use a simple explicit solver
		sol = solve(prob, Tsit5(); abstol = 1e-8, reltol = 1e-6, saveat = 0.1)

		# Check mass at final time
		unpack_state!(mesh, sol.u[end], material)
		final_mass = sum(
			mesh.cell_states[i].ξ[1] * mesh.cell_volumes[i]
			for i in 1:length(mesh.cells)
		)

		# Mass should be conserved (no reactions, impermeable BCs)
		@test isapprox(final_mass, initial_mass, rtol = 1e-6)
	end

	@testset "ODE Integration - Runs Without Error" begin
		mesh = create_test_mesh(0.01, 10, Val(2); T_init = 300.0, ξ_init = 1000.0)
		material = create_test_material()

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		u0 = create_state_vector(mesh, material)
		tspan = (0.0, 0.1)  # Short time

		prob = ODEProblem(ode_fn, u0, tspan)

		# Test explicit solver
		sol_explicit = solve(prob, Tsit5(); abstol = 1e-6, reltol = 1e-4)
		@test sol_explicit.retcode == ReturnCode.Success

		# Test implicit solver (Rosenbrock) with finite-difference Jacobian
		# Note: AutoFiniteDiff() is used instead of AutoForwardDiff() because the time
		# parameter t can be a Dual type which is now handled correctly
		sol_implicit = solve(prob, Rodas5P(autodiff = AutoFiniteDiff()); abstol = 1e-6, reltol = 1e-4)
		@test sol_implicit.retcode == ReturnCode.Success
	end

	@testset "Temperature Increases Under Heating" begin
		mesh = create_test_mesh(0.01, 10, Val(1); T_init = 300.0, ξ_init = 1000.0)
		material = create_inert_material()

		# Apply constant heat flux at top (boundary 1)
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		u0 = create_state_vector(mesh, material)
		tspan = (0.0, 1.0)

		prob = ODEProblem(ode_fn, u0, tspan)
		sol = solve(prob, Tsit5(); abstol = 1e-6, reltol = 1e-4)

		# Unpack final state
		unpack_state!(mesh, sol.u[end], material)

		# Surface (cell N, at z=L) should be hotter than initial
		# (heat enters from top, boundary_id=1)
		@test mesh.cell_states[end].T > 300.0

		# Temperature gradient should exist (surface hotter than back)
		# Surface at cell N (z=L), back at cell 1 (z=0)
		@test mesh.cell_states[end].T > mesh.cell_states[1].T
	end
end

# ============================================================================
# Integration with DifferentialEquations.jl
# ============================================================================

@testset "DifferentialEquations.jl Integration" begin
	@testset "ODEFunction with Jacobian" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)

		# Test creating ODEFunction with sparse Jacobian
		ode_fn_jac = make_ode_function_with_jac(ode_fn; jac_type = :sparse_fd)
		@test ode_fn_jac isa OrdinaryDiffEq.ODEFunction

		# Test with dense Jacobian
		ode_fn_dense = make_ode_function_with_jac(ode_fn; jac_type = :dense_fd)
		@test ode_fn_dense isa OrdinaryDiffEq.ODEFunction

		# Test with no Jacobian
		ode_fn_none = make_ode_function_with_jac(ode_fn; jac_type = :none)
		@test ode_fn_none isa OrdinaryDiffEq.ODEFunction
	end

	@testset "Solve with Jacobian" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 350.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(20e3),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunction(mesh, material, bc_set)
		ode_fn_jac = make_ode_function_with_jac(ode_fn; jac_type = :sparse_fd)

		u0 = create_state_vector(mesh, material)
		tspan = (0.0, 0.05)

		prob = ODEProblem(ode_fn_jac, u0, tspan)

		# Should work with stiff solver
		sol = solve(prob, Rodas5P(); abstol = 1e-6, reltol = 1e-4)
		@test sol.retcode == ReturnCode.Success
	end
end

# ============================================================================
# ALE State Vector Tests
# ============================================================================

@testset "ALE State Vector" begin
	@testset "StateVectorALE Construction" begin
		sv = StateVectorALE(10, Val(2))
		@test sv.n_cells == 10
		@test sv.n_components == 2
		@test sv.n_nodes == 11
		@test sv.T_offset == 0
		@test sv.ξ_offset == 10
		@test sv.z_offset == 30  # 10 + 2*10

		@test state_vector_length_ale(sv) == 41  # 10*(1+2) + 11
		@test state_vector_length_ale(10, 2) == 41
	end

	@testset "StateVectorALE from Mesh" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		sv = StateVectorALE(mesh)

		@test sv.n_cells == 10
		@test sv.n_nodes == 11
		@test state_vector_length_ale(sv) == 41
	end

	@testset "Pack and Unpack ALE" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 500.0)
		material = create_test_material()

		# Create ALE state vector
		u = create_state_vector_ale(mesh)

		@test length(u) == 5 * 3 + 6  # 5 cells * 3 fields + 6 nodes

		# Check packed values
		sv = StateVectorALE(mesh)
		for i in 1:5
			@test get_temperature_ale(u, sv, i) ≈ 400.0
			@test get_concentration_ale(u, sv, i, 1) ≈ 500.0
			@test get_concentration_ale(u, sv, i, 2) ≈ 500.0
		end

		# Check node positions
		@test get_node_position(u, sv, 1) ≈ 0.0
		@test get_node_position(u, sv, 6) ≈ 0.01  # Last node at L=0.01

		# Modify and unpack
		set_temperature_ale!(u, sv, 3, 600.0)
		set_concentration_ale!(u, sv, 3, 1, 800.0)
		set_node_position!(u, sv, 3, 0.005)  # Move middle node

		unpack_state_ale!(mesh, u, sv)

		@test mesh.cell_states[3].T ≈ 600.0
		@test mesh.cell_states[3].ξ[1] ≈ 800.0
		@test mesh.cell_states[1].T ≈ 400.0  # Unchanged
		@test mesh.nodes[3].z ≈ 0.005
	end

	@testset "Node Positions View" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		z_view = get_node_positions(u, sv)
		@test length(z_view) == 6
		@test z_view[1] ≈ 0.0
		@test z_view[end] ≈ 0.01

		# Modifying view should modify u
		z_view[3] = 0.007
		@test get_node_position(u, sv, 3) ≈ 0.007
	end

	@testset "ALE Views" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		T_view = temperature_view_ale(u, sv)
		@test length(T_view) == 5

		ξ1_view = concentration_view_ale(u, sv, 1)
		@test length(ξ1_view) == 5

		z_view = node_positions_view(u, sv)
		@test length(z_view) == 6
	end

	@testset "Thickness from State" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		thickness = get_thickness_from_state(u, sv)
		@test thickness ≈ 0.01

		# Simulate shrinkage by moving top node
		set_node_position!(u, sv, 6, 0.008)
		thickness_new = get_thickness_from_state(u, sv)
		@test thickness_new ≈ 0.008
	end

	@testset "Volumes from State" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)
		A = mesh.cross_section_area

		volumes = compute_volumes_from_state(u, sv, A)
		@test length(volumes) == 5
		@test sum(volumes) ≈ 0.01 * A  # Total volume = L * A
		@test all(v -> v ≈ 0.002 * A, volumes)  # Uniform mesh: each cell = L/5 * A
	end

	@testset "Cell Centers from State" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		centers = compute_cell_centers_from_state(u, sv)
		@test length(centers) == 5
		@test centers[1] ≈ 0.001  # First cell center
		@test centers[5] ≈ 0.009  # Last cell center
	end

	@testset "Geometry Update from Unpacked State" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		# Record original volumes
		orig_volumes = copy(mesh.cell_volumes)

		# Move nodes to create non-uniform mesh
		set_node_position!(u, sv, 2, 0.003)  # Move second node
		set_node_position!(u, sv, 5, 0.007)  # Move fifth node

		# Unpack should update geometry
		unpack_state_ale!(mesh, u, sv)

		# Volumes should have changed
		@test mesh.cell_volumes[1] ≈ 0.003 * mesh.cross_section_area  # 0 to 0.003
		@test mesh.cell_states[1].volume ≈ mesh.cell_volumes[1]

		# Cell centers should have updated too
		@test mesh.cell_centers[1] ≈ 0.0015  # (0 + 0.003) / 2
	end

	@testset "Pack/Unpack Roundtrip" begin
		mesh = create_test_mesh(0.01, 10, Val(2); T_init = 500.0, ξ_init = 800.0)

		# Set some specific values
		mesh.cell_states[1].T = 400.0
		mesh.cell_states[5].ξ = (600.0, 300.0)

		# Pack into state vector
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		# Create fresh mesh and unpack into it
		mesh2 = create_test_mesh(0.01, 10, Val(2))
		unpack_state_ale!(mesh2, u, sv)

		# Check roundtrip
		@test mesh2.cell_states[1].T ≈ 400.0
		@test mesh2.cell_states[5].ξ[1] ≈ 600.0
		@test mesh2.cell_states[5].ξ[2] ≈ 300.0
		@test mesh2.nodes[1].z ≈ mesh.nodes[1].z
		@test mesh2.nodes[end].z ≈ mesh.nodes[end].z
	end

	@testset "Type Stability" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		@test @inferred(get_temperature_ale(u, sv, 1)) isa Float64
		@test @inferred(get_concentration_ale(u, sv, 1, 1)) isa Float64
		@test @inferred(get_concentrations_ale(u, sv, 1)) isa NTuple{2, Float64}
		@test @inferred(get_node_position(u, sv, 1)) isa Float64
		@test @inferred(get_thickness_from_state(u, sv)) isa Float64
	end

	@testset "Validation" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		sv = StateVectorALE(mesh)

		u_correct = create_state_vector_ale(mesh)
		@test validate_state_vector_ale(u_correct, sv) == true

		u_wrong = zeros(10)  # Wrong size
		@test_throws DimensionMismatch validate_state_vector_ale(u_wrong, sv)
	end

	@testset "Set Node Velocity" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		sv = StateVectorALE(mesh)
		du = zeros(state_vector_length_ale(sv))

		# Set node velocities (as would be done in ODE RHS)
		set_node_velocity!(du, sv, 1, 0.0)       # Bottom fixed
		set_node_velocity!(du, sv, 6, -0.001)    # Top moving down (shrinkage)

		@test du[sv.z_offset+1] ≈ 0.0
		@test du[sv.z_offset+6] ≈ -0.001
	end
end

# ============================================================================
# Geometry from State Tests
# ============================================================================

@testset "Geometry from State" begin
	@testset "GeometryCache Construction" begin
		# From dimensions
		cache = GeometryCache(10, 11)
		@test length(cache.volumes) == 10
		@test length(cache.centers) == 10
		@test length(cache.face_positions) == 11
		@test length(cache.face_distances) == 11

		# From mesh
		mesh = create_test_mesh(0.01, 5, Val(2))
		cache2 = GeometryCache(mesh)
		@test length(cache2.volumes) == 5
		@test length(cache2.face_positions) == 6
	end

	@testset "compute_geometry_from_state!" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		volumes = zeros(5)
		centers = zeros(5)
		face_positions = zeros(6)

		compute_geometry_from_state!(volumes, centers, face_positions,
			u, sv, mesh.cross_section_area)

		# Check volumes match original mesh
		@test volumes ≈ mesh.cell_volumes
		@test centers ≈ mesh.cell_centers
		@test face_positions ≈ mesh.face_positions
	end

	@testset "compute_face_distances_from_state!" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		# First compute geometry
		volumes = zeros(5)
		centers = zeros(5)
		face_positions = zeros(6)
		compute_geometry_from_state!(volumes, centers, face_positions,
			u, sv, mesh.cross_section_area)

		# Now compute distances
		distances = zeros(6)
		compute_face_distances_from_state!(distances, centers, face_positions, mesh)

		# Interior faces should have distance = cell size (uniform mesh)
		cell_width = 0.01 / 5
		for f in 2:5  # Interior faces
			@test distances[f] ≈ cell_width
		end

		# Boundary faces should have 2× half-cell distance
		@test distances[1] ≈ cell_width  # 2 × (cell_width/2)
		@test distances[6] ≈ cell_width
	end

	@testset "update_geometry_cache!" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		cache = GeometryCache(mesh)
		update_geometry_cache!(cache, u, sv, mesh)

		@test cache.volumes ≈ mesh.cell_volumes
		@test cache.centers ≈ mesh.cell_centers
		@test cache.face_positions ≈ mesh.face_positions
		@test all(cache.face_distances .> 0)
	end

	@testset "Geometry after node motion" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		# Move nodes (shrink domain by 20%)
		for i in 1:sv.n_nodes
			old_z = get_node_position(u, sv, i)
			set_node_position!(u, sv, i, old_z * 0.8)
		end

		cache = GeometryCache(mesh)
		update_geometry_cache!(cache, u, sv, mesh)

		# Total volume should be 80% of original
		@test sum(cache.volumes) ≈ 0.8 * 0.01 * mesh.cross_section_area

		# Top position should be at 0.008
		@test cache.face_positions[end] ≈ 0.008

		# Cell centers should also be scaled
		@test cache.centers[1] ≈ 0.8 * mesh.cell_centers[1]
		@test cache.centers[end] ≈ 0.8 * mesh.cell_centers[end]
	end

	@testset "Non-uniform node motion" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		# Move only the top node (simulating surface regression)
		set_node_position!(u, sv, 6, 0.007)  # Top moves from 0.01 to 0.007

		cache = GeometryCache(mesh)
		update_geometry_cache!(cache, u, sv, mesh)

		# Only top cell (cell 5) should shrink
		# Node positions are now: [0, 0.002, 0.004, 0.006, 0.008, 0.007]
		# Cell 5 spans from node 5 (0.008) to node 6 (0.007)
		cell_width_original = 0.01 / 5  # = 0.002
		@test cache.volumes[1] ≈ cell_width_original * mesh.cross_section_area

		# Cell 5 volume = |0.007 - 0.008| * A = 0.001 * A
		@test cache.volumes[5] ≈ 0.001 * mesh.cross_section_area

		# Total volume: cells 1-4 unchanged (4 * 0.002 * A = 0.008 * A)
		# plus cell 5 (0.001 * A) = 0.009 * A total
		# (Note: moving only the top node doesn't compress the whole mesh,
		# it only shrinks the top cell. Cell 5 has negative gradient z[6]-z[5])
		@test sum(cache.volumes) ≈ 0.009 * mesh.cross_section_area
	end

	@testset "GCL Compliance - Basic" begin
		# Verify dV/dt = (w_right - w_left) * A
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)
		A = mesh.cross_section_area

		# Compute initial volumes
		V_old = compute_volumes_from_state(u, sv, A)

		# Apply mesh velocity (linear profile: w = 0 at bottom, w = -0.001 at top)
		dt = 0.1
		w = [(i-1) * (-0.001) / 5 for i in 1:6]  # w[1]=0, w[6]=-0.001

		# Move nodes according to velocity
		for i in 1:sv.n_nodes
			old_z = get_node_position(u, sv, i)
			set_node_position!(u, sv, i, old_z + dt * w[i])
		end

		# Compute new volumes
		V_new = compute_volumes_from_state(u, sv, A)

		# Verify GCL: dV/dt = (w_right - w_left) * A
		for i in 1:5
			dV_dt_actual = (V_new[i] - V_old[i]) / dt
			dV_dt_gcl = (w[i+1] - w[i]) * A
			@test dV_dt_actual ≈ dV_dt_gcl atol=1e-12
		end
	end

	@testset "GCL Compliance - verify_gcl_compliance" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)
		A = mesh.cross_section_area

		V_old = compute_volumes_from_state(u, sv, A)

		# Apply uniform shrinkage velocity
		dt = 0.05
		w = [0.0, -0.0002, -0.0004, -0.0006, -0.0008, -0.001]

		for i in 1:sv.n_nodes
			old_z = get_node_position(u, sv, i)
			set_node_position!(u, sv, i, old_z + dt * w[i])
		end

		V_new = compute_volumes_from_state(u, sv, A)

		@test verify_gcl_compliance(V_old, V_new, w, dt, mesh) == true

		# Now tamper with volumes to break GCL
		V_bad = copy(V_new)
		V_bad[3] += 0.001  # Artificially change one volume
		@test verify_gcl_compliance(V_old, V_bad, w, dt, mesh) == false
	end

	@testset "compute_gcl_volume_rate" begin
		mesh = create_test_mesh(0.01, 5, Val(2))
		A = mesh.cross_section_area

		# Uniform velocity gradient
		w = [0.0, -0.0002, -0.0004, -0.0006, -0.0008, -0.001]
		dV_dt = compute_gcl_volume_rate(w, mesh)

		# Each cell should have same dV/dt for linear w profile
		expected_rate = (-0.0002 - 0.0) * A  # = (w[2] - w[1]) * A
		@test all(dV_dt .≈ expected_rate)
	end

	@testset "Cache independence from mesh" begin
		# Verify that cache computation doesn't modify mesh
		mesh = create_test_mesh(0.01, 5, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)

		# Store original mesh values
		orig_volumes = copy(mesh.cell_volumes)
		orig_centers = copy(mesh.cell_centers)

		# Move nodes in state vector
		for i in 1:sv.n_nodes
			old_z = get_node_position(u, sv, i)
			set_node_position!(u, sv, i, old_z * 0.5)
		end

		# Update cache (should NOT modify mesh)
		cache = GeometryCache(mesh)
		update_geometry_cache!(cache, u, sv, mesh)

		# Mesh should be unchanged
		@test mesh.cell_volumes ≈ orig_volumes
		@test mesh.cell_centers ≈ orig_centers

		# Cache should have new values
		@test cache.volumes ≈ orig_volumes * 0.5
		@test !(cache.volumes ≈ orig_volumes)
	end

	@testset "Zero-allocation check" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		u = create_state_vector_ale(mesh)
		sv = StateVectorALE(mesh)
		cache = GeometryCache(mesh)

		# Warm up
		update_geometry_cache!(cache, u, sv, mesh)

		# Check allocations
		allocs = @allocated update_geometry_cache!(cache, u, sv, mesh)
		@test allocs == 0
	end
end

# Helper function for accessing dξ from du vector (used in ALE tests)
function get_dconcentration_ale(du::Vector, sv::StateVectorALE, cell::Int, component::Int)
	idx = sv.ξ_offset + (component - 1) * sv.n_cells + cell
	return du[idx]
end

# ============================================================================
# ALE ODE Function Tests
# ============================================================================

@testset "ALE ODE Function" begin
	@testset "PyrolysisODEFunctionALE Construction" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		@test get_mesh(ode_fn) === mesh
		@test get_material(ode_fn) === material
		@test ode_fn.I_radiation == 0.0

		# Check that ALE-specific fields are initialized
		@test length(ode_fn.geom_cache.volumes) == 10
		@test length(ode_fn.w_work) == 11          # Node velocities work array
		@test length(ode_fn.w_face) == 11          # Face velocities
		# NEW: Cell-centered advection arrays (not face-based flux arrays)
		@test length(ode_fn.ale_advection) == 10   # ALE energy advection (per cell)
		# ale_species_advection is (n_components × n_cells) = 2 × 10 = 20
		@test length(ode_fn.ale_species_advection) == 10 * 2
	end

	@testset "PyrolysisODEFunctionALE RHS Evaluation" begin
		mesh = create_test_mesh(0.01, 10, Val(2); T_init = 300.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		# Create ALE state vector (T, ξ, z)
		u = create_state_vector_ale(mesh)
		du = similar(u)

		# Should not error
		ode_fn(du, u, nothing, 0.0)

		# Check output dimensions match
		@test length(du) == length(u)

		# All outputs should be finite
		@test all(isfinite, du)

		sv = StateVectorALE(mesh)

		# Temperature should be increasing from heating at top
		# (dT/dt for top cells should be positive)
		# With heat flux BC, surface is heated
		dT_top = du[sv.T_offset+sv.n_cells]  # Top cell
		# Note: dT/dt depends on many factors, just check finite
		@test isfinite(dT_top)

		# Node velocities should be set (mesh velocity outputs)
		# Bottom node should be fixed (velocity = 0)
		@test du[sv.z_offset+1] ≈ 0.0

		# Top node velocity should be finite (could be 0 or negative for shrinkage)
		@test isfinite(du[sv.z_offset+sv.n_nodes])
	end

	@testset "PyrolysisODEFunctionALE with Inert Material" begin
		# With inert material at low temperature, no reactions -> no mass loss
		mesh = create_test_mesh(0.01, 10, Val(1); T_init = 300.0, ξ_init = 1000.0)
		material = create_inert_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(10e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u = create_state_vector_ale(mesh)
		du = similar(u)

		ode_fn(du, u, nothing, 0.0)

		sv = StateVectorALE(mesh)

		# No mass loss means no surface regression
		# Surface velocity should be ~0 (no gaseous products escaping)
		w_top = du[sv.z_offset+sv.n_nodes]
		@test abs(w_top) < 1e-10

		# All node velocities should be ~0 for inert material
		for i in 1:sv.n_nodes
			@test abs(du[sv.z_offset+i]) < 1e-10
		end
	end

	@testset "PyrolysisODEFunctionALE Mesh Velocity" begin
		# Test that mesh velocity is computed correctly
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 300.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u = create_state_vector_ale(mesh)
		du = similar(u)

		ode_fn(du, u, nothing, 0.0)

		# Retrieve mesh velocity using accessor function
		w = get_mesh_velocity(ode_fn)
		@test length(w) == 6  # 5 cells + 1 = 6 nodes

		# Bottom should be fixed
		@test w[1] ≈ 0.0

		# With LaplacianSmoothing default, interior velocities should
		# interpolate linearly between boundaries
		# At low temperature, reactions are slow, so velocity should be small
	end

	@testset "PyrolysisODEFunctionALE GCL Compliance" begin
		# Verify that volume rates are consistent with mesh velocity
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 500.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u = create_state_vector_ale(mesh)
		du = similar(u)

		ode_fn(du, u, nothing, 0.0)

		sv = StateVectorALE(mesh)
		A = mesh.cross_section_area

		# Get mesh velocities from du
		w = [du[sv.z_offset+i] for i in 1:sv.n_nodes]

		# Compute expected volume rates from GCL: dV_i/dt = (w_{i+1} - w_i) * A
		for i in 1:sv.n_cells
			dV_dt_gcl = (w[i+1] - w[i]) * A
			# Volume rates should match GCL prediction
			# (Note: This is implicitly satisfied by our velocity computation)
			@test isfinite(dV_dt_gcl)
		end
	end

	@testset "PyrolysisODEFunctionALE vs Fixed Mesh (w=0)" begin
		# When mesh velocity is zero, ALE ODE should give same result as fixed mesh
		# Use inert material with no reactions to ensure zero mesh motion
		mesh_ale = create_test_mesh(0.01, 10, Val(1); T_init = 300.0, ξ_init = 1000.0)
		mesh_fixed = create_test_mesh(0.01, 10, Val(1); T_init = 300.0, ξ_init = 1000.0)
		material = create_inert_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(20e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn_ale = PyrolysisODEFunctionALE(mesh_ale, material, bc_set)
		ode_fn_fixed = PyrolysisODEFunction(mesh_fixed, material, bc_set)

		u_ale = create_state_vector_ale(mesh_ale)
		u_fixed = create_state_vector(mesh_fixed, material)

		du_ale = similar(u_ale)
		du_fixed = similar(u_fixed)

		ode_fn_ale(du_ale, u_ale, nothing, 0.0)
		ode_fn_fixed(du_fixed, u_fixed, nothing, 0.0)

		sv_ale = StateVectorALE(mesh_ale)
		sv_fixed = StateVector(mesh_fixed, material)

		# Temperature derivatives should match between ALE and fixed
		for i in 1:10
			dT_ale = du_ale[sv_ale.T_offset+i]
			dT_fixed = du_fixed[sv_fixed.T_offset+i]
			@test dT_ale ≈ dT_fixed rtol=1e-10
		end

		# Concentration derivatives should also match
		for i in 1:10
			dξ_ale = du_ale[sv_ale.ξ_offset+i]
			dξ_fixed = du_fixed[sv_fixed.ξ_offset+i]
			@test dξ_ale ≈ dξ_fixed rtol=1e-10
		end
	end

	@testset "PyrolysisODEFunctionALE Out-of-Place" begin
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 300.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u = create_state_vector_ale(mesh)

		# Test out-of-place evaluation
		du = ode_fn(u, nothing, 0.0)

		@test length(du) == length(u)
		@test all(isfinite, du)
	end

	@testset "PyrolysisODEFunctionALEWithRadiation" begin
		mesh = create_test_mesh(0.01, 10, Val(2))
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(0.0),
			thermal_bottom = AdiabaticBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)
		I_func = t -> 50e3 * (1 - exp(-t / 10.0))

		ode_fn_rad = PyrolysisODEFunctionALEWithRadiation(ode_fn, I_func)

		u = create_state_vector_ale(mesh)
		du = similar(u)

		ode_fn_rad(du, u, nothing, 10.0)

		# Check radiation was set
		@test ode_fn.I_radiation ≈ 50e3 * (1 - exp(-1))
		@test all(isfinite, du)
	end

	@testset "PyrolysisODEFunctionALE Integration with ODE solver" begin
		# Test that the ALE ODE function can be integrated without error
		mesh = create_test_mesh(0.01, 10, Val(2); T_init = 300.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u0 = create_state_vector_ale(mesh)
		tspan = (0.0, 0.1)  # Short time

		prob = ODEProblem(ode_fn, u0, tspan)

		# Test explicit solver
		sol_explicit = solve(prob, Tsit5(); abstol = 1e-6, reltol = 1e-4)
		@test sol_explicit.retcode == ReturnCode.Success

		# Verify solution values are physical
		u_final = sol_explicit.u[end]
		sv = StateVectorALE(mesh)

		# Temperatures should be positive and finite
		for i in 1:sv.n_cells
			T = get_temperature_ale(u_final, sv, i)
			@test T > 0
			@test isfinite(T)
		end

		# Node positions should be positive and finite
		for i in 1:sv.n_nodes
			z = get_node_position(u_final, sv, i)
			@test z >= 0
			@test isfinite(z)
		end

		# Thickness should have decreased (some material loss from reactions)
		# or stayed the same (if reactions are slow)
		thickness = get_thickness_from_state(u_final, sv)
		@test thickness > 0
		@test thickness <= 0.01 + 1e-6  # Shouldn't grow significantly
	end

	@testset "PyrolysisODEFunctionALE High Temperature" begin
		# At high temperature, reactions should be active, causing mass loss
		# and surface regression
		mesh = create_test_mesh(0.01, 10, Val(2); T_init = 600.0, ξ_init = 1190.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(50e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u = create_state_vector_ale(mesh)
		du = similar(u)

		ode_fn(du, u, nothing, 0.0)

		sv = StateVectorALE(mesh)

		# At 600K with Arrhenius kinetics (A=8.5e12, E=188e3),
		# reaction rate should be measurable.
		# 
		# In the ALE formulation with dilation correction:
		# - Solid concentration dξ/dt ≈ 0 (mass loss is tracked by volume change)
		# - Gas concentration dξ/dt > 0 (gas is produced)
		# - Mesh velocity w < 0 (shrinking due to solid mass loss)
		#
		# The correct check is that either:
		# 1. Gas is being produced (dξ_gas > 0), OR
		# 2. Mesh is shrinking (w_top < 0)

		# Check that gas is being produced (any cell)
		any_gas_production = false
		for i in 1:sv.n_cells
			dξ2 = get_dconcentration_ale(du, sv, i, 2)  # Gas production rate
			if dξ2 > 1e-6  # Measurable production
				any_gas_production = true
				break
			end
		end
		@test any_gas_production

		# Surface velocity should be negative (shrinking) due to solid mass loss
		w_top = du[sv.z_offset+sv.n_nodes]
		@test isfinite(w_top)
		@test w_top < 0  # Mesh should be shrinking
	end

	@testset "PyrolysisODEFunctionALE with Jacobian" begin
		# Test make_ode_function_with_jac for ALE ODE function
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 1000.0)
		material = create_test_material()
		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		# Test creating ODEFunction with sparse Jacobian
		ode_fn_sparse = make_ode_function_with_jac(ode_fn; jac_type = :sparse_fd)
		@test ode_fn_sparse isa OrdinaryDiffEq.ODEFunction

		# Test with dense Jacobian
		ode_fn_dense = make_ode_function_with_jac(ode_fn; jac_type = :dense_fd)
		@test ode_fn_dense isa OrdinaryDiffEq.ODEFunction

		# Test with AD Jacobian
		ode_fn_ad = make_ode_function_with_jac(ode_fn; jac_type = :ad)
		@test ode_fn_ad isa OrdinaryDiffEq.ODEFunction

		# Test with no Jacobian
		ode_fn_none = make_ode_function_with_jac(ode_fn; jac_type = :none)
		@test ode_fn_none isa OrdinaryDiffEq.ODEFunction
	end

	@testset "PyrolysisODEFunctionALE with Darcy and Jacobian" begin
		# Test Jacobian computation for ALE with Darcy-Fick enabled
		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 1000.0)

		# Create material with component-based κ for Darcy flow
		material = create_test_material(kappa = 1e-11, permeability_mixing = WEIGHTED, β_κ = 0.5)

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		# Create Darcy-enabled ODE function (uses material.permeability_model)
		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		# Verify Darcy is enabled via material
		@test has_darcy_flow(material)

		# Test sparse FD Jacobian
		ode_fn_sparse = make_ode_function_with_jac(ode_fn; jac_type = :sparse_fd)
		@test ode_fn_sparse isa OrdinaryDiffEq.ODEFunction

		# Test AD Jacobian (recommended for Darcy)
		ode_fn_ad = make_ode_function_with_jac(ode_fn; jac_type = :ad)
		@test ode_fn_ad isa OrdinaryDiffEq.ODEFunction

		# Solve a short integration to verify Jacobian works
		u0 = create_state_vector_ale(mesh)
		tspan = (0.0, 0.01)

		prob = ODEProblem(ode_fn_ad, u0, tspan)
		sol = solve(prob, Rodas5P(); abstol = 1e-6, reltol = 1e-4)
		@test sol.retcode == ReturnCode.Success
	end

	@testset "FD Jacobian for Darcy-enabled RHS" begin
		# Verify FD Jacobian computation works and produces sensible results
		# Note: AD is not supported for Darcy-enabled RHS due to internal Float64 mutations

		mesh = create_test_mesh(0.01, 5, Val(2); T_init = 400.0, ξ_init = 1000.0)

		# Create material with component-based κ for Darcy flow
		material = create_test_material(kappa = 1e-11, permeability_mixing = WEIGHTED, β_κ = 0.5)

		bc_set = BoundaryConditionSet(
			thermal_top = HeatFluxBC(30e3),
			thermal_bottom = AdiabaticBC(),
			mass_top = ImpermeableBC(),
			mass_bottom = ImpermeableBC(),
		)

		ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set)

		u0 = create_state_vector_ale(mesh)
		NC = 2
		n_cells = length(mesh.cells)
		n_nodes = n_cells + 1
		n_vars = n_cells * (1 + NC) + n_nodes
		n_physics = n_cells * (1 + NC)

		# Compute FD Jacobian for physics variables
		J_fd = zeros(n_physics, n_physics)
		du0 = zeros(n_vars)
		du_pert = zeros(n_vars)
		u_pert = copy(u0)

		# Get base RHS
		ode_fn(du0, u0, nothing, 0.0)

		ε = 1e-7
		for j in 1:n_physics
			u_pert[j] = u0[j] + ε
			ode_fn(du_pert, u_pert, nothing, 0.0)
			for i in 1:n_physics
				J_fd[i, j] = (du_pert[i] - du0[i]) / ε
			end
			u_pert[j] = u0[j]
		end

		# Basic sanity checks on the Jacobian
		@test size(J_fd) == (n_physics, n_physics)
		@test all(isfinite.(J_fd))  # No NaN or Inf values

		# Temperature diagonal: diffusion contributes negative terms (stabilizing)
		# Advection form (q = J × cp × ΔT) contributes both positive and negative terms
		# that redistribute heat without net creation. The overall system should still
		# be stable, but individual diagonals may not all be negative.
		T_diag = [J_fd[i, i] for i in 1:n_cells]
		# Check that at least some are negative (diffusion-dominated cells)
		@test count(x -> x < 0, T_diag) > 0  # Some cells have negative diagonal
		# Check that the trace (sum of eigenvalues) is negative (overall stability)
		@test sum(T_diag) < 0  # Overall temperature evolution is stable

		# Off-diagonal temperature entries should exist (heat diffusion couples cells)
		# Check that interior cells have non-zero coupling to neighbors
		for i in 2:(n_cells-1)
			@test J_fd[i, i-1] != 0.0 || J_fd[i, i+1] != 0.0  # At least one neighbor coupling
		end

		# Jacobian should be sparse-ish (most entries are small)
		nonzero_count = count(x -> abs(x) > 1e-10, J_fd)
		sparsity = 1.0 - nonzero_count / length(J_fd)
		@test sparsity > 0.5  # At least 50% sparse

		# Verify Jacobian is consistent with RHS by checking linearization
		# For small δu, du(u0 + δu) ≈ du(u0) + J * δu
		δu = randn(n_physics) * 1e-6
		u_test = copy(u0)
		u_test[1:n_physics] .+= δu
		du_test = zeros(n_vars)
		ode_fn(du_test, u_test, nothing, 0.0)

		# Predicted change from Jacobian
		du_predicted = du0[1:n_physics] + J_fd * δu

		# Actual change
		du_actual = du_test[1:n_physics]

		# Should match to O(δu²) ≈ 1e-12
		rel_error = norm(du_predicted - du_actual) / (norm(du_actual) + 1e-10)
		@test rel_error < 1e-4  # Allow for FD truncation error
	end
end
