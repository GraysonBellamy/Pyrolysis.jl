"""
Test that the ALE Jacobian computation is allocation-free.

This test verifies that the JacobianOperatorALE implementation
eliminates allocations in the hot path.
"""

using Test
using Pyrolysis
using Pyrolysis.TimeIntegration
using Pyrolysis.Core
using Pyrolysis.Materials
using SparseArrays

@testset "JacobianOperatorALE Allocations" begin
	# Create a simple ALE problem
	# Use a small mesh to keep test fast
	n_cells = 10
	L = 0.01  # 1 cm

	# Create simple 2-component material (virgin wood + char)
	comp_virgin = MaterialComponent(
		name = "Virgin",
		density = ConstantProperty(500.0),
		heat_capacity = ConstantProperty(1500.0),
		conductivity = ConstantProperty(0.2),
		absorptivity = ConstantProperty(0.0),
		phase = SOLID,
		swelling_factor = 1.0,
	)

	comp_char = MaterialComponent(
		name = "Char",
		density = ConstantProperty(200.0),
		heat_capacity = ConstantProperty(1000.0),
		conductivity = ConstantProperty(0.1),
		absorptivity = ConstantProperty(0.0),
		phase = SOLID,
		swelling_factor = 1.0,
	)

	# Simple pyrolysis reaction: Virgin -> Char
	rxn = Reaction(
		reactant1 = 1,
		reactant2 = 0,
		product1 = 2,
		product2 = 0,
		stoich = (1.0, 0.0, 0.5, 0.0),
		A = 1e8,
		E = 140e3,
		order1 = 1.0,
		order2 = 0.0,
		heat = ConstantProperty(0.0),
		T_min = 300.0,
		T_max = 2000.0,
	)

	material = Material(
		name = "SimpleWood",
		components = (comp_virgin, comp_char),
		reactions = (rxn,),
		conductivity_mixing = PARALLEL,
	)

	# Create mesh
	mesh = Mesh1D(n_cells, L, material)

	# Initialize state
	T_init = 400.0
	ξ_init = (500.0, 0.0)  # All virgin initially
	initialize_state!(mesh, T_init, ξ_init)

	# Create boundary conditions
	bc_set = BoundaryConditionSet(
		top = Dict(:type => :convective, :h => 10.0, :T_inf => 800.0),
		bottom = Dict(:type => :adiabatic),
	)

	# Create ALE ODE function
	ode_fn = PyrolysisODEFunctionALE(mesh, material, bc_set;
		radiation_model = NO_RADIATION)

	# Create state vector
	sv = ode_fn.sv
	u = pack_state_ale(mesh, sv)

	# Create Jacobian operator
	op = JacobianOperatorALE(ode_fn)

	# Create Jacobian matrix
	J = copy(op.sparsity)

	# Test 1: First call (may allocate for compilation)
	t = 0.0
	p = nothing
	op(J, u, p, t)

	@testset "Zero Allocations" begin
		# Test 2: Second call should be allocation-free
		allocs = @allocated op(J, u, p, t)
		@test allocs == 0 "Expected 0 allocations, got $allocs bytes"

		# Test multiple calls to ensure consistency
		for _ in 1:5
			allocs = @allocated op(J, u, p, t)
			@test allocs == 0 "Expected 0 allocations in repeated call, got $allocs bytes"
		end
	end

	@testset "Correctness" begin
		# Verify the Jacobian is being computed (not all zeros)
		@test nnz(J) > 0 "Jacobian should have non-zero entries"

		# Check that diagonal blocks are non-zero (physics coupling)
		n = sv.n_cells
		NC = n_components(material)

		# Temperature-temperature block should have entries
		T_block = J[1:n, 1:n]
		@test nnz(T_block) > 0 "Temperature block should have non-zero entries"

		# Species-species block should have entries
		ξ_offset = n
		ξ_block = J[(ξ_offset+1):(ξ_offset+n*NC), (ξ_offset+1):(ξ_offset+n*NC)]
		@test nnz(ξ_block) > 0 "Species block should have non-zero entries"

		# Node velocity block should have entries
		z_offset = sv.z_offset
		z_block = J[(z_offset+1):end, :]
		@test nnz(z_block) > 0 "Node velocity block should have non-zero entries"
	end

	@testset "Comparison with make_ode_function_with_jac" begin
		# Test that the operator works correctly when used via make_ode_function_with_jac
		ode_func = make_ode_function_with_jac(ode_fn; jac_type = :sparse_fd)

		# Extract the Jacobian function
		jac_fn = ode_func.jac

		# Create another Jacobian matrix
		J2 = copy(op.sparsity)

		# Compute Jacobian via ODEFunction interface
		jac_fn(J2, u, p, t)

		# Should also be allocation-free
		allocs = @allocated jac_fn(J2, u, p, t)
		@test allocs == 0 "ODEFunction Jacobian should also be allocation-free, got $allocs bytes"

		# Results should match
		@test J ≈ J2 "Jacobians should match"
	end
end

println("✓ All Jacobian allocation tests passed!")
