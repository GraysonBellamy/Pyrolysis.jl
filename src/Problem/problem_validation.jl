"""
Problem validation utilities.

Reference: Design Document Section 6.7
"""

# ============================================================================
# Main Validation Function
# ============================================================================

"""
	validate_problem(problem::PyrolysisProblem; strict::Bool=false)

Validate a PyrolysisProblem for consistency and correctness.

# Arguments
- `problem`: The problem to validate
- `strict=false`: If true, perform additional expensive checks

# Checks Performed
1. Mesh is properly initialized
2. Material is valid
3. Boundary conditions are defined for all boundaries
4. Initial conditions are physically reasonable
5. Time span is valid
6. Component count consistency

# Returns
`true` if all checks pass, throws an error with details otherwise.

# Example
```julia
problem = PyrolysisProblem(...)
validate_problem(problem)  # Throws if invalid
```
"""
function validate_problem(problem::PyrolysisProblem; strict::Bool = false)
	errors = String[]

	# Validate mesh
	mesh_errors = validate_mesh(problem.mesh)
	append!(errors, mesh_errors)

	# Validate material
	mat_errors = validate_material_for_problem(problem.material)
	append!(errors, mat_errors)

	# Validate boundary conditions
	bc_errors = validate_boundary_conditions(problem)
	append!(errors, bc_errors)

	# Validate time span
	tspan_errors = validate_tspan(problem.tspan, problem.dt_initial)
	append!(errors, tspan_errors)

	# Validate initial conditions match component count
	ic_errors = validate_initial_condition_funcs(problem)
	append!(errors, ic_errors)

	# Additional strict checks
	if strict
		strict_errors = validate_strict(problem)
		append!(errors, strict_errors)
	end

	# Report errors
	if !isempty(errors)
		error_msg = "Problem validation failed with $(length(errors)) error(s):\n"
		for (i, err) in enumerate(errors)
			error_msg *= "  $i. $err\n"
		end
		error(error_msg)
	end

	return true
end

# ============================================================================
# Component Validation Functions
# ============================================================================

"""
	validate_mesh(mesh) -> Vector{String}

Validate mesh structure.
"""
function validate_mesh(mesh)
	errors = String[]

	# Check mesh has cells
	if isempty(mesh.cells)
		push!(errors, "Mesh has no cells")
	end

	# Check mesh has nodes
	if isempty(mesh.nodes)
		push!(errors, "Mesh has no nodes")
	end

	# Check mesh has faces
	if isempty(mesh.faces)
		push!(errors, "Mesh has no faces")
	end

	# Check cell states exist
	if isempty(mesh.cell_states)
		push!(errors, "Mesh has no cell states")
	elseif length(mesh.cell_states) != length(mesh.cells)
		push!(errors, "Mesh cell_states count ($(length(mesh.cell_states))) doesn't match cells count ($(length(mesh.cells)))")
	end

	# Check geometry cache
	if isempty(mesh.cell_centers)
		push!(errors, "Mesh cell_centers not initialized")
	end
	if isempty(mesh.cell_volumes)
		push!(errors, "Mesh cell_volumes not initialized")
	end

	# Check boundary faces
	if isempty(mesh.boundary_faces)
		push!(errors, "Mesh has no boundary face definitions")
	end

	return errors
end

"""
	validate_material_for_problem(material) -> Vector{String}

Validate material for use in a problem.
"""
function validate_material_for_problem(material)
	errors = String[]

	NC = n_components(material)

	# Must have at least one component
	if NC < 1
		push!(errors, "Material must have at least one component")
	end

	# Check reactions reference valid components
	for (i, rxn) in enumerate(material.reactions)
		if rxn.reactant1 != 0 && (rxn.reactant1 < 1 || rxn.reactant1 > NC)
			push!(errors, "Reaction $i: reactant1 index $(rxn.reactant1) out of range [1, $NC]")
		end
		if rxn.reactant2 != 0 && (rxn.reactant2 < 1 || rxn.reactant2 > NC)
			push!(errors, "Reaction $i: reactant2 index $(rxn.reactant2) out of range [1, $NC]")
		end
		if rxn.product1 != 0 && (rxn.product1 < 1 || rxn.product1 > NC)
			push!(errors, "Reaction $i: product1 index $(rxn.product1) out of range [1, $NC]")
		end
		if rxn.product2 != 0 && (rxn.product2 < 1 || rxn.product2 > NC)
			push!(errors, "Reaction $i: product2 index $(rxn.product2) out of range [1, $NC]")
		end
	end

	return errors
end

"""
	validate_boundary_conditions(problem) -> Vector{String}

Validate boundary conditions are defined for all mesh boundaries.
"""
function validate_boundary_conditions(problem)
	errors = String[]

	mesh = problem.mesh
	bc_set = problem.bc_set

	# Get all boundary IDs from mesh
	mesh_boundary_ids = Set(keys(mesh.boundary_faces))

	# Check that we have BCs for the expected boundaries (1=top, 2=bottom)
	# BoundaryConditionSet always has all 4 BCs (thermal_top, thermal_bottom, mass_top, mass_bottom)
	# so we just verify the mesh has the expected boundaries
	expected_boundaries = Set([1, 2])  # Top and bottom
	missing_boundaries = setdiff(expected_boundaries, mesh_boundary_ids)
	if !isempty(missing_boundaries)
		push!(errors, "Mesh missing expected boundaries: $missing_boundaries (expected 1=top, 2=bottom)")
	end

	# Note: BoundaryConditionSet is always fully defined with all 4 BCs,
	# so no need to check for missing BCs like we did with Dicts

	return errors
end

"""
	validate_tspan(tspan, dt_initial) -> Vector{String}

Validate time span and initial timestep.
"""
function validate_tspan(tspan, dt_initial)
	errors = String[]

	t_start, t_end = tspan

	if t_start >= t_end
		push!(errors, "Invalid time span: t_start ($t_start) >= t_end ($t_end)")
	end

	if t_start < 0
		push!(errors, "Negative start time: t_start = $t_start")
	end

	if dt_initial <= 0
		push!(errors, "Non-positive initial timestep: dt_initial = $dt_initial")
	end

	if dt_initial > (t_end - t_start)
		push!(errors, "Initial timestep ($dt_initial) larger than time span ($(t_end - t_start))")
	end

	return errors
end

"""
	validate_initial_condition_funcs(problem) -> Vector{String}

Validate initial condition functions match component count.
"""
function validate_initial_condition_funcs(problem)
	errors = String[]

	NC = n_components(problem.material)

	# Check ξ_initial length
	if length(problem.ξ_initial) != NC
		push!(errors, "Number of initial concentration functions ($(length(problem.ξ_initial))) doesn't match component count ($NC)")
	end

	return errors
end

"""
	validate_strict(problem) -> Vector{String}

Perform additional expensive validation checks.
"""
function validate_strict(problem)
	errors = String[]

	mesh = problem.mesh
	material = problem.material

	# Check reaction stoichiometry mass balance
	for (i, rxn) in enumerate(material.reactions)
		if !Core.verify_stoichiometry(rxn; tol = 1e-10)
			ν = rxn.stoich
			reactant_mass = ν[1] + ν[2]
			product_mass = ν[3] + ν[4]
			err = abs(reactant_mass - product_mass)
			push!(errors, "Reaction $i ($(rxn.name)) violates mass balance: " *
						  "reactants=$(reactant_mass), products=$(product_mass), error=$(err)")
		end
	end

	# Test initial condition functions are callable
	try
		z_test = 0.0
		problem.T_initial(z_test)
	catch e
		push!(errors, "T_initial function failed at z=0: $e")
	end

	for (j, ξ_func) in enumerate(problem.ξ_initial)
		try
			ξ_func(0.0)
		catch e
			push!(errors, "ξ_initial[$j] function failed at z=0: $e")
		end
	end

	# Test external flux function
	try
		problem.external_flux(0.0)
	catch e
		push!(errors, "external_flux function failed at t=0: $e")
	end

	# Test ambient temperature function
	try
		problem.ambient_temperature(0.0)
	catch e
		push!(errors, "ambient_temperature function failed at t=0: $e")
	end

	# Check cell volumes are positive
	for (i, vol) in enumerate(mesh.cell_volumes)
		if vol <= 0
			push!(errors, "Cell $i has non-positive volume: $vol")
		end
	end

	return errors
end

# ============================================================================
# Quick Validation Helpers
# ============================================================================

"""
	is_valid_problem(problem::PyrolysisProblem) -> Bool

Quick check if a problem is valid (doesn't throw).
"""
function is_valid_problem(problem::PyrolysisProblem)
	try
		validate_problem(problem; strict = false)
		return true
	catch
		return false
	end
end

"""
	check_consistency(problem::PyrolysisProblem)

Check consistency between problem components (mesh, material, BCs).
Returns a named tuple with status and messages.
"""
function check_consistency(problem::PyrolysisProblem)
	messages = String[]
	is_consistent = true

	# Check NC consistency
	NC_material = n_components(problem.material)
	if !isempty(problem.mesh.cell_states)
		NC_mesh = length(problem.mesh.cell_states[1].ξ)
		if NC_material != NC_mesh
			push!(messages, "Component count mismatch: material has $NC_material, mesh has $NC_mesh")
			is_consistent = false
		end
	end

	# Check ξ_initial count
	if length(problem.ξ_initial) != NC_material
		push!(messages, "ξ_initial count ($(length(problem.ξ_initial))) doesn't match component count ($NC_material)")
		is_consistent = false
	end

	return (valid = is_consistent, messages = messages)
end
