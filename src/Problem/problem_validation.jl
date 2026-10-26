"""
Problem validation utilities.
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
7. All gas-bearing problems start with physical pore volume and complete boundary data

Setup additionally checks the pore-volume uncertainty margin and matrix support.

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

    # Pore reconstruction depends on valid geometry and initial fields.
    # Check every gas-bearing material once those prerequisites are valid.
    isempty(errors) && append!(errors, validate_initial_pore_state(problem))


    # Additional strict checks
    if strict && isempty(errors)
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
function validate_mesh(mesh::Meshes.Mesh)
    errors = String[]
    topo = mesh.topology
    if length(mesh.positions) != Meshes.n_nodes(topo)
        push!(
            errors,
            "Mesh has $(length(mesh.positions)) node positions for $(Meshes.n_nodes(topo)) nodes",
        )
    end
    all(p -> all(isfinite, p), mesh.positions) ||
        push!(errors, "Mesh node positions are not all finite")
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
        length(rxn.reactants) == 1 || push!(
            errors,
            "Reaction $i ($(rxn.name)) requires exactly one reactant; " *
            "additional reactant mass stoichiometries are not supported",
        )
        for (idx, _order) in rxn.reactants
            (idx < 1 || idx > NC) && push!(
                errors,
                "Reaction $i ($(rxn.name)): reactant index $idx out of range [1, $NC]",
            )
        end
        for (idx, _ν) in rxn.products
            (idx < 1 || idx > NC) && push!(
                errors,
                "Reaction $i ($(rxn.name)): product index $idx out of range [1, $NC]",
            )
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
    bc_set = problem.bc_set
    _validate_boundaries!(errors, problem.mesh, bc_set)

    for bc in values(bc_set.mass)
        try
            validate_gas_boundary(bc, problem.material, problem.tspan[1])
        catch err
            push!(errors, sprint(showerror, err))
        end
    end

    return errors
end

# The mesh's boundary tags: `:top` (the exposed surface) and `:bottom` among
# them, and a thermal condition for each.
function _validate_boundaries!(errors, mesh::Meshes.Mesh, bc_set)
    tags = Meshes.boundary_tags(mesh.topology)
    missing = setdiff([:top, :bottom], tags)
    isempty(missing) || push!(
        errors,
        "Mesh missing expected boundary tags: $(missing) (expected :top, the exposed " *
        "surface, and :bottom)",
    )
    for tag in tags
        haskey(bc_set.thermal, tag) ||
            push!(errors, "No thermal condition for the mesh's boundary tag :$tag")
    end
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
        push!(
            errors,
            "Initial timestep ($dt_initial) larger than time span ($(t_end - t_start))",
        )
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
        push!(
            errors,
            "Number of initial concentration functions ($(length(problem.ξ_initial))) doesn't match component count ($NC)",
        )
    end

    return errors
end

"""Check actual initialized gas-bearing cells, independently of mobility."""
function validate_initial_pore_state(problem)
    errors=String[]
    material=problem.material
    has_gas(material) || return errors
    NC=n_components(material)
    length(problem.ξ_initial)==NC || return errors
    for (i, z) in enumerate(cell_centers(problem.mesh))
        T=Float64(problem.T_initial(z))
        ξ=ntuple(j->Float64(problem.ξ_initial[j](z)), NC)
        φ=1-Properties.skeletal_volume_fraction(ξ, T, material)
        if !(φ>=2Physics.POROSITY_FLOOR)
            push!(
                errors,
                "cell $i: gas-bearing material needs positive initial physical pore volume in the guard identity region (raw φ=$φ)",
            )
        end
        domain=Properties.transport_domain(ξ, T, material)
        domain.admissible || push!(
            errors,
            "cell $i: $(domain.outcome), matrix occupancy $(domain.q), skeletal occupancy $(domain.skeletal); supply applicable endpoint data or stay in the declared matrix domain",
        )
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
        if !verify_stoichiometry(rxn; tol = 1e-10)
            err = mass_balance_error(rxn)
            push!(
                errors,
                isfinite(err) ?
                "Reaction $i ($(rxn.name)) violates mass balance: |residual|=$err" :
                "Reaction $i ($(rxn.name)) has unverifiable mass stoichiometry",
            )
        end
    end

    # Test initial condition functions are callable at the first cell centre
    x_test = first(cell_centers(mesh))
    try
        problem.T_initial(x_test)
    catch e
        push!(errors, "T_initial function failed at $x_test: $e")
    end

    for (j, ξ_func) in enumerate(problem.ξ_initial)
        try
            ξ_func(x_test)
        catch e
            push!(errors, "ξ_initial[$j] function failed at $x_test: $e")
        end
    end

    # Check cell volumes are positive
    for (i, vol) in enumerate(cell_volumes(mesh))
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
