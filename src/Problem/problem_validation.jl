"""
Problem validation utilities.
"""

# ============================================================================
# Main Validation Function
# ============================================================================

"""
	validate_problem(problem::PyrolysisProblem) -> true

Validate a PyrolysisProblem for consistency, as `solve` does first. Throws
an error listing every check that failed:

1. The mesh's node positions are finite and match its nodes
2. Every reaction names components the material has
3. Boundary conditions are defined for all boundaries
4. The time span and the first step are valid
5. There is one initial concentration per component
6. Gas-bearing problems start with physical pore volume and complete boundary data

Setup additionally checks the pore-volume uncertainty margin and matrix support.
"""
function validate_problem(problem::PyrolysisProblem)
    errors = String[]
    append!(errors, validate_mesh(problem.mesh))
    append!(errors, validate_material_for_problem(problem.material))
    append!(errors, validate_boundary_conditions(problem))
    append!(errors, validate_tspan(problem.tspan, problem.dt_initial))
    append!(errors, validate_initial_condition_funcs(problem))
    # Pore reconstruction depends on valid geometry and initial fields.
    isempty(errors) && append!(errors, validate_initial_pore_state(problem))

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
        idx = rxn.reactant
        (idx < 1 || idx > NC) && push!(
            errors,
            "Reaction $i ($(rxn.name)): reactant index $idx out of range [1, $NC]",
        )
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
function _validate_boundaries!(errors, mesh, bc_set)
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
