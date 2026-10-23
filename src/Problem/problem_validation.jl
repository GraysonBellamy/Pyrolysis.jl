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
7. DarcyFick problems start with consistent porosity/permeability
   (see [`validate_darcy_initial_porosity`](@ref))

It also warns, without failing, when Fickian gas transport degenerates where
the gas fraction vanishes (see [`warn_degenerate_gas_transport`](@ref)).

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

    # Validate Darcy-mode problems have a physically consistent initial
    # porosity: a permeable material with zero pore space makes the ideal-gas
    # pressure inversion blow up on the first implicit trial states.
    darcy_errors = validate_darcy_initial_porosity(problem)
    append!(errors, darcy_errors)

    warn_degenerate_gas_transport(problem.material)

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

function validate_mesh(mesh::Mesh1D)
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
        push!(
            errors,
            "Mesh cell_states count ($(length(mesh.cell_states))) doesn't match cells count ($(length(mesh.cells)))",
        )
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
	warn_degenerate_gas_transport(material)

Warn when a gas-producing material with Fickian transport (no permeable
component) has a condensed component whose gas-transfer coefficient λ is zero.
The cell coefficient is the volume-fraction mix of the component values, so a
cell of that component holding no gas has no gas transport at all. The gas
field then ends in a sharp front, and a Newton iterate that drives the front
cell's gas concentration below zero closes both of its faces and flattens its
residual: the implicit solve fails repeatedly on the front. A small λ on the
condensed components (`SolidComponent(...; λ = 1e-6)`) keeps the transport
non-degenerate.
"""
function warn_degenerate_gas_transport(material)
    has_darcy_flow(material) && return nothing
    comps = material.components
    produces_gas = any(material.reactions) do rxn
        any(p -> is_gas(comps[first(p)]), rxn.products)
    end
    produces_gas || return nothing
    names = [c.name for c in comps if !is_gas(c) && iszero(c.gas_transfer_coeff(T_REF))]
    isempty(names) && return nothing
    @warn "Material $(repr(material.name)): the condensed component(s) " *
          "$(join(repr.(names), ", ")) have gas-transfer coefficient λ = 0 under " *
          "Fickian gas transport. Gas transport vanishes where the gas fraction " *
          "vanishes, the gas field forms a sharp front, and the implicit Newton " *
          "fails repeatedly on it. Give the condensed components a small λ " *
          "(for example `λ = 1e-6` m²/s)." maxlog = 1
    return nothing
end

"""
	validate_boundary_conditions(problem) -> Vector{String}

Validate boundary conditions are defined for all mesh boundaries.
"""
function validate_boundary_conditions(problem)
    errors = String[]
    bc_set = problem.bc_set
    _validate_boundaries!(errors, problem.mesh, bc_set)

    # Per-species ambient concentrations (ConvectiveMassBC with tuple ξ_inf)
    # must match the material's component count — the tuple is indexed by
    # component index at flux evaluation.
    NC = Materials.n_components(problem.material)
    for (side, bc) in pairs(bc_set.mass)
        _validate_xi_inf_length!(errors, bc, NC, side)
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

function _validate_boundaries!(errors, mesh::Mesh1D, bc_set)
    # Get all boundary IDs from mesh
    mesh_boundary_ids = Set(keys(mesh.boundary_faces))

    # Check that we have BCs for the expected boundaries (1=top, 2=bottom)
    # BoundaryConditionSet always has thermal/mass BCs for both :top and :bottom,
    # so we just verify the mesh has the expected boundaries
    expected_boundaries = Set([1, 2])  # Top and bottom
    missing_boundaries = setdiff(expected_boundaries, mesh_boundary_ids)
    if !isempty(missing_boundaries)
        push!(
            errors,
            "Mesh missing expected boundaries: $missing_boundaries (expected 1=top, 2=bottom)",
        )
    end

    return errors
end

_validate_xi_inf_length!(errors, bc, NC, side) = nothing
function _validate_xi_inf_length!(errors, bc::ConvectiveMassBC, NC, side)
    if bc.ξ_inf isa Tuple && length(bc.ξ_inf) != NC
        push!(
            errors,
            "ConvectiveMassBC at :$side has a per-species ξ_inf tuple of length " *
            "$(length(bc.ξ_inf)) but the material has $NC components — the tuple " *
            "is indexed by component index and must have length $NC",
        )
    end
    return nothing
end
function _validate_xi_inf_length!(errors, bc::CombinedMassBC, NC, side)
    for sub in bc.bcs
        _validate_xi_inf_length!(errors, sub, NC, side)
    end
    return nothing
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

# Minimum initial porosity accepted for a cell whose effective permeability
# is nonzero in a DarcyFick-mode problem. The pore-pressure formula
# `ideal_gas_pressure` divides by `max(φ, 1e-12)`, so φ at (or within
# floating-point dust of) zero amplifies any generated gas pressure by up to
# ×10¹², producing 1e15–1e30 Darcy fluxes that diverge the implicit solver's
# Newton iterations. The threshold is deliberately tiny: it rejects the
# "forgot ρ_intrinsic" misconfiguration (φ exactly 0) without judging
# legitimately low but physical porosities.
const _DARCY_MIN_INITIAL_POROSITY = 1e-6

"""
	validate_darcy_initial_porosity(problem) -> Vector{String}

For DarcyFick-mode problems (material with `has_darcy_flow`), check that no
cell starts with (near-)zero porosity while its effective permeability is
nonzero. That state is physically inconsistent — a permeable matrix (κ > 0)
with no pore space for gas to occupy — and numerically fatal: the pore
pressure `P = Σⱼ ξⱼRT/Mⱼ / max(φ, 1e-12)` is amplified ×10¹² the moment any
gas appears, so the implicit solver's Newton iterations diverge (typically
surfacing as a `DomainError` from a negative trial temperature deep in the
Darcy property update).

The usual cause is a solid component constructed without `ρ_intrinsic` (or
`φ`): the intrinsic density then falls back to the bulk density and the
solid fills the cell exactly, giving φ₀ = 0.

Initial porosity and permeability are evaluated from `T_initial`/`ξ_initial`
at every cell center — the same points `apply_initial_conditions!` uses.
The porosity is the unclamped pore fraction `1 − skeletal_volume_fraction`.
Cells whose effective permeability is not positive are not flagged (zero
mobility keeps the pressure field inert there). Returns no errors for Fickian
materials. Bails out silently when the IC functions are not evaluable or
their count mismatches the material — those problems are reported by the
dedicated checks.
"""
function validate_darcy_initial_porosity(problem)
    errors = String[]

    material = problem.material
    has_darcy_flow(material) || return errors

    NC = n_components(material)
    length(problem.ξ_initial) == NC || return errors

    centers = cell_centers(problem.mesh)
    isempty(centers) && return errors

    n_bad = 0
    first_bad = 0
    φ_bad = 0.0
    κ_bad = 0.0
    try
        for (i, z) in enumerate(centers)
            T = Float64(problem.T_initial(z))
            ξ = ntuple(j -> Float64(problem.ξ_initial[j](z)), NC)
            φ = 1 - Properties.skeletal_volume_fraction(ξ, T, material)
            φ >= _DARCY_MIN_INITIAL_POROSITY && continue
            κ = Properties.effective_permeability(ξ, T, material, φ)
            κ > 0.0 || continue
            n_bad += 1
            if n_bad == 1
                first_bad, φ_bad, κ_bad = i, φ, κ
            end
        end
    catch
        # IC functions not evaluable at cell centers — reported by
        # `validate_strict`; nothing meaningful to check here.
        return errors
    end

    if n_bad > 0
        push!(
            errors,
            "DarcyFick transport with (near-)zero initial porosity: φ₀ < " *
            "$(_DARCY_MIN_INITIAL_POROSITY) in $n_bad of $(length(centers)) cells " *
            "(first: cell $first_bad, φ₀=$(φ_bad), κ_eff=$(κ_bad) m²) while the " *
            "effective permeability is nonzero. A permeable matrix with no pore " *
            "space is physically inconsistent, and the pore-pressure floor " *
            "P = Σⱼ ξⱼRT/Mⱼ / max(φ, 1e-12) amplifies any generated gas by up to " *
            "×10¹², diverging the implicit solver. Give solid/liquid components " *
            "their skeletal density via `ρ_intrinsic` (or a `φ` kwarg) so φ₀ is " *
            "the real pore fraction — e.g. wood: ρ = 500 bulk, ρ_intrinsic = 1500 " *
            "→ φ₀ ≈ 0.67 — and consider initializing the pore gas at ambient " *
            "pressure with `background_gas_concentration`. To bypass this check, " *
            "pass `validate = false` to `solve`.",
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
                "Reaction $i ($(rxn.name)) violates mass balance: |residual|=$err",
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
