"""
Initial condition specification and application.
"""

# ============================================================================
# Apply Initial Conditions
# ============================================================================

"""
	apply_initial_conditions!(mesh, material; T=300.0, ξ=nothing)

Apply uniform initial conditions.

Mesh-level utility: writes only `mesh.cell_states` and does **not** touch
`mesh.m_dry_solid_initial` (the χ progress denominator). For setting up a
problem for the residual/Jacobian composition, use
`apply_initial_conditions!(problem)` (or [`initialize_state`](@ref)), which
populates both.

# Arguments
- `mesh`: The computational mesh
- `material`: Material definition
- `T=300.0`: Uniform initial temperature [K]
- `ξ=nothing`: Initial concentrations (scalar for all, tuple/vector per component,
			   or nothing for defaults)
"""
function apply_initial_conditions!(
    mesh::Mesh1D{NC},
    material::Material{NC};
    T::Real = 300.0,
    ξ = nothing,
) where {NC}
    n_cells_val = length(mesh.cells)

    # Process concentration argument
    if ξ === nothing
        # Default: first solid at reference density, others zero
        ξ_values = ntuple(NC) do j
            comp = material.components[j]
            if j == 1 && is_solid(comp)
                Float64(comp.density(T))
            else
                0.0
            end
        end
    elseif ξ isa Real
        # Uniform value for all components
        ξ_values = ntuple(_ -> Float64(ξ), NC)
    else
        # Tuple or vector of values
        ξ_values = ntuple(j -> Float64(ξ[j]), NC)
    end

    T_float = Float64(T)

    for i = 1:n_cells_val
        state = mesh.cell_states[i]
        state.T = T_float
        state.ξ = ξ_values
    end

    return mesh
end

# ============================================================================
# Initial Condition Helpers
# ============================================================================

"""
	step_profile(z; z_interface, value_left, value_right)

Create a step function profile for initial conditions.

# Arguments
- `z`: Position
- `z_interface`: Location of the step
- `value_left`: Value for z < z_interface
- `value_right`: Value for z ≥ z_interface

# Returns
Scalar value at position z.
"""
function step_profile(z; z_interface::Real, value_left::Real, value_right::Real)
    z < z_interface ? Float64(value_left) : Float64(value_right)
end

"""
	linear_profile(z; z_min, z_max, value_min, value_max)

Create a linear profile for initial conditions.

# Arguments
- `z`: Position
- `z_min`: Start of the linear region
- `z_max`: End of the linear region
- `value_min`: Value at z_min
- `value_max`: Value at z_max

# Returns
Linearly interpolated value at position z.
"""
function linear_profile(z; z_min::Real, z_max::Real, value_min::Real, value_max::Real)
    if z <= z_min
        return Float64(value_min)
    elseif z >= z_max
        return Float64(value_max)
    else
        t = (z - z_min) / (z_max - z_min)
        return Float64(value_min + t * (value_max - value_min))
    end
end

"""
	gaussian_profile(z; z_center, width, amplitude, baseline=0.0)

Create a Gaussian profile for initial conditions.

# Arguments
- `z`: Position
- `z_center`: Center of the Gaussian
- `width`: Standard deviation (width parameter)
- `amplitude`: Peak amplitude above baseline
- `baseline=0.0`: Baseline value

# Returns
Gaussian-distributed value at position z.
"""
function gaussian_profile(
    z;
    z_center::Real,
    width::Real,
    amplitude::Real,
    baseline::Real = 0.0,
)
    Float64(baseline + amplitude * exp(-0.5 * ((z - z_center) / width)^2))
end

"""
	tanh_profile(z; z_center, width, value_left, value_right)

Create a smooth step (tanh) profile for initial conditions.

# Arguments
- `z`: Position
- `z_center`: Center of the transition
- `width`: Transition width (smaller = sharper)
- `value_left`: Value far to the left
- `value_right`: Value far to the right

# Returns
Smoothly transitioning value at position z.
"""
function tanh_profile(z; z_center::Real, width::Real, value_left::Real, value_right::Real)
    t = 0.5 * (1.0 + tanh((z - z_center) / width))
    Float64(value_left + t * (value_right - value_left))
end

# ============================================================================
# Set Initial State from Problem
# ============================================================================

"""
	apply_initial_conditions!(problem::PyrolysisProblem)

Apply the problem's initial conditions to its mesh: evaluate the T/ξ
callables at every cell center AND populate the variable cross-section
progress denominator `mesh.m_dry_solid_initial` (the χ denominator read
under a `lateral_shrinkage_law`). This is the single problem-level IC
method — [`initialize_state`](@ref) calls it as its first step.

Returns the modified mesh.
"""
function apply_initial_conditions!(problem::PyrolysisProblem)
    mesh = problem.mesh
    NC = n_components(problem.material)
    # Cell states first: the dry-solid populator reads them.
    for i in eachindex(mesh.cells)
        z = mesh.cell_centers[i]
        state = mesh.cell_states[i]
        state.T = Float64(problem.T_initial(z))
        state.ξ = ntuple(j -> Float64(problem.ξ_initial[j](z)), NC)
    end
    populate_m_dry_solid_initial!(mesh, problem.material)
    return mesh
end

"""
	populate_m_dry_solid_initial!(mesh::Mesh1D{NC}, material::Material{NC}) where NC

Populate `mesh.m_dry_solid_initial[i]` from initial cell state.

`m_dry_solid_initial[i] = V_i(0) · Σ_{j ∈ solid} ξⱼ(0)`. This is the constant
denominator of the per-cell pyrolysis-progress variable
`χ_i = m_released_dry,i / m_dry_solid_initial,i`. Under the `WithChi` mode the
residual evolves χ (dχ_i/dt = S_dry,i·V_i / m_dry_solid_initial,i) and the
lateral-shrinkage law consumes the column mean χ̄; under `NoChi` the slot is
stored but χ stays at zero.
"""
function populate_m_dry_solid_initial!(mesh::Mesh1D{NC}, material::Material{NC}) where {NC}
    n = length(mesh.cells)
    resize!(mesh.m_dry_solid_initial, n)

    @inbounds for i = 1:n
        V_i = mesh.cell_volumes[i]
        ρ_solid = 0.0
        for j = 1:NC
            if is_solid(material.components[j])
                ρ_solid += mesh.cell_states[i].ξ[j]
            end
        end
        mesh.m_dry_solid_initial[i] = V_i * ρ_solid
    end

    return mesh
end
