"""
The initial-condition profiles and the initial fields at the cell centres.
"""

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
# The initial fields at the cell centres
# ============================================================================

"""
	initial_fields(problem, mesh::Meshes.Mesh) -> (T0, ξ0)

The problem's initial temperature [K] and concentrations [kg/m³] at the cell
centres of `mesh` (its reference positions): `T0[i] = T_initial(xᵢ)` and
`ξ0[j, i] = ξ_initial[j](xᵢ)`. On a line `xᵢ` is the scalar midpoint of the
cell; on more dimensions the centroid, an `SVector`.
"""
function initial_fields(problem::PyrolysisProblem, mesh::Meshes.Mesh)
    points = cell_centers(mesh)
    NC = length(problem.ξ_initial)
    T0 = [Float64(problem.T_initial(x)) for x in points]
    ξ0 = [Float64(problem.ξ_initial[j](x)) for j = 1:NC, x in points]
    return T0, ξ0
end

# The reference geometry of a mesh.
function _geometry(mesh::Meshes.Mesh)
    topo = mesh.topology
    g = Meshes.MeshGeometry{Float64}(topo)
    Meshes.geometry!(g, reduce(vcat, mesh.positions), topo.cross_section, topo)
    return g
end

"""
	cell_centers(mesh) -> Vector

The points at which the initial conditions are evaluated: the cell centres at
the reference positions, scalars on a line and `SVector`s on more dimensions.
"""
cell_centers(mesh::Meshes.Mesh{1}) = [c[1] for c in _geometry(mesh).cell_center]
cell_centers(mesh::Meshes.Mesh) = _geometry(mesh).cell_center

"""
	cell_volumes(mesh) -> Vector{Float64}

The cell volumes at the reference positions [m³], signed.
"""
cell_volumes(mesh::Meshes.Mesh) = _geometry(mesh).cell_volume
