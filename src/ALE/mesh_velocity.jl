"""
Mesh velocity computation methods.

Computes the velocity field for mesh motion in ALE formulation.
Different methods trade off between:
- Lagrangian accuracy (tracking material motion)
- Mesh quality maintenance (preventing tangling)
- Interface tracking (following reaction fronts)

Reference: Design Document Section 3.6
"""

using ..Core: AbstractMeshVelocityMethod, GAS, SOLID
using ..Geometry: Mesh1D, cell_volume

# ============================================================================
# Mesh Velocity Method Types
# ============================================================================

"""
	LaplacianSmoothing <: AbstractMeshVelocityMethod

Mesh velocity determined by solving ∇²w = 0 with boundary conditions.

For 1D, this gives a linear interpolation between boundary velocities.

# Fields
- `iterations::Int`: Number of smoothing iterations (for iterative solve)
- `relaxation::Float64`: Relaxation factor ∈ (0, 1]
"""
struct LaplacianSmoothing <: AbstractMeshVelocityMethod
	iterations::Int
	relaxation::Float64
end
LaplacianSmoothing() = LaplacianSmoothing(10, 0.5)

"""
	SpringAnalogy <: AbstractMeshVelocityMethod

Mesh velocity from spring network connecting nodes.

Each edge acts as a spring with stiffness k ∝ 1/L^power where L is edge length.
Shorter edges are stiffer, preventing collapse.

# Fields
- `stiffness_power::Float64`: Power for stiffness-length relationship
"""
struct SpringAnalogy <: AbstractMeshVelocityMethod
	stiffness_power::Float64
end
SpringAnalogy() = SpringAnalogy(1.0)

"""
	LinearElasticity <: AbstractMeshVelocityMethod

Mesh velocity from linear elasticity analogy.

Treats mesh as elastic solid, solves for displacement field.

# Fields
- `E_base::Float64`: Base Young's modulus
- `ν::Float64`: Poisson's ratio
- `stiffening_near_boundary::Float64`: Factor for boundary stiffening
"""
struct LinearElasticity <: AbstractMeshVelocityMethod
	E_base::Float64
	ν::Float64
	stiffening_near_boundary::Float64
end
LinearElasticity() = LinearElasticity(1.0, 0.3, 2.0)

"""
	WeightedALE <: AbstractMeshVelocityMethod

Weighted combination of Lagrangian and smoothing velocities.

w = α_lag·w_material + α_smooth·w_smooth + α_interface·w_interface

# Fields
- `α_lagrangian::Float64`: Weight for material-following velocity
- `α_smoothing::Float64`: Weight for mesh quality maintenance
- `α_interface::Float64`: Weight for interface tracking
- `inner_method::AbstractMeshVelocityMethod`: Method for smoothing component
"""
struct WeightedALE <: AbstractMeshVelocityMethod
	α_lagrangian::Float64
	α_smoothing::Float64
	α_interface::Float64
	inner_method::AbstractMeshVelocityMethod
end
WeightedALE() = WeightedALE(0.3, 0.6, 0.1, LaplacianSmoothing())
WeightedALE(α_lag, α_smooth) = WeightedALE(α_lag, α_smooth, 0.0, LaplacianSmoothing())

# ============================================================================
# Main Interface
# ============================================================================

"""
	compute_mesh_velocity!(w::Vector{Float64}, mesh::Mesh1D, 
						   w_top::Float64, w_bottom::Float64,
						   method::AbstractMeshVelocityMethod)

Compute mesh velocity field given boundary velocities.

# Arguments
- `w`: Output vector of mesh velocities at each node [m/s]
- `mesh`: Current mesh
- `w_top`: Mesh velocity at top boundary (boundary_id=1, z=L) [m/s]
- `w_bottom`: Mesh velocity at bottom boundary (boundary_id=2, z=0) [m/s]
- `method`: Mesh velocity computation method

# Convention
- Top surface at z=L (boundary_id=1): exposed surface, may regress (w_top < 0)
- Bottom surface at z=0 (boundary_id=2): substrate, typically fixed (w_bottom = 0)
- Surface regression: as material pyrolyzes, top moves toward bottom (w_top < 0)

# Notes
The mesh velocity field must satisfy the GCL for conservation.
"""
function compute_mesh_velocity!(w::Vector{Float64}, mesh::Mesh1D,
	w_top::Float64, w_bottom::Float64,
	method::AbstractMeshVelocityMethod)
	_compute_mesh_velocity_impl!(w, mesh, w_top, w_bottom, method)
end

# ============================================================================
# Laplacian Smoothing Implementation
# ============================================================================

"""
	_compute_mesh_velocity_impl!(w, mesh, w_top, w_bottom, method::LaplacianSmoothing)

Laplacian smoothing: solve ∇²w = 0 with Dirichlet BCs.

In 1D, the solution is simply linear interpolation.
"""
function _compute_mesh_velocity_impl!(w::Vector{Float64}, mesh::Mesh1D,
	w_top::Float64, w_bottom::Float64,
	method::LaplacianSmoothing)
	n_nodes = length(mesh.nodes)

	# Get domain bounds
	z_min = Inf
	z_max = -Inf
	top_node_idx = 0
	bottom_node_idx = 0

	@inbounds for i in 1:n_nodes
		node = mesh.nodes[i]
		if node.boundary_id == 1  # Top
			z_max = node.z
			top_node_idx = i
		elseif node.boundary_id == 2  # Bottom
			z_min = node.z
			bottom_node_idx = i
		end
		z_min = min(z_min, node.z)
		z_max = max(z_max, node.z)
	end

	L = z_max - z_min

	# In 1D, Laplacian smoothing gives linear interpolation
	@inbounds for i in 1:n_nodes
		node = mesh.nodes[i]

		if node.boundary_id == 1  # Top boundary
			w[i] = w_top
		elseif node.boundary_id == 2  # Bottom boundary
			w[i] = w_bottom
		else
			# Linear interpolation
			# t = 0 at bottom, t = 1 at top
			t = (node.z - z_min) / L
			w[i] = w_bottom + t * (w_top - w_bottom)
		end
	end

	return nothing
end

# ============================================================================
# Spring Analogy Implementation
# ============================================================================

"""
	_compute_mesh_velocity_impl!(w, mesh, w_top, w_bottom, method::SpringAnalogy)

Spring analogy: iteratively solve spring network equilibrium.

In 1D, this is equivalent to weighted averaging where weights depend on 
inverse edge lengths.
"""
function _compute_mesh_velocity_impl!(w::Vector{Float64}, mesh::Mesh1D,
	w_top::Float64, w_bottom::Float64,
	method::SpringAnalogy)
	n_nodes = length(mesh.nodes)
	power = method.stiffness_power

	# Initialize with boundary values, zero elsewhere
	fill!(w, 0.0)
	@inbounds for i in 1:n_nodes
		node = mesh.nodes[i]
		if node.boundary_id == 1
			w[i] = w_top
		elseif node.boundary_id == 2
			w[i] = w_bottom
		end
	end

	# Build neighbor information and spring stiffnesses
	# In 1D, each interior node has exactly 2 neighbors
	neighbors = [Int[] for _ in 1:n_nodes]
	stiffnesses = [Float64[] for _ in 1:n_nodes]

	@inbounds for cell in mesh.cells
		n1, n2 = cell.node_ids
		edge_length = abs(mesh.nodes[n2].z - mesh.nodes[n1].z)
		k = 1.0 / (edge_length^power + eps())  # Avoid division by zero

		# Add neighbor relationships
		push!(neighbors[n1], n2)
		push!(stiffnesses[n1], k)
		push!(neighbors[n2], n1)
		push!(stiffnesses[n2], k)
	end

	# Gauss-Seidel iteration
	max_iter = 100
	tol = 1e-10

	for iter in 1:max_iter
		max_change = 0.0

		@inbounds for i in 1:n_nodes
			node = mesh.nodes[i]

			# Skip boundary nodes
			if node.is_boundary
				continue
			end

			# Weighted average of neighbors
			w_sum = 0.0
			k_sum = 0.0
			for (j, k) in zip(neighbors[i], stiffnesses[i])
				w_sum += k * w[j]
				k_sum += k
			end

			if k_sum > 0.0
				w_new = w_sum / k_sum
				max_change = max(max_change, abs(w_new - w[i]))
				w[i] = w_new
			end
		end

		if max_change < tol
			break
		end
	end

	return nothing
end

# ============================================================================
# Linear Elasticity Implementation (Simplified 1D)
# ============================================================================

"""
	_compute_mesh_velocity_impl!(w, mesh, w_top, w_bottom, method::LinearElasticity)

Linear elasticity: in 1D, reduces to solving 1D bar problem.

With stiffening near boundaries, interior nodes move more smoothly.
"""
function _compute_mesh_velocity_impl!(w::Vector{Float64}, mesh::Mesh1D,
	w_top::Float64, w_bottom::Float64,
	method::LinearElasticity)
	n_nodes = length(mesh.nodes)

	# Get domain bounds for stiffening calculation
	z_min = minimum(n.z for n in mesh.nodes)
	z_max = maximum(n.z for n in mesh.nodes)
	L = z_max - z_min

	# Build stiffness-weighted system
	# In 1D, this is similar to spring analogy with position-dependent stiffness
	neighbors = [Int[] for _ in 1:n_nodes]
	stiffnesses = [Float64[] for _ in 1:n_nodes]

	E_base = method.E_base
	stiff_factor = method.stiffening_near_boundary

	@inbounds for cell in mesh.cells
		n1, n2 = cell.node_ids
		z1, z2 = mesh.nodes[n1].z, mesh.nodes[n2].z
		edge_length = abs(z2 - z1)

		# Position-dependent stiffness (higher near boundaries)
		z_mid = 0.5 * (z1 + z2)
		dist_to_boundary = min(z_mid - z_min, z_max - z_mid)
		boundary_factor = 1.0 + (stiff_factor - 1.0) * exp(-4.0 * dist_to_boundary / L)

		k = E_base * boundary_factor / (edge_length + eps())

		push!(neighbors[n1], n2)
		push!(stiffnesses[n1], k)
		push!(neighbors[n2], n1)
		push!(stiffnesses[n2], k)
	end

	# Initialize
	fill!(w, 0.0)
	@inbounds for i in 1:n_nodes
		node = mesh.nodes[i]
		if node.boundary_id == 1
			w[i] = w_top
		elseif node.boundary_id == 2
			w[i] = w_bottom
		end
	end

	# Gauss-Seidel iteration
	max_iter = 100
	tol = 1e-10

	for iter in 1:max_iter
		max_change = 0.0

		@inbounds for i in 1:n_nodes
			node = mesh.nodes[i]
			if node.is_boundary
				continue
			end

			w_sum = 0.0
			k_sum = 0.0
			for (j, k) in zip(neighbors[i], stiffnesses[i])
				w_sum += k * w[j]
				k_sum += k
			end

			if k_sum > 0.0
				w_new = w_sum / k_sum
				max_change = max(max_change, abs(w_new - w[i]))
				w[i] = w_new
			end
		end

		if max_change < tol
			break
		end
	end

	return nothing
end

# ============================================================================
# Weighted ALE Implementation
# ============================================================================

"""
	_compute_mesh_velocity_impl!(w, mesh, w_top, w_bottom, method::WeightedALE)

Weighted combination of velocities.

First computes smoothing velocity, then blends with Lagrangian component.
"""
function _compute_mesh_velocity_impl!(w::Vector{Float64}, mesh::Mesh1D,
	w_top::Float64, w_bottom::Float64,
	method::WeightedALE)
	n_nodes = length(mesh.nodes)

	# Compute smoothed velocity using inner method
	w_smooth = zeros(n_nodes)
	_compute_mesh_velocity_impl!(w_smooth, mesh, w_top, w_bottom, method.inner_method)

	# Normalize weights
	total_weight = method.α_lagrangian + method.α_smoothing + method.α_interface
	α_smooth_norm = method.α_smoothing / total_weight
	α_lag_norm = method.α_lagrangian / total_weight

	# For now, Lagrangian velocity is approximated as linear interpolation
	# (In full implementation, this would come from material volume change rates)
	w_lag = zeros(n_nodes)
	_compute_mesh_velocity_impl!(w_lag, mesh, w_top, w_bottom, LaplacianSmoothing())

	# Blend
	@inbounds for i in 1:n_nodes
		node = mesh.nodes[i]
		if node.is_boundary
			# Use boundary values directly
			if node.boundary_id == 1
				w[i] = w_top
			elseif node.boundary_id == 2
				w[i] = w_bottom
			end
		else
			w[i] = α_smooth_norm * w_smooth[i] + α_lag_norm * w_lag[i]
		end
	end

	return nothing
end

# ============================================================================
# Higher-Level Functions
# ============================================================================

"""
	compute_mesh_velocity_from_volume_change!(w::Vector{Float64}, mesh::Mesh1D,
											  dξ_dt::Matrix{Float64}, material,
											  method::AbstractMeshVelocityMethod)

Compute mesh velocity from species rate of change (volume change).

# Arguments
- `w`: Output mesh velocity [m/s]
- `mesh`: Current mesh
- `dξ_dt`: Rate of change of concentrations [NC × n_cells] [kg/(m³·s)]
- `material`: Material definition
- `method`: Mesh velocity method

# Notes
Volume change rate: dV/dt = V·Σⱼ(1/ρⱼ)(dξⱼ/dt)
Surface velocity comes from integrated volume change.
"""
function compute_mesh_velocity_from_volume_change!(w::Vector{Float64}, mesh::Mesh1D,
	dξ_dt::Matrix{Float64}, material,
	method::AbstractMeshVelocityMethod)
	NC, n_cells = size(dξ_dt)
	A = mesh.cross_section_area

	# Compute volume change rate for each cell
	dV_dt = zeros(n_cells)
	@inbounds for i in 1:n_cells
		V_i = mesh.cell_volumes[i]
		for j in 1:NC
			comp = material.components[j]
			ρ_j = comp.density(mesh.cell_states[i].T)
			if ρ_j > 0.0
				# Include swelling factor: γ=0 means rigid matrix (no volume change)
				# γ=1 means full Lagrangian tracking (volume changes with mass)
				dV_dt[i] += V_i * comp.swelling_factor * dξ_dt[j, i] / ρ_j
			end
		end
	end

	# Integrate to get boundary velocity
	# Surface regression rate = -∫(dV/dt)/A integrated from bottom
	cumulative_dz_dt = 0.0
	@inbounds for i in n_cells:-1:1  # From bottom toward top
		cumulative_dz_dt += dV_dt[i] / A
	end

	# Top boundary moves with integrated volume change
	w_top = cumulative_dz_dt
	w_bottom = 0.0  # Fixed bottom

	# Compute interior velocities
	compute_mesh_velocity!(w, mesh, w_top, w_bottom, method)

	return nothing
end

"""
	surface_regression_rate(mesh::Mesh1D, mass_flux::Float64, ρ_surface::Float64) -> Float64

Compute surface regression rate from mass flux.

dz_s/dt = -ṁ / ρ_surface

# Arguments
- `mesh`: Current mesh
- `mass_flux`: Mass flux leaving through surface [kg/(m²·s)]
- `ρ_surface`: Density at surface [kg/m³]

# Returns
Surface velocity [m/s] (negative = regression)
"""
function surface_regression_rate(mesh::Mesh1D, mass_flux::Float64, ρ_surface::Float64)
	if ρ_surface > 0.0
		return -mass_flux / ρ_surface
	else
		return 0.0
	end
end
