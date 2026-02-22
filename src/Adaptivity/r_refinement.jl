"""
r-refinement (node relocation).

Reference: Design Document Section 4.5, 7.5

r-refinement moves nodes (without changing mesh topology) to better resolve
solution features. The equidistribution principle is used: distribute nodes
such that each cell has approximately equal contribution to a monitor function.

## Equidistribution Principle
For a monitor function ω(z), we want:
  ∫₀^{z_n} ω(z) dz = (n/N) ∫₀^L ω(z) dz

This concentrates cells where ω is large.
"""

using ..Geometry: Mesh1D, Node1D, cell_size, cell_center, update_geometry!, n_cells, n_nodes
using ..Geometry: domain_bounds
using ..ALE: remap_solution!, FirstOrderRemap

# ============================================================================
# Monitor Functions
# ============================================================================

"""
	MonitorFunction

Abstract type for monitor functions used in r-refinement.
"""
abstract type MonitorFunction end

"""
	GradientMonitor <: MonitorFunction

Monitor based on solution gradient: ω = √(1 + α|∇φ|²)

# Fields
- `α::Float64`: Scaling parameter for gradient sensitivity
- `variable::Symbol`: Variable to monitor (:T or component index)
"""
struct GradientMonitor <: MonitorFunction
	α::Float64
	variable::Symbol
end
GradientMonitor(; α = 1.0, variable = :T) = GradientMonitor(α, variable)

"""
	InterfaceMonitor <: MonitorFunction

Monitor that concentrates nodes near an interface.

ω = ω_base + ω_peak · exp(-(z - z_interface)² / (2σ²))

# Fields
- `z_interface::Float64`: Interface position
- `σ::Float64`: Width of concentration zone
- `ω_base::Float64`: Base monitor value
- `ω_peak::Float64`: Peak monitor value at interface
"""
struct InterfaceMonitor <: MonitorFunction
	z_interface::Float64
	σ::Float64
	ω_base::Float64
	ω_peak::Float64
end
InterfaceMonitor(z::Float64; σ = 0.001, ω_base = 1.0, ω_peak = 10.0) =
	InterfaceMonitor(z, σ, ω_base, ω_peak)

"""
	CompositeMonitor <: MonitorFunction

Combines multiple monitor functions.

# Fields  
- `monitors::Vector{MonitorFunction}`
- `combination::Symbol`: :max, :sum, or :product
"""
struct CompositeMonitor <: MonitorFunction
	monitors::Vector{MonitorFunction}
	combination::Symbol
end

# ============================================================================
# Monitor Evaluation
# ============================================================================

"""
	evaluate_monitor(mon::MonitorFunction, mesh::Mesh1D, 
					 T::Vector{Float64}, ξ::Matrix{Float64}) -> Vector{Float64}

Evaluate monitor function at cell centers.
"""
function evaluate_monitor end

function evaluate_monitor(mon::GradientMonitor, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64}) where NC
	n = n_cells(mesh)
	ω = zeros(n)

	# Select variable
	φ = mon.variable == :T ? T : ξ[1, :]

	@inbounds for i in 1:n
		if !mesh.active_cell_mask[i]
			continue
		end

		# Compute gradient
		grad = compute_cell_gradient(φ, mesh, i)
		ω[i] = sqrt(1.0 + mon.α * grad^2)
	end

	return ω
end

function evaluate_monitor(mon::InterfaceMonitor, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64}) where NC
	n = n_cells(mesh)
	ω = zeros(n)

	@inbounds for i in 1:n
		if !mesh.active_cell_mask[i]
			continue
		end

		z = mesh.cell_centers[i]
		d = z - mon.z_interface
		ω[i] = mon.ω_base + mon.ω_peak * exp(-d^2 / (2 * mon.σ^2))
	end

	return ω
end

function evaluate_monitor(mon::CompositeMonitor, mesh::Mesh1D{NC},
	T::Vector{Float64}, ξ::Matrix{Float64}) where NC
	n = n_cells(mesh)

	# Evaluate all component monitors
	component_ωs = [evaluate_monitor(sub, mesh, T, ξ) for sub in mon.monitors]

	ω = zeros(n)
	if mon.combination == :max
		@inbounds for i in 1:n
			ω[i] = maximum(comp[i] for comp in component_ωs)
		end
	elseif mon.combination == :sum
		@inbounds for i in 1:n
			ω[i] = sum(comp[i] for comp in component_ωs)
		end
	elseif mon.combination == :product
		@inbounds for i in 1:n
			ω[i] = prod(comp[i] for comp in component_ωs)
		end
	end

	return ω
end

"""
	compute_cell_gradient(φ::Vector{Float64}, mesh::Mesh1D, i::Int) -> Float64

Compute gradient of field φ at cell i.
"""
function compute_cell_gradient(φ::AbstractVector{Float64}, mesh::Mesh1D, i::Int)
	cell = mesh.cells[i]
	left_face = mesh.faces[cell.face_ids[1]]
	right_face = mesh.faces[cell.face_ids[2]]

	left_id = left_face.cell_ids[1]
	right_id = right_face.cell_ids[2]

	if left_id == 0 && right_id == 0
		return 0.0
	elseif left_id == 0
		dz = mesh.cell_centers[right_id] - mesh.cell_centers[i]
		return abs(dz) > 0 ? (φ[right_id] - φ[i]) / dz : 0.0
	elseif right_id == 0
		dz = mesh.cell_centers[i] - mesh.cell_centers[left_id]
		return abs(dz) > 0 ? (φ[i] - φ[left_id]) / dz : 0.0
	else
		dz = mesh.cell_centers[right_id] - mesh.cell_centers[left_id]
		return abs(dz) > 0 ? (φ[right_id] - φ[left_id]) / dz : 0.0
	end
end

# ============================================================================
# Equidistribution Algorithm
# ============================================================================

"""
	compute_target_positions(mesh::Mesh1D, ω::Vector{Float64}) -> Vector{Float64}

Compute target node positions using equidistribution.

# Arguments
- `mesh`: Current mesh
- `ω`: Monitor function values at cell centers

# Returns
- Vector of target z-positions for interior nodes
"""
function compute_target_positions(mesh::Mesh1D, ω::Vector{Float64})
	z_min, z_max = domain_bounds(mesh)
	n_interior_nodes = n_nodes(mesh) - 2  # Exclude boundary nodes

	if n_interior_nodes ≤ 0
		return Float64[]
	end

	# Build sorted list of active cells with their monitor values
	active_cells = [(i, mesh.cell_centers[i], ω[i], cell_size(mesh, i))
					for i in 1:n_cells(mesh) if mesh.active_cell_mask[i]]
	sort!(active_cells, by = x -> x[2])  # Sort by z position

	# Compute cumulative integral of monitor function
	# ∫ω dz ≈ Σ ω_i * h_i
	Ω_total = sum(cell[3] * cell[4] for cell in active_cells)

	if Ω_total ≤ 0
		return Float64[]
	end

	# Target: each segment should have equal Ω_total / N contribution
	n_segments = length(active_cells)
	target_per_segment = Ω_total / n_segments

	# Compute target positions by equidistribution
	# For each interior node n, find z such that ∫₀^z ω dz = (n/N) * Ω_total
	target_positions = zeros(n_interior_nodes)

	# Build cumulative distribution
	cumulative = zeros(length(active_cells) + 1)
	z_positions = zeros(length(active_cells) + 1)
	z_positions[1] = z_min

	for (idx, (_, z_c, ω_i, h_i)) in enumerate(active_cells)
		cumulative[idx+1] = cumulative[idx] + ω_i * h_i
		# Estimate cell right boundary
		if idx < length(active_cells)
			next_z = active_cells[idx+1][2]
			z_positions[idx+1] = 0.5 * (z_c + h_i / 2 + next_z - active_cells[idx+1][4] / 2)
		else
			z_positions[idx+1] = z_max
		end
	end

	# Find target positions by interpolating inverse of cumulative distribution
	for n in 1:n_interior_nodes
		target_Ω = (n / (n_interior_nodes + 1)) * Ω_total

		# Find interval containing target_Ω
		idx = searchsortedlast(cumulative, target_Ω)
		idx = clamp(idx, 1, length(cumulative) - 1)

		# Linear interpolation within interval
		if cumulative[idx+1] > cumulative[idx]
			frac = (target_Ω - cumulative[idx]) / (cumulative[idx+1] - cumulative[idx])
		else
			frac = 0.5
		end

		target_positions[n] = z_positions[idx] + frac * (z_positions[idx+1] - z_positions[idx])
	end

	return target_positions
end

"""
	relocate_nodes!(mesh::Mesh1D{NC}; monitor::MonitorFunction=GradientMonitor(),
					max_iter::Int=10, relaxation::Float64=0.5, 
					tol::Float64=1e-3, h_min::Float64=1e-8) where NC

Relocate interior nodes using equidistribution.

# Arguments
- `mesh`: Mesh to modify (in-place)
- `monitor`: Monitor function for equidistribution
- `max_iter`: Maximum iterations
- `relaxation`: Relaxation factor (0 < β ≤ 1)
- `tol`: Convergence tolerance (max relative node movement)
- `h_min`: Minimum allowed cell size

# Returns
- Number of iterations performed
"""
function relocate_nodes!(mesh::Mesh1D{NC};
	monitor::MonitorFunction = GradientMonitor(),
	max_iter::Int = 10,
	relaxation::Float64 = 0.5,
	tol::Float64 = 1e-3,
	h_min::Float64 = 1e-8) where NC

	# Extract current state for monitor evaluation
	T, ξ = extract_state_arrays(mesh)

	z_min, z_max = domain_bounds(mesh)
	L = z_max - z_min

	# Get interior node indices (not on boundaries)
	interior_node_ids = [i for i in 1:n_nodes(mesh)
							   if !mesh.nodes[i].is_boundary]

	if isempty(interior_node_ids)
		return 0
	end

	# Sort interior nodes by position
	sort!(interior_node_ids, by = i -> mesh.nodes[i].z)

	for iter in 1:max_iter
		# Evaluate monitor function
		ω = evaluate_monitor(monitor, mesh, T, ξ)

		# Compute target positions
		targets = compute_target_positions(mesh, ω)

		if length(targets) != length(interior_node_ids)
			# Mismatch - skip r-refinement
			break
		end

		# Compute node movements
		max_movement = 0.0

		for (idx, node_id) in enumerate(interior_node_ids)
			z_old = mesh.nodes[node_id].z
			z_target = targets[idx]

			# Relaxed update
			z_new = z_old + relaxation * (z_target - z_old)

			# Quality constraints: ensure minimum cell size
			# Check neighbors
			if idx > 1
				prev_node_id = interior_node_ids[idx-1]
				z_prev = mesh.nodes[prev_node_id].z
				z_new = max(z_new, z_prev + h_min)
			else
				z_new = max(z_new, z_min + h_min)
			end

			if idx < length(interior_node_ids)
				next_node_id = interior_node_ids[idx+1]
				z_next = mesh.nodes[next_node_id].z
				z_new = min(z_new, z_next - h_min)
			else
				z_new = min(z_new, z_max - h_min)
			end

			# Update node
			old_node = mesh.nodes[node_id]
			mesh.nodes[node_id] = Node1D(
				old_node.id, z_new, old_node.z_ref,
				old_node.is_boundary, old_node.boundary_id,
			)

			movement = abs(z_new - z_old) / L
			max_movement = max(max_movement, movement)
		end

		# Update geometry
		update_geometry!(mesh)

		# Check convergence
		if max_movement < tol
			return iter
		end
	end

	return max_iter
end

"""
	extract_state_arrays(mesh::Mesh1D{NC}) -> (T, ξ)

Extract temperature and concentration arrays from mesh.
"""
function extract_state_arrays(mesh::Mesh1D{NC}) where NC
	n = n_cells(mesh)
	T = zeros(n)
	ξ = zeros(NC, n)

	@inbounds for i in 1:n
		state = mesh.cell_states[i]
		T[i] = state.T
		for j in 1:NC
			ξ[j, i] = state.ξ[j]
		end
	end

	return T, ξ
end

# ============================================================================
# Combined hr-Refinement
# ============================================================================

"""
	hr_adapt!(mesh::Mesh1D{NC}; 
			  h_indicator::AbstractErrorIndicator=GradientIndicator(:T),
			  r_monitor::MonitorFunction=GradientMonitor(),
			  θ_refine::Float64=0.7, θ_coarsen::Float64=0.1,
			  h_min::Float64=1e-6, h_max::Float64=1e-2,
			  max_level::Int=4, r_iterations::Int=5) where NC

Perform combined h-refinement and r-refinement.

# Strategy (from Design Doc Section 4.6)
1. Perform h-adaptivity first (change topology)
2. Perform r-adaptivity second (smooth and optimize)
"""
function hr_adapt!(mesh::Mesh1D{NC};
	h_indicator::AbstractErrorIndicator = GradientIndicator(:T),
	r_monitor::MonitorFunction = GradientMonitor(),
	θ_refine::Float64 = 0.7,
	θ_coarsen::Float64 = 0.1,
	h_min::Float64 = 1e-6,
	h_max::Float64 = 1e-2,
	max_level::Int = 4,
	r_iterations::Int = 5) where NC

	# 1. h-adaptivity
	T, ξ = extract_state_arrays(mesh)
	η = compute_indicator(h_indicator, mesh, T, ξ)

	# Refine
	marked_refine = mark_cells_for_refinement(mesh, η, θ_refine, h_min, max_level)
	marked_refine = enforce_2to1_balance!(mesh, marked_refine, max_level)
	perform_h_refinement!(mesh, marked_refine)

	# Coarsen
	T, ξ = extract_state_arrays(mesh)  # Re-extract after refinement
	η = compute_indicator(h_indicator, mesh, T, ξ)
	marked_coarsen = mark_cells_for_coarsening(mesh, η, θ_coarsen, h_max)
	perform_h_coarsening!(mesh, marked_coarsen)

	# 2. r-adaptivity
	if r_iterations > 0
		relocate_nodes!(mesh; monitor = r_monitor, max_iter = r_iterations, h_min = h_min)
	end

	return nothing
end
