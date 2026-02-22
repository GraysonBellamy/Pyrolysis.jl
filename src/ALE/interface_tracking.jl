"""
Interface tracking methods.

Track sharp interfaces (e.g., char/virgin boundary, reaction fronts)
through the mesh during pyrolysis simulations.

Reference: Design Document Section 6.5
"""

using ..Core: AbstractInterfaceTracker
using ..Geometry: Mesh1D

# ============================================================================
# Interface Tracker Type
# ============================================================================

"""
	InterfaceTracker <: AbstractInterfaceTracker

Tracks the position of interfaces defined by threshold crossing.

# Fields
- `interface_positions::Vector{Float64}`: Current interface positions [m]
- `interface_temperatures::Vector{Float64}`: Temperatures at interfaces [K]
- `tracking_variable::Symbol`: Variable to track (:T for temperature, or component name)
- `threshold::Float64`: Threshold value defining interface
- `use_level_set::Bool`: If true, use level-set method
- `φ::Vector{Float64}`: Level-set field (if used)
"""
mutable struct InterfaceTracker <: AbstractInterfaceTracker
	interface_positions::Vector{Float64}
	interface_temperatures::Vector{Float64}
	tracking_variable::Symbol
	threshold::Float64
	use_level_set::Bool
	φ::Vector{Float64}
end

"""
	InterfaceTracker(; variable=:T, threshold=0.5, use_level_set=false)

Create an interface tracker with specified parameters.

# Arguments
- `variable`: Variable to track (:T or component symbol)
- `threshold`: Value defining interface (e.g., 0.5 for 50% conversion)
- `use_level_set`: Whether to use level-set method for interface tracking
"""
function InterfaceTracker(; variable::Symbol = :T, threshold::Float64 = 0.5,
	use_level_set::Bool = false)
	InterfaceTracker(Float64[], Float64[], variable, threshold, use_level_set, Float64[])
end

# ============================================================================
# Interface Detection
# ============================================================================

"""
	find_interface_position(mesh::Mesh1D, values::Vector{Float64}, 
						   threshold::Float64) -> Vector{Float64}

Find positions where values cross the threshold.

Returns all crossing positions (there may be multiple interfaces).

# Arguments
- `mesh`: Current mesh
- `values`: Cell-centered values of the tracked variable
- `threshold`: Threshold value defining interface

# Returns
Vector of interface positions [m]
"""
function find_interface_position(mesh::Mesh1D, values::Vector{Float64},
	threshold::Float64)
	interfaces = Float64[]
	n_cells = length(mesh.cells)

	@inbounds for i in 1:(n_cells-1)
		if !mesh.active_cell_mask[i] || !mesh.active_cell_mask[i+1]
			continue
		end

		v_i = values[i]
		v_next = values[i+1]

		# Check for crossing
		if (v_i - threshold) * (v_next - threshold) < 0.0
			# Linear interpolation to find crossing point
			z_i = mesh.cell_centers[i]
			z_next = mesh.cell_centers[i+1]

			# t such that v_i + t*(v_next - v_i) = threshold
			t = (threshold - v_i) / (v_next - v_i)
			z_interface = z_i + t * (z_next - z_i)

			push!(interfaces, z_interface)
		end
	end

	return interfaces
end

"""
	update_interface!(tracker::InterfaceTracker, mesh::Mesh1D, 
					  u::Vector{Float64}, NC::Int)

Update interface positions based on current solution.

# Arguments
- `tracker`: Interface tracker to update
- `mesh`: Current mesh
- `u`: Solution vector
- `NC`: Number of components
"""
function update_interface!(tracker::InterfaceTracker, mesh::Mesh1D,
	u::Vector{Float64}, NC::Int)
	n_cells = length(mesh.cells)
	n_dof_per_cell = 1 + NC

	# Extract the tracked variable
	values = zeros(n_cells)

	if tracker.tracking_variable == :T
		# Track temperature
		@inbounds for i in 1:n_cells
			idx = (i - 1) * n_dof_per_cell + 1
			values[i] = u[idx]
		end
	else
		# Track a specific component (by index or name)
		# Assume tracking_variable is something like :ξ1, :ξ2, etc.
		var_str = string(tracker.tracking_variable)

		if startswith(var_str, "ξ") || startswith(var_str, "xi")
			# Extract component index
			comp_idx = parse(Int, filter(isdigit, var_str))
			if comp_idx >= 1 && comp_idx <= NC
				@inbounds for i in 1:n_cells
					idx = (i - 1) * n_dof_per_cell + 1 + comp_idx
					values[i] = u[idx]
				end
			end
		end
	end

	# Find interface positions
	tracker.interface_positions = find_interface_position(mesh, values, tracker.threshold)

	# Also store temperatures at interfaces
	if !isempty(tracker.interface_positions)
		tracker.interface_temperatures = zeros(length(tracker.interface_positions))

		for (k, z_int) in enumerate(tracker.interface_positions)
			# Interpolate temperature at interface
			tracker.interface_temperatures[k] = interpolate_at_position(
				mesh, u, NC, 1, z_int)  # var_idx=1 for temperature
		end
	end

	return nothing
end

"""
	update_interface_from_state!(tracker::InterfaceTracker, mesh::Mesh1D)

Update interface positions from mesh cell states directly.
"""
function update_interface_from_state!(tracker::InterfaceTracker, mesh::Mesh1D)
	n_cells = length(mesh.cells)
	NC = length(mesh.cell_states[1].ξ)
	values = zeros(n_cells)

	if tracker.tracking_variable == :T
		@inbounds for i in 1:n_cells
			values[i] = mesh.cell_states[i].T
		end
	else
		# Component tracking
		var_str = string(tracker.tracking_variable)
		if startswith(var_str, "ξ") || startswith(var_str, "xi")
			comp_idx = parse(Int, filter(isdigit, var_str))
			if comp_idx >= 1 && comp_idx <= NC
				@inbounds for i in 1:n_cells
					values[i] = mesh.cell_states[i].ξ[comp_idx]
				end
			end
		end
	end

	tracker.interface_positions = find_interface_position(mesh, values, tracker.threshold)

	# Store temperatures at interfaces
	if !isempty(tracker.interface_positions)
		tracker.interface_temperatures = zeros(length(tracker.interface_positions))
		for (k, z_int) in enumerate(tracker.interface_positions)
			tracker.interface_temperatures[k] = interpolate_temperature_at_position(mesh, z_int)
		end
	end

	return nothing
end

# ============================================================================
# Interface Velocity
# ============================================================================

"""
	interface_velocity(tracker::InterfaceTracker, 
					   positions_old::Vector{Float64}, dt::Float64) -> Vector{Float64}

Compute interface velocity from position change.

# Arguments
- `tracker`: Interface tracker with current positions
- `positions_old`: Previous interface positions
- `dt`: Time step

# Returns
Interface velocities [m/s]
"""
function interface_velocity(tracker::InterfaceTracker,
	positions_old::Vector{Float64}, dt::Float64)
	n_interfaces = min(length(tracker.interface_positions), length(positions_old))

	if n_interfaces == 0 || dt <= 0.0
		return Float64[]
	end

	velocities = zeros(n_interfaces)

	@inbounds for i in 1:n_interfaces
		velocities[i] = (tracker.interface_positions[i] - positions_old[i]) / dt
	end

	return velocities
end

"""
	interface_velocity_from_flux(mesh::Mesh1D, reaction_rates::Vector{Float64},
								  interface_position::Float64,
								  Δρ::Float64) -> Float64

Compute interface velocity from reaction rate using Stefan condition.

v_interface = -r / Δρ

where r is the reaction rate at the interface and Δρ is the density jump.

# Arguments
- `mesh`: Current mesh
- `reaction_rates`: Reaction rate at each cell [kg/(m³·s)]
- `interface_position`: Position of interface [m]
- `Δρ`: Density jump across interface [kg/m³]

# Returns
Interface velocity [m/s]
"""
function interface_velocity_from_flux(mesh::Mesh1D, reaction_rates::Vector{Float64},
	interface_position::Float64,
	Δρ::Float64)
	if abs(Δρ) < eps()
		return 0.0
	end

	# Interpolate reaction rate at interface
	r_interface = interpolate_value_at_position(mesh, reaction_rates, interface_position)

	return -r_interface / Δρ
end

# ============================================================================
# Helper Functions
# ============================================================================

"""
	interpolate_at_position(mesh, u, NC, var_idx, z) -> Float64

Interpolate a variable from solution vector at position z.
"""
function interpolate_at_position(mesh::Mesh1D, u::Vector{Float64},
	NC::Int, var_idx::Int, z::Float64)
	n_cells = length(mesh.cells)
	n_dof_per_cell = 1 + NC

	# Find bracketing cells
	i_left = 0
	i_right = 0

	@inbounds for i in 1:n_cells
		if mesh.active_cell_mask[i]
			if mesh.cell_centers[i] <= z
				i_left = i
			end
			if mesh.cell_centers[i] >= z && i_right == 0
				i_right = i
			end
		end
	end

	if i_left == 0 && i_right > 0
		# Extrapolate from right
		idx = (i_right - 1) * n_dof_per_cell + var_idx
		return u[idx]
	elseif i_right == 0 && i_left > 0
		# Extrapolate from left
		idx = (i_left - 1) * n_dof_per_cell + var_idx
		return u[idx]
	elseif i_left == i_right
		idx = (i_left - 1) * n_dof_per_cell + var_idx
		return u[idx]
	else
		# Linear interpolation
		z_left = mesh.cell_centers[i_left]
		z_right = mesh.cell_centers[i_right]

		idx_left = (i_left - 1) * n_dof_per_cell + var_idx
		idx_right = (i_right - 1) * n_dof_per_cell + var_idx

		v_left = u[idx_left]
		v_right = u[idx_right]

		t = (z - z_left) / (z_right - z_left)
		return v_left + t * (v_right - v_left)
	end
end

"""
	interpolate_temperature_at_position(mesh, z) -> Float64

Interpolate temperature from cell states at position z.
"""
function interpolate_temperature_at_position(mesh::Mesh1D, z::Float64)
	n_cells = length(mesh.cells)

	i_left = 0
	i_right = 0

	@inbounds for i in 1:n_cells
		if mesh.active_cell_mask[i]
			if mesh.cell_centers[i] <= z
				i_left = i
			end
			if mesh.cell_centers[i] >= z && i_right == 0
				i_right = i
			end
		end
	end

	if i_left == 0 && i_right > 0
		return mesh.cell_states[i_right].T
	elseif i_right == 0 && i_left > 0
		return mesh.cell_states[i_left].T
	elseif i_left == i_right
		return mesh.cell_states[i_left].T
	else
		z_left = mesh.cell_centers[i_left]
		z_right = mesh.cell_centers[i_right]
		T_left = mesh.cell_states[i_left].T
		T_right = mesh.cell_states[i_right].T

		t = (z - z_left) / (z_right - z_left)
		return T_left + t * (T_right - T_left)
	end
end

"""
	interpolate_value_at_position(mesh, values, z) -> Float64

Interpolate a cell-centered value at position z.
"""
function interpolate_value_at_position(mesh::Mesh1D, values::Vector{Float64},
	z::Float64)
	n_cells = length(mesh.cells)

	i_left = 0
	i_right = 0

	@inbounds for i in 1:n_cells
		if mesh.active_cell_mask[i]
			if mesh.cell_centers[i] <= z
				i_left = i
			end
			if mesh.cell_centers[i] >= z && i_right == 0
				i_right = i
			end
		end
	end

	if i_left == 0 && i_right > 0
		return values[i_right]
	elseif i_right == 0 && i_left > 0
		return values[i_left]
	elseif i_left == i_right
		return values[i_left]
	else
		z_left = mesh.cell_centers[i_left]
		z_right = mesh.cell_centers[i_right]
		v_left = values[i_left]
		v_right = values[i_right]

		t = (z - z_left) / (z_right - z_left)
		return v_left + t * (v_right - v_left)
	end
end

# ============================================================================
# Mesh Velocity for Interface Tracking
# ============================================================================

"""
	interface_tracking_velocity!(w::Vector{Float64}, mesh::Mesh1D,
								  tracker::InterfaceTracker, v_interface::Float64;
								  width::Float64=0.0)

Compute mesh velocity contribution for interface tracking.

Near the interface, mesh velocity is set to track the interface.
Away from interface, velocity transitions smoothly to background.

# Arguments
- `w`: Output mesh velocity to add to [m/s]
- `mesh`: Current mesh
- `tracker`: Interface tracker with current position
- `v_interface`: Target interface velocity [m/s]
- `width`: Width of transition region (0 = sharp transition)
"""
function interface_tracking_velocity!(w::Vector{Float64}, mesh::Mesh1D,
	tracker::InterfaceTracker,
	v_interface::Float64;
	width::Float64 = 0.0)
	if isempty(tracker.interface_positions)
		return nothing
	end

	# Use first interface if multiple
	z_interface = tracker.interface_positions[1]

	n_nodes = length(mesh.nodes)

	if width <= 0.0
		# Sharp tracking: set mesh velocity = interface velocity at interface location
		# This is typically blended with other velocity components
		@inbounds for i in 1:n_nodes
			z = mesh.nodes[i].z
			dist = abs(z - z_interface)

			# Simple linear blend over nearest cells
			h_local = mesh.cell_volumes[max(1, min(i, length(mesh.cells)))] / mesh.cross_section_area
			if dist < h_local
				blend = 1.0 - dist / h_local
				w[i] += blend * v_interface
			end
		end
	else
		# Smooth tracking with Gaussian blending
		@inbounds for i in 1:n_nodes
			z = mesh.nodes[i].z
			dist = abs(z - z_interface)

			blend = exp(-dist^2 / (2.0 * width^2))
			w[i] += blend * v_interface
		end
	end

	return nothing
end
