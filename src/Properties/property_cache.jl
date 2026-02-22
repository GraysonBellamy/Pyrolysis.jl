"""
Property caching for performance.

Pre-allocated storage for evaluated properties to avoid allocations in inner loops.
"""

"""
	PropertyCache

Pre-allocated storage for evaluated material properties.

# Fields
## Per-cell properties
- `ρ::Vector{Float64}`: Mixture density [kg/m³]
- `ρcp::Vector{Float64}`: Volumetric heat capacity [J/(m³·K)]
- `cp::Vector{Float64}`: Specific heat capacity [J/(kg·K)]
- `k::Vector{Float64}`: Thermal conductivity [W/(m·K)]
- `λ::Vector{Float64}`: Gas transfer coefficient [m²/s]
- `ε::Vector{Float64}`: Emissivity [-]
- `α::Vector{Float64}`: Absorption coefficient [m⁻¹]

## Per-face properties
- `k_face::Vector{Float64}`: Face-centered conductivity [W/(m·K)]
- `λρ_face::Vector{Float64}`: Face-centered gas transport coefficient λ·ρ [kg/(m·s)]

## Darcy flow properties (optional, for pressure-driven gas transport)
- `P::Vector{Float64}`: Pressure field [Pa]
- `φ::Vector{Float64}`: Porosity [-]
- `κ::Vector{Float64}`: Permeability [m²]
- `μ::Vector{Float64}`: Gas viscosity [Pa·s]
- `v_gas::Vector{Float64}`: Darcy velocity at faces [m/s]
- `∇P::Vector{Float64}`: Pressure gradient at faces [Pa/m]
"""
mutable struct PropertyCache
	# Per-cell properties
	ρ::Vector{Float64}
	ρcp::Vector{Float64}
	cp::Vector{Float64}
	k::Vector{Float64}
	λ::Vector{Float64}
	ε::Vector{Float64}
	α::Vector{Float64}

	# Per-face properties (for flux computation)
	k_face::Vector{Float64}
	λρ_face::Vector{Float64}

	# Darcy flow properties (for pressure-driven gas transport)
	P::Vector{Float64}           # Pressure field [Pa]
	φ::Vector{Float64}           # Porosity [-]
	κ::Vector{Float64}           # Permeability [m²]
	μ::Vector{Float64}           # Gas viscosity [Pa·s]
	v_gas::Vector{Float64}       # Darcy velocity at faces [m/s]
	∇P::Vector{Float64}          # Pressure gradient at faces [Pa/m]
end

"""
	PropertyCache(n_cells::Int, n_faces::Int)

Create a PropertyCache with pre-allocated arrays.
"""
function PropertyCache(n_cells::Int, n_faces::Int)
	PropertyCache(
		zeros(n_cells),  # ρ
		zeros(n_cells),  # ρcp
		zeros(n_cells),  # cp
		zeros(n_cells),  # k
		zeros(n_cells),  # λ
		zeros(n_cells),  # ε
		zeros(n_cells),  # α
		zeros(n_faces),  # k_face
		zeros(n_faces),  # λρ_face
		# Darcy flow arrays
		zeros(n_cells),  # P
		zeros(n_cells),  # φ
		zeros(n_cells),  # κ
		zeros(n_cells),  # μ
		zeros(n_faces),  # v_gas
		zeros(n_faces),  # ∇P
	)
end

"""
	PropertyCache(mesh::Mesh1D)

Create a PropertyCache sized for the given mesh.
"""
PropertyCache(mesh) = PropertyCache(length(mesh.cells), length(mesh.faces))

"""
	resize_cache!(cache::PropertyCache, n_cells::Int, n_faces::Int)

Resize cache arrays (e.g., after mesh refinement).
"""
function resize_cache!(cache::PropertyCache, n_cells::Int, n_faces::Int)
	resize!(cache.ρ, n_cells)
	resize!(cache.ρcp, n_cells)
	resize!(cache.cp, n_cells)
	resize!(cache.k, n_cells)
	resize!(cache.λ, n_cells)
	resize!(cache.ε, n_cells)
	resize!(cache.α, n_cells)
	resize!(cache.k_face, n_faces)
	resize!(cache.λρ_face, n_faces)
	# Darcy flow arrays
	resize!(cache.P, n_cells)
	resize!(cache.φ, n_cells)
	resize!(cache.κ, n_cells)
	resize!(cache.μ, n_cells)
	resize!(cache.v_gas, n_faces)
	resize!(cache.∇P, n_faces)
	return cache
end

"""
	update_properties!(cache::PropertyCache, mesh::Mesh1D{NC}, material::Material{NC}) where NC

Evaluate all material properties at the current state and store in cache.

This function is designed to be allocation-free after the cache is created.
The `NC` type parameter enables the compiler to specialize for the specific
number of components, eliminating type instability and allocations.
"""
function update_properties!(cache::PropertyCache, mesh::Mesh1D{NC}, material::Material{NC}) where NC
	n_cells_val = length(mesh.cells)
	n_faces_val = length(mesh.faces)

	# Ensure cache is properly sized
	if length(cache.ρ) != n_cells_val
		resize_cache!(cache, n_cells_val, n_faces_val)
	end

	# Evaluate cell-centered properties
	@inbounds for i in 1:n_cells_val
		state = mesh.cell_states[i]
		T = state.T
		ξ = state.ξ

		cache.ρ[i] = mixture_density(ξ, T, material)
		cache.ρcp[i] = effective_heat_capacity(ξ, T, material)
		cache.cp[i] = effective_specific_heat(ξ, T, material)
		cache.k[i] = effective_conductivity(ξ, T, material)
		cache.λ[i] = effective_gas_transfer(ξ, T, material)
		cache.ε[i] = effective_emissivity(ξ, T, material)
		cache.α[i] = effective_absorption(ξ, material)
	end

	# Compute face-centered properties using distance-weighted harmonic mean
	# - k_face: harmonic mean of conductivity (for heat flux)
	# - λρ_face: harmonic mean of λ·ρ product (for gas flux)
	#
	# The distance-weighted harmonic mean gives: φ_face = (d_L + d_R) / (d_L/φ_L + d_R/φ_R)
	# This is the correct resistance-in-series interpolation for non-uniform meshes.
	# The λρ product formulation ensures physically correct flux at interfaces
	# with sharp property jumps (e.g., solid/gas pore boundaries).
	@inbounds for f in 1:n_faces_val
		face = mesh.faces[f]
		left_id, right_id = face.cell_ids
		face_pos = mesh.face_positions[f]

		if left_id == 0
			# Left boundary: use right cell value
			cache.k_face[f] = cache.k[right_id]
			cache.λρ_face[f] = cache.λ[right_id] * cache.ρ[right_id]
		elseif right_id == 0
			# Right boundary: use left cell value
			cache.k_face[f] = cache.k[left_id]
			cache.λρ_face[f] = cache.λ[left_id] * cache.ρ[left_id]
		else
			# Interior face: distance-weighted harmonic mean for resistance-in-series
			# Compute distances from cell centers to face
			d_L = abs(face_pos - mesh.cell_centers[left_id])
			d_R = abs(mesh.cell_centers[right_id] - face_pos)

			cache.k_face[f] = harmonic_mean_weighted(
				cache.k[left_id], cache.k[right_id], d_L, d_R)
			# Distance-weighted harmonic mean of λ·ρ product for gas transport
			cache.λρ_face[f] = harmonic_mean_product_weighted(
				cache.λ[left_id], cache.ρ[left_id],
				cache.λ[right_id], cache.ρ[right_id],
				d_L, d_R,
			)
		end
	end

	return cache
end

"""
	get_cell_properties(cache::PropertyCache, i::Int)

Get all properties for cell i as a named tuple.
"""
@inline function get_cell_properties(cache::PropertyCache, i::Int)
	(ρ = cache.ρ[i], ρcp = cache.ρcp[i], cp = cache.cp[i],
		k = cache.k[i], λ = cache.λ[i], ε = cache.ε[i], α = cache.α[i])
end
