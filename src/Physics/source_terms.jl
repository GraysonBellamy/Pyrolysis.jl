"""
Combined source terms.

Reference: Design Document Section 4.5
"""

"""
	compute_energy_sources!(Q::Vector{Float64}, mesh::Mesh1D{NC},
						   rates::Vector{Float64}, heats::Vector{Float64}) where NC

Compute total energy source term for each cell [W/m³].

Q = -Σᵢ hᵢ rᵢ

Sign convention:
- h > 0: endothermic (absorbs heat) → Q < 0 (cools)
- h < 0: exothermic (releases heat) → Q > 0 (heats)

# Arguments
- `Q`: Output vector of heat sources [W/m³]
- `mesh`: Mesh (for sizing)
- `rates`: Reaction rates [kg/(m³·s)]
- `heats`: Heats of reaction [J/kg]
"""
function compute_energy_sources!(Q::Vector{Float64}, mesh::Mesh1D{NC},
	rates::Vector{Float64}, heats::Vector{Float64}) where NC
	n_cells_val = length(mesh.cells)

	@inbounds for i in 1:n_cells_val
		Q[i] = reaction_heat_source(rates, heats)
	end

	return nothing
end

"""
	SourceTermCache

Pre-allocated storage for source term computations.
"""
mutable struct SourceTermCache
	# Reaction data (per cell)
	rates::Vector{Float64}         # Reaction rates
	heats::Vector{Float64}         # Heats of reaction
	Q_rxn::Vector{Float64}         # Reaction heat source

	# Species sources (NC × n_cells)
	dξ_rxn::Matrix{Float64}        # Species source from reactions

	# Radiation
	Q_rad::Vector{Float64}         # Radiation heat source

	# Heat fluxes
	q_cond::Vector{Float64}        # Conductive fluxes at faces
	q_conv::Vector{Float64}        # Convective fluxes at faces

	# Gas fluxes (NC × n_faces)
	J_gas::Matrix{Float64}         # Gas mass fluxes
end

"""
	SourceTermCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int)

Create a SourceTermCache with pre-allocated arrays.
"""
function SourceTermCache(n_cells::Int, n_faces::Int, NC::Int, NR::Int)
	SourceTermCache(
		zeros(NR),           # rates
		zeros(NR),           # heats
		zeros(n_cells),      # Q_rxn
		zeros(NC, n_cells),  # dξ_rxn
		zeros(n_cells),      # Q_rad
		zeros(n_faces),      # q_cond
		zeros(n_faces),      # q_conv
		zeros(NC, n_faces),  # J_gas
	)
end

"""
	SourceTermCache(mesh::Mesh1D{NC}, material::Material{NC}) where NC

Create a SourceTermCache sized for the given mesh and material.
"""
function SourceTermCache(mesh::Mesh1D{NC}, material::Material{NC}) where NC
	n_cells_val = length(mesh.cells)
	n_faces_val = length(mesh.faces)
	NR = n_reactions(material)
	return SourceTermCache(n_cells_val, n_faces_val, NC, NR)
end

"""
	resize_source_cache!(cache::SourceTermCache, n_cells::Int, n_faces::Int, NC::Int, NR::Int)

Resize source cache arrays.
"""
function resize_source_cache!(cache::SourceTermCache, n_cells::Int, n_faces::Int, NC::Int, NR::Int)
	resize!(cache.rates, NR)
	resize!(cache.heats, NR)
	resize!(cache.Q_rxn, n_cells)
	cache.dξ_rxn = zeros(NC, n_cells)
	resize!(cache.Q_rad, n_cells)
	resize!(cache.q_cond, n_faces)
	resize!(cache.q_conv, n_faces)
	cache.J_gas = zeros(NC, n_faces)
	return cache
end
