"""
Combined source terms.
"""

"""
	compute_energy_sources!(Q::Vector{Float64}, mesh::Mesh1D{NC},
						   material::Material{NC}, prop_cache) where NC

Compute total energy source term for each cell [W/m³].

Q = -Σᵢ hᵢ rᵢ

Evaluates reactions per-cell using each cell's local state (T, ξ).

Sign convention:
- h > 0: endothermic (absorbs heat) → Q < 0 (cools)
- h < 0: exothermic (releases heat) → Q > 0 (heats)

# Arguments
- `Q`: Output vector of heat sources [W/m³]
- `mesh`: Mesh with current cell states
- `material`: Material definition with reactions
- `prop_cache`: Property cache (unused, for API consistency)

# Note
For production use, prefer the inline per-cell evaluation in `compute_rhs!`
which avoids allocating temporary rate/heat vectors. This function is provided
for testing and post-processing convenience.
"""
function compute_energy_sources!(
    Q::AbstractVector,
    state,
    material::Material{NC},
    prop_cache = nothing,
) where {NC}
    n_cells_val = length(state.T)
    NR = n_reactions(material)
    Tu = eltype(state.T)
    rates = zeros(Tu, NR)
    heats = zeros(Tu, NR)

    @inbounds for i = 1:n_cells_val
        T_i = state.T[i]
        ξ_i = ntuple(j -> state.ξ[j, i], Val(NC))
        evaluate_reactions!(rates, heats, ξ_i, T_i, material)
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
    S_conv::Vector{Float64}        # Gas convective volumetric source per cell

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
        zeros(n_cells),      # S_conv
        zeros(NC, n_faces),  # J_gas
    )
end

"""
	SourceTermCache(mesh::Mesh1D{NC}, material::Material{NC}) where NC

Create a SourceTermCache sized for the given mesh and material.
"""
function SourceTermCache(mesh::Mesh1D{NC}, material::Material{NC}) where {NC}
    n_cells_val = length(mesh.cells)
    n_faces_val = length(mesh.faces)
    NR = n_reactions(material)
    return SourceTermCache(n_cells_val, n_faces_val, NC, NR)
end

"""
	resize_source_cache!(cache::SourceTermCache, n_cells::Int, n_faces::Int, NC::Int, NR::Int)

Resize source cache arrays.
"""
function resize_source_cache!(
    cache::SourceTermCache,
    n_cells::Int,
    n_faces::Int,
    NC::Int,
    NR::Int,
)
    resize!(cache.rates, NR)
    resize!(cache.heats, NR)
    resize!(cache.Q_rxn, n_cells)
    cache.dξ_rxn = zeros(NC, n_cells)
    resize!(cache.Q_rad, n_cells)
    resize!(cache.q_cond, n_faces)
    resize!(cache.S_conv, n_cells)
    cache.J_gas = zeros(NC, n_faces)
    return cache
end
