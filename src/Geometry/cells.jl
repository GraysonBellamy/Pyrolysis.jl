"""
Cell type definitions for mesh elements.
"""

# ============================================================================
# 1D Cell
# ============================================================================

"""
	Cell1D <: AbstractCell

A cell (element) in a 1D mesh.

# Fields
- `id::Int`: Unique cell identifier
- `node_ids::Tuple{Int,Int}`: (left, right) node indices
- `face_ids::Tuple{Int,Int}`: (left, right) face indices
- `level::Int`: Refinement level (0 = coarsest)
- `parent::Int`: Parent cell index (0 = original/root cell)
- `children::Tuple{Int,Int}`: Child cell indices (0,0 = leaf cell)
"""
struct Cell1D <: AbstractCell
    id::Int
    node_ids::Tuple{Int,Int}     # (left, right) node indices
    face_ids::Tuple{Int,Int}     # (left, right) face indices
    level::Int                   # refinement level (0 = coarsest)
    parent::Int                  # parent cell (0 = original)
    children::Tuple{Int,Int}     # child cells (0,0 = leaf)
end

# Convenience constructor for leaf cells
Cell1D(id::Int, node_ids::Tuple{Int,Int}, face_ids::Tuple{Int,Int}; level::Int = 0) =
    Cell1D(id, node_ids, face_ids, level, 0, (0, 0))

# ============================================================================
# Cell State
# ============================================================================

"""
	CellState{NC}

State variables stored per cell.

# Type Parameters
- `NC`: Number of components

# Fields
- `T::Float64`: Temperature [K]
- `ξ::NTuple{NC,Float64}`: Mass concentrations [kg/m³]
- `volume::Float64`: Current cell volume [m³]
- `volume_ref::Float64`: Reference cell volume [m³]
"""
mutable struct CellState{NC}
    T::Float64                   # temperature [K]
    ξ::NTuple{NC,Float64}        # mass concentrations [kg/m³]
    volume::Float64              # current volume [m³]
    volume_ref::Float64          # reference volume [m³]
end

# Constructor with initial values
function CellState{NC}(T::Real, ξ::NTuple{NC,<:Real}, volume::Real) where {NC}
    CellState{NC}(Float64(T), map(Float64, ξ), Float64(volume), Float64(volume))
end

# Constructor with uniform concentration
function CellState{NC}(T::Real, ξ_uniform::Real, volume::Real) where {NC}
    ξ_tuple = ntuple(_ -> Float64(ξ_uniform), NC)
    CellState{NC}(Float64(T), ξ_tuple, Float64(volume), Float64(volume))
end
