"""
Cell type definitions for mesh elements.

Reference: Design Document Section 6.3
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
	node_ids::Tuple{Int, Int}     # (left, right) node indices
	face_ids::Tuple{Int, Int}     # (left, right) face indices
	level::Int                   # refinement level (0 = coarsest)
	parent::Int                  # parent cell (0 = original)
	children::Tuple{Int, Int}     # child cells (0,0 = leaf)
end

# Convenience constructor for leaf cells
Cell1D(id::Int, node_ids::Tuple{Int, Int}, face_ids::Tuple{Int, Int}; level::Int = 0) =
	Cell1D(id, node_ids, face_ids, level, 0, (0, 0))

# ============================================================================
# 2D Cell (Quadrilateral)
# ============================================================================

"""
	Cell2D <: AbstractCell

A quadrilateral cell in a 2D mesh.

# Fields
- `id::Int`: Unique cell identifier
- `node_ids::NTuple{4,Int}`: Node indices (counter-clockwise from bottom-left)
- `face_ids::NTuple{4,Int}`: Face indices (bottom, right, top, left)
- `level::Int`: Refinement level (0 = coarsest)
- `parent::Int`: Parent cell index (0 = original/root cell)
- `children::NTuple{4,Int}`: Child cell indices (all 0 = leaf cell)
"""
struct Cell2D <: AbstractCell
	id::Int
	node_ids::NTuple{4, Int}      # counter-clockwise from bottom-left
	face_ids::NTuple{4, Int}      # (bottom, right, top, left)
	level::Int
	parent::Int
	children::NTuple{4, Int}
end

# Convenience constructor for leaf cells
Cell2D(id::Int, node_ids::NTuple{4, Int}, face_ids::NTuple{4, Int}; level::Int = 0) =
	Cell2D(id, node_ids, face_ids, level, 0, (0, 0, 0, 0))

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
	ξ::NTuple{NC, Float64}        # mass concentrations [kg/m³]
	volume::Float64              # current volume [m³]
	volume_ref::Float64          # reference volume [m³]
end

# Constructor with initial values
function CellState{NC}(T::Real, ξ::NTuple{NC, <:Real}, volume::Real) where NC
	CellState{NC}(Float64(T), map(Float64, ξ), Float64(volume), Float64(volume))
end

# Constructor with uniform concentration
function CellState{NC}(T::Real, ξ_uniform::Real, volume::Real) where NC
	ξ_tuple = ntuple(_ -> Float64(ξ_uniform), NC)
	CellState{NC}(Float64(T), ξ_tuple, Float64(volume), Float64(volume))
end

# ============================================================================
# Cell Accessors
# ============================================================================

"""
	is_leaf(cell::Cell1D)

Check if cell is a leaf (has no children).
"""
@inline is_leaf(cell::Cell1D) = cell.children == (0, 0)

"""
	is_leaf(cell::Cell2D)

Check if cell is a leaf (has no children).
"""
@inline is_leaf(cell::Cell2D) = cell.children == (0, 0, 0, 0)

"""
	is_root(cell::Union{Cell1D,Cell2D})

Check if cell is a root (has no parent).
"""
@inline is_root(cell::Union{Cell1D, Cell2D}) = cell.parent == 0

"""
	left_node(cell::Cell1D)

Get left node index of a 1D cell.
"""
@inline left_node(cell::Cell1D) = cell.node_ids[1]

"""
	right_node(cell::Cell1D)

Get right node index of a 1D cell.
"""
@inline right_node(cell::Cell1D) = cell.node_ids[2]

"""
	left_face(cell::Cell1D)

Get left face index of a 1D cell.
"""
@inline left_face(cell::Cell1D) = cell.face_ids[1]

"""
	right_face(cell::Cell1D)

Get right face index of a 1D cell.
"""
@inline right_face(cell::Cell1D) = cell.face_ids[2]

# ============================================================================
# Cell State Accessors
# ============================================================================

"""
	total_density(state::CellState)

Compute total density from component concentrations: ρ = Σξⱼ
"""
@inline function total_density(state::CellState{NC}) where NC
	ρ = 0.0
	@inbounds for j in 1:NC
		ρ += state.ξ[j]
	end
	return ρ
end

"""
	mass_fractions(state::CellState)

Compute mass fractions: Yⱼ = ξⱼ / Σξⱼ

Type-stable implementation using Val(NC).
"""
@inline function mass_fractions(state::CellState{NC}) where NC
	ρ = total_density(state)
	if ρ > 0
		inv_ρ = 1.0 / ρ
		return ntuple(Val(NC)) do j
			state.ξ[j] * inv_ρ
		end
	else
		return ntuple(_ -> 0.0, Val(NC))
	end
end

"""
	total_mass(state::CellState)

Compute total mass in cell: m = V * Σξⱼ
"""
@inline total_mass(state::CellState) = state.volume * total_density(state)
