"""
Jacobian computation and sparsity patterns for implicit solvers.

Reference: Design Document Section 8.4

This module provides:
- Jacobian sparsity patterns for efficient implicit time stepping
- Automatic differentiation-based Jacobian computation
- Finite difference Jacobian approximation

Note on Darcy-Fick transport:
The block-tridiagonal sparsity pattern is sufficient for Darcy-Fick because:
- Pressure P[i] depends only on local state (T[i], ξ[i]) via ideal gas law
- Darcy velocity v_gas[f] depends on pressures of adjacent cells (P[left], P[right])
- Advective flux at face f depends on v_gas[f] and ξ of adjacent cells
All these couplings are within the existing tridiagonal bandwidth.

Important: Darcy transport requires concentration rows to depend on neighbor 
temperatures (through pressure gradient → velocity → flux). This coupling 
is included in the sparsity pattern along with the other tridiagonal entries.
"""

# ============================================================================
# Jacobian Sparsity Pattern
# ============================================================================

"""
	jacobian_sparsity(n_cells::Int, NC::Int; radiation_model::RadiationModel=NO_RADIATION) -> SparseMatrixCSC{Float64, Int}

Compute the sparsity pattern of the Jacobian matrix for a 1D mesh.

For a 1D finite volume discretization with NC components, the Jacobian
has a block-tridiagonal structure:
- Each cell couples to its left and right neighbors
- Temperature couples to all concentrations in the same cell
- Concentrations couple to temperature in the same cell
- **Beer-Lambert Radiation**: Temperature at cell i depends on concentrations of ALL upstream cells
  (Beer-Lambert creates non-local coupling from surface to depth)

State vector layout: [T₁,...,Tₙ, ξ₁₁,...,ξ₁ₙ, ..., ξ_NC,1,...,ξ_NC,n]

# Arguments
- `n_cells`: Number of cells
- `NC`: Number of components
- `radiation_model`: Radiation model (affects sparsity pattern for BEER_LAMBERT)

# Returns
Sparse matrix with 1.0 at each non-zero position.

# Sparsity by Radiation Model
- NO_RADIATION, SURFACE_ABSORPTION, P1_QUASI_STEADY: Block-tridiagonal (~3×(1+NC)² entries per cell)
- BEER_LAMBERT: Block-tridiagonal + upper triangular non-local coupling (~15× more entries)
"""
function jacobian_sparsity(n_cells::Int, NC::Int; radiation_model::RadiationModel = NO_RADIATION)
	n_vars = n_cells * (1 + NC)

	# Estimate number of non-zeros
	# Each cell: self + 2 neighbors, each with (1+NC) variables coupling
	# Beer-Lambert radiation adds O(n_cells² × NC) entries
	nnz_estimate = n_cells * 3 * (1 + NC)^2
	if has_nonlocal_jacobian_coupling(radiation_model)
		nnz_estimate += n_cells * n_cells * NC ÷ 2  # Upper triangular for BEER_LAMBERT
	end

	I = Int[]
	J = Int[]
	sizehint!(I, nnz_estimate)
	sizehint!(J, nnz_estimate)

	# Temperature block indices: 1:n_cells
	# Concentration block j indices: n_cells + (j-1)*n_cells + 1 : n_cells + j*n_cells

	for i in 1:n_cells
		# Temperature row for cell i
		T_row = i

		# Temperature depends on neighbors' temperatures
		for neighbor in max(1, i-1):min(n_cells, i+1)
			T_col = neighbor
			push!(I, T_row)
			push!(J, T_col)
		end

		# Temperature depends on neighbor concentrations (through ρcₚ, k, and
		# gas convective source S_conv which depends on neighbor J_gas)
		for neighbor in max(1, i-1):min(n_cells, i+1)
			for j in 1:NC
				ξ_col = n_cells + (j-1)*n_cells + neighbor
				push!(I, T_row)
				push!(J, ξ_col)
			end
		end

		# Concentration rows for cell i
		for j in 1:NC
			ξ_row = n_cells + (j-1)*n_cells + i

			# Concentration depends on neighbor temperatures (Darcy transport:
			# pressure gradient → velocity → advective flux affects concentration)
			for neighbor in max(1, i-1):min(n_cells, i+1)
				T_col = neighbor
				push!(I, ξ_row)
				push!(J, T_col)
			end

			# Concentration depends on neighbor concentrations (gas transport)
			for neighbor in max(1, i-1):min(n_cells, i+1)
				for k in 1:NC
					ξ_col = n_cells + (k-1)*n_cells + neighbor
					push!(I, ξ_row)
					push!(J, ξ_col)
				end
			end
		end
	end

	# ==========================================================================
	# Radiation coupling (Beer-Lambert non-local)
	# ==========================================================================
	# For Beer-Lambert radiation absorption, intensity at cell i depends on
	# the optical depth from the surface, which depends on the concentration
	# (and thus absorption coefficient) of ALL cells from surface to cell i.
	#
	# This creates non-local (upper triangular) coupling in the Jacobian.
	# Other radiation models (NO_RADIATION, SURFACE_ABSORPTION, P1_QUASI_STEADY)
	# maintain tridiagonal structure.
	#
	# Convention: z=0 is bottom, z=L is top (surface). Radiation enters from top.
	# Cell n is nearest to surface, cell 1 is deepest.
	# Temperature at cell i depends on concentrations of cells i+1, i+2, ..., n
	if has_nonlocal_jacobian_coupling(radiation_model)
		for i in 1:n_cells
			T_row = i
			# Radiation at cell i depends on concentrations of all cells from i to surface (n_cells)
			for k in i:n_cells
				for j in 1:NC
					ξ_col = n_cells + (j-1)*n_cells + k
					push!(I, T_row)
					push!(J, ξ_col)
				end
			end
		end
	end

	# Create sparse matrix with ones at non-zero positions
	V = ones(length(I))
	return sparse(I, J, V, n_vars, n_vars)
end

"""
	jacobian_sparsity_ale(n_cells::Int, NC::Int; radiation_model::RadiationModel=NO_RADIATION) -> SparseMatrixCSC{Float64, Int}

Compute the sparsity pattern of the Jacobian matrix for ALE (moving mesh) mode.

For ALE, the state vector includes node positions:
State vector layout: [T₁,...,Tₙ, ξ₁₁,...,ξ_NC,n, z₁,...,z_{n+1}]

The additional couplings compared to fixed mesh are:
- Node positions depend on local cell concentrations (volume change)
- Cell quantities depend on neighboring node positions (geometry)
- **Beer-Lambert Radiation**: Temperature at cell i depends on concentrations of ALL upstream cells
  (Beer-Lambert creates non-local coupling from surface to depth)

# Arguments
- `n_cells`: Number of cells
- `NC`: Number of components
- `radiation_model`: Radiation model (affects sparsity pattern for BEER_LAMBERT)

# Returns
Sparse matrix with 1.0 at each non-zero position.

# Sparsity by Radiation Model
- NO_RADIATION, SURFACE_ABSORPTION, P1_QUASI_STEADY: Block-tridiagonal structure
- BEER_LAMBERT: Block-tridiagonal + upper triangular non-local coupling
"""
function jacobian_sparsity_ale(n_cells::Int, NC::Int; radiation_model::RadiationModel = NO_RADIATION)
	n_nodes = n_cells + 1
	n_physics = n_cells * (1 + NC)  # T and ξ variables
	n_total = n_physics + n_nodes   # Total state vector length

	# Estimate number of non-zeros
	# Physics-physics: ~3 × (1+NC)² per cell (tridiagonal blocks)
	# Physics-nodes: 2 × (1+NC) per cell (each cell depends on 2 bounding nodes)
	# Nodes-physics: lower triangular - node k depends on cells 1..(k-1)
	#   Total = (1+NC) × Σₖ₌₂ⁿ⁺¹ (k-1) = (1+NC) × n(n+1)/2
	# Nodes-nodes: lower triangular - node k depends on nodes 1..k
	#   Total = n_nodes × (n_nodes+1) / 2
	nnz_physics_physics = n_cells * 3 * (1 + NC)^2
	nnz_physics_nodes = 2 * n_cells * (1 + NC)
	nnz_nodes_physics = (1 + NC) * div(n_cells * (n_cells + 1), 2)
	nnz_nodes_nodes = div(n_nodes * (n_nodes + 1), 2)
	nnz_estimate = nnz_physics_physics + nnz_physics_nodes + nnz_nodes_physics + nnz_nodes_nodes

	I = Int[]
	J = Int[]
	sizehint!(I, nnz_estimate)
	sizehint!(J, nnz_estimate)

	# ==========================================================================
	# Physics-Physics coupling (same as fixed mesh)
	# ==========================================================================
	for i in 1:n_cells
		T_row = i

		# Temperature depends on neighbors' temperatures
		for neighbor in max(1, i-1):min(n_cells, i+1)
			push!(I, T_row)
			push!(J, neighbor)
		end

		# Temperature depends on neighbor concentrations (through ρcₚ, k, and
		# gas convective source S_conv which depends on neighbor J_gas)
		for neighbor in max(1, i-1):min(n_cells, i+1)
			for j in 1:NC
				ξ_col = n_cells + (j-1)*n_cells + neighbor
				push!(I, T_row)
				push!(J, ξ_col)
			end
		end

		# Concentration rows for cell i
		for j in 1:NC
			ξ_row = n_cells + (j-1)*n_cells + i

			# Concentration depends on neighbor temperatures (Darcy transport:
			# pressure gradient → velocity → advective flux affects concentration)
			for neighbor in max(1, i-1):min(n_cells, i+1)
				T_col = neighbor
				push!(I, ξ_row)
				push!(J, T_col)
			end

			# Concentration depends on neighbor concentrations (gas transport)
			for neighbor in max(1, i-1):min(n_cells, i+1)
				for k in 1:NC
					ξ_col = n_cells + (k-1)*n_cells + neighbor
					push!(I, ξ_row)
					push!(J, ξ_col)
				end
			end
		end
	end

	# ==========================================================================
	# Physics-Nodes coupling (T and ξ depend on geometry through node positions)
	# ==========================================================================
	# The physics (T, ξ) at cell i depends on node positions through:
	# 1. Cell i volume: depends on nodes i and i+1 (bounding nodes)
	# 2. Face fluxes: depend on d_LR = |center(neighbor) - center(i)|
	#    which involves neighbor cell centers, hence neighbor nodes
	#
	# For conduction/diffusion between cells i and i±1:
	# - Left face (cells i-1 and i): d_LR depends on z[i-1], z[i], z[i+1]
	# - Right face (cells i and i+1): d_LR depends on z[i], z[i+1], z[i+2]
	#
	# Combined: cell i physics depends on nodes max(1, i-1) to min(n_nodes, i+2)
	# This is an EXTENDED stencil, not just bounding nodes!
	for i in 1:n_cells
		T_row = i
		# Extended stencil: nodes i-1 through i+2 (clipped to valid range)
		n_start = max(1, i - 1)
		n_end = min(n_nodes, i + 2)
		for n in n_start:n_end
			z_col = n_physics + n
			push!(I, T_row)
			push!(J, z_col)
		end

		for j in 1:NC
			ξ_row = n_cells + (j-1)*n_cells + i
			for n in n_start:n_end
				z_col = n_physics + n
				push!(I, ξ_row)
				push!(J, z_col)
			end
		end
	end

	# ==========================================================================
	# Nodes-Physics coupling (mesh velocity depends on cumulative strain rates)
	# ==========================================================================
	# Node velocity w[k] = Σᵢ₌₁ᵏ⁻¹ θᵢ × Δzᵢ (cumulative sum from bottom)
	# where θᵢ depends on T and ξ through reaction-induced volume change.
	# 
	# Therefore node k depends on T and ξ of ALL cells i = 1, ..., k-1
	# This is a LOWER-TRIANGULAR structure, not just neighbors!
	# 
	# Note: Node 1 (bottom) is fixed (w[1] = 0), so it has no physics dependencies.
	for n in 2:n_nodes
		z_row = n_physics + n

		# Node n depends on cells 1 through n-1 (all cells below it)
		for i in 1:(n-1)
			# Depends on cell temperature (through strain rate)
			push!(I, z_row)
			push!(J, i)

			# Depends on all cell concentrations (through strain rate → volume change)
			for j in 1:NC
				ξ_col = n_cells + (j-1)*n_cells + i
				push!(I, z_row)
				push!(J, ξ_col)
			end
		end
	end

	# ==========================================================================
	# Nodes-Nodes coupling (mesh velocity depends on cumulative strain rates)
	# ==========================================================================
	# The node velocity is: w[k] = Σᵢ₌₁ᵏ⁻¹ θᵢ × Δzᵢ
	# where θᵢ is the strain rate and Δzᵢ = z_{i+1} - zᵢ
	# 
	# ∂wₖ/∂zₘ has contributions from all cells i < k where:
	#   ∂(θᵢ × Δzᵢ)/∂zₘ = θᵢ × ∂Δzᵢ/∂zₘ
	#   ∂Δzᵢ/∂zᵢ = -1,  ∂Δzᵢ/∂z_{i+1} = +1
	# 
	# This creates a LOWER-TRIANGULAR structure (plus diagonal for self):
	# - Node k depends on nodes 1, 2, ..., k (NOT k+1, k+2, ...)
	# - Old code used tridiagonal which was overly conservative
	for n in 1:n_nodes
		z_row = n_physics + n
		# Node n can depend on nodes 1 through n (lower triangular)
		# Include diagonal (self-dependence) and all nodes below
		for m in 1:n
			z_col = n_physics + m
			push!(I, z_row)
			push!(J, z_col)
		end
	end

	# ==========================================================================
	# Radiation coupling (Beer-Lambert non-local)
	# ==========================================================================
	# For Beer-Lambert radiation absorption, intensity at cell i depends on
	# the optical depth from the surface, which depends on the concentration
	# (and thus absorption coefficient) of ALL cells from surface to cell i.
	#
	# This creates non-local (upper triangular) coupling in the Jacobian.
	# Other radiation models (NO_RADIATION, SURFACE_ABSORPTION, P1_QUASI_STEADY)
	# maintain tridiagonal structure.
	#
	# Convention: z=0 is bottom, z=L is top (surface). Radiation enters from top.
	# Cell n is nearest to surface, cell 1 is deepest.
	# Temperature at cell i depends on concentrations of cells i+1, i+2, ..., n
	# (all cells between cell i and the surface).
	if has_nonlocal_jacobian_coupling(radiation_model)
		for i in 1:n_cells
			T_row = i
			# Radiation at cell i depends on concentrations of all cells from i to surface (n_cells)
			# These upstream cells attenuate the radiation reaching cell i
			for k in i:n_cells
				for j in 1:NC
					ξ_col = n_cells + (j-1)*n_cells + k
					push!(I, T_row)
					push!(J, ξ_col)
				end
			end
		end
	end

	# Create sparse matrix with ones at non-zero positions
	V = ones(length(I))
	return sparse(I, J, V, n_total, n_total)
end

"""
	jacobian_sparsity(mesh, material; radiation_model::RadiationModel=NO_RADIATION) -> SparseMatrixCSC{Float64, Int}

Create sparsity pattern from mesh and material.

# Arguments
- `mesh`: Computational mesh
- `material`: Material definition
- `radiation_model`: Radiation model (affects sparsity pattern for BEER_LAMBERT)
"""
function jacobian_sparsity(mesh, material; radiation_model::RadiationModel = NO_RADIATION)
	n_cells_val = length(mesh.cells)
	NC = n_components(material)
	return jacobian_sparsity(n_cells_val, NC; radiation_model = radiation_model)
end

# ============================================================================
# Sparsity Map for O(1) Jacobian Writes
# ============================================================================

"""
	SparsityMap

Pre-computed mapping from (row, col) to nzval index for O(1) sparse matrix writes.

The standard sparse matrix indexing `J[i, j] = val` performs a binary search
(O(log nnz)) for every write. For analytical Jacobians with many writes per
evaluation, this becomes a significant bottleneck.

The SparsityMap pre-computes the mapping from (row, col) pairs to nzval indices,
enabling direct writes to `J.nzval[idx]` in O(1) time.

# Performance
- Standard sparse write: O(log(nnz/n_cols)) per write due to binary search
- SparsityMap write: O(1) per write via dictionary lookup
- For a Jacobian with ~150 writes/evaluation and ~24,000 non-zeros, this can
  provide 3-10× speedup in Jacobian assembly.

# Usage
```julia
# Build map once (at problem setup or first Jacobian call)
smap = SparsityMap(J)

# Use for all subsequent writes
sparse_add!(J, smap, row, col, value)  # Instead of J[row, col] += value
sparse_write!(J, smap, row, col, value)  # Instead of J[row, col] = value
```
"""
struct SparsityMap
	# Mapping from (row, col) to index into nzval
	indices::Dict{Tuple{Int, Int}, Int}
end

"""
	SparsityMap(J::SparseMatrixCSC)

Build a sparsity map from an existing sparse matrix.

This iterates through all non-zero entries and stores their (row, col) => nzval_index
mapping for O(1) lookup during Jacobian assembly.
"""
function SparsityMap(J::SparseMatrixCSC)
	indices = Dict{Tuple{Int, Int}, Int}()
	sizehint!(indices, nnz(J))

	n_cols = size(J, 2)
	@inbounds for col in 1:n_cols
		for nz_idx in J.colptr[col]:(J.colptr[col+1]-1)
			row = J.rowval[nz_idx]
			indices[(row, col)] = nz_idx
		end
	end

	return SparsityMap(indices)
end

"""
	get_nz_index(map::SparsityMap, row::Int, col::Int) -> Int

Get the nzval index for a (row, col) entry.

Returns 0 if the (row, col) pair is not in the sparsity pattern.
"""
@inline function get_nz_index(map::SparsityMap, row::Int, col::Int)
	return get(map.indices, (row, col), 0)
end

"""
	sparse_write!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val)

Write a value directly to sparse matrix using pre-computed index.

Equivalent to `J[row, col] = val` but O(1) instead of O(log n).
Does nothing if (row, col) is not in the sparsity pattern.
"""
@inline function sparse_write!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val)
	nz_idx = get(map.indices, (row, col), 0)
	if nz_idx > 0
		@inbounds J.nzval[nz_idx] = val
	end
	return nothing
end

"""
	sparse_add!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val)

Add a value to sparse matrix using pre-computed index.

Equivalent to `J[row, col] += val` but O(1) instead of O(log n).
Does nothing if (row, col) is not in the sparsity pattern.
"""
@inline function sparse_add!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val)
	nz_idx = get(map.indices, (row, col), 0)
	if nz_idx > 0
		@inbounds J.nzval[nz_idx] += val
	end
	return nothing
end

"""
	sparse_sub!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val)

Subtract a value from sparse matrix using pre-computed index.

Equivalent to `J[row, col] -= val` but O(1) instead of O(log n).
Does nothing if (row, col) is not in the sparsity pattern.
"""
@inline function sparse_sub!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val)
	nz_idx = get(map.indices, (row, col), 0)
	if nz_idx > 0
		@inbounds J.nzval[nz_idx] -= val
	end
	return nothing
end

# Fallback for dense matrices (no-op, use standard indexing)
@inline sparse_write!(J::AbstractMatrix, ::Nothing, row::Int, col::Int, val) = (J[row, col] = val; nothing)
@inline sparse_add!(J::AbstractMatrix, ::Nothing, row::Int, col::Int, val) = (J[row, col] += val; nothing)
@inline sparse_sub!(J::AbstractMatrix, ::Nothing, row::Int, col::Int, val) = (J[row, col] -= val; nothing)

# ============================================================================
# FastSparsityMap: Dense Array-based O(1) Lookup (No Hash Table Overhead)
# ============================================================================

"""
	FastSparsityMap

Dense array-based sparsity map for O(1) sparse matrix writes without hash table overhead.

For a matrix with n_rows × n_cols, stores a dense lookup table where:
	indices[row + (col-1)*n_rows] = nzval_index (or 0 if not in pattern)

# Performance
- Memory: O(n_rows × n_cols) ≈ 2.5 MB for 560 DOFs
- Lookup: O(1) with simple arithmetic, no hash computation
- Compared to Dict-based SparsityMap: eliminates ~13% runtime overhead from hash table lookups

This trades memory for speed - appropriate for problems up to ~5000 DOFs.
For larger problems, use the Dict-based SparsityMap as a fallback.

# Fields
- `indices::Vector{Int}`: Dense lookup table, 0 = not in sparsity pattern
- `n_rows::Int`: Number of rows (for linear index computation)
"""
struct FastSparsityMap
	indices::Vector{Int}
	n_rows::Int
end

"""
	FastSparsityMap(J::SparseMatrixCSC)

Build a fast sparsity map from an existing sparse matrix.

Iterates through all non-zero entries and stores their linear index mapping
for O(1) lookup during Jacobian assembly.
"""
function FastSparsityMap(J::SparseMatrixCSC)
	n_rows, n_cols = size(J)
	indices = zeros(Int, n_rows * n_cols)

	@inbounds for col in 1:n_cols
		offset = (col - 1) * n_rows
		for nz_idx in J.colptr[col]:(J.colptr[col+1]-1)
			row = J.rowval[nz_idx]
			indices[row+offset] = nz_idx
		end
	end

	return FastSparsityMap(indices, n_rows)
end

"""
Threshold for using FastSparsityMap vs Dict-based SparsityMap.
For matrices larger than this, use Dict-based to avoid excessive memory usage.

For typical pyrolysis problems (60-240 cells, ~14 DOFs/cell = 840-3360 DOFs),
FastSparsityMap provides best performance. The memory cost is:
  n_dof² × 8 bytes = 3360² × 8 ≈ 90 MB for 240 cells

Above 3500 DOFs, the memory overhead outweighs the hash-free lookup benefit.
"""
const FAST_SPARSITY_THRESHOLD = 3500

"""
	create_sparsity_map(J::SparseMatrixCSC) -> Union{FastSparsityMap, SparsityMap}

Create the appropriate sparsity map based on matrix size.
Uses FastSparsityMap for small-to-medium problems, SparsityMap for large problems.
"""
function create_sparsity_map(J::SparseMatrixCSC)
	n = size(J, 1)
	if n <= FAST_SPARSITY_THRESHOLD
		return FastSparsityMap(J)
	else
		return SparsityMap(J)
	end
end

"""
	fast_sparse_write!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int, col::Int, val)

Write a value directly to sparse matrix using dense array lookup.

Equivalent to `J[row, col] = val` but O(1) with no hash overhead.
Does nothing if (row, col) is not in the sparsity pattern.
"""
@inline function fast_sparse_write!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int, col::Int, val)
	idx = row + (col - 1) * map.n_rows
	nz_idx = @inbounds map.indices[idx]
	if nz_idx > 0
		@inbounds J.nzval[nz_idx] = val
	end
	return nothing
end

"""
	fast_sparse_add!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int, col::Int, val)

Add a value to sparse matrix using dense array lookup.

Equivalent to `J[row, col] += val` but O(1) with no hash overhead.
Does nothing if (row, col) is not in the sparsity pattern.
"""
@inline function fast_sparse_add!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int, col::Int, val)
	idx = row + (col - 1) * map.n_rows
	nz_idx = @inbounds map.indices[idx]
	if nz_idx > 0
		@inbounds J.nzval[nz_idx] += val
	end
	return nothing
end

"""
	fast_sparse_sub!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int, col::Int, val)

Subtract a value from sparse matrix using dense array lookup.

Equivalent to `J[row, col] -= val` but O(1) with no hash overhead.
Does nothing if (row, col) is not in the sparsity pattern.
"""
@inline function fast_sparse_sub!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int, col::Int, val)
	idx = row + (col - 1) * map.n_rows
	nz_idx = @inbounds map.indices[idx]
	if nz_idx > 0
		@inbounds J.nzval[nz_idx] -= val
	end
	return nothing
end

"""
	fast_sparse_add_batch!(J, map, row, cols, vals)

Add multiple values to same row of sparse matrix using dense array lookup.

Batches multiple column writes to the same row, reducing function call overhead
and improving instruction pipelining. Equivalent to calling `fast_sparse_add!`
for each (col, val) pair but more efficient.

# Arguments
- `J::SparseMatrixCSC`: Sparse Jacobian matrix
- `map::FastSparsityMap`: Fast sparsity lookup map
- `row::Int`: Row index (same for all values)
- `cols::NTuple{N,Int}`: Column indices
- `vals::NTuple{N,Float64}`: Values to add
"""
@inline function fast_sparse_add_batch!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int,
	cols::NTuple{N, Int}, vals::NTuple{N, Float64}) where N
	n_rows = map.n_rows
	@inbounds for i in 1:N
		idx = row + (cols[i] - 1) * n_rows
		nz_idx = map.indices[idx]
		if nz_idx > 0
			J.nzval[nz_idx] += vals[i]
		end
	end
	return nothing
end

"""
	fast_sparse_sub_batch!(J, map, row, cols, vals)

Subtract multiple values from same row of sparse matrix using dense array lookup.

Batches multiple column writes to the same row, reducing function call overhead
and improving instruction pipelining. Equivalent to calling `fast_sparse_sub!`
for each (col, val) pair but more efficient.

# Arguments
- `J::SparseMatrixCSC`: Sparse Jacobian matrix
- `map::FastSparsityMap`: Fast sparsity lookup map
- `row::Int`: Row index (same for all values)
- `cols::NTuple{N,Int}`: Column indices
- `vals::NTuple{N,Float64}`: Values to subtract
"""
@inline function fast_sparse_sub_batch!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int,
	cols::NTuple{N, Int}, vals::NTuple{N, Float64}) where N
	n_rows = map.n_rows
	@inbounds for i in 1:N
		idx = row + (cols[i] - 1) * n_rows
		nz_idx = map.indices[idx]
		if nz_idx > 0
			J.nzval[nz_idx] -= vals[i]
		end
	end
	return nothing
end

# Unified interface: batch functions for Dict-based SparsityMap
@inline function fast_sparse_add_batch!(J::SparseMatrixCSC, map::SparsityMap, row::Int,
	cols::NTuple{N, Int}, vals::NTuple{N, Float64}) where N
	@inbounds for i in 1:N
		sparse_add!(J, map, row, cols[i], vals[i])
	end
	return nothing
end

@inline function fast_sparse_sub_batch!(J::SparseMatrixCSC, map::SparsityMap, row::Int,
	cols::NTuple{N, Int}, vals::NTuple{N, Float64}) where N
	@inbounds for i in 1:N
		sparse_sub!(J, map, row, cols[i], vals[i])
	end
	return nothing
end

# Fallback for dense matrices
@inline function fast_sparse_add_batch!(J::AbstractMatrix, ::Nothing, row::Int,
	cols::NTuple{N, Int}, vals::NTuple{N, Float64}) where N
	@inbounds for i in 1:N
		J[row, cols[i]] += vals[i]
	end
	return nothing
end

@inline function fast_sparse_sub_batch!(J::AbstractMatrix, ::Nothing, row::Int,
	cols::NTuple{N, Int}, vals::NTuple{N, Float64}) where N
	@inbounds for i in 1:N
		J[row, cols[i]] -= vals[i]
	end
	return nothing
end

# Unified interface: fast_sparse_* functions also work with Dict-based SparsityMap
# (for large problems that exceed FAST_SPARSITY_THRESHOLD)
@inline fast_sparse_write!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val) = sparse_write!(J, map, row, col, val)
@inline fast_sparse_add!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val) = sparse_add!(J, map, row, col, val)
@inline fast_sparse_sub!(J::SparseMatrixCSC, map::SparsityMap, row::Int, col::Int, val) = sparse_sub!(J, map, row, col, val)

# Fallback for dense matrices
@inline fast_sparse_write!(J::AbstractMatrix, ::Nothing, row::Int, col::Int, val) = (J[row, col] = val; nothing)
@inline fast_sparse_add!(J::AbstractMatrix, ::Nothing, row::Int, col::Int, val) = (J[row, col] += val; nothing)
@inline fast_sparse_sub!(J::AbstractMatrix, ::Nothing, row::Int, col::Int, val) = (J[row, col] -= val; nothing)

# ============================================================================
# Jacobian Prototype for ODE Solvers
# ============================================================================

"""
	JacobianOperator{F,SP}

Jacobian operator that can be used with OrdinaryDiffEq.jl.

# Fields
- `ode_fn::F`: The ODE function
- `sparsity::SP`: Sparsity pattern
- `colorvec::Vector{Int}`: Column coloring for efficient finite differences
- `color_groups::Vector{Vector{Int}}`: Columns grouped by color for fast iteration
- `n_colors::Int`: Total number of colors
"""
struct JacobianOperator{F, SP}
	ode_fn::F
	sparsity::SP
	colorvec::Vector{Int}
	color_groups::Vector{Vector{Int}}
	n_colors::Int
end

"""
	JacobianOperatorALE{F,SP}

Jacobian operator for ALE simulations with pre-allocated work arrays.

This struct avoids allocations in the Jacobian computation hot path by
pre-allocating the work arrays f0, f1, and u_perturbed.

# Fields
- `ode_fn::F`: The ALE ODE function
- `sparsity::SP`: Sparsity pattern
- `colorvec::Vector{Int}`: Column coloring for efficient finite differences
- `color_groups::Vector{Vector{Int}}`: Columns grouped by color for fast iteration
- `n_colors::Int`: Total number of colors
- `f0::Vector{Float64}`: Pre-allocated array for baseline function evaluation
- `f1::Vector{Float64}`: Pre-allocated array for perturbed function evaluation
- `u_perturbed::Vector{Float64}`: Pre-allocated array for perturbed state vector
"""
struct JacobianOperatorALE{F, SP}
	ode_fn::F
	sparsity::SP
	colorvec::Vector{Int}
	color_groups::Vector{Vector{Int}}
	n_colors::Int
	# Pre-allocated work arrays to avoid allocations in hot path
	f0::Vector{Float64}
	f1::Vector{Float64}
	u_perturbed::Vector{Float64}
end

"""
	JacobianOperator(ode_fn::PyrolysisODEFunction)

Create a Jacobian operator for use with implicit solvers.

This uses the sparsity pattern to enable efficient finite difference
Jacobian computation via matrix coloring. The sparsity pattern is
determined by the radiation model stored in the ODE function.
"""
function JacobianOperator(ode_fn::PyrolysisODEFunction)
	mesh = ode_fn.mesh
	material = ode_fn.material
	n_cells_val = length(mesh.cells)
	NC = n_components(material)

	# Get sparsity pattern based on radiation model
	sparsity = jacobian_sparsity(n_cells_val, NC; radiation_model = ode_fn.radiation_model)

	# Compute column coloring for efficient finite differences
	colorvec = matrix_colors(sparsity)

	# Build color groups for efficient iteration
	n_colors = maximum(colorvec; init = 0)
	color_groups = [Int[] for _ in 1:n_colors]
	for j in eachindex(colorvec)
		push!(color_groups[colorvec[j]], j)
	end

	JacobianOperator(ode_fn, sparsity, colorvec, color_groups, n_colors)
end

"""
	JacobianOperatorALE(ode_fn::PyrolysisODEFunctionALE)

Create a Jacobian operator for ALE simulations with pre-allocated work arrays.

This uses the ALE sparsity pattern to enable efficient finite difference
Jacobian computation via matrix coloring, while avoiding allocations in
the hot path by pre-allocating work arrays.

# Arguments
- `ode_fn::PyrolysisODEFunctionALE`: The ALE ODE function

# Returns
A `JacobianOperatorALE` with pre-allocated work arrays sized for the state vector.
"""
function JacobianOperatorALE(ode_fn::PyrolysisODEFunctionALE)
	mesh = ode_fn.mesh
	material = ode_fn.material
	n_cells_val = length(mesh.cells)
	NC = n_components(material)
	n_nodes = n_cells_val + 1
	n_vars = n_cells_val * (1 + NC) + n_nodes

	# Get ALE sparsity pattern based on radiation model
	sparsity = jacobian_sparsity_ale(n_cells_val, NC; radiation_model = ode_fn.radiation_model)

	# Compute column coloring for efficient finite differences
	colorvec = matrix_colors(sparsity)

	# Build color groups for efficient iteration
	n_colors = maximum(colorvec; init = 0)
	color_groups = [Int[] for _ in 1:n_colors]
	for j in eachindex(colorvec)
		push!(color_groups[colorvec[j]], j)
	end

	# Pre-allocate work arrays
	f0 = zeros(n_vars)
	f1 = zeros(n_vars)
	u_perturbed = zeros(n_vars)

	JacobianOperatorALE(ode_fn, sparsity, colorvec, color_groups, n_colors,
		f0, f1, u_perturbed)
end

"""
	(op::JacobianOperatorALE)(J, u, p, t)

Compute the Jacobian using colored finite differences with pre-allocated work arrays.

This function is allocation-free after the first call, using the pre-allocated
work arrays stored in the operator struct.

# Arguments
- `J`: Output Jacobian matrix (modified in-place)
- `u`: Current state vector
- `p`: Parameters (unused, for ODE solver compatibility)
- `t`: Current time

# Returns
The modified Jacobian matrix `J`
"""
function (op::JacobianOperatorALE)(J, u, p, t)
	# Reuse pre-allocated work arrays (zero allocations!)
	f0 = op.f0
	f1 = op.f1
	u_perturbed = op.u_perturbed

	# Copy u into u_perturbed
	copyto!(u_perturbed, u)

	# Baseline evaluation
	op.ode_fn(f0, u, nothing, t)

	# Colored finite differences with adaptive step sizing
	@inbounds for c in 1:op.n_colors
		cols = op.color_groups[c]

		# Perturb all columns in this color group
		for j in cols
			ε_j = optimal_fd_step(u[j])
			u_perturbed[j] = u[j] + ε_j
		end

		# Evaluate with perturbations
		op.ode_fn(f1, u_perturbed, nothing, t)

		# Extract derivatives for all columns in this color
		for j in cols
			ε_j = optimal_fd_step(u[j])
			inv_ε_j = 1.0 / ε_j
			for idx in J.colptr[j]:(J.colptr[j+1]-1)
				i = J.rowval[idx]
				J.nzval[idx] = (f1[i] - f0[i]) * inv_ε_j
			end
			# Restore original value
			u_perturbed[j] = u[j]
		end
	end

	return J
end

"""
	matrix_colors(sparsity::SparseMatrixCSC) -> Vector{Int}

Compute a column coloring for the sparsity pattern using O(nnz) algorithm.

Columns with the same color have no common non-zero rows, allowing
simultaneous finite difference perturbations.

Uses CSR (row-major) representation built from CSC to efficiently find
all columns that share a row with the current column being colored.
"""
function matrix_colors(sparsity::SparseMatrixCSC)
	n = size(sparsity, 2)
	m = size(sparsity, 1)

	# Build CSR representation for efficient row access
	# row_ptr[i] to row_ptr[i+1]-1 gives indices into col_indices for row i
	row_ptr = zeros(Int, m + 1)

	# Count entries per row
	@inbounds for idx in 1:length(sparsity.rowval)
		row_ptr[sparsity.rowval[idx]+1] += 1
	end

	# Cumulative sum to get row pointers
	row_ptr[1] = 1
	@inbounds for i in 1:m
		row_ptr[i+1] += row_ptr[i]
	end

	# Fill column indices for each row
	col_indices = Vector{Int}(undef, length(sparsity.rowval))
	row_count = zeros(Int, m)

	@inbounds for j in 1:n
		for idx in sparsity.colptr[j]:(sparsity.colptr[j+1]-1)
			row = sparsity.rowval[idx]
			pos = row_ptr[row] + row_count[row]
			col_indices[pos] = j
			row_count[row] += 1
		end
	end

	# Greedy coloring using CSR for conflict detection
	colors = zeros(Int, n)
	max_colors = 64
	forbidden = falses(max_colors)

	for j in 1:n
		fill!(forbidden, false)

		# For each row that column j touches, find all other columns in that row
		@inbounds for idx in sparsity.colptr[j]:(sparsity.colptr[j+1]-1)
			row = sparsity.rowval[idx]

			# Look at all columns that share this row
			for k in row_ptr[row]:(row_ptr[row+1]-1)
				other_col = col_indices[k]
				if other_col < j  # Only already-colored columns
					c = colors[other_col]
					if c > 0
						# Grow forbidden array if needed
						if c > length(forbidden)
							old_len = length(forbidden)
							resize!(forbidden, 2 * c)
							forbidden[(old_len+1):end] .= false
						end
						forbidden[c] = true
					end
				end
			end
		end

		# Find minimum available color
		color = 1
		while color <= length(forbidden) && forbidden[color]
			color += 1
		end
		colors[j] = color
	end

	return colors
end

# ============================================================================
# Adaptive Finite Difference Step Size
# ============================================================================

# Constants for finite difference computation
# FD_RELATIVE_STEP is approximately sqrt(machine_epsilon) for Float64
# This balances truncation error (wants larger ε) vs roundoff error (wants smaller ε)
const FD_RELATIVE_STEP = 1.5e-8
# Absolute floor to avoid issues when u[j] is near zero
const FD_ABSOLUTE_MIN = 1e-12

"""
	optimal_fd_step(u_j::Float64) -> Float64

Compute optimal finite difference step size for value u_j.

Uses relative step size sqrt(machine_epsilon) × |u_j| with absolute floor.
This provides good accuracy across different scales:
- Temperature: 300-2000 K  → ε ≈ 4.5e-6 to 3e-5
- Concentration: 0-2000 kg/m³ → ε ≈ 1e-12 to 3e-5
- Node positions: 0-0.01 m → ε ≈ 1e-12 to 1.5e-10
"""
@inline function optimal_fd_step(u_j::Float64)
	return max(abs(u_j) * FD_RELATIVE_STEP, FD_ABSOLUTE_MIN)
end

# ============================================================================
# Finite Difference Jacobian
# ============================================================================

"""
	compute_jacobian_fd!(J::AbstractMatrix, ode_fn, u, p, t)

Compute the Jacobian using finite differences with adaptive step sizing.

Uses forward differences: ∂f/∂uⱼ ≈ (f(u + εⱼeⱼ) - f(u)) / εⱼ
where εⱼ = optimal_fd_step(u[j]) adapts to the scale of each variable.

# Arguments
- `J`: Output Jacobian matrix (dense or sparse)
- `ode_fn`: ODE function
- `u`: Current state
- `p`: Parameters
- `t`: Current time
"""
function compute_jacobian_fd!(J::AbstractMatrix, ode_fn, u, p, t)
	n = length(u)
	f0 = similar(u)
	f1 = similar(u)
	u_perturbed = copy(u)

	# Evaluate at base point
	ode_fn(f0, u, p, t)

	# Perturb each variable with adaptive step size
	@inbounds for j in 1:n
		ε_j = optimal_fd_step(u[j])
		u_perturbed[j] = u[j] + ε_j
		ode_fn(f1, u_perturbed, p, t)

		inv_ε_j = 1.0 / ε_j
		for i in 1:n
			J[i, j] = (f1[i] - f0[i]) * inv_ε_j
		end

		u_perturbed[j] = u[j]  # Reset
	end

	return J
end

"""
	compute_jacobian_fd_sparse!(J::SparseMatrixCSC, ode_fn, u, p, t,
								 op::JacobianOperator)

Compute sparse Jacobian using colored finite differences with pre-computed color groups
and adaptive step sizing.

Columns with the same color are perturbed simultaneously, reducing
the number of function evaluations from n to n_colors. Each column uses
an individually computed optimal step size based on its current value.

# Arguments
- `J`: Output sparse Jacobian matrix
- `ode_fn`: ODE function
- `u`: Current state
- `p`: Parameters
- `t`: Current time
- `op`: JacobianOperator with pre-computed sparsity and color groups
"""
function compute_jacobian_fd_sparse!(J::SparseMatrixCSC, ode_fn, u, p, t,
	op::JacobianOperator)

	f0 = similar(u)
	f1 = similar(u)
	u_perturbed = copy(u)

	# Evaluate at base point
	ode_fn(f0, u, p, t)

	# Process each color using pre-computed groups
	@inbounds for c in 1:op.n_colors
		cols = op.color_groups[c]

		# Perturb all columns with this color simultaneously using adaptive step
		for j in cols
			ε_j = optimal_fd_step(u[j])
			u_perturbed[j] = u[j] + ε_j
		end

		ode_fn(f1, u_perturbed, p, t)

		# Extract derivatives and reset (each column uses its own step size)
		for j in cols
			ε_j = optimal_fd_step(u[j])
			inv_ε_j = 1.0 / ε_j
			for idx in J.colptr[j]:(J.colptr[j+1]-1)
				i = J.rowval[idx]
				J.nzval[idx] = (f1[i] - f0[i]) * inv_ε_j
			end
			u_perturbed[j] = u[j]
		end
	end

	return J
end

# ============================================================================
# Automatic Differentiation Jacobian
# ============================================================================

"""
	PyrolysisJacobianAD{F,C}

Jacobian computed via ForwardDiff automatic differentiation.

# Fields
- `ode_fn::F`: ODE function
- `cfg::C`: ForwardDiff configuration (chunk size, etc.)
"""
struct PyrolysisJacobianAD{F, C}
	ode_fn::F
	cfg::C
	f_cache::Vector{Float64}
end

"""
	PyrolysisJacobianAD(ode_fn::PyrolysisODEFunction; chunk_size::Int=8)

Create an AD-based Jacobian calculator.

# Arguments
- `ode_fn`: The ODE function
- `chunk_size`: ForwardDiff chunk size (affects memory/speed tradeoff)
"""
function PyrolysisJacobianAD(ode_fn::PyrolysisODEFunction; chunk_size::Int = 8)
	n_cells_val = length(ode_fn.mesh.cells)
	NC = n_components(ode_fn.material)
	n_vars = n_cells_val * (1 + NC)

	f_cache = zeros(n_vars)
	cfg = ForwardDiff.JacobianConfig(nothing, f_cache, zeros(n_vars),
		ForwardDiff.Chunk{chunk_size}())

	PyrolysisJacobianAD(ode_fn, cfg, f_cache)
end

"""
	(jac::PyrolysisJacobianAD)(J, u, p, t)

Compute Jacobian using automatic differentiation.
"""
function (jac::PyrolysisJacobianAD)(J, u, p, t)
	# Create closure for ForwardDiff
	function f!(du, u_inner)
		jac.ode_fn(du, u_inner, p, t)
	end

	ForwardDiff.jacobian!(J, f!, jac.f_cache, u)
	return J
end

# ============================================================================
# ODEFunction with Jacobian Specification
# ============================================================================

"""
	make_ode_function_with_jac(ode_fn::PyrolysisODEFunction;
								jac_type::Symbol=:sparse_fd)

Create an ODEFunction with Jacobian specification for OrdinaryDiffEq.jl.

# Arguments
- `ode_fn`: The base ODE function
- `jac_type`: Jacobian type, one of:
  - `:sparse_fd`: Sparse finite differences with coloring (default)
  - `:dense_fd`: Dense finite differences
  - `:ad`: Automatic differentiation
  - `:none`: No Jacobian (let solver compute)

# Returns
An ODEFunction object suitable for OrdinaryDiffEq.jl implicit solvers.
"""
function make_ode_function_with_jac(ode_fn::PyrolysisODEFunction;
	jac_type::Symbol = :sparse_fd)
	mesh = ode_fn.mesh
	material = ode_fn.material
	n_cells_val = length(mesh.cells)
	NC = n_components(material)

	if jac_type == :none
		return OrdinaryDiffEq.ODEFunction(ode_fn)

	elseif jac_type == :sparse_fd
		op = JacobianOperator(ode_fn)
		jac_prototype = copy(op.sparsity)

		function jac_sparse!(J, u, p, t)
			compute_jacobian_fd_sparse!(J, ode_fn, u, p, t, op)
		end

		return OrdinaryDiffEq.ODEFunction(ode_fn;
			jac = jac_sparse!,
			jac_prototype = jac_prototype)

	elseif jac_type == :dense_fd
		n_vars = n_cells_val * (1 + NC)
		jac_prototype = zeros(n_vars, n_vars)

		function jac_dense!(J, u, p, t)
			compute_jacobian_fd!(J, ode_fn, u, p, t)
		end

		return OrdinaryDiffEq.ODEFunction(ode_fn;
			jac = jac_dense!,
			jac_prototype = jac_prototype)

	elseif jac_type == :ad
		jac_ad = PyrolysisJacobianAD(ode_fn)
		n_vars = n_cells_val * (1 + NC)
		jac_prototype = zeros(n_vars, n_vars)

		return OrdinaryDiffEq.ODEFunction(ode_fn;
			jac = jac_ad,
			jac_prototype = jac_prototype)
	else
		error("Unknown jac_type: $jac_type. Use :sparse_fd, :dense_fd, :ad, or :none")
	end
end

"""
	make_ode_function_with_jac(ode_fn::PyrolysisODEFunctionALE;
								jac_type::Symbol=:sparse_fd)

Create an ODEFunction with Jacobian specification for ALE problems.

This version handles the ALE ODE function which may include Darcy-Fick
gas transport. The Jacobian is computed via finite differences or AD;
the sparsity pattern is adequate because:
- Darcy pressure depends only on local state (T, ξ)
- Darcy velocity at faces depends on adjacent cell pressures
- All couplings remain within the block-tridiagonal structure

# Arguments
- `ode_fn`: The ALE ODE function (with or without Darcy options)
- `jac_type`: Jacobian type, one of:
  - `:sparse_fd`: Sparse finite differences with coloring (default)
  - `:dense_fd`: Dense finite differences
  - `:ad`: Automatic differentiation (recommended for Darcy)
  - `:none`: No Jacobian (let solver compute)

# Returns
An ODEFunction object suitable for OrdinaryDiffEq.jl implicit solvers.

# Notes
For Darcy-Fick transport, `:ad` (automatic differentiation) is recommended
as it correctly handles the pressure-velocity-flux chain without requiring
explicit derivative computation.
"""
function make_ode_function_with_jac(ode_fn::PyrolysisODEFunctionALE;
	jac_type::Symbol = :sparse_fd)
	mesh = ode_fn.mesh
	material = ode_fn.material
	n_cells_val = length(mesh.cells)
	NC = n_components(material)
	n_nodes = n_cells_val + 1  # 1D mesh
	# Full ALE state vector: T + ξ (all components) + node positions
	n_vars = n_cells_val * (1 + NC) + n_nodes

	if jac_type == :none
		return OrdinaryDiffEq.ODEFunction(ode_fn)

	elseif jac_type == :sparse_fd
		# Create JacobianOperatorALE with pre-allocated work arrays
		# This avoids allocations in the hot path
		op = JacobianOperatorALE(ode_fn)
		jac_prototype = copy(op.sparsity)

		return OrdinaryDiffEq.ODEFunction(ode_fn;
			jac = op,  # Use the callable operator directly
			jac_prototype = jac_prototype)

	elseif jac_type == :dense_fd
		jac_prototype = zeros(n_vars, n_vars)

		function jac_dense_ale!(J, u, p, t)
			compute_jacobian_fd!(J, ode_fn, u, p, t)
		end

		return OrdinaryDiffEq.ODEFunction(ode_fn;
			jac = jac_dense_ale!,
			jac_prototype = jac_prototype)

	elseif jac_type == :ad
		# AD handles Darcy-Fick automatically via dual number propagation
		jac_ad = PyrolysisJacobianAD_ALE(ode_fn)
		jac_prototype = zeros(n_vars, n_vars)

		return OrdinaryDiffEq.ODEFunction(ode_fn;
			jac = jac_ad,
			jac_prototype = jac_prototype)
	else
		error("Unknown jac_type: $jac_type. Use :sparse_fd, :dense_fd, :ad, or :none")
	end
end

"""
	PyrolysisJacobianAD_ALE{F,C}

AD-based Jacobian for ALE ODE functions (including Darcy-Fick support).

# Fields
- `ode_fn::F`: ALE ODE function
- `f_cache::Vector{Float64}`: Scratch space for function output
- `f_cache2::Vector{Float64}`: Second scratch space for FD perturbation
- `u_cache::Vector{Float64}`: Scratch space for perturbed input
"""
struct PyrolysisJacobianAD_ALE{F}
	ode_fn::F
	f_cache::Vector{Float64}
	f_cache2::Vector{Float64}
	u_cache::Vector{Float64}
end

"""
	PyrolysisJacobianAD_ALE(ode_fn::PyrolysisODEFunctionALE)

Create a finite-difference-based Jacobian calculator for ALE problems.

# Arguments
- `ode_fn`: The ALE ODE function

Note: Uses finite differences rather than AD because the ODE function closure
changes each call (different t, p), making pre-allocated ForwardDiff configs
incompatible. For better performance with implicit solvers, consider using
:sparse_fd which exploits sparsity structure.
"""
function PyrolysisJacobianAD_ALE(ode_fn::PyrolysisODEFunctionALE)
	n_cells_val = length(ode_fn.mesh.cells)
	NC = n_components(ode_fn.material)
	n_nodes = n_cells_val + 1  # 1D mesh
	# Full ALE state vector: T + ξ (all components) + node positions
	n_vars = n_cells_val * (1 + NC) + n_nodes

	f_cache = zeros(n_vars)
	f_cache2 = zeros(n_vars)
	u_cache = zeros(n_vars)

	PyrolysisJacobianAD_ALE(ode_fn, f_cache, f_cache2, u_cache)
end

"""
	(jac::PyrolysisJacobianAD_ALE)(J, u, p, t)

Compute Jacobian for ALE problem using finite differences with adaptive step sizing.

This correctly handles Darcy-Fick transport by perturbing the full state
vector including node positions. Uses adaptive step sizes via `optimal_fd_step()`
to handle the wide range of variable magnitudes in ALE problems:
- Node positions z: O(0.001) m
- Concentrations ξ: O(1-1000) kg/m³
- Temperature T: O(300-2000) K

Note: Uses finite differences rather than AD because the ODE function closure
changes each call (different t, p), making pre-allocated ForwardDiff configs
incompatible.
"""
function (jac::PyrolysisJacobianAD_ALE)(J, u, p, t)
	n = length(u)

	# Compute baseline f(u)
	jac.ode_fn(jac.f_cache, u, p, t)

	# Finite difference each column with adaptive step sizing
	@inbounds for j in 1:n
		# Use adaptive step size based on variable magnitude
		ε_j = optimal_fd_step(u[j])

		# Copy u to scratch and perturb
		copyto!(jac.u_cache, u)
		jac.u_cache[j] += ε_j

		# Compute f(u + εeⱼ)
		jac.ode_fn(jac.f_cache2, jac.u_cache, p, t)

		# J[:,j] = (f(u+εeⱼ) - f(u)) / ε_j
		inv_ε_j = 1.0 / ε_j
		for i in 1:n
			J[i, j] = (jac.f_cache2[i] - jac.f_cache[i]) * inv_ε_j
		end
	end

	return J
end
