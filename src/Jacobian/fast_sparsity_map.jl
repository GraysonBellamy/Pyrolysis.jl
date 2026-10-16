# ============================================================================
# FastSparsityMap: Dense Array-based O(1) Lookup (No Hash Table Overhead)
# ============================================================================
#
# Canonical sparse-write helper used by the `Structured` backend's
# operator scatter paths.

"""
	FastSparsityMap

Dense array-based sparsity map for O(1) sparse matrix writes without hash table overhead.

For a matrix with n_rows × n_cols, stores a dense lookup table where:
	indices[row + (col-1)*n_rows] = nzval_index (or 0 if not in pattern)

# Performance
- Memory: O(n_rows × n_cols) ≈ 2.5 MB for 560 DOFs
- Lookup: O(1) with simple arithmetic, no hash computation

This trades memory for speed. The dense lookup table grows quadratically with
DOF count (~90 MB at 3360 DOFs / 240 cells); very large problems pay this cost.

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

    @inbounds for col = 1:n_cols
        offset = (col - 1) * n_rows
        for nz_idx = J.colptr[col]:(J.colptr[col+1]-1)
            row = J.rowval[nz_idx]
            indices[row+offset] = nz_idx
        end
    end

    return FastSparsityMap(indices, n_rows)
end

"""
	fast_sparse_add!(J::SparseMatrixCSC, map::FastSparsityMap, row::Int, col::Int, val)

Add a value to sparse matrix using dense array lookup.

Equivalent to `J[row, col] += val` but O(1) with no hash overhead.
Does nothing if (row, col) is not in the sparsity pattern.
"""
@inline function fast_sparse_add!(
    J::SparseMatrixCSC,
    map::FastSparsityMap,
    row::Int,
    col::Int,
    val,
)
    idx = row + (col - 1) * map.n_rows
    nz_idx = @inbounds map.indices[idx]
    if nz_idx > 0
        @inbounds J.nzval[nz_idx] += val
    end
    return nothing
end
