# ============================================================================
# Dense automatic differentiation: the reference Jacobian
# ============================================================================

"""
	dense_jacobian(def, ws, u, t) -> Matrix{Float64}

The residual's Jacobian at `(u, t)` by dense forward-mode automatic
differentiation: the reference every structured Jacobian is checked against.
"""
function dense_jacobian(def::CoreDefinition, ws::CoreWorkspace, u::AbstractVector, t)
    f!(du, x) = residual!(du, x, def, t, ws)
    du = zeros(length(u))
    return ForwardDiff.jacobian(f!, du, Vector{Float64}(u))
end

"""
	missing_entries(J_ref, pattern) -> Vector{Tuple{Int,Int}}

The structurally nonzero entries of the reference Jacobian `J_ref` (entries
with a nonzero value) that `pattern` does not contain. Empty when the pattern
contains the reference.
"""
function missing_entries(J_ref::AbstractMatrix, pattern::AbstractMatrix{Bool})
    size(J_ref) == size(pattern) || throw(DimensionMismatch("J_ref and pattern"))
    out = Tuple{Int,Int}[]
    for j in axes(J_ref, 2), i in axes(J_ref, 1)
        J_ref[i, j] != 0 && !pattern[i, j] && push!(out, (i, j))
    end
    return out
end

"""
	compare_entries(J, J_ref; rtol = 1e-12, atol = 0.0)
	    -> (passed, max_abs_dev, worst_ratio, worst_index)

Per-entry comparison with a row-magnitude floor: every entry must satisfy
`|J[i,j] − J_ref[i,j]| ≤ atol + rtol · max(|J_ref[i,j]|, maxⱼ |J_ref[i,:]|)`.
Rows of a pyrolysis Jacobian span many orders of magnitude, so a single global
normalization would let a wholly missing small-magnitude block pass; with the
row floor a missing entry cannot hide in its own row.
"""
function compare_entries(
    J::AbstractMatrix,
    J_ref::AbstractMatrix;
    rtol::Float64 = 1e-12,
    atol::Float64 = 0.0,
)
    size(J) == size(J_ref) || throw(DimensionMismatch("J and J_ref"))
    max_dev = 0.0
    worst_ratio = 0.0
    worst = (1, 1)
    for i in axes(J_ref, 1)
        row_max = maximum(abs, view(J_ref, i, :); init = 0.0)
        for j in axes(J_ref, 2)
            dev = abs(J[i, j] - J_ref[i, j])
            max_dev = max(max_dev, dev)
            bound = atol + rtol * max(abs(J_ref[i, j]), row_max)
            ratio = bound > 0 ? dev / bound : (dev > 0 ? Inf : 0.0)
            if ratio > worst_ratio
                worst_ratio = ratio
                worst = (i, j)
            end
        end
    end
    return (
        passed = worst_ratio <= 1,
        max_abs_dev = max_dev,
        worst_ratio,
        worst_index = worst,
    )
end
