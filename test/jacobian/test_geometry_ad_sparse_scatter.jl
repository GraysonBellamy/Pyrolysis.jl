"""
Dense-vs-sparse scatter equivalence for `apply_geometry_jacobian!`.

The dense method overwrites every entry in the (n_total × n_nodes) z-column
block. The sparse method walks the existing
`colptr`/`rowval` and writes only the structural entries. This test
asserts that on the production sparsity pattern (`S_phys ∪ S_geom` from
`build_sparsity`), the two scatters agree bit-for-bit on the structural
entries, and the entries the sparse method "skips" are already zero in the
dense output (i.e., the sparsity pattern is not missing any structural
contribution).

Coverage: both ALE canonical configs from `perf_helpers.jl`. Eulerian
configs are skipped — they have no z block.
"""

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis.Jacobian: GeometryADCache, apply_geometry_jacobian!, build_sparsity
const _state_length = Pyrolysis.Problem.state_length

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

function _build_empty_sparse_csc_from_pattern(S_pattern::SparseMatrixCSC{Bool,Int})
    # Same colptr/rowval as the Bool pattern, but Float64 nzval initialized
    # to zero. We write into this matrix and compare against the dense.
    return SparseMatrixCSC{Float64,Int}(
        size(S_pattern, 1),
        size(S_pattern, 2),
        copy(S_pattern.colptr),
        copy(S_pattern.rowval),
        zeros(Float64, length(S_pattern.nzval)),
    )
end

@testset "apply_geometry_jacobian! sparse vs dense scatter" begin
    for (name, builder) in (
        ("ale_no_radiation", build_ale_no_radiation),
        ("ale_beer_lambert_darcy", build_ale_beer_lambert_darcy),
    )
        @testset "$(name)" begin
            def, ws, u, _ = builder()
            cache = GeometryADCache(def, ws)
            n_total = _state_length(ws.layout)
            n_nodes = ws.layout.n_nodes
            NC = 2  # both canonical configs are 2-component
            z_off = ws.layout.n_cells * (1 + NC)
            z_range = (z_off+1):(z_off+n_nodes)

            # Dense scatter: full overwrite of the z-column block.
            J_dense = zeros(Float64, n_total, n_total)
            apply_geometry_jacobian!(J_dense, cache, u, def, 0.0)

            # Sparse scatter: production pattern (S_phys ∪ S_geom).
            S_phys, S_geom = build_sparsity(ws.layout, ws.layout.mode)
            S_union = S_phys .| S_geom
            J_sparse = _build_empty_sparse_csc_from_pattern(S_union)
            apply_geometry_jacobian!(J_sparse, cache, u, def, 0.0)

            # 1. Structural entries must agree bit-for-bit (same closure, same
            #    cache, same Float64 J_dense buffer). No ≈ tolerance: this is
            #    a memory-copy equivalence, not a numerics question.
            @testset "structural entries agree exactly" begin
                colptr = J_sparse.colptr
                rowval = J_sparse.rowval
                nzval = J_sparse.nzval
                mismatches = 0
                for col = 1:n_nodes
                    J_col = z_off + col
                    for idx = colptr[J_col]:(colptr[J_col+1]-1)
                        row = rowval[idx]
                        mismatches += (nzval[idx] != J_dense[row, J_col])
                    end
                end
                @test mismatches == 0
            end

            # 2. Entries the sparse pattern omits in the z-column block must
            #    already be zero in the dense output. If this fails, the
            #    sparsity pattern (`_emit_geometry_pattern!` in sparsity.jl)
            #    is missing a structural entry — fix the pattern, not the
            #    scatter.
            @testset "non-structural z-column entries are zero in dense" begin
                structural = falses(n_total, n_nodes)
                colptr = J_sparse.colptr
                rowval = J_sparse.rowval
                for col = 1:n_nodes
                    J_col = z_off + col
                    for idx = colptr[J_col]:(colptr[J_col+1]-1)
                        structural[rowval[idx], col] = true
                    end
                end
                max_dropped = 0.0
                for col = 1:n_nodes, row = 1:n_total
                    if !structural[row, col]
                        v = abs(J_dense[row, z_off+col])
                        v > max_dropped && (max_dropped = v)
                    end
                end
                # 1e-12 absolute is well below any meaningful Jacobian entry
                # in the canonical problems; if a dropped entry is larger,
                # it's a real missed coupling.
                @test max_dropped < 1e-12
            end

            # 3. As a single combined contract: the sparse-as-dense view of
            #    the z columns equals the dense block within bit-precision.
            @testset "z-column block equality (combined)" begin
                @test Matrix(J_sparse[:, z_range]) ≈ J_dense[:, z_range] atol=1e-12 rtol=eps(
                    Float64,
                )*16
            end
        end
    end
end
