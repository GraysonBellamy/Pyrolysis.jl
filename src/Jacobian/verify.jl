"""
Differential-correctness verification.

Compares the active backend's Jacobian against a reference
(`DenseAD_ForwardDiff` by default) at one or more sample points. The
purpose is to catch backend drift mid-development without committing the
team to a particular backend in production: if any backend's Jacobian
diverges from the AD reference beyond tolerance, this surfaces it before
the linear solver does.

This file ships the standalone comparison function `verify_jacobian_at`;
the `solve(...; verify_jacobian = ...)` kwarg that installs a
`DiscreteCallback` lives in `Solver/solve.jl`. The CI matrix in
`test/perf/test_verify_jacobian.jl` drives this function over
`(canonical_config × backend)`.

# Example

```julia
def, ws, u, _ = build_ale_no_radiation()
spec_active = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)
report = verify_jacobian_at(spec_active, def, ws, u, 0.0;
    reference = DenseAD_ForwardDiff(), rtol = 1e-6, atol = 1e-10)

report.passed              # ::Bool
report.max_abs_dev         # ::Float64
report.max_rel_dev         # ::Float64
report.worst_index         # ::Tuple{Int,Int}  – (row, col) of largest dev
```
"""

using SparseArrays: SparseMatrixCSC

"""
	verify_jacobian_at(spec_active::JacobianSpec, def, ws, u, t;
	        reference::JacobianBackend = DenseAD_ForwardDiff(),
	        rtol::Float64 = 1e-6,
	        atol::Float64 = 1e-10) -> NamedTuple

Compute both `spec_active`'s Jacobian and the `reference`-backend Jacobian
at state `(u, t)` and report the maximum absolute and relative deviations.

# Comparison modes (chosen by the [`is_approximate`](@ref) trait)

- **Strict / support-complete** (default — `is_approximate(backend) ==
  false`). Every entry must satisfy the **per-entry** tolerance

  ```
  |ΔJ_ij| ≤ atol + rtol · max(|J_ref[i,j]|, maxⱼ|J_ref[i,:]|)
  ```

  The row-magnitude floor `maxⱼ|J_ref[i,:]|` keeps roundoff at true zeros
  of a stiff row from failing the check, while a wholly missing block in
  a *small-magnitude* row (e.g. the ALE z-rows, ~10 orders below the
  Arrhenius ξ-rows) still fails — the missing entry contributes to its
  own row's scale, so it cannot hide behind the global matrix maximum.
  Entries outside the backend's declared support are compared against an
  implicit zero, so a missing sparsity slot with a non-zero reference
  entry fails at that `(row, col)`. This is the regime every "exact"
  backend must satisfy. The verifier additionally exposes
  `max_abs_dev_outside_support` so a missing-entry failure is
  distinguishable from an in-support disagreement.
- **Masked** (`is_approximate(backend) == true`). Entries outside the
  active backend's support are zeroed in the diff, and the pass criterion
  is the coarse global-scale rule
  `max_abs_dev ≤ atol || max_abs_dev / max|J_ref| ≤ rtol`. Approximate
  backends *deliberately* drop coupling contributions, so their in-support
  values legitimately deviate by the dropped magnitude; the masked compare
  is a sanity check, not an exactness gate. A diagonal slot the profile stores
  but never writes (`structural_zero_diagonal`) counts as outside the
  support. Used by
  `Structured(strategy = ApproxSolve())` and any future backend tagged
  approximate.

The reference backend is built fresh each call (lazy AD workspace
allocation on first hit, reused across calls within the same Julia
session via `JacobianSpec` construction-time caching). For repeated
verification, build the reference spec once and pass it via the
`reference_spec` keyword (low-level form, see below).

# Returned NamedTuple
- `passed::Bool` — strict: every entry within its per-entry tolerance;
  masked: global-scale criterion (see above).
- `max_abs_dev::Float64` — largest absolute element-wise difference under
  the mode's comparison rule (strict: full diff, masked: in-support diff).
- `max_rel_dev::Float64` — strict: `|ΔJ| / max(|J_ref[i,j]|, row_max_i,
  1e-12)` at the worst (most tolerance-violating) entry; masked:
  `max_abs_dev / max(maximum(abs, J_ref), 1e-12)`.
- `worst_index::Tuple{Int,Int}` — `(row, col)` of the worst entry (strict:
  largest violation ratio; masked: largest absolute deviation).
- `max_abs_dev_in_support::Float64` — in-support disagreement.
- `max_abs_dev_outside_support::Float64` — magnitude of reference entries
  the backend's support omits. For exact backends this should be within
  the per-entry tolerance; for approximate backends it is reported but
  not asserted.
- `mode::Symbol` — `:strict` or `:masked`.
- `J_active::AbstractMatrix` — the active backend's Jacobian (for
  inspection; do not retain across calls — backends may reuse buffers).
- `J_ref::Matrix{Float64}` — the dense reference Jacobian.
"""
function verify_jacobian_at(
    spec_active::JacobianSpec,
    def,
    ws,
    u::AbstractVector,
    t::Real;
    reference::JacobianBackend = DenseAD_ForwardDiff(),
    rtol::Real = 1e-6,
    atol::Real = 1e-10,
    reference_spec::Union{Nothing,JacobianSpec} = nothing,
)
    n = length(u)

    J_ref = zeros(Float64, n, n)
    ref_spec =
        reference_spec === nothing ? JacobianSpec(reference, def, ws) : reference_spec
    compute_jacobian!(J_ref, ref_spec, u, def, Float64(t), ws)

    J_active = _alloc_target(spec_active, n)
    compute_jacobian!(J_active, spec_active, u, def, Float64(t), ws)

    J_active_dense = J_active isa SparseMatrixCSC ? Matrix(J_active) : J_active
    mask = if spec_active.sparsity_phys === nothing
        ones(Bool, n, n)
    else
        Matrix(spec_active.sparsity_phys .| spec_active.sparsity_geom)
    end
    approximate = is_approximate(spec_active.backend)
    # A diagonal slot the approximate profile stores but never writes holds a
    # dropped coupling, like every entry outside its support.
    if approximate
        for r in _structural_zero_rows(spec_active.backend, def, ws)
            mask[r, r] = false
        end
    end

    # In-support disagreement — backend wrote a value that disagrees with
    # the reference at an entry inside its declared support.
    in_diff = abs.(J_active_dense .- J_ref) .* mask
    max_abs_dev_in_support = maximum(in_diff)

    # Outside-support reference magnitude — entries the backend's support
    # omits. For exact backends, this must be ≤ atol; for approximate
    # backends it is informational.
    out_diff = abs.(J_ref) .* (.!mask)
    max_abs_dev_outside_support = maximum(out_diff)

    mode = approximate ? :masked : :strict

    if approximate
        # Masked compare — only in-support disagreement is asserted, under
        # the coarse global-scale criterion (approximate backends drop
        # couplings by design; per-entry exactness is not their contract).
        max_abs_dev = max_abs_dev_in_support
        worst = argmax(in_diff)
        worst_idx = (worst[1], worst[2])
        scale = max(maximum(abs, J_ref), 1e-12)
        max_rel_dev = max_abs_dev / scale
        passed = (max_abs_dev <= atol) || (max_rel_dev <= rtol)
    else
        # Strict / support-complete compare — per-entry tolerance over the
        # full diff (in-support disagreement and missing entries alike).
        passed, max_abs_dev, max_rel_dev, worst_idx =
            _strict_per_entry_compare(J_active_dense, J_ref, Float64(rtol), Float64(atol))
    end

    return (
        passed = passed,
        max_abs_dev = max_abs_dev,
        max_rel_dev = max_rel_dev,
        worst_index = worst_idx,
        max_abs_dev_in_support = max_abs_dev_in_support,
        max_abs_dev_outside_support = max_abs_dev_outside_support,
        mode = mode,
        J_active = J_active,
        J_ref = J_ref,
    )
end

"""
	_strict_per_entry_compare(J_active, J_ref, rtol, atol)
	    -> (passed, max_abs_dev, max_rel_dev, worst_index)

Per-entry strict comparison: every entry must satisfy

```
|J_active[i,j] − J_ref[i,j]| ≤ atol + rtol · max(|J_ref[i,j]|, row_max_i)
```

with `row_max_i = maximum(abs, J_ref[i, :])`. Pyrolysis Jacobians span
~10 orders of magnitude across rows (Arrhenius ξ-rows ~1e6 vs ALE z-rows
~1e-2), so a single global normalization lets a wholly missing
small-magnitude block pass — pass/fail then depends on problem size, not
correctness. Normalizing per entry with a
*row*-magnitude floor fixes that: a missing entry contributes to its own
row's reference scale and cannot hide, while honest roundoff at true
zeros of a stiff row stays within `rtol · row_max`.

`worst_index` is the entry with the largest violation ratio
`|ΔJ| / (atol + rtol·scale)`; `max_rel_dev` is `|ΔJ| / max(scale, 1e-12)`
at that entry.
"""
function _strict_per_entry_compare(
    J_active::AbstractMatrix{Float64},
    J_ref::AbstractMatrix{Float64},
    rtol::Float64,
    atol::Float64,
)
    size(J_active) == size(J_ref) ||
        throw(DimensionMismatch("_strict_per_entry_compare: size mismatch"))
    n_rows, n_cols = size(J_ref)

    max_abs_dev = 0.0
    worst_ratio = 0.0
    worst_idx = (1, 1)
    worst_scale = 1e-12
    @inbounds for i = 1:n_rows
        row_max = 0.0
        for j = 1:n_cols
            a = abs(J_ref[i, j])
            a > row_max && (row_max = a)
        end
        for j = 1:n_cols
            dev = abs(J_active[i, j] - J_ref[i, j])
            dev > max_abs_dev && (max_abs_dev = dev)
            scale = max(abs(J_ref[i, j]), row_max)
            ratio = dev / (atol + rtol * scale)
            if ratio > worst_ratio
                worst_ratio = ratio
                worst_idx = (i, j)
                worst_scale = scale
            end
        end
    end

    passed = worst_ratio <= 1.0
    max_rel_dev =
        abs(J_active[worst_idx[1], worst_idx[2]] - J_ref[worst_idx[1], worst_idx[2]]) /
        max(worst_scale, 1e-12)
    return passed, max_abs_dev, max_rel_dev, worst_idx
end

# Internal: allocate a backend's Jacobian target — sparse on its pattern when
# it carries one, dense otherwise (the `DenseAD_ForwardDiff` reference).
function _alloc_target(spec::JacobianSpec, n::Int)
    if spec.sparsity_phys === nothing
        return zeros(Float64, n, n)
    else
        S = spec.sparsity_phys .| spec.sparsity_geom
        return SparseMatrixCSC{Float64,Int}(
            S.m,
            S.n,
            copy(S.colptr),
            copy(S.rowval),
            zeros(Float64, length(S.nzval)),
        )
    end
end

# Rows whose stored diagonal the backend leaves as a structural zero (see
# `structural_zero_diagonal`); none unless the backend declares them.
_structural_zero_rows(::JacobianBackend, def, ws) = Int[]
_structural_zero_rows(backend::Structured, def, ws) =
    structural_zero_diagonal(ws.layout, def.mode, sparsity_profile(backend.strategy))
