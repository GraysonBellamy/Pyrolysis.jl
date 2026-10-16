"""
Benchmark harness.

Time and accuracy comparison across `JacobianBackend` implementations on a
single problem. **No `BenchmarkTools` dep** — the harness uses `Base.@elapsed`
and takes the minimum over `n_repeats` calls. For higher-precision results
(e.g. nanosecond-scale timing for the cell-local AD kernel),
`benchmark/Project.toml` scripts can wrap their own
`BenchmarkTools.@belapsed` around `compute_jacobian!` directly; this
harness is the convenience layer for end-to-end backend comparisons.

# Usage

```julia
using Pyrolysis
include("test/perf/perf_helpers.jl")
def, ws, u, _ = build_ale_no_radiation()

results = benchmark_jacobian_backends(def, ws, u;
    backends = [Structured(strategy = DirectSolve()),
                Structured(strategy = ApproxSolve())],
    reference = DenseAD_ForwardDiff(),
    n_repeats = 10)

# results is a Vector{NamedTuple} with fields:
#   name::String, assembly_seconds::Float64, bytes::Int,
#   max_abs_err::Float64, max_rel_err::Float64
```
"""

using SparseArrays: SparseMatrixCSC

"""
	benchmark_jacobian_backends(def, ws, u;
	        backends = [Structured(strategy = DirectSolve()),
	                    Structured(strategy = ApproxSolve())],
	        reference = DenseAD_ForwardDiff(),
	        n_repeats::Int = 10,
	        t::Float64 = 0.0) -> Vector{NamedTuple}

Run each backend on the given `(def, ws, u)` and report timing,
allocations, and accuracy against the reference backend. Returns one
`NamedTuple` per backend.

# Arguments
- `def, ws, u`: as returned by canonical builders in `perf_helpers.jl`.
- `backends`: vector of `JacobianBackend` instances.
- `reference`: a `JacobianBackend` for the accuracy comparison.
  Defaults to `DenseAD_ForwardDiff()`.
- `n_repeats`: number of timed repetitions; the minimum elapsed time is
  reported (filters scheduling jitter).
- `t`: simulation time at which to assemble.

# Output schema
Per backend:
- `name::String` — backend type name.
- `assembly_seconds::Float64` — minimum-over-`n_repeats` `@elapsed`.
- `bytes::Int` — `@allocated` after warm-up (zero is the goal for
  production backends).
- `max_abs_err::Float64` — vs reference Jacobian on this backend's
  structural support.
- `max_rel_err::Float64` — `max_abs_err` divided by
  `max(maximum(abs, J_ref), 1e-12)`.
"""
function benchmark_jacobian_backends(
    def,
    ws,
    u;
    backends = [Structured(strategy = DirectSolve()), Structured(strategy = ApproxSolve())],
    reference::JacobianBackend = DenseAD_ForwardDiff(),
    n_repeats::Int = 10,
    t::Float64 = 0.0,
)
    n = length(u)

    # Reference Jacobian (dense Float64) — built once.
    J_ref = zeros(n, n)
    ref_spec = JacobianSpec(reference, def, ws)
    compute_jacobian!(J_ref, ref_spec, u, def, t, ws)
    scale = max(maximum(abs, J_ref), 1e-12)

    rows = NamedTuple[]
    for backend in backends
        spec = JacobianSpec(backend, def, ws)
        J = _alloc_target(spec, n)

        # Warm-up: lazy AD workspace, FastSparsityMap initialization, etc.
        compute_jacobian!(J, spec, u, def, t, ws)

        # Accuracy.
        max_abs = _max_abs_diff(J, J_ref, spec)

        # Timing: minimum over n_repeats.
        t_min = Inf
        for _ = 1:n_repeats
            elapsed = @elapsed compute_jacobian!(J, spec, u, def, t, ws)
            t_min = min(t_min, elapsed)
        end

        # Allocations after warm-up.
        bytes = @allocated compute_jacobian!(J, spec, u, def, t, ws)

        push!(
            rows,
            (
                name = string(typeof(backend).name.name),
                assembly_seconds = t_min,
                bytes = bytes,
                max_abs_err = max_abs,
                max_rel_err = max_abs / scale,
            ),
        )
    end
    return rows
end

function benchmark_jacobian_backends(
    problem::Problem.PyrolysisProblem;
    use_ale::Bool = false,
    radiation_model = NO_RADIATION,
    radiation_intensity::Real = 0.0,
    kwargs...,
)
    def, ws, u = Problem.initialize_state(
        problem;
        use_ale = use_ale,
        radiation_model = radiation_model,
        radiation_intensity = Float64(radiation_intensity),
    )
    return benchmark_jacobian_backends(def, ws, u; kwargs...)
end

# Internal: allocate the Jacobian target matching this spec — sparse if the
# spec carries a structural pattern, dense otherwise.
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

# Internal: compute max abs error of the backend's Jacobian against the
# dense reference. Sparse backends are compared on their structural support.
function _max_abs_diff(J::Matrix{Float64}, J_ref::Matrix{Float64}, ::JacobianSpec)
    return maximum(abs, J .- J_ref)
end

function _max_abs_diff(
    J::SparseMatrixCSC{Float64,Int},
    J_ref::Matrix{Float64},
    spec::JacobianSpec,
)
    S = spec.sparsity_phys .| spec.sparsity_geom
    J_dense = Matrix(J)
    J_ref_masked = J_ref .* Matrix(S)
    return maximum(abs, J_dense .- J_ref_masked)
end
