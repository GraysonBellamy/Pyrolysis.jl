"""
The unified Jacobian interface.

Every Jacobian backend dispatches through `_compute_jacobian!`. Users never
call `_compute_jacobian!` directly — they construct a `JacobianSpec` and pass
it to `solve(..., jacobian = spec)`, which routes through `compute_jacobian!`.

Adding a new backend is a three-step contract:

  1. Define `struct MyBackend <: JacobianBackend end`
  2. Define `_build_jacobian_spec(::MyBackend, problem; kwargs...) -> JacobianSpec`
  3. Define `_compute_jacobian!(J, ::MyBackend, cache, sparsity_phys, sparsity_geom,
                                u, def, t, ws)`

Backend implementations live in the Jacobian module or in package
extensions. `MyBackend` lives in `Pyrolysis.Jacobian` (so `using Pyrolysis`
exposes the type); a package extension may provide the implementation if it
requires an optional dep.
"""

# ============================================================================
# Backend tag
# ============================================================================

"""
	JacobianBackend

Abstract supertype of every Jacobian-assembly strategy. The production
subtype is `Structured{S<:LinearSolveStrategy}`; `DenseAD_ForwardDiff`
is preserved as the verifier reference path.

Concrete subtypes should be lightweight bitstypes whenever possible — heavy
state belongs in the cache field of [`JacobianSpec`](@ref), not on the tag.
"""
abstract type JacobianBackend end

# ============================================================================
# Approximation trait
# ============================================================================

"""
	is_approximate(backend::JacobianBackend) -> Bool

Trait declaring whether `backend` deliberately omits known-non-zero entries
from its declared support pattern. Default: `false`.

Set to `true` for backends that approximate the Jacobian by *dropping*
couplings rather than computing them — e.g. `Structured(strategy = ApproxSolve())`,
which drops the structured ALE mesh-velocity coupling. Used by
[`verify_jacobian_at`](@ref) to decide whether to run the
**support-complete** strict check (assert that the reference Jacobian
is also zero outside the active backend's declared support) or fall
back to the masked check (compare only inside support).

The default `false` makes the strict check the norm: any backend that
"forgets" entries inside its declared support pattern AND has a sparsity
pattern that omits the same entries will fail verification, instead of
silently passing as both shortcuts cancel. Approximate backends opt in
to the masked compare via this trait.

# Adding the trait

```julia
struct MyApproximateBackend <: JacobianBackend end
Pyrolysis.Jacobian.is_approximate(::MyApproximateBackend) = true
```

"""
is_approximate(::JacobianBackend) = false

# ============================================================================
# JacobianSpec — single user-facing handle
# ============================================================================

"""
	JacobianSpec{B<:JacobianBackend, S, C}

The single user-facing handle for Jacobian configuration. Carried as a
keyword argument to `solve(...)`.

# Type parameters
- `B`: backend tag (concrete subtype of [`JacobianBackend`](@ref))
- `S`: sparsity-pattern type. `SparseMatrixCSC{Bool,Int}` for the
       `Structured` backend, `Nothing` for the verifier-only
       `DenseAD_ForwardDiff`. The two sparsity fields share this type
       — they must both be present or both absent for a given backend.
- `C`: backend-specific cache type. Holds Dual workspaces, JacobianConfigs,
       reusable scratch arrays, etc. Each backend declares its own cache
       struct; sharing is forbidden (different backends would race over
       a shared `build_ad_workspace` Dual workspace).

# Fields
- `backend`: instance of `B`
- `sparsity_phys`: physics-block (T, ξ, χ-row) sparsity. May be `nothing`.
- `sparsity_geom`: geometry-block (z-column) sparsity. ALE only; `nothing`
  otherwise.
- `cache`: backend-private workspace
- `colorvec_phys`: column coloring of `sparsity_phys` (for FD/sparse-AD).
                   `Int[]` when not applicable.
- `colorvec_geom`: column coloring of `sparsity_geom`.
                   `Int[]` when not applicable.
"""
struct JacobianSpec{B<:JacobianBackend,S,C}
    backend::B
    sparsity_phys::S
    sparsity_geom::S
    cache::C
    colorvec_phys::Vector{Int}
    colorvec_geom::Vector{Int}
end

# ============================================================================
# Dispatcher (the only public entry point besides the constructor)
# ============================================================================

"""
	compute_jacobian!(J, spec::JacobianSpec, u, def, t, ws)

Assemble the Jacobian into `J` using the backend carried by `spec`. This is
the single dispatch entry point: backends specialize on the backend argument
of [`_compute_jacobian!`](@ref), never on `compute_jacobian!` itself.

`def` is expected to be a [`ProblemDef`](@ref) and `ws` a [`Workspace`](@ref),
but the dispatcher itself does not assert these types — it forwards to the
backend, which is free to specialize on whatever problem-shape it understands.
"""
@inline compute_jacobian!(J, spec::JacobianSpec, u, def, t, ws) = _compute_jacobian!(
    J,
    spec.backend,
    spec.cache,
    spec.sparsity_phys,
    spec.sparsity_geom,
    u,
    def,
    t,
    ws,
)

"""
	_compute_jacobian!(J, backend::JacobianBackend, cache, sparsity_phys, sparsity_geom,
	                   u, def, t, ws)

Backend-private hook. **Backends must specialize on `backend`** with their
concrete subtype. The default fallback throws a clean `ErrorException`: if
you reach this method, you tried to use a backend whose implementation
hasn't been loaded (typically an extension that requires a missing
optional dep).
"""
function _compute_jacobian!(
    J,
    backend::JacobianBackend,
    cache,
    sparsity_phys,
    sparsity_geom,
    u,
    def,
    t,
    ws,
)
    error(
        "Jacobian.jl: no `_compute_jacobian!` method for backend " *
        "$(typeof(backend)). The backend is either part of an unloaded " *
        "package extension, or no implementation has shipped for this backend.",
    )
end

# ============================================================================
# Construction
# ============================================================================

"""
	_build_jacobian_spec(backend::JacobianBackend, problem; kwargs...) -> JacobianSpec

Backend-private constructor hook. Each backend implements this to build the
sparsity patterns, column colorings, and cache appropriate for its assembly
strategy. Default fallback raises a clean error mirroring
[`_compute_jacobian!`](@ref).
"""
function _build_jacobian_spec(backend::JacobianBackend, problem; kwargs...)
    error(
        "Jacobian.jl: no `_build_jacobian_spec` method for backend " *
        "$(typeof(backend)). The backend is either part of an unloaded " *
        "package extension, or no implementation has shipped for this backend.",
    )
end

"""
	JacobianSpec(backend::JacobianBackend, problem; kwargs...)
	JacobianSpec(backend::JacobianBackend, def::ProblemDef, ws::Workspace; kwargs...)

Two-form user-facing constructor. Both delegate to
[`_build_jacobian_spec`](@ref); each backend is responsible for building
its sparsity patterns and cache once, at construction time. **Do not call
inside the integrator hot path** — construction may build sparsity patterns
and run a warm-up Jacobian solve.

Use the form that matches what's already in scope:

- `JacobianSpec(backend, problem)` accepts a [`PyrolysisProblem`](@ref) and
  builds the canonical `(def, ws)` pair internally. This is the form most
  user code reaches for.
- `JacobianSpec(backend, def, ws)` accepts an already-built `(ProblemDef,
  Workspace)` pair. The solver and tests use this form because they
  already have the pair on hand. Each backend ultimately specializes on
  `_build_jacobian_spec(::B, ::ProblemDef, ::Workspace; ...)`.
"""
JacobianSpec(backend::JacobianBackend, problem; kwargs...) =
    _build_jacobian_spec(backend, problem; kwargs...)

JacobianSpec(backend::JacobianBackend, def::ProblemDef, ws::Workspace; kwargs...) =
    _build_jacobian_spec(backend, def, ws; kwargs...)

"""
	_build_jacobian_spec(backend::JacobianBackend, def::ProblemDef, ws::Workspace; kwargs...)

Three-argument backend hook. Mirrors the two-argument form but receives the
already-built `(def, ws)` pair directly — no `PyrolysisProblem` extraction.
Default fallback raises a clean error.
"""
function _build_jacobian_spec(
    backend::JacobianBackend,
    def::ProblemDef,
    ws::Workspace;
    kwargs...,
)
    error(
        "Jacobian.jl: no `_build_jacobian_spec` method for backend " *
        "$(typeof(backend)) with (ProblemDef, Workspace) arguments. The " *
        "backend is either part of an unloaded package extension, or no " *
        "implementation has shipped for this backend.",
    )
end

# ============================================================================
# Generic (PyrolysisProblem)-form constructor
# ============================================================================
#
# Bridges from the user-facing `JacobianSpec(backend, problem)` shape to
# each backend's `(def, ws)` implementation. The mode/radiation kwargs
# duplicate what `solve()` passes — they default to Eulerian / no-radiation
# so the most common case "Eulerian PMMA-style problem" needs no extra
# arguments. ALE callers must pass `use_ale = true`.
#
# For routine production work, prefer `solve(...; jacobian = nothing)` so
# the auto-default machinery picks the right backend for the active mode
# without requiring callers to thread mode flags through twice.

function _build_jacobian_spec(
    backend::JacobianBackend,
    problem::Problem.PyrolysisProblem;
    use_ale::Bool = false,
    radiation_model = nothing,
    radiation_intensity::Real = 0.0,
    energy_form::Symbol = :standard,
    use_transpiration_bc::Bool = false,
    use_blowing_correction::Bool = false,
    min_thickness::Real = 1e-6,
    kwargs...,
)
    rm = radiation_model === nothing ? NO_RADIATION : radiation_model
    def = Problem.to_problem_def(
        problem;
        use_ale = use_ale,
        radiation_model = rm,
        radiation_intensity = Float64(radiation_intensity),
        energy_form = energy_form,
        use_transpiration_bc = use_transpiration_bc,
        use_blowing_correction = use_blowing_correction,
    )
    ws = Problem.build_workspace(problem, def; min_thickness = min_thickness)
    return _build_jacobian_spec(backend, def, ws; kwargs...)
end
