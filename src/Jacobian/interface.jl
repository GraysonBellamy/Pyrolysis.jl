"""
The unified Jacobian interface.

Every Jacobian backend dispatches through `_compute_jacobian!`. Users never
call `_compute_jacobian!` directly: `solve(..., jacobian = backend)` builds the
`JacobianSpec` for its own `(def, ws)` and routes through `compute_jacobian!`.

Adding a new backend is a three-step contract:

  1. Define `struct MyBackend <: JacobianBackend end`
  2. Define `_build_jacobian_spec(::MyBackend, def::ProblemDef, ws::Workspace;
                                  kwargs...) -> JacobianSpec`
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
is the verifier reference path.

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

A backend bundled with the sparsity patterns and caches it assembles into,
built for one `(ProblemDef, Workspace)` pair. `solve` builds its own from its
`jacobian` backend; build one directly for [`compute_jacobian!`](@ref) or
[`verify_jacobian_at`](@ref).

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
"""
struct JacobianSpec{B<:JacobianBackend,S,C}
    backend::B
    sparsity_phys::S
    sparsity_geom::S
    cache::C
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
	JacobianSpec(backend::JacobianBackend, def::ProblemDef, ws::Workspace; kwargs...)

Build the spec for an already-built `(ProblemDef, Workspace)` pair by
delegating to [`_build_jacobian_spec`](@ref); each backend builds its
sparsity patterns and cache once, at construction time. **Do not call inside
the integrator hot path.** The spec binds `ws` (and its mesh), so it is valid
only for residual and Jacobian evaluations on that workspace; `solve` builds
its own.
"""
JacobianSpec(backend::JacobianBackend, def::ProblemDef, ws::Workspace; kwargs...) =
    _build_jacobian_spec(backend, def, ws; kwargs...)

"""
	_build_jacobian_spec(backend::JacobianBackend, def::ProblemDef, ws::Workspace; kwargs...)

Backend-private constructor hook. Each backend implements this to build the
sparsity patterns and cache appropriate for its assembly strategy. The
default fallback raises a clean error.
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
