# 11. Solver and Jacobian Configuration; Extensions

This chapter is the reference for tuning *how* Pyrolysis.jl solves a problem,
as opposed to *what* problem it solves (Chapters 3–7) or *how to read the
result* (Chapter 9). It covers choosing a Jacobian backend and linear-solve
strategy by problem size and simulation mode, configuring the integrator's
linear solver (sparse-direct or ILU-preconditioned Krylov), loading the optional
vendor and sensitivity extensions, debugging Jacobian correctness with
`verify_jacobian_at`, and dropping down to the pure-residual layer
(`ProblemDef` / `Workspace` / `residual!`) when you need to build a custom
solver loop. The mathematics behind the structured Jacobian and the implicit-RK
`W` matrix is derived in Technical Reference §15; this chapter is about the
knobs.

Two distinct linear systems appear in a stiff implicit solve, and Pyrolysis.jl
lets you configure each one separately. Keeping them straight is the key to this
chapter:

1. **The Jacobian backend** (`jacobian = ...` keyword) decides *how the matrix*
   `J = ∂(du/dt)/∂u` *is assembled* — by structured operator assembly (the
   default: `solve` auto-builds `Structured(strategy = DirectSolve())`, or
   `ApproxSolve` under ALE with depletion handling, when `jacobian = nothing`),
   or by dense automatic differentiation (verification only).
2. **The integrator's linear solver** (`linear_solver = ...` keyword) decides
   *how the implicit-stage system* `W x = b`, with `W = I/(γ_RK Δt) − J`, *is
   solved* once `J` exists — by sparse-direct factorization (KLU, UMFPACK,
   Sparspak, Pardiso) or by a preconditioned Krylov iteration (ILU-GMRES, ILU-BiCGSTAB).

You can mix and match: a pluggable structured Jacobian with a sparse-direct
linear solver, an integrator-built finite-difference Jacobian with KLU, and so
on. The defaults (auto-built structured Jacobian; KLU or UMFPACK chosen from its pattern) work well for
the typical 1D cone-calorimeter problem and you can ignore this chapter until a
solve is too slow or you need sensitivities.

---

## 11.1 The two layers, and what the defaults do

When you call `solve(problem)` with no Jacobian or linear-solver keywords, the
following happens (see `src/Solver/solve.jl`):

- The integrator defaults to `KenCarp4(autodiff = AutoFiniteDiff())` with the
  default linear solver below. `KenCarp4` is an SDIRK (singly
  diagonally implicit Runge–Kutta) method well suited to the moderate stiffness
  of Arrhenius pyrolysis. `solve` always supplies the Jacobian, so the
  `AutoFiniteDiff()` setting never builds `J`; it only affects methods that also
  differentiate the residual in time (the Rosenbrock `∂f/∂t`).
- The Jacobian defaults to the **structured analytical backend**: with
  `jacobian = nothing` (the default), `solve` constructs
  `Structured(strategy = DirectSolve())` internally — for both Eulerian and
  ALE problems — and hands the integrator an exact sparse Jacobian plus its
  sparsity prototype. Under ALE with `handle_depletion = true` the default is
  `Structured(strategy = ApproxSolve())` (§10.8).
- The linear solver follows the Jacobian's sparsity pattern:
  `UMFPACKFactorization()` for the ALE geometry block that `DirectSolve`
  materializes under ALE, `KLUFactorization()` for every other sparse
  Jacobian (§11.4.1).

In short: by default, **Layer 1 is the structured analytical Jacobian**
(`Structured(strategy = DirectSolve())`, built for you) and **Layer 2 is
KLU or UMFPACK**, chosen from that Jacobian's pattern. The `jacobian` keyword is for *overriding* the strategy
(`ApproxSolve`, `StructuredSolve`, `MatrixFreeSolve`) or supplying a custom
spec — not for switching the structured machinery on.

```julia
using Pyrolysis

# The default solve: auto-built structured Jacobian + KLU/UMFPACK linear solve.
sol = solve(problem)
```

> Defaults at a glance: `integrator = nothing` → `KenCarp4`,
> `linear_solver = nothing` → `KLUFactorization()` (`UMFPACKFactorization()`
> for the exact ALE Jacobian),
> `jacobian = nothing` → auto-built `Structured(strategy = DirectSolve())`
> (`ApproxSolve` under ALE with `handle_depletion = true`),
> `abstol = 1e-8`, `reltol = 1e-6`.

---

## 11.2 Choosing a Jacobian backend

A **Jacobian backend** is supplied as a `JacobianSpec` to the `jacobian`
keyword of `solve`. Construct it from a backend tag plus your problem:

```julia
spec = JacobianSpec(Structured(strategy = DirectSolve()), problem)
sol  = solve(problem; jacobian = spec)
```

There are two production-relevant backend tags, both exported at top level:

| Backend tag | Purpose | Use in production? |
|---|---|---|
| `Structured(strategy = ...)` | Operator-based structured assembly (fast, exact or approximate per strategy) | Yes |
| `DenseAD_ForwardDiff()` | Whole-residual ForwardDiff; dense Jacobian | No — verification reference only |

`Structured` is parameterized by a **linear-solve strategy** (§11.3) that
determines how the assembled structured Jacobian is inverted. The default,
selected when you pass `jacobian = nothing`, is `Structured(strategy =
DirectSolve())`.

`DenseAD_ForwardDiff()` runs ForwardDiff over the entire residual vector and
returns a dense matrix. It is correct by construction but slow; its only role is
as the reference against which the structured backend is verified (see §11.7).
Do not use it for production solves.

### 11.2.1 The `JacobianSpec` constructor

`JacobianSpec` has two forms (`src/Jacobian/interface.jl`):

```julia
# From a PyrolysisProblem (the common form). Pass use_ale = true for ALE.
JacobianSpec(backend, problem; use_ale = false, radiation_model = nothing,
             radiation_intensity = 0.0, energy_form = :standard,
             use_transpiration_bc = false, min_thickness = 1e-6)

# From an already-built (ProblemDef, Workspace) pair (low-level form).
JacobianSpec(backend, def, ws)
```

Important: a `JacobianSpec` builds its sparsity pattern and AD workspaces once,
at construction time, and it bakes in the simulation mode. **If your solve uses
ALE, you must tell the spec by passing `use_ale = true`** — otherwise the
sparsity pattern omits the geometry (z-column) block and the Jacobian will be
wrong:

```julia
# WRONG for an ALE solve — the spec is built for the Eulerian mode.
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem)
solve(problem; jacobian = spec, use_ale = true)   # mode mismatch

# RIGHT — the spec's mode matches the solve.
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
solve(problem; jacobian = spec, use_ale = true)
```

If a problem also uses radiation, pass the matching `radiation_model` (e.g.
`radiation_model = BEER_LAMBERT`) so the spec includes the Beer–Lambert
non-local sparsity. Similarly, pass `energy_form = :with_generation_sink` if you
solve with that energy form. For most users the simplest correct approach is to
let `solve` build the spec for you — see §11.2.3.

### 11.2.2 Building the spec for the right mode

A `JacobianSpec` is built for one state layout. `JacobianSpec(backend,
problem)` defaults to the Eulerian layout with no radiation; an ALE solve or a
Beer–Lambert solve needs the same `use_ale` / `radiation_model` keywords the
solve uses, and `solve` refuses a spec whose state length disagrees with the
problem's:

```julia
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem;
                    use_ale = true, radiation_model = BEER_LAMBERT)
solve(problem; use_ale = true, radiation_model = BEER_LAMBERT, jacobian = spec)
```

The simpler form is to pass the bare backend and let `solve` build the spec
for its own mode:

```julia
solve(problem; use_ale = true, jacobian = Structured(strategy = ApproxSolve()))
```

### 11.2.3 Letting `solve` pick the backend

For routine work, the cleanest pattern is to pass `jacobian = nothing` (the
default) and let `solve` construct the default structured Jacobian for you, *or* to
opt into a non-default backend explicitly only when you need it. The auto-default when
`jacobian` is left at its default is `Structured(strategy = DirectSolve())`
constructed internally with the correct mode flags — you do not have to thread
`use_ale` / `radiation_model` through `JacobianSpec` yourself. Supply your own
spec only when you want a non-default strategy (`ApproxSolve`, `StructuredSolve`,
`MatrixFreeSolve`) or a custom reference for verification.

---

## 11.3 Linear-solve strategies for the structured Jacobian

The `strategy` argument to `Structured(...)` is a `LinearSolveStrategy`. It
governs how the structured Jacobian — a sparse local block `A` plus a list of
fast-apply non-local couplings (prefix-sum, suffix-sum, rank-1) — is inverted
inside the Newton step. All four strategies are top-level exports.

| Strategy | What it does | Exact? | Best for |
|---|---|---|---|
| `DirectSolve()` | Densify the structured couplings into `A`, factor the full sparse matrix (KLU/UMFPACK) | Exact | Small/typical problems; the verification anchor |
| `StructuredSolve()` | Sparse-LU on `A` + preconditioned Richardson over the couplings (2–5 iters) | Exact | Large problems where dense assembly is costly but exactness is required |
| `ApproxSolve()` | Drop the structured couplings; solve `A⁻¹ b` only (Schur-style approximation) | Approximate | Large ALE problems where speed dominates (~7.8× faster assembly) |
| `MatrixFreeSolve()` | GMRES over the operator action without densifying | Exact (iterative) | Research / very large `N`; benchmark before relying on it |

Construct any of these as a parameterless instance:

```julia
JacobianSpec(Structured(strategy = DirectSolve()),     problem)               # default
JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
JacobianSpec(Structured(strategy = ApproxSolve()),     problem; use_ale = true)
JacobianSpec(Structured(strategy = MatrixFreeSolve()), problem; use_ale = true)
```

### 11.3.1 `DirectSolve` — the default and the anchor

`DirectSolve` materializes `A + Σⱼ densify(couplingⱼ)` and factors it directly.
It is exact and predictable. For Eulerian problems and modest ALE meshes it is
the right choice and the one you get when you do nothing. It is also the backend
you should verify against first when proving a configuration is correct.

```julia
spec = JacobianSpec(Structured(strategy = DirectSolve()), problem)
sol  = solve(problem; jacobian = spec)
```

### 11.3.2 `StructuredSolve` — exact, cheaper assembly at scale

`StructuredSolve` keeps the solve exact but avoids densifying the non-local
ALE/radiation couplings. It factors only the sparse local block `A`, then runs a
preconditioned Richardson iteration that applies each coupling through its O(N)
`apply!` method:

```
F = lu(A); y = F \ b
repeat:  r = b − A·y − Σⱼ Cⱼ·y;  Δy = F \ r;  y ← y + Δy   (2–5 iterations)
```

Convergence is fast (2–5 iterations to machine precision) when the coupling
magnitudes are small relative to `A`'s diagonal — which is the usual ALE regime.
Choose `StructuredSolve` when an ALE mesh is large enough that the dense
assembly of `DirectSolve` dominates, but you still need an exact Newton step.

### 11.3.3 `ApproxSolve` — fastest, for large ALE

`ApproxSolve` drops the structured couplings entirely and solves `A⁻¹ b`. The
Newton step is then a Schur-style approximation: it omits the dense ALE
mesh-velocity coupling block, takes a few extra Newton iterations, but each
linear solve and the assembly itself are far cheaper (the assembly profile drops
to `:ale_reduced`, skipping the cumulative ALE widening, the geometry AD pass,
the Beer–Lambert suffix, and the χ-column rank-1 term). This is the strategy the
Douglas-fir ALE example uses for roughly an 8× speedup over the full Jacobian.

```julia
# Large ALE cone problem — fast approximate Jacobian.
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true)
```

`ApproxSolve` is the **only** strategy that deliberately omits known-nonzero
entries; it carries the `is_approximate` trait, which switches
`verify_jacobian_at` (§11.7) into masked-comparison mode so the dropped entries
are reported informationally rather than counted as failures.

### 11.3.4 `MatrixFreeSolve` — research / very large N

`MatrixFreeSolve` runs GMRES over the structured-Jacobian operator action
without ever forming `A + C`. It is exact in the limit of GMRES convergence and
useful when even materializing the coupled operator is expensive. It is marked
experimental in the source; benchmark a specific problem before using it in
production, and prefer `DirectSolve`/`StructuredSolve` unless the matrix-free
path measurably wins.

### 11.3.5 A sizing rule of thumb

- Eulerian (fixed mesh), any size, no depletion → `jacobian = nothing`
  (default) is fine; it already gives you
  `Structured(strategy = DirectSolve())`.
- ALE, small/medium mesh → `DirectSolve`.
- ALE, large mesh, exactness required → `StructuredSolve`.
- ALE, large mesh, speed over exactness → `ApproxSolve`.
- Depletion tracking → no pluggable Jacobian (`jacobian = nothing`).

---

## 11.4 Configuring the integrator's linear solver

Independently of the Jacobian backend, you choose how the per-stage implicit
system `W x = b` (`W = I/(γ_RK Δt) − J`, see Technical Reference §15.5) is
solved. This is the `linear_solver` keyword. There are two families:
sparse-direct factorizations and ILU-preconditioned Krylov iterations.

### 11.4.1 Sparse-direct factorizations

With `linear_solver = nothing` the factorization follows the Jacobian's
sparsity pattern:

| Jacobian | Default factorization |
|---|---|
| `DirectSolve` under ALE (the exact ALE Jacobian, with its dense geometry block) | `UMFPACKFactorization()` |
| Every other `Structured` Jacobian: `ApproxSolve` (the ALE default with `handle_depletion = true`) and `DirectSolve` on an Eulerian mesh | `KLUFactorization()` |
| A dense backend (`DenseAD_ForwardDiff`) | `LUFactorization()` |

On 100-cell PMMA and charring-wood cone decks and on the V-BO-01 burnout case,
KLU runs 2.6–4.4× faster end to end than Sparspak on the local patterns, with
step counts equal to within 1%, and UMFPACK runs the exact ALE Jacobian
1.5–1.9× faster than Sparspak (Technical Reference §15.5.1). The Eulerian
`BEER_LAMBERT` pattern of `DirectSolve` stays on KLU: it is the fastest there,
and UMFPACK's threshold pivoting returned inaccurate solves on it. Every
backend stays selectable by pairing it with a *constructor* integrator:

```julia
using LinearSolve: KLUFactorization, UMFPACKFactorization, SparspakFactorization

sol = solve(problem; linear_solver = UMFPACKFactorization())
sol = solve(problem; linear_solver = SparspakFactorization())   # pure-Julia sparse LU
```

> Pairing rule: `linear_solver` is combined with `integrator` only when
> `integrator` is a *constructor* such as `KenCarp4` (or left at its default).
> If you want a fully constructed algorithm with a non-default linear solver,
> build it yourself and pass it as the `integrator`:
> `solve(problem; integrator = KenCarp4(linsolve = KLUFactorization()))`.
> Passing both a fully constructed `integrator` and a separate `linear_solver`
> raises an error.

A note on dense Jacobian backends: `DenseAD_ForwardDiff` has no sparse
`jac_prototype`, so no sparse factorization can run. `solve` then defaults
the linear solver to a dense `LUFactorization()` — and if you explicitly
request a sparse factorization it warns and rebuilds your algorithm with
`LUFactorization` instead. Every `Structured` Jacobian, including the
`ApproxSolve` default under depletion handling, carries the sparse
prototype.

### 11.4.2 ILU-preconditioned GMRES

For large ALE problems where direct factorization of `W` dominates the cost,
Pyrolysis.jl ships an ILU-preconditioned GMRES solver,
`ILUGMRESFactorization` (exported by the `Solver` submodule; qualify as
`Pyrolysis.ILUGMRESFactorization` or `using Pyrolysis.Solver`). It builds an incomplete-LU
preconditioner from the concrete sparse `W` and iterates GMRES against it,
falling back to diagonal scaling if the ILU factorization fails on a stiff
transient.

```julia
# ILUGMRESFactorization(; τ = 0.1, gmres_rtol = 1e-8, restart = 30)
alg = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1))

sol = solve(problem;
    integrator = alg,
    jacobian   = JacobianSpec(Structured(strategy = DirectSolve()), problem; use_ale = true),
    use_ale    = true)
```

Tuning parameters and their defaults (the field name `τ` is the ILU drop
tolerance, written `τ_ILU` in the nomenclature):

| Parameter | Default | Meaning |
|---|---|---|
| `τ` | `0.1` | ILU drop tolerance. **Smaller** → denser ILU → better preconditioner, slower setup. |
| `gmres_rtol` | `1e-8` | GMRES relative convergence tolerance. |
| `restart` | `30` | GMRES restart frequency (memory parameter). |

Caveat: inexact linear solves can increase Newton rejections on Arrhenius
transients, which sometimes *raises* total wall time despite cheaper individual
solves. Always benchmark an ILU-GMRES configuration end-to-end against the
default direct factorization before adopting it for a production study.

### 11.4.3 ILU-preconditioned BiCGSTAB

`ILUBiCGSTABFactorization` (exported by the `Solver` submodule; qualify or
`using Pyrolysis.Solver`) is the BiCGSTAB analogue. It uses
the same ILU preconditioner but a fixed-memory BiCGSTAB iteration (no restart
parameter), which can be faster for some matrix structures:

```julia
# ILUBiCGSTABFactorization(; τ = 0.1, rtol = 1e-8)
alg = KenCarp4(linsolve = ILUBiCGSTABFactorization(τ = 0.1, rtol = 1e-8))
sol = solve(problem; integrator = alg,
    jacobian = JacobianSpec(Structured(strategy = DirectSolve()), problem; use_ale = true),
    use_ale  = true)
```

| Parameter | Default | Meaning |
|---|---|---|
| `τ` | `0.1` | ILU drop tolerance. |
| `rtol` | `1e-8` | BiCGSTAB relative convergence tolerance. |

### 11.4.4 Linear solver and depletion

Depletion handling keeps the state size, so the sparse factorization keeps
its pattern through every park event and nothing about the linear solver
changes with `handle_depletion`. The dense LU is used only with the dense
`DenseAD_ForwardDiff` backend, which attaches no sparse prototype; an explicitly
requested sparse factorization is then rebuilt to `LUFactorization` inside the
user's algorithm, with a warning.

---

## 11.5 Tolerances and integrator selection

The ODE error tolerances are independent of both Jacobian and linear-solver
choices:

| Keyword | Default | Meaning |
|---|---|---|
| `abstol` | `1e-8` | Absolute error tolerance. |
| `reltol` | `1e-6` | Relative error tolerance. |
| `maxiters` | `1000000` | Maximum integration steps. |
| `dtmin` | `1e-12` | Minimum time step [s]. |
| `dtmax` | `Inf` | Maximum time step [s]. |

To use a different stiff integrator, pass its constructor (the linear solver and
`AutoFiniteDiff` autodiff are wired in for you) or a fully built instance:

```julia
using OrdinaryDiffEq: Rodas5P, FBDF

# Constructor form — linear_solver is paired in automatically:
sol = solve(problem; integrator = Rodas5P, linear_solver = SparspakFactorization())

# Fully constructed form (do not also pass linear_solver):
sol = solve(problem; integrator = FBDF(linsolve = KLUFactorization()))
```

Loosening tolerances (`reltol = 1e-4`) speeds up exploratory runs; tighten them
for convergence and validation studies (Technical Reference §18.5). Tolerance
choice does not affect mass/energy *bookkeeping* closure — that closes to the
integrator tolerance by construction (Technical Reference §16) — but it does
affect physical accuracy.

---

## 11.6 Loading and using extensions

Optional functionality is delivered through Julia's package-extension system
(Julia 1.9+). The hooks (`klu_factorization`, `pardiso_factorization`,
`amg_preconditioner`, `forward_sensitivity`, `adjoint_sensitivity`,
`verify_mixing_derivatives`) are **exported stubs** in core Pyrolysis.jl. Each
stub errors with a load hint until you `using` the corresponding optional
package, at which point the extension attaches a concrete method. No extension
depends on another; all are independent peers.

| Extension | Trigger (`using ...`) | Attaches | Optional dep(s) |
|---|---|---|---|
| `PyrolysisKLUExt` | `KLU` | `klu_factorization()` | `KLU` |
| `PyrolysisPardisoExt` | `Pardiso` | `pardiso_factorization()` | `Pardiso` |
| `PyrolysisAMGExt` | `AlgebraicMultigrid` | `amg_preconditioner(A)` | `AlgebraicMultigrid` |
| `PyrolysisSciMLSensitivityExt` | `SciMLSensitivity, Zygote` | `forward_sensitivity`, `adjoint_sensitivity` | `SciMLSensitivity`, `Zygote` |
| `PyrolysisSymbolicsExt` | `Symbolics` | `verify_mixing_derivatives` | `Symbolics` |
| `PyrolysisPlotsExt` | `RecipesBase`/`Plots` | `plot(sol; kind = ...)` recipes | `RecipesBase` |

### 11.6.1 Linear-solver extensions: KLU, Pardiso, AMG

These three extensions each expose a vendor solver under a stable
Pyrolysis-level name so the benchmark harness and out-of-tree code can request
it without core Pyrolysis taking a hard dependency.

```julia
using Pyrolysis
using KLU                              # triggers PyrolysisKLUExt

ls  = Pyrolysis.klu_factorization()    # -> LinearSolve.KLUFactorization()
sol = solve(problem; linear_solver = ls)
```

```julia
using Pyrolysis
using Pardiso                          # triggers PyrolysisPardisoExt

# Returns LinearSolve.MKLPardisoFactorize(); for users with an Intel MKL license.
ls  = Pyrolysis.pardiso_factorization()
sol = solve(problem; integrator = KenCarp4(linsolve = ls))
```

```julia
using Pyrolysis
using AlgebraicMultigrid               # triggers PyrolysisAMGExt
using SparseArrays, LinearAlgebra

# Smoothed-aggregation AMG preconditioner for a sparse matrix A, for use as a
# preconditioner inside a Krylov linsolve.
A    = sprand(100, 100, 0.05) + 100I   # illustrative sparse matrix
prec = Pyrolysis.amg_preconditioner(A)
```

Note that `KLUFactorization` is already reachable as
`LinearSolve.KLUFactorization` (LinearSolve is a hard dependency); the KLU
extension simply confirms KLU itself is loaded and exposes it under the stable
name `Pyrolysis.klu_factorization`. The AMG preconditioner is meant for the
preconditioner slot of a Krylov solver, not as a `linear_solver` on its own.

**Stub behavior without the trigger.** If you call any of these before loading
the optional package, you get a `MethodError` pointing at the unloaded
extension:

```julia
using Pyrolysis
Pyrolysis.klu_factorization()   # MethodError — `using KLU` not yet loaded
```

The fix is always the same: add the corresponding `using` line.

### 11.6.2 Sensitivity extension

`forward_sensitivity` and `adjoint_sensitivity` are stubs until you load
`SciMLSensitivity` and `Zygote`, which activate `PyrolysisSciMLSensitivityExt`.
Both take a `p_inject` closure that builds a fresh problem from a parameter
vector `θ` (it must not mutate the base problem). The closure rebuilds the
problem from the base problem's mesh and BCs with a material that carries the
injected parameters; the exact way you parameterize the material is up to you
(Chapter 13 develops realistic `p_inject` closures in full):

```julia
using Pyrolysis
using SciMLSensitivity, Zygote          # triggers PyrolysisSciMLSensitivityExt

# p_inject(base, θ) -> fresh PyrolysisProblem with θ injected.
# `make_material(θ)` is your own function that builds a Material from the
# parameter vector (e.g. setting reaction A, E). The fresh problem reuses the
# base problem's mesh, BCs, time span, and initial conditions.
function p_inject(base, θ)
    mat = make_material(θ)
    return PyrolysisProblem(base.mesh, mat, base.bc_set, base.tspan;
                            T_initial = base.T_initial,
                            output_times = base.output_times)
end

# Forward mode (ForwardDiff) — small parameter counts:
θ = [1e13, 2.0e5]
output, ∂out_∂θ = forward_sensitivity(base_problem, θ, p_inject;
                                      output_fn = sol -> maximum(sol.surface_T))

# Adjoint mode (Zygote + SciMLSensitivity) — large parameter counts, scalar loss:
loss, grad = adjoint_sensitivity(base_problem, θ, p_inject, sol -> sum(sol.surface_T);
                                 sensealg = InterpolatingAdjoint(autojacvec = ReverseDiffVJP()))
```

This chapter only covers the *loading* mechanics; the full sensitivity workflow
(forward vs. adjoint choice, inverse analysis, Sobol global sensitivity) is the
subject of Chapter 13, and the theory is in Technical Reference §17.

### 11.6.3 Symbolics verification extension (test-only)

`verify_mixing_derivatives` is a test-only stub activated by `using Symbolics`.
It compares the package's analytical mixing-rule derivatives against
Symbolics-derived ground truth. Production code never calls it; it exists to
catch drift bugs in the mixing-rule Jacobian terms.

```julia
using Pyrolysis
using Symbolics                         # triggers PyrolysisSymbolicsExt

# Returns a per-rule max relative error (NaN for BRUGGEMAN, which has no closed form).
report = verify_mixing_derivatives(material; rtol = 1e-8, T_test = 600.0)
```

---

## 11.7 Debugging Jacobian correctness with `verify_jacobian_at`

When a custom or approximate Jacobian configuration produces suspicious Newton
behavior (excessive rejections, divergence, or wrong steady state), compare it
against the dense ForwardDiff reference at a state point. There are two ways: a
one-shot call, and an in-solve callback.

### 11.7.1 One-shot verification

`verify_jacobian_at` takes an active spec, a `(def, ws)` pair, a state `u`, and
a time `t`, and reports the deviation from the `DenseAD_ForwardDiff` reference.
Its signature is `verify_jacobian_at(spec, def, ws, u, t; reference =
DenseAD_ForwardDiff(), rtol = 1e-6, atol = 1e-10)`:

```julia
using Pyrolysis

# Build the (def, ws, u) triple for the mode you intend to solve in
# (initialize_state applies the problem's ICs before packing u).
def, ws, u = initialize_state(problem; use_ale = true)

spec = JacobianSpec(Structured(strategy = ApproxSolve()), def, ws)

report = verify_jacobian_at(spec, def, ws, u, 0.0;
                            reference = DenseAD_ForwardDiff(),
                            rtol = 1e-6, atol = 1e-10)

report.passed                      # Bool — strict: per-entry |ΔJ_ij| ≤ atol + rtol·max(|J_ref[i,j]|, row max);
                                   #        masked: max_abs_dev ≤ atol || max_rel_dev ≤ rtol
report.max_abs_dev                 # largest element-wise deviation under the mode's rule
report.max_rel_dev
report.worst_index                 # (row, col) of the worst deviation
report.mode                        # :strict (exact backends) or :masked (approximate)
report.max_abs_dev_in_support      # in-support disagreement
report.max_abs_dev_outside_support # reference magnitude the backend's support omits
```

Comparison mode is chosen automatically by the backend's `is_approximate`
trait (§11.3.3):

- **Strict / support-complete** (exact backends, e.g. `DirectSolve`,
  `StructuredSolve`): every entry of the dense reference is compared against the
  active backend; a non-zero reference entry that the backend's support omits
  surfaces as a deviation, and its magnitude appears in
  `max_abs_dev_outside_support`. The strict criterion is per-entry with a
  row-magnitude floor, so a missing coupling in a small-magnitude row cannot
  hide behind the global matrix maximum. This is the regime exact backends
  must satisfy.
- **Masked** (`ApproxSolve`): entries outside the backend's declared support are
  zeroed in the diff; the omitted-entry magnitudes are reported in
  `max_abs_dev_outside_support` but not asserted. Use `is_approximate(spec.backend)`
  to confirm which mode you are in.

`verify_jacobian_at` reports deviations; it does not throw. You inspect
`report.passed` and the per-entry fields to localize a discrepancy.

### 11.7.2 In-solve verification callback

To sample Jacobian agreement at several times during an actual solve, pass the
`verify_jacobian` keyword as a `NamedTuple` with keys drawn from `{:rtol, :atol,
:sample_times}`. This installs a callback that compares the active backend
against the dense reference each time the integrator crosses one of the listed
times and **warns** (does not abort) on a mismatch:

```julia
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
sol  = solve(problem;
    jacobian        = spec,
    use_ale         = true,
    verify_jacobian = (rtol = 1e-6, atol = 1e-10, sample_times = [0.0, 1.0, 10.0]))
```

The callback only runs when a pluggable `jacobian` is supplied (there is nothing
to verify against the integrator's own finite-difference Jacobian). It warns
rather than throws so a single transient outlier does not abort an otherwise
healthy production run.

---

## 11.8 The pure-residual layer: building a custom solver loop

Beneath `solve`, Pyrolysis.jl exposes a decoupled residual API so you can drive
the time integration yourself, build a non-standard solver, or feed the residual
to external tooling. The three pieces are all top-level exports:

- `ProblemDef` — an **immutable** frozen problem definition (material, mesh
  topology, BCs, initial conditions, simulation-mode traits, options). Safe to
  share across threads and sweep workers.
- `Workspace` — a **mutable** bundle of per-step caches, work arrays, and the
  geometry-from-state machinery. Allocate **one per worker**.
- `residual!(du, u, def, t, ws)` — the pure ODE right-hand side. Note the
  argument order: `du`, `u`, `def`, `t`, `ws`.

The canonical builder is `initialize_state`, which turns a `PyrolysisProblem`
into this lower-level form in one call — applying the problem's initial
conditions, freezing the definition, allocating the workspace, and packing
the initial state (the same sequence `solve()` runs internally):

```julia
using Pyrolysis
import OrdinaryDiffEq                                  # from your environment

# 1. ICs → immutable def → mutable workspace → packed u0 [T, ξ, z (ALE), χ (WithChi)].
def, ws, u0 = initialize_state(problem;
    use_ale         = false,
    radiation_model = NO_RADIATION,
    energy_form     = :standard)

# 2. Wrap residual! in the (du, u, p, t) shape OrdinaryDiffEq expects.
#    Note the residual! argument order is (du, u, def, t, ws).
rhs!(du, u, p, t) = residual!(du, u, def, t, ws)

# 3. Hand-roll your own integration.
prob = OrdinaryDiffEq.ODEProblem(rhs!, u0, problem.tspan)
raw  = OrdinaryDiffEq.solve(prob, OrdinaryDiffEq.KenCarp4(autodiff = OrdinaryDiffEq.AutoFiniteDiff()))
```

The underlying primitives `to_problem_def`, `build_workspace`, and
`create_state_vector` remain public as the advanced layer (AD workspaces via
`build_ad_workspace`, mid-simulation repacking via `pack_state!`, warm
starts). If you compose them by hand, the IC contract is on you:
`create_state_vector` packs the mesh **as-is**, so call
`Pyrolysis.Problem.apply_initial_conditions!(problem)` first — otherwise `u0`
silently carries the mesh constructor's placeholder state instead of the
problem's `T_initial`/`ξ_initial`, an error no retcode will flag.

`StateLayout` (also exported) describes the packing and provides accessors so
you can read fields out of a raw state vector without re-deriving the offsets:

```julia
layout = ws.layout

# Index helpers and block accessors:
state_length(layout)                       # total state-vector length
get_temperature(u, layout, i)              # T in cell i
get_concentration(u, layout, i, j)         # ξ of component j in cell i (cell index first)
get_node_positions(u, layout)              # z nodes (ALE only)
get_thickness_from_state(u, layout)        # current sample thickness
block_view(u, layout, ...)                 # a view of one block

# For variable-cross-section runs:
column_chi_bar_from_state(u, layout, material)        # column-average χ̄
cross_section_area_from_state(mesh, material, u, layout)
```

You can also attach a structured Jacobian to your hand-rolled loop. Build a spec
from the same `(def, ws)` pair and call `compute_jacobian!`:

```julia
spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
J    = zeros(state_length(layout), state_length(layout))
compute_jacobian!(J, spec, u0, def, 0.0, ws)
```

When driving the solver yourself, remember the contract that `solve` normally
enforces for you: **the residual and Jacobian assembly are mesh-read-only**. The
mutable mesh state is only meant to be synchronized to the integrator's `u`
between accepted steps (this is what the built-in `StepUpdateCallback` does
inside `solve`). If your custom loop relies on the mesh geometry being current
for post-processing or AMR, you must call
`update_mesh_from_state!(mesh, u, layout)` between steps yourself.

---

## 11.9 The `Internal` umbrella and importing trait subtypes

The top-level `Pyrolysis` namespace deliberately exports only a curated public
surface. Two escape hatches give you the rest.

### 11.9.1 `Pyrolysis.Internal`

`Pyrolysis.Internal` re-exports **every** symbol exported by the core
submodules (Materials, BoundaryConditions, Geometry, Properties, Physics,
Adaptivity, Discretization, Problem, Residual, Jacobian, Diagnostics, Solver) in
one namespace:

```julia
using Pyrolysis.Internal   # everything the core submodules export, in one place
```

Use this when you need an internal building block — a mixing function, a flux
kernel, a cache type — that is not on the curated top-level list. Caveat:
symbols reached through `Internal` are **not** part of the stable public API and
may move or change between minor releases; prefer top-level exports in
production code. (The `Experimental` namespace — P1 radiation, h-AMR — is *not*
included in `Internal`; reach it explicitly as `Pyrolysis.Experimental.X`.)

### 11.9.2 Importing the simulation-mode trait subtypes

The `SimulationMode` *abstract* trait is exported, but the concrete subtypes —
`Eulerian`, `ALE`, `Fickian`, `DarcyFick`, `NoChi`, `WithChi`, `NoRadiation`,
`SurfaceAbsorption`, `BeerLambert` — are **not** re-exported at top level, to
keep the public surface clean and avoid shadowing generic names like `ALE`. If
you need to dispatch on or construct one (e.g. for a custom residual
specialization), reach it through the `Problem` submodule:

```julia
# Fully qualified:
mode_tag = Pyrolysis.Problem.Eulerian()

# Or import the ones you need explicitly:
using Pyrolysis.Problem: Eulerian, ALE, Fickian, DarcyFick, WithChi
```

In normal usage you never construct these directly — `solve` and
`to_problem_def` select the right mode bundle from your keyword flags
(`use_ale`, `radiation_model`, the presence of permeability for `DarcyFick`, a
lateral-shrinkage law for `WithChi`). You only reach for the trait subtypes when
building custom dispatch at the residual layer (see Chapter 7 for how the trait
axes map to the four residual paths, and Technical Reference §15.1).

> Note: P1 radiation (`P1_QUASI_STEADY`) is **not** part of the unified
> residual; it lives in `Pyrolysis.Experimental` with its own residual wrapper
> and is rejected by `solve` (`radiation_model = P1_QUASI_STEADY` raises an
> error). The supported radiation models for `solve` are `NO_RADIATION`,
> `SURFACE_ABSORPTION`, and `BEER_LAMBERT` (Technical Reference §8.12).

---

## 11.10 Quick configuration recipes

```julia
# (a) Default, simplest correct solve — Eulerian, no depletion.
sol = solve(problem)

# (b) Explicit exact structured Jacobian (Eulerian).
spec = JacobianSpec(Structured(strategy = DirectSolve()), problem)
sol  = solve(problem; jacobian = spec)

# (c) Large ALE, fast approximate Jacobian (~8× faster assembly).
spec = JacobianSpec(Structured(strategy = ApproxSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true)

# (d) Large ALE, exact but cheap assembly via Richardson.
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true)

# (e) ILU-GMRES linear solve for a large W matrix.
sol = solve(problem;
    integrator = KenCarp4(linsolve = ILUGMRESFactorization(τ = 0.1)),
    jacobian   = JacobianSpec(Structured(strategy = DirectSolve()), problem; use_ale = true),
    use_ale    = true)

# (f) Depletion under ALE — the default Jacobian is Structured(strategy = ApproxSolve()).
sol = solve(problem;
    use_ale          = true,
    handle_depletion = true)

# (g) Vendor solver via extension.
using KLU
sol = solve(problem; linear_solver = Pyrolysis.klu_factorization())

# (h) Verify a custom backend against the dense reference during a solve.
spec = JacobianSpec(Structured(strategy = StructuredSolve()), problem; use_ale = true)
sol  = solve(problem; jacobian = spec, use_ale = true,
             verify_jacobian = (rtol = 1e-6, atol = 1e-10, sample_times = [0.0, 5.0]))
```

---

## 11.11 Cross-references

- Structured-Jacobian architecture, the implicit-RK `W = I/(γ_RK Δt) − J`
  matrix, per-operator Jacobian formulas, and the prefix/suffix/rank-1
  couplings: **Technical Reference §15**.
- Surface-depletion handling and the depletion ↔ pluggable-Jacobian
  incompatibility: **Chapter 10** and **Technical Reference §14.7**.
- Building a problem and selecting simulation modes: **Chapter 7**.
- The `solve` keyword reference and post-processing: **Chapter 8** and
  **Chapter 9**.
- Sensitivity workflows that build on the sensitivity extension: **Chapter 13**
  and **Technical Reference §17**.
- The complete API and default-value tables: **Chapter 16**.
