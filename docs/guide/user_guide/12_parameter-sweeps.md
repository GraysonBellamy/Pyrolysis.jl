# 12. Parameter Sweeps

A *parameter sweep* runs the same base problem many times, each time with one or
more inputs changed, and reduces every run to a small summary. Pyrolysis.jl ships
a single harness for this — `parameter_sweep(problem, grid; modify, extract)` —
that materializes a fresh problem at each grid point, solves it, and collects the
result. This chapter explains how to build the grid, write the `modify` and
`extract` closures correctly (Pyrolysis types are immutable, so you *rebuild*
rather than mutate), run the sweep in parallel with `Distributed.jl`, reuse a
warm start between grid points, and tabulate the output. Sweeps give you discrete
maps of how outputs respond to inputs; for *derivatives* of outputs with respect
to inputs (gradients), use the sensitivity machinery in Chapter 13 instead. For
the theory behind global/variance-based (Sobol) sensitivity built on top of this
harness, see Technical Reference §17.

> **Coordinate and symbol conventions.** `z = 0` is the bottom/substrate
> (boundary id 2), `z = L` is the top/exposed surface (boundary id 1); heat
> enters from the top. Arrhenius parameters are the pre-exponential `A_i` [s⁻¹],
> the activation energy `E_i` [J·mol⁻¹], and the heat of reaction `h` [J·kg⁻¹]
> with the **`h > 0` endothermic** convention (see Technical Reference §2 and §6).

---

## 12.1 The `parameter_sweep` function

The public entry point is exported from `Pyrolysis` and has this signature:

```julia
parameter_sweep(
    problem,
    grid;
    modify::Function = _default_sweep_modify,   # (base_problem, point) -> new_problem
    extract::Function = identity,               # solution -> summary
    parallel::Bool = nworkers() > 1,            # use Distributed.pmap when true
    solve_kwargs...,                            # forwarded to every solve() call
) -> Vector
```

The body is intentionally tiny — the harness only orchestrates; you supply the
two closures that define *what changes* (`modify`) and *what you keep*
(`extract`):

```julia
work = function (point)
    new_problem = modify(problem, point)
    sol = solve(new_problem; progress = false, solve_kwargs...)
    return extract(sol)
end

pts = collect(grid)
return parallel ? pmap(work, pts) : map(work, pts)
```

Key facts that follow directly from this implementation:

| Aspect | Behavior |
|---|---|
| Return value | A `Vector`, one entry per grid point, in `collect(grid)` order. Each entry is whatever `extract` returns (the full `PyrolysisSolution` by default). |
| `progress` | Always forced to `false` inside the sweep (a progress bar per grid point would be noise); you cannot re-enable it via `solve_kwargs`. |
| Base problem | **Never mutated.** Each grid point gets a fresh problem from `modify`. |
| Grid materialization | `collect(grid)` runs first, so lazy iterators (e.g. `Iterators.product`) are realized into a concrete collection before dispatch. |
| Failures | Not retried. A `solve` that errors propagates out of `pmap`/`map` unless you guard inside `extract` (see §12.7). |
| `solve_kwargs` | Forwarded verbatim to each `solve()`. This is where `abstol`, `reltol`, `use_ale`, `radiation_model`, `dtmin`, etc. go (see Chapter 8 for the full `solve` option list). |

The default `modify` is a guard that throws a helpful `ArgumentError` — you must
supply your own:

```julia
function _default_sweep_modify(_, _)
    throw(ArgumentError(
        "parameter_sweep requires a `modify::(problem, point) -> problem` " *
        "closure that materialises a fresh PyrolysisProblem from each grid " *
        "point. See the docstring for an example."))
end
```

The default `extract` is `identity`, i.e. each entry of the returned vector is the
entire `PyrolysisSolution`. For anything beyond a handful of points you should
pass an `extract` that returns a small `NamedTuple` (see §12.4) — full solutions
are expensive to ship across worker boundaries and hold a lot of memory.

---

## 12.2 Building the grid

The `grid` argument is **any iterable of parameter points**. A "point" is whatever
your `modify` closure knows how to consume — a `NamedTuple`, a `Vector`, a `Dict`,
a scalar — the harness does not interpret it. Common shapes:

**1-D sweep (a vector of scalars).** Sweep a single quantity such as the incident
heat flux:

```julia
fluxes = [15e3, 25e3, 35e3, 50e3, 75e3]   # W/m², a Vector{Float64}
```

**N-D Cartesian grid (a vector of `NamedTuple`s).** A comprehension over two axes
produces a `Matrix{NamedTuple}`; `collect`/`pmap` flatten it column-major:

```julia
grid = [(A = a, E = e)
        for a in (1e12, 1e13, 1e14),   # rows
            e in (1.5e5, 2.0e5)]        # columns
# grid isa Matrix{@NamedTuple{A::Float64, E::Float64}}, size (3, 2)
```

**Lazy product (materialized by the harness).** For higher dimensions a
comprehension allocates the whole array up front; `Iterators.product` defers that
and is `collect`-ed internally:

```julia
grid = Iterators.product((1e12, 1e13, 1e14), (1.5e5, 2.0e5))
# each point is a Tuple (a, e); your modify must index it as point[1], point[2]
```

**Named axes via `Dict`.** Useful when `modify` looks up parameters by name and
not every point varies the same set:

```julia
grid = [Dict(:A => 1e13, :E => 1.8e5),
        Dict(:A => 5e13, :E => 1.9e5)]
```

Pick the shape that makes `modify` simplest to write. Because the result vector is
in `collect(grid)` order, keep the grid object around so you can pair each result
back with the input that produced it (see §12.8).

---

## 12.3 The `modify` closure — rebuild, do not mutate

`modify(base_problem, point)` must **return a fresh `PyrolysisProblem`** for the
given point. Pyrolysis's `Material`, `Reaction`, component, and `PyrolysisProblem`
types are immutable `struct`s with no in-place setters and **no copy
constructors** — you rebuild the pieces you want to change and reuse the rest by
reading the struct's fields. The relevant fields:

- `PyrolysisProblem`: `mesh`, `material`, `bc_set`, `T_initial`, `ξ_initial`,
  `tspan`, `dt_initial`, `output_times`.
- `Material`: `name`, `components` (a tuple), `reactions` (a tuple), `mixing`,
  `lateral_shrinkage_law`.
- `Reaction`: `name`, `reactants`, `products`, `A`, `E`, `heat`, `T_min`, `T_max`.

The constructors you rebuild with are the standard keyword constructors (Chapters
3–5):

- `PyrolysisProblem(; mesh, material, bc_set, T_initial, ξ_initial, tspan,
  dt_initial, output_times)`.
- `Material(; name, components, reactions, mixing, lateral_shrinkage_law)`.
- `Reaction(name, reactant => product; A, E, h, n, T_min, T_max)` and
  `Reaction(name, reactant => (p1, p2); A, E, h, yields, n, T_min, T_max)`.

### 12.3.1 Sweeping Arrhenius parameters (A, E)

The cleanest pattern: rebuild the single reaction you are sweeping, splice it back
into a fresh reactions tuple, rebuild the material, and rebuild the problem. The
example below assumes the swept reaction is the first one and reuses everything
else from the base material. Note the `Reaction` fields: `reactants` is a tuple of
`(component_idx, order)` pairs and `products` is a tuple of
`(component_idx, mass_stoich)` pairs (see Chapter 4):

```julia
function modify_AE(prob, point)
    mat = prob.material
    old_rxn = mat.reactions[1]               # the reaction we are sweeping

    # Rebuild reaction 1 with the swept A, E; reuse its other parameters.
    reactant_idx = old_rxn.reactants[1][1]   # (component_idx, order)
    order        = old_rxn.reactants[1][2]
    product_idx  = old_rxn.products[1][1]     # (component_idx, mass_stoich)
    new_rxn = Reaction(
        old_rxn.name,
        reactant_idx => product_idx;
        A     = point.A,                     # swept pre-exponential [s⁻¹]
        E     = point.E,                     # swept activation energy [J/mol]
        h     = old_rxn.heat,                # reuse heat-of-reaction property
        n     = order,
        T_min = old_rxn.T_min,
        T_max = old_rxn.T_max,
    )

    # Splice the rebuilt reaction back into the tuple. Tuples are immutable, so
    # construct a new one (here: replace index 1, keep the rest):
    new_reactions = (new_rxn, mat.reactions[2:end]...)

    new_mat = Material(
        name                  = mat.name,
        components            = mat.components,        # reused unchanged
        reactions             = new_reactions,
        mixing                = mat.mixing,
        lateral_shrinkage_law = mat.lateral_shrinkage_law,
    )

    return PyrolysisProblem(
        mesh         = prob.mesh,
        material     = new_mat,
        bc_set       = prob.bc_set,
        T_initial    = prob.T_initial,
        ξ_initial    = prob.ξ_initial,
        tspan        = prob.tspan,
        dt_initial   = prob.dt_initial,
        output_times = prob.output_times,
    )
end
```

A few notes:

- `old_rxn.heat` is already a property function (an `AbstractPropertyFunction`),
  so passing it straight back as `h` is fine — the `Reaction` constructor accepts
  either a number or a property function and runs `to_property` on it internally.
- If the reaction has **two products**, use the two-product `Reaction`
  constructor and pass the yields you read back from the old reaction:
  `reactant_idx => (old_rxn.products[1][1], old_rxn.products[2][1])` with
  `yields = (old_rxn.products[1][2], old_rxn.products[2][2])`.
- Mass balance is re-validated when the reaction is reconstructed
  (`validate_mass = true` by default), so a malformed point fails fast.

### 12.3.2 Sweeping the heat of reaction (h)

Heat of reaction is a property of the reaction, swept the same way. Remember the
sign convention: **`h > 0` is endothermic** (see Technical Reference §6):

```julia
function modify_h(prob, point)
    mat = prob.material
    old = mat.reactions[1]
    new_rxn = Reaction(
        old.name,
        old.reactants[1][1] => old.products[1][1];
        A = old.A, E = old.E,
        h = point.h,                         # [J/kg], h > 0 endothermic
        n = old.reactants[1][2],
        T_min = old.T_min, T_max = old.T_max,
    )
    new_mat = Material(
        name = mat.name, components = mat.components,
        reactions = (new_rxn, mat.reactions[2:end]...),
        mixing = mat.mixing,
        lateral_shrinkage_law = mat.lateral_shrinkage_law,
    )
    return PyrolysisProblem(
        mesh = prob.mesh, material = new_mat, bc_set = prob.bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end
```

### 12.3.3 Sweeping a boundary condition (incident heat flux)

To sweep an *applied* quantity rather than a material property, rebuild only the
`bc_set`. A `BoundaryConditionSet` exposes two `NamedTuple` fields, `thermal` and
`mass`, each keyed by `:top` and `:bottom`. The example below sweeps the incident
radiant flux delivered by a `RadiativeFluxBC` at the top surface (the cone-heater
flux), keeping the bottom thermal BC and both mass BCs unchanged:

```julia
function modify_flux(prob, q_incident)            # point is a scalar flux [W/m²]
    bcs = prob.bc_set
    new_top_thermal =
        RadiativeFluxBC(flux = q_incident, absorptivity = 0.96) +
        RadiativeBC(ε = 0.96, T_inf = 330.0) +
        ConvectiveBC(h = 7.2, T_inf = 330.0)

    new_bc_set = BoundaryConditionSet(
        thermal = (top = new_top_thermal, bottom = bcs.thermal.bottom),
        mass    = (top = bcs.mass.top,    bottom = bcs.mass.bottom),
    )

    return PyrolysisProblem(
        mesh = prob.mesh, material = prob.material, bc_set = new_bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end
```

The `+` operator on thermal BCs builds a `CombinedThermalBC` (Chapter 5). Each
grid point therefore gets its own combined boundary condition with the swept
incident flux while re-radiation and convection stay fixed. Because the swept
quantity is the incident flux of a `RadiativeFluxBC`, run the sweep with a
matching `radiation_model` in `solve_kwargs` (e.g. `SURFACE_ABSORPTION`); see
Chapter 8 and Technical Reference §8 for how the radiation model interprets that
flux.

> **Caveat — keep `modify` self-contained for parallel runs.** When you run with
> workers (§12.5) the closure body executes on a remote process. Anything it
> references — helper functions, constants — must be defined on the workers too
> (use `@everywhere`). Reading fields off the `base_problem` is safe because the
> base problem is shipped to each worker by `pmap`.

---

## 12.4 The `extract` closure — keep results small

`extract(solution)` reduces a `PyrolysisSolution` to whatever you want to keep.
The default is `identity` (the full solution). For real sweeps, return a small
`NamedTuple` of scalars or short vectors. The fields available on a
`PyrolysisSolution` (Chapter 9) include:

| Field | Meaning |
|---|---|
| `sol.t` | Saved times [s] |
| `sol.T` | Temperature field, `(n_cells, n_steps)` [K] |
| `sol.ξ` | Concentrations, `(NC, n_cells, n_steps)` [kg/m³] |
| `sol.surface_T` | Surface temperature at each saved time [K] |
| `sol.mass_loss_rate` | MLR per *current* area `A(t)` (top + permeable-bottom) [kg/(m²·s)] |
| `sol.retcode` | Integrator return code (`:Success`, `:Terminated`, …) |
| `sol.solve_time` | Wall-clock solve time [s] |
| `sol.extras.total_mass` | Domain mass per initial area [kg/m²] |
| `sol.extras.thickness_history` | Sample thickness over time [m] |
| `sol.extras.MLR_gauge` | Gauge MLR (`MLR_local · A(t)/A_0`) [kg/(m²·s)] |
| `sol.extras.gas_generation_rate` | Volumetric gas production rate [kg/(m²·s)] |

A typical scalar-summary `extract`:

```julia
extract_summary(sol) = (
    retcode    = sol.retcode,
    peak_T     = isempty(sol.surface_T)      ? NaN : maximum(sol.surface_T),
    peak_MLR   = isempty(sol.mass_loss_rate) ? NaN : maximum(sol.mass_loss_rate),
    final_mass = isempty(sol.extras.total_mass) ? NaN : sol.extras.total_mass[end],
    solve_time = sol.solve_time,
)
```

If you need a *curve* (say the surface-temperature trace) but not the full field,
return just that vector plus its time base:

```julia
extract_curve(sol) = (t = copy(sol.t), surface_T = copy(sol.surface_T))
```

> Use `copy` when returning arrays from `extract` so the summary does not alias
> the (about-to-be-garbage-collected) solution's internal buffers.

Always carry `sol.retcode` into your summary so you can detect grid points that
did not finish cleanly — `:Success` means the integrator reached `t_end`,
`:Terminated` means a callback stopped it (e.g. ALE thickness termination), and
codes like `:MaxIters` or `:DtLessThanMin` flag a failed solve.

---

## 12.5 Parallel execution with `Distributed.jl`

`parallel` defaults to `true` whenever more than one worker is available
(`nworkers() > 1`); otherwise the sweep runs sequentially with `map`. When
parallel, the harness uses `Distributed.pmap`, which load-balances grid points
across workers. The thread-safety model is built into the design:

- The pure residual `residual!(du, u, def, t, ws)` is **immutable** and reused;
  the mutable `Workspace` is allocated lazily inside each worker's `solve()`, so
  each grid point's solve has its own scratch state — no contention.
- The immutable `ProblemDef`/`Material` data is shipped to workers by `pmap`.

A complete parallel sweep:

```julia
using Distributed
addprocs(4)                       # start 4 worker processes
@everywhere using Pyrolysis       # load the package on every worker

# Make the closures (and any helpers/constants they use) available on workers:
@everywhere function modify_AE(prob, point)
    mat = prob.material
    old = mat.reactions[1]
    new_rxn = Reaction(old.name, old.reactants[1][1] => old.products[1][1];
                       A = point.A, E = point.E, h = old.heat,
                       n = old.reactants[1][2], T_min = old.T_min, T_max = old.T_max)
    new_mat = Material(name = mat.name, components = mat.components,
                       reactions = (new_rxn, mat.reactions[2:end]...),
                       mixing = mat.mixing,
                       lateral_shrinkage_law = mat.lateral_shrinkage_law)
    return PyrolysisProblem(mesh = prob.mesh, material = new_mat, bc_set = prob.bc_set,
                            T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
                            tspan = prob.tspan, dt_initial = prob.dt_initial,
                            output_times = prob.output_times)
end

@everywhere extract_summary(sol) = (
    retcode = sol.retcode,
    peak_T  = isempty(sol.surface_T) ? NaN : maximum(sol.surface_T),
)

grid = [(A = a, E = e) for a in (1e12, 1e13, 1e14), e in (1.5e5, 2.0e5)]

results = parameter_sweep(
    problem, grid;
    modify  = modify_AE,
    extract = extract_summary,
    abstol  = 1e-8, reltol = 1e-6,    # solve_kwargs forwarded to every solve()
)
```

Notes on the parallel path:

- Define `modify` and `extract` with `@everywhere` (or in a module loaded
  `@everywhere`). A bare `function` defined only on the head process will not be
  found on the workers and the sweep will fail with an `UndefVarError`.
- `parallel` is computed once at call time from `nworkers()`. If you call
  `parameter_sweep` before `addprocs`, it runs sequentially even if you add
  workers afterward. Force the path explicitly with `parallel = true` /
  `parallel = false` if you want to override the default.
- There is no benefit to more workers than grid points; `pmap` simply leaves the
  extra workers idle.
- The order of the returned vector matches `collect(grid)`, *not* the order in
  which workers finished — `pmap` preserves input order.

For a sequential sweep (debugging, single core, or small grids where worker
start-up cost dominates) either omit `addprocs` or pass `parallel = false`:

```julia
results = parameter_sweep(problem, fluxes; modify = modify_flux,
                          extract = extract_summary, parallel = false)
```

---

## 12.6 Warm starting across grid points

`solve` accepts `warm_start::PyrolysisSolution`, which seeds the integrator's
initial state with a previous solution's *final* state. On a sweep whose points
form a smooth continuation (e.g. a monotone heat-flux ramp), warm starting can cut
per-point solve time by giving the implicit integrator a better starting guess.

The compatibility rules are strict and enforced by a state-vector length check:
the component count `NC`, the mesh size, the ALE mode (`use_ale`), and the χ mode
must all be identical between the prior solve and the new one — otherwise `solve`
errors. (`tspan[1]` still sets the integrator's start time; warm starting just
reuses the prior final state as fresh initial data.) The prior solution must also
still carry its `raw_solution`, so warm-starting from a heavily reduced `extract`
output will not work — keep the full `PyrolysisSolution` for the run you chain
from.

Because `pmap` runs grid points concurrently and out of any meaningful order,
warm starting is **not compatible with `parallel = true`** — there is no
"previous" point on a parallel worker. Do warm-started sweeps **sequentially**,
chaining each solve to the last:

```julia
function warm_chained_sweep(problem, points; modify, extract, solve_kwargs...)
    results = Vector{Any}(undef, length(points))
    prev = nothing
    for (i, pt) in enumerate(points)
        prob_i = modify(problem, pt)
        sol = prev === nothing ?
              solve(prob_i; progress = false, solve_kwargs...) :
              solve(prob_i; progress = false, warm_start = prev, solve_kwargs...)
        results[i] = extract(sol)
        prev = sol      # chain the FULL solution into the next point
    end
    return results
end

# Sweep a heat-flux ramp, warm-starting each point from the previous one:
fluxes = collect(15e3:5e3:75e3)
results = warm_chained_sweep(problem, fluxes;
                             modify = modify_flux, extract = extract_summary,
                             abstol = 1e-8, reltol = 1e-6)
```

> **Accuracy caveat.** Warm starting changes only the *initial guess*, not the
> equations, so the converged trajectory of each point is unaffected to integrator
> tolerance. But it does couple grid points: a point solved from a warm start may
> take a slightly different step sequence than the same point solved cold. For
> sweeps feeding a sensitivity or calibration study where reproducibility of each
> point in isolation matters, prefer cold solves. See the warm-start discussion in
> Technical Reference §17 (sensitivity pitfalls).

---

## 12.7 Fault tolerance

A `solve` that throws aborts the whole sweep (the exception propagates out of
`pmap`/`map`). For long grids you usually want a failed point to record a marker
and let the rest continue. The harness always calls `solve` itself, so to own the
error handling completely, take over the per-point body: make `modify` a
pass-through and do the rebuild, solve, and reduction inside one `extract`-side
function guarded by `try`/`catch`. The cleanest way is to fold the work into
`extract` while letting `modify` build the problem normally and catching only the
solve inside a small wrapper of your own:

```julia
@everywhere function point_result(prob, point)
    try
        new_prob = modify_AE(prob, point)
        sol = solve(new_prob; progress = false, abstol = 1e-8, reltol = 1e-6)
        return (ok = sol.retcode === :Success, retcode = sol.retcode,
                peak_T = isempty(sol.surface_T) ? NaN : maximum(sol.surface_T))
    catch err
        return (ok = false, retcode = :Error, peak_T = NaN, error = string(err))
    end
end

# Drive it through your own map/pmap rather than the harness, so the try/catch
# wraps the solve:
using Distributed
pts = collect(grid)
results = (nworkers() > 1 ? pmap : map)(p -> point_result(problem, p), pts)
```

If you stick with `parameter_sweep` itself, remember its `extract` runs *after*
`solve` returns — so an `extract`-side `try` catches problems in your reduction
(e.g. `maximum` on an empty vector), but **not** an exception thrown inside
`solve`. To catch solve-time exceptions you must move the `solve` into your own
loop, as above.

---

## 12.8 Collecting and tabulating results

The result is a flat `Vector` in `collect(grid)` order. To turn it into something
you can plot or write out, zip it back against the grid.

**1-D sweep → paired vectors:**

```julia
fluxes  = [15e3, 25e3, 35e3, 50e3, 75e3]
results = parameter_sweep(problem, fluxes; modify = modify_flux,
                          extract = extract_summary, parallel = false)

peak_T = [r.peak_T for r in results]
for (q, pT) in zip(fluxes, peak_T)
    println("q = $(q/1000) kW/m²  ->  peak T_s = $(round(pT, digits=1)) K")
end
```

**2-D grid → reshape into a matrix.** Because a comprehension grid is a `Matrix`
and `pmap`/`map` flatten column-major, reshape the result back to the grid's
shape:

```julia
A_vals = (1e12, 1e13, 1e14)
E_vals = (1.5e5, 2.0e5)
grid   = [(A = a, E = e) for a in A_vals, e in E_vals]   # 3×2 Matrix

results = parameter_sweep(problem, grid; modify = modify_AE,
                          extract = extract_summary)

peak_T = reshape([r.peak_T for r in results], size(grid))  # 3×2 again
# peak_T[i, j] corresponds to (A = A_vals[i], E = E_vals[j])
```

**Write to a delimited file** for external post-processing:

```julia
using DelimitedFiles
open("sweep_results.csv", "w") do io
    println(io, "A,E,peak_T,final_mass,retcode")
    for (p, r) in zip(vec(grid), results)
        println(io, "$(p.A),$(p.E),$(r.peak_T),$(r.final_mass),$(r.retcode)")
    end
end
```

(If you prefer a `DataFrame`, build one from the same `zip` of grid points and
result `NamedTuple`s; that requires the `DataFrames.jl` package, which is not a
dependency of Pyrolysis.jl.)

---

## 12.9 Worked example: heat-flux series on a cone-calorimeter problem

This puts the pieces together on a self-contained PMMA cone problem (the setup is
the Chapter 2 / example 06 problem in condensed form). The sweep varies the
incident cone flux and records peak surface temperature, peak MLR, and final mass.

```julia
using Pyrolysis

# --- Base material: PMMA (virgin -> char + MMA gas) ---
pmma_cp = CallableProperty(T -> T < 378.0 ? (-2290.0 + 11.2T) : (1040.0 + 3.08T))
pmma_k  = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4T))

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96),
        SolidComponent(:char,   ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96),
        GaseousComponent(:MMA,  M = 0.10012, c = 1500.0, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:decomposition, 1 => (2, 3);
                 A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0,
                 yields = (0.002, 0.998)),     # h > 0 endothermic
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
)

# --- Mesh, BCs, and base problem ---
L, n_cells = 5.8e-3, 72
mesh = create_uniform_mesh(L, n_cells, Val(3);
                           T_initial = 300.0, ξ_initial = (1210.0, 0.0, 0.0))

base_bc = BoundaryConditionSet(
    thermal = (
        top = RadiativeFluxBC(flux = 25e3, absorptivity = 0.96) +
              RadiativeBC(ε = 0.96, T_inf = 330.0) +
              ConvectiveBC(h = 7.2, T_inf = 330.0),
        bottom = RadiativeBC(ε = 0.95, T_inf = 307.0) +
                 ConvectiveBC(h = 4.0, T_inf = 307.0),
    ),
    mass = (
        top    = ConvectiveMassBC(h_m = 0.05, ξ_inf = 0.0),
        bottom = ImpermeableBC(),
    ),
)

problem = PyrolysisProblem(
    mesh = mesh, material = pmma, bc_set = base_bc,
    T_initial = z -> 300.0,
    ξ_initial = [z -> 1210.0, z -> 0.0, z -> 0.0],
    tspan = (0.0, 600.0),
    output_times = collect(range(0.0, 600.0, length = 121)),
)

# --- The sweep: incident cone flux, 15..75 kW/m² ---
function modify_flux(prob, q_incident)
    bcs = prob.bc_set
    new_top =
        RadiativeFluxBC(flux = q_incident, absorptivity = 0.96) +
        RadiativeBC(ε = 0.96, T_inf = 330.0) +
        ConvectiveBC(h = 7.2, T_inf = 330.0)
    new_bc_set = BoundaryConditionSet(
        thermal = (top = new_top, bottom = bcs.thermal.bottom),
        mass    = (top = bcs.mass.top, bottom = bcs.mass.bottom),
    )
    return PyrolysisProblem(
        mesh = prob.mesh, material = prob.material, bc_set = new_bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end

extract_summary(sol) = (
    retcode    = sol.retcode,
    peak_T     = maximum(sol.surface_T),
    peak_MLR   = maximum(sol.mass_loss_rate),
    final_mass = sol.extras.total_mass[end],
)

fluxes  = [15e3, 25e3, 35e3, 50e3, 75e3]
results = parameter_sweep(problem, fluxes;
                          modify = modify_flux, extract = extract_summary,
                          parallel = false,
                          radiation_model = SURFACE_ABSORPTION,  # interprets the cone flux
                          abstol = 1e-8, reltol = 1e-6)

for (q, r) in zip(fluxes, results)
    println("q = $(q/1000) kW/m²: retcode=$(r.retcode), ",
            "peak T_s=$(round(r.peak_T, digits=1)) K, ",
            "peak MLR=$(round(r.peak_MLR*1000, digits=3)) g/(m²·s)")
end
```

This sweep runs in fixed-mesh (Eulerian) mode because no `use_ale = true` is
passed in `solve_kwargs`. To sweep an ALE problem instead — for example to study
how the swept flux changes the final sample thickness — add the ALE options to
the sweep call and read the ALE-specific extras:

```julia
extract_ale(sol) = (
    retcode         = sol.retcode,
    peak_T          = maximum(sol.surface_T),
    final_thickness = sol.extras.thickness_history[end],
    peak_MLR_gauge  = maximum(sol.extras.MLR_gauge),  # gauge-normalized MLR
)

results_ale = parameter_sweep(problem, fluxes;
    modify = modify_flux, extract = extract_ale, parallel = false,
    radiation_model = SURFACE_ABSORPTION,
    use_ale = true, handle_depletion = true,
    min_thickness = 0.01e-3, depletion_threshold = 0.05, min_cells = 1,
    dtmin = 1e-14)
```

In variable-cross-section/ALE runs use `sol.extras.MLR_gauge` rather than
`sol.mass_loss_rate` when comparing across grid points: `mass_loss_rate` is per
*current* area `A(t)`, whereas `MLR_gauge = MLR_local · A(t)/A_0` is normalized to
the initial gauge area `A_0`, which is the apples-to-apples quantity for cone
comparisons (see Chapter 8 and Technical Reference §10).

---

## 12.10 Worked example: an activation-energy / pre-exponential grid

A 2-D kinetics grid maps how peak surface temperature and the timing of peak gas
release respond to the Arrhenius pair `(A_1, E_1)`. Reusing `modify_AE` from
§12.3.1:

```julia
A_vals = (5e12, 1e13, 5e13, 1e14)
E_vals = (1.7e5, 1.9e5, 2.1e5)
grid   = [(A = a, E = e) for a in A_vals, e in E_vals]   # 4×3 Matrix

extract_kinetics(sol) = (
    retcode = sol.retcode,
    peak_T  = maximum(sol.surface_T),
    # time of peak gas-generation rate as a crude "ignition/burnout" proxy:
    t_peak  = sol.t[argmax(sol.extras.gas_generation_rate)],
)

results = parameter_sweep(problem, grid;
                          modify = modify_AE, extract = extract_kinetics,
                          radiation_model = SURFACE_ABSORPTION)

peak_T = reshape([r.peak_T for r in results], size(grid))
t_peak = reshape([r.t_peak for r in results], size(grid))
# peak_T[i, j], t_peak[i, j] correspond to (A = A_vals[i], E = E_vals[j])
```

Run this one in parallel for the 12 points by `addprocs`-ing first and defining
`modify_AE`/`extract_kinetics` with `@everywhere` (§12.5). Because the grid is a
`Matrix`, the reshape recovers the `(A, E)` layout for heatmaps or contour plots.

---

## 12.11 From sweeps to sensitivities

A parameter sweep gives you a **discrete map** of output-vs-input: you sample the
input space and read off outputs. When you instead need the **derivative** of an
output with respect to inputs — for gradient-based calibration, uncertainty
propagation, or optimization — use the sensitivity entry points covered in
Chapter 13:

- `forward_sensitivity(base_problem, θ, p_inject; output_fn)` — ForwardDiff
  Jacobian of an output with respect to a small parameter vector `θ`.
- `adjoint_sensitivity(base_problem, θ, p_inject, loss_fn; sensealg)` — Zygote
  gradient of a scalar loss for large parameter counts / inverse problems.

These require loading the `PyrolysisSciMLSensitivityExt` extension
(`using SciMLSensitivity, Zygote`); without it the exported names are stubs that
error with a hint to load the extension. The `p_inject` closure plays exactly the
role `modify` plays here — it rebuilds a fresh problem from a parameter vector —
so a `modify` you write for a sweep is usually one rename away from a `p_inject`
you can hand to `forward_sensitivity` (the only difference is the parameter shape:
`p_inject` takes a numeric vector `θ`, `modify` takes an arbitrary grid point).

The sweep harness is also the substrate for **global / variance-based (Sobol)**
sensitivity: you sample the parameter space (the grid), run the sweep, and feed
the input/output table into a Sobol estimator. See the
`examples/sensitivity_analysis/sensitivity_sobol.jl` study and Technical Reference
§17.5 for that workflow.

---

## 12.12 Quick reference

| Item | Value / note |
|---|---|
| Function | `parameter_sweep(problem, grid; modify, extract, parallel, solve_kwargs...)` |
| `modify` default | guard that throws `ArgumentError` — you must supply one |
| `extract` default | `identity` (returns the full `PyrolysisSolution`) |
| `parallel` default | `nworkers() > 1` |
| `progress` | forced `false` inside the sweep |
| Return | `Vector` in `collect(grid)` order |
| Grid shapes | `Vector{<:NamedTuple}`, `Iterators.product`, `Vector`/`Dict`, comprehension `Matrix` |
| Rebuild, don't mutate | no copy constructors; reuse struct fields and call keyword constructors |
| Parallel safety | immutable residual reused; mutable `Workspace` allocated per worker `solve` |
| Define closures on workers | `@everywhere` for `modify`/`extract` and their helpers |
| Warm start | sequential only; `NC` / mesh / ALE / χ mode must match; prior `raw_solution` required (see §12.6) |
| Cross-grid MLR (ALE) | compare `sol.extras.MLR_gauge`, not `sol.mass_loss_rate` |
| Gradients instead of maps | Chapter 13 (`forward_sensitivity` / `adjoint_sensitivity`) |
