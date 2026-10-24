# 12. Parameter Sweeps

A *parameter sweep* runs the same base problem many times, each time with one or
more inputs changed, and reduces every run to a small summary. Pyrolysis.jl ships
a single harness for this, `parameter_sweep(problem, grid; modify, extract)`,
which builds a fresh problem at each grid point, solves it with `solve`, and
collects the result. This chapter explains how to build the grid, write the
`modify` and `extract` closures correctly (Pyrolysis types are immutable, so you
*rebuild* rather than mutate), run the sweep in parallel with `Distributed.jl`,
sweep the mesh and the solver settings, and tabulate the output. Sweeps give you
discrete maps of how outputs respond to inputs; for *derivatives* of outputs with
respect to inputs, take finite differences of whole solves (Chapter 13). For the
theory behind global/variance-based (Sobol) sensitivity built on top of this
harness, see Technical Reference §17.

> **Coordinate and symbol conventions.** `z = 0` is the back face (tag
> `:bottom`), `z = L` the exposed surface (tag `:top`); heat enters from the top.
> Arrhenius parameters are the pre-exponential `A_i` [s⁻¹], the activation energy
> `E_i` [J·mol⁻¹], and the heat of reaction `h` [J·kg⁻¹] with the **`h > 0`
> endothermic** convention (see Technical Reference §2 and §6).

---

## 12.1 The `parameter_sweep` function

The public entry point is exported from `Pyrolysis` and has this signature:

```julia
parameter_sweep(
    problem,
    grid;
    modify,                       # (base_problem, point) -> problem; required
    extract = identity,           # solution -> summary
    parallel = nworkers() > 1,    # Distributed.pmap when true, map otherwise
    solve_kwargs...,              # forwarded to every solve() call
) -> Vector
```

The harness only orchestrates; you supply the two closures that define *what
changes* (`modify`) and *what you keep* (`extract`). Its body is, in effect:

```julia
work = point -> extract(solve(modify(problem, point); solve_kwargs...))
points = collect(grid)
return parallel ? pmap(work, points) : map(work, points)
```

Key facts that follow directly from this implementation:

| Aspect | Behavior |
|---|---|
| Return value | A `Vector`, one entry per grid point, in `collect(grid)` order. Each entry is whatever `extract` returns (the full `PyrolysisSolution` by default). |
| Base problem | Never mutated. Each grid point gets its own problem from `modify`; it may share the base problem's mesh, which a solve never changes. |
| Grid materialization | `collect(grid)` runs first, so lazy iterators (e.g. `Iterators.product`) are realized into a concrete collection before dispatch. |
| `solve_kwargs` | Forwarded verbatim to each `solve()`. This is where `abstol`, `reltol`, `use_ale`, `radiation_model`, `deactivation`, `min_thickness`, `dtmin`, etc. go (Chapter 8 lists every `solve` keyword). An unknown keyword is refused by `solve` with an `ArgumentError` at the first point. |
| `progress` | Not set by the harness: `solve`'s default (`progress = false`) applies, and `progress = true` in `solve_kwargs` draws one bar per point. |
| Failures | Not retried. A `solve` that throws propagates out of `map`/`pmap`; a solve that ends without throwing (for example with return code `:Failure`) is passed to `extract` like any other (§12.7). |

`modify` has no usable default: leaving it out throws

```text
ArgumentError: parameter_sweep needs `modify = (problem, point) -> problem`, which builds the problem of each grid point
```

The default `extract` is `identity`, i.e. each entry of the returned vector is the
entire `PyrolysisSolution`. For anything beyond a handful of points pass an
`extract` that returns a small `NamedTuple` (§12.4): full solutions hold every
saved field, the integrator's solution and the solve's definition, which is a lot
of memory and expensive to ship back from workers.

---

## 12.2 Building the grid

The `grid` argument is **any iterable of parameter points**. A "point" is whatever
your `modify` closure knows how to consume (a `NamedTuple`, a `Vector`, a `Dict`,
a scalar); the harness does not interpret it. Common shapes:

**1-D sweep (a vector of scalars).** Sweep a single quantity such as the incident
heat flux:

```julia
fluxes = [15e3, 25e3, 35e3, 50e3, 75e3]   # W/m², a Vector{Float64}
```

**N-D Cartesian grid (a matrix of `NamedTuple`s).** A comprehension over two axes
produces a `Matrix{NamedTuple}`; `collect`, `map` and `pmap` traverse it
column-major:

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
back with the input that produced it (§12.8).

---

## 12.3 The `modify` closure: rebuild, do not mutate

`modify(base_problem, point)` must **return a `PyrolysisProblem`** for the given
point. Pyrolysis's `Material`, `Reaction`, component and `PyrolysisProblem` types
are immutable `struct`s with no in-place setters and **no copy constructors**:
you rebuild the pieces you want to change and reuse the rest by reading the
struct's fields. The relevant fields:

- `PyrolysisProblem`: `mesh`, `material`, `bc_set`, `T_initial`, `ξ_initial`,
  `tspan`, `dt_initial`, `output_times`.
- `Material`: `name`, `components` (a tuple), `reactions` (a tuple), `mixing`,
  `lateral_shrinkage_law`, `depletion_limiter`, `gas_viscosity`.
- `Reaction`: `name`, `reactants`, `products`, `A`, `E`, `heat`, `T_min`,
  `T_max`, `gate_width`.

The constructors you rebuild with are the standard keyword constructors (Chapters
3–5):

- `PyrolysisProblem(; mesh, material, bc_set, T_initial, ξ_initial, tspan,
  dt_initial, output_times)`.
- `Material(; name, components, reactions, mixing, lateral_shrinkage_law,
  depletion_limiter, gas_viscosity)`. Pass every field you do not change
  explicitly: a field left out takes its default (for example the air viscosity
  curve `SutherlandProperty()` for `gas_viscosity`), not the base material's value.
- `Reaction(name, reactant => product; A, E, h, n, T_min, T_max, gate_width)` and
  `Reaction(name, reactant => (p1, p2); A, E, h, yields, n, T_min, T_max,
  gate_width)`.

Reusing `prob.mesh` in every point is safe, on a moving mesh included: `solve`
works on its own copy of the mesh, so each point starts from the base geometry.
A `modify` may also build a mesh of its own (§12.6). The initial fields
`T_initial` and `ξ_initial` are functions of position that `setup` evaluates at
the cell centres of whatever mesh the problem carries.

### 12.3.1 Sweeping Arrhenius parameters (A, E)

The cleanest pattern: rebuild the single reaction you are sweeping, splice it back
into a reactions tuple, rebuild the material, and rebuild the problem. The
example below assumes the swept reaction is the first one and reuses everything
else from the base material. Note the `Reaction` fields: `reactants` is a tuple of
`(component_idx, order)` pairs and `products` a tuple of
`(component_idx, mass_stoich)` pairs (Chapter 4):

```julia
function modify_AE(prob, point)
    mat = prob.material
    old_rxn = mat.reactions[1]               # the reaction we are sweeping

    # Rebuild reaction 1 with the swept A, E; reuse its other parameters.
    reactant_idx = old_rxn.reactants[1][1]   # (component_idx, order)
    order        = old_rxn.reactants[1][2]
    product_idx  = old_rxn.products[1][1]    # (component_idx, mass_stoich)
    new_rxn = Reaction(
        old_rxn.name,
        reactant_idx => product_idx;
        A          = point.A,                # swept pre-exponential [s⁻¹]
        E          = point.E,                # swept activation energy [J/mol]
        h          = old_rxn.heat,           # reuse the heat-of-reaction property
        n          = order,
        T_min      = old_rxn.T_min,
        T_max      = old_rxn.T_max,
        gate_width = old_rxn.gate_width,
    )

    # Splice the rebuilt reaction back into the tuple. Tuples are immutable, so
    # construct another one (here: replace index 1, keep the rest):
    new_reactions = (new_rxn, Base.tail(mat.reactions)...)

    new_mat = Material(
        name                  = mat.name,
        components            = mat.components,        # reused unchanged
        reactions             = new_reactions,
        mixing                = mat.mixing,
        lateral_shrinkage_law = mat.lateral_shrinkage_law,
        depletion_limiter     = mat.depletion_limiter,
        gas_viscosity         = mat.gas_viscosity,
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
  so passing it straight back as `h` is fine: the `Reaction` constructor accepts
  either a number or a property function and runs `to_property` on it.
- If the reaction has **two products**, use the two-product `Reaction`
  constructor and pass the yields you read back from the old reaction:
  `reactant_idx => (old_rxn.products[1][1], old_rxn.products[2][1])` with
  `yields = (old_rxn.products[1][2], old_rxn.products[2][2])`.
- Mass balance is re-validated when the reaction is reconstructed
  (`validate_mass = true` by default), so a malformed point fails fast.

### 12.3.2 Sweeping the heat of reaction (h)

Heat of reaction is a property of the reaction, swept the same way. Remember the
sign convention: **`h > 0` is endothermic** (Technical Reference §6):

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
        T_min = old.T_min, T_max = old.T_max, gate_width = old.gate_width,
    )
    new_mat = Material(
        name = mat.name, components = mat.components,
        reactions = (new_rxn, Base.tail(mat.reactions)...),
        mixing = mat.mixing,
        lateral_shrinkage_law = mat.lateral_shrinkage_law,
        depletion_limiter = mat.depletion_limiter,
        gas_viscosity = mat.gas_viscosity,
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
`bc_set`. A `BoundaryConditionSet` has two `NamedTuple` fields, `thermal` and
`mass`, each keyed by the mesh's boundary tags (`:top` and `:bottom` on a line).
The example below sweeps the incident radiant flux delivered by a
`RadiativeFluxBC` at the top surface (the cone-heater flux), keeping the bottom
thermal condition and both mass conditions unchanged:

```julia
function modify_flux(prob, q_incident)            # point is a scalar flux [W/m²]
    bcs = prob.bc_set
    new_top_thermal =
        RadiativeFluxBC(flux = q_incident, absorptivity = 0.96) +
        RadiativeBC(ε = 0.96, T_inf = 330.0) +
        ConvectiveBC(h = 7.2, T_inf = 330.0)

    new_bc_set = BoundaryConditionSet(
        thermal = (top = new_top_thermal, bottom = bcs.thermal.bottom),
        mass    = bcs.mass,
    )

    return PyrolysisProblem(
        mesh = prob.mesh, material = prob.material, bc_set = new_bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end
```

The `+` operator on thermal conditions builds a `CombinedThermalBC` (Chapter 5).
Each grid point therefore gets its own combined condition with the swept incident
flux while re-radiation and convection stay fixed. The incident flux of a
`RadiativeFluxBC` enters the surface balance under `NO_RADIATION` and
`SURFACE_ABSORPTION` and is absorbed in depth under `BEER_LAMBERT`; pass the
`radiation_model` you want in `solve_kwargs` (Chapter 8, Technical Reference §8).

> **Caveat: keep `modify` self-contained for parallel runs.** When you run with
> workers (§12.5) the closure body executes on a remote process. Anything it
> references (helper functions, constants) must be defined on the workers too
> (use `@everywhere`). Reading fields off the `base_problem` is safe because the
> base problem is shipped to each worker by `pmap`.

---

## 12.4 The `extract` closure: keep results small

`extract(solution)` reduces a `PyrolysisSolution` to whatever you want to keep.
The default is `identity` (the full solution). For real sweeps, return a small
`NamedTuple` of scalars or short vectors. The fields available on a
`PyrolysisSolution` (Chapter 9) include:

| Field | Meaning |
|---|---|
| `sol.t` | Saved times [s] |
| `sol.T` | Cell temperatures, `n_cells × n_times` [K], `NaN` in a deactivated cell |
| `sol.ξ` | Concentrations, `NC × n_cells × n_times` [kg/m³], `NaN` in a deactivated cell |
| `sol.surface_T` | Surface temperature of the `:top` faces at each saved time [K] |
| `sol.mass_loss_rate` | Mass leaving through every boundary face, per unit *current* area of the `:top` faces [kg/(m²·s)] |
| `sol.retcode` | Integrator return code (`:Success`, `:Terminated`, `:Failure`, …) |
| `sol.outcome` | `:running` unless deactivation or the thickness stop ended the solve (`:burned_out`, `:interior`, `:min_thickness`) |
| `sol.statistics` | Solver counters (`accepted`, `rejected`, `domain_rejections`, `invocations`, …) |
| `sol.solve_time` | Wall-clock solve time [s] |
| `sol.extras.total_mass` | Summed component masses per unit initial `:top` area [kg/m²] |
| `sol.extras.thickness_history` | Distance from the first node to the outermost live one [m] |
| `sol.extras.MLR_gauge` | Mass-loss rate per unit *initial* `:top` area [kg/(m²·s)] |
| `sol.extras.gas_generation_rate` | Gas the reactions of the live cells produce, per unit initial `:top` area [kg/(m²·s)] |
| `sol.extras.live_cells` | Number of live cells at each saved time |

A typical scalar-summary `extract`:

```julia
extract_summary(sol) = (
    retcode    = sol.retcode,
    ok         = succeeded(sol),
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
> the solution's buffers, which then become garbage.

Always carry `sol.retcode` (or `succeeded(sol)`) into your summary so you can
detect grid points that did not finish: `succeeded(sol)` is true for `:Success`
and for a `:Terminated` solve whose outcome is a complete burnout
(`:burned_out`) or the `min_thickness` stop; `completed(sol)` is true for
`:Success` only. `:Failure` (among its causes an event that leaves the state
inadmissible or the node velocities inconsistent) and the outcome `:interior` (a
cell below the threshold inward of a live one) flag a solve that did not finish.

A per-point conservation check costs nothing extra: `conservation_audit(sol)`
closes mass and energy by the trapezoid over the saved series (Chapter 9), so an
`extract` can return `conservation_audit(sol).mass.relative` beside the
quantities of interest. A `ConservationRecorder` passed through `solve_kwargs`
would be shared by every point of the sweep, so audit through the saved series
inside `extract` instead.

---

## 12.5 Parallel execution with `Distributed.jl`

`parallel` defaults to `true` whenever more than one worker is available
(`nworkers() > 1`); otherwise the sweep runs sequentially with `map`. When
parallel, the harness uses `Distributed.pmap`, which load-balances grid points
across workers. No state is shared between points: every `solve` sets up its
own problem (the definition, the workspace, the state and the colored Jacobian
cache, Chapter 11) from the immutable problem `modify` returned, and the base
problem and its material are shipped to the workers by `pmap`.

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
                       n = old.reactants[1][2], T_min = old.T_min, T_max = old.T_max,
                       gate_width = old.gate_width)
    new_mat = Material(name = mat.name, components = mat.components,
                       reactions = (new_rxn, Base.tail(mat.reactions)...),
                       mixing = mat.mixing,
                       lateral_shrinkage_law = mat.lateral_shrinkage_law,
                       depletion_limiter = mat.depletion_limiter,
                       gas_viscosity = mat.gas_viscosity)
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
  `@everywhere`). A bare `function` defined only on the head process is not
  found on the workers and the sweep fails with an `UndefVarError`.
- `parallel` is computed once at call time from `nworkers()`. If you call
  `parameter_sweep` before `addprocs`, it runs sequentially even if you add
  workers afterward. Force the path explicitly with `parallel = true` /
  `parallel = false` to override the default.
- A parallel sweep returns the same values as the serial one: each point's solve
  is deterministic and independent of where it runs.
- There is no benefit to more workers than grid points; `pmap` leaves the extra
  workers idle.
- The order of the returned vector matches `collect(grid)`, *not* the order in
  which workers finished: `pmap` preserves input order.
- A whole `PyrolysisSolution` can be sent back from a worker (with the default
  `extract = identity`), but it is large; reduce on the worker instead.

For a sequential sweep (debugging, single core, or small grids where worker
start-up cost dominates) either omit `addprocs` or pass `parallel = false`:

```julia
results = parameter_sweep(problem, fluxes; modify = modify_flux,
                          extract = extract_summary, parallel = false)
```

Solves also run concurrently on threads: the inverse-analysis and MaCFP examples
(Chapter 14) start one `Threads.@spawn` task per scenario. `parameter_sweep`
itself distributes over processes only.

---

## 12.6 Sweeping resolution and solver settings

`modify` may change anything the problem holds, the mesh included, so a
grid-convergence study is a sweep whose points are cell counts:

```julia
const L_SAMPLE = 5.8e-3   # the thickness of the base problem's mesh [m]

function modify_cells(prob, n_cells)
    return PyrolysisProblem(
        mesh = create_uniform_mesh(L_SAMPLE, n_cells),
        material = prob.material, bc_set = prob.bc_set,
        T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
        tspan = prob.tspan, dt_initial = prob.dt_initial,
        output_times = prob.output_times,
    )
end

results = parameter_sweep(problem, [120, 60, 30, 15];
                          modify = modify_cells, extract = extract_summary,
                          radiation_model = SURFACE_ABSORPTION)
```

The fixed `output_times` make the saved series of every resolution share one time
base, so curves compare element by element (`examples/convergence_douglas_fir.jl`,
§14.3, builds its problem the same way inside its own `build_and_solve(n_cells)`).
Since `T_initial` and `ξ_initial` are functions of position, they carry over to
any mesh unchanged.

The `solve_kwargs` of a sweep are the same at every point. To vary a solver
setting (the tolerances, the integrator policy, the linear solver) per point,
map over the settings yourself:

```julia
tolerances = [1e-5, 1e-6, 1e-7, 1e-8]
results = map(tolerances) do rt
    extract_summary(solve(problem; reltol = rt, abstol = 1e-8 * rt / 1e-6,
                          radiation_model = SURFACE_ABSORPTION))
end
```

(use `pmap` with `@everywhere` definitions for a parallel version). Settings that
`integrate` takes (`integrator_policy`, `linear_solver`, `dtmax`, `saveat`, …)
can be swept on one set-up problem: `cp = setup(problem; …)` once, then
`integrate(cp; integrator_policy = …)` per setting. Every `integrate` call starts
from the set-up initial state with every cell live, so repeated integrations of
one `cp` are independent of one another (Chapter 11). The tolerances, the
radiation model, `use_ale` and deactivation are `setup`'s and need a `setup` per
value.

---

## 12.7 Fault tolerance

A `solve` that throws aborts the whole sweep (the exception propagates out of
`pmap`/`map`). Most failures do not throw: a solve that cannot continue ends with
return code `:Failure` (or `:Terminated` with outcome `:interior`) and
returns a solution, which reaches `extract`; record
`sol.retcode`, `sol.outcome` and `succeeded(sol)` there (§12.4). Errors that do
throw are invalid input: an unknown `solve` keyword, a problem that fails
validation, an inadmissible initial state, a `modify` that cannot build the
point's problem. For long grids you usually want such a point to record a marker
and let the rest continue. The harness calls `modify` and `solve` itself, so to
own the error handling completely, take over the per-point body and drive it
through your own `map`/`pmap`:

```julia
@everywhere function point_result(prob, point)
    try
        sol = solve(modify_AE(prob, point); abstol = 1e-8, reltol = 1e-6)
        return (ok = succeeded(sol), retcode = sol.retcode, outcome = sol.outcome,
                peak_T = isempty(sol.surface_T) ? NaN : maximum(sol.surface_T),
                error = "")
    catch err
        return (ok = false, retcode = :Error, outcome = :none, peak_T = NaN,
                error = sprint(showerror, err))
    end
end

using Distributed
pts = collect(grid)
results = (nworkers() > 1 ? pmap : map)(p -> point_result(problem, p), pts)
```

If you stick with `parameter_sweep` itself, remember that its `extract` runs
*after* `solve` returns, so an `extract`-side `try` catches problems in your
reduction (e.g. `maximum` on an empty vector) but **not** an exception thrown
inside `modify` or `solve`.

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
and `map`/`pmap` traverse it column-major, reshape the result back to the grid's
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

This puts the pieces together on a self-contained PMMA cone problem (the deck of
Chapter 2 and `examples/06_cone_calorimeter.jl` in condensed form, with a small
gas-transfer coefficient on the solids so that Fickian gas transport does not
degenerate, §13.9). The sweep varies the incident cone flux and records peak
surface temperature, peak mass-loss rate and final mass.

```julia
using Pyrolysis

# --- Base material: PMMA (virgin -> char + MMA gas) ---
# Piecewise fits held constant outside their ~300–700 K range.
pmma_cp = CallableProperty(
    T -> T < 378.0 ? (-2290.0 + 11.2 * max(T, 300.0)) : (1040.0 + 3.08 * T),
)
pmma_k = CallableProperty(T -> T < 378.0 ? 0.15 : (0.34 - 3.9e-4 * min(T, 700.0)))

pmma = Material(
    name = :PMMA,
    components = (
        SolidComponent(:virgin, ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96, λ = 1e-6),
        SolidComponent(:char,   ρ = 1210.0, c = pmma_cp, k = pmma_k, ε = 0.96, λ = 1e-6),
        GaseousComponent(:MMA,  M = 0.10012, c = 1500.0, k = 0.02, λ = 1e-5),
    ),
    reactions = (
        Reaction(:decomposition, 1 => (2, 3);
                 A = 1.5e14, E = 2.03e5, h = 820e3, n = 1.0,
                 yields = (0.002, 0.998)),     # h > 0 endothermic
    ),
    mixing = MaterialMixing(k = MixingSpec(WEIGHTED, 0.5), λ = MixingSpec(PARALLEL)),
)

# --- Mesh, boundary conditions and base problem ---
mesh = create_uniform_mesh(5.8e-3, 72)

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
        mass    = bcs.mass,
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
                          radiation_model = SURFACE_ABSORPTION,
                          abstol = 1e-8, reltol = 1e-6)

for (q, r) in zip(fluxes, results)
    println("q = $(q/1000) kW/m²: retcode=$(r.retcode), ",
            "peak T_s=$(round(r.peak_T, digits=1)) K, ",
            "peak MLR=$(round(r.peak_MLR*1000, digits=3)) g/(m²·s)")
end
```

This sweep runs on a fixed mesh because no `use_ale = true` is passed in
`solve_kwargs`. To sweep a moving-mesh problem instead (for example to study how
the swept flux changes the final sample thickness), add the moving-mesh and
burnout options to the sweep call and read the corresponding extras:

```julia
extract_ale(sol) = (
    retcode         = sol.retcode,
    outcome         = sol.outcome,
    live_cells      = sol.extras.live_cells[end],
    peak_T          = maximum(sol.surface_T),
    final_thickness = sol.extras.thickness_history[end],
    peak_MLR_gauge  = maximum(sol.extras.MLR_gauge),  # per initial area
)

results_ale = parameter_sweep(problem, fluxes;
    modify = modify_flux, extract = extract_ale, parallel = false,
    radiation_model = SURFACE_ABSORPTION,
    use_ale = true, deactivation = true, min_thickness = 0.01e-3,
    dtmin = 1e-14)
```

With `deactivation = true` the cells whose condensed mass falls below the
threshold (1 % of their initial condensed mass by default) are deactivated into
the live cell beneath (Chapter 10). This PMMA leaves a 0.2 % char residue, so at
the higher fluxes the layer deactivates down to one cell, which holds the residue
of the whole layer, and the solve runs to the end of its time span
(`:Success`, outcome `:running`); `live_cells` and `final_thickness` show it.

On a variable cross section and a moving mesh, use `sol.extras.MLR_gauge` rather
than `sol.mass_loss_rate` when comparing across grid points: `mass_loss_rate` is
per *current* `:top` area, whereas `MLR_gauge` is per *initial* `:top` area
`A_0`, the like-for-like quantity for cone comparisons (Chapter 9, Technical
Reference §10).

---

## 12.10 Worked example: an activation-energy / pre-exponential grid

A 2-D kinetics grid maps how peak surface temperature and the timing of peak gas
release respond to the Arrhenius pair `(A_1, E_1)`. The PMMA reaction of §12.9 has
two products, so its `modify` uses the two-product constructor:

```julia
function modify_AE2(prob, point)
    mat = prob.material
    old = mat.reactions[1]
    (p1, y1), (p2, y2) = old.products
    new_rxn = Reaction(old.name, old.reactants[1][1] => (p1, p2);
                       A = point.A, E = point.E, h = old.heat, yields = (y1, y2),
                       n = old.reactants[1][2], T_min = old.T_min, T_max = old.T_max,
                       gate_width = old.gate_width)
    new_mat = Material(name = mat.name, components = mat.components,
                       reactions = (new_rxn, Base.tail(mat.reactions)...),
                       mixing = mat.mixing,
                       lateral_shrinkage_law = mat.lateral_shrinkage_law,
                       depletion_limiter = mat.depletion_limiter,
                       gas_viscosity = mat.gas_viscosity)
    return PyrolysisProblem(mesh = prob.mesh, material = new_mat, bc_set = prob.bc_set,
                            T_initial = prob.T_initial, ξ_initial = prob.ξ_initial,
                            tspan = prob.tspan, dt_initial = prob.dt_initial,
                            output_times = prob.output_times)
end

A_vals = (5e13, 1e14, 1.5e14, 3e14)
E_vals = (1.9e5, 2.03e5, 2.15e5)
grid   = [(A = a, E = e) for a in A_vals, e in E_vals]   # 4×3 Matrix

extract_kinetics(sol) = (
    retcode = sol.retcode,
    peak_T  = maximum(sol.surface_T),
    # time of peak gas-generation rate as a crude "ignition" proxy:
    t_peak  = sol.t[argmax(sol.extras.gas_generation_rate)],
)

results = parameter_sweep(problem, grid;
                          modify = modify_AE2, extract = extract_kinetics,
                          radiation_model = SURFACE_ABSORPTION)

peak_T = reshape([r.peak_T for r in results], size(grid))
t_peak = reshape([r.t_peak for r in results], size(grid))
# peak_T[i, j], t_peak[i, j] correspond to (A = A_vals[i], E = E_vals[j])
```

Run this one in parallel for the 12 points by `addprocs`-ing first and defining
the base problem's pieces, `modify_AE2` and `extract_kinetics` with `@everywhere`
(§12.5). Because the grid is a `Matrix`, the reshape recovers the `(A, E)` layout
for heatmaps or contour plots.

---

## 12.11 From sweeps to sensitivities

A parameter sweep gives you a **discrete map** of output against input: you
sample the input space and read off outputs. When you instead need the
**derivative** of an output with respect to inputs (for gradient-based
calibration, uncertainty propagation or optimization), take finite differences of
the same solves (Chapter 13); no derivative is taken through the solve. The
`p_inject` closure that Chapter 13 perturbs plays exactly the role `modify` plays
here: it rebuilds a problem from a parameter vector. A `modify` you write for a
sweep is usually one rename away from a `p_inject`; the only difference is the
parameter shape (`p_inject` takes a numeric vector `θ`, `modify` an arbitrary
grid point). A sweep over the perturbed parameter vectors runs the
finite-difference solves in parallel (§13.5).

The sweep harness is also the substrate for **global / variance-based (Sobol)**
sensitivity: you sample the parameter space (the grid), run the sweep, and feed
the input/output table into a Sobol estimator. See
`examples/sensitivity_analysis/sensitivity_douglas_fir.jl` and
`sensitivity_sobol.jl` (§14.8) and Technical Reference §17.5 for that workflow.

---

## 12.12 Quick reference

| Item | Value / note |
|---|---|
| Function | `parameter_sweep(problem, grid; modify, extract, parallel, solve_kwargs...)` |
| `modify` | required; leaving it out throws an `ArgumentError` |
| `extract` default | `identity` (returns the full `PyrolysisSolution`) |
| `parallel` default | `nworkers() > 1` (`pmap`), otherwise `map` |
| `progress` | `solve`'s default (`false`); pass `progress = true` in `solve_kwargs` for a bar per point |
| Return | `Vector` in `collect(grid)` order |
| Grid shapes | `Vector`, comprehension `Matrix`, `Iterators.product`, `Vector{Dict}` |
| Rebuild, don't mutate | no copy constructors; reuse struct fields (all of them, `depletion_limiter` and `gas_viscosity` included) and call the keyword constructors |
| Parallel safety | each solve sets up its own definition, workspace and Jacobian cache; nothing is shared between points |
| Define closures on workers | `@everywhere` for `modify`/`extract` and their helpers |
| Per-point solver settings | your own `map`/`pmap` over `solve`, or `setup` once and `integrate` per setting (§12.6) |
| Finished? | `succeeded(sol)` (`:Success`, or burnout / `min_thickness` termination); `completed(sol)` (`:Success` only) |
| Cross-grid MLR (moving mesh, variable area) | compare `sol.extras.MLR_gauge`, not `sol.mass_loss_rate` |
| Gradients instead of maps | Chapter 13 (finite differences of `solve`) |
