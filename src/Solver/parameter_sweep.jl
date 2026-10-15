"""
Parameter-sweep harness — `pmap`-style parallel evaluation of a base
problem at every point in a parameter grid.

Designed for sensitivity studies, uncertainty propagation, and
inverse-analysis grids. The grid is an iterable of parameter sets; the
caller supplies a `modify` closure that materialises a fresh
`PyrolysisProblem` from the base problem and a parameter point, plus an
optional `extract` closure that reduces each `PyrolysisSolution` to a
small summary so the worker -> head transfer cost stays bounded.

Each worker process gets its own `Workspace` (allocated lazily inside
the `solve()` it runs), so per-worker state is isolated and parameter
points run in parallel without contention.
"""

using Distributed: pmap, nworkers

"""
	parameter_sweep(problem, grid;
	                modify::Function = default_modify,
	                extract::Function = identity,
	                parallel::Bool = nworkers() > 1,
	                solve_kwargs...) -> Vector

Sweep `problem` over `grid` and return a vector of results, one per grid
point. Each grid point `p ∈ grid` is materialised into a fresh problem
via `modify(problem, p)`, solved with `solve_kwargs`, then reduced via
`extract(solution)` before being returned.

# Arguments
- `problem::PyrolysisProblem`: The base problem. **Not** mutated.
- `grid`: Any iterable of parameter points. Common shapes:
  `Vector{<:NamedTuple}`, `Iterators.product(...)`, `Vector{Vector{Float64}}`.
- `modify`: `(base_problem, point) -> new_problem`. Must return a fresh
  `PyrolysisProblem` for each call; do **not** mutate `base_problem`. The new
  problem may reuse the base problem's mesh: `solve` works on its own copy of
  the mesh, so every point starts from the base geometry.
  The default raises a helpful error pointing to this kwarg.
- `extract`: `solution -> summary`. Defaults to `identity`. For large
  grids prefer a small `NamedTuple` — full `PyrolysisSolution` instances
  are expensive to ship across worker boundaries.
- `parallel`: Use `Distributed.pmap` when `true`, otherwise fall back to
  a serial `map`. Defaults to parallel when more than one worker is
  available (`nworkers() > 1`).
- `solve_kwargs`: Forwarded to `solve()`.

# Notes
- Use `Distributed.addprocs()` (and `@everywhere using Pyrolysis`) before
  calling this function to make workers available.
- Failed solves are not retried. Wrap the body of `extract` in a `try`
  if you need fault tolerance.

# Example

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13, 1e14), e in (1.5e5, 2.0e5)]

function modify(prob, p)
	# Rebuild the material with the swept Arrhenius parameters.
	mat = prob.material
	new_rxn = Reaction(:decomp, 1 => 2; A = p.A, E = p.E, h = 200e3)
	new_mat = Material(name = mat.name, components = mat.components,
	                   reactions = (new_rxn,), mixing = mat.mixing,
	                   lateral_shrinkage_law = mat.lateral_shrinkage_law)
	return PyrolysisProblem(mesh = prob.mesh, material = new_mat,
	                        bc_set = prob.bc_set, T_initial = prob.T_initial,
	                        ξ_initial = prob.ξ_initial, tspan = prob.tspan,
	                        dt_initial = prob.dt_initial,
	                        output_times = prob.output_times)
end

extract(sol) = (peak_T = maximum(sol.surface_T),
                final_mass = sol.extras.total_mass[end])

results = parameter_sweep(problem, grid; modify, extract)
```
"""
function parameter_sweep(
    problem,
    grid;
    modify::Function = _default_sweep_modify,
    extract::Function = identity,
    parallel::Bool = nworkers() > 1,
    solve_kwargs...,
)

    work = function (point)
        new_problem = modify(problem, point)
        sol = solve(new_problem; progress = false, solve_kwargs...)
        return extract(sol)
    end

    # Materialise iterators (e.g. `Iterators.product`) before dispatch so
    # pmap sees a concrete collection. `collect` over `Vector` is a no-op.
    pts = collect(grid)
    return parallel ? pmap(work, pts) : map(work, pts)
end

function _default_sweep_modify(_, _)
    throw(
        ArgumentError(
            "parameter_sweep requires a `modify::(problem, point) -> problem` " *
            "closure that materialises a fresh PyrolysisProblem from each grid " *
            "point. See the docstring for an example.",
        ),
    )
end
