# ============================================================================
# Parameter sweeps
# ============================================================================

"""
	parameter_sweep(problem, grid; modify, extract = identity,
	                parallel = nworkers() > 1, solve_kwargs...) -> Vector

Solve `problem` at every point of `grid` and return one result per point, in
the grid's order: each point `p` becomes a problem by `modify(problem, p)`,
which is solved with `solve_kwargs` ([`solve`](@ref)) and reduced by
`extract(solution)`.

- `grid`: any iterable of parameter points (a vector of `NamedTuple`s,
  `Iterators.product(...)`, …); it is collected before the sweep.
- `modify`: `(problem, point) -> problem`, returning a new problem and leaving
  the base problem alone; the new problem may share the base problem's mesh,
  which a solve never changes. Required.
- `extract`: `solution -> result`. With workers, prefer a small result (a
  `NamedTuple` of the quantities wanted): a whole solution is costly to send
  back.
- `parallel`: `Distributed.pmap` over the workers, by default whenever there
  is more than one; otherwise a serial `map`. The workers need the package
  and the definitions `modify` and `extract` use (`@everywhere`).

# Example

```julia
using Distributed
addprocs(4)
@everywhere using Pyrolysis

grid = [(A = a, E = e) for a in (1e12, 1e13), e in (1.5e5, 2.0e5)]
@everywhere function modify(problem, p)
	m = problem.material
	reaction = Reaction(:decomp, 1 => 2; A = p.A, E = p.E, h = 2e5)
	material = Material(; name = m.name, components = m.components,
	                    reactions = (reaction,), mixing = m.mixing)
	return PyrolysisProblem(; mesh = problem.mesh, material, bc_set = problem.bc_set,
	                        T_initial = problem.T_initial, ξ_initial = problem.ξ_initial,
	                        tspan = problem.tspan, output_times = problem.output_times)
end
@everywhere extract(sol) = (peak_T = maximum(sol.surface_T),)

results = parameter_sweep(problem, grid; modify, extract)
```
"""
function parameter_sweep(
    problem::PyrolysisProblem{<:Mesh},
    grid;
    modify::Function = _default_sweep_modify,
    extract::Function = identity,
    parallel::Bool = nworkers() > 1,
    solve_kwargs...,
)
    work = point -> extract(solve(modify(problem, point); solve_kwargs...))
    points = collect(grid)
    return parallel ? pmap(work, points) : map(work, points)
end

_default_sweep_modify(_, _) = throw(
    ArgumentError(
        "parameter_sweep needs `modify = (problem, point) -> problem`, which builds " *
        "the problem of each grid point",
    ),
)
