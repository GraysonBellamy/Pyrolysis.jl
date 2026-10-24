"""
	Pyrolysis.Internal

Umbrella namespace re-exporting every Pyrolysis submodule's public symbol in
one place. The top-level `Pyrolysis` module exports a curated stable surface;
`Internal` is the escape hatch for advanced users who want the union of
submodule exports without reaching into individual submodules: the residual
(`residual!`), the colored Jacobian (`jacobian!`, `dense_jacobian`), the state
layout and its indices, the kernels and the closures of a set-up problem.

```julia
using Pyrolysis.Internal   # everything submodules export
```

Symbols loaded here are not part of the stable public API and may move or
change between minor releases. Prefer the top-level `Pyrolysis` exports for
production code; reach into `Pyrolysis.Internal` only when an internal building
block is genuinely needed. A name two submodules export as different objects
is not re-exported; reach it through its submodule.
"""
module Internal

# Resolve all submodule references through the parent (`Pyrolysis`) instead of
# bare names, which could shadow Base names.
const _SUBMODULES = let p = parentmodule(@__MODULE__)
    (
        p.Materials,
        p.Properties,
        p.BoundaryConditions,
        p.Physics,
        p.Meshes,
        p.States,
        p.Kernels,
        p.Closures,
        p.Assembly,
        p.Deactivation,
        p.ColoredJacobian,
        p.Problem,
        p.Outputs,
        p.Integration,
    )
end

let owners = Dict{Symbol,Vector{Any}}()
    for m in _SUBMODULES, n in names(m; all = false)
        n === nameof(m) && continue
        push!(get!(owners, n, Any[]), getfield(m, n))
    end
    for (n, objects) in owners
        all(o -> o === first(objects), objects) || continue
        Core.eval(@__MODULE__, Expr(:const, Expr(:(=), n, first(objects))))
        Core.eval(@__MODULE__, Expr(:export, n))
    end
end

end # module Internal
