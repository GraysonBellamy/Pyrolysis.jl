"""
	Pyrolysis.Internal

Umbrella namespace re-exporting every Pyrolysis submodule's public symbol in
one place. The top-level `Pyrolysis` module exports a curated stable surface;
`Internal` is the escape hatch for advanced users who want the union of
submodule exports without reaching into individual submodules.

```julia
using Pyrolysis.Internal   # everything submodules export
```

Symbols loaded here are not part of the stable public API and may move or
change between minor releases. Prefer the top-level `Pyrolysis` exports for
production code; reach into `Pyrolysis.Internal` only when an internal building
block is genuinely needed.
"""
module Internal

using ..Materials
using ..BoundaryConditions
using ..Geometry
using ..Properties
using ..Physics
using ..Adaptivity
using ..Discretization
using ..Problem
using ..Residual
using ..Jacobian
using ..Diagnostics
using ..Solver

# Resolve all submodule references through the parent (`Pyrolysis`) instead of
# bare names (some, e.g. a former `Core`, would shadow Base under Julia 1.12).
const _SUBMODULES = let p = parentmodule(@__MODULE__)
    (
        p.Materials,
        p.BoundaryConditions,
        p.Geometry,
        p.Properties,
        p.Physics,
        p.Adaptivity,
        p.Discretization,
        p.Problem,
        p.Residual,
        p.Jacobian,
        p.Diagnostics,
        p.Solver,
    )
end

for m in _SUBMODULES
    for n in names(m; all = false)
        n === nameof(m) && continue
        @eval export $n
    end
end

end # module Internal
