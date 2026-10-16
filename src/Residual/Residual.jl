"""
	Residual

The unified, pure ODE right-hand side (`residual!`) for all non-P1 pyrolysis
simulations.

`traits.jl` deliberately **stays in `Problem`**: the trait types
(`SimulationMode`, `Eulerian`, `Fickian`, …) are *upstream* of `Problem`'s
state/problem definitions (they parameterize `StateLayout`/`ProblemDef`), while
`residual!` is *downstream* of them. Co-locating both in one module would form a
dependency cycle, so this module is loaded **after** `Problem` and reaches the
trait structs via explicit imports below.

P1 radiation has its own residual in `Pyrolysis.Experimental.RadiationP1`.
"""
module Residual

using ..Materials
using ..Materials: n_components, NO_RADIATION, Material
using ..Geometry
using ..Properties
using ..Physics
using ..Discretization
using ForwardDiff

# Problem surface — exported types/accessors plus the unexported internals the
# residual kernel dispatches on (trait structs) and reuses (chi/geometry helpers).
using ..Problem
using ..Problem:
    ProblemDef,
    Workspace,
    StateLayout,
    get_option,
    chi_range,
    mesh_mode,
    chi_mode,
    populate!,
    update_cross_section_area!,
    update_state_geometry!,
    _chi_bar_from_state,
    Eulerian,
    ALE,
    NoChi,
    WithChi

include("residual.jl")

export residual!, column_chi_bar_from_state, cross_section_area_from_state

end # module Residual
