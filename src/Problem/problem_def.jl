"""
ProblemDef — immutable problem definition.

A `ProblemDef` carries the *unchanging* contract of a pyrolysis problem:
material, topology, boundary conditions, initial conditions, time span,
options. Every mutable cache lives in [`Workspace`](@ref) instead.

This split is what makes the residual a pure function of `(du, u, def, t, ws)`
and unlocks parameter sweeps, sensitivity analysis, and pluggable Jacobian
backends in later phases.

Field access is direct: read `def.material`, `def.topology`, `def.bcs`,
`def.tspan`, `def.mode`, `def.options`. The single helper is
[`get_option`](@ref) for fallback semantics on the runtime-options
NamedTuple.

# Type Parameters
- `M`:    Material type
- `Top`:  MeshTopology type
- `B`:    BoundaryConditionSet type
- `IC`:   Initial-conditions container type
- `Mode`: bundled SimulationMode type
- `Opts`: NamedTuple of runtime options
"""
struct ProblemDef{M,Top,B,IC,Mode<:SimulationMode,Opts<:NamedTuple}
    material::M
    topology::Top
    bcs::B
    initial::IC
    tspan::Tuple{Float64,Float64}
    mode::Mode
    options::Opts
end

"""
	ProblemDef(; material, topology, bcs, initial,
	             tspan, mode, options=NamedTuple())

Direct keyword constructor for advanced users. Most callers will build a
`PyrolysisProblem` first and use [`to_problem_def`](@ref) to derive the
`ProblemDef`.
"""
function ProblemDef(;
    material,
    topology,
    bcs,
    initial,
    tspan::Tuple{<:Real,<:Real},
    mode::SimulationMode,
    options::NamedTuple = NamedTuple(),
)
    tspan_f = (Float64(tspan[1]), Float64(tspan[2]))
    ProblemDef(material, topology, bcs, initial, tspan_f, mode, options)
end

"""
	get_option(def::ProblemDef, key::Symbol, default)

Look up a runtime option by name with a fallback default. NamedTuple lookup
is compile-time when `key` is statically known. This is a true convenience
because the alternative — `hasproperty(def.options, key) ? getproperty(def.options, key) : default` —
is verbose at the call site.
"""
@inline get_option(def::ProblemDef, key::Symbol, default) =
    hasproperty(def.options, key) ? getproperty(def.options, key) : default

# ---------------------------------------------------------------------------
# Initial conditions container
# ---------------------------------------------------------------------------

"""
	InitialConditions{FT, Fξ}

Immutable initial-conditions descriptor used by `ProblemDef`. Mirrors the
`PyrolysisProblem` fields. Concrete types are stored so the descriptor is
type-stable for downstream consumers.
"""
struct InitialConditions{FT,Fξ}
    T_initial::FT
    ξ_initial::Vector{Fξ}
end
