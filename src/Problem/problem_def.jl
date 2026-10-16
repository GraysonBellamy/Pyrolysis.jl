"""
ProblemDef — immutable problem definition.

A `ProblemDef` carries the *unchanging* contract of a pyrolysis problem:
material, topology, boundary conditions, simulation mode and runtime options.
Every mutable cache lives in [`Workspace`](@ref) instead.

This split is what makes the residual a pure function of `(du, u, def, t, ws)`
and unlocks parameter sweeps, sensitivity analysis, and pluggable Jacobian
backends in later phases.

Field access is direct: read `def.material`, `def.topology`, `def.bcs`,
`def.mode`, `def.options`. The single helper is
[`get_option`](@ref) for fallback semantics on the runtime-options
NamedTuple.

# Type Parameters
- `M`:    Material type
- `Top`:  MeshTopology type
- `B`:    BoundaryConditionSet type
- `Mode`: bundled SimulationMode type
- `Opts`: NamedTuple of runtime options
"""
struct ProblemDef{M,Top,B,Mode<:SimulationMode,Opts<:NamedTuple}
    material::M
    topology::Top
    bcs::B
    mode::Mode
    options::Opts
end

"""
	ProblemDef(; material, topology, bcs, mode, options=NamedTuple())

Direct keyword constructor for advanced users. Most callers will build a
`PyrolysisProblem` first and use [`to_problem_def`](@ref) to derive the
`ProblemDef`.
"""
ProblemDef(;
    material,
    topology,
    bcs,
    mode::SimulationMode,
    options::NamedTuple = NamedTuple(),
) = ProblemDef(material, topology, bcs, mode, options)

"""
	get_option(def::ProblemDef, key::Symbol, default)

Look up a runtime option by name with a fallback default. NamedTuple lookup
is compile-time when `key` is statically known. This is a true convenience
because the alternative — `hasproperty(def.options, key) ? getproperty(def.options, key) : default` —
is verbose at the call site.
"""
@inline get_option(def::ProblemDef, key::Symbol, default) =
    hasproperty(def.options, key) ? getproperty(def.options, key) : default
