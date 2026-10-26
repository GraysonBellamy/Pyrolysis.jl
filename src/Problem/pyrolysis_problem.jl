"""
The PyrolysisProblem type.
"""

# ============================================================================
# PyrolysisProblem Type
# ============================================================================

"""
	PyrolysisProblem{M,Mat,BCS,FT,Fξ} <: AbstractPyrolysisProblem

Complete specification of a pyrolysis simulation problem.

# Type Parameters
- `M`: Mesh type
- `Mat`: Material type
- `BCS`: BoundaryConditionSet type (type-stable BC storage)
- `FT`: Initial temperature function type
- `Fξ`: Initial concentration functions type

# Fields
- `mesh::M`: The computational mesh
- `material::Mat`: Material definition with components and reactions
- `bc_set::BCS`: BoundaryConditionSet with all boundary conditions (type-stable)
- `T_initial::FT`: Initial temperature function z -> T(z)
- `ξ_initial::Vector{Fξ}`: Initial concentration functions [z -> ξⱼ(z)]
- `tspan::Tuple{Float64,Float64}`: Simulation time span (t_start, t_end)
- `dt_initial::Float64`: Initial timestep [s]
- `output_times::Vector{Float64}`: Times at which to save output

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
	mass    = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(
	mesh = create_uniform_mesh(0.01, 100),
	material = pmma,
	bc_set = bc_set,
	T_initial = z -> 300.0,
	ξ_initial = [z -> 1190.0, z -> 0.0],
	tspan = (0.0, 100.0)
)
```
"""
struct PyrolysisProblem{M,Mat,BCS,FT,Fξ} <: AbstractPyrolysisProblem
    mesh::M
    material::Mat
    bc_set::BCS
    T_initial::FT
    ξ_initial::Vector{Fξ}
    tspan::Tuple{Float64,Float64}
    dt_initial::Float64
    output_times::Vector{Float64}
end

# ============================================================================
# Constructors
# ============================================================================

"""
	PyrolysisProblem(; mesh, material, bc_set::BoundaryConditionSet,
					 T_initial=z->300.0, ξ_initial=nothing,
					 tspan, dt_initial=1e-4, output_times=Float64[])

Create a PyrolysisProblem with BoundaryConditionSet (preferred, type-stable).

Radiative heating, convective coupling, and ambient temperatures are carried
exclusively by `bc_set` (`RadiativeFluxBC`, `ConvectiveBC`, `RadiativeBC`).

# Arguments
- `mesh`: The computational mesh
- `material`: Material definition
- `bc_set`: BoundaryConditionSet with all boundary conditions (type-stable)
- `T_initial=z->300.0`: Initial temperature function
- `ξ_initial=nothing`: Initial concentration functions (defaults from material)
- `tspan`: Simulation time span (required)
- `dt_initial=1e-4`: Initial timestep
- `output_times=Float64[]`: Times to save output (auto-generated if empty)

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3) + RadiativeBC(ε = 0.9, T_inf = 300),
	           bottom = AdiabaticBC()),
	mass    = (top = AmbientGasBC(GasReservoir(P=0., T=300.)), bottom = ImpermeableBC()),
)

problem = PyrolysisProblem(
	mesh = mesh,
	material = pmma,
	bc_set = bc_set,
	T_initial = z -> 300.0,
	tspan = (0.0, 100.0)
)
```
"""
function PyrolysisProblem(;
    mesh,
    material,
    bc_set::BoundaryConditionSet,
    T_initial = z -> 300.0,
    ξ_initial = nothing,
    tspan::Tuple{Real,Real},
    dt_initial::Real = 1e-4,
    output_times::Vector{<:Real} = Float64[],
)
    NC = n_components(material)

    # Create default initial concentration functions if not provided
    if ξ_initial === nothing
        # Default: first solid component at full density, others at zero
        ξ_funcs = Vector{Function}(undef, NC)
        for j = 1:NC
            comp = material.components[j]
            if j == 1 && is_solid(comp)
                # First solid at reference density
                ρ_ref = comp.density(300.0)
                ξ_funcs[j] = z -> ρ_ref
            else
                ξ_funcs[j] = z -> 0.0
            end
        end
        ξ_initial = ξ_funcs
    end

    # Convert ξ_initial to vector of functions if needed
    if ξ_initial isa Function
        ξ_initial = [ξ_initial for _ = 1:NC]
    end

    # Auto-generate output times if empty
    output_times_float = Float64.(output_times)
    if isempty(output_times_float)
        t_start, t_end = Float64.(tspan)
        n_outputs = min(101, max(11, round(Int, (t_end - t_start) / dt_initial / 100)))
        output_times_float = collect(range(t_start, t_end, length = n_outputs))
    end

    PyrolysisProblem(
        mesh,
        material,
        bc_set,
        T_initial,
        ξ_initial,
        (Float64(tspan[1]), Float64(tspan[2])),
        Float64(dt_initial),
        output_times_float,
    )
end

"""
	PyrolysisProblem(mesh, material, bc_set::BoundaryConditionSet, tspan; kwargs...)

Simplified constructor for common use cases with type-stable BoundaryConditionSet.

# Example
```julia
bc_set = BoundaryConditionSet(
	thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
)
problem = PyrolysisProblem(mesh, pmma, bc_set, (0.0, 100.0))
```
"""
function PyrolysisProblem(
    mesh,
    material,
    bc_set::BoundaryConditionSet,
    tspan::Tuple;
    kwargs...,
)
    PyrolysisProblem(;
        mesh = mesh,
        material = material,
        bc_set = bc_set,
        tspan = tspan,
        kwargs...,
    )
end

# ============================================================================
# Display
# ============================================================================

# The number of cells of a problem's mesh.
_n_cells(mesh::Meshes.Mesh) = Meshes.n_cells(mesh.topology)

function Base.show(io::IO, prob::PyrolysisProblem)
    t_start, t_end = prob.tspan
    n_cells = _n_cells(prob.mesh)
    NC = n_components(prob.material)
    print(
        io,
        "PyrolysisProblem($(prob.material.name), $n_cells cells, $NC components, t∈[$t_start, $t_end])",
    )
end

function Base.show(io::IO, ::MIME"text/plain", prob::PyrolysisProblem)
    t_start, t_end = prob.tspan
    n_cells = _n_cells(prob.mesh)
    NC = n_components(prob.material)
    NR = n_reactions(prob.material)

    println(io, "PyrolysisProblem")
    println(io, "  Material: $(prob.material.name)")
    println(io, "  Components: $NC")
    println(io, "  Reactions: $NR")
    println(io, "  Mesh: $n_cells cells")
    println(io, "  Time span: [$t_start, $t_end] s")
    println(io, "  Initial dt: $(prob.dt_initial) s")
    println(io, "  Output times: $(length(prob.output_times))")
    println(
        io,
        "  Thermal BCs: top=$(typeof(prob.bc_set.thermal.top)), bottom=$(typeof(prob.bc_set.thermal.bottom))",
    )
    println(
        io,
        "  Mass BCs: top=$(typeof(prob.bc_set.mass.top)), bottom=$(typeof(prob.bc_set.mass.bottom))",
    )
end
