# ============================================================================
# Solution
# ============================================================================

"""
	SolverStatistics

Counters of a solve: accepted and rejected steps, right-hand-side
evaluations, Jacobian evaluations, linear solves and nonlinear-solver
failures.
"""
Base.@kwdef struct SolverStatistics
    accepted::Int
    rejected::Int
    rhs_evaluations::Int
    jacobians::Int
    linear_solves::Int
    nonlinear_failures::Int
end

"""
	CoreSolution{P}

The solution of the finite-volume core at the saved times:

- `problem`: the problem that was solved
- `t`: saved times [s]
- `T`: cell temperatures, `n_cells × n_times` [K]
- `ξ`: concentrations, `NC × n_cells × n_times` [kg/m³]
- `m`: component masses per cell, `NC × n_cells × n_times` [kg]
- `z`: node positions, `n_nodes × n_times` [m] on a moving mesh, `0 × 0` on a
  fixed one
- `cell_centers`: cell-center positions, `n_cells × n_times` [m]
- `boundary_T`: per boundary tag, the area-weighted surface temperature of
  its faces at every saved time [K]
- `surface_T`: the surface temperature of the `:top` faces (`NaN` without
  that tag) [K]
- `boundary_mass_flux`: per boundary tag, the outward mass flux of every
  component summed, area-weighted over its faces [kg/(m²·s)]
- `boundary_blowing_ratio`: per boundary tag, the area-weighted ratio
  `h/h₀` its convective condition applied at the solved face temperature (1
  without outflow, without a blowing model, or on a face its condition fixes
  in closed form)
- `extras`: the series under the names of `PyrolysisSolution.extras`:
  `surface_mass_flux`, the outward mass flux of the `:top` faces
  [kg/(m²·s)], and `surface_blowing_ratio`, their blowing ratio (`NaN`
  without that tag); `radiation_transmitted`, the intensity leaving through
  the last face of the in-depth beam, and `volumetric_absorption`, the
  intensity the cells absorbed, both per unit exposed area [W/m²] (zero
  without in-depth radiation); `thickness_history`, the distance between
  the end nodes [m]; `cross_section_area_history`, the cross-section measure
  (the area of a planar column) [m²]; `total_mass`, the summed component
  masses per unit initial area of the `:top` faces [kg/m²]
- `retcode`: the integrator's return code
- `solve_time`: wall time of the solve [s]
- `statistics`: the [`SolverStatistics`](@ref)
- `raw_solution`: the integrator's solution
"""
struct CoreSolution{P}
    problem::P
    t::Vector{Float64}
    T::Matrix{Float64}
    ξ::Array{Float64,3}
    m::Array{Float64,3}
    z::Matrix{Float64}
    cell_centers::Matrix{Float64}
    boundary_T::Dict{Symbol,Vector{Float64}}
    surface_T::Vector{Float64}
    boundary_mass_flux::Dict{Symbol,Vector{Float64}}
    boundary_blowing_ratio::Dict{Symbol,Vector{Float64}}
    extras::@NamedTuple{
        surface_mass_flux::Vector{Float64},
        surface_blowing_ratio::Vector{Float64},
        radiation_transmitted::Vector{Float64},
        volumetric_absorption::Vector{Float64},
        thickness_history::Vector{Float64},
        cross_section_area_history::Vector{Float64},
        total_mass::Vector{Float64},
    }
    retcode::Symbol
    solve_time::Float64
    statistics::SolverStatistics
    raw_solution::Any
end

"""
	CoreSolution(problem, t, states, def, ws, retcode, solve_time, statistics, raw)

Reconstruct the fields from the saved state vectors `states`: one residual
evaluation per saved state fills the derived state and the kernel
contributions (the face temperatures and mass fluxes among them) in the plain
buffers of `ws`.
"""
function CoreSolution(
    problem,
    t::AbstractVector{Float64},
    states::AbstractVector,
    def::CoreDefinition,
    ws::CoreWorkspace,
    retcode::Symbol,
    solve_time::Float64,
    statistics::SolverStatistics,
    raw,
)
    layout = def.layout
    topo = def.setup.topology
    NC = n_components(layout)
    n = n_cells(layout)
    nt = length(t)
    T = zeros(n, nt)
    ξ = zeros(NC, n, nt)
    m = zeros(NC, n, nt)
    z = is_moving(layout) ? zeros(n_nodes(layout), nt) : zeros(0, 0)
    centers = zeros(n, nt)
    tags = boundary_tags(topo)
    boundary_T = Dict(tag => zeros(nt) for tag in tags)
    boundary_mass_flux = Dict(tag => zeros(nt) for tag in tags)
    boundary_blowing_ratio = Dict(tag => zeros(nt) for tag in tags)
    transmitted = zeros(nt)
    absorbed = zeros(nt)
    thickness = zeros(nt)
    area = zeros(nt)
    total_mass = zeros(nt)
    thermal = find_kernel(ThermalBoundaryKernel, def.kernels)
    beam = find_kernel(BeerLambertKernel, def.kernels)
    ds, c = buffers(ws, def, Float64)
    du = zeros(state_length(layout))
    ref = def.setup.reference
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, ref.positions, ref.cross_section, topo)
    top = haskey(topo.boundary_tags, :top) ? boundary_faces(topo, :top) : Int[]
    exposed_area = isempty(top) ? NaN : sum(g0.face_area[f] for f in top)
    for (k, u) in enumerate(states)
        residual!(du, u, def, t[k], ws)
        @views T[:, k] .= ds.T
        @views ξ[:, :, k] .= ds.ξ
        @views m[:, :, k] .= ds.m
        is_moving(layout) && (@views z[:, k] .= ds.x)
        for i = 1:n
            centers[i, k] = ds.geometry.cell_center[i][1]
        end
        thickness[k] = ds.x[end] - ds.x[1]
        area[k] = ds.A[]
        total_mass[k] = sum(ds.m) / exposed_area
        for tag in tags
            faces = boundary_faces(topo, tag)
            A = sum(ds.geometry.face_area[f] for f in faces)
            boundary_T[tag][k] =
                sum(c.T_face[f] * ds.geometry.face_area[f] for f in faces) / A
            boundary_mass_flux[tag][k] =
                sum(sum(view(c.F_mass, :, f)) * ds.geometry.face_area[f] for f in faces) / A
            boundary_blowing_ratio[tag][k] =
                sum(
                    thermal.states[f].blowing_ratio * ds.geometry.face_area[f] for
                    f in faces
                ) / A
        end
        if beam !== nothing
            e, b = beam.exposed, transmitted_face(beam)
            A_e = ds.geometry.face_area[e]
            transmitted[k] = ds.I[b] * ds.geometry.face_area[b] / A_e
            absorbed[k] = ds.I[e] - transmitted[k]
        end
    end
    surface_T = haskey(boundary_T, :top) ? copy(boundary_T[:top]) : fill(NaN, nt)
    surface_mass_flux =
        haskey(boundary_mass_flux, :top) ? copy(boundary_mass_flux[:top]) : fill(NaN, nt)
    surface_blowing_ratio =
        haskey(boundary_blowing_ratio, :top) ? copy(boundary_blowing_ratio[:top]) :
        fill(NaN, nt)
    return CoreSolution(
        problem,
        Vector{Float64}(t),
        T,
        ξ,
        m,
        z,
        centers,
        boundary_T,
        surface_T,
        boundary_mass_flux,
        boundary_blowing_ratio,
        (;
            surface_mass_flux,
            surface_blowing_ratio,
            radiation_transmitted = transmitted,
            volumetric_absorption = absorbed,
            thickness_history = thickness,
            cross_section_area_history = area,
            total_mass,
        ),
        retcode,
        solve_time,
        statistics,
        raw,
    )
end

"""
	cell_positions(sol::CoreSolution) -> Vector{Float64}

The cell-center positions at the last saved time [m].
"""
cell_positions(sol::CoreSolution) = sol.cell_centers[:, end]
