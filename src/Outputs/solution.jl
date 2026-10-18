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
    retcode::Symbol
    solve_time::Float64
    statistics::SolverStatistics
    raw_solution::Any
end

"""
	CoreSolution(problem, t, states, def, ws, retcode, solve_time, statistics, raw)

Reconstruct the fields from the saved state vectors `states`: one residual
evaluation per saved state fills the derived state and the kernel
contributions (the face temperatures among them) in the plain buffers of `ws`.
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
    ds, c = buffers(ws, def, Float64)
    du = zeros(state_length(layout))
    for (k, u) in enumerate(states)
        residual!(du, u, def, t[k], ws)
        @views T[:, k] .= ds.T
        @views ξ[:, :, k] .= ds.ξ
        @views m[:, :, k] .= ds.m
        is_moving(layout) && (@views z[:, k] .= ds.x)
        for i = 1:n
            centers[i, k] = ds.geometry.cell_center[i][1]
        end
        for tag in tags
            faces = boundary_faces(topo, tag)
            A = sum(ds.geometry.face_area[f] for f in faces)
            boundary_T[tag][k] =
                sum(c.T_face[f] * ds.geometry.face_area[f] for f in faces) / A
        end
    end
    surface_T = haskey(boundary_T, :top) ? copy(boundary_T[:top]) : fill(NaN, nt)
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
