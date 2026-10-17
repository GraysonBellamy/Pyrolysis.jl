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
    retcode::Symbol
    solve_time::Float64
    statistics::SolverStatistics
    raw_solution::Any
end

"""
	CoreSolution(problem, t, states, layout, ds, setup, retcode, solve_time, statistics, raw)

Reconstruct the fields from the saved state vectors `states` through the
state mapping into the plain derived state `ds`.
"""
function CoreSolution(
    problem,
    t::AbstractVector{Float64},
    states::AbstractVector,
    layout::StateLayout{NC},
    ds::DerivedState,
    setup::StateSetup,
    retcode::Symbol,
    solve_time::Float64,
    statistics::SolverStatistics,
    raw,
) where {NC}
    n = n_cells(layout)
    nt = length(t)
    T = zeros(n, nt)
    ξ = zeros(NC, n, nt)
    m = zeros(NC, n, nt)
    z = is_moving(layout) ? zeros(n_nodes(layout), nt) : zeros(0, 0)
    for (k, u) in enumerate(states)
        derive!(ds, u, layout, setup)
        @views T[:, k] .= ds.T
        @views ξ[:, :, k] .= ds.ξ
        @views m[:, :, k] .= ds.m
        is_moving(layout) && (@views z[:, k] .= ds.x)
    end
    return CoreSolution(
        problem,
        Vector{Float64}(t),
        T,
        ξ,
        m,
        z,
        retcode,
        solve_time,
        statistics,
        raw,
    )
end
