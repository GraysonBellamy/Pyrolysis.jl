# ============================================================================
# Accessors
# ============================================================================

"""
	get_times(sol::PyrolysisSolution) -> Vector{Float64}

The saved times [s].
"""
get_times(sol::PyrolysisSolution) = sol.t

"""
	get_temperatures(sol::PyrolysisSolution) -> Matrix{Float64}

The cell temperatures at the saved times, `n_cells × n_times` [K], `NaN` in a
deactivated cell.
"""
get_temperatures(sol::PyrolysisSolution) = sol.T

"""
	get_concentrations(sol::PyrolysisSolution) -> Array{Float64,3}

The concentrations at the saved times, `NC × n_cells × n_times` [kg/m³],
`NaN` in a deactivated cell.
"""
get_concentrations(sol::PyrolysisSolution) = sol.ξ

"""
	get_gas_generation_rate(sol::PyrolysisSolution) -> Vector{Float64}

The gas the reactions of the live cells produce, per unit initial area of the
`:top` faces, at the saved times [kg/(m²·s)]: a volumetric production, not
the mass-loss rate through the surface (`sol.mass_loss_rate`).
"""
get_gas_generation_rate(sol::PyrolysisSolution) = sol.extras.gas_generation_rate

"""
	n_timesteps(sol::PyrolysisSolution) -> Int

The number of saved times.
"""
n_timesteps(sol::PyrolysisSolution) = length(sol.t)

"""
	succeeded(sol::PyrolysisSolution) -> Bool

Whether the solve finished: with return code `:Success`, or `:Terminated`
by a complete burnout (outcome `:burned_out`) or the `min_thickness` stop
(outcome `:min_thickness`). An interior void (`:interior`) and `:Failure`
(among them an event that leaves the state inadmissible or the node
velocities inconsistent) are not; nor is a termination by a callback of the
caller's. [`completed`](@ref) is the strict test.
"""
succeeded(sol::PyrolysisSolution) =
    sol.retcode === :Success ||
    (sol.retcode === :Terminated && sol.outcome in (:burned_out, :min_thickness))

"""
	completed(sol::PyrolysisSolution) -> Bool

Whether the solve reached the end of its time span: return code `:Success`.
"""
completed(sol::PyrolysisSolution) = sol.retcode === :Success

"""
	get_cell_positions(sol::PyrolysisSolution; at = lastindex(sol.t)) -> Vector{Float64}

The cell-centre positions at the saved time `sol.t[at]` [m]: `NaN` for a
deactivated cell, and so for every cell beyond the outermost live node. The
initial positions are `get_cell_positions(sol; at = 1)`.
"""
get_cell_positions(sol::PyrolysisSolution; at::Integer = lastindex(sol.t)) =
    sol.cell_centers[:, at]

"""
	get_profile(sol::PyrolysisSolution, field, i) -> Array

The field `field` at the saved time `sol.t[i]`: `:T` (the cell
temperatures), `:ξ` (the concentrations, `NC × n_cells`) or `:ξj` for
component `j`, `:m` (the component masses, `NC × n_cells`) or `:mj`, and `:z`
(the node positions; on a fixed mesh the mesh's). Deactivated cells are
`NaN`, and so are the nodes beyond the outermost live one.
"""
function get_profile(sol::PyrolysisSolution, field::Symbol, i::Integer)
    field === :T && return sol.T[:, i]
    field === :ξ && return sol.ξ[:, :, i]
    if field === :m
        m = sol.m[:, :, i]
        m[:, isnan.(sol.T[:, i])] .= NaN
        return m
    end
    if field === :z
        isempty(sol.z) || return sol.z[:, i]
        return reduce(vcat, sol.problem.mesh.positions)
    end
    s = string(field)
    for (prefix, f) in (("ξ", :ξ), ("xi", :ξ), ("m", :m))
        startswith(s, prefix) || continue
        j = tryparse(Int, chopprefix(s, prefix))
        j === nothing && continue
        return get_profile(sol, f, i)[j, :]
    end
    throw(ArgumentError("unknown field :$field; use :T, :ξ, :ξj, :m, :mj or :z"))
end

# The mask of a state at time `t`: every event at or before `t` applied.
function _mask_at!(layout, events, t)
    fill!(layout.mask, true)
    for e in events
        e.t <= t || break
        layout.mask[e.rows] .= false
    end
    return layout
end

"""
	(sol::PyrolysisSolution)(t) -> NamedTuple
	interpolate_solution(sol::PyrolysisSolution, t) -> NamedTuple

The fields at time `t`: `(; t, T, ξ, m, z, cell_centers)` as one saved time
of the solution holds them, reconstructed from the integrator's interpolant
between the saved states through the state mapping, with the mask of the
last deactivation event at or before `t`. The interpolant between two saved
states that an event separates is not the trajectory, so a `t` strictly
inside such an interval is refused.
"""
function (sol::PyrolysisSolution)(t::Real)
    raw = sol.raw_solution
    ts = raw.t
    first(ts) <= t <= last(ts) || throw(
        ArgumentError(
            "t = $t lies outside the solution's span [$(first(ts)), $(last(ts))]",
        ),
    )
    k = searchsortedlast(ts, t)
    if ts[k] != t
        a, b = ts[k], ts[k+1]
        for e in sol.events
            a < e.t <= b && throw(
                ArgumentError(
                    "t = $t lies between the saved states at $a and $b s, which the " *
                    "deactivation of cell $(e.cell) at t = $(e.t) s (invocation " *
                    "$(e.invocation)) separates: the interpolant there is not the trajectory",
                ),
            )
        end
    end
    u = ts[k] == t ? raw.u[k] : raw(t)
    def = sol.definition
    layout = def.layout
    topo = def.setup.topology
    warm = [isfinite(x) ? x : T_REF for x in sol.T[:, k]]
    ws = CoreWorkspace(def, warm)
    _mask_at!(layout, sol.events, t)
    du = zeros(state_length(layout))
    residual!(du, u, def, t, ws)
    ds, _ = buffers(ws, def, Float64)
    n = n_cells(layout)
    T, ξ, m = copy(ds.T), copy(ds.ξ), copy(ds.m)
    centers = [ds.geometry.cell_center[i][1] for i = 1:n]
    for i = 1:n
        ds.active[i] && continue
        T[i] = NaN
        centers[i] = NaN
        ξ[:, i] .= NaN
    end
    z = Float64[]
    if is_moving(layout)
        z = copy(ds.x)
        if haskey(topo.boundary_tags, :top)
            top = live_face(topo, ds.active, only(boundary_faces(topo, :top)))
            z[(top+1):end] .= NaN
        end
    end
    _mask_at!(layout, sol.events, Inf)
    return (; t = Float64(t), T, ξ, m, z, cell_centers = centers)
end

"""
	interpolate_solution(sol::PyrolysisSolution, t) -> NamedTuple

The reconstructed fields at `t`, `sol(t)`: see [`PyrolysisSolution`](@ref)'s
call method for the interpolation and its refusal across an event's step.
"""
interpolate_solution(sol::PyrolysisSolution, t::Real) = sol(t)

function Base.show(io::IO, sol::PyrolysisSolution)
    print(
        io,
        "PyrolysisSolution(:$(sol.retcode), outcome :$(sol.outcome), $(length(sol.t)) saved times, ",
        "$(round(sol.solve_time, digits = 3)) s)",
    )
end

function Base.show(io::IO, ::MIME"text/plain", sol::PyrolysisSolution)
    println(io, "PyrolysisSolution")
    println(io, "  Return code: $(sol.retcode)")
    println(io, "  Outcome: $(sol.outcome)")
    println(io, "  Saved times: $(length(sol.t))")
    isempty(sol.t) || println(io, "  Time range: [$(sol.t[1]), $(sol.t[end])] s")
    isempty(sol.t) ||
        println(io, "  Live cells: $(sol.extras.live_cells[end]) of $(size(sol.T, 1))")
    println(io, "  Invocations: $(sol.statistics.invocations)")
    println(io, "  Accepted steps: $(sol.statistics.accepted)")
    print(io, "  Solve time: $(round(sol.solve_time, digits = 3)) s")
end
