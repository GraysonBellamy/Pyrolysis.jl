# ============================================================================
# Accessors
# ============================================================================

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
	get_profile(sol::PyrolysisSolution, field, i) -> Array

The field `field` at the saved time `sol.t[i]`: `:T` (the cell
temperatures), `:ξ` (the concentrations, `NC × n_cells`) or `:ξj` for
component `j`, `:m` (the component masses, `NC × n_cells`) or `:mj`, and `:z`
(the node positions; on a fixed mesh the mesh's). Deactivated cells are
`NaN`, and so are the nodes beyond the outermost live one.
"""
function get_profile(sol::PyrolysisSolution, field::Symbol, i::Integer)
    checkbounds(sol.t, i)
    field === :T && return sol.T[:, i]
    field === :ξ && return sol.ξ[:, :, i]
    if field === :m
        m = sol.m[:, :, i]
        m[:, isnan.(sol.T[:, i])] .= NaN
        return m
    end
    if field === :z
        isempty(sol.z) || return sol.z[:, i]
        z = reduce(vcat, sol.problem.mesh.positions)
        live = findlast(isfinite, view(sol.T, :, i))
        z[(isnothing(live) ? 1 : live+2):end] .= NaN
        return z
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

"""
	(sol::PyrolysisSolution)(t) -> NamedTuple

The fields at time `t`: `(; t, T, ξ, m, z, cell_centers)` as one saved time
of the solution holds them, reconstructed from the integrator's interpolant
between the saved states through the state mapping, with the mask and the
tables of the last deactivation event at or before `t`. The interpolant between two saved
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
    # Replay on a copy: a query must not change the solution's mask or tables,
    # including when a property evaluation throws.
    def, cells = deepcopy((sol.definition, sol.deactivation))
    layout = def.layout
    topo = def.setup.topology
    cells === nothing || replay_events!(cells, layout, t)
    warm = [isfinite(x) ? x : T_REF for x in sol.T[:, k]]
    ws = ResidualWorkspace(def, warm)
    ds, _ = buffers(ws, def, Float64)
    derive!(ds, u, layout, def.setup)
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
    return (; t = Float64(t), T, ξ, m, z, cell_centers = centers)
end

"""
	interior_events(sol::PyrolysisSolution) -> Vector{DeactivationEvent}

The events of `sol` whose cell left from under live cells, which were spent: a
void of the transfer's making that the pass closed ([`deactivate!`](@ref)). A void
beneath a skin is never among them: it ends the solve with outcome `:interior`.
"""
interior_events(sol::PyrolysisSolution) = filter(e -> e.interior, sol.events)

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
    println(io, "  Interior events: $(length(interior_events(sol)))")
    println(io, "  Accepted steps: $(sol.statistics.accepted)")
    print(io, "  Solve time: $(round(sol.solve_time, digits = 3)) s")
end

"""Gas fractions in `sol.gas.species` order; basis is `:mole` or `:mass`, zero at vacuum."""
function get_gas_composition(sol::PyrolysisSolution; basis = :mole)
    basis===:mole && return sol.gas.mole_fraction
    basis===:mass && return sol.gas.mass_fraction
    throw(ArgumentError("composition basis must be :mole or :mass"))
end
"""
    get_boundary_exchange(sol, tag)

Exact assembled outward species mass flux and energy flux, inward reservoir
heat flux, carried enthalpy, surface partial pressure, Darcy velocity, and
maximum scaled surface residual/iteration count at each saved time. Fluxes
are area-weighted over live faces [kg/(m² s), W/m²]; see `GasDiagnostics`.
"""
get_boundary_exchange(sol::PyrolysisSolution, tag::Symbol) =
    (; species = sol.gas.species, sol.gas.boundary[tag]...)
