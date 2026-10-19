# ============================================================================
# Solution
# ============================================================================

"""
	SolverStatistics

Counters of a solve: accepted and rejected steps, right-hand-side
evaluations, Jacobian evaluations, linear solves, nonlinear-solver failures,
the steps the domain test rejected, and the largest scaled constraint residual
of the accepted states (one or below passes; zero for a fixed-step solve).
"""
Base.@kwdef struct SolverStatistics
    accepted::Int
    rejected::Int
    rhs_evaluations::Int
    jacobians::Int
    linear_solves::Int
    nonlinear_failures::Int
    domain_rejections::Int
    largest_constraint_residual::Float64
end

"""
	CoreExtras

The secondary series of a [`CoreSolution`](@ref), under the names of
`PyrolysisSolution.extras`, one value per saved time unless noted. Quantities
of the `:top` (surface) or `:bottom` (back) faces are `NaN` without that tag
or once no live cell borders it; per-area quantities divide by the current
area of the faces they name unless noted.

- `gas_generation_rate`: the gas the reactions of the live cells produce, per
  unit initial area of the `:top` faces [kg/(m²·s)]
- `back_temperature`: the temperature of the `:bottom` faces [K]
- `total_mass`: the summed component masses, per unit initial area of the
  `:top` faces [kg/m²]
- `thickness_history`: the distance from the first node to the outermost live
  one [m]
- `surface_mass_flux`: the outward mass flux of the `:top` faces [kg/(m²·s)]
- `surface_heat_flux_net`: the heat flux the `:top` condition applies at the
  face, positive into the domain (under a prescribed temperature, the
  conducted flux) [W/m²]; under `BEER_LAMBERT` the incident flux enters in
  depth and is not part of it (`volumetric_absorption`)
- `surface_heat_flux_absorbed`, `surface_heat_flux_reradiation`,
  `surface_heat_flux_convection`: its incident radiation absorbed at the face
  and its re-radiative and convective losses (positive leaving) [W/m²]
- `surface_heat_flux_conduction`: the flux conducted from the `:top` faces
  into their cells, `k (T_s − T_cell) / d` [W/m²]
- `surface_blowing_ratio`, `surface_pooling_factor`, `surface_face_fraction`,
  `surface_face_buoyancy`: the blowing diagnostics of the `:top` surface
  balance (1, 1, 0 and 1 without outflow, without a blowing model, or on a
  face its condition fixes in closed form)
- `volumetric_absorption`, `radiation_transmitted`: the intensity the cells
  absorbed and the intensity leaving through the last face of the in-depth
  beam, per unit exposed area [W/m²] (zero without in-depth radiation)
- `surface_emissivity`, `surface_conductivity`: the emissivity and the
  conductivity [W/(m·K)] of the cell under the `:top` faces
- `back_heat_flux`: the heat flux the `:bottom` condition applies, positive
  into the domain [W/m²]
- `MLR_gauge`: the mass-loss rate per unit initial area of the `:top` faces
  [kg/(m²·s)]
- `cross_section_area_history`: the cross-section measure (the area of a
  planar column) [m²]
- `live_cells`: the number of live cells
- `park_times`: the time of every deactivation, one entry per deactivated
  cell, in order [s]
"""
const CoreExtras = @NamedTuple{
    gas_generation_rate::Vector{Float64},
    back_temperature::Vector{Float64},
    total_mass::Vector{Float64},
    thickness_history::Vector{Float64},
    surface_mass_flux::Vector{Float64},
    surface_heat_flux_net::Vector{Float64},
    surface_heat_flux_absorbed::Vector{Float64},
    surface_heat_flux_reradiation::Vector{Float64},
    surface_heat_flux_convection::Vector{Float64},
    surface_heat_flux_conduction::Vector{Float64},
    surface_blowing_ratio::Vector{Float64},
    surface_pooling_factor::Vector{Float64},
    surface_face_fraction::Vector{Float64},
    surface_face_buoyancy::Vector{Float64},
    volumetric_absorption::Vector{Float64},
    radiation_transmitted::Vector{Float64},
    surface_emissivity::Vector{Float64},
    surface_conductivity::Vector{Float64},
    back_heat_flux::Vector{Float64},
    MLR_gauge::Vector{Float64},
    cross_section_area_history::Vector{Float64},
    live_cells::Vector{Int},
    park_times::Vector{Float64},
}

"""
	CoreSolution{P}

The solution of the finite-volume core at the saved times:

- `problem`: the problem that was solved
- `t`: saved times [s]
- `T`: cell temperatures, `n_cells × n_times` [K], `NaN` in a deactivated cell
- `ξ`: concentrations, `NC × n_cells × n_times` [kg/m³]
- `m`: component masses per cell, `NC × n_cells × n_times` [kg]
- `z`: node positions, `n_nodes × n_times` [m] on a moving mesh, `0 × 0` on a
  fixed one; `NaN` beyond the outermost live node
- `cell_centers`: cell-center positions, `n_cells × n_times` [m], `NaN` in a
  deactivated cell
- `boundary_T`: per boundary tag, the area-weighted surface temperature of
  its faces at every saved time [K]; after cells beyond a tag's faces are
  deactivated, of the face bounding the live cells on that side
- `surface_T`: the surface temperature of the `:top` faces (`NaN` without
  that tag) [K]
- `mass_loss_rate`: the mass leaving through every boundary face, per unit
  current area of the `:top` faces [kg/(m²·s)]
- `boundary_mass_flux`: per boundary tag, the outward mass flux of every
  component summed, area-weighted over its faces [kg/(m²·s)]
- `boundary_blowing_ratio`: per boundary tag, the area-weighted ratio
  `h/h₀` its convective condition applied at the solved face temperature (1
  without outflow, without a blowing model, or on a face its condition fixes
  in closed form)
- `extras`: the [`CoreExtras`](@ref) series
- `events`: the [`DeactivationEvent`](@ref)s of the solve, in order
- `outcome`: `:running` unless the deactivation ended the solve
  (`:burned_out`, `:interior`)
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
    mass_loss_rate::Vector{Float64}
    boundary_mass_flux::Dict{Symbol,Vector{Float64}}
    boundary_blowing_ratio::Dict{Symbol,Vector{Float64}}
    extras::CoreExtras
    events::Vector{DeactivationEvent}
    outcome::Symbol
    retcode::Symbol
    solve_time::Float64
    statistics::SolverStatistics
    raw_solution::Any
end

# The extras of `nt` saved times: every series zero except the blowing
# diagnostics, which start at their no-outflow values.
_extras(nt::Int) = CoreExtras((
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    ones(nt),
    ones(nt),
    zeros(nt),
    ones(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(nt),
    zeros(Int, nt),
    Float64[],
))

# The live faces of boundary tag `tag`.
_live_faces(topo, active, tag) =
    filter!(!=(0), [live_face(topo, active, f) for f in boundary_faces(topo, tag)])

# The outward mass flux of every component through face `f` [kg/(m²·s)].
_face_mass_fluxes(c::Contributions{NC}, f) where {NC} = ntuple(j -> c.F_mass[j, f], Val(NC))

# The area-weighted heat fluxes of the faces `faces` under their thermal
# condition: the applied (net), absorbed, re-radiated and convected fluxes of
# the condition at the solved face temperatures, and the conducted flux.
function _face_heat_fluxes(thermal, tag, faces, c, ds, topo, material, t)
    bc = boundary_condition(thermal, tag)
    net, absorbed, reradiation, convection, conduction = 0.0, 0.0, 0.0, 0.0, 0.0
    A = 0.0
    for f in faces
        i = topo.faces[f].cells[1]
        T_s = c.T_face[f]
        J = _face_mass_fluxes(c, f)
        q = compute_boundary_heat_fluxes(
            bc,
            T_s,
            Float64(t),
            ds.ε[i],
            thermal.radiation_model,
            BlowingFace(J, T_s, material);
            a_eff = ds.a[i],
        )
        q_cond = conducted_flux(thermal, c, ds, topo, f)
        a = ds.geometry.face_area[f]
        net += (bc isa DirichletThermalBC ? q_cond : q.net) * a
        absorbed += q.absorbed * a
        reradiation += q.reradiation * a
        convection += q.convection * a
        conduction += q_cond * a
        A += a
    end
    return (;
        net = net / A,
        absorbed = absorbed / A,
        reradiation = reradiation / A,
        convection = convection / A,
        conduction = conduction / A,
    )
end

"""
	CoreSolution(problem, t, states, def, ws, retcode, solve_time, statistics, raw;
	             deactivation = nothing)

Reconstruct the fields from the saved state vectors `states`: one residual
evaluation per saved state, under the mask and reference masses its time had
(the events of the [`CellDeactivation`](@ref) `deactivation` before it), fills
the derived state and the kernel contributions (the face temperatures and
mass fluxes among them) in the plain buffers of `ws`. The mask and the
references are left as after the last event.
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
    raw;
    deactivation::Union{Nothing,CellDeactivation} = nothing,
)
    layout = def.layout
    topo = def.setup.topology
    material = def.setup.material
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
    mass_loss_rate = zeros(nt)
    x = _extras(nt)
    gas = map(is_gas, material.components)
    thermal = find_kernel(ThermalBoundaryKernel, def.kernels)
    beam = find_kernel(BeerLambertKernel, def.kernels)
    ds, c = buffers(ws, def, Float64)
    du = zeros(state_length(layout))
    ref = def.setup.reference
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, ref.positions, ref.cross_section, topo)
    has_top = haskey(topo.boundary_tags, :top)
    has_bottom = haskey(topo.boundary_tags, :bottom)
    top = has_top ? boundary_faces(topo, :top) : Int[]
    exposed_area = has_top ? sum(g0.face_area[f] for f in top) : NaN
    for (k, u) in enumerate(states)
        deactivation === nothing || replay_events!(deactivation, layout, t[k])
        residual!(du, u, def, t[k], ws)
        @views T[:, k] .= ds.T
        @views ξ[:, :, k] .= ds.ξ
        @views m[:, :, k] .= ds.m
        for i = 1:n
            centers[i, k] = ds.geometry.cell_center[i][1]
            if !ds.active[i]
                T[i, k] = NaN
                centers[i, k] = NaN
            end
            ds.active[i] || continue
            for j = 1:NC
                gas[j] && (x.gas_generation_rate[k] += c.S_mass[j, i])
            end
        end
        x.gas_generation_rate[k] /= exposed_area
        x.live_cells[k] = count(ds.active)
        top_node = has_top ? live_face(topo, ds.active, only(top)) : n_nodes(layout)
        if is_moving(layout)
            @views z[:, k] .= ds.x
            @views z[(top_node+1):end, k] .= NaN
        end
        x.thickness_history[k] = ds.x[max(top_node, 1)] - ds.x[1]
        x.cross_section_area_history[k] = ds.A[]
        x.total_mass[k] = sum(ds.m) / exposed_area
        outflow = 0.0
        for tag in tags
            faces = _live_faces(topo, ds.active, tag)
            if isempty(faces)
                boundary_T[tag][k] = NaN
                boundary_mass_flux[tag][k] = NaN
                boundary_blowing_ratio[tag][k] = NaN
                continue
            end
            A = sum(ds.geometry.face_area[f] for f in faces)
            boundary_T[tag][k] =
                sum(c.T_face[f] * ds.geometry.face_area[f] for f in faces) / A
            out = sum(sum(view(c.F_mass, :, f)) * ds.geometry.face_area[f] for f in faces)
            outflow += out
            boundary_mass_flux[tag][k] = out / A
            boundary_blowing_ratio[tag][k] =
                sum(
                    thermal.states[f].blowing_ratio * ds.geometry.face_area[f] for
                    f in faces
                ) / A
        end
        top_faces = has_top ? _live_faces(topo, ds.active, :top) : Int[]
        if isempty(top_faces)
            mass_loss_rate[k] = NaN
            x.MLR_gauge[k] = NaN
            for s in (
                x.surface_heat_flux_net,
                x.surface_heat_flux_absorbed,
                x.surface_heat_flux_reradiation,
                x.surface_heat_flux_convection,
                x.surface_heat_flux_conduction,
                x.surface_blowing_ratio,
                x.surface_pooling_factor,
                x.surface_face_fraction,
                x.surface_face_buoyancy,
                x.surface_emissivity,
                x.surface_conductivity,
            )
                s[k] = NaN
            end
        else
            A_top = sum(ds.geometry.face_area[f] for f in top_faces)
            mass_loss_rate[k] = outflow / A_top
            x.MLR_gauge[k] = outflow / exposed_area
            q = _face_heat_fluxes(thermal, :top, top_faces, c, ds, topo, material, t[k])
            x.surface_heat_flux_net[k] = q.net
            x.surface_heat_flux_absorbed[k] = q.absorbed
            x.surface_heat_flux_reradiation[k] = q.reradiation
            x.surface_heat_flux_convection[k] = q.convection
            x.surface_heat_flux_conduction[k] = q.conduction
            f = first(top_faces)
            i = topo.faces[f].cells[1]
            state = thermal.states[f]
            x.surface_blowing_ratio[k] = boundary_blowing_ratio[:top][k]
            x.surface_pooling_factor[k] = state.pooling_factor
            x.surface_face_fraction[k] = state.face_fraction
            x.surface_face_buoyancy[k] = state.face_buoyancy
            x.surface_emissivity[k] = ds.ε[i]
            x.surface_conductivity[k] = ds.k[i]
        end
        bottom_faces = has_bottom ? _live_faces(topo, ds.active, :bottom) : Int[]
        if isempty(bottom_faces)
            x.back_temperature[k] = NaN
            x.back_heat_flux[k] = NaN
        else
            x.back_temperature[k] = boundary_T[:bottom][k]
            x.back_heat_flux[k] = _face_heat_fluxes(
                thermal,
                :bottom,
                bottom_faces,
                c,
                ds,
                topo,
                material,
                t[k],
            ).net
        end
        if beam !== nothing
            e, b = max(exposed_face(beam, ds), 1), transmitted_face(beam)
            A_e = ds.geometry.face_area[e]
            x.radiation_transmitted[k] = ds.I[b] * ds.geometry.face_area[b] / A_e
            x.volumetric_absorption[k] = ds.I[e] - x.radiation_transmitted[k]
        end
    end
    deactivation === nothing || replay_events!(deactivation, layout, Inf)
    deactivation === nothing || append!(x.park_times, e.t for e in deactivation.events)
    surface_T = haskey(boundary_T, :top) ? copy(boundary_T[:top]) : fill(NaN, nt)
    if haskey(boundary_mass_flux, :top)
        x.surface_mass_flux .= boundary_mass_flux[:top]
    else
        fill!(x.surface_mass_flux, NaN)
    end
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
        mass_loss_rate,
        boundary_mass_flux,
        boundary_blowing_ratio,
        x,
        deactivation === nothing ? DeactivationEvent[] : copy(deactivation.events),
        deactivation === nothing ? :running : deactivation.outcome[],
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
