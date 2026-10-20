# ============================================================================
# Deactivation policies and events
# ============================================================================

"""
	DeactivationPolicy

How a deactivated cell's inventories and volume leave the live mesh.
"""
abstract type DeactivationPolicy end

"""
	ConservativeTransfer <: DeactivationPolicy

The dead cells' inventories (masses, sensible enthalpy, history pair) go to
the live cell inward of them, the receiver, whose temperature follows from
the summed sensible enthalpy of all of them; the energy state is set from it.
The caloric extension's regularization, which only a cell far below the
heat-capacity floor carries, is not transferred: what the cells held of it
beyond the receiver's after the merge is written to the event. On a moving
mesh the receiver also absorbs the dead cells' volume: its far face moves out
to enclose the merged volume, where every node of the dead cells then
coincides (the moving-boundary collapse), and the merge never destroys pore
space: the merged volume is the larger of the summed volumes and the merged
skeletal volume plus the summed pre-merge pore volumes, the difference
written to the event as created volume. On a fixed mesh the receiver keeps
its volume and takes the condensed inventories only: its pores could hold
the dead cell's gas only at a higher pressure, so that gas is vented, its
masses and sensible enthalpy written to the event as lost.
"""
struct ConservativeTransfer <: DeactivationPolicy end

"""
	DeactivationEvent

One deactivated cell: its time, the callback invocation it belongs to
(numbered from one over the solve, shared by every cell one invocation
deactivates), the dead cell, the receiver (zero when none was left, the dead
cell's inventories then lost), the dead cell's condensed mass as a fraction
of its reference when it went, the rows masked, the volume the merge created
[m³], the masses [kg] and sensible enthalpy [J] lost, and the caloric
regularization [J] the cells held that the receiver does not (the extended
energies of the merged cells before the event less the receiver's after it
and the enthalpy lost). A merge of several cells into one receiver records
one event per dead cell, the outermost first; the merge's created volume,
losses and regularization are written to that first event and the others
carry zeros, so the sums over the events are the merges'.
"""
struct DeactivationEvent
    t::Float64
    invocation::Int
    cell::Int
    receiver::Int
    fraction::Float64
    rows::Vector{Int}
    created_volume::Float64
    lost_mass::Vector{Float64}
    lost_enthalpy::Float64
    regularization::Float64
end

"""
	CellDeactivation{P,NC,E}

The deactivation state of a solve: the policy, the threshold `ε_m`, the
near-threshold band `ε_band` (a fraction of a cell's reference, read on a
moving mesh only, see [`deactivate!`](@ref)), the cells' initial condensed
masses (`reference`; a receiver keeps its own, so a merged cell deactivates
once its whole merged inventory is below `ε_m` of what it started with; the
reaction kernel's switch reads the same masses), the condensed component
flags, the caloric extension merged temperatures are found with, the events
so far, and the outcome: `:running`; `:burned_out` once the last live cell
went; `:interior` when a cell below the threshold lies inward of a live one
that [`deactivate!`](@ref) leaves, which the policies here do not handle.
"""
struct CellDeactivation{P<:DeactivationPolicy,NC,E<:CaloricExtension{NC}}
    policy::P
    threshold::Float64
    band::Float64
    reference::Vector{Float64}
    solid::NTuple{NC,Bool}
    caloric::E
    events::Vector{DeactivationEvent}
    outcome::Base.RefValue{Symbol}
end

"""
	CellDeactivation(policy, threshold, band, reference, solid, caloric)

The deactivation state before a solve, `reference` the initial condensed
masses; `threshold ≤ band < 1`.
"""
function CellDeactivation(
    policy::DeactivationPolicy,
    threshold::Real,
    band::Real,
    reference::Vector{Float64},
    solid::NTuple{NC,Bool},
    caloric::CaloricExtension{NC},
) where {NC}
    threshold <= band < 1 || throw(
        ArgumentError("the band must lie in [threshold, 1) = [$threshold, 1), got $band"),
    )
    return CellDeactivation(
        policy,
        Float64(threshold),
        Float64(band),
        reference,
        solid,
        caloric,
        DeactivationEvent[],
        Ref(:running),
    )
end

"""
	reset_cells!(d, layout) -> d

Every row live and no events: the state a solve starts from.
"""
function reset_cells!(d::CellDeactivation, layout::StateLayout)
    fill!(layout.mask, true)
    empty!(d.events)
    d.outcome[] = :running
    return d
end

"""
	condensed_mass(d, layout, u, i)

The condensed (solid and liquid) mass of cell `i` in the state `u` [kg].
"""
@inline function condensed_mass(
    d::CellDeactivation{P,NC},
    layout::StateLayout{NC},
    u::AbstractVector,
    i::Int,
) where {P,NC}
    m = zero(eltype(u))
    @inbounds for j = 1:NC
        d.solid[j] && (m += u[mass_index(layout, i, j)])
    end
    return m
end

"""
	below_threshold(d, layout, u, i) -> Bool

Whether live cell `i` holds less than `ε_m` of its initial condensed mass. A
cell that started with none never deactivates.
"""
@inline below_threshold(d::CellDeactivation, layout::StateLayout, u, i::Int) =
    layout.mask[energy_index(layout, i)] &&
    d.reference[i] > 0 &&
    condensed_mass(d, layout, u, i) < d.threshold * d.reference[i]

"""
	any_below(d, layout, u) -> Bool

Whether any live cell is below the threshold: the event's condition.
"""
function any_below(d::CellDeactivation, layout::StateLayout, u::AbstractVector)
    for i = 1:n_cells(layout)
        below_threshold(d, layout, u, i) && return true
    end
    return false
end

"""
	top_live_cell(layout) -> Int

The outermost live cell, zero when none is live. Cells die from the exposed
(`:top`) side only, so the live cells are `1:top_live_cell(layout)`.
"""
function top_live_cell(layout::StateLayout)
    for i = n_cells(layout):-1:1
        layout.mask[energy_index(layout, i)] && return i
    end
    return 0
end

# Whether live cell `i` holds less than `ε_band` of its initial condensed mass.
@inline _in_band(d::CellDeactivation, layout::StateLayout, u, i::Int) =
    layout.mask[energy_index(layout, i)] &&
    d.reference[i] > 0 &&
    condensed_mass(d, layout, u, i) < d.band * d.reference[i]

# The inward end of the run of cells from `top` inward that are in the band.
function _run_start(d::CellDeactivation, layout::StateLayout, u, top::Int)
    lo = top
    while lo >= 1 && _in_band(d, layout, u, lo)
        lo -= 1
    end
    return lo + 1
end

"""
	deactivate!(u, d, def, t) -> Symbol

Deactivate, in one pass, the cells below the threshold at the exposed side.
On a moving mesh the run of live cells from the outermost one inward whose
condensed mass is below `ε_band` of their reference goes, when it holds a
cell below `ε_m`, in one merge into the live cell beneath it
([`ConservativeTransfer`](@ref)); a run that reaches the innermost cell
leaves that cell, holding the run's inventories, as the survivor. On a fixed
mesh the outermost live cell goes while it is below `ε_m`, each merge
re-testing the receiver with the dead cell's condensed inventories added.
The last live cell, once below `ε_m`, goes with no receiver, its inventories
written to the event as lost. The events of one call share their invocation.

Returns `:transferred`, `:burned_out` when the last live cell went, or
`:interior` when a cell below the threshold remains inward of a live one.
"""
function deactivate!(u::AbstractVector, d::CellDeactivation, def::CoreDefinition, t::Real)
    layout = def.layout
    invocation = isempty(d.events) ? 1 : last(d.events).invocation + 1
    while true
        top = top_live_cell(layout)
        top == 0 && break
        lo = is_moving(layout) ? _run_start(d, layout, u, top) : top
        any(i -> below_threshold(d, layout, u, i), lo:top) || break
        if top == 1
            _lose!(u, d, def, top, Float64(t), invocation)
            d.outcome[] = :burned_out
            return :burned_out
        end
        r = max(lo - 1, 1)
        _merge!(u, d, d.policy, def, r, (r+1):top, Float64(t), invocation)
    end
    if any_below(d, layout, u)
        d.outcome[] = :interior
        return :interior
    end
    return :transferred
end

"""
	replay_events!(d, layout, t) -> d

The mask of a state saved at `t`: every event of `d` at or before `t`
applied, starting from every row live. The integrator applies an event before
it saves the state of the step that triggered it, so a state saved at an
event's own time is the one after the event.
"""
function replay_events!(d::CellDeactivation, layout::StateLayout, t::Real)
    fill!(layout.mask, true)
    for e in d.events
        e.t <= t || break
        layout.mask[e.rows] .= false
    end
    return d
end

# ----------------------------------------------------------------------------
# Transfer
# ----------------------------------------------------------------------------

# The masses of cell `i`.
@inline _masses(u, layout::StateLayout{NC}, i) where {NC} =
    ntuple(j -> u[mass_index(layout, i, j)], Val(NC))

# The temperature of cell `i` from its energy state.
@inline _temperature(::TemperatureState, u, d, laws, layout, i, m) =
    u[energy_index(layout, i)]
@inline _temperature(::EnthalpyState, u, d, laws, layout, i, m) =
    first(invert_caloric(d.caloric, laws, u[energy_index(layout, i)], m, i, T_REF))

# Sensible enthalpy of the masses `m` at `T` [J].
_sensible_enthalpy(laws, T, m) = sum(map(*, m, component_enthalpies(laws, T)))

# The extended energy of cell `i` (the inventory itself in the enthalpy form).
@inline _energy(::EnthalpyState, u, d, laws, layout, i, m, T) = u[energy_index(layout, i)]
@inline _energy(::TemperatureState, u, d, laws, layout, i, m, T) =
    extended_energy(d.caloric, laws, T, m, i)

# The receiver's and the dead cell's nodes: the receiver's inward node, the
# node they share and the dead cell's far node.
function _nodes(topo::MeshTopology{1}, r::Int, dead::Int)
    rn, dn = topo.cells[r].nodes, topo.cells[dead].nodes
    shared = rn[1] in dn ? rn[1] : rn[2]
    inner = rn[1] == shared ? rn[2] : rn[1]
    far = dn[1] == shared ? dn[2] : dn[1]
    return inner, shared, far
end

# The volume between nodes `a` and `b` at the positions of `u`.
function _volume(u, layout::StateLayout, topo::MeshTopology{1}, a::Int, b::Int)
    A = has_chi(layout) ? u[area_index(layout)] : topo.cross_section
    xa, xb = u[position_index(layout, a)], u[position_index(layout, b)]
    return shell_volume(topo.coordinates, min(xa, xb), max(xa, xb), A)
end

# Merge the cells `dead`, which lie outward of `r` and next to it, into `r`.
function _merge!(
    u,
    d::CellDeactivation{P,NC},
    ::ConservativeTransfer,
    def::CoreDefinition,
    r::Int,
    dead::UnitRange{Int},
    t::Float64,
    invocation::Int,
) where {P,NC}
    layout, setup = def.layout, def.setup
    laws, topo, material = setup.laws, setup.topology, setup.material
    energy = energy_state(layout)
    moving = is_moving(layout)
    outermost_first = reverse(dead)
    fractions = [condensed_mass(d, layout, u, i) / d.reference[i] for i in outermost_first]
    m_r = _masses(u, layout, r)
    T_r = _temperature(energy, u, d, laws, layout, r, m_r)
    E = _energy(energy, u, d, laws, layout, r, m_r, T_r)
    # The merged temperature carries the summed sensible enthalpy of what the
    # receiver takes; on a fixed mesh the dead cells' gas is vented. On a
    # moving mesh the volumes and the pore volumes are summed, Σⱼ mⱼ/ρⱼ(T)
    # being the skeletal volume of the masses.
    H = _sensible_enthalpy(laws, T_r, m_r)
    m, vented, H_vented = m_r, ntuple(j -> 0.0, Val(NC)), 0.0
    skeletal(mm, T) = skeletal_volume_fraction(mm, T, material)
    inner, shared, _ = _nodes(topo, r, first(dead))
    V = moving ? _volume(u, layout, topo, inner, shared) : 0.0
    pores = moving ? V - skeletal(m_r, T_r) : 0.0
    for i in outermost_first
        m_d = _masses(u, layout, i)
        T_d = _temperature(energy, u, d, laws, layout, i, m_d)
        E += _energy(energy, u, d, laws, layout, i, m_d, T_d)
        v =
            moving ? ntuple(j -> 0.0, Val(NC)) :
            ntuple(j -> d.solid[j] ? 0.0 : m_d[j], Val(NC))
        H_v = _sensible_enthalpy(laws, T_d, v)
        m = ntuple(j -> m[j] + (m_d[j] - v[j]), Val(NC))
        H += _sensible_enthalpy(laws, T_d, m_d) - H_v
        vented = map(+, vented, v)
        H_vented += H_v
        if moving
            _, near, far = _nodes(topo, i - 1, i)
            V_d = _volume(u, layout, topo, near, far)
            V += V_d
            pores += V_d - skeletal(m_d, T_d)
        end
    end
    T_m = first(invert_caloric(d.caloric, laws, H, m, r, T_r))
    created =
        moving ?
        _enclose!(u, layout, topo, inner, shared, dead, max(V, skeletal(m, T_m) + pores)) -
        V : 0.0
    for j = 1:NC
        u[mass_index(layout, r, j)] = m[j]
    end
    E_m = extended_energy(d.caloric, laws, T_m, m, r)
    u[energy_index(layout, r)] = energy isa EnthalpyState ? E_m : T_m
    for (k, i) in enumerate(outermost_first)
        for j = 1:NC
            u[mass_index(layout, i, j)] = 0
        end
        u[energy_index(layout, i)] = energy isa EnthalpyState ? 0.0 : T_m
        if has_chi(layout)
            for h = 1:2
                u[history_index(layout, r, h)] += u[history_index(layout, i, h)]
                u[history_index(layout, i, h)] = 0
            end
        end
        _, _, far = _nodes(topo, i - 1, i)
        first_event = k == 1
        push!(
            d.events,
            DeactivationEvent(
                t,
                invocation,
                i,
                r,
                fractions[k],
                _mask!(layout, i, far),
                first_event ? created : 0.0,
                first_event ? collect(Float64, vented) : zeros(NC),
                first_event ? H_vented : 0.0,
                first_event ? E - E_m - H_vented : 0.0,
            ),
        )
    end
    return nothing
end

# The receiver's far node moves out to enclose the volume `V`, every far
# node of the dead cells onto it; returns `V`.
function _enclose!(u, layout, topo, inner, shared, dead, V)
    A = has_chi(layout) ? u[area_index(layout)] : topo.cross_section
    x = outer_position(topo.coordinates, u[position_index(layout, inner)], V, A)
    u[position_index(layout, shared)] = x
    for i in dead
        _, _, far = _nodes(topo, i - 1, i)
        u[position_index(layout, far)] = x
        u[velocity_index(layout, far)] = 0
    end
    return V
end

# The last live cell goes with no receiver: its inventories are written to
# the event as lost.
function _lose!(
    u,
    d::CellDeactivation{P,NC},
    def::CoreDefinition,
    i::Int,
    t::Float64,
    invocation::Int,
) where {P,NC}
    layout, setup = def.layout, def.setup
    laws = setup.laws
    energy = energy_state(layout)
    fraction = condensed_mass(d, layout, u, i) / d.reference[i]
    m = _masses(u, layout, i)
    T = _temperature(energy, u, d, laws, layout, i, m)
    H = _sensible_enthalpy(laws, T, m)
    E = _energy(energy, u, d, laws, layout, i, m, T)
    for j = 1:NC
        u[mass_index(layout, i, j)] = 0
    end
    energy isa EnthalpyState && (u[energy_index(layout, i)] = 0)
    if has_chi(layout)
        for k = 1:2
            u[history_index(layout, i, k)] = 0
        end
    end
    rows = _mask!(layout, i, maximum(setup.topology.cells[i].nodes))
    push!(
        d.events,
        DeactivationEvent(t, invocation, i, 0, fraction, rows, 0.0, collect(m), H, E - H),
    )
    return nothing
end

# Mask the rows of cell `i` and of node `far`; returns them.
function _mask!(layout::StateLayout, i::Int, far::Int)
    rows = vcat(collect(cell_rows(layout, i)), collect(node_rows(layout, far)))
    layout.mask[rows] .= false
    return rows
end
