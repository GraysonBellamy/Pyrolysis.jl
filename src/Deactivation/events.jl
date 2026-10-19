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

The dead cell's inventories (masses, sensible enthalpy, history pair) go to
the live cell inward of it, the receiver, whose temperature follows from the
summed sensible enthalpy of both; the energy state is set from it (the
caloric extension's regularization, which only a cell far below the
heat-capacity floor carries, is not transferred). On a moving mesh the receiver also absorbs the dead cell's
volume: its far face moves to the dead cell's far face, where both of the dead
cell's nodes then coincide (the moving-boundary collapse), and the merge never
destroys pore space: the merged volume is the larger of the summed volumes
and the merged skeletal volume plus the summed pre-merge pore volumes, the
difference written to the event as created volume. On a fixed mesh the
receiver keeps its volume and takes the inventories only.
"""
struct ConservativeTransfer <: DeactivationPolicy end

"""
	DeactivationEvent

One deactivation: its time, the dead cell, the receiver (zero when none was
left, the dead cell's inventories then lost), the rows masked, the volume the
merge created [m³], and the masses [kg] and sensible enthalpy [J] lost.
"""
struct DeactivationEvent
    t::Float64
    cell::Int
    receiver::Int
    rows::Vector{Int}
    created_volume::Float64
    lost_mass::Vector{Float64}
    lost_enthalpy::Float64
end

"""
	CellDeactivation{P,NC,E}

The deactivation state of a solve: the policy, the threshold `ε_m`, the cells'
initial condensed masses (`reference`; a receiver keeps its own, so a merged
cell deactivates once its whole merged inventory is below `ε_m` of what it
started with; the reaction kernel's switch reads the same masses), the
condensed component flags, the caloric extension merged temperatures are
found with, the events so far, and the outcome: `:running`; `:burned_out` once
the last live cell went; `:interior` when a cell below the threshold lies
inward of a live one, which the policies here do not handle.
"""
struct CellDeactivation{P<:DeactivationPolicy,NC,E<:CaloricExtension{NC}}
    policy::P
    threshold::Float64
    reference::Vector{Float64}
    solid::NTuple{NC,Bool}
    caloric::E
    events::Vector{DeactivationEvent}
    outcome::Base.RefValue{Symbol}
end

"""
	CellDeactivation(policy, threshold, reference, solid, caloric)

The deactivation state before a solve, `reference` the initial condensed
masses.
"""
CellDeactivation(
    policy::DeactivationPolicy,
    threshold::Real,
    reference::Vector{Float64},
    solid::NTuple{NC,Bool},
    caloric::CaloricExtension{NC},
) where {NC} = CellDeactivation(
    policy,
    Float64(threshold),
    reference,
    solid,
    caloric,
    DeactivationEvent[],
    Ref(:running),
)

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

Whether live cell `i` holds less than `ε_m` of its initial condensed mass.
"""
@inline below_threshold(d::CellDeactivation, layout::StateLayout, u, i::Int) =
    layout.mask[energy_index(layout, i)] &&
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

"""
	deactivate!(u, d, def, t) -> Symbol

Deactivate, in one pass, every cell below the threshold from the exposed
side: the outermost live cell while it is below, each transfer re-testing the
receiver with the dead cell's inventory added. Returns
`:transferred`, `:burned_out` when the last live cell went (its inventories
written to the event as lost), or `:interior` when a cell below the threshold
remains inward of a live one.
"""
function deactivate!(u::AbstractVector, d::CellDeactivation, def::CoreDefinition, t::Real)
    layout = def.layout
    while true
        top = top_live_cell(layout)
        (top == 0 || !below_threshold(d, layout, u, top)) && break
        if top == 1
            _lose!(u, d, def, top, Float64(t))
            d.outcome[] = :burned_out
            return :burned_out
        end
        _transfer!(u, d, d.policy, def, top, top - 1, Float64(t))
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

# The receiver's and the dead cell's nodes: the receiver's inward node, the
# node they share and the dead cell's far node.
function _nodes(topo::MeshTopology{1}, r::Int, dead::Int)
    rn, dn = topo.cells[r].nodes, topo.cells[dead].nodes
    shared = rn[1] in dn ? rn[1] : rn[2]
    inner = rn[1] == shared ? rn[2] : rn[1]
    far = dn[1] == shared ? dn[2] : dn[1]
    return inner, shared, far
end

function _transfer!(
    u,
    d::CellDeactivation{P,NC},
    ::ConservativeTransfer,
    def::CoreDefinition,
    dead::Int,
    r::Int,
    t::Float64,
) where {P,NC}
    layout, setup = def.layout, def.setup
    laws = setup.laws
    energy = energy_state(layout)
    m_d, m_r = _masses(u, layout, dead), _masses(u, layout, r)
    m = map(+, m_r, m_d)
    T_d = _temperature(energy, u, d, laws, layout, dead, m_d)
    T_r = _temperature(energy, u, d, laws, layout, r, m_r)
    # The merged temperature carries the summed sensible enthalpy.
    H = _sensible_enthalpy(laws, T_r, m_r) + _sensible_enthalpy(laws, T_d, m_d)
    T_m = first(invert_caloric(d.caloric, laws, H, m, r, T_r))
    inner, shared, far = _nodes(setup.topology, r, dead)
    created =
        is_moving(layout) ?
        _merge_volume!(u, def, inner, shared, far, m_d, m_r, m, T_d, T_r, T_m) : 0.0
    for j = 1:NC
        u[mass_index(layout, r, j)] = m[j]
        u[mass_index(layout, dead, j)] = 0
    end
    u[energy_index(layout, r)] =
        energy isa EnthalpyState ? extended_energy(d.caloric, laws, T_m, m, r) : T_m
    u[energy_index(layout, dead)] = energy isa EnthalpyState ? 0.0 : T_m
    if has_chi(layout)
        for k = 1:2
            u[history_index(layout, r, k)] += u[history_index(layout, dead, k)]
            u[history_index(layout, dead, k)] = 0
        end
    end
    rows = _mask!(layout, dead, far)
    push!(d.events, DeactivationEvent(t, dead, r, rows, created, zeros(NC), 0.0))
    return nothing
end

# The receiver's far node moves out to enclose the merged volume, the dead
# cell's far node onto it; returns the created volume.
function _merge_volume!(u, def, inner, shared, far, m_d, m_r, m, T_d, T_r, T_m)
    layout, setup = def.layout, def.setup
    topo, material = setup.topology, setup.material
    coords = topo.coordinates
    A = has_chi(layout) ? u[area_index(layout)] : topo.cross_section
    x0 = u[position_index(layout, inner)]
    x1 = u[position_index(layout, shared)]
    x2 = u[position_index(layout, far)]
    V_r = shell_volume(coords, min(x0, x1), max(x0, x1), A)
    V_d = shell_volume(coords, min(x1, x2), max(x1, x2), A)
    # Σⱼ mⱼ/ρⱼ(T): the skeletal volume of the masses.
    skeletal(mm, T) = skeletal_volume_fraction(mm, T, material)
    pores = (V_r - skeletal(m_r, T_r)) + (V_d - skeletal(m_d, T_d))
    V = max(V_r + V_d, skeletal(m, T_m) + pores)
    x = outer_position(coords, x0, V, A)
    u[position_index(layout, shared)] = x
    u[position_index(layout, far)] = x
    u[velocity_index(layout, far)] = 0
    return V - V_r - V_d
end

# The last live cell goes with no receiver: its inventories are written to
# the event as lost.
function _lose!(
    u,
    d::CellDeactivation{P,NC},
    def::CoreDefinition,
    i::Int,
    t::Float64,
) where {P,NC}
    layout, setup = def.layout, def.setup
    laws = setup.laws
    energy = energy_state(layout)
    m = _masses(u, layout, i)
    T = _temperature(energy, u, d, laws, layout, i, m)
    H = _sensible_enthalpy(laws, T, m)
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
    push!(d.events, DeactivationEvent(t, i, 0, rows, 0.0, collect(m), H))
    return nothing
end

# Mask the rows of cell `i` and of node `far`; returns them.
function _mask!(layout::StateLayout, i::Int, far::Int)
    rows = vcat(collect(cell_rows(layout, i)), collect(node_rows(layout, far)))
    layout.mask[rows] .= false
    return rows
end
