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
the live cell beside them, the receiver, whose temperature follows from the
summed sensible enthalpy of all of them; the energy state is set from it. The
caloric extension's regularization, which only a cell far below the
heat-capacity floor carries, is not transferred: what the cells held of it
beyond the receiver's after the merge is written to the event. The receiver is
the live cell inward of the dead run, or, for a run from the innermost cell
under live cells, the live cell outward of it. On a moving mesh the receiver
also absorbs the dead cells' volume, and the merge never destroys pore space:
the merged volume is the larger of the summed volumes and the merged skeletal
volume plus the summed pre-merge pore volumes, the difference written to the
event as created volume. Inward, the receiver's far face moves out to enclose
the merged volume, where every node of the dead cells then coincides (the
moving-boundary collapse); outward, the merged volume is what lies between the
innermost node and the receiver's far face. In both the created volume moves
that face and every live node beyond it, each live cell keeping its volume. On
a fixed mesh the receiver keeps its volume and takes the condensed inventories
only: its pores could hold the dead cell's gas only at a higher pressure, so
that gas is vented, its masses and sensible enthalpy written to the event as
lost.
"""
struct ConservativeTransfer <: DeactivationPolicy end

"""
	DeactivationEvent

One deactivated cell: its time, the callback invocation it belongs to
(numbered from one over the solve, shared by every cell one invocation
deactivates), the dead cell and the receiver as cells of the initial mesh (the
receiver zero when none was left, the dead cell's inventories then lost), the
dead cell's condensed mass as a fraction of its reference when it went, the
rows masked, the outermost slots once the live cells have moved down, the volume
the merge created
[m³], the masses [kg] and sensible enthalpy [J] lost, the caloric
regularization [J] the cells held that the receiver does not (the extended
energies of the merged cells before the event less the receiver's after it
and the enthalpy lost), the initial cell in each slot after the event
([`SlotTables`](@ref)), and the reference updates of the merge, each a
reaction, an initial cell and that cell's new reference, and whether live cells
lay outward of the dead cell (`interior`): the cell left from under them,
through a void of the tail's making ([`deactivate!`](@ref)). A merge of several
cells into one receiver records one event per dead cell, the outermost first;
the merge's created volume, losses, regularization and reference updates are
written to that first event and the others carry zeros, so the sums over the
events are the merges'.
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
    origin::Vector{Int}
    reference_updates::Vector{Tuple{Int,Int,Float64}}
    interior::Bool
end

"""
	CellDeactivation{P,NC,E}

The deactivation state of a solve: the policy, the threshold `ε_m`, the
near-threshold band `ε_band` (a fraction of a cell's reference, read on a
moving mesh only, see [`deactivate!`](@ref)), the cells' initial condensed
masses (`reference`; a receiver keeps its own, so a merged cell deactivates
once its whole merged inventory is below `ε_m` of what it started with), the
condensed component flags, the caloric extension merged temperatures are
found with, the events so far, the outcome: `:running`; `:burned_out` once the last live cell
went; `:interior` when a cell below the threshold lies inward of a live one on
a fixed mesh, whose geometry stays with its slots, or beneath a skin on a moving
one; the [`SlotTables`](@ref) of the solve; and, at a reaction order other than one,
what the threshold reads there ([`below_threshold`](@ref)): the components
`reactants`, one per reaction with a condensed reactant, the cells' `reactive`
reference, and `lineage`, the reaction kernel's references by slot.
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
    tables::SlotTables
    reactants::Vector{Int}
    reactive::Vector{Float64}
    lineage::Matrix{Float64}
end

"""
	CellDeactivation(policy, threshold, band, reference, solid, caloric;
	                 origin = collect(eachindex(reference)), reactants = Int[],
	                 reactive = Float64[], lineage = zeros(0, 0))

The deactivation state before a solve, `reference` the initial condensed
masses; `threshold ≤ band < 1`. `origin` is the array of the
[`SlotTables`](@ref), with no table registered yet.
"""
function CellDeactivation(
    policy::DeactivationPolicy,
    threshold::Real,
    band::Real,
    reference::Vector{Float64},
    solid::NTuple{NC,Bool},
    caloric::CaloricExtension{NC};
    origin::Vector{Int} = collect(eachindex(reference)),
    reactants::Vector{Int} = Int[],
    reactive::Vector{Float64} = Float64[],
    lineage::Matrix{Float64} = zeros(0, 0),
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
        SlotTables(origin),
        reactants,
        reactive,
        lineage,
    )
end

"""
	reset_cells!(d, layout) -> d

Every row live, every cell in its initial slot with the tables as set up, and
no events: the state a solve starts from.
"""
function reset_cells!(d::CellDeactivation, layout::StateLayout)
    fill!(layout.mask, true)
    empty!(d.events)
    d.outcome[] = :running
    arrange_slots!(d.tables, layout, 1:n_cells(layout), Tuple{Int,Int,Float64}[])
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
	reactive_mass(d, layout, u, i)

The reactive inventory of cell `i` in the state `u` [kg]: the mass of each
reaction's condensed reactant, summed over the reactions.
"""
@inline function reactive_mass(
    d::CellDeactivation,
    layout::StateLayout,
    u::AbstractVector,
    i::Int,
)
    m = zero(eltype(u))
    @inbounds for j in d.reactants
        m += u[mass_index(layout, i, j)]
    end
    return m
end

# Whether cell `i` is the last live one.
@inline _last_live(layout::StateLayout, i::Int) =
    i == 1 && (n_cells(layout) == 1 || !layout.mask[energy_index(layout, 2)])

# What the threshold reads in live cell `i`: its inventory and that inventory's
# reference.
@inline function _measure(d::CellDeactivation, layout::StateLayout, u, i::Int)
    (isempty(d.reactants) || _last_live(layout, i)) &&
        return condensed_mass(d, layout, u, i), d.reference[i]
    return reactive_mass(d, layout, u, i), d.reactive[i]
end

# Whether live cell `i` holds less than `1 + slack` times the fraction `f` of
# its reference.
@inline function _below(d::CellDeactivation, layout::StateLayout, u, i::Int, f, slack)
    layout.mask[energy_index(layout, i)] || return false
    held, reference = _measure(d, layout, u, i)
    reference > 0 || return false
    return held < f * reference * (1 + slack)
end

"""
	below_threshold(d, layout, u, i; slack = 0) -> Bool

Whether live cell `i` holds less than `ε_m` of its reference, or, with a positive
`slack`, less than `(1 + slack) ε_m` of it. At order one the test reads the
condensed mass against the initial condensed mass. At an order other than one
it reads the reactive inventory ([`reactive_mass`](@ref)) against the cell's
`reactive` reference, the sum over the reactions with a condensed reactant of the
reaction's lineage reference, or for a first-order reaction its reactant's initial
mass, so that a residue never holds a spent cell live; the last live cell has no
receiver and goes only below the threshold on its condensed mass, so a column
that leaves a residue ends with that residue. A cell whose reference is zero,
which started with no condensed mass or with none of a reaction's lineage, never
deactivates.
"""
@inline below_threshold(
    d::CellDeactivation,
    layout::StateLayout,
    u,
    i::Int;
    slack::Float64 = 0.0,
) = _below(d, layout, u, i, d.threshold, slack)

"""
	any_below(d, layout, u; slack = 0) -> Bool

Whether any live cell is below the threshold ([`below_threshold`](@ref)): the
end-of-step event's condition.
"""
function any_below(
    d::CellDeactivation,
    layout::StateLayout,
    u::AbstractVector;
    slack::Float64 = 0.0,
)
    for i = 1:n_cells(layout)
        below_threshold(d, layout, u, i; slack) && return true
    end
    return false
end

"""
	threshold_margin(d, layout, u) -> Float64

The smallest margin of a live cell's inventory `m` over its threshold
`ε_m m_ref`, `(m − ε_m m_ref)/(ε_m m_ref)`, in the measure of
[`below_threshold`](@ref); `Inf` when no live cell is tested. The margin is zero
exactly when the cell holds exactly its threshold and negative exactly when it
is below it. Its first root in a step is the located event's time.
"""
function threshold_margin(d::CellDeactivation, layout::StateLayout, u::AbstractVector)
    margin = Inf
    for i = 1:n_cells(layout)
        layout.mask[energy_index(layout, i)] || continue
        held, reference = _measure(d, layout, u, i)
        reference > 0 || continue
        limit = d.threshold * reference
        margin = min(margin, (held - limit) / limit)
    end
    return margin
end

"""
	top_live_cell(layout) -> Int

The outermost live cell, zero when none is live. The dead cells leave the
outermost slots, so the live cells are `1:top_live_cell(layout)`.
"""
function top_live_cell(layout::StateLayout)
    for i = n_cells(layout):-1:1
        layout.mask[energy_index(layout, i)] && return i
    end
    return 0
end

# The outermost live cell below the threshold, zero when none.
function _outermost_below(d::CellDeactivation, layout::StateLayout, u, top::Int, slack)
    for i = top:-1:1
        below_threshold(d, layout, u, i; slack) && return i
    end
    return 0
end

"""
	VOID_SKIN

How many times its threshold a live cell outward of a crossing must hold for
the crossing to be a void beneath a skin. A tail delays its receiver past the
cell inward of it, which then crosses first under a receiver that holds
`ε (e^{kτ} + 1)`: a layer burning out that way holds up to six thresholds on the
decks measured, a skin from 29 to thousands of them.
"""
const VOID_SKIN = 10.0

# Whether a live cell outward of slot `hi`, up to `top`, is a skin: it holds
# `VOID_SKIN` thresholds of its reference, or has no reference, an inert layer or
# one that keeps its volume.
function _skin_above(d::CellDeactivation, layout::StateLayout, u, hi::Int, top::Int)
    for i = (hi+1):top
        held, reference = _measure(d, layout, u, i)
        (reference > 0 && held < VOID_SKIN * d.threshold * reference) || return true
    end
    return false
end

# The contiguous live cells around `c` that hold less than `1 + slack` times
# `ε_band` of their reference.
function _band_run(d::CellDeactivation, layout::StateLayout, u, c::Int, top::Int, slack)
    lo, hi = c, c
    while hi < top && _below(d, layout, u, hi + 1, d.band, slack)
        hi += 1
    end
    while lo > 1 && _below(d, layout, u, lo - 1, d.band, slack)
        lo -= 1
    end
    return lo, hi
end

"""
	deactivate!(u, d, def, t; slack = 0) -> Symbol

Deactivate, in one pass, the cells below the threshold. On a moving mesh the
outermost cell below `ε_m` goes, with the contiguous live cells around it below
`ε_band` of their reference, in one merge ([`ConservativeTransfer`](@ref)) into
the live cell beneath the run, or, for a run from the innermost cell under live
cells, into the live cell above it; a run of the whole column leaves the
innermost cell, holding the run's inventories, as the survivor. A run with live
cells outward of it goes only if they are spent, each under [`VOID_SKIN`](@ref)
thresholds: a tail delayed their receiver, and the void is the transfer's. A skin
over the run is a void the model cannot carry heat across, so the pass stops
there. The live cells above a run that went move down into its slots, with
every table of the [`SlotTables`](@ref), so that the dead slots are the outermost
ones. After a merge the pass starts again, the receiver holding the run's
inventories against its own reference, so that every cell below `ε_m` goes in the
same invocation. On a fixed mesh, whose geometry stays with its slots, only the
outermost live cell can go, while it is below `ε_m`, each merge re-testing the
receiver with the dead cell's condensed inventories added. The last live cell,
once below `ε_m`, goes with no receiver, its inventories written to the event as
lost. The events of one call share their invocation. A positive `slack` widens
both tests to `1 + slack` times `ε_m` and `ε_band`: an event located at a root
passes the solve's relative tolerance, so that the crossing cell, which holds
exactly its threshold there, goes, and with it the cells within the tolerance of
their thresholds, whose crossing the solve does not resolve from the located one.

Returns `:transferred`, `:burned_out` when the last live cell went, or `:interior`
when a cell below the threshold remains inward of a live one on a fixed mesh or
beneath a skin on a moving one.
"""
function deactivate!(
    u::AbstractVector,
    d::CellDeactivation,
    def::CoreDefinition,
    t::Real;
    slack::Float64 = 0.0,
)
    layout = def.layout
    invocation = isempty(d.events) ? 1 : last(d.events).invocation + 1
    while true
        top = top_live_cell(layout)
        top == 0 && break
        c = _outermost_below(d, layout, u, top, slack)
        c == 0 && break
        is_moving(layout) || c == top || break
        lo, hi = is_moving(layout) ? _band_run(d, layout, u, c, top, slack) : (c, c)
        _skin_above(d, layout, u, hi, top) && break
        if top == 1
            _lose!(u, d, def, top, Float64(t), invocation)
            d.outcome[] = :burned_out
            return :burned_out
        end
        lo == 1 && hi == top && (lo = 2)
        r = lo == 1 ? hi + 1 : lo - 1
        _merge!(u, d, d.policy, def, r, lo:hi, Float64(t), invocation)
        if hi < top
            _compact!(u, d, def, lo, hi, top)
        else
            arrange_slots!(d.tables, layout, d.tables.origin, d.tables.held)
        end
    end
    if any_below(d, layout, u; slack)
        d.outcome[] = :interior
        return :interior
    end
    return :transferred
end

"""
	replay_events!(d, layout, t) -> d

The mask and the tables of a state saved at `t`: every event of `d` at or before
`t` applied, starting from every row live and every cell in its initial slot.
The integrator applies an event before it saves the state of the step that
triggered it, so a state saved at an event's own time is the one after the
event.
"""
function replay_events!(d::CellDeactivation, layout::StateLayout, t::Real)
    fill!(layout.mask, true)
    origin = 1:n_cells(layout)
    updates = Tuple{Int,Int,Float64}[]
    for e in d.events
        e.t <= t || break
        layout.mask[e.rows] .= false
        origin = e.origin
        append!(updates, e.reference_updates)
    end
    arrange_slots!(d.tables, layout, origin, updates)
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

# The volume of cell `i` at the positions of `u`.
function _cell_volume(u, layout::StateLayout, topo::MeshTopology{1}, i::Int)
    A = has_chi(layout) ? u[area_index(layout)] : topo.cross_section
    return shell_volume(
        topo.coordinates,
        u[position_index(layout, i)],
        u[position_index(layout, i + 1)],
        A,
    )
end

# Merge the cells `dead`, a run of slots next to `r`, into `r`: inward when `r`
# lies beneath them, outward when it lies above them.
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
    origin = d.tables.origin
    outermost_first = reverse(dead)
    interior = last(dead) < top_live_cell(layout)
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
    V = moving ? _cell_volume(u, layout, topo, r) : 0.0
    pores = moving ? V - skeletal(m_r, T_r) : 0.0
    T_lo, T_hi = T_r, T_r
    for i in outermost_first
        m_d = _masses(u, layout, i)
        T_d = _temperature(energy, u, d, laws, layout, i, m_d)
        T_lo, T_hi = min(T_lo, T_d), max(T_hi, T_d)
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
            V_d = _cell_volume(u, layout, topo, i)
            V += V_d
            pores += V_d - skeletal(m_d, T_d)
        end
    end
    T_m = merged_temperature(d.caloric, laws, H, m, r, T_r, T_lo, T_hi)
    created =
        moving ?
        _enclose!(u, layout, topo, r, dead, max(V, skeletal(m, T_m) + pores), V) - V : 0.0
    updates = _receive_lineage!(d, r, dead)
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
        if has_history(layout)
            for h = 1:2
                u[history_index(layout, r, h)] += u[history_index(layout, i, h)]
                u[history_index(layout, i, h)] = 0
            end
        end
        first_event = k == 1
        push!(
            d.events,
            DeactivationEvent(
                t,
                invocation,
                origin[i],
                origin[r],
                fractions[k],
                _mask!(layout, i, i + 1),
                first_event ? created : 0.0,
                first_event ? collect(Float64, vented) : zeros(NC),
                first_event ? H_vented : 0.0,
                first_event ? E - E_m - H_vented : 0.0,
                copy(origin),
                first_event ? updates : Tuple{Int,Int,Float64}[],
                interior,
            ),
        )
    end
    return nothing
end

# A receiver with no lineage of a reaction whose lineage the dead cells hold
# takes theirs: it then holds none of that lineage but what it received, which
# keeps the dead cells' conversion. Returns the updates, written to the live and
# the initial reference tables.
function _receive_lineage!(d::CellDeactivation, r::Int, dead::UnitRange{Int})
    updates = Tuple{Int,Int,Float64}[]
    isempty(d.reactants) && return updates
    tables = d.tables
    _, initial = only(tables.lineage)
    for k in axes(d.lineage, 1)
        iszero(d.lineage[k, r]) || continue
        value = sum(d.lineage[k, i] for i in dead)
        value > 0 || continue
        cell = tables.origin[r]
        d.lineage[k, r] = value
        initial[k, cell] = value
        push!(updates, (k, cell, value))
        push!(tables.held, (k, cell, value))
    end
    return updates
end

# The merged volume `merged` of the receiver `r` and the cells `dead`, which
# held `V`; returns `merged`. Inward the receiver's far node moves out to
# enclose it, every far node of the dead cells onto it. Outward the receiver's
# far node stays unless the merge creates volume. Created volume also carries
# every live node beyond, each live cell keeping its volume.
function _enclose!(u, layout, topo, r::Int, dead::UnitRange{Int}, merged, V)
    A = has_chi(layout) ? u[area_index(layout)] : topo.cross_section
    beyond = (max(r, last(dead))+1):top_live_cell(layout)
    carried = merged > V ? [_cell_volume(u, layout, topo, i) for i in beyond] : Float64[]
    x = u[position_index(layout, r + 1)]
    if r < first(dead)
        x = outer_position(topo.coordinates, u[position_index(layout, r)], merged, A)
        u[position_index(layout, r + 1)] = x
        for i in dead
            u[position_index(layout, i + 1)] = x
            u[velocity_index(layout, i + 1)] = 0
        end
    elseif merged > V
        x = outer_position(topo.coordinates, u[position_index(layout, 1)], merged, A)
        u[position_index(layout, r + 1)] = x
    end
    for (i, V_i) in zip(beyond, carried)
        x = outer_position(topo.coordinates, x, V_i, A)
        u[position_index(layout, i + 1)] = x
    end
    return merged
end

# After the dead cells `lo:hi`, which live cells up to `top` lay outward of,
# merged: those cells move down into their slots, cell block and node block
# alike, and the tables with them, so that the dead slots are the outermost
# ones; the last events' rows and permutations are those of the result. The
# dead cells' nodes coincide with the receiver's far node, or lie inside the
# receiver when it is the cell that moves into the innermost slot.
function _compact!(
    u,
    d::CellDeactivation{P,NC},
    def::CoreDefinition,
    lo::Int,
    hi::Int,
    top::Int,
) where {P,NC}
    layout = def.layout
    energy = energy_state(layout)
    tables = d.tables
    origin = tables.origin
    nd = hi - lo + 1
    last_live = top - nd
    for c = (hi+1):top
        u[cell_rows(layout, c - nd)] .= view(u, cell_rows(layout, c))
    end
    for n = (hi+2):(top+1)
        u[node_rows(layout, n - nd)] .= view(u, node_rows(layout, n))
    end
    origin[lo:top] = vcat(origin[(hi+1):top], origin[lo:hi])
    x_top = u[position_index(layout, last_live + 1)]
    T_top = u[energy_index(layout, last_live)]
    for c = (last_live+1):top
        for j = 1:NC
            u[mass_index(layout, c, j)] = 0
        end
        u[energy_index(layout, c)] = energy isa EnthalpyState ? 0.0 : T_top
        if has_history(layout)
            for h = 1:2
                u[history_index(layout, c, h)] = 0
            end
        end
        u[node_rows(layout, c + 1)] .= 0
        u[position_index(layout, c + 1)] = x_top
    end
    for c = lo:last_live
        layout.mask[cell_rows(layout, c)] .= true
        layout.mask[node_rows(layout, c + 1)] .= true
    end
    for k = 1:nd
        slot = top - k + 1
        i = length(d.events) - nd + k
        e = d.events[i]
        d.events[i] = DeactivationEvent(
            e.t,
            e.invocation,
            e.cell,
            e.receiver,
            e.fraction,
            _mask!(layout, slot, slot + 1),
            e.created_volume,
            e.lost_mass,
            e.lost_enthalpy,
            e.regularization,
            copy(origin),
            e.reference_updates,
            e.interior,
        )
    end
    permute_slots!(tables, layout)
    return nothing
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
    if has_history(layout)
        for k = 1:2
            u[history_index(layout, i, k)] = 0
        end
    end
    rows = _mask!(layout, i, maximum(setup.topology.cells[i].nodes))
    push!(
        d.events,
        DeactivationEvent(
            t,
            invocation,
            d.tables.origin[i],
            0,
            fraction,
            rows,
            0.0,
            collect(m),
            H,
            E - H,
            copy(d.tables.origin),
            Tuple{Int,Int,Float64}[],
            false,
        ),
    )
    return nothing
end

# Mask the rows of cell `i` and of node `far`; returns them.
function _mask!(layout::StateLayout, i::Int, far::Int)
    rows = vcat(collect(cell_rows(layout, i)), collect(node_rows(layout, far)))
    layout.mask[rows] .= false
    return rows
end
