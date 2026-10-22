# ============================================================================
# Energy state
# ============================================================================

"""
	EnergyState

Which energy variable the integrator carries per cell. The state mapping
([`derive!`](@ref)) and the row assembly are the only code that knows which
one is in use.
"""
abstract type EnergyState end

"""
	TemperatureState <: EnergyState

The cell temperature `Tᵢ`. Its row is `(Σⱼ mᵢⱼ c_p,ⱼ) dTᵢ/dt = Eᵢ − Σⱼ hⱼ(Tᵢ) Ṁᵢⱼ`,
with `Eᵢ` the energy right-hand side and `Ṁᵢⱼ` the face-flux part of the mass
rates; energy is conserved to the integrator's order.
"""
struct TemperatureState <: EnergyState end

"""
	EnthalpyState <: EnergyState

The cell's sensible-enthalpy inventory `Hᵢ = Σⱼ mᵢⱼ hⱼ(Tᵢ)`. Its row is
`dHᵢ/dt = Eᵢ + Σⱼ ωᵢⱼ Vᵢ hⱼ(Tᵢ)`, a linear invariant every Runge–Kutta step
conserves to round-off; the temperature is recovered by the regularized
inversion of [`CaloricExtension`](@ref).
"""
struct EnthalpyState <: EnergyState end

# ============================================================================
# State layout
# ============================================================================

"""
	StateLayout{NC,E,Mov,Rad,Chi,D,U,Car}

Position of every unknown in the state vector of a `D`-dimensional mesh with
`NC` components and energy state `E`.

The state is a sequence of node blocks and cell blocks, then the global
scalars, in the order of [`block_order`](@ref): on a line, node 1, cell 1,
node 2, …, cell N, node N + 1, so the Jacobian is banded in every geometry;
in more dimensions each cell's block follows those of its nodes not yet
placed. A node block holds, in order, the `D` coordinates of the node
position `x` (moving mesh, `Mov`), the `D` components of its velocity `w`
when the closure's unknowns are the node velocities (`U = :velocities`), and
the intensity `I` of the face at that node (Beer–Lambert, `Rad`, 1D). A cell
block holds the component masses `m₁ … m_NC`, the energy state, the fields
the solid carries (`Car`, in order: `:history`, the conversion-history pair
`G`, `M`, which the variable cross section `Chi` needs; `:deformation`, the
`D × D` deformation gradient `F` column by column), and the closure's
multiplier `φ` when the closure's unknowns are one multiplier per cell
(`U = :multipliers`). The global block holds the cross-section measure `A`
and its lateral strain rate `θ_A` (`Chi`).

Kernels never compute an index; every access goes through the index
functions ([`mass_index`](@ref), [`energy_index`](@ref), …). `mask` marks the
live rows, a dead row being an identity row of the system; `differential`
marks the rows with a time derivative (see [`is_differential`](@ref)).
"""
struct StateLayout{NC,E<:EnergyState,Mov,Rad,Chi,D,U,Car}
    n_cells::Int
    n_nodes::Int
    node_offset::Vector{Int}
    cell_offset::Vector{Int}
    global_offset::Int
    mask::BitVector
    differential::BitVector
end

"""
	block_order(topology) -> Vector{Tuple{Symbol,Int}}

The order of the node and cell blocks in the state vector, as `(:node, n)`
and `(:cell, i)` entries: on a line node `i` then cell `i`, and the last node
after the last cell; in more dimensions, cell by cell, the nodes of each cell
not yet placed, in its node order, then the cell.
"""
function block_order(topo::MeshTopology{1})
    n_nodes(topo) == n_cells(topo) + 1 &&
    all(i -> topo.cells[i].nodes == (i, i + 1), 1:n_cells(topo)) ||
        throw(ArgumentError("a line's cells must lie between consecutive nodes"))
    return _line_order(n_cells(topo))
end

function block_order(topo::MeshTopology)
    order = Tuple{Symbol,Int}[]
    placed = falses(n_nodes(topo))
    for i = 1:n_cells(topo)
        for n in topo.cells[i].nodes
            placed[n] && continue
            placed[n] = true
            push!(order, (:node, n))
        end
        push!(order, (:cell, i))
    end
    all(placed) || throw(ArgumentError("every node must belong to a cell"))
    return order
end

# The interleaved block order of a line of `n` cells.
_line_order(n::Int) =
    push!([b for i = 1:n for b in ((:node, i), (:cell, i))], (:node, n + 1))

"""
	StateLayout(NC, n_cells, energy; moving = false, intensities = false, chi = false)

The layout for `n_cells` cells of a line with every row live.
"""
StateLayout(
    NC::Int,
    n_cells::Int,
    energy::EnergyState;
    moving::Bool = false,
    intensities::Bool = false,
    chi::Bool = false,
) = _layout(
    NC,
    n_cells,
    n_cells + 1,
    _line_order(n_cells),
    energy,
    1,
    :velocities;
    moving,
    intensities,
    chi,
    history = chi,
    deformation = false,
)

"""
	StateLayout(NC, topology, energy; moving = false, intensities = false, chi = false,
	            history = chi, deformation = false, unknowns = :velocities)

The layout for the cells and nodes of `topology` with every row live, its
blocks in [`block_order`](@ref). `unknowns` names the closure's algebraic
unknowns on a moving mesh: `:velocities`, `D` per node, or `:multipliers`,
one per cell. `history` carries the conversion-history pair, which the
variable cross section needs, and `deformation` the deformation gradient of
a moving mesh. Face intensities and the variable cross section are 1D only.
"""
function StateLayout(
    NC::Int,
    topo::MeshTopology{D},
    energy::EnergyState;
    moving::Bool = false,
    intensities::Bool = false,
    chi::Bool = false,
    history::Bool = chi,
    deformation::Bool = false,
    unknowns::Symbol = :velocities,
) where {D}
    D == 1 ||
        !(intensities || chi) ||
        throw(ArgumentError("face intensities and the variable cross section are 1D only"))
    return _layout(
        NC,
        n_cells(topo),
        n_nodes(topo),
        block_order(topo),
        energy,
        D,
        unknowns;
        moving,
        intensities,
        chi,
        history,
        deformation,
    )
end

function _layout(
    NC::Int,
    n_cells::Int,
    n_nodes::Int,
    order,
    energy::EnergyState,
    D::Int,
    unknowns::Symbol;
    moving::Bool,
    intensities::Bool,
    chi::Bool,
    history::Bool,
    deformation::Bool,
)
    NC >= 1 || throw(ArgumentError("a layout needs at least one component"))
    n_cells >= 1 || throw(ArgumentError("a layout needs at least one cell"))
    chi && !moving && throw(ArgumentError("the variable cross section needs a moving mesh"))
    chi &&
        !history &&
        throw(ArgumentError("the variable cross section needs the history pair"))
    deformation &&
        !moving &&
        throw(ArgumentError("the deformation gradient needs a moving mesh"))
    unknowns in (:velocities, :multipliers) ||
        throw(ArgumentError("unknowns must be :velocities or :multipliers, got :$unknowns"))
    carried = ((history ? (:history,) : ())..., (deformation ? (:deformation,) : ())...)
    L = StateLayout{NC,typeof(energy),moving,intensities,chi,D,unknowns,carried}
    nw, cw = _node_width(L), _cell_width(L)
    node_offset = zeros(Int, n_nodes)
    cell_offset = zeros(Int, n_cells)
    k = 0
    for (kind, id) in order
        if kind === :node
            node_offset[id] = k
            k += nw
        else
            cell_offset[id] = k
            k += cw
        end
    end
    len = k + _global_width(L)
    layout = L(n_cells, n_nodes, node_offset, cell_offset, k, trues(len), trues(len))
    layout.differential .= _differential_rows(layout)
    return layout
end

@inline n_components(::StateLayout{NC}) where {NC} = NC
@inline energy_state(::StateLayout{NC,E}) where {NC,E} = E()
@inline is_moving(::StateLayout{NC,E,Mov}) where {NC,E,Mov} = Mov
@inline has_intensities(::StateLayout{NC,E,Mov,Rad}) where {NC,E,Mov,Rad} = Rad
@inline has_chi(::StateLayout{NC,E,Mov,Rad,Chi}) where {NC,E,Mov,Rad,Chi} = Chi
@inline n_cells(l::StateLayout) = l.n_cells
@inline n_nodes(l::StateLayout) = l.n_nodes

"""
	dimension(layout) -> Int

The dimension `D` of the mesh the layout describes.
"""
@inline dimension(::StateLayout{NC,E,Mov,Rad,Chi,D}) where {NC,E,Mov,Rad,Chi,D} = D

"""
	carried_fields(layout) -> Tuple

The fields the solid carries per cell, in their order in the cell block.
"""
@inline carried_fields(
    ::StateLayout{NC,E,Mov,Rad,Chi,D,U,Car},
) where {NC,E,Mov,Rad,Chi,D,U,Car} = Car

"""
	has_history(layout), has_deformation(layout) -> Bool

Whether the cells carry the conversion-history pair, the deformation gradient.
"""
@inline has_history(l::StateLayout) = _carries(carried_fields(l), :history)
@inline has_deformation(l::StateLayout) = _carries(carried_fields(l), :deformation)

@inline _carries(::Tuple{}, f::Symbol) = false
@inline _carries(c::Tuple, f::Symbol) = first(c) === f || _carries(Base.tail(c), f)

# The rows of a carried field per cell, and the rows before field `f`.
@inline _field_width(f::Symbol, D::Int) = f === :history ? 2 : D * D
@inline _carried_width(::Tuple{}, D::Int) = 0
@inline _carried_width(c::Tuple, D::Int) =
    _field_width(first(c), D) + _carried_width(Base.tail(c), D)
@inline _carried_offset(c::Tuple, D::Int, f::Symbol) =
    first(c) === f ? 0 : _field_width(first(c), D) + _carried_offset(Base.tail(c), D, f)

"""
	closure_unknowns(layout) -> Symbol

The closure's algebraic unknowns: `:velocities` (`D` per node) or
`:multipliers` (one per cell).
"""
@inline closure_unknowns(::StateLayout{NC,E,Mov,Rad,Chi,D,U}) where {NC,E,Mov,Rad,Chi,D,U} =
    U

@inline _node_width(
    ::Type{<:StateLayout{NC,E,Mov,Rad,Chi,D,U}},
) where {NC,E,Mov,Rad,Chi,D,U} = (Mov ? (U === :velocities ? 2D : D) : 0) + Rad
@inline _cell_width(
    ::Type{<:StateLayout{NC,E,Mov,Rad,Chi,D,U,Car}},
) where {NC,E,Mov,Rad,Chi,D,U,Car} =
    NC + 1 + _carried_width(Car, D) + (Mov && U === :multipliers)
@inline _global_width(::Type{<:StateLayout{NC,E,Mov,Rad,Chi}}) where {NC,E,Mov,Rad,Chi} =
    2 * Chi

@inline node_width(l::StateLayout) = _node_width(typeof(l))
@inline cell_width(l::StateLayout) = _cell_width(typeof(l))

"""
	state_length(layout) -> Int

Length of the state vector.
"""
@inline state_length(l::StateLayout) = l.global_offset + _global_width(typeof(l))

@inline _node_offset(l::StateLayout, n::Int) = @inbounds l.node_offset[n]
@inline _cell_offset(l::StateLayout, i::Int) = @inbounds l.cell_offset[i]
@inline _global_offset(l::StateLayout) = l.global_offset

"""
	cell_rows(layout, i) -> UnitRange

Every row of the cell block of cell `i`.
"""
@inline cell_rows(l::StateLayout, i::Int) =
    (_cell_offset(l, i)+1):(_cell_offset(l, i)+cell_width(l))

"""
	node_rows(layout, n) -> UnitRange

Every row of the node block of node `n` (empty on a fixed mesh without
intensities).
"""
@inline node_rows(l::StateLayout, n::Int) =
    (_node_offset(l, n)+1):(_node_offset(l, n)+node_width(l))

"""
	mass_index(layout, i, j)

Row of the mass of component `j` in cell `i`.
"""
@inline mass_index(l::StateLayout, i::Int, j::Int) = _cell_offset(l, i) + j

"""
	energy_index(layout, i)

Row of the energy state of cell `i`.
"""
@inline energy_index(l::StateLayout, i::Int) = _cell_offset(l, i) + n_components(l) + 1

"""
	history_index(layout, i, k)

Row of the conversion-history inventory of cell `i`: `k = 1` the cumulative
dry gas produced `G`, `k = 2` the initial dry solid mass `M`.
"""
@inline function history_index(l::StateLayout, i::Int, k::Int)
    @boundscheck has_history(l) || throw(ArgumentError("layout carries no history pair"))
    return _cell_offset(l, i) +
           n_components(l) +
           1 +
           _carried_offset(carried_fields(l), dimension(l), :history) +
           k
end

"""
	deformation_index(layout, i, a, b)

Row of the entry `F[a, b]` of the deformation gradient of cell `i`.
"""
@inline function deformation_index(l::StateLayout, i::Int, a::Int, b::Int)
    @boundscheck has_deformation(l) ||
                 throw(ArgumentError("layout carries no deformation gradient"))
    D = dimension(l)
    return _cell_offset(l, i) +
           n_components(l) +
           1 +
           _carried_offset(carried_fields(l), D, :deformation) +
           D * (b - 1) +
           a
end

"""
	multiplier_index(layout, i)

Row of the closure's multiplier `φᵢ` of cell `i` (moving mesh, multipliers).
"""
@inline function multiplier_index(l::StateLayout, i::Int)
    @boundscheck (is_moving(l) && closure_unknowns(l) === :multipliers) ||
                 throw(ArgumentError("layout carries no multipliers"))
    return _cell_offset(l, i) +
           n_components(l) +
           1 +
           _carried_width(carried_fields(l), dimension(l)) +
           1
end

"""
	position_index(layout, n, k = 1)

Row of coordinate `k` of the position of node `n` (moving mesh).
"""
@inline function position_index(l::StateLayout, n::Int, k::Int = 1)
    @boundscheck is_moving(l) || throw(ArgumentError("layout carries no node positions"))
    return _node_offset(l, n) + k
end

"""
	velocity_index(layout, n, k = 1)

Row of component `k` of the velocity of node `n` (moving mesh, node
velocities).
"""
@inline function velocity_index(l::StateLayout, n::Int, k::Int = 1)
    @boundscheck (is_moving(l) && closure_unknowns(l) === :velocities) ||
                 throw(ArgumentError("layout carries no node velocities"))
    return _node_offset(l, n) + dimension(l) + k
end

"""
	intensity_index(layout, f)

Row of the radiant intensity of face `f` (Beer–Lambert). In 1D face `f` sits at
node `f`.
"""
@inline function intensity_index(l::StateLayout, f::Int)
    @boundscheck has_intensities(l) || throw(ArgumentError("layout carries no intensities"))
    return _node_offset(l, f) + node_width(l)
end

"""
	area_index(layout), lateral_rate_index(layout)

Rows of the cross-section measure `A` and its lateral strain rate `θ_A`
(variable cross section).
"""
@inline function area_index(l::StateLayout)
    @boundscheck has_chi(l) || throw(ArgumentError("layout carries no cross-section rows"))
    return _global_offset(l) + 1
end
@inline function lateral_rate_index(l::StateLayout)
    @boundscheck has_chi(l) || throw(ArgumentError("layout carries no cross-section rows"))
    return _global_offset(l) + 2
end

"""
	is_differential(layout) -> BitVector

`true` on the rows with a time derivative (masses, energy, node positions,
carried fields), `false` on the algebraic rows (node velocities,
multipliers, intensities, cross-section rows).
"""
is_differential(l::StateLayout) = l.differential

function _differential_rows(l::StateLayout)
    d = trues(state_length(l))
    if is_moving(l)
        if closure_unknowns(l) === :velocities
            for n = 1:n_nodes(l), k = 1:dimension(l)
                d[velocity_index(l, n, k)] = false
            end
        else
            for i = 1:n_cells(l)
                d[multiplier_index(l, i)] = false
            end
        end
    end
    if has_intensities(l)
        for f = 1:n_nodes(l)
            d[intensity_index(l, f)] = false
        end
    end
    if has_chi(l)
        d[area_index(l)] = false
        d[lateral_rate_index(l)] = false
    end
    return d
end

"""
	mass_matrix(layout) -> Diagonal

The integrator's mass matrix: one on the differential rows, zero on the
algebraic rows.
"""
mass_matrix(l::StateLayout) = Diagonal(Float64.(is_differential(l)))
