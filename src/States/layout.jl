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
	StateLayout{NC,E,Mov,Rad,Chi}

Position of every unknown in the state vector of a 1D mesh with `NC`
components and energy state `E`.

The ordering is interleaved by position — node 1, cell 1, node 2, …, cell N,
node N + 1, then the global scalars — so the Jacobian is banded in every
geometry. A node block holds, in order, the node position `x` and velocity `w`
(moving mesh, `Mov`) and the intensity `I` of the face at that node
(Beer–Lambert, `Rad`). A cell block holds the component masses `m₁ … m_NC`,
the energy state, and the conversion-history pair `G`, `M` (variable cross
section, `Chi`). The global block holds the cross-section measure `A` and its
lateral strain rate `θ_A` (`Chi`).

Kernels never compute an index; every access goes through the index
functions ([`mass_index`](@ref), [`energy_index`](@ref), …). `mask` marks the
live rows, a dead row being an identity row of the system; `differential`
marks the rows with a time derivative (see [`is_differential`](@ref)).
"""
struct StateLayout{NC,E<:EnergyState,Mov,Rad,Chi}
    n_cells::Int
    mask::BitVector
    differential::BitVector
end

"""
	StateLayout(NC, n_cells, energy; moving = false, intensities = false, chi = false)

The layout for `n_cells` cells of a 1D mesh with every row live.
"""
function StateLayout(
    NC::Int,
    n_cells::Int,
    energy::EnergyState;
    moving::Bool = false,
    intensities::Bool = false,
    chi::Bool = false,
)
    NC >= 1 || throw(ArgumentError("a layout needs at least one component"))
    n_cells >= 1 || throw(ArgumentError("a layout needs at least one cell"))
    chi && !moving && throw(ArgumentError("the variable cross section needs a moving mesh"))
    L = StateLayout{NC,typeof(energy),moving,intensities,chi}
    len = _state_length(L, n_cells)
    layout = L(n_cells, trues(len), trues(len))
    layout.differential .= _differential_rows(layout)
    return layout
end

@inline n_components(::StateLayout{NC}) where {NC} = NC
@inline energy_state(::StateLayout{NC,E}) where {NC,E} = E()
@inline is_moving(::StateLayout{NC,E,Mov}) where {NC,E,Mov} = Mov
@inline has_intensities(::StateLayout{NC,E,Mov,Rad}) where {NC,E,Mov,Rad} = Rad
@inline has_chi(::StateLayout{NC,E,Mov,Rad,Chi}) where {NC,E,Mov,Rad,Chi} = Chi
@inline n_cells(l::StateLayout) = l.n_cells
@inline n_nodes(l::StateLayout) = l.n_cells + 1

@inline _node_width(::Type{<:StateLayout{NC,E,Mov,Rad}}) where {NC,E,Mov,Rad} =
    2 * Mov + Rad
@inline _cell_width(::Type{<:StateLayout{NC,E,Mov,Rad,Chi}}) where {NC,E,Mov,Rad,Chi} =
    NC + 1 + 2 * Chi
@inline _global_width(::Type{<:StateLayout{NC,E,Mov,Rad,Chi}}) where {NC,E,Mov,Rad,Chi} =
    2 * Chi
@inline _state_length(L::Type{<:StateLayout}, n::Int) =
    n * (_node_width(L) + _cell_width(L)) + _node_width(L) + _global_width(L)

@inline node_width(l::StateLayout) = _node_width(typeof(l))
@inline cell_width(l::StateLayout) = _cell_width(typeof(l))

"""
	state_length(layout) -> Int

Length of the state vector.
"""
@inline state_length(l::StateLayout) = _state_length(typeof(l), l.n_cells)

@inline _stride(l::StateLayout) = node_width(l) + cell_width(l)
@inline _node_offset(l::StateLayout, n::Int) = (n - 1) * _stride(l)
@inline _cell_offset(l::StateLayout, i::Int) = (i - 1) * _stride(l) + node_width(l)
@inline _global_offset(l::StateLayout) = l.n_cells * _stride(l) + node_width(l)

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
    @boundscheck has_chi(l) || throw(ArgumentError("layout carries no history pair"))
    return _cell_offset(l, i) + n_components(l) + 1 + k
end

"""
	position_index(layout, n)

Row of the position of node `n` (moving mesh).
"""
@inline function position_index(l::StateLayout, n::Int)
    @boundscheck is_moving(l) || throw(ArgumentError("layout carries no node positions"))
    return _node_offset(l, n) + 1
end

"""
	velocity_index(layout, n)

Row of the velocity of node `n` (moving mesh).
"""
@inline function velocity_index(l::StateLayout, n::Int)
    @boundscheck is_moving(l) || throw(ArgumentError("layout carries no node velocities"))
    return _node_offset(l, n) + 2
end

"""
	intensity_index(layout, f)

Row of the radiant intensity of face `f` (Beer–Lambert). In 1D face `f` sits at
node `f`.
"""
@inline function intensity_index(l::StateLayout, f::Int)
    @boundscheck has_intensities(l) || throw(ArgumentError("layout carries no intensities"))
    return _node_offset(l, f) + 2 * is_moving(l) + 1
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
history pair), `false` on the algebraic rows (node velocities, intensities,
cross-section rows).
"""
is_differential(l::StateLayout) = l.differential

function _differential_rows(l::StateLayout)
    d = trues(state_length(l))
    if is_moving(l)
        for n = 1:n_nodes(l)
            d[velocity_index(l, n)] = false
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
