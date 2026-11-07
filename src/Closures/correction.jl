# ============================================================================
# The correction: the mesh velocity relative to the solid's
# ============================================================================

"""
	Correction

How the node velocities of a moving mesh depart from the solid's: the
per-node correction `c = w − vₛ` that [`correction!`](@ref) fills into the
derived state before the kernels run. The relative-flux kernel carries every
inventory across the faces at `−c`, and the balance rows hold the solid
velocity `w − c` to the cells' volume balance.
"""
abstract type Correction end

"""
	NoCorrection <: Correction

The mesh moves with the solid, `c ≡ 0`, and nothing crosses its faces.
"""
struct NoCorrection <: Correction end

"""
	PrescribedCorrection{F} <: Correction

A prescribed mesh velocity `w(x, t)` over a solid at rest, for verification:
every velocity row is `0 = wₙ − w(xₙ, t)`, and the correction is the state's
node velocity, `c = w`.
"""
struct PrescribedCorrection{F} <: Correction
    w::F
end

"""
	EQUIDISTRIBUTION_TIME

A fixed equidistribution time `τ` [s] for `EquidistributionCorrection(τ = …)`.
"""
const EQUIDISTRIBUTION_TIME = 720.0

"""
	ACTIVATION_TIME_FACTOR

The factor `λ` of [`ActivationTime`](@ref): a layer's equidistribution time
is `λ` times the time its thickness would take to vanish at its thinning rate
when the time is estimated. A layer synchronizes in about a fifth of that
estimated life.
"""
const ACTIVATION_TIME_FACTOR = 1.8

"""
	ActivationTime(; λ = ACTIVATION_TIME_FACTOR)

A layer's equidistribution time from its own thinning, estimated twice: when
its first cell activates (falls below `r_hi` of its initial thickness) and
when its first cell falls below `r_lo`, each time from the layer's thickness
`h` and its thinning rate `ḣ` then,

	τ = λ h / |ḣ|.

The first estimate holds the activated cells together while the layer's rate
still rises; the second, taken nearer the layer's end, is blended in over the
correction's second switch ([`EquidistributionCorrection`](@ref)), so neither
estimate changes the right-hand side where it is taken. The first estimate
needs a cell to fall below `r_hi`, so the activation's `r_hi` must be below 1.
A layer not thinning at its first estimate has no correction until its
second; one not thinning at its second keeps its first. The rate scales with the physical rates and `λ` is
a pure number, so a deck whose rates are all multiplied by `s` gets `τ/s`; no
estimate reads the solve's horizon. Before the first estimate the layer has no
time and no gain, and a cell that re-activates after the second changes
nothing.
"""
struct ActivationTime
    λ::Float64
    function ActivationTime(; λ::Real = ACTIVATION_TIME_FACTOR)
        λ > 0 ||
            throw(ArgumentError("the activation time's factor must be positive, got $λ"))
        return new(Float64(λ))
    end
end

"""
	EquidistributionCorrection(; τ = ActivationTime(), activation = (0.3, 0.1),
	                           layers = :vanishing) <: Correction

The relaxed equidistribution of a vanishing layer's cells. At an interior
node `n` of a layer of `N` cells, between cell `n − 1` (inward) and cell `n`,

	cₙ = aₙ N² (hₙ − hₙ₋₁) / τ,

which moves the node into the thicker cell, so `h` relaxes
toward a common value over the layer in about `τ/π²` whatever `N`: `τ` [s] is
the layer's equidistribution time, by default estimated for each layer from
its thinning ([`ActivationTime`](@ref)), or one value for every layer. The
activation `aₙ = 1 − (1 − aₙ₋₁)(1 − aₙ)` is the smooth "or" of the two cells'
switches `aᵢ = S(xᵢ)`, `xᵢ = (ln r_hi − ln(hᵢ/hᵢ(0))) / (ln r_hi − ln r_lo)`,
`S(x) = x²(3 − 2x)` clamped to `[0, 1]`, with `activation = (r_hi, r_lo)`:
zero while a cell keeps at least `r_hi` of its initial thickness, one at and
below `r_lo`, `C¹` in `ln h` between. Under [`ActivationTime`](@ref) `r_hi`
must be below 1, while a fixed `τ` accepts any `r_hi`, and the layer's first
estimate `τ₁` gives the time above and its second `τ₂` enters
through the cells' second switches `bᵢ = S(xᵢ − 1)`, zero at and above `r_lo`
and one at and below `r_lo²/r_hi`:

	cₙ = (aₙ/τ₁ + bₙ (1/τ₂ − 1/τ₁)) N² (hₙ − hₙ₋₁),

`bₙ` the smooth "or" of the two cells' `bᵢ`.

A layer is bounded by the pinned node, the mesh's Lagrangian nodes (the layer
interfaces) and the end nodes; `c = 0` at each of them, at every node beside
a dead cell (the outermost live node among them) and wherever both cells keep
at least `r_hi` of their initial thickness. With `layers = :vanishing` only
the layers deactivation can remove move ([`vanishes`](@ref)) and the others
keep their mesh with the solid; `:all` moves every layer, for decks that
exercise the correction on layers that cannot vanish. `setup` binds the correction to the
mesh ([`bind_correction`](@ref)).
"""
struct EquidistributionCorrection{T<:Union{Float64,ActivationTime}} <: Correction
    τ::T
    r_hi::Float64
    r_lo::Float64
    layers::Symbol
end

function EquidistributionCorrection(;
    τ::Union{Real,ActivationTime} = ActivationTime(),
    activation::Tuple{Real,Real} = (0.3, 0.1),
    layers::Symbol = :vanishing,
)
    r_hi, r_lo = activation
    τ isa ActivationTime ||
        τ > 0 ||
        throw(ArgumentError("the equidistribution time must be positive, got $τ"))
    0 < r_lo < r_hi || throw(
        ArgumentError(
            "the activation needs 0 < r_lo < r_hi, got (r_hi, r_lo) = $activation",
        ),
    )
    τ isa ActivationTime &&
        r_hi >= 1 &&
        throw(
            ArgumentError(
                "an ActivationTime takes a layer's first estimate when one of its cells " *
                "first falls below r_hi of its initial thickness, so it needs r_hi < 1, " *
                "got (r_hi, r_lo) = $activation; a fixed τ accepts any activation",
            ),
        )
    layers in (:vanishing, :all) ||
        throw(ArgumentError("layers must be :vanishing or :all, got $(repr(layers))"))
    return EquidistributionCorrection(
        τ isa ActivationTime ? τ : Float64(τ),
        r_hi,
        r_lo,
        layers,
    )
end

"""
	LayerCorrection <: Correction

An [`EquidistributionCorrection`](@ref) bound to a mesh: per node the gain
`N²/τ₁` of its layer and the shift `N²/τ₂ − N²/τ₁` of its second estimate
(both zero at the layers' bounding nodes, the shift zero under one fixed time),
by slot, and the same by initial node; the activation's bounds `r_hi` and `r_lo`
with `ln r_hi` and `ln r_hi − ln r_lo`; per moving layer its bounding nodes on
the initial mesh, its equidistribution time (`Inf` before the first estimate,
then the latest), the time of its latest estimate (`NaN` before the first), the
stage of its estimate (zero before the first activation, one after the first
estimate, two after the second) and its first estimate, with the
[`ActivationTime`](@ref) rule (`nothing` for one fixed time); per initial cell
the layer it starts in, and per slot the initial cell there.
"""
struct LayerCorrection <: Correction
    gain::Vector{Float64}
    initial_gain::Vector{Float64}
    shift::Vector{Float64}
    initial_shift::Vector{Float64}
    r_hi::Float64
    r_lo::Float64
    log_hi::Float64
    log_span::Float64
    layers::Vector{Tuple{Int,Int}}
    time::Vector{Float64}
    set_at::Vector{Float64}
    stage::Vector{Int}
    first_estimate::Vector{Float64}
    estimate::Union{Nothing,ActivationTime}
    layer::Vector{Int}
    origin::Vector{Int}
end

"""
	residues(material) -> NTuple{NC,Float64}

The condensed mass each condensed component of `material` leaves, as a
fraction of its own, once the reactions have converted it completely: the
condensed yields of a reaction that consumes it, each weighted by that
product's own residue, the largest over the reactions that consume it (they
compete, so the least gasifying one bounds what is left). Zero for a
component every consuming reaction gasifies; one for an inert component, a
residue no reaction consumes and the components of a cycle that yields
condensed mass only. A reaction that yields part of its own reactant back is
read by the rest of its products, renormalized; one that yields the whole of
it back consumes nothing. A gas component reads zero.
"""
residues(material::Material{NC}) where {NC} = _left_by_conversion(
    material,
    ntuple(j -> is_gas(material.components[j]) ? 0.0 : 1.0, Val(NC)),
)

"""
	retained_volumes(material, T) -> NTuple{NC,Float64}

The volume [m³/kg] each condensed component of `material` leaves per unit of its
mass once the reactions have converted it completely, at the densities at `T`:
the [`residues`](@ref) rule with each component that no reaction consumes
counted by its volume. A low-density char keeps far more of a cell than its mass
suggests. A gas component reads zero.
"""
retained_volumes(material::Material{NC}, T::Real) where {NC} = _left_by_conversion(
    material,
    ntuple(
        j ->
            is_gas(material.components[j]) ? 0.0 : 1 / material.components[j].density(T),
        Val(NC),
    ),
)

# What each condensed component leaves once converted, a component no reaction
# consumes leaving `left_alone[j]` per unit of its mass.
function _left_by_conversion(material::Material{NC}, left_alone::NTuple{NC}) where {NC}
    condensed = ntuple(j -> !is_gas(material.components[j]), Val(NC))
    returned(rxn, j) = sum((last(p) for p in rxn.products if first(p) == j); init = 0.0)
    consumes(rxn, j) = rxn.reactant == j && returned(rxn, j) < 1
    left(rxn, j, r) =
        sum((last(p) * r[first(p)] for p in rxn.products if first(p) != j); init = 0.0) /
        (1 - returned(rxn, j))
    consumers = [[rxn for rxn in material.reactions if consumes(rxn, j)] for j = 1:NC]
    largest = maximum(left_alone)
    r = [
        condensed[j] ? (isempty(consumers[j]) ? left_alone[j] : largest) : 0.0 for j = 1:NC
    ]
    # From the largest value, they fall to the greatest fixed point: exactly along
    # a chain, geometrically around a cycle that loses mass to the gas.
    for _ = 1:RESIDUE_PASSES
        change = 0.0
        for j = 1:NC
            (condensed[j] && !isempty(consumers[j])) || continue
            x = min(r[j], maximum(rxn -> left(rxn, j, r), consumers[j]))
            change = max(change, r[j] - x)
            r[j] = x
        end
        change <= eps() && break
    end
    return ntuple(j -> r[j], Val(NC))
end

# The passes `residues` takes at most.
const RESIDUE_PASSES = 10_000

"""
	vanishes(residues, m0, cells, condensed, threshold) -> Bool

Whether deactivation can remove the cells `cells`, with initial masses `m0`
(`NC × n`): some of them hold condensed mass (`condensed` flags the condensed
components), and each that does keeps less than the fraction `threshold` of
it once its components are converted completely ([`residues`](@ref)).
"""
function vanishes(
    r::NTuple{NC,Float64},
    m0::AbstractMatrix,
    cells,
    condensed,
    threshold,
) where {NC}
    held = false
    for i in cells
        total = sum((m0[j, i] for j = 1:NC if condensed[j]); init = 0.0)
        total > 0 || continue
        held = true
        kept = sum((r[j] * m0[j, i] for j = 1:NC if condensed[j]); init = 0.0)
        kept < threshold * total || return false
    end
    return held
end

"""
	bind_correction(correction, topology, reference, m0, material, threshold) -> Correction

The correction a solve evaluates: an [`EquidistributionCorrection`](@ref)
bound, as a [`LayerCorrection`](@ref), to the layers of the line mesh
`topology` and to its initial masses `m0` (`NC × n`), with `layers =
:vanishing` a layer that deactivation at `threshold` cannot remove
([`vanishes`](@ref)) with zero gain; any other correction as given.
"""
bind_correction(corr::Correction, topo, reference, m0, material, threshold) = corr

function bind_correction(
    corr::EquidistributionCorrection,
    topo::MeshTopology{1},
    reference::CellReference,
    m0::AbstractMatrix{Float64},
    material::Material{NC},
    threshold::Real,
) where {NC}
    n = n_cells(topo)
    nn = n_nodes(topo)
    all(i -> topo.cells[i].nodes == (i, i + 1), 1:n) || throw(
        ArgumentError(
            "the equidistribution correction needs a line mesh whose cell i lies " *
            "between nodes i and i + 1",
        ),
    )
    condensed = ntuple(j -> !is_gas(material.components[j]), Val(NC))
    bounds = sort!(unique([1; nn; topo.pinned; topo.lagrangian]))
    v = residues(material)
    layers = Tuple{Int,Int}[]
    layer = zeros(Int, n)
    for (k, (lo, hi)) in enumerate(zip(bounds[1:(end-1)], bounds[2:end]))
        moves = corr.layers === :all || vanishes(v, m0, lo:(hi-1), condensed, threshold)
        moves && push!(layers, (lo, hi))
        layer[lo:(hi-1)] .= k
    end
    estimate = corr.τ isa ActivationTime ? corr.τ : nothing
    nl = length(layers)
    lc = LayerCorrection(
        zeros(nn),
        zeros(nn),
        zeros(nn),
        zeros(nn),
        corr.r_hi,
        corr.r_lo,
        log(corr.r_hi),
        log(corr.r_hi) - log(corr.r_lo),
        layers,
        fill(Inf, nl),
        fill(NaN, nl),
        zeros(Int, nl),
        fill(Inf, nl),
        estimate,
        layer,
        collect(1:n),
    )
    return _bind_times!(lc, corr.τ)
end

# The bound correction's layer times: one fixed time from the start, or every
# layer before its first activation.
function _bind_times!(corr::LayerCorrection, τ::Float64)
    for k in eachindex(corr.layers)
        _layer_gains!(corr, k, τ, τ)
        corr.time[k] = τ
        corr.set_at[k] = 0.0
        corr.stage[k] = 2
        corr.first_estimate[k] = τ
    end
    return corr
end
_bind_times!(corr::LayerCorrection, ::ActivationTime) = reset_layer_times!(corr)

# Layer `k`'s gain `N²/τ₁` and shift `N²/τ₂ − N²/τ₁` at its nodes, by initial
# node and by slot (an infinite time gives no gain).
function _layer_gains!(corr::LayerCorrection, k::Int, τ₁::Float64, τ₂::Float64)
    lo, hi = corr.layers[k]
    g₁, g₂ = (hi - lo)^2 / τ₁, (hi - lo)^2 / τ₂
    corr.initial_gain[(lo+1):(hi-1)] .= g₁
    corr.initial_shift[(lo+1):(hi-1)] .= g₂ - g₁
    slot_gains!(corr)
    return corr
end

"""
	slot_gains!(correction::LayerCorrection) -> correction

The gains and shifts by slot: a node between two cells of one layer keeps
those of the outer cell's inward node, and a node between cells of two layers
has none, so a correction never moves a node across the layers it bounds,
whichever cells have left the column between them.
"""
function slot_gains!(corr::LayerCorrection)
    origin = corr.origin
    for v in (corr.gain, corr.shift)
        v[1] = 0.0
        v[end] = 0.0
    end
    for s = 2:length(origin)
        below, above = origin[s-1], origin[s]
        same = corr.layer[below] == corr.layer[above]
        corr.gain[s] = same ? corr.initial_gain[above] : 0.0
        corr.shift[s] = same ? corr.initial_shift[above] : 0.0
    end
    return corr
end

"""
	estimates_layer_times(correction) -> Bool

Whether `correction` sets some layer's equidistribution time during a solve
([`ActivationTime`](@ref)).
"""
estimates_layer_times(::Correction) = false
estimates_layer_times(corr::LayerCorrection) =
    corr.estimate !== nothing && !isempty(corr.layers)

"""
	reset_layer_times!(correction::LayerCorrection) -> correction

Every layer of an [`ActivationTime`](@ref) correction back before its first
activation, with no time and no gain; a correction with one fixed time is
unchanged.
"""
function reset_layer_times!(corr::LayerCorrection)
    corr.estimate === nothing && return corr
    for k in eachindex(corr.layers)
        _layer_gains!(corr, k, Inf, Inf)
        corr.time[k] = Inf
        corr.set_at[k] = NaN
        corr.stage[k] = 0
        corr.first_estimate[k] = Inf
    end
    return corr
end

# The live slots of layer `k`, contiguous since the live cells keep their order.
function _live_cells(corr::LayerCorrection, layout::StateLayout, k::Int)
    lo, hi = corr.layers[k]
    first_slot, last_slot = 0, -1
    for s in eachindex(corr.origin)
        layout.mask[energy_index(layout, s)] || break
        lo <= corr.origin[s] < hi || continue
        first_slot == 0 && (first_slot = s)
        last_slot = s
    end
    return first_slot == 0 ? (lo:(lo-1)) : (first_slot:last_slot)
end

"""
	layer_time_margin(correction, layout, reference, u) -> Float64

The smallest margin `hᵢ/(r hᵢ(0)) − 1` over the live cells of the layers
still to be estimated, `r` the activation's `r_hi` before a layer's first
estimate and `r_lo` before its second; `Inf` when none is waiting. It crosses
zero downward where [`estimate_layer_times!`](@ref) runs.
"""
function layer_time_margin(
    corr::LayerCorrection,
    layout::StateLayout,
    reference::CellReference,
    u::AbstractVector,
)
    margin = Inf
    for k in eachindex(corr.layers)
        stage = corr.stage[k]
        stage == 2 && continue
        r = stage == 0 ? corr.r_hi : corr.r_lo
        for i in _live_cells(corr, layout, k)
            h = u[position_index(layout, i + 1)] - u[position_index(layout, i)]
            margin = min(margin, h / (r * reference.thickness[i]) - 1)
        end
    end
    return margin
end

"""
	estimate_layer_times!(correction, layout, reference, u, t; slack) -> correction

At the state `u` at time `t`, every waiting layer whose margin
([`layer_time_margin`](@ref)) is at most `slack` takes its next estimate of
its equidistribution time from its thickness `h` (the distance between its
bounding nodes, the outer one its outermost live node) and its thinning rate
`ḣ` (their velocity difference) as [`ActivationTime`](@ref) states: the first
sets its gain, the second its shift.
"""
function estimate_layer_times!(
    corr::LayerCorrection,
    layout::StateLayout,
    reference::CellReference,
    u::AbstractVector,
    t::Real;
    slack::Real = 0.0,
)
    est = corr.estimate
    est === nothing && return corr
    for k in eachindex(corr.layers)
        stage = corr.stage[k]
        stage == 2 && continue
        r = stage == 0 ? corr.r_hi : corr.r_lo
        cells = _live_cells(corr, layout, k)
        isempty(cells) && continue
        crossed = any(cells) do i
            h = u[position_index(layout, i + 1)] - u[position_index(layout, i)]
            h / (r * reference.thickness[i]) - 1 <= slack
        end
        crossed || continue
        lo, hi = first(cells), last(cells) + 1
        h = u[position_index(layout, hi)] - u[position_index(layout, lo)]
        ḣ = u[velocity_index(layout, hi)] - u[velocity_index(layout, lo)]
        τ = ḣ < 0 ? est.λ * h / -ḣ : stage == 0 ? Inf : corr.first_estimate[k]
        if stage == 0
            corr.first_estimate[k] = τ
            _layer_gains!(corr, k, τ, τ)
        else
            _layer_gains!(corr, k, corr.first_estimate[k], τ)
        end
        corr.time[k] = τ
        corr.set_at[k] = Float64(t)
        corr.stage[k] = stage + 1
    end
    return corr
end

"""
	activation(correction::LayerCorrection, h, h0)

The switch of a cell of thickness `h` and initial thickness `h0`: zero at and
above `r_hi h0`, one at and below `r_lo h0` (zero and negative thicknesses
included, the thickness floored smoothly at `DISTANCE_FLOOR_FRACTION · h0`
before the logarithm), and `C¹` in `ln h` between.
"""
@inline activation(corr::LayerCorrection, h, h0) = first(switches(corr, h, h0))

"""
	switches(correction::LayerCorrection, h, h0) -> (a, b)

The two switches of a cell of thickness `h` and initial thickness `h0`: its
[`activation`](@ref) `a = S(x)` and its second switch `b = S(x − 1)`, zero at
and above `r_lo h0` and one at and below `(r_lo²/r_hi) h0`, `C¹` in `ln h`
between, `x` as in [`EquidistributionCorrection`](@ref).
"""
@inline function switches(corr::LayerCorrection, h, h0)
    h >= corr.r_hi * h0 && return (zero(h), zero(h))
    r = smooth_floor(h, DISTANCE_FLOOR_FRACTION * h0) / h0
    x = (corr.log_hi - log(r)) / corr.log_span
    x >= 2 && return (one(x), one(x))
    x >= 1 && return (one(x), _smoothstep(x - 1))
    return (_smoothstep(x), zero(x))
end

@inline _smoothstep(x) = x * x * (3 - 2x)

# The node corrections `c` at the derived state `ds` (`correction!`).
@inline _correction!(c, ::NoCorrection, ds, reference) = nothing

@inline function _correction!(c, ::PrescribedCorrection, ds, reference)
    copyto!(c, ds.w)
    return nothing
end

function _correction!(c, corr::LayerCorrection, ds::DerivedState, reference)
    x = ds.x
    h0 = reference.thickness
    nn = length(c)
    T = eltype(c)
    c[1] = zero(T)
    c[nn] = zero(T)
    h_in = x[2] - x[1]
    a_in, b_in = switches(corr, h_in, h0[1])
    @inbounds for nd = 2:(nn-1)
        h_out = x[nd+1] - x[nd]
        a_out, b_out = switches(corr, h_out, h0[nd])
        a = 1 - (1 - a_in) * (1 - a_out)
        g, s = corr.gain[nd], corr.shift[nd]
        c[nd] =
            if (iszero(g) && iszero(s)) || !(ds.active[nd-1] && ds.active[nd]) || iszero(a)
                zero(T)
            elseif iszero(s)
                a * g * (h_out - h_in)
            else
                (a * g + (1 - (1 - b_in) * (1 - b_out)) * s) * (h_out - h_in)
            end
        h_in, a_in, b_in = h_out, a_out, b_out
    end
    return nothing
end
