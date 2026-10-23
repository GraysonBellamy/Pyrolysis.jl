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

The fallback equidistribution time `τ` [s] of [`ActivationTime`](@ref), for
a layer whose thinning gives no estimate: with it the 5.9 mm PMMA layer of
V-BO-01 at 50 kW/m² ends in one deactivation at 30, 60 and 120 cells across
it.
"""
const EQUIDISTRIBUTION_TIME = 720.0

"""
	ACTIVATION_TIME_FACTOR

The factor `λ` of [`ActivationTime`](@ref): a layer's equidistribution time
is `λ` times the time its thickness would take to vanish at the thinning rate
of its first activation. With it the 5.9 mm PMMA layer of V-BO-01 at 50 kW/m²
gets about [`EQUIDISTRIBUTION_TIME`](@ref), and a layer synchronizes in about
a fifth of that estimated life.
"""
const ACTIVATION_TIME_FACTOR = 1.8

"""
	ActivationTime(; λ = ACTIVATION_TIME_FACTOR, fallback = EQUIDISTRIBUTION_TIME,
	               floor = 1.0)

A layer's equidistribution time set once, when its first cell activates
(falls below `r_hi` of its initial thickness), from the layer's thickness
`h` and its thinning rate `ḣ` at that time:

	τ = λ h / |ḣ|.

The rate scales with the physical rates and `λ` is a pure number, so a deck
whose rates are all multiplied by `s` gets `τ/s`. Before the estimate the
layer's correction is zero, since no cell of it is active. A layer thinning
no faster than `floor` times its initial thickness over the solve's horizon
keeps `fallback` [s] until a cell crosses `r_lo` of its initial thickness,
where the estimate is taken once more; a layer that thickens keeps
`fallback`. A layer is estimated once: a cell that re-activates later changes
nothing.
"""
struct ActivationTime
    λ::Float64
    fallback::Float64
    floor::Float64
    function ActivationTime(;
        λ::Real = ACTIVATION_TIME_FACTOR,
        fallback::Real = EQUIDISTRIBUTION_TIME,
        floor::Real = 1.0,
    )
        λ > 0 ||
            throw(ArgumentError("the activation time's factor must be positive, got $λ"))
        fallback > 0 ||
            throw(ArgumentError("the fallback time must be positive, got $fallback"))
        floor >= 0 ||
            throw(ArgumentError("the rate floor must not be negative, got $floor"))
        return new(Float64(λ), Float64(fallback), Float64(floor))
    end
end

"""
	EquidistributionCorrection(; τ = ActivationTime(), activation = (0.3, 0.1),
	                           monitor = (1.0, 0.0), layers = :vanishing) <: Correction

The relaxed equidistribution of a vanishing layer's cells. At an interior
node `n` of a layer of `N` cells, between cell `n − 1` (inward) and cell `n`,

	cₙ = aₙ N² (M̃ₙ hₙ − M̃ₙ₋₁ hₙ₋₁) / τ,

which moves the node into the thicker (or heavier) cell, so `M̃ h` relaxes
toward a common value over the layer in about `τ/π²` whatever `N`: `τ` [s] is
the layer's equidistribution time, by default set for each layer at its first
activation from its thinning ([`ActivationTime`](@ref)), or one value for
every layer. The activation `aₙ = 1 − (1 − aₙ₋₁)(1 − aₙ)`
is the smooth "or" of the two cells' switches
`aᵢ = S((ln r_hi − ln(hᵢ/hᵢ(0))) / (ln r_hi − ln r_lo))`, `S(x) = x²(3 − 2x)`
clamped to `[0, 1]`, with `activation = (r_hi, r_lo)`: zero while a cell keeps
at least `r_hi` of its initial thickness, one at and below `r_lo`, `C¹` in
`ln h` between. The monitor `M = α + β ρ/ρ_ref`, with `monitor = (α, β)`, `ρ`
a cell's condensed mass per volume and `ρ_ref` the largest initial one over
the mesh, is divided by its mean over the layer's cells at the initial state,
so `α` and `β` set its shape and not the relaxation time; with `β = 0`,
`M̃ = 1`.

A layer is bounded by the pinned node, the mesh's Lagrangian nodes (the layer
interfaces) and the end nodes; `c = 0` at each of them, at every node beside
a dead cell (the outermost live node among them) and wherever both cells keep
at least `r_hi` of their initial thickness. With `layers = :vanishing` only
the layers that can vanish move ([`vanishes`](@ref)) and the others keep their
mesh with the solid; `:all` moves every layer, for decks that exercise the
correction on layers that cannot vanish. `setup` binds the correction to the
mesh ([`bind_correction`](@ref)).
"""
struct EquidistributionCorrection{T<:Union{Float64,ActivationTime}} <: Correction
    τ::T
    r_hi::Float64
    r_lo::Float64
    α::Float64
    β::Float64
    layers::Symbol
end

function EquidistributionCorrection(;
    τ::Union{Real,ActivationTime} = ActivationTime(),
    activation::Tuple{Real,Real} = (0.3, 0.1),
    monitor::Tuple{Real,Real} = (1.0, 0.0),
    layers::Symbol = :vanishing,
)
    r_hi, r_lo = activation
    α, β = monitor
    τ isa ActivationTime ||
        τ > 0 ||
        throw(ArgumentError("the equidistribution time must be positive, got $τ"))
    0 < r_lo < r_hi || throw(
        ArgumentError(
            "the activation needs 0 < r_lo < r_hi, got (r_hi, r_lo) = $activation",
        ),
    )
    (α > 0 && β >= 0) ||
        throw(ArgumentError("the monitor needs α > 0 and β ≥ 0, got (α, β) = $monitor"))
    layers in (:vanishing, :all) ||
        throw(ArgumentError("layers must be :vanishing or :all, got $(repr(layers))"))
    return EquidistributionCorrection(
        τ isa ActivationTime ? τ : Float64(τ),
        r_hi,
        r_lo,
        α,
        β,
        layers,
    )
end

"""
	LayerCorrection{NC} <: Correction

An [`EquidistributionCorrection`](@ref) bound to a mesh: per node the gain
`N²/τ` of its layer (zero at the layers' bounding nodes), per cell the
reciprocal of its layer's mean initial monitor, the activation's bounds
`r_hi` and `r_lo` with `ln r_hi` and `ln r_hi − ln r_lo`, the monitor's `α`
and `β/ρ_ref`, and the condensed component flags; per moving layer its
bounding nodes, its equidistribution time, the time that time was set at
(`NaN` while it is not) and the stage of its estimate (zero before the first
activation, one while a floored rate waits for the second, two once set),
with the [`ActivationTime`](@ref) rule (`nothing` for one fixed time).
"""
struct LayerCorrection{NC} <: Correction
    gain::Vector{Float64}
    monitor_scale::Vector{Float64}
    r_hi::Float64
    r_lo::Float64
    log_hi::Float64
    log_span::Float64
    α::Float64
    β::Float64
    condensed::NTuple{NC,Bool}
    layers::Vector{Tuple{Int,Int}}
    time::Vector{Float64}
    set_at::Vector{Float64}
    stage::Vector{Int}
    estimate::Union{Nothing,ActivationTime}
end

"""
	vanishing(material) -> NTuple{NC,Bool}

Whether each condensed component of `material` gasifies completely: some
reaction consumes it, and every reaction that does yields no condensed
component that does not itself gasify completely. The least such set, so an
inert component, a residue that no reaction consumes and the components of a
cycle do not vanish; a gas component reads `false`. A reaction that yields
part of its own reactant back is read by the rest of its products, and one
that yields the whole of it back consumes nothing.
"""
function vanishing(material::Material{NC}) where {NC}
    condensed = ntuple(j -> !is_gas(material.components[j]), Val(NC))
    returned(rxn, j) = sum((last(p) for p in rxn.products if first(p) == j); init = 0.0)
    consumes(rxn, j) = any(r -> first(r) == j, rxn.reactants) && returned(rxn, j) < 1
    v = fill(false, NC)
    changed = true
    while changed
        changed = false
        for j = 1:NC
            (condensed[j] && !v[j]) || continue
            consumers = [rxn for rxn in material.reactions if consumes(rxn, j)]
            isempty(consumers) && continue
            gasified =
                p -> first(p) == j || !condensed[first(p)] || iszero(last(p)) || v[first(p)]
            all(rxn -> all(gasified, rxn.products), consumers) || continue
            v[j] = true
            changed = true
        end
    end
    return ntuple(j -> v[j], Val(NC))
end

"""
	vanishes(vanishing, m0, cells, condensed) -> Bool

Whether the cells `cells`, with initial masses `m0` (`NC × n`), can vanish:
they hold condensed mass (`condensed` flags the condensed components), and
every condensed component present in any of them vanishes (`vanishing`,
[`vanishing`](@ref)'s flags).
"""
function vanishes(v::NTuple{NC,Bool}, m0::AbstractMatrix, cells, condensed) where {NC}
    present(j) = condensed[j] && any(i -> m0[j, i] > 0, cells)
    return any(present, 1:NC) && all(j -> !present(j) || v[j], 1:NC)
end

"""
	bind_correction(correction, topology, reference, m0, material) -> Correction

The correction a solve evaluates: an [`EquidistributionCorrection`](@ref)
bound, as a [`LayerCorrection`](@ref), to the layers of the line mesh
`topology` and to its initial masses `m0` (`NC × n`), with `layers =
:vanishing` a layer that cannot vanish ([`vanishes`](@ref)) with zero gain;
any other correction as given.
"""
bind_correction(corr::Correction, topo, reference, m0, material) = corr

function bind_correction(
    corr::EquidistributionCorrection,
    topo::MeshTopology{1},
    reference::CellReference,
    m0::AbstractMatrix{Float64},
    material::Material{NC},
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
    ρ0 = [
        sum(m0[j, i] for j = 1:NC if condensed[j]; init = 0.0) / reference.volume[i] for
        i = 1:n
    ]
    ρ_ref = maximum(ρ0)
    β = ρ_ref > 0 ? corr.β / ρ_ref : 0.0
    bounds = sort!(unique([1; nn; topo.pinned; topo.lagrangian]))
    v = vanishing(material)
    gain = zeros(nn)
    scale = zeros(n)
    layers = Tuple{Int,Int}[]
    for (lo, hi) in zip(bounds[1:(end-1)], bounds[2:end])
        N = hi - lo
        moves = corr.layers === :all || vanishes(v, m0, lo:(hi-1), condensed)
        moves && push!(layers, (lo, hi))
        M̄ = sum(corr.α + β * ρ0[i] for i = lo:(hi-1)) / N
        scale[lo:(hi-1)] .= 1 / M̄
    end
    estimate = corr.τ isa ActivationTime ? corr.τ : nothing
    lc = LayerCorrection(
        gain,
        scale,
        corr.r_hi,
        corr.r_lo,
        log(corr.r_hi),
        log(corr.r_hi) - log(corr.r_lo),
        corr.α,
        β,
        condensed,
        layers,
        zeros(length(layers)),
        zeros(length(layers)),
        zeros(Int, length(layers)),
        estimate,
    )
    return _bind_times!(lc, corr.τ)
end

# The bound correction's layer times: one fixed time set from the start, or
# every layer before its first activation.
function _bind_times!(corr::LayerCorrection, τ::Float64)
    for k in eachindex(corr.layers)
        _set_layer_time!(corr, k, τ, 0.0)
    end
    return corr
end
_bind_times!(corr::LayerCorrection, ::ActivationTime) = reset_layer_times!(corr)

# Layer `k`'s equidistribution time `τ`, set at `t`, and its nodes' gain.
function _set_layer_time!(corr::LayerCorrection, k::Int, τ::Float64, t::Float64)
    lo, hi = corr.layers[k]
    corr.gain[(lo+1):(hi-1)] .= (hi - lo)^2 / τ
    corr.time[k] = τ
    corr.set_at[k] = t
    corr.stage[k] = 2
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
activation, at the fallback time; a correction with one fixed time is
unchanged.
"""
function reset_layer_times!(corr::LayerCorrection)
    est = corr.estimate
    est === nothing && return corr
    for k in eachindex(corr.layers)
        _set_layer_time!(corr, k, est.fallback, NaN)
        corr.stage[k] = 0
    end
    return corr
end

# The live cells of layer `k`, at most `lo:(hi − 1)`.
function _live_cells(corr::LayerCorrection, layout::StateLayout, k::Int)
    lo, hi = corr.layers[k]
    top = lo - 1
    for i = lo:(hi-1)
        layout.mask[energy_index(layout, i)] || break
        top = i
    end
    return lo:top
end

"""
	layer_time_margin(correction, layout, reference, u) -> Float64

The smallest margin `hᵢ/(r hᵢ(0)) − 1` over the live cells of the layers
whose equidistribution time is still to be estimated, `r` the activation's
`r_hi` before a layer's first activation and `r_lo` while a floored rate waits
for the second; `Inf` when none is waiting. It crosses zero downward where
[`estimate_layer_times!`](@ref) runs.
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
	estimate_layer_times!(correction, layout, reference, u, t, horizon; slack) -> correction

At the state `u` at time `t`, every waiting layer whose margin
([`layer_time_margin`](@ref)) is at most `slack` takes its equidistribution
time from its thickness `h` (the distance between its bounding nodes, the
outer one its outermost live node) and its thinning rate `ḣ` (their velocity
difference) as [`ActivationTime`](@ref) states, the floor being the
[`ActivationTime`](@ref)'s `floor` times the layer's initial thickness over
`horizon` [s].
"""
function estimate_layer_times!(
    corr::LayerCorrection,
    layout::StateLayout,
    reference::CellReference,
    u::AbstractVector,
    t::Real,
    horizon::Real;
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
        lo0, hi0 = corr.layers[k]
        floor = est.floor * (reference.positions[hi0] - reference.positions[lo0]) / horizon
        if ḣ < -floor
            _set_layer_time!(corr, k, est.λ * h / -ḣ, Float64(t))
        elseif ḣ <= 0 && stage == 0
            corr.stage[k] = 1
        else
            _set_layer_time!(corr, k, est.fallback, Float64(t))
        end
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
@inline function activation(corr::LayerCorrection, h, h0)
    h >= corr.r_hi * h0 && return zero(h)
    r = smooth_floor(h, DISTANCE_FLOOR_FRACTION * h0) / h0
    x = (corr.log_hi - log(r)) / corr.log_span
    x >= 1 && return one(x)
    return x * x * (3 - 2x)
end

# The normalized monitor of cell `i`; one without the mass term.
@inline function _monitor(corr::LayerCorrection{NC}, ds::DerivedState{NC}, i) where {NC}
    iszero(corr.β) && return one(eltype(ds.m))
    ρ = zero(eltype(ds.m))
    @inbounds for j = 1:NC
        corr.condensed[j] && (ρ += ds.m[j, i])
    end
    return (corr.α + corr.β * ρ / ds.V_eval[i]) * corr.monitor_scale[i]
end

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
    a_in = activation(corr, h_in, h0[1])
    @inbounds for nd = 2:(nn-1)
        h_out = x[nd+1] - x[nd]
        a_out = activation(corr, h_out, h0[nd])
        a = 1 - (1 - a_in) * (1 - a_out)
        g = corr.gain[nd]
        c[nd] =
            iszero(g) || !(ds.active[nd-1] && ds.active[nd]) || iszero(a) ? zero(T) :
            a * g * (_monitor(corr, ds, nd) * h_out - _monitor(corr, ds, nd - 1) * h_in)
        h_in, a_in = h_out, a_out
    end
    return nothing
end
