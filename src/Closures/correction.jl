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

The default equidistribution time `τ` [s] of a layer: with it the 5.9 mm
PMMA layer of V-BO-01 at 50 kW/m² ends in one deactivation at 30, 60 and 120
cells across it.
"""
const EQUIDISTRIBUTION_TIME = 720.0

"""
	EquidistributionCorrection(; τ = EQUIDISTRIBUTION_TIME, activation = (0.3, 0.1),
	                           monitor = (1.0, 0.0), layers = :vanishing) <: Correction

The relaxed equidistribution of a vanishing layer's cells. At an interior
node `n` of a layer of `N` cells, between cell `n − 1` (inward) and cell `n`,

	cₙ = aₙ N² (M̃ₙ hₙ − M̃ₙ₋₁ hₙ₋₁) / τ,

which moves the node into the thicker (or heavier) cell, so `M̃ h` relaxes
toward a common value over the layer in about `τ/π²` whatever `N`: `τ` [s] is
the layer's equidistribution time. The activation `aₙ = 1 − (1 − aₙ₋₁)(1 − aₙ)`
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
struct EquidistributionCorrection <: Correction
    τ::Float64
    r_hi::Float64
    r_lo::Float64
    α::Float64
    β::Float64
    layers::Symbol
end

function EquidistributionCorrection(;
    τ::Real = EQUIDISTRIBUTION_TIME,
    activation::Tuple{Real,Real} = (0.3, 0.1),
    monitor::Tuple{Real,Real} = (1.0, 0.0),
    layers::Symbol = :vanishing,
)
    r_hi, r_lo = activation
    α, β = monitor
    τ > 0 || throw(ArgumentError("the equidistribution time must be positive, got $τ"))
    0 < r_lo < r_hi || throw(
        ArgumentError(
            "the activation needs 0 < r_lo < r_hi, got (r_hi, r_lo) = $activation",
        ),
    )
    (α > 0 && β >= 0) ||
        throw(ArgumentError("the monitor needs α > 0 and β ≥ 0, got (α, β) = $monitor"))
    layers in (:vanishing, :all) ||
        throw(ArgumentError("layers must be :vanishing or :all, got $(repr(layers))"))
    return EquidistributionCorrection(τ, r_hi, r_lo, α, β, layers)
end

"""
	LayerCorrection{NC} <: Correction

An [`EquidistributionCorrection`](@ref) bound to a mesh: per node the gain
`N²/τ` of its layer (zero at the layers' bounding nodes), per cell the
reciprocal of its layer's mean initial monitor, the activation's upper bound
`r_hi` with `ln r_hi` and `ln r_hi − ln r_lo`, the monitor's `α` and `β/ρ_ref`,
and the condensed component flags.
"""
struct LayerCorrection{NC} <: Correction
    gain::Vector{Float64}
    monitor_scale::Vector{Float64}
    r_hi::Float64
    log_hi::Float64
    log_span::Float64
    α::Float64
    β::Float64
    condensed::NTuple{NC,Bool}
end

"""
	vanishing(material) -> NTuple{NC,Bool}

Whether each condensed component of `material` gasifies completely: some
reaction consumes it, and every reaction that does yields no condensed
component that does not itself gasify completely. The least such set, so an
inert component, a residue that no reaction consumes and the components of a
cycle do not vanish; a gas component reads `false`.
"""
function vanishing(material::Material{NC}) where {NC}
    condensed = ntuple(j -> !is_gas(material.components[j]), Val(NC))
    consumes(rxn, j) = any(r -> first(r) == j, rxn.reactants)
    v = fill(false, NC)
    changed = true
    while changed
        changed = false
        for j = 1:NC
            (condensed[j] && !v[j]) || continue
            consumers = [rxn for rxn in material.reactions if consumes(rxn, j)]
            isempty(consumers) && continue
            gasified = p -> !condensed[first(p)] || iszero(last(p)) || v[first(p)]
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
    for (lo, hi) in zip(bounds[1:(end-1)], bounds[2:end])
        N = hi - lo
        moves = corr.layers === :all || vanishes(v, m0, lo:(hi-1), condensed)
        moves && (gain[(lo+1):(hi-1)] .= N^2 / corr.τ)
        M̄ = sum(corr.α + β * ρ0[i] for i = lo:(hi-1)) / N
        scale[lo:(hi-1)] .= 1 / M̄
    end
    return LayerCorrection(
        gain,
        scale,
        corr.r_hi,
        log(corr.r_hi),
        log(corr.r_hi) - log(corr.r_lo),
        corr.α,
        β,
        condensed,
    )
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
