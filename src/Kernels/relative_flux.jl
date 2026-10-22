# ============================================================================
# Inventory fluxes relative to a moving mesh
# ============================================================================

"""
	van_albada(a, b, ε)

The smooth van Albada slope of two one-sided gradients `a` (upwind) and `b`
(downwind): `(a (b² + ε) + b (a² + ε)) / (a² + b² + 2ε)`. It equals `a` when
`a = b`, is zero on a uniform field and at `a = −b`, stays below
`0.21 max(|a|, |b|)` in magnitude when they differ in sign, and is C^∞ for
`ε > 0`.
"""
@inline van_albada(a, b, ε) = (a * (b * b + ε) + b * (a * a + ε)) / (a * a + b * b + 2ε)

"""
	RelativeFluxKernel{L,NC} <: Kernel

The inventories a moving mesh carries across its faces when the mesh velocity
`w` differs from the solid velocity `vₛ`. Through each face, along its normal,
every component moves with the solid at the relative velocity
`u = (vₛ − w)·n = −c·n`, with `c` the node correction of the derived state
(`ds.c`), its mass flux `Fⱼ = u ξⱼ,f`; the conversion-history pair
moves at the solid mass flux times its per-solid-mass values; and the energy
flux gains `Σⱼ Fⱼ hⱼ(T_f)`, the sensible enthalpy of the same mass fluxes at
the face temperature `T_f` reconstructed from the same cells, so an
isothermal composition profile advected across the mesh stays isothermal.

Face states are second-order upwind reconstructions: the upwind cell's value
plus its smooth [`van_albada`](@ref) slope times its distance to the face,
the slope built from the gradient from the cell beyond it and the gradient
toward the downwind cell, the gradient from beyond taken as zero where no
live cell lies there: beside a resolved gradient the face state is then the
upwind cell's own value, so a boundary cell holding almost none of a
component is not drained. A boundary face takes the one-sided linear
extrapolation of its cell whichever way the mesh crosses it, so its flux has
no switch. The face state
of an interior face changes side with the sign of `u`; the flux `u ξ_f` is
continuous there, with a kink. The side follows the value of `u` alone, so
at `u = 0` every derivative pass takes the derivative of the `u ≥ 0` side,
whatever its seeds.

A concentration's face state is bounded by the upwind cell's inventory: with
`p` the cell's concentration, smoothly floored at zero below its
concentration tolerance `η` ([`positive_part`](@ref)), the reconstruction's
increment `δ` over the cell value becomes `δ p / √(p² + δ² + η²)`. The face
state lies in `[0, 2p]`, differs from the reconstruction by `O(δ³/p²)` where
`|δ| ≪ p`, and is exactly zero for an empty cell, so no face draws a
component out of a cell that holds none of it, whatever the slopes around
it: the inventories' semi-discrete equations keep them nonnegative. The
temperature and the history values are not inventories and keep the
reconstruction.

A face whose relative velocity is exactly zero carries nothing and is
skipped. The kernel runs after the thermal boundary kernel, so the surface
balance and its blowing models see the gas-transport fluxes alone.

`solid` flags the components that carry the history pair; `ε` holds, per
cell, the slope scales of the concentrations, the temperature and the
per-solid-mass history values (each the square of a tolerance over the
reference cell thickness); `η` the concentration tolerance of each cell;
`far` and `beyond` give, per face and side, the
cell beyond the face's cell on that side and the face between them (zero
where there is none); `d_floor` the floors of the center-to-face distances
of each face's owner and neighbour.
"""
struct RelativeFluxKernel{L<:CaloricLaws,NC} <: Kernel
    laws::L
    solid::NTuple{NC,Bool}
    history::Bool
    ε::Matrix{Float64}
    η::Vector{Float64}
    far::Matrix{Int}
    beyond::Matrix{Int}
    d_floor::Matrix{Float64}
    mass_scale::Vector{Float64}
end

"""
	RelativeFluxKernel(topology, geometry, material, laws, reference;
	                   history, temperature_scale, history_scale)

The kernel with its slope scales: the concentration tolerance of each cell
(`reference.mass_scale / reference.volume`), `temperature_scale` [K] and
`history_scale` (dimensionless), each over the cell's reference thickness;
`history` carries the conversion-history pair.
"""
function RelativeFluxKernel(
    topo::MeshTopology{1},
    g::MeshGeometry,
    material::Material,
    laws::CaloricLaws,
    reference::CellReference;
    history::Bool,
    temperature_scale::Float64,
    history_scale::Float64,
)
    n = n_cells(topo)
    nf = n_faces(topo)
    ε = zeros(3, n)
    η = reference.mass_scale ./ reference.volume
    for i = 1:n
        f1, f2 = topo.cells[i].faces
        h =
            center_to_face_distance(g, topo, f1, i) +
            center_to_face_distance(g, topo, f2, i)
        ε[1, i] = (η[i] / h)^2
        ε[2, i] = (temperature_scale / h)^2
        ε[3, i] = (history_scale / h)^2
    end
    far = zeros(Int, 2, nf)
    beyond = zeros(Int, 2, nf)
    d = zeros(2, nf)
    for f = 1:nf
        for side = 1:2
            c = topo.faces[f].cells[side]
            c == 0 && continue
            d[side, f] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, f, c)
            cf = topo.cells[c].faces
            other = cf[1] == f ? cf[2] : cf[1]
            o, nb = topo.faces[other].cells
            beyond[side, f] = other
            far[side, f] = o == c ? nb : o
        end
    end
    return RelativeFluxKernel(
        laws,
        map(is_solid, material.components),
        history,
        ε,
        η,
        far,
        beyond,
        d,
        copy(reference.mass_scale),
    )
end

# Smoothly floored distance from the center of the cell on `side` of face `f`
# to the face.
@inline function _side_distance(kernel::RelativeFluxKernel, g, topo, f, side)
    c = topo.faces[f].cells[side]
    return smooth_floor(center_to_face_distance(g, topo, f, c), kernel.d_floor[side, f])
end

# Distance between the centers of the two cells of interior face `f`.
@inline _across(kernel::RelativeFluxKernel, g, topo, f) =
    _side_distance(kernel, g, topo, f, 1) + _side_distance(kernel, g, topo, f, 2)

# The cells and distances of a face state: the upwind cell `U`, the downwind
# cell `D` (zero on a boundary face), the cell `X` whose gradient into `U`
# sets the slope (zero where there is none or where it is dead, `active` being
# the live cells), `d_f` from `U` to the face, `d_UD` and `d_XU`. On a
# boundary face `X` is the boundary cell's inner neighbour and the state is
# its one-sided extrapolation. The side is read from the value of `u`: a dual
# number of zero value compares by its partials, which differ between seeds.
@inline function _stencil(kernel::RelativeFluxKernel, g, topo, f, u, active)
    o, nb = topo.faces[f].cells
    side = nb == 0 || value(u) >= 0 ? 1 : 2
    U = side == 1 ? o : nb
    D = nb == 0 ? 0 : (side == 1 ? nb : o)
    X = kernel.far[side, f]
    X = X != 0 && !active[X] ? 0 : X
    d_f = _side_distance(kernel, g, topo, f, side)
    d_UD = D == 0 ? zero(d_f) : _across(kernel, g, topo, f)
    d_XU = X == 0 ? zero(d_f) : _across(kernel, g, topo, kernel.beyond[side, f])
    return U, D, X, d_f, d_UD, d_XU
end

# The face value of a cell quantity `q(i)` on the stencil, with slope scale ε;
# with no cell beyond the upwind cell the gradient from beyond is zero.
@inline function _face_value(q::F, U, D, X, d_f, d_UD, d_XU, ε) where {F}
    qU = q(U)
    if D == 0
        X == 0 && return qU
        return qU + d_f * (qU - q(X)) / d_XU
    end
    b = (q(D) - qU) / d_UD
    X == 0 && return qU + d_f * van_albada(zero(b), b, ε)
    a = (qU - q(X)) / d_XU
    return qU + d_f * van_albada(a, b, ε)
end

# The face value of a concentration `q(i)`: the face value's increment over
# the upwind cell's value, scaled by the cell's inventory `p` against the
# concentration tolerance η, so the state lies in [0, 2p].
@inline function _inventory_face_value(q::F, U, D, X, d_f, d_UD, d_XU, ε, η) where {F}
    qU = q(U)
    δ = _face_value(q, U, D, X, d_f, d_UD, d_XU, ε) - qU
    p = positive_part(qU, η)
    return p + δ * p / sqrt(p * p + δ * δ + η * η)
end

# Solid mass of cell `i`, smoothly floored at its mass scale, the carrier of
# the per-solid-mass history values.
@inline function _solid_mass(kernel::RelativeFluxKernel{L,NC}, ds, i) where {L,NC}
    m = zero(eltype(ds.m))
    @inbounds for j = 1:NC
        kernel.solid[j] && (m += ds.m[j, i])
    end
    return smooth_floor(m, kernel.mass_scale[i])
end

function contribute!(
    c::Contributions{NC},
    kernel::RelativeFluxKernel{L,NC},
    ds::DerivedState{NC},
    topo::MeshTopology{1},
    material,
    t,
) where {L,NC}
    g = ds.geometry
    @inbounds for f = 1:n_faces(topo)
        u = -ds.c[topo.faces[f].nodes[1]] * g.face_normal[f][1]
        iszero(u) && continue
        s = _stencil(kernel, g, topo, f, u, ds.active)
        U = s[1]
        T_f = _face_value(i -> ds.T[i], s..., kernel.ε[2, U])
        h = component_enthalpies(kernel.laws, T_f)
        q = zero(eltype(c.F_energy))
        F_solid = zero(q)
        for j = 1:NC
            ξ_f = _inventory_face_value(i -> ds.ξ[j, i], s..., kernel.ε[1, U], kernel.η[U])
            F = u * ξ_f
            c.F_mass[j, f] += F
            q += F * h[j]
            kernel.solid[j] && (F_solid += F)
        end
        c.F_energy[f] += q
        kernel.history || continue
        ε = kernel.ε[3, U]
        g_f = _face_value(i -> ds.G[i] / _solid_mass(kernel, ds, i), s..., ε)
        μ_f = _face_value(i -> ds.M[i] / _solid_mass(kernel, ds, i), s..., ε)
        c.F_carried[1, f] += F_solid * g_f
        c.F_carried[2, f] += F_solid * μ_f
    end
    return c
end

"""
	RelativeFluxKernel(topology::MeshTopology{D}, geometry, material, laws, reference;
	                   history, temperature_scale, history_scale)    (D ≥ 2)

The relative fluxes of a mesh of more than one dimension, first-order upwind:
through each face every inventory moves at the solid's velocity relative to
the face, `u A_f = −Σ_q s_{q,f} · c_q` with `c` the node corrections and `s`
the face's sweep vectors, carrying the upwind cell's concentrations, their
enthalpy at its temperature and its per-solid-mass history values; a boundary
face carries its own cell's whichever way the mesh crosses it. A
concentration's face state is the upwind cell's positive part at its
concentration tolerance, so no face draws a component out of a cell that
holds none of it. The slope scales do not enter.
"""
function RelativeFluxKernel(
    topo::MeshTopology,
    g::MeshGeometry,
    material::Material,
    laws::CaloricLaws,
    reference::CellReference;
    history::Bool,
    temperature_scale::Float64,
    history_scale::Float64,
)
    return RelativeFluxKernel(
        laws,
        map(is_solid, material.components),
        history,
        zeros(3, 0),
        reference.mass_scale ./ reference.volume,
        zeros(Int, 2, 0),
        zeros(Int, 2, 0),
        zeros(2, 0),
        copy(reference.mass_scale),
    )
end

function contribute!(
    c::Contributions{NC},
    kernel::RelativeFluxKernel{L,NC},
    ds::DerivedState{NC},
    topo::MeshTopology{D},
    material,
    t,
) where {L,NC,D}
    g = ds.geometry
    @inbounds for f = 1:n_faces(topo)
        nodes = topo.faces[f].nodes
        sweep = zero(eltype(c.F_energy))
        for q in eachindex(nodes)
            cq = SVector{D}(ntuple(k -> ds.c[D*(nodes[q]-1)+k], Val(D)))
            sweep += dot(g.face_sweep[q, f], cq)
        end
        iszero(sweep) && continue
        u = -sweep / g.face_area[f]
        o, nb = topo.faces[f].cells
        U = nb == 0 || value(u) >= 0 ? o : nb
        h = component_enthalpies(kernel.laws, ds.T[U])
        q = zero(eltype(c.F_energy))
        F_solid = zero(q)
        for j = 1:NC
            F = u * positive_part(ds.ξ[j, U], kernel.η[U])
            c.F_mass[j, f] += F
            q += F * h[j]
            kernel.solid[j] && (F_solid += F)
        end
        c.F_energy[f] += q
        kernel.history || continue
        m = _solid_mass(kernel, ds, U)
        c.F_carried[1, f] += F_solid * ds.G[U] / m
        c.F_carried[2, f] += F_solid * ds.M[U] / m
    end
    return c
end
