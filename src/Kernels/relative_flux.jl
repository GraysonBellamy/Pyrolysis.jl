# ============================================================================
# Inventory fluxes relative to a moving mesh
# ============================================================================

"""
	van_albada(a, b, ε)

The smooth van Albada slope of two one-sided gradients `a` (upwind) and `b`
(downwind): `(a (b² + ε) + b (a² + ε)) / (a² + b² + 2ε)`. It equals `a` when
`a = b`, nearly vanishes when they differ in sign, is zero on a uniform field
and is C^∞ for `ε > 0`.
"""
@inline van_albada(a, b, ε) = (a * (b * b + ε) + b * (a * a + ε)) / (a * a + b * b + 2ε)

"""
	RelativeFluxKernel{L,NC} <: Kernel

The inventories a moving mesh carries across its faces when the mesh velocity
`w` differs from the solid velocity `vₛ`. Through each face, along its normal,
every component moves with the solid at the relative velocity
`u = (vₛ − w)·n`, its mass flux `Fⱼ = u ξⱼ,f`; the conversion-history pair
moves at the solid mass flux times its per-solid-mass values; and the energy
flux gains `Σⱼ Fⱼ hⱼ(T_f)`, the sensible enthalpy of the same mass fluxes at
the face temperature `T_f` reconstructed from the same cells, so an
isothermal composition profile advected across the mesh stays isothermal.

Face states are second-order upwind reconstructions: the upwind cell's value
plus its smooth [`van_albada`](@ref) slope times its distance to the face,
the slope built from the gradient from the cell beyond it and the gradient
toward the downwind cell (the downwind gradient alone where no cell lies
beyond). A boundary face takes the one-sided linear extrapolation of its cell
whichever way the mesh crosses it, so its flux has no switch. The face state
of an interior face changes side with the sign of `u`; the flux `u ξ_f` is
continuous there, with a kink.

The solid is at rest (`vₛ = 0`): the kernel runs under a prescribed mesh
velocity, on a material whose reactions do not change its volume. It runs
after the thermal boundary kernel, so the surface balance and its blowing
models see the gas-transport fluxes alone.

`solid` flags the components that carry the history pair; `ε` holds, per
cell, the slope scales of the concentrations, the temperature and the
per-solid-mass history values (each the square of a tolerance over the
reference cell thickness); `far` and `beyond` give, per face and side, the
cell beyond the face's cell on that side and the face between them (zero
where there is none); `d_floor` the floors of the center-to-face distances
of each face's owner and neighbour.
"""
struct RelativeFluxKernel{L<:CaloricLaws,NC} <: Kernel
    laws::L
    solid::NTuple{NC,Bool}
    history::Bool
    ε::Matrix{Float64}
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
    for i = 1:n
        f1, f2 = topo.cells[i].faces
        h =
            center_to_face_distance(g, topo, f1, i) +
            center_to_face_distance(g, topo, f2, i)
        ε[1, i] = (reference.mass_scale[i] / reference.volume[i] / h)^2
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
# sets the slope (zero where there is none), `d_f` from `U` to the face,
# `d_UD` and `d_XU`. On a boundary face `X` is the boundary cell's inner
# neighbour and the state is its one-sided extrapolation.
@inline function _stencil(kernel::RelativeFluxKernel, g, topo, f, u)
    o, nb = topo.faces[f].cells
    side = nb == 0 || u >= 0 ? 1 : 2
    U = side == 1 ? o : nb
    D = nb == 0 ? 0 : (side == 1 ? nb : o)
    X = kernel.far[side, f]
    d_f = _side_distance(kernel, g, topo, f, side)
    d_UD = D == 0 ? zero(d_f) : _across(kernel, g, topo, f)
    d_XU = X == 0 ? zero(d_f) : _across(kernel, g, topo, kernel.beyond[side, f])
    return U, D, X, d_f, d_UD, d_XU
end

# The face value of a cell quantity `q(i)` on the stencil, with slope scale ε.
@inline function _face_value(q::F, U, D, X, d_f, d_UD, d_XU, ε) where {F}
    qU = q(U)
    if D == 0
        X == 0 && return qU
        return qU + d_f * (qU - q(X)) / d_XU
    end
    b = (q(D) - qU) / d_UD
    X == 0 && return qU + d_f * b
    a = (qU - q(X)) / d_XU
    return qU + d_f * van_albada(a, b, ε)
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
        u = -ds.w[topo.faces[f].nodes[1]] * g.face_normal[f][1]
        s = _stencil(kernel, g, topo, f, u)
        U = s[1]
        T_f = _face_value(i -> ds.T[i], s..., kernel.ε[2, U])
        h = component_enthalpies(kernel.laws, T_f)
        q = zero(eltype(c.F_energy))
        F_solid = zero(q)
        for j = 1:NC
            F = u * _face_value(i -> ds.ξ[j, i], s..., kernel.ε[1, U])
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
