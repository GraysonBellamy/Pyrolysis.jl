# ============================================================================
# Conduction
# ============================================================================

# Smooth floor on the center-to-face distances, as a fraction of their
# reference values, so an inverted trial mesh keeps finite fluxes.
const DISTANCE_FLOOR_FRACTION = 1e-8

# The Tikhonov term of a cell's least-squares gradient, relative to the trace
# of its normal matrix, where its live neighbours do not span every direction
# (a column one cell wide, a cell beside dead ones): the unresolved component
# of the gradient is then zero instead of undetermined. The neighbours span
# when the normal matrix's determinant exceeds GRADIENT_SPAN times its mean
# eigenvalue to the power D.
const GRADIENT_REGULARIZATION = 1e-8
const GRADIENT_SPAN = 1e-6

"""
	ConductionKernel{A} <: Kernel

Fourier conduction across the interior faces between live cells. The face's
conductivity tensor `K_f` is the distance-weighted harmonic mean `k_f` of the
cell conductivities, with the smooth face-mean guard, times the
distance-weighted mean of the cells' dimensionless tensors (`tensor`, one
`SMatrix{D,D}` per cell; `nothing` on a line, where the tensor is the scalar
conductivity). With `d` the vector between the cell centers and `d_o`, `d_n`
the distances of the owner's and neighbour's centers from the face along its
normal, the flux along the normal splits into a two-point normal part and a
cross-diffusion part,

	q_f = −(n·K_f·n) (T_nb − T_o) / (d_o + d_n) − t·∇T_f,
	t = K_f n − (n·K_f·n) d / (d_o + d_n),

`t` tangential to the face and `∇T_f` the distance-weighted mean of the two
cells' least-squares gradients over their live face neighbours, so a linear
temperature field's flux is exact on any mesh where the neighbours span
every direction; where they do not, the gradient's unresolved component is
zero. Both parts are evaluated in every residual, so the Jacobian carries
both. On a line `d` is along the normal, the cross part vanishes identically
and is not evaluated.

`d_floor` holds the per-face floors of the two center-to-face distances.
"""
struct ConductionKernel{A} <: Kernel
    d_floor::Vector{Float64}
    tensor::A
end

"""
	ConductionKernel(topology, geometry; tensor = nothing)

The conduction kernel with distance floors taken from the reference
geometry; `tensor` gives each cell's dimensionless conductivity tensor on a
mesh of more than one dimension (the identity by default).
"""
function ConductionKernel(
    topo::MeshTopology{D},
    g::MeshGeometry;
    tensor = nothing,
) where {D}
    d = zeros(n_faces(topo))
    for f = 1:n_faces(topo)
        o = topo.faces[f].cells[1]
        d[f] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, f, o)
    end
    D == 1 && tensor === nothing && return ConductionKernel(d, nothing)
    D == 1 && throw(ArgumentError("a line's conductivity is the scalar of its cells"))
    K = tensor === nothing ? fill(SMatrix{D,D,Float64}(I), n_cells(topo)) : tensor
    length(K) == n_cells(topo) || throw(DimensionMismatch("one tensor per cell"))
    for (i, Ki) in enumerate(K)
        Ki isa SMatrix{D,D} && issymmetric(Ki) && isposdef(Matrix(Ki)) || throw(
            ArgumentError(
                "the tensor of cell $i must be a symmetric positive-definite SMatrix",
            ),
        )
    end
    return ConductionKernel(d, Vector{SMatrix{D,D,Float64,D * D}}(K))
end

# Normal component n·K·n of a conductivity tensor at a face.
@inline normal_conductivity(K::SMatrix{1,1}, n::SVector{1}) = K[1, 1] * n[1] * n[1]
@inline normal_conductivity(K::SMatrix{D,D}, n::SVector{D}) where {D} = dot(n, K * n)

"""
	cell_gradient(c, Val(D), i) -> SVector{D}

The least-squares temperature gradient of cell `i` that the conduction kernel
wrote into the contributions `c` (a mesh of more than one dimension).
"""
@inline cell_gradient(c::Contributions, ::Val{D}, i::Int) where {D} =
    SVector{D}(ntuple(k -> c.gradient[k, i], Val(D)))

function contribute!(
    c::Contributions,
    kernel::ConductionKernel,
    ds::DerivedState,
    topo::MeshTopology{D},
    material,
    t,
) where {D}
    g = ds.geometry
    _gradients!(c, kernel.tensor, kernel, ds, topo)
    @inbounds for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        (nb == 0 || !(ds.active[o] && ds.active[nb])) && continue
        floor = kernel.d_floor[f]
        d_o = smooth_floor(center_to_face_distance(g, topo, f, o), floor)
        d_n = smooth_floor(center_to_face_distance(g, topo, f, nb), floor)
        k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
        K = _face_tensor(kernel.tensor, k_f, o, nb, d_o, d_n)
        c.F_energy[f] -=
            normal_conductivity(K, g.face_normal[f]) * (ds.T[nb] - ds.T[o]) / (d_o + d_n)
        _cross_flux!(c, kernel.tensor, K, ds, g, f, o, nb, d_o, d_n)
    end
    return c
end

# The face's conductivity tensor: the scalar face conductivity on a line, times
# the distance-weighted mean of the cells' dimensionless tensors otherwise.
@inline _face_tensor(::Nothing, k_f, o, nb, d_o, d_n) = SMatrix{1,1}(k_f)
@inline _face_tensor(A::AbstractVector, k_f, o, nb, d_o, d_n) =
    k_f * (d_n * A[o] + d_o * A[nb]) / (d_o + d_n)

# The cross-diffusion part of the face flux, `t·∇T_f`: nothing on a line.
@inline _cross_flux!(c, ::Nothing, K, ds, g, f, o, nb, d_o, d_n) = nothing
@inline function _cross_flux!(
    c,
    ::AbstractVector,
    K::SMatrix{D,D},
    ds,
    g,
    f,
    o,
    nb,
    d_o,
    d_n,
) where {D}
    n = g.face_normal[f]
    d = g.cell_center[nb] - g.cell_center[o]
    tf = K * n - normal_conductivity(K, n) * d / (d_o + d_n)
    ∇T =
        (d_n * cell_gradient(c, Val(D), o) + d_o * cell_gradient(c, Val(D), nb)) /
        (d_o + d_n)
    c.F_energy[f] -= dot(tf, ∇T)
    return nothing
end

# The least-squares gradient of every live cell over its live face neighbours,
# each difference weighted by the inverse squared distance:
# (Σ d dᵀ/|d|² + R) ∇T = Σ d (T_j − T_i)/|d|², with R the Tikhonov term where
# the neighbours do not span every direction.
@inline _gradients!(c, ::Nothing, kernel, ds, topo) = nothing

# The Tikhonov term of a normal matrix whose neighbours do not span every
# direction, decided on its plain values: a property of the geometry.
@inline function _regularization(M::SMatrix{D,D}) where {D}
    Mv = value.(M)
    t = tr(Mv)
    (t <= 0 || det(Mv) > GRADIENT_SPAN * (t / D)^D) && return zero(M)
    return GRADIENT_REGULARIZATION * tr(M) * I
end
function _gradients!(c, ::AbstractVector, kernel, ds, topo::MeshTopology{D}) where {D}
    g = ds.geometry
    T = eltype(c.gradient)
    @inbounds for i = 1:n_cells(topo)
        M = zero(SMatrix{D,D,T})
        b = zero(SVector{D,T})
        if ds.active[i]
            for f in topo.cells[i].faces
                o, nb = topo.faces[f].cells
                j = o == i ? nb : o
                (j == 0 || !ds.active[j]) && continue
                d = g.cell_center[j] - g.cell_center[i]
                w = 1 / dot(d, d)
                M += w * (d * d')
                b += w * (ds.T[j] - ds.T[i]) * d
            end
            M += _regularization(M)
        end
        ∇T = abs(det(M)) > 0 ? M \ b : zero(SVector{D,T})
        for k = 1:D
            c.gradient[k, i] = ∇T[k]
        end
    end
    return nothing
end
