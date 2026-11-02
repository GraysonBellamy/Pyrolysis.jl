# ============================================================================
# Conduction
# ============================================================================

# Smooth floor on the center-to-face distances, as a fraction of their
# reference values, so an inverted trial mesh keeps finite fluxes.
const DISTANCE_FLOOR_FRACTION = 1e-8

# The Tikhonov term of a cell's least-squares gradient where its neighbours'
# images do not span every direction (a column one cell wide, a cell beside
# dead ones, neighbour tensors that make the images collinear): the
# unresolved component of the gradient is then zero instead of undetermined.
# Its strength is GRADIENT_REGULARIZATION times the trace of the normal matrix
# M times a cubic smooth step of the normalized determinant
# x = det M / (tr M / D)^D, full at and below GRADIENT_SPAN and zero at and
# above GRADIENT_SUPPORT, so the gradient is continuous in the
# conductivities.
const GRADIENT_REGULARIZATION = 1e-8
const GRADIENT_SPAN = 1e-6
const GRADIENT_SUPPORT = 1e-4

# The normal conductivity of a face neighbour, relative to the cell's, at and
# below which the neighbour is insulating: its datum in the cell's gradient is
# the limit row (K_i n)·∇T = 0 of its transmission image.
const INSULATING_NEIGHBOUR = 1e-12

"""
	ConductionKernel{A,F} <: Kernel

Fourier conduction across the interior faces between live cells, with `d_o`,
`d_n` the distances of the owner's and the neighbour's centers from the face
along its normal and `δ = d_o + d_n`. On a line the face conductivity `k_f` is
the distance-weighted harmonic mean of the cell conductivities with the smooth
face-mean guard and `q_f = −k_f (T_nb − T_o) / δ`.

On a mesh of more than one dimension every cell has the complete tensor
`K_i = k_i A_i`, `A_i` its dimensionless conductivity tensor (`tensor`, one
`SMatrix{D,D}` per cell) and `k_i` its unfloored effective conductivity. The
flux is the transmission flux of the face, which carries a common temperature
trace and a common normal flux across it: with `a_i = n·K_i n`,
`b_i = P K_i n`, `P = I − n nᵀ` and `d` the vector between the cell centers,

	q_f = −a_f (T_nb − T_o) / δ − (b_f − a_f P d / δ)·∇T_f,
	a_f = δ a_o a_n / (d_o a_n + d_n a_o + η δ),
	b_f = (d_o a_n b_o + d_n a_o b_n) / (d_o a_n + d_n a_o + η δ),

`η` the face-mean guard and `∇T_f = (d_n ∇T_o + d_o ∇T_nb) / δ`, of which only
the tangential part enters. A cell's gradient `∇T_i` is the least-squares fit
over its live face neighbours, each taken at its transmission image through
the shared face `f`,

	d′ = (x_j − x_i) + (n·(x_j − x_f)) (K_i − K_j) n / (n·K_j n),

where a field affine on either side of the face, continuous in temperature and
normal flux across it, is affine as seen from cell `i`. A neighbour whose
normal conductivity is at most `INSULATING_NEIGHBOUR` times the cell's
contributes the limit of its image, the insulated-face row `(K_i n)·∇T = 0`
with unit weight.

The operator is exact for piecewise-affine fields across grid-fitted planar
interfaces wherever no cell's fit is regularized, its normalized determinant
at or above `GRADIENT_SUPPORT` (a uniform tensor's linear fields among them),
and second order on layered coefficients. It is not monotone for rotated
anisotropy: on unit squares at unit volumetric heat capacity an interior
cell's coefficients on the two diagonal neighbours across the cross coupling
are `−K₁₂/2`. Where the images miss a direction or barely span it, a Tikhonov
term that fades smoothly with the fit's conditioning makes the unresolved
component of the gradient small or zero ([`gradient_support`](@ref) counts
those cells). The flux goes to zero as either side's normal conductivity
vanishes. Every part is evaluated in every residual, so the Jacobian carries
the gradients and their images.

`d_floor` holds the per-face floors of the two center-to-face distances. `F`
is whether the mesh has more than one dimension and every node pinned: its
metrics then do not depend on the state, and the kernel reads their plain
values under every derivative pass.
"""
struct ConductionKernel{A,F} <: Kernel
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
    D == 1 && tensor === nothing && return ConductionKernel{Nothing,false}(d, nothing)
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
    A = Vector{SMatrix{D,D,Float64,D * D}}(K)
    return ConductionKernel{typeof(A),is_fixed(topo)}(d, A)
end

# A metric of the mesh as the kernel reads it: its plain value where it does
# not depend on the state.
@inline _metric(::ConductionKernel{A,true}, x) where {A} = primal.(x)
@inline _metric(::ConductionKernel, x) = x

# Normal component n·K·n of a conductivity tensor at a face.
@inline normal_conductivity(K::SMatrix{D,D}, n::SVector{D}) where {D} = dot(n, K * n)

# The complete conductivity tensor of cell `i`, and its product with the face
# normal `n`.
@inline cell_tensor(A::AbstractVector, ds, i) = ds.k_physical[i] * A[i]
@inline _flux_direction(A, ds, i, n) = ds.k_physical[i] * (A[i] * n)

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
    _gradients!(c, kernel, ds, topo)
    @inbounds for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        (nb == 0 || !(ds.active[o] && ds.active[nb])) && continue
        floor = kernel.d_floor[f]
        d_o = smooth_floor(_metric(kernel, center_to_face_distance(g, topo, f, o)), floor)
        d_n = smooth_floor(_metric(kernel, center_to_face_distance(g, topo, f, nb)), floor)
        c.F_energy[f] += _face_flux(kernel, kernel.tensor, c, ds, g, f, o, nb, d_o, d_n)
    end
    return c
end

# The two-point flux of the harmonic face conductivity on a line.
@inline function _face_flux(kernel, ::Nothing, c, ds, g, f, o, nb, d_o, d_n)
    n = g.face_normal[f]
    k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
    return -(k_f * n[1] * n[1] * (ds.T[nb] - ds.T[o]) / (d_o + d_n))
end

# The transmission flux of the complete tensors on either side of the face.
@inline function _face_flux(
    kernel,
    A::AbstractVector{<:SMatrix{D,D}},
    c,
    ds,
    g,
    f,
    o,
    nb,
    d_o,
    d_n,
) where {D}
    n = _metric(kernel, g.face_normal[f])
    Kn_o, Kn_n = _flux_direction(A, ds, o, n), _flux_direction(A, ds, nb, n)
    a_o, a_n = dot(n, Kn_o), dot(n, Kn_n)
    b_o, b_n = Kn_o - a_o * n, Kn_n - a_n * n
    δ = d_o + d_n
    a_f = guarded_harmonic_mean(a_o, a_n, d_o, d_n)
    b_f =
        (d_o * a_n * b_o + d_n * a_o * b_n) / (d_o * a_n + d_n * a_o + FACE_MEAN_GUARD * δ)
    d = _metric(kernel, g.cell_center[nb]) - _metric(kernel, g.cell_center[o])
    r = d - dot(n, d) * n
    ∇T = (d_n * cell_gradient(c, Val(D), o) + d_o * cell_gradient(c, Val(D), nb)) / δ
    return -a_f * (ds.T[nb] - ds.T[o]) / δ - dot(b_f - a_f * r / δ, ∇T)
end

# The least-squares gradient of every live cell over its live face
# neighbours' transmission images: (Σ w d′ d′ᵀ + R) ∇T = Σ w (T_j − T_i) d′,
# w = 1/|d′|², R the Tikhonov term.
@inline _gradients!(c, ::ConductionKernel{Nothing}, ds, topo) = nothing
function _gradients!(
    c,
    kernel::ConductionKernel{<:AbstractVector},
    ds,
    topo::MeshTopology{D},
) where {D}
    T = eltype(c.gradient)
    @inbounds for i = 1:n_cells(topo)
        M = zero(SMatrix{D,D,T})
        b = zero(SVector{D,T})
        if ds.active[i]
            M, b = _normal_equations(kernel, ds, topo, i, T)
            M += _regularization(M)
        end
        ∇T = abs(primal(det(M))) > 0 ? M \ b : zero(SVector{D,T})
        for k = 1:D
            c.gradient[k, i] = ∇T[k]
        end
    end
    return nothing
end

# The normal matrix and right-hand side of live cell `i`'s fit, element type `T`.
@inline function _normal_equations(
    kernel,
    ds,
    topo::MeshTopology{D},
    i,
    ::Type{T},
) where {D,T}
    M = zero(SMatrix{D,D,T})
    b = zero(SVector{D,T})
    for f in topo.cells[i].faces
        o, nb = topo.faces[f].cells
        j = o == i ? nb : o
        (j == 0 || !ds.active[j]) && continue
        Mj, bj = _datum(kernel, ds, f, i, j)
        M += Mj
        b += bj
    end
    return M, b
end

# The term of face neighbour `j` across face `f` in cell `i`'s fit: its
# transmission image, or the insulated-face row where it does not conduct
# normal to the face (none if neither cell does), the branch chosen on the
# plain values so that every derivative pass differentiates the plain one.
@inline function _datum(kernel, ds, f, i, j)
    g, A = ds.geometry, kernel.tensor
    n = _metric(kernel, g.face_normal[f])
    Kn_i, Kn_j = _flux_direction(A, ds, i, n), _flux_direction(A, ds, j, n)
    a_j = dot(n, Kn_j)
    if primal(a_j) <= INSULATING_NEIGHBOUR * primal(dot(n, Kn_i))
        ee = dot(Kn_i, Kn_i)
        return (primal(ee) > 0 ? (Kn_i * Kn_i') / ee : zero(Kn_i * Kn_i')), zero(Kn_i)
    end
    x_i, x_j = _metric(kernel, g.cell_center[i]), _metric(kernel, g.cell_center[j])
    s = dot(n, x_j - _metric(kernel, g.face_center[f]))
    d = x_j - x_i + s * (Kn_i - Kn_j) / a_j
    w = 1 / dot(d, d)
    return w * (d * d'), w * (ds.T[j] - ds.T[i]) * d
end

# The normalized determinant det M / (tr M / D)^D of a normal matrix: 1 for
# neighbours spread evenly over every direction, 0 where they miss one.
@inline function _conditioning(M::SMatrix{D,D}) where {D}
    t = tr(M)
    return primal(t) > 0 ? det(M) / (t / D)^D : zero(t)
end

# The Tikhonov term of a normal matrix, on its own values and derivatives; the
# step's branch from the plain values.
@inline function _regularization(M::SMatrix)
    x = _conditioning(M)
    (primal(tr(M)) <= 0 || primal(x) >= GRADIENT_SUPPORT) && return zero(M)
    σ =
        primal(x) <= GRADIENT_SPAN ? one(x) :
        _smooth_step((GRADIENT_SUPPORT - x) / (GRADIENT_SUPPORT - GRADIENT_SPAN))
    return GRADIENT_REGULARIZATION * tr(M) * σ * I
end
@inline _smooth_step(s) = s * s * (3 - 2s)

"""
	gradient_support(kernel::ConductionKernel, ds, topology) -> (; cells, worst)

The live cells whose least-squares gradient is a regularized fit, and so not
exact for piecewise-affine fields: `cells` counts those whose normal matrix
has a normalized determinant `det M / (tr M / D)^D` below `GRADIENT_SUPPORT`,
where the Tikhonov term starts, and `worst` is the smallest normalized
determinant over the live cells (`NaN` on a line or without live cells),
at the plain derived state `ds`.
"""
gradient_support(::ConductionKernel{Nothing}, ds, topo::MeshTopology) =
    (cells = 0, worst = NaN)
function gradient_support(
    kernel::ConductionKernel{<:AbstractVector},
    ds,
    topo::MeshTopology,
)
    cells, worst = 0, NaN
    for i = 1:n_cells(topo)
        ds.active[i] || continue
        M, _ = _normal_equations(kernel, ds, topo, i, eltype(ds.T))
        x = primal(_conditioning(M))
        cells += x < GRADIENT_SUPPORT
        worst = isnan(worst) ? x : min(worst, x)
    end
    return (; cells, worst)
end
