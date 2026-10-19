# ============================================================================
# Conduction
# ============================================================================

# Smooth floor on the center-to-face distances, as a fraction of their
# reference values, so an inverted trial mesh keeps finite fluxes.
const DISTANCE_FLOOR_FRACTION = 1e-8

"""
	ConductionKernel <: Kernel

Fourier conduction across the interior faces between live cells: the two-point flux
`q_f = −(n·K_f·n) (T_nb − T_o) / d` along the face normal, with `d` the
distance between the cell centers and `K_f` the conductivity tensor of the
face, built from the distance-weighted harmonic mean of the cell values with
the smooth face-mean guard. In 1D the tensor is the scalar conductivity and
the cross-diffusion part vanishes identically.

`d_floor` holds the per-face floors of the two center-to-face distances.
"""
struct ConductionKernel <: Kernel
    d_floor::Vector{Float64}
end

"""
	ConductionKernel(topology, geometry)

The conduction kernel with distance floors taken from the reference geometry.
"""
function ConductionKernel(topo::MeshTopology, g::MeshGeometry)
    d = zeros(n_faces(topo))
    for f = 1:n_faces(topo)
        o = topo.faces[f].cells[1]
        d[f] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, f, o)
    end
    return ConductionKernel(d)
end

# Normal component n·K·n of a conductivity tensor at a face.
@inline normal_conductivity(K::SMatrix{1,1}, n::SVector{1}) = K[1, 1] * n[1] * n[1]

function contribute!(
    c::Contributions,
    kernel::ConductionKernel,
    ds::DerivedState,
    topo::MeshTopology{1},
    material,
    t,
)
    g = ds.geometry
    @inbounds for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        (nb == 0 || !(ds.active[o] && ds.active[nb])) && continue
        floor = kernel.d_floor[f]
        d_o = smooth_floor(center_to_face_distance(g, topo, f, o), floor)
        d_n = smooth_floor(center_to_face_distance(g, topo, f, nb), floor)
        k_f = guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
        K = SMatrix{1,1}(k_f)
        c.F_energy[f] -=
            normal_conductivity(K, g.face_normal[f]) * (ds.T[nb] - ds.T[o]) / (d_o + d_n)
    end
    return c
end
