"""
    GasTransportKernel

One conservative species-flux vector per interior face, relative to the
skeleton. Total flux is matched through the two half cells. The flux carries
the gas enthalpies at the distance-weighted face temperature
`T_f = (d_n T_o + d_o T_n)/(d_o + d_n)`, plus the conduction
`D (T_o − T_n)/(d_o + d_n)` that keeps each cell's coupling to the other,
through the face's conduction and gas transport together, nonnegative:
`D = max(0, a_o − k_f, a_n − k_f)`, with `k_f` the face conductivity of
conduction and `a_i` the heat-capacity flux leaving cell `i` through the face
times its center-to-face distance, from the enthalpy secants between `T_f` and
the cell temperatures. Where conduction alone keeps both couplings
nonnegative `D` is zero and the flux is second-order accurate (with constant
heat capacities: up to a cell Péclet number of 2 on equal half cells); where
the gas flux dominates, the face becomes upwind and first order. Conduction and gas transport together then
create no temperature extremum, wherever every changing inventory lies above
the caloric relation's mass scale (see [`face_peclet`](@ref)). Boundary
exchange is assembled by the coupled thermal/mass boundary kernel.
"""
struct GasTransportKernel{M,NC,NG,I,L,S} <: Kernel
    model::M
    gas::NTuple{NC,Bool}
    indices::I
    d_floor::Vector{Float64}
    laws::L
    masses::S
end
function GasTransportKernel(
    model::GasTransport,
    topo::MeshTopology,
    g::MeshGeometry,
    bc_set::BoundaryConditionSet,
    material::Material{NC},
    laws::CaloricLaws,
) where {NC}
    indices=gas_component_indices(material.components)
    d=[
        DISTANCE_FLOOR_FRACTION*center_to_face_distance(g, topo, f, topo.faces[f].cells[1])
        for f = 1:n_faces(topo)
    ]
    masses=SVector{length(indices)}(map(j->molar_mass(material.components[j]), indices))
    return GasTransportKernel{
        typeof(model),
        NC,
        length(indices),
        typeof(indices),
        typeof(laws),
        typeof(masses),
    }(
        model,
        map(is_gas, material.components),
        indices,
        d,
        laws,
        masses,
    )
end

@inline function _cell_pore(ds, i, indices::NTuple{NG}) where {NG}
    p=SVector{NG}(ntuple(k->ds.p[indices[k], i], Val(NG)))
    K=SVector{NG}(ntuple(k->ds.K[indices[k], i], Val(NG)))
    return p, (K = K, t = ds.molecular_resistance[i], mob = ds.mob[i])
end

# The partial pressures an interior face reads. The lumped law is linear in
# them, so it takes the signed inventories: exact at every nonnegative mass,
# with slope at vacuum, where the guarded pressures have none and Newton's
# method fails; a negative inventory draws gas back across the face. The
# dusty-gas interface solve is posed on nonnegative pressures, so it keeps the
# guarded ones, as do the boundary solves.
@inline function _face_pressures(
    ::LumpedTransfer,
    kernel::GasTransportKernel{M,NC,NG},
    ds,
    i,
    guarded,
) where {M,NC,NG}
    return SVector{NG}(
        ntuple(
            k->ds.ξ[kernel.indices[k], i]/ds.φ[i]*R_GAS*ds.T_trial[i]/kernel.masses[k],
            Val(NG),
        ),
    )
end
@inline _face_pressures(::GasTransport, kernel, ds, i, guarded) = guarded

# The kernel's species mass fluxes through the interior face `f` between the
# live cells `o` and `nb` [kg/(m² s)], with the two floored center-to-face
# distances.
@inline function _face_mass_flux(kernel::GasTransportKernel, ds, topo, f, o, nb)
    d_o=smooth_floor(center_to_face_distance(ds.geometry, topo, f, o), kernel.d_floor[f])
    d_n=smooth_floor(center_to_face_distance(ds.geometry, topo, f, nb), kernel.d_floor[f])
    guarded_l, cl=_cell_pore(ds, o, kernel.indices)
    guarded_r, cr=_cell_pore(ds, nb, kernel.indices)
    pl=_face_pressures(kernel.model, kernel, ds, o, guarded_l)
    pr=_face_pressures(kernel.model, kernel, ds, nb, guarded_r)
    flux=interface_flux(pl, pr, ds.T_trial[o], ds.T_trial[nb], d_o, d_n, cl, cr)
    return flux.n .* kernel.masses, d_o, d_n
end

# Below this fraction of the larger cell temperature the enthalpy secant takes
# its limit, the heat capacity at the mean temperature. An enthalpy carries a
# rounding error of about ε c T, so above it the difference quotient keeps
# about four digits. Below it the limit misses the secant by at most c′ ΔT,
# with c′ the heat-capacity slope, and the face coupling by about
# |F| c′ ΔT²: under the flux's own rounding ε |F| c T for slopes up to about
# 1e5 J/(kg K²), four orders above any physical law.
const SECANT_RESOLUTION = 1e-11

# The conductivity `D` the face adds to its flux, and the face conductivity
# `k_f` of conduction. Per unit area the face moves heat into the owner's
# temperature row at (k_f + D − a_o)(T_n − T_o)/(d_o + d_n) and into the
# neighbour's at (k_f + D − a_n)(T_o − T_n)/(d_o + d_n), where
# a_i = (d_o + d_n) Σⱼ Fⱼ (hⱼ(T_f) − hⱼ(T_i))/(T_n − T_o) is the heat-capacity
# flux leaving cell i times its center-to-face distance; `D` is the least
# conductivity that keeps both coefficients nonnegative. `h` holds the gas
# enthalpies at the face temperature `T_f`, `z` the zero of the flux's type.
@inline function _added_conductivity(
    kernel::GasTransportKernel{M,NC,NG},
    ds,
    F,
    h,
    T_f,
    o,
    nb,
    d_o,
    d_n,
    z,
) where {M,NC,NG}
    T_o, T_n=ds.T[o], ds.T[nb]
    Δ=T_n-T_o
    d=d_o+d_n
    k_f=guarded_harmonic_mean(ds.k[o], ds.k[nb], d_o, d_n)
    a_o, a_n=z, z
    if abs(primal(Δ))>SECANT_RESOLUTION*max(abs(primal(T_o)), abs(primal(T_n)))
        for k = 1:NG
            j=kernel.indices[k]
            a_o+=F[k]*(h[j]-ds.h[j, o])
            a_n+=F[k]*(h[j]-ds.h[j, nb])
        end
        a_o*=d/Δ
        a_n*=d/Δ
    else
        c_o=map(c->c+z, component_enthalpy_slopes(kernel.laws, (T_o+T_f)/2))
        c_n=map(c->c+z, component_enthalpy_slopes(kernel.laws, (T_n+T_f)/2))
        for k = 1:NG
            j=kernel.indices[k]
            a_o+=F[k]*c_o[j]
            a_n-=F[k]*c_n[j]
        end
        a_o*=d_o
        a_n*=d_n
    end
    excess=(primal(a_o)>=primal(a_n) ? a_o : a_n)-k_f
    return (primal(excess)>0 ? excess : zero(excess)), k_f
end

function contribute!(
    c::Contributions{NC},
    kernel::GasTransportKernel{M,NC,NG},
    ds::DerivedState{NC},
    topo::MeshTopology{1},
    material,
    t,
) where {M,NC,NG}
    @inbounds for f = 1:n_faces(topo)
        o, nb=topo.faces[f].cells
        (nb==0 || !(ds.active[o] && ds.active[nb])) && continue
        F, d_o, d_n=_face_mass_flux(kernel, ds, topo, f, o, nb)
        Ts=(d_n*ds.T[o]+d_o*ds.T[nb])/(d_o+d_n)
        h=component_enthalpies(kernel.laws, Ts, kernel.gas)
        q=zero(eltype(c.F_energy))
        for k = 1:NG
            j=kernel.indices[k]
            c.F_mass[j, f]+=F[k]
            q+=F[k]*h[j]
        end
        D, _=_added_conductivity(kernel, ds, F, h, Ts, o, nb, d_o, d_n, zero(q))
        primal(D)>0 && (q+=D*(ds.T[o]-ds.T[nb])/(d_o+d_n))
        c.F_energy[f]+=q
    end
    return c
end

"""
    face_peclet(kernel::GasTransportKernel, ds, topo, f) -> (; peclet, diffusion)

The cell Péclet number of the interior face `f` and the conductivity the face
adds to its gas flux, relative to the face conductivity. With the species mass
fluxes `Fⱼ`, the enthalpy slopes `cⱼ` at the face temperature, the
center-to-face distances `d_o`, `d_n` and the harmonic face conductivity `k_f`
of conduction,

    peclet    = (d_o + d_n) |Σⱼ Fⱼ cⱼ| / k_f
    diffusion = D / k_f

with `D` the conductivity of [`GasTransportKernel`](@ref). `diffusion` is zero
where conduction alone keeps the face's couplings nonnegative and the face is
second-order accurate; with constant heat capacities on equal half cells it is
`max(0, peclet/2 − 1)`. A positive value marks a face that is upwind and first
order. NaN on a boundary face or beside a dead cell.
"""
function face_peclet(
    kernel::GasTransportKernel{M,NC,NG},
    ds,
    topo::MeshTopology{1},
    f,
) where {M,NC,NG}
    o, nb=topo.faces[f].cells
    (nb==0 || !(ds.active[o] && ds.active[nb])) && return (peclet = NaN, diffusion = NaN)
    F, d_o, d_n=_face_mass_flux(kernel, ds, topo, f, o, nb)
    T_f=(d_n*ds.T[o]+d_o*ds.T[nb])/(d_o+d_n)
    h=component_enthalpies(kernel.laws, T_f, kernel.gas)
    D, k_f=_added_conductivity(kernel, ds, F, h, T_f, o, nb, d_o, d_n, zero(T_f))
    s=component_enthalpy_slopes(kernel.laws, T_f)
    C=zero(T_f)
    for k = 1:NG
        C+=F[k]*s[kernel.indices[k]]
    end
    return (peclet = (d_o+d_n)*abs(C)/k_f, diffusion = D/k_f)
end
