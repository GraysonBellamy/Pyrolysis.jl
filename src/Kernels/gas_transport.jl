"""
    GasTransportKernel

One conservative species-flux vector per interior face, relative to the
skeleton. Total flux is matched through the two half cells; its thermal
interface state also supplies transported enthalpy. Boundary exchange is
assembled by the coupled thermal/mass boundary kernel.
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
        d_o=smooth_floor(
            center_to_face_distance(ds.geometry, topo, f, o),
            kernel.d_floor[f],
        )
        d_n=smooth_floor(
            center_to_face_distance(ds.geometry, topo, f, nb),
            kernel.d_floor[f],
        )
        pl, cl=_cell_pore(ds, o, kernel.indices)
        pr, cr=_cell_pore(ds, nb, kernel.indices)
        flux=interface_flux(pl, pr, ds.T_trial[o], ds.T_trial[nb], d_o, d_n, cl, cr)
        Ts=(d_n*ds.T[o]+d_o*ds.T[nb])/(d_o+d_n)
        h=component_enthalpies(kernel.laws, Ts, kernel.gas)
        q=zero(eltype(c.F_energy))
        for k = 1:NG
            j=kernel.indices[k]
            F=flux.n[k]*kernel.masses[k]
            c.F_mass[j, f]+=F
            q+=F*h[j]
        end
        c.F_energy[f]+=q
    end
    return c
end
