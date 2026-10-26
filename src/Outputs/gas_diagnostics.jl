"""
    GasDiagnostics

Saved pore states and exact assembled boundary exchange. Gas arrays follow
`species`, in material order. Cell arrays have dimensions species × cells ×
times (porosity and total pressure: cells × times). Dead cells and absent live
boundaries report NaN. Vacuum fractions are zero. Boundary fluxes are outward,
area-weighted over the current live faces of each tag; Darcy velocity is [m/s],
species mass flux [kg/(m² s)], partial pressure [Pa], and energy flux [W/m²].
`heat` is positive inward at the reservoir control surface; `enthalpy` and
`energy` are positive outward, with energy = enthalpy − heat.
"""
struct GasDiagnostics
    species::Vector{Symbol}
    porosity::Matrix{Float64}
    pressure::Matrix{Float64}
    partial_pressure::Array{Float64,3}
    intrinsic_density::Array{Float64,3}
    mole_fraction::Array{Float64,3}
    mass_fraction::Array{Float64,3}
    boundary::Dict{Symbol,NamedTuple}
    inventories::Matrix{Float64}
    production::Matrix{Float64}
    escape::Matrix{Float64}
    ingress::Matrix{Float64}
end
function GasDiagnostics(material, n, nt, tags)
    species=[c.name for c in material.components if is_gas(c)]
    ng=length(species)
    boundary=Dict{Symbol,NamedTuple}(
        tag=>(;
            fluxes = fill(NaN, ng, nt),
            partial_pressure = fill(NaN, ng, nt),
            darcy_velocity = fill(NaN, nt),
            enthalpy = fill(NaN, nt),
            heat = fill(NaN, nt),
            energy = fill(NaN, nt),
            residual = fill(NaN, nt),
            iterations = zeros(Int, nt),
        ) for tag in tags
    )
    GasDiagnostics(
        species,
        fill(NaN, n, nt),
        fill(NaN, n, nt),
        (fill(NaN, ng, n, nt) for _ = 1:4)...,
        boundary,
        (zeros(ng, nt) for _ = 1:4)...,
    )
end
function _gas_diagnostics!(out::GasDiagnostics, k, def, ds, c)
    mat=def.setup.material
    topo=def.setup.topology
    indices=gas_component_indices(mat.components)
    for i = 1:n_cells(topo)
        ds.active[i] || continue
        p=pore_state(Tuple(ds.ξ⁺[:, i]), ds.T_trial[i], mat)
        out.porosity[i, k]=p.φ
        out.pressure[i, k]=p.P
        for (g, j) in enumerate(indices)
            out.partial_pressure[g, i, k]=p.p[j]
            out.intrinsic_density[g, i, k]=p.ρ[j]
            out.mole_fraction[g, i, k]=p.x[j]
            out.mass_fraction[g, i, k]=p.Y[j]
            out.inventories[g, k]+=ds.m[j, i]
            out.production[g, k]+=c.S_mass[j, i]
        end
    end
    for (tag, b) in out.boundary
        faces=_live_faces(topo, ds.active, tag)
        isempty(faces) && continue
        A=sum(ds.geometry.face_area[f] for f in faces)
        wt(f) = A>0 ? ds.geometry.face_area[f]/A : 1/length(faces)
        for (g, j) in enumerate(indices)
            b.fluxes[g, k]=sum(
                _outward_sign(topo, ds.active, f)*c.F_mass[j, f]*wt(f) for f in faces
            )
            b.partial_pressure[g, k]=sum(c.p_surface[j, f]*wt(f) for f in faces)
            for f in faces
                N=_outward_sign(topo, ds.active, f)*c.F_mass[j, f]*ds.geometry.face_area[f]
                out.escape[g, k]+=max(N, 0.0)
                out.ingress[g, k]+=max(-N, 0.0)
            end
        end
        b.darcy_velocity[k]=sum(c.darcy_velocity[f]*wt(f) for f in faces)
        b.enthalpy[k]=sum(c.boundary_enthalpy[f]*wt(f) for f in faces)
        b.heat[k]=sum(c.boundary_heat[f]*wt(f) for f in faces)
        b.energy[k]=sum(
            _outward_sign(topo, ds.active, f)*c.F_energy[f]*wt(f) for f in faces
        )
        b.residual[k]=maximum(c.surface_residual[f] for f in faces)
        b.iterations[k]=maximum(c.surface_iterations[f] for f in faces)
    end
    return out
end
