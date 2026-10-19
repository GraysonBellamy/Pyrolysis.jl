# ============================================================================
# Thermal boundary conditions
# ============================================================================

"""
	BoundaryPatch{BC}

The boundary faces carrying one tag and the condition applied on them.
"""
struct BoundaryPatch{BC}
    tag::Symbol
    faces::Vector{Int}
    bc::BC
end

"""
	ThermalBoundaryKernel{P<:Tuple,R,L} <: Kernel

The energy flux through every boundary face, one [`BoundaryPatch`](@ref) per
boundary tag, written along the outward face normal, and the face temperature.
A patch applies on the face that bounds the live cells on its side
(`Meshes.live_face`): after cells beyond it are deactivated, the face of the
outermost live cell.

A prescribed temperature or heat flux, or an adiabatic face, gives the flux in
closed form. Every other condition goes through the surface energy balance
`q_bc(T_s) = k (T_s − T_cell) / d` of the face, solved by the Newton iteration
of `BoundaryConditions.solve_boundary_temperature!` started from the clamped
cell temperature; the outward gas mass flux of the face feeds the blowing
models of its convective members, and the cell's emissivity and absorptivity
its radiative members. `states` holds the per-face diagnostics of the balance.

With gas transport (`laws` the caloric laws, `nothing` without), the face's
energy flux also carries the sensible enthalpy `Σⱼ Jⱼ hⱼ(T_s)` of its outward
mass fluxes at the face temperature: gas leaves at the surface temperature,
and the surface balance itself keeps no gas-enthalpy term.
"""
struct ThermalBoundaryKernel{P<:Tuple,R,L} <: Kernel
    patches::P
    d_floor::Vector{Float64}
    states::Vector{BoundaryState}
    radiation_model::R
    laws::L
end

"""
	ThermalBoundaryKernel(topology, geometry, bc_set; radiation_model = NO_RADIATION,
	                      laws = nothing)

One patch per boundary tag of the mesh, holding the thermal condition the
boundary-condition set assigns to that tag; `laws` are the caloric laws when
gas crosses the boundary faces.
"""
function ThermalBoundaryKernel(
    topo::MeshTopology,
    g::MeshGeometry,
    bc_set::BoundaryConditionSet;
    radiation_model::RadiationModel = NO_RADIATION,
    laws::Union{Nothing,CaloricLaws} = nothing,
)
    tags = boundary_tags(topo)
    for tag in keys(bc_set.thermal)
        tag in tags || throw(
            ArgumentError(
                "the boundary conditions name the tag :$tag, which the mesh does not have " *
                "(mesh tags: $tags)",
            ),
        )
    end
    patches = Tuple(
        BoundaryPatch(tag, boundary_faces(topo, tag), bc_set.thermal[tag]) for tag in tags
    )
    d = zeros(n_faces(topo))
    for f = 1:n_faces(topo)
        o = topo.faces[f].cells[1]
        d[f] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, f, o)
    end
    states = [BoundaryState() for _ = 1:n_faces(topo)]
    return ThermalBoundaryKernel(patches, d, states, radiation_model, laws)
end

function contribute!(
    c::Contributions,
    kernel::ThermalBoundaryKernel,
    ds::DerivedState,
    topo::MeshTopology,
    material,
    t,
)
    _patches!(c, kernel.patches, kernel, ds, topo, material, t)
    return c
end

# Peel the patch tuple so each patch's condition is a concrete type.
@inline _patches!(c, ::Tuple{}, kernel, ds, topo, material, t) = nothing
@inline function _patches!(c, patches::Tuple, kernel, ds, topo, material, t)
    p = first(patches)
    @inbounds for f0 in p.faces
        f = live_face(topo, ds.active, f0)
        f == 0 && continue
        q, T_s = boundary_energy_flux(p.bc, f, kernel, c, ds, topo, material, t)
        c.F_energy[f] += q + _carried_enthalpy(kernel.laws, c, f, T_s)
        c.T_face[f] = T_s
    end
    return _patches!(c, Base.tail(patches), kernel, ds, topo, material, t)
end

# Sensible enthalpy carried by the outward mass fluxes of face `f` at `T_s`.
@inline _carried_enthalpy(::Nothing, c, f, T_s) = zero(T_s)
@inline function _carried_enthalpy(
    laws::CaloricLaws,
    c::Contributions{NC},
    f,
    T_s,
) where {NC}
    h = component_enthalpies(laws, T_s)
    q = zero(promote_type(eltype(c.F_mass), typeof(T_s)))
    @inbounds for j = 1:NC
        q += c.F_mass[j, f] * h[j]
    end
    return q
end

# The boundary cell of face `f` and its smoothly floored center-to-face distance.
@inline function _boundary_cell(kernel, ds, topo, f)
    i = topo.faces[f].cells[1]
    d = smooth_floor(center_to_face_distance(ds.geometry, topo, f, i), kernel.d_floor[f])
    return i, d
end

"""
	boundary_energy_flux(bc, f, kernel, c, ds, topology, material, t) -> (q, T_s)

Outward energy flux per unit area [W/m²] through boundary face `f` under the
thermal condition `bc`, and the face temperature [K].
"""
@inline function boundary_energy_flux(
    ::AdiabaticBC,
    f,
    kernel,
    c::Contributions,
    ds,
    topo,
    material,
    t,
)
    i, _ = _boundary_cell(kernel, ds, topo, f)
    return zero(eltype(ds.T)), ds.T[i]
end

@inline function boundary_energy_flux(
    bc::HeatFluxBC,
    f,
    kernel,
    c::Contributions,
    ds,
    topo,
    material,
    t,
)
    i, d = _boundary_cell(kernel, ds, topo, f)
    q_in = evaluate_param(bc.flux, t) + zero(eltype(ds.T))
    return -q_in, ds.T[i] + q_in * d / ds.k[i]
end

@inline function boundary_energy_flux(
    bc::DirichletThermalBC,
    f,
    kernel,
    c::Contributions,
    ds,
    topo,
    material,
    t,
)
    i, d = _boundary_cell(kernel, ds, topo, f)
    T_bc = evaluate_param(bc.T, t) + zero(eltype(ds.T))
    return ds.k[i] * (ds.T[i] - T_bc) / d, T_bc
end

@inline function boundary_energy_flux(
    bc::AbstractThermalBC,
    f,
    kernel,
    c::Contributions{NC},
    ds,
    topo,
    material,
    t,
) where {NC}
    i, d = _boundary_cell(kernel, ds, topo, f)
    J_out = ntuple(j -> c.F_mass[j, f], Val(NC))
    r = solve_boundary_temperature!(
        kernel.states[f],
        ds.T[i],
        ds.k[i],
        d,
        bc,
        material,
        J_out,
        t,
        ds.ε[i];
        a_eff = ds.a[i],
        radiation_model = kernel.radiation_model,
    )
    return -r.q_conducted, r.T_s
end
