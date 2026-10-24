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
outermost live cell. On a mesh of more than one dimension, where cells leave
by element death, a face between a live cell and a dead one is a boundary face
of the live cell, its normal turned out of it, and carries the condition of
the exposed tag (`exposed`, the condition itself, or `nothing` on a line).

A prescribed temperature or heat flux, or an adiabatic face, gives the flux in
closed form. Every other condition goes through the surface energy balance
`q_bc(T_s) = k (T_s − T_cell) / d` of the face, solved by the Newton iteration
of `BoundaryConditions.solve_boundary_temperature!` started from the cell's
trial temperature, with the condition evaluated at the trial temperature of
the iterate ([`trial_temperature`](@ref)), so the balance has a root at any
trial state and its solve never clamps; the outward gas mass flux of the face
feeds the blowing models of its convective members, and the cell's emissivity
and absorptivity its radiative members. `states` holds the per-face
diagnostics of the balance.

With gas transport (`laws` the caloric laws, `nothing` without), the face's
energy flux also carries the sensible enthalpy `Σⱼ Jⱼ hⱼ(T_s)` of its outward
mass fluxes at the face temperature: gas leaves at the surface temperature,
and the surface balance itself keeps no gas-enthalpy term.

On a mesh of more than one dimension the conductivity is a tensor `K = k A`
(`tensor` holds the cells' dimensionless `A`, shared with the
[`ConductionKernel`](@ref)), and the flux conducted from the face into its
cell is `(n·K·n)(T_s − T_cell)/d + t·∇T_cell`, `t = K n − (n·K·n)(x_f − x_cell)/d`
tangential to the face and `∇T_cell` the conduction kernel's least-squares
gradient, so a linear field's flux is exact; the surface balance solves for
`T_s` with the cross part folded into the cell temperature. On a line the
cross part vanishes and `n·K·n` is the cell's conductivity.
"""
struct ThermalBoundaryKernel{P<:Tuple,R,L,A,X} <: Kernel
    patches::P
    d_floor::Vector{Float64}
    states::Vector{BoundaryState}
    radiation_model::R
    laws::L
    tensor::A
    exposed::X
end

"""
	ThermalBoundaryKernel(topology, geometry, bc_set; radiation_model = NO_RADIATION,
	                      laws = nothing, tensor = nothing, exposed = :top)

One patch per boundary tag of the mesh, holding the thermal condition the
boundary-condition set assigns to that tag; `laws` are the caloric laws when
gas crosses the boundary faces; `tensor` the cells' dimensionless
conductivity tensors on a mesh of more than one dimension, those of the
[`ConductionKernel`](@ref), and `exposed` the tag whose condition the faces
element death exposes take there.
"""
function ThermalBoundaryKernel(
    topo::MeshTopology{D},
    g::MeshGeometry,
    bc_set::BoundaryConditionSet;
    radiation_model::RadiationModel = NO_RADIATION,
    laws::Union{Nothing,CaloricLaws} = nothing,
    tensor = nothing,
    exposed::Symbol = :top,
) where {D}
    (D == 1) == (tensor === nothing) || throw(
        ArgumentError(
            "the cells' tensors are given on, and only on, more than one dimension",
        ),
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
    condition =
        D == 1 || !(exposed in tags) ? nothing :
        only(p.bc for p in patches if p.tag === exposed)
    return ThermalBoundaryKernel(
        patches,
        d,
        states,
        radiation_model,
        laws,
        tensor,
        condition,
    )
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
    _exposed!(c, kernel.exposed, kernel, ds, topo, material, t)
    return c
end

# The faces between a live cell and a dead one, under the exposed condition.
@inline _exposed!(c, ::Nothing, kernel, ds, topo, material, t) = nothing
function _exposed!(c, bc, kernel, ds, topo, material, t)
    @inbounds for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        (nb == 0 || ds.active[o] == ds.active[nb]) && continue
        q, T_s = boundary_energy_flux(bc, f, kernel, c, ds, topo, material, t)
        σ = ds.active[o] ? 1 : -1
        c.F_energy[f] += σ * q + _carried_enthalpy(kernel.laws, c, f, T_s)
        c.T_face[f] = T_s
    end
    return nothing
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

# The boundary cell of face `f`, its live side, and its smoothly floored
# center-to-face distance.
@inline function _boundary_cell(kernel, ds, topo, f)
    o, nb = topo.faces[f].cells
    i = nb == 0 || ds.active[o] ? o : nb
    d = smooth_floor(center_to_face_distance(ds.geometry, topo, f, i), kernel.d_floor[f])
    return i, d
end

# The sign of face `f`'s normal out of its cell `i`.
@inline _outward(topo, f, i) = topo.faces[f].cells[1] == i ? 1 : -1

# The normal conductivity `n·K·n` of boundary face `f`'s cell `i` and the cross
# part `t·∇T_i` of the flux conducted from the face into the cell: the cell's
# conductivity and nothing on a line.
@inline _boundary_conduction(::Nothing, c, ds, topo, f, i, d) =
    (ds.k[i], zero(eltype(ds.T)))
@inline function _boundary_conduction(
    A::AbstractVector,
    c,
    ds,
    topo::MeshTopology{D},
    f,
    i,
    d,
) where {D}
    g = ds.geometry
    n = _outward(topo, f, i) * g.face_normal[f]
    K = ds.k[i] * A[i]
    k_n = normal_conductivity(K, n)
    t = K * n - k_n * (g.face_center[f] - g.cell_center[i]) / d
    return k_n, dot(t, cell_gradient(c, Val(D), i))
end

"""
	boundary_condition(kernel::ThermalBoundaryKernel, tag) -> AbstractThermalBC

The thermal condition the kernel applies on the faces of boundary tag `tag`.
"""
function boundary_condition(kernel::ThermalBoundaryKernel, tag::Symbol)
    for p in kernel.patches
        p.tag === tag && return p.bc
    end
    throw(ArgumentError("the thermal boundary kernel has no patch for the tag :$tag"))
end

"""
	conducted_flux(kernel::ThermalBoundaryKernel, c, ds, topology, f) -> Real

The heat flux [W/m²] conducted from boundary face `f` into its cell,
`k (T_s − T_cell) / d`, at the face temperature of the contributions `c`.
"""
function conducted_flux(kernel::ThermalBoundaryKernel, c::Contributions, ds, topo, f)
    i, d = _boundary_cell(kernel, ds, topo, f)
    k, χ = _boundary_conduction(kernel.tensor, c, ds, topo, f, i, d)
    return k * (c.T_face[f] - ds.T[i]) / d + χ
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
    k, χ = _boundary_conduction(kernel.tensor, c, ds, topo, f, i, d)
    q_in = evaluate_param(bc.flux, t) + zero(eltype(ds.T))
    return -q_in, ds.T[i] + (q_in - χ) * d / k
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
    k, χ = _boundary_conduction(kernel.tensor, c, ds, topo, f, i, d)
    T_bc = evaluate_param(bc.T, t) + zero(eltype(ds.T))
    return k * (ds.T[i] - T_bc) / d - χ, T_bc
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
    k, χ = _boundary_conduction(kernel.tensor, c, ds, topo, f, i, d)
    σ = _outward(topo, f, i)
    J_out = ntuple(j -> σ * c.F_mass[j, f], Val(NC))
    r = solve_boundary_temperature!(
        kernel.states[f],
        ds.T_trial[i] - χ * d / k,
        k,
        d,
        bc,
        material,
        J_out,
        t,
        ds.ε[i];
        a_eff = ds.a[i],
        radiation_model = kernel.radiation_model,
        trial = trial_temperature,
    )
    return -r.q_conducted, r.T_s
end
