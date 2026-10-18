# ============================================================================
# Gas transport in flux form
# ============================================================================

# Guard of the smooth face mean of the Darcy mobility κ/μ [m²/(Pa·s)]: far
# below the mobility of any permeable solid (κ = 1e-18 m², μ = 5e-5 Pa·s gives
# 2e-14).
const MOBILITY_GUARD = 1e-30

"""
	GasTransportKernel{M,NC,P<:Tuple,L} <: Kernel

The mass flux of every gas component through every face, per unit area along
the face normal, and the sensible enthalpy it carries.

Interior faces: the volume-fraction-gradient flux
`J = −(λ_f/T̄) (ξ_nb T_nb − ξ_o T_o) / d` (`Physics.gas_flux`) at the trial
temperatures, with `λ_f` the guarded distance-weighted harmonic mean of the
cell gas-transfer coefficients and `d` the center distance; under Darcy–Fick
transport (`M = DarcyFickTransport`) also the advective flux `v ξ_up / φ_f`
with the superficial velocity `v = −m_f (P_nb − P_o) / d` from the guarded face
mean `m_f` of the cell mobilities, the face porosity `φ_f` the mean of the
smoothly floored cell porosities, and `ξ_up` upwind of `v`. The energy flux
of the face gains `Σⱼ Jⱼ hⱼ(T_f)` with `T_f` the distance-weighted face
temperature, the same mass flux the mass rows carry.

Boundary faces, one [`BoundaryPatch`](@ref) of mass conditions per tag: the
film flux of the condition (`BoundaryConditions.mass_bc_flux`, oriented
outward) and, under Darcy–Fick transport, the outflow the pressure condition
drives at the outward Darcy velocity of [`outward_darcy_velocity`](@ref)
against an ambient that carries none of the tracked gases. The enthalpy those
fluxes carry leaves at the face temperature, which the thermal boundary
kernel adds once the surface balance has fixed it.

`gas` flags the gas components; `laws` gives the sensible enthalpies.
"""
struct GasTransportKernel{M<:GasTransport,NC,P<:Tuple,L<:CaloricLaws} <: Kernel
    model::M
    gas::NTuple{NC,Bool}
    patches::P
    d_floor::Vector{Float64}
    laws::L
end

"""
	GasTransportKernel(model, topology, geometry, bc_set, material, laws)

The gas-transport kernel with one mass-condition patch per boundary tag of
the mesh and distance floors from the reference geometry.
"""
function GasTransportKernel(
    model::GasTransport,
    topo::MeshTopology,
    g::MeshGeometry,
    bc_set::BoundaryConditionSet,
    material::Material,
    laws::CaloricLaws,
)
    tags = boundary_tags(topo)
    for tag in keys(bc_set.mass)
        tag in tags || throw(
            ArgumentError(
                "the mass conditions name the tag :$tag, which the mesh does not have " *
                "(mesh tags: $tags)",
            ),
        )
    end
    patches = Tuple(
        BoundaryPatch(tag, boundary_faces(topo, tag), bc_set.mass[tag]) for tag in tags
    )
    d = zeros(n_faces(topo))
    for f = 1:n_faces(topo)
        o = topo.faces[f].cells[1]
        d[f] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, f, o)
    end
    return GasTransportKernel(model, map(is_gas, material.components), patches, d, laws)
end

function contribute!(
    c::Contributions{NC},
    kernel::GasTransportKernel{M,NC},
    ds::DerivedState{NC},
    topo::MeshTopology{1},
    material,
    t,
) where {M,NC}
    g = ds.geometry
    gas = kernel.gas
    @inbounds for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        nb == 0 && continue
        floor = kernel.d_floor[f]
        d_o = smooth_floor(center_to_face_distance(g, topo, f, o), floor)
        d_n = smooth_floor(center_to_face_distance(g, topo, f, nb), floor)
        d = d_o + d_n
        λ_f = guarded_harmonic_mean(ds.λ[o], ds.λ[nb], d_o, d_n)
        T_o, T_n = ds.T_trial[o], ds.T_trial[nb]
        v, φ_f = _darcy_face(kernel.model, ds, o, nb, d_o, d_n)
        h = component_enthalpies(kernel.laws, (d_n * ds.T[o] + d_o * ds.T[nb]) / d)
        q = zero(eltype(c.F_energy))
        for j = 1:NC
            gas[j] || continue
            ξ_o, ξ_n = ds.ξ[j, o], ds.ξ[j, nb]
            F = gas_flux(λ_f, ξ_o, ξ_n, T_o, T_n, d) + _advective_flux(v, φ_f, ξ_o, ξ_n)
            c.F_mass[j, f] += F
            q += F * h[j]
        end
        c.F_energy[f] += q
    end
    _mass_patches!(c, kernel.patches, kernel, ds, topo, material, t)
    return c
end

# Superficial Darcy velocity along the normal of an interior face, and the
# face porosity; none under Fickian transport.
@inline _darcy_face(::FickTransport, ds, o, nb, d_o, d_n) = (nothing, nothing)
@inline function _darcy_face(::DarcyFickTransport, ds, o, nb, d_o, d_n)
    m_f = guarded_harmonic_mean(ds.mob[o], ds.mob[nb], d_o, d_n, MOBILITY_GUARD)
    v = -m_f * (ds.P[nb] - ds.P[o]) / (d_o + d_n)
    return v, (ds.φ[o] + ds.φ[nb]) / 2
end

# First-order upwind advection of the intrinsic gas density ξ/φ.
@inline _advective_flux(::Nothing, ::Nothing, ξ_o, ξ_n) = zero(ξ_o)
@inline _advective_flux(v, φ_f, ξ_o, ξ_n) = v * (v >= 0 ? ξ_o : ξ_n) / φ_f

# Peel the patch tuple so each patch's condition is a concrete type.
@inline _mass_patches!(c, ::Tuple{}, kernel, ds, topo, material, t) = nothing
@inline function _mass_patches!(c, patches::Tuple, kernel, ds, topo, material, t)
    p = first(patches)
    @inbounds for f in p.faces
        _boundary_mass_flux!(c, p.bc, f, kernel, ds, topo, material, t)
    end
    return _mass_patches!(c, Base.tail(patches), kernel, ds, topo, material, t)
end

# Outward gas mass fluxes through boundary face `f` under the mass condition `bc`.
@inline function _boundary_mass_flux!(
    c::Contributions{NC},
    bc::AbstractMassBC,
    f,
    kernel::GasTransportKernel{M,NC},
    ds,
    topo,
    material,
    t,
) where {M,NC}
    i = topo.faces[f].cells[1]
    d = smooth_floor(center_to_face_distance(ds.geometry, topo, f, i), kernel.d_floor[f])
    ξ⁺ = ntuple(j -> ds.ξ⁺[j, i], Val(NC))
    ρ = mixture_density(ξ⁺, ds.T_trial[i], material)
    λρ = ds.λ[i] * ρ
    v_over_φ = _vent_rate(kernel.model, bc, ds, i, d)
    @inbounds for j = 1:NC
        kernel.gas[j] || continue
        ξ = ds.ξ[j, i]
        c.F_mass[j, f] += v_over_φ * ξ - mass_bc_flux(bc, ξ, λρ, ρ, d, t, j)
    end
    return nothing
end

# Outflow per unit intrinsic gas density, v_out/φ, that the pressure condition
# drives through the face; inflow carries none of the tracked gases.
@inline _vent_rate(::FickTransport, bc, ds, i, d) = zero(eltype(ds.T))
@inline function _vent_rate(::DarcyFickTransport, bc, ds, i, d)
    v = outward_darcy_velocity(bc, ds.P[i], ds.mob[i], d)
    return v > 0 ? v / ds.φ[i] : zero(v)
end

"""
	outward_darcy_velocity(bc, P, mobility, d) -> v

Superficial Darcy velocity [m/s] along the outward normal of a boundary face
under the mass condition `bc`, from the pressure `P` [Pa] and the mobility
`κ/μ` [m²/(Pa·s)] of the boundary cell, its center lying `d` [m] from the face:

- `DirichletPressureBC`: the one-sided gradient, `v = m (P − P_bc) / d`;
- `NeumannPressureBC`: `∂P∂n` is the pressure gradient along the outward
  normal, `v = −m ∂P∂n`;
- `ConvectivePressureBC`: the film law `v = h_P (P − P_ext)` where the cell is
  permeable, zero where it is not;
- every other condition drives none; a `CombinedMassBC` sums its members.
"""
@inline outward_darcy_velocity(bc::DirichletPressureBC, P, mob, d) = mob * (P - bc.P) / d
@inline outward_darcy_velocity(bc::NeumannPressureBC, P, mob, d) = -mob * bc.∂P∂n + zero(P)
@inline outward_darcy_velocity(bc::ConvectivePressureBC, P, mob, d) =
    mob > 0 ? bc.h_P * (P - bc.P_ext) : zero(P * mob)
@inline outward_darcy_velocity(::AbstractMassBC, P, mob, d) = zero(P * mob / d)
@inline outward_darcy_velocity(bc::CombinedMassBC, P, mob, d) =
    _sum_velocity(bc.bcs, P, mob, d)
@inline _sum_velocity(::Tuple{}, P, mob, d) = zero(P * mob / d)
@inline _sum_velocity(bcs::Tuple, P, mob, d) =
    outward_darcy_velocity(first(bcs), P, mob, d) + _sum_velocity(Base.tail(bcs), P, mob, d)
