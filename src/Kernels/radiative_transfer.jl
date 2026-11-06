# ============================================================================
# Radiative transfer in the solid: gray two-stream rows
# ============================================================================

"""
Diffusivity factor of the two-stream rows: the ratio of a hemispherical flux's
mean path to the slab's thickness. With 3/2 the flux of the continuous
two-stream equations in an optically thick medium is the Rosseland flux.
"""
const DIFFUSIVITY = 3 / 2

"""
Smooth floor on the exposed face's transmissivity: the thermal rows stay
nonsingular however transparent the column and reflective the faces.
"""
const INTERFACE_FLOOR = 1e-3

# The radiation bands: one gray band.
const RADIATION_BANDS = 1

# The thermal families of a band: q⁻ (toward the back) and q⁺ (toward the
# exposed face).
const Q_DOWN = 1
const Q_UP = 2

"""
	MirrorBack()

A back face that exchanges no radiation: it reflects everything that reaches it.
"""
struct MirrorBack end

"""
	OpenBack(condition)

A back face that is a gray interface at the emissivity of the back cell, open
to the diffuse incoming radiation of its condition's radiative members.
"""
struct OpenBack{C<:AbstractThermalBC}
    condition::C
end

"""
	RadiativeTransferKernel{T,B} <: Kernel

Radiative transfer through the solid, gray, with two hemispherical fluxes per
face as algebraic unknowns: `q⁻` toward the back and `q⁺` toward the exposed
face [W/m²]. Every live cell `i` absorbs and emits at its absorption
coefficient `α` and temperature: with `τᵢ = α Δᵢ`, `tᵢ = exp(−D τᵢ)`
(`D = 3/2`, [`DIFFUSIVITY`](@ref)) and `Bᵢ = σ Tᵢ⁴`, the rows of a cell with
upstream face `u` (toward the exposed face) and downstream face `d` are the
exact solution of the two-stream equations across it:

	0 = q⁻_d − [Bᵢ + (q⁻_u − Bᵢ) tᵢ]
	0 = q⁺_u − [Bᵢ + (q⁺_d − Bᵢ) tᵢ]

and the cell absorbs `A_u (q⁻_u − q⁺_u) − A_d (q⁻_d − q⁺_d)`.

The exposed face `e` (the upstream face of the first live cell `i₁`) is a gray
interface with the transmissivity `tₑ = ε(i₁)`, smoothly floored at
[`INTERFACE_FLOOR`](@ref), both ways:

	0 = q⁻_e − [tₑ Q + (1 − tₑ) q⁺_e]

with `Q` the diffuse radiation arriving from outside: `F σ T∞⁴` per
`RadiativeBC` member and the flux per `RadiativeFluxBC` member of the `:top`
condition. The back face `b` either reflects everything (`q⁺_b = q⁻_b`, a
[`MirrorBack`](@ref)) or is a gray interface at the back cell's emissivity
`ε_l` ([`OpenBack`](@ref)): `q⁺_b = ε_l Q_b + (1 − ε_l) q⁻_b`.

An isothermal medium facing radiation at its own temperature satisfies every
row with `q⁻ = q⁺ = B` at any optical thickness, and the cell sources sum to
the net flux through the exposed face minus the net flux through the back.

`up` and `down` hold the upstream and downstream face of every cell, `order`
the cells from the exposed face into the depth, `d_floor` the floors of their
center-to-face distances. The faces upstream of the first live cell carry no
rows.
"""
struct RadiativeTransferKernel{T<:AbstractThermalBC,B} <: Kernel
    top::T
    back::B
    up::Vector{Int}
    down::Vector{Int}
    order::Vector{Int}
    d_floor::Matrix{Float64}
end

# A condition's radiative members.
_radiates(::AbstractThermalBC) = false
_radiates(::Union{RadiativeBC,RadiativeFluxBC}) = true
_radiates(bc::CombinedThermalBC) = any(_radiates, bc.bcs)

# Members the radiation rows refuse: a `RadiativeBC` with its own emissivity,
# a `RadiativeFluxBC` with its own absorptivity.
_given_emissivity(::AbstractThermalBC) = false
_given_emissivity(bc::RadiativeBC) = bc.ε !== nothing
_given_emissivity(bc::CombinedThermalBC) = any(_given_emissivity, bc.bcs)
_given_absorptivity(::AbstractThermalBC) = false
_given_absorptivity(bc::RadiativeFluxBC) = bc.absorptivity !== nothing
_given_absorptivity(bc::CombinedThermalBC) = any(_given_absorptivity, bc.bcs)
_has_flux(::AbstractThermalBC) = false
_has_flux(::RadiativeFluxBC) = true
_has_flux(bc::CombinedThermalBC) = any(_has_flux, bc.bcs)

_sets_absorptivity(component) =
    hasproperty(component, :surface_absorptivity) &&
    component.surface_absorptivity !== nothing

"""
	incoming_radiation(condition, t)

The diffuse radiation arriving at a face from outside [W/m²]: `F σ T∞⁴` per
`RadiativeBC` member of `condition` and the flux per `RadiativeFluxBC` member,
zero for every other member.
"""
@inline incoming_radiation(::AbstractThermalBC, t) = 0.0
@inline incoming_radiation(bc::RadiativeBC, t) =
    evaluate_param(bc.F, t) * STEFAN_BOLTZMANN * evaluate_param(bc.T_inf, t)^4
@inline incoming_radiation(bc::RadiativeFluxBC, t) = evaluate_param(bc.flux, t)
@inline incoming_radiation(bc::CombinedThermalBC, t) = _incoming_sum(bc.bcs, t)
@inline _incoming_sum(::Tuple{}, t) = 0.0
@inline _incoming_sum(b::Tuple, t) =
    incoming_radiation(first(b), t) + _incoming_sum(Base.tail(b), t)

# The largest surroundings temperature of a condition's `RadiativeBC` members
# at time `t`, zero without one.
_surroundings_temperature(::AbstractThermalBC, t) = 0.0
_surroundings_temperature(bc::RadiativeBC, t) = Float64(evaluate_param(bc.T_inf, t))
_surroundings_temperature(bc::CombinedThermalBC, t) =
    maximum(b -> _surroundings_temperature(b, t), bc.bcs; init = 0.0)

"""
	radiative_transfer_rows(bc_set) -> Int

The radiation rows per face under `RADIATIVE_TRANSFER`: the two thermal
families of every band.
"""
radiative_transfer_rows(bc_set::BoundaryConditionSet) = 2 * RADIATION_BANDS

# The row of family `j` of band `b` in a face's radiation rows.
@inline _family_row(b, j) = 2 * (b - 1) + j

# The cells from the face `exposed` into the depth, each cell's upstream and
# downstream face, and the smooth floors of their center-to-face distances.
function _column_tables(topo::MeshTopology{1}, g::MeshGeometry, exposed::Int, tag::Symbol)
    n = n_cells(topo)
    up = zeros(Int, n)
    down = zeros(Int, n)
    order = Int[]
    f, c = exposed, topo.faces[exposed].cells[1]
    while c != 0
        push!(order, c)
        cf = topo.cells[c].faces
        other = cf[1] == f ? cf[2] : cf[1]
        up[c], down[c] = f, other
        o, nb = topo.faces[other].cells
        f, c = other, (o == c ? nb : o)
    end
    length(order) == n ||
        throw(ArgumentError("the cells do not form one column from :$tag"))
    d = zeros(2, n)
    for c = 1:n
        d[1, c] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, up[c], c)
        d[2, c] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, down[c], c)
    end
    return up, down, order, d
end

"""
	RadiativeTransferKernel(topology, geometry, bc_set, material)

The kernel for radiation entering and leaving through the face tagged `:top`,
the exposed face, with the face tagged `:bottom` the back. Refused: a `:top`
condition without a `RadiativeBC` or `RadiativeFluxBC` member (radiation
enters and leaves the solid through the exposed face), a radiative member on
any other tag than `:top` and `:bottom`, a `RadiativeBC` with its own `ε` (the
interface's transmissivity is the emissivity of its cell, both ways), and a
`RadiativeFluxBC` with its own `absorptivity` or on a material where a
component sets one (diffuse radiation enters at the emissivity).
"""
function RadiativeTransferKernel(
    topo::MeshTopology{1},
    g::MeshGeometry,
    bc_set::BoundaryConditionSet,
    material,
)
    thermal = bc_set.thermal
    haskey(thermal, :top) && _radiates(thermal[:top]) || throw(
        ArgumentError(
            "radiative transfer needs a RadiativeBC or a RadiativeFluxBC on :top: " *
            "radiation enters and leaves the solid through the exposed face",
        ),
    )
    for (tag, bc) in pairs(thermal)
        _radiates(bc) &&
            !(tag in (:top, :bottom)) &&
            throw(
                ArgumentError(
                    "radiative transfer exchanges radiation through :top and :bottom only; " *
                    "the condition on :$tag has a radiative member",
                ),
            )
        _given_emissivity(bc) && throw(
            ArgumentError(
                "under radiative transfer a face's emissivity is that of its cell: " *
                "give the RadiativeBC on :$tag no ε",
            ),
        )
        _given_absorptivity(bc) && throw(
            ArgumentError(
                "diffuse radiation enters at the emissivity of the cell: give the " *
                "RadiativeFluxBC on :$tag no absorptivity",
            ),
        )
    end
    if any(_has_flux, values(thermal))
        for c in material.components
            _sets_absorptivity(c) && throw(
                ArgumentError(
                    "diffuse radiation enters at the emissivity of the cell; " *
                    "component :$(c.name) sets its own absorptivity",
                ),
            )
        end
    end
    faces = boundary_faces(topo, :top)
    length(faces) == 1 || throw(
        ArgumentError(
            "radiative transfer needs one exposed face (:top has $(length(faces)))",
        ),
    )
    up, down, order, d = _column_tables(topo, g, only(faces), :top)
    back =
        haskey(thermal, :bottom) && _radiates(thermal[:bottom]) ?
        OpenBack(thermal[:bottom]) : MirrorBack()
    return RadiativeTransferKernel(thermal[:top], back, up, down, order, d)
end

# The position in `order` of the first live cell, zero when none is live.
@inline function _first_live(kernel::RadiativeTransferKernel, ds)
    @inbounds for s in eachindex(kernel.order)
        ds.active[kernel.order[s]] && return s
    end
    return 0
end

@inline function exposed_face(kernel::RadiativeTransferKernel, ds)
    s = _first_live(kernel, ds)
    return s == 0 ? 0 : kernel.up[kernel.order[s]]
end

# The transmissivity `t` of cell `c` (`exp(−D τ)`) and its emissive power `B`.
@inline function _cell_transfer(kernel::RadiativeTransferKernel, ds, topo, c)
    g = ds.geometry
    Δ =
        smooth_floor(
            center_to_face_distance(g, topo, kernel.up[c], c),
            kernel.d_floor[1, c],
        ) + smooth_floor(
            center_to_face_distance(g, topo, kernel.down[c], c),
            kernel.d_floor[2, c],
        )
    return exp(-DIFFUSIVITY * ds.α[c] * Δ), STEFAN_BOLTZMANN * ds.T_trial[c]^4
end

# The back face's reflectivity `ρ` and injection `J`: `q⁺_b = J + ρ q⁻_b`.
@inline _back_row(::MirrorBack, ε, t) = (one(ε), zero(ε))
@inline _back_row(back::OpenBack, ε, t) = (1 - ε, ε * incoming_radiation(back.condition, t))

function contribute!(
    c::Contributions,
    kernel::RadiativeTransferKernel,
    ds::DerivedState,
    topo::MeshTopology{1},
    material,
    t,
)
    s = _first_live(kernel, ds)
    s == 0 && return c
    g = ds.geometry
    q, r = ds.q, c.r_radiation
    order = kernel.order
    i₁ = order[s]
    e = kernel.up[i₁]
    l = last(order)
    b_face = kernel.down[l]
    tₑ = smooth_floor(ds.ε[i₁], INTERFACE_FLOOR)
    Q = incoming_radiation(kernel.top, t)
    ρ_b, J_b = _back_row(kernel.back, ds.ε[l], t)
    for b = 1:RADIATION_BANDS
        dn, up = _family_row(b, Q_DOWN), _family_row(b, Q_UP)
        r[dn, e] = q[dn, e] - (tₑ * Q + (1 - tₑ) * q[up, e])
        @inbounds for k = s:length(order)
            i = order[k]
            u, d = kernel.up[i], kernel.down[i]
            tr, B = _cell_transfer(kernel, ds, topo, i)
            r[dn, d] = q[dn, d] - (B + (q[dn, u] - B) * tr)
            r[up, u] = q[up, u] - (B + (q[up, d] - B) * tr)
            c.S_energy[i] +=
                (q[dn, u] - q[up, u]) * g.face_area[u] -
                (q[dn, d] - q[up, d]) * g.face_area[d]
        end
        r[up, b_face] = q[up, b_face] - (J_b + ρ_b * q[dn, b_face])
    end
    return c
end

"""
	consistent_radiation!(q, kernel, ds, topology, t) -> q

The face fluxes `q` (radiation rows × faces) that satisfy every radiation row
at the derived state `ds`, and zero on the faces upstream of the first live
cell. The rows are linear in the fluxes: with `q⁻_e = a` as a parameter, a
downward pass gives `q⁻_b`, the back row `q⁺_b`, an upward pass
`q⁺_e = A a + C`, the exposed row `a`; a downward and an upward pass write the
values. No iteration, no allocation.
"""
function consistent_radiation!(
    q::AbstractMatrix,
    kernel::RadiativeTransferKernel,
    ds::DerivedState,
    topo::MeshTopology{1},
    t,
)
    order = kernel.order
    n = length(order)
    s = _first_live(kernel, ds)
    for k = 1:(s==0 ? n : s-1), j in axes(q, 1)
        q[j, kernel.up[order[k]]] = 0
    end
    s == 0 && return q
    i₁ = order[s]
    e = kernel.up[i₁]
    l = order[n]
    b_face = kernel.down[l]
    tₑ = smooth_floor(ds.ε[i₁], INTERFACE_FLOOR)
    Q = incoming_radiation(kernel.top, t)
    ρ_b, J_b = _back_row(kernel.back, ds.ε[l], t)
    for b = 1:RADIATION_BANDS
        dn, up = _family_row(b, Q_DOWN), _family_row(b, Q_UP)
        A, C = one(tₑ), zero(tₑ)
        for k = s:n
            tr, B = _cell_transfer(kernel, ds, topo, order[k])
            A, C = A * tr, B + (C - B) * tr
        end
        A, C = ρ_b * A, J_b + ρ_b * C
        for k = n:-1:s
            tr, B = _cell_transfer(kernel, ds, topo, order[k])
            A, C = A * tr, B + (C - B) * tr
        end
        q[dn, e] = (tₑ * Q + (1 - tₑ) * C) / (1 - (1 - tₑ) * A)
        for k = s:n
            i = order[k]
            tr, B = _cell_transfer(kernel, ds, topo, i)
            q[dn, kernel.down[i]] = B + (q[dn, kernel.up[i]] - B) * tr
        end
        q[up, b_face] = J_b + ρ_b * q[dn, b_face]
        for k = n:-1:s
            i = order[k]
            tr, B = _cell_transfer(kernel, ds, topo, i)
            q[up, kernel.up[i]] = B + (q[up, kernel.down[i]] - B) * tr
        end
    end
    return q
end
