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
# exposed face). The beam families follow them.
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
	BeamSource{F,A}

The beams of one `RadiativeFluxBC` member with a [`Collimated`](@ref) or
[`Cone`](@ref) distribution: its flux and absorptivity, the direction cosine
`μ[h]` of each beam and the share `w[h]` of the flux it carries (`Σ w = 1`),
and `row`, the radiation row of its first beam in a band.
"""
struct BeamSource{F,A}
    flux::F
    absorptivity::A
    μ::Vector{Float64}
    w::Vector{Float64}
    row::Int
end

"""
	RadiativeTransferKernel{T,B,S} <: Kernel

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
`RadiativeBC` member and the flux per diffuse `RadiativeFluxBC` member of the
`:top` condition. The back face `b` either reflects everything
([`MirrorBack`](@ref)) or is a gray interface at the back cell's emissivity
`ε_l` ([`OpenBack`](@ref)), open to the diffuse radiation `Q_b` of its
condition.

A `RadiativeFluxBC` member of the `:top` condition with a [`Collimated`](@ref)
or [`Cone`](@ref) distribution adds one beam family per direction, its flux
along the normal `I_h` an algebraic unknown of every face
([`BeamSource`](@ref)). A beam enters with its share `w_h` of the member's
flux `q` at the member's absorptivity `a`, or else at the effective
absorptivity of the first live cell, and is attenuated along its direction
cosine `μ_h`:

	0 = I_{h,e} − a w_h q
	0 = I_{h,d} − I_{h,u} exp(−τᵢ/μ_h)

and the cell absorbs `A_u I_{h,u} − A_d I_{h,d}`. The beams reaching the
back face join the diffuse flux there:

	0 = q⁺_b − [ε_l Q_b + (1 − ε_l)(q⁻_b + Σ_h I_{h,b})]     (open back)
	0 = q⁺_b − (q⁻_b + Σ_h I_{h,b})                          (mirror)

An isothermal medium facing radiation at its own temperature satisfies every
thermal row with `q⁻ = q⁺ = B` at any optical thickness, and the cell sources
sum to the net flux through the exposed face minus the net flux through the
back.

`beams` holds the beam sources, `up` and `down` the upstream and downstream
face of every cell, `order` the cells from the exposed face into the depth,
`d_floor` the floors of their center-to-face distances. The faces upstream
of the first live cell carry no rows.
"""
struct RadiativeTransferKernel{T<:AbstractThermalBC,B,S<:Tuple} <: Kernel
    top::T
    back::B
    beams::S
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
# a diffuse `RadiativeFluxBC` with its own absorptivity.
_given_emissivity(::AbstractThermalBC) = false
_given_emissivity(bc::RadiativeBC) = bc.ε !== nothing
_given_emissivity(bc::CombinedThermalBC) = any(_given_emissivity, bc.bcs)
_given_absorptivity(::AbstractThermalBC) = false
_given_absorptivity(bc::RadiativeFluxBC) =
    bc.distribution isa Diffuse && bc.absorptivity !== nothing
_given_absorptivity(bc::CombinedThermalBC) = any(_given_absorptivity, bc.bcs)
_diffuse_flux(::AbstractThermalBC) = false
_diffuse_flux(bc::RadiativeFluxBC) = bc.distribution isa Diffuse
_diffuse_flux(bc::CombinedThermalBC) = any(_diffuse_flux, bc.bcs)

_sets_absorptivity(component) =
    hasproperty(component, :surface_absorptivity) &&
    component.surface_absorptivity !== nothing

# The beams of a distribution, of a condition, of a kernel's sources.
_direction_count(::Diffuse) = 0
_direction_count(::Collimated) = 1
_direction_count(c::Cone) = c.directions
_beam_count(::AbstractThermalBC) = 0
_beam_count(bc::RadiativeFluxBC) = _direction_count(bc.distribution)
_beam_count(bc::CombinedThermalBC) = sum(_beam_count, bc.bcs; init = 0)
_beam_count(sources::Tuple) = sum(src -> length(src.μ), sources; init = 0)

# The `RadiativeFluxBC` members of a condition that arrive as beams.
_beam_members(::AbstractThermalBC) = ()
_beam_members(bc::RadiativeFluxBC) = _beam_count(bc) > 0 ? (bc,) : ()
_beam_members(bc::CombinedThermalBC) =
    foldl((members, b) -> (members..., _beam_members(b)...), bc.bcs; init = ())

# Golub–Welsch Gauss–Legendre nodes and weights on [−1, 1].
function _gauss_legendre(n::Int)
    β = [k / sqrt(4.0 * k^2 - 1.0) for k = 1:(n-1)]
    x, V = eigen(SymTridiagonal(zeros(n), β))
    return x, 2 .* abs2.(V[1, :])
end

# The direction cosines of a distribution's beams and the share of the flux
# each carries. Under a cone of uniform intensity the flux per unit μ is
# proportional to μ, so a node's share is its weight times μ over the
# quadrature's ∫ μ dμ, Σ ω μ, which Gauss–Legendre integrates exactly.
_beam_directions(::Collimated) = ([1.0], [1.0])
function _beam_directions(c::Cone)
    μ₁, μ₂ = cos(c.θ₁), cos(c.θ₂)
    x, ω = _gauss_legendre(c.directions)
    μ = (μ₁ + μ₂) / 2 .+ (μ₁ - μ₂) / 2 .* x
    ωμ = ω .* μ
    return μ, ωμ ./ sum(ωμ)
end

# The beam sources of `members`, their rows numbered from `row`.
_beam_sources(::Tuple{}, row) = ()
function _beam_sources(members::Tuple, row)
    m = first(members)
    μ, w = _beam_directions(m.distribution)
    source = BeamSource(m.flux, m.absorptivity, μ, w, row)
    return (source, _beam_sources(Base.tail(members), row + length(μ))...)
end

"""
	incoming_radiation(condition, t)

The diffuse radiation arriving at a face from outside [W/m²]: `F σ T∞⁴` per
`RadiativeBC` member of `condition` and the flux per [`Diffuse`](@ref)
`RadiativeFluxBC` member, zero for every other member.
"""
@inline incoming_radiation(::AbstractThermalBC, t) = 0.0
@inline incoming_radiation(bc::RadiativeBC, t) =
    evaluate_param(bc.F, t) * STEFAN_BOLTZMANN * evaluate_param(bc.T_inf, t)^4
@inline incoming_radiation(bc::RadiativeFluxBC{<:Any,<:Any,Diffuse}, t) =
    evaluate_param(bc.flux, t)
@inline incoming_radiation(bc::CombinedThermalBC, t) = _incoming_sum(bc.bcs, t)
@inline _incoming_sum(::Tuple{}, t) = 0.0
@inline _incoming_sum(b::Tuple, t) =
    incoming_radiation(first(b), t) + _incoming_sum(Base.tail(b), t)

"""
	arriving_radiation(condition, t)

All the radiation arriving at a face from outside [W/m²]: the diffuse
radiation of [`incoming_radiation`](@ref) and the flux of every beam.
"""
arriving_radiation(bc::AbstractThermalBC, t) = incoming_radiation(bc, t)
arriving_radiation(bc::RadiativeFluxBC, t) = evaluate_param(bc.flux, t)
arriving_radiation(bc::CombinedThermalBC, t) = _arriving_sum(bc.bcs, t)
_arriving_sum(::Tuple{}, t) = 0.0
_arriving_sum(b::Tuple, t) =
    arriving_radiation(first(b), t) + _arriving_sum(Base.tail(b), t)

# The largest surroundings temperature of a condition's `RadiativeBC` members
# at time `t`, zero without one.
_surroundings_temperature(::AbstractThermalBC, t) = 0.0
_surroundings_temperature(bc::RadiativeBC, t) = Float64(evaluate_param(bc.T_inf, t))
_surroundings_temperature(bc::CombinedThermalBC, t) =
    maximum(b -> _surroundings_temperature(b, t), bc.bcs; init = 0.0)

"""
	radiative_transfer_rows(bc_set) -> Int

The radiation rows per face under `RADIATIVE_TRANSFER`: in every band, the two
thermal families and one beam family per direction of the beams on `:top`.
"""
function radiative_transfer_rows(bc_set::BoundaryConditionSet)
    thermal = bc_set.thermal
    beams = haskey(thermal, :top) ? _beam_count(thermal[:top]) : 0
    return RADIATION_BANDS * (2 + beams)
end

# The families of a band: q⁻, q⁺ and the beams.
_families(kernel::RadiativeTransferKernel) = 2 + _beam_count(kernel.beams)

# The row of family `j` of band `b` in a face's radiation rows, with `F`
# families per band.
@inline _family_row(F, b, j) = F * (b - 1) + j

"""
	radiation_reach(kernel::RadiativeTransferKernel) -> NamedTuple

The radiation rows each radiation row reads, for the Jacobian's pattern. A
row of family `k` at face `f` reads its own family at `f` and at
`f + along[k]`, the face upstream of it in its direction of travel (toward
the back for `q⁻` and the beams, toward the exposed face for `q⁺`), and the
rows `same_node[k]` at `f`: `q⁻` and `q⁺` read each other there, since any
face can become the exposed face and the back face reflects `q⁻` into `q⁺`.
At the face `back_face` each pair `(k, j)` of `back` is a row `k` that also
reads the row `j`: `q⁺` reads every beam arriving there.
"""
function radiation_reach(kernel::RadiativeTransferKernel)
    c = first(kernel.order)
    toward_exposed = sign(kernel.up[c] - kernel.down[c])
    F = _families(kernel)
    family(k) = mod1(k, F)
    R = RADIATION_BANDS * F
    along = [family(k) == Q_UP ? -toward_exposed : toward_exposed for k = 1:R]
    same_node =
        [family(k) == Q_DOWN ? [k + 1] : family(k) == Q_UP ? [k - 1] : Int[] for k = 1:R]
    back = [
        (_family_row(F, b, Q_UP), _family_row(F, b, j)) for b = 1:RADIATION_BANDS for
        j = (Q_UP+1):F
    ]
    return (; along, same_node, back_face = kernel.down[last(kernel.order)], back)
end

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
any other tag than `:top` and `:bottom`, a beam on any other tag than `:top`,
a `RadiativeBC` with its own `ε` (the interface's transmissivity is the
emissivity of its cell, both ways), and a diffuse `RadiativeFluxBC` with its
own `absorptivity` or on a material where a component sets one (diffuse
radiation enters at the emissivity).
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
        _beam_count(bc) > 0 &&
            tag !== :top &&
            throw(
                ArgumentError(
                    "radiative transfer takes beams through :top only: give the " *
                    "RadiativeFluxBC on :$tag a Diffuse() distribution",
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
                "diffuse RadiativeFluxBC on :$tag no absorptivity",
            ),
        )
    end
    if any(_diffuse_flux, values(thermal))
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
    beams = _beam_sources(_beam_members(thermal[:top]), Q_UP + 1)
    return RadiativeTransferKernel(thermal[:top], back, beams, up, down, order, d)
end

# The position in `order` of the first live cell, zero when none is live.
@inline function _first_live(kernel::RadiativeTransferKernel, ds)
    @inbounds for s in eachindex(kernel.order)
        ds.active[kernel.order[s]] && return s
    end
    return 0
end

"""
	exposed_face(kernel, ds) -> Int

The face the radiation enters through: the upstream face of the first live
cell from the exposed boundary (the boundary face itself while no cell is
dead), zero when no cell is live.
"""
@inline function exposed_face(kernel::RadiativeTransferKernel, ds)
    s = _first_live(kernel, ds)
    return s == 0 ? 0 : kernel.up[kernel.order[s]]
end

# The thickness of cell `c` between its two faces.
@inline function _cell_thickness(kernel::RadiativeTransferKernel, ds, topo, c)
    g = ds.geometry
    return smooth_floor(
        center_to_face_distance(g, topo, kernel.up[c], c),
        kernel.d_floor[1, c],
    ) + smooth_floor(
        center_to_face_distance(g, topo, kernel.down[c], c),
        kernel.d_floor[2, c],
    )
end

# The transmissivity `t` of cell `c` (`exp(−D τ)`), its emissive power `B` and
# its thickness `Δ`.
@inline function _cell_transfer(kernel::RadiativeTransferKernel, ds, topo, c)
    Δ = _cell_thickness(kernel, ds, topo, c)
    return exp(-DIFFUSIVITY * ds.α[c] * Δ), STEFAN_BOLTZMANN * ds.T_trial[c]^4, Δ
end

# The back face's reflectivity `ρ` and injection `J`:
# `q⁺_b = J + ρ (q⁻_b + Σ_h I_{h,b})`.
@inline _back_row(::MirrorBack, ε, t) = (one(ε), zero(ε))
@inline _back_row(back::OpenBack, ε, t) = (1 - ε, ε * incoming_radiation(back.condition, t))

# The fraction of a beam of direction cosine `μ` that crosses the optical
# thickness `τ`.
@inline _beam_transmission(τ, μ) = exp(-τ / μ)

# The absorptivity a beam source enters at: its own, else the first live
# cell's `a`.
@inline _entering_absorptivity(::Nothing, a, t) = a
@inline _entering_absorptivity(absorptivity, a, t) = evaluate_param(absorptivity, t)

# The rows of the beams entering at the face `e`, the band's rows offset by `o`.
@inline _enter_beams!(r, q, ::Tuple{}, o, e, a, t) = r
@inline function _enter_beams!(r, q, sources::Tuple, o, e, a, t)
    src = first(sources)
    entering = _entering_absorptivity(src.absorptivity, a, t) * evaluate_param(src.flux, t)
    @inbounds for h in eachindex(src.μ)
        k = o + src.row + h - 1
        r[k, e] = q[k, e] - src.w[h] * entering
    end
    return _enter_beams!(r, q, Base.tail(sources), o, e, a, t)
end

# The rows of the beams across a cell of optical thickness `τ` from its face
# `u` to its face `d`, and the power the cell absorbs from them.
@inline _cross_beams!(r, q, ::Tuple{}, o, u, d, τ, A_u, A_d) = zero(eltype(q))
@inline function _cross_beams!(r, q, sources::Tuple, o, u, d, τ, A_u, A_d)
    src = first(sources)
    absorbed = _cross_beams!(r, q, Base.tail(sources), o, u, d, τ, A_u, A_d)
    @inbounds for h in eachindex(src.μ)
        k = o + src.row + h - 1
        r[k, d] = q[k, d] - q[k, u] * _beam_transmission(τ, src.μ[h])
        absorbed += q[k, u] * A_u - q[k, d] * A_d
    end
    return absorbed
end

# The flux of every beam at the face `f`.
@inline _beams_at(q, ::Tuple{}, o, f) = zero(eltype(q))
@inline function _beams_at(q, sources::Tuple, o, f)
    src = first(sources)
    arriving = _beams_at(q, Base.tail(sources), o, f)
    @inbounds for h in eachindex(src.μ)
        arriving += q[o+src.row+h-1, f]
    end
    return arriving
end

# Each beam from the exposed face of the first live cell `order[s]` into the
# depth, exactly: the values that satisfy its rows.
_sweep_beams!(q, ::Tuple{}, kernel, ds, topo, o, s, a, t) = q
function _sweep_beams!(q, sources::Tuple, kernel, ds, topo, o, s, a, t)
    src = first(sources)
    order = kernel.order
    e = kernel.up[order[s]]
    entering = _entering_absorptivity(src.absorptivity, a, t) * evaluate_param(src.flux, t)
    for h in eachindex(src.μ)
        k = o + src.row + h - 1
        q[k, e] = src.w[h] * entering
        for m = s:length(order)
            i = order[m]
            τ = ds.α[i] * _cell_thickness(kernel, ds, topo, i)
            q[k, kernel.down[i]] = q[k, kernel.up[i]] * _beam_transmission(τ, src.μ[h])
        end
    end
    return _sweep_beams!(q, Base.tail(sources), kernel, ds, topo, o, s, a, t)
end

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
    beams = kernel.beams
    i₁ = order[s]
    e = kernel.up[i₁]
    l = last(order)
    b_face = kernel.down[l]
    tₑ = smooth_floor(ds.ε[i₁], INTERFACE_FLOOR)
    Q = incoming_radiation(kernel.top, t)
    ρ_b, J_b = _back_row(kernel.back, ds.ε[l], t)
    F = _families(kernel)
    for b = 1:RADIATION_BANDS
        o = _family_row(F, b, 0)
        dn, up = o + Q_DOWN, o + Q_UP
        r[dn, e] = q[dn, e] - (tₑ * Q + (1 - tₑ) * q[up, e])
        _enter_beams!(r, q, beams, o, e, ds.a[i₁], t)
        @inbounds for k = s:length(order)
            i = order[k]
            u, d = kernel.up[i], kernel.down[i]
            tr, B, Δ = _cell_transfer(kernel, ds, topo, i)
            r[dn, d] = q[dn, d] - (B + (q[dn, u] - B) * tr)
            r[up, u] = q[up, u] - (B + (q[up, d] - B) * tr)
            A_u, A_d = g.face_area[u], g.face_area[d]
            c.S_energy[i] += (q[dn, u] - q[up, u]) * A_u - (q[dn, d] - q[up, d]) * A_d
            c.S_energy[i] += _cross_beams!(r, q, beams, o, u, d, ds.α[i] * Δ, A_u, A_d)
        end
        arriving = q[dn, b_face] + _beams_at(q, beams, o, b_face)
        r[up, b_face] = q[up, b_face] - (J_b + ρ_b * arriving)
    end
    return c
end

"""
	radiation_crossing(kernel::RadiativeTransferKernel, ds, t) -> (; absorbed, emitted, back)

The radiation crossing the exposed and back faces at the derived state `ds`
[W/m²]: `absorbed`, the radiation entering through the exposed face, the
diffuse radiation arriving there at the face's transmissivity `tₑ` and every
beam entering; `emitted`, the radiation leaving through it, `absorbed` less
the net intake `q⁻_e − q⁺_e + Σ_h I_{h,e}`, which is `tₑ q⁺_e` where the
exposed row holds; and `back`, the net radiation leaving through the back
face, `q⁻_b + Σ_h I_{h,b} − q⁺_b`, zero at a mirror where its row holds. The
cell sources sum to `A_e (absorbed − emitted) − A_b back` at any fluxes. `NaN`
when no cell is live.
"""
function radiation_crossing(kernel::RadiativeTransferKernel, ds, t)
    q = ds.q
    s = _first_live(kernel, ds)
    s == 0 && return (absorbed = NaN, emitted = NaN, back = NaN)
    i₁ = kernel.order[s]
    e = kernel.up[i₁]
    b_face = kernel.down[last(kernel.order)]
    tₑ = smooth_floor(ds.ε[i₁], INTERFACE_FLOOR)
    Q = incoming_radiation(kernel.top, t)
    F = _families(kernel)
    absorbed, intake, back = zero(eltype(q)), zero(eltype(q)), zero(eltype(q))
    for b = 1:RADIATION_BANDS
        o = _family_row(F, b, 0)
        dn, up = o + Q_DOWN, o + Q_UP
        entering = _beams_at(q, kernel.beams, o, e)
        absorbed += tₑ * Q + entering
        intake += q[dn, e] - q[up, e] + entering
        back += q[dn, b_face] + _beams_at(q, kernel.beams, o, b_face) - q[up, b_face]
    end
    return (; absorbed, emitted = absorbed - intake, back)
end

"""
	consistent_radiation!(q, kernel, ds, topology, t) -> q

The face fluxes `q` (radiation rows × faces) that satisfy every radiation row
at the derived state `ds`, and zero on the faces upstream of the first live
cell. Each beam is carried down from the exposed face first, exactly. The
thermal rows are then linear in the fluxes: with `q⁻_e = a` as a parameter, a
downward pass gives `q⁻_b`, the back row `q⁺_b` (with the beams arriving
there), an upward pass `q⁺_e = A a + C`, the exposed row `a`; a downward and
an upward pass write the values. No iteration, no allocation.
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
    F = _families(kernel)
    for b = 1:RADIATION_BANDS
        o = _family_row(F, b, 0)
        dn, up = o + Q_DOWN, o + Q_UP
        _sweep_beams!(q, kernel.beams, kernel, ds, topo, o, s, ds.a[i₁], t)
        beams_b = _beams_at(q, kernel.beams, o, b_face)
        A, C = one(tₑ), zero(tₑ)
        for k = s:n
            tr, B = _cell_transfer(kernel, ds, topo, order[k])
            A, C = A * tr, B + (C - B) * tr
        end
        A, C = ρ_b * A, J_b + ρ_b * (C + beams_b)
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
        q[up, b_face] = J_b + ρ_b * (q[dn, b_face] + beams_b)
        for k = n:-1:s
            i = order[k]
            tr, B = _cell_transfer(kernel, ds, topo, i)
            q[up, kernel.up[i]] = B + (q[up, kernel.down[i]] - B) * tr
        end
    end
    return q
end
