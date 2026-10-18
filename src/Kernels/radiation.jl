# ============================================================================
# In-depth radiation: Beer–Lambert with lifted face intensities
# ============================================================================

"""
	BeerLambertKernel{R<:RadiativeFluxBC} <: Kernel

In-depth absorption of the radiative flux incident on the exposed face, with
the intensity of every face an algebraic unknown of the state.

The rows, one per face:

- exposed face: `0 = I_e − a I_inc(t)`, with `I_inc` the condition's incident
  flux, `a` its explicit absorptivity folded into `I_inc` or, without one,
  the effective absorptivity of the exposed cell (Kirchhoff's law unless a
  component sets its own);
- every other face, the face of its cell farther from the exposed face:
  `0 = I_down − I_up exp(−α Δ)`, with `α` the cell's absorption coefficient at
  the positive-part concentrations and `Δ` its thickness.

The cell absorbs `I_up A_up − I_down A_down`, the exact cell average of the
attenuated beam; the flux leaving the last face is transmitted. Each row
reads one cell and its two faces, so the Jacobian stays banded where the
eliminated recurrence would couple every cell to every cell above it.

`up` and `down` hold the upstream and downstream face of every cell,
`d_floor` the floors of their center-to-face distances.
"""
struct BeerLambertKernel{R<:RadiativeFluxBC} <: Kernel
    source::R
    exposed::Int
    up::Vector{Int}
    down::Vector{Int}
    order::Vector{Int}
    d_floor::Matrix{Float64}
end

"""
	BeerLambertKernel(topology, geometry, bc_set; tag = :top)

The kernel for the radiative flux condition on the boundary tagged `tag`,
whose single face is the exposed face; the beam then crosses the cells in
order of their distance from it.
"""
function BeerLambertKernel(
    topo::MeshTopology{1},
    g::MeshGeometry,
    bc_set::BoundaryConditionSet;
    tag::Symbol = :top,
)
    for other in keys(bc_set.thermal)
        other === tag && continue
        get_radiative_flux_bc(bc_set.thermal[other]) === nothing || throw(
            ArgumentError(
                "in-depth radiation enters through the :$tag faces only; the condition " *
                "on :$other carries an incident radiative flux",
            ),
        )
    end
    source = get_radiative_flux_bc(bc_set.thermal[tag])
    source === nothing && (source = RadiativeFluxBC(0.0))
    faces = boundary_faces(topo, tag)
    length(faces) == 1 || throw(
        ArgumentError(
            "in-depth radiation needs one exposed face (:$tag has $(length(faces)))",
        ),
    )
    exposed = only(faces)
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
        throw(ArgumentError("the beam from :$tag does not cross every cell"))
    d = zeros(2, n)
    for c = 1:n
        d[1, c] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, up[c], c)
        d[2, c] = DISTANCE_FLOOR_FRACTION * center_to_face_distance(g, topo, down[c], c)
    end
    return BeerLambertKernel(source, exposed, up, down, order, d)
end

"""
	incident_intensity(kernel, ds, t)

The intensity entering through the exposed face [W/m²]: the incident flux,
times the condition's explicit absorptivity or else the exposed cell's
effective absorptivity.
"""
@inline function incident_intensity(kernel::BeerLambertKernel, ds, t)
    src = kernel.source
    q = evaluate_param(src.flux, t) + zero(eltype(ds.a))
    a =
        src.absorptivity === nothing ? ds.a[first(kernel.order)] :
        evaluate_param(src.absorptivity, t) + zero(eltype(ds.a))
    return a * q
end

# Optical thickness of cell `c`.
@inline function _optical_thickness(kernel::BeerLambertKernel, ds, topo, c)
    g = ds.geometry
    Δ =
        smooth_floor(
            center_to_face_distance(g, topo, kernel.up[c], c),
            kernel.d_floor[1, c],
        ) + smooth_floor(
            center_to_face_distance(g, topo, kernel.down[c], c),
            kernel.d_floor[2, c],
        )
    return ds.α[c] * Δ
end

function contribute!(
    c::Contributions,
    kernel::BeerLambertKernel,
    ds::DerivedState,
    topo::MeshTopology{1},
    material,
    t,
)
    g = ds.geometry
    e = kernel.exposed
    c.r_intensity[e] = ds.I[e] - incident_intensity(kernel, ds, t)
    @inbounds for i in kernel.order
        up, down = kernel.up[i], kernel.down[i]
        I_up, I_down = ds.I[up], ds.I[down]
        c.r_intensity[down] = I_down - I_up * exp(-_optical_thickness(kernel, ds, topo, i))
        c.S_energy[i] += I_up * g.face_area[up] - I_down * g.face_area[down]
    end
    return c
end

"""
	sweep_intensities!(I, kernel, ds, topology, t) -> I

The face intensities that satisfy every intensity row at the derived state
`ds`: the incident intensity at the exposed face, attenuated cell by cell
into the depth. The consistent initial values of the algebraic rows.
"""
function sweep_intensities!(
    I::AbstractVector,
    kernel::BeerLambertKernel,
    ds::DerivedState,
    topo::MeshTopology{1},
    t,
)
    I[kernel.exposed] = incident_intensity(kernel, ds, t)
    for i in kernel.order
        I[kernel.down[i]] = I[kernel.up[i]] * exp(-_optical_thickness(kernel, ds, topo, i))
    end
    return I
end

"""
	transmitted_face(kernel) -> Int

The face the beam leaves through after crossing every cell.
"""
transmitted_face(kernel::BeerLambertKernel) = kernel.down[last(kernel.order)]
