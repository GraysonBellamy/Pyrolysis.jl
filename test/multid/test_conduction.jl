using Test
using LinearAlgebra
using StaticArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.Meshes: MeshGeometry
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels
using Pyrolysis.Kernels: Contributions, reset!, contribute!
using Pyrolysis.Properties: guarded_harmonic_mean

_rot(θ) = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
_rot3(θ, ϕ) =
    (@SMatrix [cos(θ) -sin(θ) 0; sin(θ) cos(θ) 0; 0 0 1]) *
    (@SMatrix [1 0 0; 0 cos(ϕ) -sin(ϕ); 0 sin(ϕ) cos(ϕ)])
_tensor(θ) = _rot(θ) * SMatrix{2,2}(Diagonal(SVector(4.0, 1.0))) * _rot(θ)'
_tensor3(θ, ϕ) = _rot3(θ, ϕ) * SMatrix{3,3}(Diagonal(SVector(4.0, 1.0, 2.0))) * _rot3(θ, ϕ)'

_shear(x::SVector{2}) = SVector(x[1] + 0.3x[2], x[2])
_shear(x::SVector{3}) = SVector(x[1] + 0.3x[3], x[2] - 0.2x[3], x[3])
_wave(x::SVector{2}) =
    SVector(x[1] + 0.08sin(π * x[1]) * sin(2π * x[2]), x[2] + 0.06sin(2π * x[1]) * x[2])

function _bent()
    m = box_mesh(0:0.5:1.0, 0:0.5:1.0, 0:0.5:1.0; transform = _shear)
    p = copy(m.positions)
    p[14] += SVector(0.1, -0.07, 0.05)
    cells = [m.topology.cells[i].nodes for i = 1:n_cells(m.topology)]
    return cell_mesh(
        p,
        cells;
        tag = (i, k, s) -> ((:left, :right), (:front, :back), (:bottom, :top))[k][s+1],
    )
end

# A single solid at uniform density on `mesh` with cell temperatures
# `T(center)`: the derived state, the conduction kernel with every cell's
# tensor `A` (one for all, or a function of the center), the contributions
# after it ran, and the thermal boundary kernel sharing its tensors. Given
# `factor`, a function of the center, the cells' conductivities are set to it
# before the kernel runs.
function _conduction(mesh, T, A; k = 2.0, factor = nothing)
    s = SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = k)
    material = Material(name = :solid, components = (s,), reactions = ())
    topo = mesh.topology
    D = length(first(mesh.positions))
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, reduce(vcat, mesh.positions), 1.0, topo)
    n = n_cells(topo)
    m0 = reshape(1000.0 .* g0.cell_volume, 1, n)
    reference = CellReference(mesh, m0, material; reltol = 1e-6, threshold = 1e-2)
    setup = StateSetup(topo, material, reference, nothing)
    layout = StateLayout(1, topo, TemperatureState())
    T0 = [T(g0.cell_center[i]) for i = 1:n]
    ds = DerivedState{Float64}(layout, setup, copy(T0))
    u = zeros(state_length(layout))
    for i = 1:n
        u[mass_index(layout, i, 1)] = m0[1, i]
        u[energy_index(layout, i)] = T0[i]
    end
    derive!(ds, u, layout, setup)
    factor === nothing || (ds.k_physical .= [factor(g0.cell_center[i]) for i = 1:n])
    tensor = A isa Function ? [A(g0.cell_center[i]) for i = 1:n] : fill(A, n)
    kernel = ConductionKernel(topo, g0; tensor)
    c = Contributions{Float64}(1, 0, n, n_faces(topo), n_nodes(topo), D)
    reset!(c)
    contribute!(c, kernel, ds, topo, material, 0.0)
    tags = boundary_tags(topo)
    bc_set = BoundaryConditionSet(
        thermal = NamedTuple{Tuple(tags)}(Tuple(AdiabaticBC() for _ in tags)),
    )
    boundary = ThermalBoundaryKernel(topo, g0, bc_set; tensor)
    return (; material, topo, ds, c, kernel, boundary, g = ds.geometry)
end

# Interior and boundary fluxes of the linear field T = T₀ + b·x against the
# exact −n·K·b, and the cells' gradients against b.
function _linear_errors(mesh, A, b; k = 2.0)
    T(x) = 600.0 + dot(b, x)
    r = _conduction(mesh, T, A; k)
    (; topo, ds, c, g) = r
    D = length(b)
    scale = k * opnorm(Matrix(A)) * norm(b)
    interior, boundary, temperature, gradient = 0.0, 0.0, 0.0, 0.0
    for i = 1:n_cells(topo)
        gradient = max(gradient, norm(cell_gradient(c, Val(D), i) - b))
    end
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        n = g.face_normal[f]
        if nb != 0
            d_o = center_to_face_distance(g, topo, f, o)
            d_n = center_to_face_distance(g, topo, f, nb)
            K = guarded_harmonic_mean(k, k, d_o, d_n) * A
            interior = max(interior, abs(c.F_energy[f] + dot(n, K * b)))
        else
            exact = -dot(n, k * A * b)
            Tf = T(g.face_center[f])
            q, Ts = boundary_energy_flux(
                DirichletThermalBC(Tf),
                f,
                r.boundary,
                c,
                ds,
                topo,
                r.material,
                0.0,
            )
            boundary = max(boundary, abs(q - exact))
            @test Ts == Tf
            q, Ts = boundary_energy_flux(
                HeatFluxBC(-exact),
                f,
                r.boundary,
                c,
                ds,
                topo,
                r.material,
                0.0,
            )
            @test q == exact
            temperature = max(temperature, abs(Ts - Tf))
        end
    end
    return (
        interior = interior / scale,
        boundary = boundary / scale,
        temperature,
        gradient = gradient / norm(b),
    )
end

@testset "linear fields: exact fluxes with a rotated diag(4, 1) on $name" for (
    name,
    mesh,
) in (
    ("a sheared quadrilateral mesh", box_mesh(0:0.25:1.0, 0:0.2:1.0; transform = _shear)),
    ("a distorted quadrilateral mesh", box_mesh(0:0.2:1.0, 0:0.25:1.0; transform = _wave)),
)
    b = SVector(300.0, -200.0)
    for θ in (π / 6, π / 2)
        e = _linear_errors(mesh, _tensor(θ), b)
        @test e.interior < 1e-13
        @test e.boundary < 1e-13
        @test e.temperature < 1e-10
        @test e.gradient < 1e-12
    end
end

@testset "linear fields: exact fluxes with a rotated diag(4, 1, 2) on $name" for (
    name,
    mesh,
) in (
    ("sheared hexahedra", box_mesh(0:0.5:1.0, 0:0.5:1.0, 0:0.25:1.0; transform = _shear)),
    ("hexahedra with bent faces", _bent()),
)
    b = SVector(300.0, -200.0, 150.0)
    for (θ, ϕ) in ((π / 6, π / 5), (π / 2, 0.0))
        e = _linear_errors(mesh, _tensor3(θ, ϕ), b)
        @test e.interior < 1e-13
        @test e.boundary < 1e-13
        @test e.temperature < 1e-10
        @test e.gradient < 1e-12
    end
end

@testset "the cross part carries what the two-point part misses" begin
    mesh = box_mesh(0:0.25:1.0, 0:0.2:1.0; transform = _shear)
    b = SVector(300.0, -200.0)
    A = _tensor(π / 6)
    r = _conduction(mesh, x -> 600.0 + dot(b, x), A)
    (; topo, ds, g) = r
    worst = 0.0
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        nb == 0 && continue
        n = g.face_normal[f]
        d_o = center_to_face_distance(g, topo, f, o)
        d_n = center_to_face_distance(g, topo, f, nb)
        K = guarded_harmonic_mean(2.0, 2.0, d_o, d_n) * A
        two_point = -dot(n, K * n) * (ds.T[nb] - ds.T[o]) / (d_o + d_n)
        worst = max(worst, abs(two_point + dot(n, K * b)) / abs(dot(n, K * b)))
    end
    @test worst > 0.1
end

@testset "a column one cell wide: finite gradients, no cross flux on rectangles" begin
    mesh = box_mesh(0:1.0, 0:0.2:1.0)
    b = SVector(0.0, -200.0)
    r = _conduction(mesh, x -> 600.0 + dot(b, x), SMatrix{2,2}(1.0I))
    (; topo, c) = r
    for i = 1:n_cells(topo)
        ∇T = cell_gradient(c, Val(2), i)
        @test all(isfinite, ∇T)
        @test ∇T[2] ≈ b[2] rtol = 1e-7
    end
    @test all(
        f -> topo.faces[f].cells[2] == 0 || c.F_energy[f] ≈ 2.0 * 200.0,
        1:n_faces(topo),
    )
end

# ============================================================================
# Material interfaces
# ============================================================================

_sym(M) = (M + M') / 2
_shear2(s) = x -> SVector(x[1] + s * x[2], x[2])
_shear3(s) = x -> SVector(x[1] + s * x[2] + 0.3s * x[3], x[2], x[3])
const _Ra = _rot(27π / 180) * SMatrix{2,2}(Diagonal(SVector(4.0, 1.0))) * _rot(27π / 180)'
const _Rb = _rot(-19π / 180) * SMatrix{2,2}(Diagonal(SVector(3.0, 1.0))) * _rot(-19π / 180)'
const _A1 = _sym(_tensor3(π / 6, π / 5))
const _A2 = _sym(
    _rot3(-π / 7, π / 9) *
    SMatrix{3,3}(Diagonal(SVector(3.0, 1.0, 1.5))) *
    _rot3(-π / 7, π / 9)',
)

# A field affine on either side of the plane `n·x = c₀`, with the tangential
# gradient `t` and the normal flux `q` common to both sides: continuous in
# temperature and normal flux, `K₁` where `n·x < c₀` and `K₂` beyond. Returns
# the field, the gradient under a tensor, and the side test.
function _transmission_field(n, c0, K1, K2; t, q = 5.0)
    t = t - dot(t, n) * n
    grad(K) = t + n * ((-q - dot(n, K * t)) / dot(n, K * n))
    first_side(x) = dot(n, x) < c0
    T(x) = 500.0 + dot(t, x) + (dot(n, x) - c0) * dot(n, grad(first_side(x) ? K1 : K2))
    return T, grad, first_side
end

# The largest interior-face flux error and cell-gradient error of the kernel
# on the transmission field across the plane through `origin` with normal `n`,
# the conductivity `k₁ A₁` on its first side and `k₂ A₂` beyond, relative to
# the largest exact flux and gradient.
function _interface_errors(mesh, n, origin, (A1, k1), (A2, k2); t)
    T, grad, first_side = _transmission_field(n, dot(n, origin), k1 * A1, k2 * A2; t)
    r = _conduction(
        mesh,
        T,
        x -> first_side(x) ? A1 : A2;
        factor = x -> first_side(x) ? k1 : k2,
    )
    (; topo, c, g) = r
    K(i) = first_side(g.cell_center[i]) ? k1 * A1 : k2 * A2
    D = length(n)
    flux, flux_scale, gradient, gradient_scale = 0.0, 0.0, 0.0, 0.0
    for i = 1:n_cells(topo)
        exact = grad(K(i))
        gradient = max(gradient, norm(cell_gradient(c, Val(D), i) - exact))
        gradient_scale = max(gradient_scale, norm(exact))
    end
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        nb == 0 && continue
        exact = -dot(g.face_normal[f], K(o) * grad(K(o)))
        flux = max(flux, abs(c.F_energy[f] - exact))
        flux_scale = max(flux_scale, abs(exact))
    end
    return (flux = flux / flux_scale, gradient = gradient / gradient_scale)
end

@testset "piecewise-affine fields across a planar interface: exact on $(name), shear $(s)" for (
        name,
        pair1,
        pair2,
    ) in (
        ("a 100× scalar jump with one rotated tensor", (_Ra, 1.0), (_Ra, 100.0)),
        ("a rotated-tensor jump", (_Ra, 1.0), (_Rb, 100.0)),
        (
            "I against diag(100, 1)",
            (SMatrix{2,2}(1.0I), 1.0),
            (SMatrix{2,2}(Diagonal(SVector(100.0, 1.0))), 1.0),
        ),
    ),
    s in (0.0, 0.35)

    mesh = box_mesh(0:1.0:4.0, 0:1.0:4.0; transform = _shear2(s))
    e = _interface_errors(
        mesh,
        normalize(SVector(1.0, -s)),
        SVector(2.0, 0.0),
        pair1,
        pair2;
        t = SVector(1.0, 1.0),
    )
    @test e.flux < 1e-12
    @test e.gradient < 1e-12
end

@testset "piecewise-affine fields across a planar interface: exact on sheared hexahedra with $(name), shear $(s)" for (
        name,
        pair1,
        pair2,
    ) in (
        ("a 100× scalar jump", (_A1, 1.0), (_A1, 100.0)),
        ("a rotated-tensor jump", (_A1, 1.0), (_A2, 100.0)),
    ),
    s in (0.0, 0.35)

    mesh = box_mesh(0:1.0:4.0, 0:1.0:4.0, 0:1.0:2.0; transform = _shear3(s))
    e = _interface_errors(
        mesh,
        normalize(SVector(1.0, -s, -0.3s)),
        SVector(2.0, 0.0, 0.0),
        pair1,
        pair2;
        t = SVector(1.0, 1.0, -0.5),
    )
    @test e.flux < 1e-12
    @test e.gradient < 1e-12
end

# Two interfaces meeting: K∇T = (0, 1) everywhere, the flux along every face
# normal −n₂, on 4×4 unit squares.
@testset "the $(name): exact" for (name, K, T) in (
    (
        "ridge between two sheared tensors",
        x -> x[1] < 2 ? (@SMatrix [1.0 -2.0; -2.0 5.0]) : (@SMatrix [1.0 2.0; 2.0 5.0]),
        x -> 500.0 + x[2] - 2abs(x[1] - 2),
    ),
    (
        "junction of four materials",
        x -> begin
            a = x[1] < 2 ? 2.0 : -2.0
            b = x[2] < 2 ? 1.0 : 3.0
            @SMatrix [1.0 -a/b; -a/b 1/b+(a/b)^2]
        end,
        x -> 500.0 - 2abs(x[1] - 2) + (x[2] < 2 ? x[2] : 2 + 3(x[2] - 2)),
    ),
)
    r = _conduction(box_mesh(0:1.0:4.0, 0:1.0:4.0), T, K; factor = x -> 1.0)
    (; topo, c, g) = r
    worst = 0.0
    for f = 1:n_faces(topo)
        topo.faces[f].cells[2] == 0 && continue
        worst = max(worst, abs(c.F_energy[f] + g.face_normal[f][2]))
    end
    @test worst < 1e-12
end

# Whether `x` is the point `y`, up to the round-off of computed centroids.
_at(x, y) = norm(x - y) < 1e-12

# The cell centered at `x`.
_cell(r, x) = only(i for i = 1:n_cells(r.topo) if _at(r.g.cell_center[i], x))

# The face between the center of 3×3 unit squares and its right neighbour,
# whose conductivity is `ε` times the others'.
function _insulating(ε)
    T(x) = 600.0 + 30x[1] - 20x[2] + 5x[1] * x[2]
    r = _conduction(
        box_mesh(0:1.0:3.0, 0:1.0:3.0),
        T,
        _Ra;
        factor = x -> _at(x, SVector(2.5, 1.5)) ? ε : 1.0,
    )
    (; topo, c) = r
    i, j = _cell(r, SVector(1.5, 1.5)), _cell(r, SVector(2.5, 1.5))
    f = only(f for f in topo.cells[i].faces if j in topo.faces[f].cells)
    return (; r, i, j, f, q = c.F_energy[f], ∇T = cell_gradient(c, Val(2), i))
end

@testset "an insulating neighbour: the face flux vanishes with it, the gradient is continuous through the limit row" begin
    fluxes = [abs(_insulating(ε).q) / ε for ε in (1e-3, 1e-6, 1e-9)]
    @test maximum(fluxes) < 2 * minimum(fluxes)
    zero_limit = _insulating(0.0)
    @test zero_limit.q == 0
    for k in (zero_limit.r.c.F_energy, zero_limit.r.c.gradient)
        @test all(isfinite, k)
    end
    # The two cells share one dimensionless tensor, so the neighbour's normal
    # conductivity is ε times the cell's: through the switch to the limit row,
    # and from above the switch to zero.
    τ = Pyrolysis.Kernels.INSULATING_NEIGHBOUR
    above, below = _insulating(τ * (1 + 1e-6)), _insulating(τ * (1 - 1e-6))
    scale = norm(above.∇T)
    @test norm(above.∇T - below.∇T) < 1e-9 * scale
    @test norm(_insulating(1e-11).∇T - zero_limit.∇T) < 1e-9 * scale
end

@testset "an insulating cell at the boundary: no conducted heat, its trace pinned" begin
    T(x) = 600.0 + 30x[1] - 20x[2]
    r = _conduction(
        box_mesh(0:1.0:3.0, 0:1.0:3.0),
        T,
        _Ra;
        factor = x -> _at(x, SVector(0.5, 0.5)) ? 0.0 : 1.0,
    )
    (; topo, ds, c, boundary, material) = r
    i = _cell(r, SVector(0.5, 0.5))
    @test all(isfinite, c.F_energy) && all(isfinite, c.gradient)
    for f in topo.cells[i].faces
        if topo.faces[f].cells[2] == 0
            @test boundary_conductivity(boundary, c, ds, topo, f) == 0
            for (bc, T_s) in (
                (HeatFluxBC(1e4), ds.T[i]),
                (DirichletThermalBC(350.0), 350.0),
                (AdiabaticBC(), ds.T[i]),
                (RadiativeFluxBC(5e4) + RadiativeBC(T_inf = 300.0), ds.T[i]),
            )
                q, Ts = boundary_energy_flux(bc, f, boundary, c, ds, topo, material, 0.0)
                @test q == 0
                @test Ts == T_s
                c.T_face[f] = Ts
                @test conducted_flux(boundary, c, ds, topo, f) == 0
            end
        else
            @test c.F_energy[f] == 0
        end
    end
end

# The center of 3×3 unit squares with K = I, its x-neighbours under
# [1 c; c 1+c²] and its y-neighbours under [1+c² c; c 1]: the neighbours'
# images are collinear at c = 2, where no fit over them determines the
# gradient.
function _collinear(c)
    Kx = @SMatrix [1.0 c; c 1+c^2]
    Ky = @SMatrix [1+c^2 c; c 1.0]
    K(x) =
        _at(x, SVector(0.5, 1.5)) || _at(x, SVector(2.5, 1.5)) ? Kx :
        _at(x, SVector(1.5, 0.5)) || _at(x, SVector(1.5, 2.5)) ? Ky : SMatrix{2,2}(1.0I)
    r = _conduction(
        box_mesh(0:1.0:3.0, 0:1.0:3.0),
        x -> 600.0 + 30x[1] - 20x[2] + 5x[1] * x[2],
        K;
        factor = x -> 1.0,
    )
    ∇T = cell_gradient(r.c, Val(2), _cell(r, SVector(1.5, 1.5)))
    return (; ∇T, support = gradient_support(r.kernel, r.ds, r.topo))
end

# The center's regularized fit in BigFloat, from the closed-form images of its
# neighbours, (1, −c/2), (−1, c/2), (−c/2, 1) and (c/2, −1), and the
# temperature differences of the field of `_collinear`.
function _collinear_reference(c)
    c = BigFloat(c)
    data = (
        (SVector(1, -c / 2), 37.5),
        (SVector(-1, c / 2), -37.5),
        (SVector(-c / 2, 1), -12.5),
        (SVector(c / 2, -1), 12.5),
    )
    M = sum(d * d' / dot(d, d) for (d, _) in data)
    b = sum(ΔT * d / dot(d, d) for (d, ΔT) in data)
    x = det(M) / (tr(M) / 2)^2
    span, support = Pyrolysis.Kernels.GRADIENT_SPAN, Pyrolysis.Kernels.GRADIENT_SUPPORT
    s = clamp((support - x) / (support - span), 0, 1)
    λ = Pyrolysis.Kernels.GRADIENT_REGULARIZATION * tr(M) * s^2 * (3 - 2s)
    return (M + λ * I) \ b
end

# At c = 2, and where the normalized determinant of the center's fit crosses
# GRADIENT_SPAN: the gradient is the regularized fit, and it converges to its
# value there from either side, at a one-sided rate that matches the centered
# one; a step in the regularization or a wrong value at the point itself would
# grow the quotients as the inverse halfwidth.
@testset "neighbour images that turn collinear: a continuous gradient, and the cell is reported" begin
    for c in (2.0, 2sqrt((1 - 1e-3) / (1 + 1e-3)))
        at = _collinear(c).∇T
        @test at ≈ Float64.(_collinear_reference(c)) rtol = 1e-6
        for h in (1e-5, 1e-7, 1e-9)
            above, below = _collinear(c + h).∇T, _collinear(c - h).∇T
            centered = norm(above - below) / 2h
            @test norm(above - at) / h ≈ centered rtol = 0.02
            @test norm(at - below) / h ≈ centered rtol = 0.02
        end
    end
    at = _collinear(2.0)
    @test all(isfinite, at.∇T)
    @test at.support.cells >= 1
    @test at.support.worst < 1e-12
    @test _collinear(0.0).support.cells == 0
end

# The frozen rate matrix of the kernel at unit heat capacity per volume.
function _rate_matrix(r)
    (; topo, ds, c, g, kernel, material) = r
    n = n_cells(topo)
    R = zeros(n, n)
    for j = 1:n
        ds.T .= 0.0
        ds.T[j] = 1.0
        reset!(c)
        contribute!(c, kernel, ds, topo, material, 0.0)
        for f = 1:n_faces(topo)
            o, nb = topo.faces[f].cells
            nb == 0 && continue
            q = c.F_energy[f] * g.face_area[f]
            R[o, j] -= q / g.cell_volume[o]
            R[nb, j] += q / g.cell_volume[nb]
        end
    end
    return R
end

# Bound preservation is a separate property: rotated anisotropy gives negative
# coupling coefficients. The record below documents them; it asserts only
# conservation.
@testset "monotonicity record: $(name)" for (name, mesh, K) in (
    (
        "a uniform diag(4, 1) rotated by 30° on 5×5 unit squares",
        box_mesh(0:1.0:5.0, 0:1.0:5.0),
        x -> _rot(π / 6) * SMatrix{2,2}(Diagonal(SVector(4.0, 1.0))) * _rot(π / 6)',
    ),
    (
        "an SPD checkerboard on 4×4 unit squares",
        box_mesh(0:1.0:4.0, 0:1.0:4.0),
        x ->
            isodd(floor(Int, x[1]) + floor(Int, x[2])) ? (@SMatrix [1.0 -1.5; -1.5 4.0]) :
            (@SMatrix [4.0 1.5; 1.5 1.0]),
    ),
)
    R = _rate_matrix(_conduction(mesh, x -> 600.0, K; factor = x -> 1.0))
    off = [R[i, j] for i in axes(R, 1) for j in axes(R, 2) if i != j]
    @info "monotonicity record: $(name)" negative = count(<(-1e-12), off) of = length(off) most_negative =
        minimum(off)
    @test maximum(abs, sum(R; dims = 1)) < 1e-12
    @test maximum(abs, sum(R; dims = 2)) < 1e-12
end

# ============================================================================
# A steady manufactured solution
# ============================================================================

using SparseArrays
using Pyrolysis.Assembly: residual!, VolumetricSources
using Pyrolysis.ColoredJacobian: jacobian!, jac_prototype
using Pyrolysis.Integration: setup

# T = 300 + 100 sin(πx) sin(πy) on the unit square, held at 300 K on every
# face, under the conductivity k A, A = R diag(4, 1) Rᵀ rotated by 30°: the
# source −∇·(kA∇T) / (ρc), evaluated at the cell centroids, makes it steady.
# ρc is well above the heat-capacity floor.
const _ρ, _c, _k = 1000.0, 1000.0, 1000.0
const _K = _tensor(π / 6)
_exact(x) = 300.0 + 100.0 * sin(π * x[1]) * sin(π * x[2])
function _source(x, t)
    s1, s2 = sin(π * x[1]), sin(π * x[2])
    c1, c2 = cos(π * x[1]), cos(π * x[2])
    div = 100π^2 * (-(_K[1, 1] + _K[2, 2]) * s1 * s2 + 2 * _K[1, 2] * c1 * c2)
    return -_k * div / (_ρ * _c)
end
# A smooth distortion that keeps the square's outline.
_square(x::SVector{2}) = SVector(
    x[1] + 0.04sin(π * x[1]) * sin(2π * x[2]),
    x[2] + 0.04sin(π * x[2]) * sin(2π * x[1]),
)

function _poisson_error(n)
    s = SolidComponent(:s, ρ = _ρ, c = _c, k = _k)
    held = DirichletThermalBC(300.0)
    problem = PyrolysisProblem(
        mesh = box_mesh(
            range(0.0, 1.0; length = n + 1),
            range(0.0, 1.0; length = n + 1);
            transform = _square,
            fixed = true,
        ),
        material = Material(name = :conductor, components = (s,), reactions = ()),
        bc_set = BoundaryConditionSet(
            thermal = (left = held, right = held, bottom = held, top = held),
        ),
        T_initial = x -> 300.0,
        ξ_initial = [x -> _ρ],
        tspan = (0.0, 1.0),
    )
    cp = setup(problem; tensor = _K, sources = VolumetricSources(nothing, _source))
    layout = cp.def.layout
    rows = [energy_index(layout, i) for i = 1:n_cells(layout)]
    u = copy(cp.u0)
    du = zeros(length(u))
    J = jac_prototype(cp.jacobian)
    for _ = 1:2
        residual!(du, u, cp.def, 0.0, cp.ws)
        jacobian!(J, cp.jacobian, u, 0.0)
        u[rows] .-= J[rows, rows] \ du[rows]
    end
    residual!(du, u, cp.def, 0.0, cp.ws)
    @test maximum(abs, du[rows]) <= 1e-8 * maximum(abs, _source(SVector(0.5, 0.5), 0.0))
    g = cp.ws.derived.geometry
    e = [u[rows[i]] - _exact(g.cell_center[i]) for i = 1:n_cells(layout)]
    return sqrt(sum(e .^ 2 .* g.cell_volume) / sum(g.cell_volume)), maximum(abs, e)
end

@testset "a steady manufactured solution with a rotated tensor on a distorted mesh: second order" begin
    errors = [_poisson_error(n) for n in (8, 16, 32, 64)]
    l2 = first.(errors)
    linf = last.(errors)
    orders = log2.(l2[1:(end-1)] ./ l2[2:end])
    @info "steady manufactured solution at 8, 16, 32 and 64 cells a side" l2 =
        join(l2, ", ") linf = join(linf, ", ") orders = join(orders, ", ")
    @test all(>=(1.9), orders)
    # The largest error is one-sided at the coarsest pair and second order once
    # the mesh resolves the distortion.
    @test log2(linf[end-1] / linf[end]) >= 1.9
end

# ============================================================================
# Material interfaces through setup
# ============================================================================

using Pyrolysis.ColoredJacobian: dense_jacobian, missing_entries, compare_entries
using Pyrolysis.Kernels: find_kernel
using Pyrolysis.Outputs: PyrolysisSolution, SolverStatistics
using ForwardDiff

const _mixed = (
    left = DirichletThermalBC(350.0),
    right = HeatFluxBC(1e3),
    bottom = AdiabaticBC(),
    top = RadiativeFluxBC(5e4) + RadiativeBC(T_inf = 300.0),
)
const _convective = (
    left = DirichletThermalBC(350.0),
    right = HeatFluxBC(1e3),
    bottom = AdiabaticBC(),
    top = ConvectiveBC(h = 10.0, T_inf = 300.0),
)
const _insulated = (
    left = AdiabaticBC(),
    right = AdiabaticBC(),
    bottom = AdiabaticBC(),
    top = AdiabaticBC(),
)

# A solid of conductivity `k` (a number or a function of temperature) under
# the per-cell tensors `tensor` on `mesh` at the temperatures `T`, under the
# thermal conditions `thermal` (by default held at 350 K on the left, given
# 1 kW/m² on the right, insulated below and heated radiatively above), set up
# with `options`.
function _interface_setup(
    mesh,
    tensor;
    k = 1.0,
    T = x -> 500.0 + 30x[1] - 20x[2] + 5x[1] * x[2],
    thermal = _mixed,
    energy = :temperature,
    options...,
)
    s = SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = k)
    problem = PyrolysisProblem(
        mesh = mesh,
        material = Material(name = :solid, components = (s,), reactions = ()),
        bc_set = BoundaryConditionSet(; thermal),
        T_initial = T,
        ξ_initial = [x -> 1000.0],
        tspan = (0.0, 1.0),
    )
    return setup(problem; energy_state = energy, tensor, options...)
end

# The residual of a set-up problem at its initial state, with its
# contributions.
function _residual(cp; t = 0.0)
    du = zeros(length(cp.u0))
    residual!(du, cp.u0, cp.def, t, cp.ws)
    return du, cp.ws.contributions
end

@testset "two unit squares, I against diag(100, 1): the series-resistance flux" begin
    mesh = box_mesh([0.0, 1.0, 2.0], [0.0, 1.0])
    tags = boundary_tags(mesh.topology)
    problem = PyrolysisProblem(
        mesh = mesh,
        material = Material(
            name = :solid,
            components = (SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 1.0),),
            reactions = (),
        ),
        bc_set = BoundaryConditionSet(
            thermal = NamedTuple{Tuple(tags)}(Tuple(AdiabaticBC() for _ in tags)),
        ),
        T_initial = x -> x[1] < 1 ? 300.0 : 400.0,
        ξ_initial = [x -> 1000.0],
        tspan = (0.0, 1.0),
    )
    cp = setup(
        problem;
        tensor = x ->
            x[1] < 1 ? SMatrix{2,2}(1.0I) : SMatrix{2,2}(Diagonal(SVector(100.0, 1.0))),
    )
    _, c = _residual(cp)
    topo = mesh.topology
    f = only(f for f = 1:n_faces(topo) if topo.faces[f].cells[2] != 0)
    @test c.F_energy[f] ≈ -100 / (0.5 / 1 + 0.5 / 100) rtol = 1e-12
    @test c.F_energy[f] ≈ -198.019801980 rtol = 1e-11
end

# (k, A) and (s k, A/s) are one complete tensor: the residual and every face
# flux are the same, below the conductivity floor of the line as above it.
@testset "the residual reads the complete tensor only: k = $k" for k in (1e-7, 1e-5, 1.0)
    mesh = box_mesh(0:0.25:1.0, 0:0.25:1.0; transform = _shear2(0.35))
    A(x) = x[1] - 0.35x[2] < 0.5 ? _Ra : _Rb
    interior = [f for f = 1:n_faces(mesh.topology) if mesh.topology.faces[f].cells[2] != 0]
    du_ref, c_ref = _residual(_interface_setup(mesh, A; k, thermal = _convective))
    F_ref = copy(c_ref.F_energy)
    for s in (1e-2, 1e2)
        du, c = _residual(
            _interface_setup(mesh, x -> A(x) / s; k = s * k, thermal = _convective),
        )
        @test maximum(abs, du - du_ref) <= 1e-12 * maximum(abs, du_ref)
        @test maximum(abs, c.F_energy[interior] - F_ref[interior]) <=
              1e-12 * maximum(abs, F_ref[interior])
        @test maximum(abs, c.F_energy - F_ref) <= 1e-12 * maximum(abs, F_ref)
    end
end

# Central differences of the residual over every column.
function _central_differences(cp, u; t = 0.0)
    J = zeros(length(u), length(u))
    du_p, du_m, v = similar(u), similar(u), copy(u)
    for j in eachindex(u)
        h = max(1e-7 * abs(u[j]), 1e-9)
        v[j] = u[j] + h
        residual!(du_p, v, cp.def, t, cp.ws)
        v[j] = u[j] - h
        residual!(du_m, v, cp.def, t, cp.ws)
        J[:, j] .= (du_p .- du_m) ./ (2h)
        v[j] = u[j]
    end
    return J
end

# The dense Jacobian holds its nonzeros inside the pattern, equals the colored
# Jacobian per entry and central differences to their truncation error.
function _derivatives(cp; t = 0.0)
    u = cp.u0
    Jd = dense_jacobian(cp.def, cp.ws, u, t)
    J = jacobian!(jac_prototype(cp.jacobian), cp.jacobian, u, t)
    @test all(isfinite, Jd)
    @test isempty(missing_entries(Jd, cp.jacobian.pattern))
    @test compare_entries(Matrix(J), Jd; rtol = 1e-12).passed
    @test compare_entries(Jd, _central_differences(cp, u; t); rtol = 1e-6).passed
    return Jd
end

# Allocations of one residual and one colored Jacobian after warm-up, measured
# inside a function.
function _allocations(cp; t = 0.0)
    u = copy(cp.u0)
    du = similar(u)
    J = jac_prototype(cp.jacobian)
    residual!(du, u, cp.def, t, cp.ws)
    jacobian!(J, cp.jacobian, u, t)
    jacobian!(J, cp.jacobian, u, t)
    return (
        @allocated(residual!(du, u, cp.def, t, cp.ws)),
        @allocated(jacobian!(J, cp.jacobian, u, t))
    )
end

const _kT = T -> 1.0 + 0.002 * (T - 500.0)

@testset "a rotated-tensor jump with k(T), $(name), $(energy) form: derivatives through the images, the pattern, no allocations" for (
        name,
        moving,
    ) in (
        ("fixed mesh", false),
        ("potential closure", true),
    ),
    energy in (:temperature, :enthalpy)

    mesh = box_mesh(0:0.25:1.0, 0:0.25:1.0; transform = _shear2(0.35), fixed = !moving)
    cp = _interface_setup(
        mesh,
        x -> x[1] - 0.35x[2] < 0.5 ? _Ra : 100 * _Rb;
        k = _kT,
        energy,
        (moving ? (use_ale = true, slip = (:left => 1, :right => 1)) : (;))...,
    )
    _derivatives(cp)
    r, j = _allocations(cp)
    @test r == 0
    @test j == 0
end

# The center of 3×3 unit squares in the transition of its regularization:
# its fit's normalized determinant lies between GRADIENT_SPAN and
# GRADIENT_SUPPORT and moves with the temperatures through k(T).
@testset "derivatives through the regularization's transition: $(energy) form" for energy in
                                                                                   (
    :temperature,
    :enthalpy,
)
    Kx = @SMatrix [1.0 2.0; 2.0 5.0]
    Ky = @SMatrix [5.0 2.0; 2.0 1.0]
    K(x) =
        _at(x, SVector(0.5, 1.5)) || _at(x, SVector(2.5, 1.5)) ? Kx :
        _at(x, SVector(1.5, 0.5)) || _at(x, SVector(1.5, 2.5)) ? Ky : SMatrix{2,2}(1.0I)
    cp = _interface_setup(
        box_mesh(0:1.0:3.0, 0:1.0:3.0; fixed = true),
        K;
        k = _kT,
        T = x -> 500.0 + 3x[1] - 2x[2] + 0.5x[1] * x[2],
        energy,
    )
    support = gradient_support(
        find_kernel(ConductionKernel, cp.def.kernels),
        cp.ws.derived,
        cp.def.setup.topology,
    )
    @test support.cells == 1
    @test Pyrolysis.Kernels.GRADIENT_SPAN <
          support.worst <
          Pyrolysis.Kernels.GRADIENT_SUPPORT
    _derivatives(cp)
end

# The right neighbour of the center of 3×3 insulated unit squares conducts `ε`
# times as much: its datum in the center's fit switches to the limit row at
# ε = INSULATING_NEIGHBOUR. Across the switch the Jacobian changes only as
# the neighbour's conductivity does.
@testset "derivatives through the insulating neighbour's limit row" begin
    insulating(ε) = _interface_setup(
        box_mesh(0:1.0:3.0, 0:1.0:3.0; fixed = true),
        x -> _at(x, SVector(2.5, 1.5)) ? ε * _Ra : _Ra;
        thermal = _insulated,
    )
    τ = Pyrolysis.Kernels.INSULATING_NEIGHBOUR
    above, below = insulating(τ * (1 + 1e-6)), insulating(τ * (1 - 1e-6))
    J_above = dense_jacobian(above.def, above.ws, above.u0, 0.0)
    J_below = dense_jacobian(below.def, below.ws, below.u0, 0.0)
    @test compare_entries(J_above, J_below; rtol = 3e-6).passed
    _derivatives(insulating(1e-9))
    _derivatives(insulating(1e-13))
end

# A dual number with one partial, `seed`.
_Seeded(x, seed) = ForwardDiff.Dual{Nothing}(x, seed)

# The center's gradient on 3×3 unit squares under A = I, the cells at the
# temperatures of `_insulating`, with the cell conductivities `k` (plain or
# dual) set before the kernel runs: the right neighbour's normal conductivity
# is then exactly `k[right]/k[center]` times the center's.
function _seeded_gradient(k_center, k_right)
    T = promote_type(typeof(k_center), typeof(k_right))
    mesh = box_mesh(0:1.0:3.0, 0:1.0:3.0)
    s = SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 1.0)
    material = Material(name = :solid, components = (s,), reactions = ())
    topo = mesh.topology
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, reduce(vcat, mesh.positions), 1.0, topo)
    n = n_cells(topo)
    m0 = reshape(1000.0 .* g0.cell_volume, 1, n)
    setup_ = StateSetup(
        topo,
        material,
        CellReference(mesh, m0, material; reltol = 1e-6, threshold = 1e-2),
        nothing,
    )
    layout = StateLayout(1, topo, TemperatureState())
    T0 = [600.0 + 30x[1] - 20x[2] + 5x[1] * x[2] for x in g0.cell_center]
    ds = DerivedState{T}(layout, setup_, copy(T0))
    u = zeros(T, state_length(layout))
    for i = 1:n
        u[mass_index(layout, i, 1)] = m0[1, i]
        u[energy_index(layout, i)] = T0[i]
    end
    derive!(ds, u, layout, setup_)
    center = only(i for i = 1:n if _at(g0.cell_center[i], SVector(1.5, 1.5)))
    right = only(i for i = 1:n if _at(g0.cell_center[i], SVector(2.5, 1.5)))
    ds.k_physical .= one(T)
    ds.k_physical[center], ds.k_physical[right] = k_center, k_right
    kernel = ConductionKernel(topo, g0; tensor = fill(SMatrix{2,2}(1.0I), n))
    c = Contributions{T}(1, 0, n, n_faces(topo), n_nodes(topo), 2)
    contribute!(c, kernel, ds, topo, material, 0.0)
    return cell_gradient(c, Val(2), center)
end

# Every branch of the fit is chosen on the plain values: a derivative seed of
# either sign, at exactly the insulating cutoff or with both cells insulating,
# differentiates the branch the plain residual takes. At the cutoff that is
# the limit row, which does not read the neighbour's conductivity; above it,
# the image does.
@testset "derivative seeds take the plain branch: $(name)" for (name, k_center, k_right) in
                                                               (
    ("at the insulating cutoff", 1.0, Pyrolysis.Kernels.INSULATING_NEIGHBOUR),
    ("both cells insulating", 0.0, 0.0),
)
    plain = _seeded_gradient(k_center, k_right)
    @test all(isfinite, plain)
    for seed in (-1.0, 0.0, 1.0)
        for (kc, kr) in (
            (_Seeded(k_center, 0.0), _Seeded(k_right, seed)),
            (_Seeded(k_center, seed), _Seeded(k_right, 0.0)),
        )
            seeded = _seeded_gradient(kc, kr)
            @test ForwardDiff.value.(seeded) == plain
            @test all(isfinite, ForwardDiff.partials.(seeded)[k][1] for k = 1:2)
        end
        at_cutoff = _seeded_gradient(_Seeded(k_center, 0.0), _Seeded(k_right, seed))
        @test all(iszero, ForwardDiff.partials.(at_cutoff)[k][1] for k = 1:2)
    end
    if k_center > 0
        above = _seeded_gradient(1.0, _Seeded(k_right * (1 + 1e-3), 1.0))
        @test any(!iszero, ForwardDiff.partials.(above)[k][1] for k = 1:2)
    end
end

# The fields saved at two states of a column one cell wide, the initial one
# and a hotter one: the conductivity normal to the top face is the top cell's
# k(T) times n·A n, unfloored, and every cell's gradient is a regularized fit,
# its neighbours lying along one direction only.
@testset "saved outputs read the assembled quantities: the surface conductivity normal to the face, the gradient support" begin
    A = _rot(π / 3) * SMatrix{2,2}(Diagonal(SVector(4.0, 1.0))) * _rot(π / 3)'
    cp = _interface_setup(box_mesh([0.0, 0.25], 0:0.25:1.0; fixed = true), x -> A; k = _kT)
    layout = cp.def.layout
    hot = copy(cp.u0)
    for i = 1:n_cells(layout)
        hot[energy_index(layout, i)] += 200.0
    end
    statistics = SolverStatistics(;
        accepted = 0,
        rejected = 0,
        rhs_evaluations = 0,
        jacobians = 0,
        linear_solves = 0,
        nonlinear_failures = 0,
        domain_rejections = 0,
        largest_constraint_residual = 0.0,
        invocations = 0,
    )
    sol = PyrolysisSolution(
        cp.problem,
        [0.0, 1.0],
        [cp.u0, hot],
        cp.def,
        cp.ws,
        :Success,
        0.0,
        statistics,
        nothing,
    )
    x = sol.extras
    topo = cp.def.setup.topology
    i = topo.faces[first(boundary_faces(topo, :top))].cells[1]
    @test x.surface_conductivity ≈ [_kT(T) * A[2, 2] for T in sol.T[i, :]] rtol = 1e-12
    @test x.regularized_gradients == [4, 4]
    @test all(<(1e-12), abs.(x.gradient_support))
end

# ============================================================================
# A layered manufactured solution
# ============================================================================

# The plane s = n·x − c₀ = 0 through grid nodes splits the domain into K₁ and
# K₂ (k = 1). With b a bubble vanishing on the outline,
#   T = 300 + 100 [b² + s (γ b − 2 b n·K∇b) / (n·K n)]
# is continuous with continuous normal flux 100 γ b across the plane and 300 K
# on the outline; the source is −tr(K H) / (ρc) from its Hessian H at the cell
# centroid.
struct _Layered{D,L,F,G}
    K1::SMatrix{D,D,Float64,L}
    K2::SMatrix{D,D,Float64,L}
    n::SVector{D,Float64}
    c0::Float64
    bubble::F
    gradb::G
end
_layered_tensor(L::_Layered, x) = dot(L.n, x) < L.c0 ? L.K1 : L.K2
function _layered_exact(L::_Layered, x, K)
    b, s = L.bubble(x), dot(L.n, x) - L.c0
    return 300 + 100 * (b^2 + s * (3b - 2b * dot(L.n, K * L.gradb(x))) / dot(L.n, K * L.n))
end
_layered_exact(L::_Layered, x) = _layered_exact(L, x, _layered_tensor(L, x))
function _layered_source(L::_Layered{D}, x) where {D}
    K = _layered_tensor(L, x)
    H = ForwardDiff.hessian(y -> _layered_exact(L, SVector{D}(y), K), Vector(x))
    return -tr(K * H) / 1e6
end

# The bubble Π sin(π ξₖ) of the reference coordinates ξ = B x, and its gradient.
_bubble(B, x) = prod(sin.(π * (B * x)))
function _bubble_gradient(B, x::SVector{D}) where {D}
    s, c = sin.(π * (B * x)), cos.(π * (B * x))
    return B' * SVector(ntuple(k -> π * c[k] * prod(s[m] for m = 1:D if m != k), Val(D)))
end

# The unit square or cube sheared by `s`, its bubble, and the interface plane
# through the nodes of its reference half.
function _sheared_domain(n, s; D = 2)
    axis = range(0.0, 1.0; length = n + 1)
    shear =
        D == 2 ? (@SMatrix [1.0 s; 0.0 1.0]) :
        (@SMatrix [1.0 s 0.3s; 0.0 1.0 0.0; 0.0 0.0 1.0])
    mesh =
        D == 2 ? box_mesh(axis, axis; transform = _shear2(s), fixed = true) :
        box_mesh(axis, axis, axis; transform = _shear3(s), fixed = true)
    B = inv(shear)
    n = normalize(D == 2 ? SVector(1.0, -s) : SVector(1.0, -s, -0.3s))
    return mesh, x -> _bubble(B, x), x -> _bubble_gradient(B, x), n, n[1] * 0.5
end

# The unit square distorted smoothly inside its outline, keeping the plane
# x₁ = 1/2.
function _distorted_domain(n)
    axis = range(0.0, 1.0; length = n + 1)
    transform(x) = SVector(
        x[1] + 0.04sin(2π * x[1]) * sin(2π * x[2]),
        x[2] + 0.04sin(π * x[2]) * sin(2π * x[1]),
    )
    bubble(x) = sin(π * x[1]) * sin(π * x[2])
    gradb(x) = SVector(π * cos(π * x[1]) * sin(π * x[2]), π * sin(π * x[1]) * cos(π * x[2]))
    mesh = box_mesh(axis, axis; transform, fixed = true)
    return mesh, bubble, gradb, SVector(1.0, 0.0), 0.5
end

# The cell-temperature L² and L∞ errors of the steady state, held at 300 K on
# every face, under the complete tensors K₁ and K₂.
function _layered_error(domain, K1, K2)
    mesh, bubble, gradb, n, c0 = domain
    L = _Layered(K1, K2, n, c0, bubble, gradb)
    tags = Tuple(boundary_tags(mesh.topology))
    held = DirichletThermalBC(300.0)
    problem = PyrolysisProblem(
        mesh = mesh,
        material = Material(
            name = :conductor,
            components = (SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 1.0),),
            reactions = (),
        ),
        bc_set = BoundaryConditionSet(
            thermal = NamedTuple{tags}(Tuple(held for _ in tags)),
        ),
        T_initial = x -> 300.0,
        ξ_initial = [x -> 1000.0],
        tspan = (0.0, 1.0),
    )
    cp = setup(
        problem;
        tensor = x -> _layered_tensor(L, x),
        sources = VolumetricSources(
            nothing,
            (x, t) -> _layered_source(L, ForwardDiff.value.(x)),
        ),
    )
    layout = cp.def.layout
    rows = [energy_index(layout, i) for i = 1:n_cells(layout)]
    u = copy(cp.u0)
    du = zeros(length(u))
    J = jac_prototype(cp.jacobian)
    residual!(du, u, cp.def, 0.0, cp.ws)
    scale = maximum(abs, du[rows])
    for _ = 1:2
        residual!(du, u, cp.def, 0.0, cp.ws)
        jacobian!(J, cp.jacobian, u, 0.0)
        u[rows] .-= J[rows, rows] \ du[rows]
    end
    residual!(du, u, cp.def, 0.0, cp.ws)
    @test maximum(abs, du[rows]) <= 1e-10 * scale
    g = cp.ws.derived.geometry
    e = [u[rows[i]] - _layered_exact(L, g.cell_center[i]) for i = 1:n_cells(layout)]
    return sqrt(sum(e .^ 2 .* g.cell_volume) / sum(g.cell_volume)), maximum(abs, e)
end

@testset "a layered manufactured solution: second order with $(name) on $(family)" for (
    family,
    domain,
    name,
    K1,
    K2,
    sizes,
) in (
    (
        "sheared quadrilaterals",
        n -> _sheared_domain(n, 0.35),
        "a 100× scalar jump",
        _Ra,
        100 * _Ra,
        (8, 16, 32, 64),
    ),
    (
        "sheared quadrilaterals",
        n -> _sheared_domain(n, 0.35),
        "a rotated-tensor jump",
        _Ra,
        100 * _Rb,
        (8, 16, 32, 64),
    ),
    (
        "distorted quadrilaterals",
        _distorted_domain,
        "a 100× scalar jump",
        _Ra,
        100 * _Ra,
        (8, 16, 32, 64),
    ),
    (
        "distorted quadrilaterals",
        _distorted_domain,
        "a rotated-tensor jump",
        _Ra,
        100 * _Rb,
        (8, 16, 32, 64, 128),
    ),
    (
        "orthogonal quadrilaterals",
        n -> _sheared_domain(n, 0.0),
        "I against diag(1, 100), aligned with the grid",
        SMatrix{2,2}(1.0I),
        SMatrix{2,2}(Diagonal(SVector(1.0, 100.0))),
        (8, 16, 32, 64, 128, 256),
    ),
    (
        "sheared hexahedra",
        n -> _sheared_domain(n, 0.35; D = 3),
        "a 100× scalar jump",
        _A1,
        100 * _A1,
        (4, 8, 16, 32),
    ),
    (
        "sheared hexahedra",
        n -> _sheared_domain(n, 0.35; D = 3),
        "a rotated-tensor jump",
        _A1,
        100 * _A2,
        (4, 8, 16, 32),
    ),
)
    errors = [_layered_error(domain(n), K1, K2) for n in sizes]
    l2, linf = first.(errors), last.(errors)
    orders = log2.(l2[1:(end-1)] ./ l2[2:end])
    finest = log2(linf[end-1] / linf[end])
    @info "layered manufactured solution, $(name) on $(family), $(join(sizes, ", ")) cells a side" l2 =
        join(l2, ", ") orders = join(orders, ", ") finest_linf_order = finest
    @test all(>=(1.9), orders)
    @test finest >= 1.9
end
