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
# tensor `A`, the contributions after it ran, and the thermal boundary kernel
# sharing its tensors.
function _conduction(mesh, T, A; k = 2.0)
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
    tensor = fill(A, n)
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
