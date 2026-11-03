using Test
using ForwardDiff
using LinearAlgebra
using StaticArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.Meshes: MeshGeometry

# Every node's coordinates in turn: the flat position vector of a mesh.
_flat(mesh) = reduce(vcat, mesh.positions)

function _geometry(mesh, x = _flat(mesh))
    topo = mesh.topology
    g = MeshGeometry{eltype(x)}(topo)
    geometry!(g, x, topo.cross_section, topo)
    return g
end

# The sign of face `f` in the balance of cell `i`: +1 where `i` owns it.
_σ(topo, f, i) = topo.faces[f].cells[1] == i ? 1 : -1

# Meshes whose cells are not rectangles: a uniform shear (every cell the same
# parallelogram or parallelepiped), a smooth distortion that leaves no two
# cells alike, and in 3D one displaced corner, which bends its faces.
_shear(x::SVector{2}) = SVector(x[1] + 0.3x[2], x[2])
_shear(x::SVector{3}) = SVector(x[1] + 0.3x[3], x[2] - 0.2x[3], x[3])
_wave(x::SVector{2}) =
    SVector(x[1] + 0.08sin(π * x[1]) * sin(2π * x[2]), x[2] + 0.06sin(2π * x[1]) * x[2])
_wave(x::SVector{3}) = SVector(
    x[1] + 0.08sin(π * x[1]) * sin(2π * x[3]),
    x[2] + 0.05sin(π * x[2]) * x[3],
    x[3] + 0.06sin(2π * x[1]) * x[2] * x[3],
)

_meshes_2d() = (
    box = box_mesh(0:0.25:1.0, [0.0, 0.2, 0.5, 1.0]),
    shear = box_mesh(0:0.25:1.0, 0:0.25:1.0; transform = _shear),
    wave = box_mesh(0:0.2:1.0, 0:0.25:1.0; transform = _wave),
)

function _meshes_3d()
    bent = box_mesh(0:0.5:1.0, 0:0.5:1.0, 0:0.5:1.0; transform = _shear)
    p = copy(bent.positions)
    p[14] += SVector(0.1, -0.07, 0.05)
    cells = [bent.topology.cells[i].nodes for i = 1:n_cells(bent.topology)]
    return (
        box = box_mesh(0:0.5:1.0, [0.0, 0.3, 1.0], 0:0.25:1.0),
        shear = box_mesh(0:0.5:1.0, 0:0.5:1.0, 0:0.5:1.0; transform = _shear),
        wave = box_mesh(0:0.25:1.0, 0:0.5:1.0, 0:0.25:1.0; transform = _wave),
        bent = cell_mesh(p, cells; tag = (i, k, s) -> (:a, :b, :c)[k], pinned = ()),
    )
end

@testset "box meshes: connectivity, tags and pinning" begin
    m = box_mesh(0:0.25:1.0, 0:0.5:1.0)
    topo = m.topology
    @test n_cells(topo) == 8
    @test n_nodes(topo) == 15
    @test n_faces(topo) == 5 * 2 + 4 * 3
    @test boundary_tags(topo) == [:bottom, :left, :right, :top]
    @test length.(boundary_faces.(Ref(topo), [:bottom, :left, :right, :top])) ==
          [4, 2, 2, 4]
    @test all(f -> topo.faces[f].cells[2] == 0, boundary_faces(topo, :top))
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        @test (nb == 0) == (topo.faces[f].tag !== :interior)
        nb == 0 || @test o < nb
    end
    # Every cell lists its faces as −x, +x, −y, +y and its nodes in tensor order.
    c = topo.cells[6]
    @test c.nodes == (7, 8, 12, 13)
    @test [m.positions[n] for n in c.nodes] == [SVector(0.25, 0.5), SVector(0.5, 0.5), SVector(0.25, 1.0), SVector(0.5, 1.0)]
    @test topo.pinned == [1, 2, 3, 4, 5]
    @test !is_fixed(topo)
    @test is_fixed(box_mesh(0:1.0, 0:1.0; fixed = true).topology)
    @test box_mesh(0:1.0, 0:1.0; pinned = (:left, :bottom)).topology.pinned == [1, 2, 3]
    @test box_mesh(0:1.0, 0:1.0; pinned = ()).topology.pinned == Int[]

    m3 = box_mesh(0:0.5:1.0, 0:1.0, 0:0.25:1.0)
    t3 = m3.topology
    @test (n_cells(t3), n_nodes(t3)) == (8, 3 * 2 * 5)
    @test n_faces(t3) == 3 * 1 * 4 + 2 * 2 * 4 + 2 * 1 * 5
    @test boundary_tags(t3) == [:back, :bottom, :front, :left, :right, :top]
    @test length(boundary_faces(t3, :top)) == 2
    @test length(t3.pinned) == 6

    @test_throws ArgumentError box_mesh(0:1.0, [0.0, 0.5, 0.5])
    @test_throws ArgumentError box_mesh(0:1.0, [0.0])
    @test_throws ArgumentError box_mesh(0:1.0, 0:1.0; pinned = (:side,))
    @test_throws ArgumentError box_mesh(0:1.0, 0:1.0; transform = x -> SVector(-x[1], x[2]))
    square = [SVector(0.0, 0.0), SVector(1.0, 0.0), SVector(0.0, 1.0), SVector(1.0, 1.0)]
    @test_throws ArgumentError cell_mesh(square, [(1, 2, 3, 3)]; tag = (i, k, s) -> :b)
    @test_throws ArgumentError cell_mesh(square, [(1, 2, 3, 5)]; tag = (i, k, s) -> :b)
    @test_throws ArgumentError cell_mesh(square, [(1, 2, 3)]; tag = (i, k, s) -> :b)
    @test_throws ArgumentError cell_mesh(
        square,
        [(1, 2, 3, 4)];
        tag = (i, k, s) -> :interior,
    )
end

@testset "normals point from owner to neighbour, out of the domain: $(D)D $name" for (
        D,
        meshes,
    ) in (
        (2, _meshes_2d()),
        (3, _meshes_3d()),
    ),
    name in keys(meshes)

    mesh = meshes[name]
    topo = mesh.topology
    g = _geometry(mesh)
    @test all(>(0), g.face_area)
    @test all(n -> norm(n) ≈ 1, g.face_normal)
    for f = 1:n_faces(topo)
        o, nb = topo.faces[f].cells
        n = g.face_normal[f]
        @test dot(g.face_center[f] - g.cell_center[o], n) > 0
        nb == 0 || @test dot(g.cell_center[nb] - g.face_center[f], n) > 0
        nb == 0 || @test cell_distance(g, topo, f) > 0
    end
end

@testset "metric invariants: $(D)D $name" for (D, meshes) in
                                              ((2, _meshes_2d()), (3, _meshes_3d())),
    name in keys(meshes)

    mesh = meshes[name]
    topo = mesh.topology
    g = _geometry(mesh)
    for i = 1:n_cells(topo)
        # Every cell closes: Σ_f σ A_f n_f = 0.
        S = sum(
            _σ(topo, f, i) * g.face_area[f] * g.face_normal[f] for f in topo.cells[i].faces
        )
        @test norm(S) <= 1e-15
        # The divergence identity over the sweep vectors against ∫ det J.
        p = ntuple(a -> mesh.positions[topo.cells[i].nodes[a]], 1 << D)
        @test g.cell_volume[i] ≈ cell_volume_quadrature(p) rtol = 1e-13
        @test g.cell_volume[i] > 0
    end
    # A face's sweep vectors sum to its area vector.
    for f = 1:n_faces(topo)
        @test sum(g.face_sweep[:, f]) ≈ g.face_area[f] * g.face_normal[f] atol = 1e-15
    end
    @test vec(sum(g.cell_node_volume; dims = 1)) ≈ g.cell_volume rtol = 1e-14
    @test all(>(0), g.cell_node_volume)

    # ∂Vᵢ/∂x_n by AD equals the sweep operator Σ_f σ s_{n,f}.
    x = _flat(mesh)
    J = ForwardDiff.jacobian(y -> _geometry(mesh, y).cell_volume, x)
    Dop = zeros(n_cells(topo), length(x))
    for i = 1:n_cells(topo), f in topo.cells[i].faces
        for (q, n) in enumerate(topo.faces[f].nodes)
            Dop[i, (D*(n-1)+1):(D*n)] .+= _σ(topo, f, i) .* g.face_sweep[q, f]
        end
    end
    @test J ≈ Dop atol = 1e-14
end

@testset "centroids and lumped volumes of affine cells" begin
    # Two unit squares side by side: a quarter of each square per corner.
    g = _geometry(box_mesh(0:1.0:2.0, 0:1.0))
    @test g.cell_node_volume ≈ fill(0.25, 4, 2)
    @test g.cell_center ≈ [SVector(0.5, 0.5), SVector(1.5, 0.5)]
    # An affine map takes centroids to centroids and scales every volume alike.
    A = @SMatrix [1.2 0.4 0.1; -0.3 0.9 0.2; 0.1 0.0 1.1]
    b = SVector(0.3, -0.1, 0.2)
    ref = box_mesh(0:0.5:1.0, 0:1.0, 0:0.5:1.0)
    map_ = box_mesh(0:0.5:1.0, 0:1.0, 0:0.5:1.0; transform = x -> A * x + b)
    g0, g1 = _geometry(ref), _geometry(map_)
    @test g1.cell_center ≈ [A * c + b for c in g0.cell_center]
    @test g1.cell_volume ≈ det(A) .* g0.cell_volume
    @test g1.cell_node_volume ≈ det(A) .* g0.cell_node_volume
    @test g1.face_center ≈ [A * c + b for c in g0.face_center]
end

@testset "dual-number positions pass through the multi-D metrics" begin
    for mesh in (_meshes_2d().wave, _meshes_3d().bent)
        f(y) =
            (g = _geometry(mesh, y); [g.cell_volume; g.face_area; vec(g.cell_node_volume)])
        J = ForwardDiff.jacobian(f, _flat(mesh))
        @test all(isfinite, J)
    end
end

@testset "1D sweep vectors and lumped volumes" begin
    x = [0.1, 0.25, 0.45, 0.7, 1.0]
    for coords in (Planar(), Cylindrical(), Spherical())
        mesh = line_mesh(x; coordinates = coords)
        topo = mesh.topology
        g = _geometry(mesh, x)
        @test all(
            g.face_sweep[1, f] == g.face_area[f] * g.face_normal[f] for f = 1:n_faces(topo)
        )
        @test vec(sum(g.cell_node_volume; dims = 1)) ≈ g.cell_volume rtol = 1e-14
        coords isa Planar &&
            @test g.cell_node_volume ≈ [g.cell_volume'; g.cell_volume'] ./ 2
    end
end
