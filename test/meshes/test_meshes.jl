using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.Meshes: MeshGeometry, MeshTopology

const _coords = (Planar(), Cylindrical(), Spherical())

_radii(n; r0 = 0.0, r1 = 1.0, stretch = 1.1) =
    r0 .+
    (r1 - r0) .* cumsum([0.0; [stretch^k for k = 0:(n-1)]]) ./ sum(stretch^k for k = 0:(n-1))

function _metrics(coords, x; A = nothing)
    mesh =
        A === nothing ? line_mesh(x; coordinates = coords) :
        line_mesh(x; coordinates = coords, cross_section = A)
    topo = mesh.topology
    g = MeshGeometry{eltype(x)}(topo)
    geometry!(g, x, topo.cross_section, topo)
    return topo, g
end

@testset "line meshes: connectivity and tags" begin
    x = _radii(6; r0 = 0.1)
    mesh = line_mesh(x)
    topo = mesh.topology
    @test n_cells(topo) == 6
    @test n_nodes(topo) == 7
    @test n_faces(topo) == 7
    @test boundary_tags(topo) == [:bottom, :top]
    @test boundary_faces(topo, :bottom) == [1]
    @test boundary_faces(topo, :top) == [7]
    @test topo.faces[1].cells == (1, 0)
    @test topo.faces[7].cells == (6, 0)
    @test all(topo.faces[f].cells == (f - 1, f) for f = 2:6)
    @test all(topo.faces[f].tag === :interior for f = 2:6)
    @test topo.pinned == [1]
    @test !is_fixed(topo)
    @test is_fixed(line_mesh(x; fixed = true).topology)
    @test [p[1] for p in mesh.positions] == x
    @test_throws ArgumentError line_mesh([0.0, 1.0, 1.0])
    @test_throws ArgumentError line_mesh([-0.1, 1.0]; coordinates = Cylindrical())
    @test_throws ArgumentError line_mesh([0.0, 1.0]; bottom = :top)
    @test_throws ArgumentError line_mesh([0.0, 1.0]; top = :interior)
    @test_throws ArgumentError boundary_faces(topo, :side)
end

@testset "metrics against the closed-form face areas and volumes" begin
    x = _radii(8; r0 = 0.05)
    for (coords, area, volume) in (
        (Planar(), r -> 2.5, (a, b) -> 2.5 * (b - a)),
        (Cylindrical(), r -> 2π * r * 2.5, (a, b) -> π * (b^2 - a^2) * 2.5),
        (Spherical(), r -> 4π * r^2, (a, b) -> 4π / 3 * (b^3 - a^3)),
    )
        topo, g = coords isa Spherical ? _metrics(coords, x) : _metrics(coords, x; A = 2.5)
        @test all(g.face_area[f] ≈ area(x[f]) for f = 1:n_faces(topo))
        @test all(g.cell_volume[i] ≈ volume(x[i], x[i+1]) for i = 1:n_cells(topo))
        @test all(g.cell_center[i][1] ≈ (x[i] + x[i+1]) / 2 for i = 1:n_cells(topo))
        # V = (1/d) Σ_f A_f (x_f · n_f) over each cell's faces.
        d = geometric_dimension(coords)
        for i = 1:n_cells(topo)
            s = sum(
                g.face_area[f] *
                g.face_center[f][1] *
                g.face_normal[f][1] *
                (topo.faces[f].cells[1] == i ? 1 : -1) for f in topo.cells[i].faces
            )
            @test s / d ≈ g.cell_volume[i] rtol = 1e-13
        end
    end
end

@testset "normals point from owner to neighbour, out of the domain" begin
    topo, g = _metrics(Planar(), _radii(5))
    @test g.face_normal[1][1] == -1
    @test all(g.face_normal[f][1] == 1 for f = 2:n_faces(topo))
    # Planar cells close: Σ_f A_f n_f = 0 with the outward orientation.
    for i = 1:n_cells(topo)
        @test sum(
            g.face_area[f] * g.face_normal[f][1] * (topo.faces[f].cells[1] == i ? 1 : -1)
            for f in topo.cells[i].faces
        ) == 0
    end
    for f = 2:(n_faces(topo)-1)
        o, nb = topo.faces[f].cells
        @test cell_distance(g, topo, f) ≈ g.cell_center[nb][1] - g.cell_center[o][1]
        @test center_to_face_distance(g, topo, f, o) > 0
        @test center_to_face_distance(g, topo, f, nb) > 0
    end
    @test center_to_face_distance(g, topo, 1, 1) > 0
end

@testset "volume derivative by AD equals the face sweep operator" begin
    # dVᵢ/dt = Σ_f A_f (w·n)_f: ∂Vᵢ/∂x_{i+1} = A_{i+1}, ∂Vᵢ/∂x_i = −A_i, all else 0.
    for coords in _coords
        x = _radii(7; r0 = 0.2)
        topo = line_mesh(x; coordinates = coords).topology
        V(y) = begin
            g = MeshGeometry{eltype(y)}(topo)
            geometry!(g, y, topo.cross_section, topo)
            g.cell_volume
        end
        J = ForwardDiff.jacobian(V, x)
        g = MeshGeometry{Float64}(topo)
        geometry!(g, x, topo.cross_section, topo)
        D = zeros(n_cells(topo), n_nodes(topo))
        for i = 1:n_cells(topo), f in topo.cells[i].faces
            σ = topo.faces[f].cells[1] == i ? 1 : -1
            D[i, topo.faces[f].nodes[1]] += σ * g.face_area[f] * g.face_normal[f][1]
        end
        @test J ≈ D rtol = 1e-13 atol = 1e-15
    end
end

@testset "dual-number positions and cross section pass through" begin
    x = _radii(4; r0 = 0.1)
    topo = line_mesh(x; coordinates = Cylindrical()).topology
    f(p) = begin
        y = p[1:(end-1)]
        g = MeshGeometry{eltype(p)}(topo)
        geometry!(g, y, p[end], topo)
        [g.cell_volume; g.face_area]
    end
    J = ForwardDiff.jacobian(f, [x; 2.0])
    @test all(isfinite, J)
    # Volumes and areas are linear in the cross-section measure.
    g = MeshGeometry{Float64}(topo)
    geometry!(g, x, 1.0, topo)
    @test J[:, end] ≈ [g.cell_volume; g.face_area]
end

_nodes(mesh) = [p[1] for p in mesh.positions]

function _volumes(mesh)
    topo = mesh.topology
    g = MeshGeometry{Float64}(topo)
    geometry!(g, _nodes(mesh), topo.cross_section, topo)
    return g.cell_volume
end

@testset "uniform meshes: the nodes, the tags and the volumes" begin
    for (L, n) in ((0.01, 100), (0.023, 7), (1.0, 3)), coords in _coords
        mesh = create_uniform_mesh(L, n; coordinates = coords)
        topo = mesh.topology
        @test mesh isa Mesh{1}
        @test topo.coordinates === coords
        @test _nodes(mesh) == [(i - 1) * (L / n) for i = 1:(n+1)]
        @test boundary_tags(topo) == [:bottom, :top]
        @test boundary_faces(topo, :top) == [n + 1]
        @test topo.pinned == [1] && isempty(topo.lagrangian)
        A = coords isa Spherical ? 4π : 1.0
        @test topo.cross_section == A
        z = _nodes(mesh)
        exact = [shell_volume(coords, z[i], z[i+1], A) for i = 1:n]
        @test _volumes(mesh) == exact
    end
    @test create_uniform_mesh(0.01, 4; cross_section = 2.5).topology.cross_section == 2.5
    @test _volumes(create_uniform_mesh(0.01, 4; cross_section = 2.5)) ≈
          fill(2.5 * 0.0025, 4)
    sphere = create_uniform_mesh(0.01, 4; coordinates = Spherical(), cross_section = π)
    @test sphere.topology.cross_section == Float64(π)
    @test_throws ArgumentError create_uniform_mesh(0.01, 0)
    @test_throws ArgumentError create_uniform_mesh(-0.01, 4)
    @test_throws ArgumentError create_uniform_mesh(0.01, 4; cross_section = 0.0)
    @test_throws ArgumentError create_uniform_mesh(0.01, 4; cross_section = Inf)
end

@testset "stretched meshes: the ratio, the refined side and the ends" begin
    for (L, n, r) in ((0.01, 50, 1.1), (0.02, 9, 1.2), (3.0, 4, 2.0))
        top = _nodes(create_stretched_mesh(L, n, r))
        bottom = _nodes(create_stretched_mesh(L, n, r; stretch_from = :bottom))
        @test top[1] == 0 && top[end] == L
        @test bottom[1] == 0 && bottom[end] == L
        # Each cell is `r` times its neighbour on the refined side.
        h_top, h_bottom = diff(top), diff(bottom)
        @test all(h_top[i] / h_top[i+1] ≈ r for i = 1:(n-1))
        @test all(h_bottom[i+1] / h_bottom[i] ≈ r for i = 1:(n-1))
        @test h_top[end] ≈ L * (r - 1) / (r^n - 1)
        @test h_bottom ≈ reverse(h_top)
    end
    # A unit ratio is the uniform range, reflected for the refined top.
    uniform = collect(range(0.0, 0.01; length = 7))
    @test _nodes(create_stretched_mesh(0.01, 6, 1.0; stretch_from = :bottom)) == uniform
    @test _nodes(create_stretched_mesh(0.01, 6, 1.0)) == 0.01 .- reverse(uniform)
    @test _nodes(create_stretched_mesh(0.01, 6, 1.0 + 1e-12)) == 0.01 .- reverse(uniform)
    for coords in _coords
        mesh = create_stretched_mesh(0.02, 5, 1.3; coordinates = coords)
        @test mesh.topology.coordinates === coords
        z, A = _nodes(mesh), mesh.topology.cross_section
        @test _volumes(mesh) == [shell_volume(coords, z[i], z[i+1], A) for i = 1:5]
    end
    @test_throws ArgumentError create_stretched_mesh(0.01, 6, 0.9)
    @test_throws ArgumentError create_stretched_mesh(0.01, 6, 1.1; stretch_from = :side)
end

@testset "meshes from node lists" begin
    z = [0.0, 0.001, 0.0025, 0.004, 0.01]
    for coords in _coords
        mesh = create_mesh_from_nodes(z; coordinates = coords)
        @test _nodes(mesh) == z
        @test mesh.topology.coordinates === coords
    end
    @test _nodes(create_mesh_from_nodes([0, 1, 3])) == [0.0, 1.0, 3.0]
    @test create_mesh_from_nodes(z; cross_section = 0.3).topology.cross_section == 0.3
    @test_throws ArgumentError create_mesh_from_nodes([0.001, 0.002])
    @test_throws ArgumentError create_mesh_from_nodes([0.0])
    @test_throws ArgumentError create_mesh_from_nodes(Float64[])
    @test_throws ArgumentError create_mesh_from_nodes([0.0, 0.002, 0.001])
    @test_throws ArgumentError create_mesh_from_nodes([0.0, NaN])
end
