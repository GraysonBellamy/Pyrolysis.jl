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
