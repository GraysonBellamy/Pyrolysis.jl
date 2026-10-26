using Test
using LinearAlgebra
using StaticArrays
using Pyrolysis
using Pyrolysis.Meshes
using Pyrolysis.Meshes: MeshGeometry
using Pyrolysis.States
using Pyrolysis.States: StateLayout, state_length
using Pyrolysis.Kernels: Contributions
using Pyrolysis.Closures
using Pyrolysis.Closures: PotentialClosure, consistent_multipliers, multiplier_rows!

# A single inert solid on a moving `mesh`: its derived state at the mesh's
# positions and empty contributions, the closure's unknowns `unknowns`.
function _state(mesh; unknowns = :velocities)
    s = SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 1.0)
    material = Material(name = :solid, components = (s,), reactions = ())
    topo = mesh.topology
    D = length(first(mesh.positions))
    x = reduce(vcat, mesh.positions)
    g0 = MeshGeometry{Float64}(topo)
    geometry!(g0, x, topo.cross_section, topo)
    n = n_cells(topo)
    m0 = reshape(1000.0 .* g0.cell_volume, 1, n)
    reference = CellReference(mesh, m0, material; reltol = 1e-6, threshold = 1e-2)
    setup = StateSetup(topo, material, reference, nothing)
    layout = StateLayout(1, topo, TemperatureState(); moving = true, unknowns)
    ds = DerivedState{Float64}(layout, setup, fill(300.0, n))
    u = zeros(state_length(layout))
    for i = 1:n
        u[mass_index(layout, i, 1)] = m0[1, i]
        u[energy_index(layout, i)] = 300.0
    end
    for nd = 1:n_nodes(topo), k = 1:D
        u[position_index(layout, nd, k)] = x[D*(nd-1)+k]
    end
    derive!(ds, u, layout, setup)
    c = Contributions{Float64}(1, 0, n, n_faces(topo), n_nodes(topo), D)
    return (; ds, c, setup, topo, D, n)
end

# The potential closure's velocities and balance residuals at strain rates θ.
function _potential(mesh, θ; slip = (), weighted = true)
    s = _state(mesh; unknowns = :multipliers)
    cl = PotentialClosure(s.topo; slip, weighted)
    s.c.θ .= θ
    s.ds.multiplier .= consistent_multipliers(cl, s.c, s.ds, s.topo)
    correction!(s.ds, cl, s.setup)
    multiplier_rows!(s.c, cl, s.ds, s.topo)
    return (; s..., cl, w = reshape(copy(s.ds.w), s.D, :))
end

@testset "two unit squares: the lumped weights keep a uniform expansion uniform" begin
    mesh = box_mesh(0:1.0:2.0, 0:1.0)
    # Bottom pinned, lateral motion suppressed: only the top nodes rise.
    slip = (:left => 1, :right => 1, :top => 1)
    r = _potential(mesh, 1.0; slip)
    @test r.w[2, 4:6] ≈ [1.0, 1.0, 1.0] atol = 1e-12
    @test norm(r.w[1, :]) == 0
    @test norm(r.w[:, 1:3]) == 0
    @test maximum(abs, r.c.r_multiplier) < 1e-14
    bad = _potential(mesh, 1.0; slip, weighted = false)
    @test bad.w[2, 4:6] ≈ [2 / 3, 4 / 3, 2 / 3] atol = 1e-12
    @test maximum(abs, bad.c.r_multiplier) < 1e-14
end

# Bottom pinned, the side walls slipping, the top free: the potential flow of
# a uniform expansion is the axial `θ z`.
@testset "uniform motion on a $(D)D block between slip walls" for (D, mesh, slip) in (
    (2, box_mesh(0:0.25:1.0, 0:0.2:0.6), [:left => 1, :right => 1]),
    (
        3,
        box_mesh(0:0.5:1.0, 0:0.5:1.0, 0:0.5:1.0),
        [:left => 1, :right => 1, :front => 2, :back => 2],
    ),
)
    θ = 0.7
    r = _potential(mesh, θ; slip)
    z = [p[D] for p in mesh.positions]
    @test r.w[D, :] ≈ θ .* z atol = 1e-12
    @test norm(r.w[1:(D-1), :]) < 1e-12
    @test maximum(abs, r.c.r_multiplier) < 1e-13
    @test rank_check(r.cl, r.topo, r.ds.geometry)
end

@testset "the rank check: a fully pinned body cannot expand" begin
    mesh = box_mesh(0:0.5:1.0, 0:0.5:1.0; pinned = (:bottom, :top, :left, :right))
    s = _state(mesh; unknowns = :multipliers)
    cl = PotentialClosure(s.topo)
    @test !rank_check(cl, s.topo, s.ds.geometry)
    s.c.θ .= 1.0
    @test_throws ArgumentError consistent_multipliers(cl, s.c, s.ds, s.topo)
    open = box_mesh(0:0.5:1.0, 0:0.5:1.0)
    so = _state(open; unknowns = :multipliers)
    @test rank_check(PotentialClosure(so.topo), so.topo, so.ds.geometry)
    @test_throws ArgumentError PotentialClosure(so.topo; slip = (:left => 3,))
    @test_throws ArgumentError PotentialClosure(
        box_mesh(0:1.0, 0:1.0; fixed = true).topology,
    )
end

@testset "on a line the multipliers reproduce the balance rows: $(nameof(typeof(coords)))" for coords in
                                                                                               (
    Planar(),
    Cylindrical(),
    Spherical(),
)
    x = [0.1, 0.25, 0.45, 0.7, 1.0]
    mesh = line_mesh(x; coordinates = coords)
    θ = [0.3, -0.1, 0.5, 0.2]
    r = _potential(mesh, θ)
    s = _state(mesh)
    cl = BalanceRowClosure(s.topo)
    s.c.θ .= θ
    w = consistent_velocities!(zeros(n_nodes(s.topo)), cl, s.c, s.ds, s.topo, 0.0)
    @test vec(r.w) ≈ w rtol = 1e-12
    @test maximum(abs, r.c.r_multiplier) < 1e-14
end

@testset "a column of more than one dimension: the balance rows equal the line's: $(D)D" for D in
                                                                                             (
    2,
    3,
)
    z = [0.0, 0.1, 0.25, 0.45, 0.7, 1.0]
    width = 0.3
    mesh = D == 2 ? box_mesh([0.0, width], z) : box_mesh([0.0, width], [0.0, 1.0], z)
    θ = [0.3, -0.1, 0.5, 0.2, 0.4]
    s = _state(mesh)
    cl = BalanceRowClosure(s.topo)
    @test rank_check(cl, s.topo, s.ds.geometry)
    @test cl.order == collect(1:5)
    s.c.θ .= θ
    w = consistent_velocities!(zeros(D * n_nodes(s.topo)), cl, s.c, s.ds, s.topo, 0.0)
    line = _state(line_mesh(z; cross_section = D == 2 ? width : width * 1.0))
    lc = BalanceRowClosure(line.topo)
    line.c.θ .= θ
    w1 = consistent_velocities!(zeros(6), lc, line.c, line.ds, line.topo, 0.0)
    W = reshape(w, D, :)
    @test norm(W[1:(D-1), :]) == 0
    for (nd, p) in enumerate(mesh.positions)
        level = findfirst(==(p[D]), z)
        @test W[D, nd] ≈ w1[level] rtol = 1e-13 atol = 1e-16
    end
    # The rows vanish at the consistent velocities.
    s.ds.w .= w
    Closures._velocity_rows!(s.c, cl, cl.correction, s.ds, s.topo, 0.0)
    @test maximum(abs, s.c.r_velocity) < 1e-15
    # A mesh of two columns is refused.
    @test_throws ArgumentError BalanceRowClosure(box_mesh(0:0.5:1.0, z).topology)
end
