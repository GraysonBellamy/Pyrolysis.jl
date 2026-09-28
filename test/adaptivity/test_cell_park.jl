"""
Parking a depleted cell: `park_cell!` (mesh) and `park!` (mesh + state)
conserve component mass, condensed enthalpy, χ-weighted dry mass and live
geometry; the live range stays `1:n_live`, face `n_live+1` is the exposed
boundary, ghosts carry zero concentration and (ALE) zero thickness.
"""

using Test
using Pyrolysis
using Pyrolysis.Materials: T_REF
using Pyrolysis.Adaptivity: park_cell!
using Pyrolysis.Solver:
    DepletionCallbackState,
    park!,
    depletion_candidate,
    interior_void_cell,
    find_thin_cells,
    n_live_cells
using Pyrolysis.Properties: MIXING_BLEND_LO, matrix_volume_fraction
using Pyrolysis.Problem: update_mesh_from_state!, get_node_positions, z_offset, chi_range

include(joinpath(@__DIR__, "..", "helpers.jl"))

# Make the state non-trivial: a temperature ramp, partial conversion and some
# gas, then mirror it into the mesh (the callback runs after StepUpdateCallback,
# so the mesh always mirrors u when a park happens).
function _prepared_state(; n_cells = 6, use_ale = true, lateral = false)
    problem = fixture_problem(n_cells = n_cells, lateral = lateral)
    def, ws, u, _ = fixture_state(problem; use_ale = use_ale)
    layout = ws.layout
    n = layout.n_cells
    for i = 1:n
        u[i] = 300.0 + 80.0 * i                    # T
        u[n+i] = 1000.0 * (1.0 - 0.1 * i)          # solid
        u[2n+i] = 2.0 * i                          # gas
    end
    rng = chi_range(layout)
    if !isempty(rng)
        for i = 1:n
            u[first(rng)+i-1] = 0.05 * i
        end
    end
    update_mesh_from_state!(ws.mesh, u, layout)
    return problem, def, ws, u, layout
end

# Extensive invariants over the live cells (fixture: c_solid = c_gas = 1000).
function _invariants(u, layout, mesh, material)
    n_live = mesh.n_active_cells
    n = layout.n_cells
    V = mesh.cell_volumes
    m_solid = sum(u[n+i] * V[i] for i = 1:n_live)
    m_gas = sum(u[2n+i] * V[i] for i = 1:n_live)
    H = sum(V[i] * u[n+i] * 1000.0 * (u[i] - T_REF) for i = 1:n_live)
    rng = chi_range(layout)
    χm =
        isempty(rng) ? 0.0 :
        sum(u[first(rng)+i-1] * mesh.m_dry_solid_initial[i] for i = 1:n_live)
    V_live = sum(V[i] for i = 1:n_live)
    return (; m_solid, m_gas, H, χm, V_live)
end

function _check_mesh_invariants(mesh, layout, u, n_live)
    n = layout.n_cells
    @test mesh.n_active_cells == n_live
    @test all(mesh.active_cell_mask[1:n_live])
    @test !any(mesh.active_cell_mask[(n_live+1):n])
    @test mesh.boundary_faces[1] == [n_live + 1]
    top = mesh.faces[n_live+1]
    @test top.cell_ids == (n_live, 0) && top.is_boundary && top.boundary_id == 1
    @test count(f -> f.boundary_id == 1, mesh.faces) == 1
    @test count(nd -> nd.boundary_id == 1, mesh.nodes) == 1
    @test mesh.nodes[n_live+1].boundary_id == 1
    for i = 1:n_live
        @test mesh.faces[i+1].cell_ids[1] == i
        @test mesh.cells[i].node_ids == (i, i + 1)
    end
    # live state mirrors the mesh; ghosts are inert
    for i = 1:n
        @test u[i] == mesh.cell_states[i].T
        @test u[n+i] == mesh.cell_states[i].ξ[1]
        @test u[2n+i] == mesh.cell_states[i].ξ[2]
    end
    for i = (n_live+1):n
        @test u[n+i] == 0.0 && u[2n+i] == 0.0
        @test mesh.m_dry_solid_initial[i] == 0.0
    end
end

@testset "park_cell! / park! — ALE" begin
    problem, def, ws, u, layout = _prepared_state()
    mesh = ws.mesh
    material = problem.material
    n = layout.n_cells
    cb = DepletionCallbackState(layout, mesh, material, u)
    @test cb.n_live == n
    @test cb.reference ≈ fill(0.01 / n, n)

    inv0 = _invariants(u, layout, mesh, material)
    z0 = copy(get_node_positions(u, layout))
    L0 = z0[end] - z0[1]

    @testset "surface cell" begin
        n_live = park!(u, layout, mesh, material, cb, n, 1.5)
        @test n_live == n - 1
        _check_mesh_invariants(mesh, layout, u, n_live)
        inv = _invariants(u, layout, mesh, material)
        @test inv.m_solid ≈ inv0.m_solid rtol = 1e-13
        @test inv.m_gas ≈ inv0.m_gas rtol = 1e-13
        @test inv.H ≈ inv0.H rtol = 1e-10       # T-merge Newton stops at |ΔT| < 1e-9 K
        @test inv.V_live ≈ inv0.V_live rtol = 1e-13
        z = get_node_positions(u, layout)
        @test z[n_live+1] - z[1] ≈ L0
        @test z[n+1] == z[n_live+1]                  # ghost node on the surface
        @test mesh.cell_volumes[n] == 0.0
        @test issorted(z[1:(n_live+1)])
        @test all(diff(z[1:(n_live+1)]) .> 0)
        @test z[1:(n_live)] == z0[1:(n_live)]        # nodes below the merge untouched
        @test cb.park_times == [1.5] && cb.parked_cells == [n]
        @test n_live_cells(cb, 1.5) == n             # saved state at t_park is pre-park
        @test n_live_cells(cb, 1.5 + eps(1.5)) == n - 1
    end

    @testset "interior cell shifts the cells above it" begin
        n_live0 = mesh.n_active_cells
        T_above = [u[i] for i = 4:n_live0]
        z_before = copy(get_node_positions(u, layout))
        inv_a = _invariants(u, layout, mesh, material)
        n_live = park!(u, layout, mesh, material, cb, 3, 2.0)
        @test n_live == n_live0 - 1
        _check_mesh_invariants(mesh, layout, u, n_live)
        inv = _invariants(u, layout, mesh, material)
        @test inv.m_solid ≈ inv_a.m_solid rtol = 1e-13
        @test inv.m_gas ≈ inv_a.m_gas rtol = 1e-13
        @test inv.H ≈ inv_a.H rtol = 1e-10
        @test inv.V_live ≈ inv_a.V_live rtol = 1e-13
        # cells 4.. moved down one slot, keeper 2 spans old nodes 2..4
        @test [u[i] for i = 3:n_live] == T_above
        z = get_node_positions(u, layout)
        @test z[1:2] == z_before[1:2]
        @test z[3:(n_live+1)] == z_before[4:(n_live0+1)]
        @test mesh.cell_volumes[2] ≈ (z_before[4] - z_before[2]) * mesh.cross_section_area
        @test z[n_live+1] - z[1] ≈ L0
        @test cb.parked_cells == [n, 3]
        @test cb.reference[1:n_live] ≈ fill(0.01 / n, n_live)
    end

    @testset "bottom cell parks into the cell above" begin
        n_live0 = mesh.n_active_cells
        inv_a = _invariants(u, layout, mesh, material)
        n_live = park!(u, layout, mesh, material, cb, 1, 2.5)
        @test n_live == n_live0 - 1
        _check_mesh_invariants(mesh, layout, u, n_live)
        inv = _invariants(u, layout, mesh, material)
        @test inv.m_solid ≈ inv_a.m_solid rtol = 1e-13
        @test inv.H ≈ inv_a.H rtol = 1e-10
        z = get_node_positions(u, layout)
        @test z[1] == 0.0
        @test z[n_live+1] - z[1] ≈ L0
    end

    @testset "termination floor" begin
        @test_throws ErrorException park_cell!(
            create_uniform_mesh(0.01, 1, Val(2)),
            1,
            2;
            material = material,
            ale = true,
        )
    end
end

@testset "park! — WithChi conserves Σ χᵢ·m_dry_init,i" begin
    problem, def, ws, u, layout = _prepared_state(lateral = true)
    mesh = ws.mesh
    material = problem.material
    n = layout.n_cells
    @test !isempty(chi_range(layout))
    cb = DepletionCallbackState(layout, mesh, material, u)
    inv0 = _invariants(u, layout, mesh, material)
    M0 = sum(mesh.m_dry_solid_initial)
    park!(u, layout, mesh, material, cb, 4, 0.1)
    park!(u, layout, mesh, material, cb, mesh.n_active_cells, 0.2)
    inv = _invariants(u, layout, mesh, material)
    @test inv.χm ≈ inv0.χm rtol = 1e-12
    @test inv.m_solid ≈ inv0.m_solid rtol = 1e-13
    @test sum(mesh.m_dry_solid_initial) ≈ M0 rtol = 1e-13
    rng = chi_range(layout)
    @test u[first(rng)+n-1] == 0.0 && u[first(rng)+n-2] == 0.0
end

@testset "park_cell! / park! — fixed mesh" begin
    problem, def, ws, u, layout = _prepared_state(use_ale = false)
    mesh = ws.mesh
    material = problem.material
    n = layout.n_cells
    cb = DepletionCallbackState(layout, mesh, material, u)
    @test cb.reference ≈ [u[n+i] for i = 1:n]        # initial condensed concentration
    V0 = copy(mesh.cell_volumes)
    inv0 = _invariants(u, layout, mesh, material)

    @test_throws ErrorException park_cell!(mesh, 3, 2; material = material, ale = false)

    n_live = park!(u, layout, mesh, material, cb, n, 3.0)
    @test n_live == n - 1
    _check_mesh_invariants(mesh, layout, u, n_live)
    @test mesh.cell_volumes == V0                     # nothing moves
    @test mesh.cell_centers[n] == 0.5 * (mesh.nodes[n].z + mesh.nodes[n+1].z)
    inv = _invariants(u, layout, mesh, material)
    @test inv.m_solid ≈ inv0.m_solid rtol = 1e-13
    @test inv.m_gas ≈ inv0.m_gas rtol = 1e-13
    @test inv.H ≈ inv0.H rtol = 1e-10
    @test cb.reference[1:n_live] ≈ [1000.0 * (1.0 - 0.1 * i) for i = 1:n_live]
end

@testset "depletion_candidate and find_thin_cells on the live range" begin
    problem, def, ws, u, layout = _prepared_state()
    mesh = ws.mesh
    material = problem.material
    n = layout.n_cells
    cb = DepletionCallbackState(layout, mesh, material, u)
    @test depletion_candidate(u, layout, mesh, material, cb, 0.05) == 0
    zoff = z_offset(layout)
    Δz = 0.01 / n
    # cell 2 at 4% of its thickness, cell 5 at 3%: surface-most first
    u_thin = copy(u)
    u_thin[zoff+3] = u_thin[zoff+2] + 0.04 * Δz
    u_thin[zoff+6] = u_thin[zoff+5] + 0.03 * Δz
    @test find_thin_cells(u_thin, layout, cb.reference, 0.05, n) == [5, 2]
    @test depletion_candidate(u_thin, layout, mesh, material, cb, 0.05) == 5
    # ghosts are never scanned
    @test find_thin_cells(u_thin, layout, cb.reference, 0.05, 4) == [2]
    @test_throws DimensionMismatch find_thin_cells(
        u,
        layout,
        cb.reference[1:(n-1)],
        0.05,
        n,
    )
    u_bad = copy(u)
    u_bad[zoff+3] = u_bad[zoff+4] + 5e-4
    @test_throws ErrorException find_thin_cells(u_bad, layout, cb.reference, 0.05, n)

    # fixed mesh: only the top live cell, only below the threshold
    problem_e, _, ws_e, u_e, layout_e = _prepared_state(use_ale = false)
    cb_e = DepletionCallbackState(layout_e, ws_e.mesh, problem_e.material, u_e)
    @test depletion_candidate(u_e, layout_e, ws_e.mesh, problem_e.material, cb_e, 0.05) == 0
    u_e[n+n] = 0.04 * cb_e.reference[n]
    @test depletion_candidate(u_e, layout_e, ws_e.mesh, problem_e.material, cb_e, 0.05) == n
    u_e[n+n-1] = 0.0                                   # interior emptied cell: not a candidate
    @test depletion_candidate(u_e, layout_e, ws_e.mesh, problem_e.material, cb_e, 0.05) == n
end

@testset "interior_void_cell: a void below the top live cell (fixed mesh)" begin
    problem, _, ws, u, layout = _prepared_state(use_ale = false)
    mesh = ws.mesh
    material = problem.material
    n = layout.n_cells
    cb = DepletionCallbackState(layout, mesh, material, u)
    @test interior_void_cell(u, layout, material, cb) == 0

    # Cell 3 holds a trace of solid and no gas (the fixture's gas has γ = 1 and
    # counts as matrix volume): φ_m = ξ_s/ρ_s, just below the void edge.
    φ_m(v, i) = matrix_volume_fraction((v[n+i], v[2n+i]), v[i], material)
    ρ_s = 1000.0
    u_void = copy(u)
    u_void[2n+3] = 0.0
    u_void[n+3] = 0.5 * MIXING_BLEND_LO * ρ_s
    @test φ_m(u_void, 3) ≈ 0.5 * MIXING_BLEND_LO
    @test interior_void_cell(u_void, layout, material, cb) == 3

    # Twice the edge: a thin porous cell, not a void.
    u_thin = copy(u_void)
    u_thin[n+3] = 2.0 * MIXING_BLEND_LO * ρ_s
    @test interior_void_cell(u_thin, layout, material, cb) == 0

    # The lowest void is reported.
    u_two = copy(u_void)
    u_two[2n+5] = 0.0
    u_two[n+5] = 0.0
    @test interior_void_cell(u_two, layout, material, cb) == 3

    # An empty top live cell is the park path, not an interior void.
    u_top = copy(u)
    u_top[2n+n] = 0.0
    u_top[n+n] = 0.0
    @test interior_void_cell(u_top, layout, material, cb) == 0

    # A cell that started with no condensed mass is not tested.
    u_gap = copy(u_void)
    u_gap[n+3] = 0.0
    cb_gap = DepletionCallbackState(layout, mesh, material, u_gap)
    @test cb_gap.reference[3] == 0.0
    @test interior_void_cell(u_gap, layout, material, cb_gap) == 0

    # After a park the new top live cell is excluded and ghosts are not scanned.
    u_park = copy(u)
    u_park[n+n] = 0.04 * cb.reference[n]
    update_mesh_from_state!(mesh, u_park, layout)
    park!(u_park, layout, mesh, material, cb, n, 1.0)
    @test cb.n_live == n - 1
    u_park[2n+n-1] = 0.0
    u_park[n+n-1] = 0.0
    @test interior_void_cell(u_park, layout, material, cb) == 0
    u_park[2n+n-2] = 0.0
    u_park[n+n-2] = 0.0
    @test interior_void_cell(u_park, layout, material, cb) == n - 2
end
