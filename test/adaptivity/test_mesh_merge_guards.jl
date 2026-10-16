"""
Guards around the conservative merge and the depletion callback.

- The depletion callback's *condition* keeps firing at `n_live ≤ min_cells`
  (the affect! terminates there).
- The energy-conserving temperature merge conserves the condensed-phase
  sensible enthalpy Σ m_j ∫c_p dT — exact for temperature-dependent c_p,
  and gas heat capacity is excluded.
- Inverted (tangled) cells are detected from signed widths; the abs()-based
  cached volumes can never see them, and `find_thin_cells` refuses to park
  an inverted cell.
"""

using Test
using Pyrolysis
using Pyrolysis.Geometry: CellState, update_geometry!
using Pyrolysis.Adaptivity: merge_temperatures_conservatively
using Pyrolysis.Solver: find_thin_cells, create_depletion_callback, DepletionCallbackState
using Pyrolysis.Problem: z_offset, get_node_positions

include(joinpath(@__DIR__, "..", "helpers.jl"))

# ============================================================================
# Enthalpy-conserving merge (every component)
# ============================================================================
@testset "merge conserves Σm∫c_p dT over every component" begin
    # (a) Constant c_p: reduces to the ρcp-weighted mean — the 720 K worked
    #     example from the docstring/guide (no gas present).
    mat_const = Material(
        name = :const_cp,
        components = (
            SolidComponent(:s, ρ = 1000.0, c = 1000.0, k = 0.2),
            GaseousComponent(:g, M = 0.05, c = 5000.0, k = 0.02, λ = 1e-5),
        ),
        reactions = (),
    )
    a = CellState{2}(800.0, (2.0, 0.0), 1.0)   # ρcp = 2000 J/(m³K)
    b = CellState{2}(400.0, (0.5, 0.0), 1.0)   # ρcp = 500
    T_m = merge_temperatures_conservatively(1.0, a, 1.0, b, mat_const)
    @test T_m ≈ (2000.0 * 800.0 + 500.0 * 400.0) / 2500.0  # = 720 K

    # (b) The pore gas stores sensible heat like every other component:
    #     loading both cells with gas moves the merged temperature to the
    #     all-component ρcp-weighted mean.
    a_gas = CellState{2}(800.0, (2.0, 50.0), 1.0)   # ρcp = 2000 + 250000
    b_gas = CellState{2}(400.0, (0.5, 50.0), 1.0)   # ρcp = 500 + 250000
    T_gas = (252000.0 * 800.0 + 250500.0 * 400.0) / (252000.0 + 250500.0)
    @test merge_temperatures_conservatively(1.0, a_gas, 1.0, b_gas, mat_const) ≈ T_gas

    # (c) Linear c_p(T) = a + bT: the conserved functional is
    #     h(T) = a(T−T₀) + b/2(T²−T₀²); pre-fix the product form c_p(T)·T
    #     carried the b-term with coefficient 1 instead of ½.
    cp_a, cp_b = 500.0, 2.0
    mat_lin = Material(
        name = :lin_cp,
        components = (
            SolidComponent(:s, ρ = 1000.0, c = T -> cp_a + cp_b * T, k = 0.2),
            GaseousComponent(:g, M = 0.05, c = 1000.0, k = 0.02, λ = 1e-5),
        ),
        reactions = (),
    )
    T_a, T_b = 800.0, 400.0
    ξ_s = 600.0
    a2 = CellState{2}(T_a, (ξ_s, 0.0), 1.0)
    b2 = CellState{2}(T_b, (ξ_s, 0.0), 1.0)
    T_lin = merge_temperatures_conservatively(1.0, a2, 1.0, b2, mat_lin)

    # Equal masses: the exact root solves h(T) = (h(T_a)+h(T_b))/2 with
    # h(T) = cp_a·T + cp_b/2·T² (datum cancels):
    h(T) = cp_a * T + 0.5 * cp_b * T^2
    rhs = 0.5 * (h(T_a) + h(T_b))
    T_exact = (-cp_a + sqrt(cp_a^2 + 2 * cp_b * rhs)) / cp_b
    @test T_lin ≈ T_exact atol = 1e-6

    # The pre-fix product-form answer is measurably different (a few K) —
    # pin that the new result does NOT reproduce it.
    ρcp(T) = ξ_s * (cp_a + cp_b * T)
    E_old = ρcp(T_a) * T_a + ρcp(T_b) * T_b
    T_est = E_old / (ρcp(T_a) + ρcp(T_b))
    T_old = E_old / (2 * ξ_s * (cp_a + cp_b * T_est))
    @test abs(T_old - T_exact) > 1.0      # the old bias is > 1 K here
    @test abs(T_lin - T_old) > 1.0

    # (d) Both cells essentially pure gas: condensed balance is degenerate —
    #     volume-weighted fallback.
    a3 = CellState{2}(800.0, (0.0, 10.0), 1.0)
    b3 = CellState{2}(400.0, (0.0, 10.0), 3.0)
    @test merge_temperatures_conservatively(1.0, a3, 3.0, b3, mat_const) ≈
          (1.0 * 800.0 + 3.0 * 400.0) / 4.0
end

# ============================================================================
# Signed-width tangling detection
# ============================================================================
@testset "inverted cells are detected (signed widths)" begin
    mesh = create_uniform_mesh(0.01, 4, Val(2))
    signed_widths(m) =
        [m.nodes[c.node_ids[2]].z - m.nodes[c.node_ids[1]].z for c in m.cells]
    @test all(>(0.0), signed_widths(mesh))

    # Push node 3 past node 4: cell 3 inverts. The cached volumes pass
    # through abs(), so only the signed width sees the inversion.
    z4 = mesh.nodes[4].z
    mesh.nodes[3].z = z4 + 5e-4
    update_geometry!(mesh)
    @test all(>(0.0), mesh.cell_volumes)   # the abs() blind spot, by design
    @test signed_widths(mesh)[3] < 0.0
end

@testset "find_thin_cells refuses inverted cells" begin
    problem = fixture_problem(n_cells = 4)
    def, ws, u, _ = fixture_state(problem; use_ale = true)
    layout = ws.layout

    # Healthy mesh: nothing thin.
    init = fill(0.01 / 4, 4)
    @test isempty(find_thin_cells(u, layout, init, 0.05, 4))

    # Invert cell 3 in the state vector's z-block: a signed thickness below
    # threshold·init is a tangled mesh, not a depleted cell.
    zoff = z_offset(layout)
    u_bad = copy(u)
    u_bad[zoff+3] = u_bad[zoff+4] + 5e-4
    @test_throws ErrorException find_thin_cells(u_bad, layout, init, 0.05, 4)
end

# ============================================================================
# Depletion condition keeps firing at n_live ≤ min_cells
# ============================================================================
@testset "condition fires at the min_cells floor" begin
    n = 4
    problem = fixture_problem(n_cells = n)
    def, ws, u, _ = fixture_state(problem; use_ale = true)

    # min_cells = n: the callback can only terminate, never park; the
    # condition must still fire so the affect!'s terminate branch runs.
    cb_state = DepletionCallbackState(ws.layout, ws.mesh, problem.material, u)
    cb = create_depletion_callback(
        ws,
        ws.mesh,
        problem.material,
        cb_state;
        threshold = 0.05,
        min_cells = n,
    )

    @test cb.condition(u, 0.0, nothing) == false   # healthy mesh: no event

    # Shrink the surface cell below 5% of its initial thickness.
    zoff = z_offset(ws.layout)
    u_thin = copy(u)
    z = get_node_positions(u_thin, ws.layout)
    u_thin[zoff+n+1] = z[n] + 0.01 * (0.01 / n)
    @test cb.condition(u_thin, 0.0, nothing) == true
end
