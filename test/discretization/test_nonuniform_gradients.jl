"""
Regression tests for nonuniform-mesh gradients: cell-centered gradients and
the limited reconstruction were first-order on exactly the nonuniform/AMR
meshes the package targets.

- The wide stencil (φ_{i+1} − φ_{i−1})/(d_L + d_R) and the unweighted mean of
  one-sided slopes both carry an O(|d_R − d_L|·φ″) error. The fix is
  opposite-distance weighting (`opposite_distance_gradient`), which is exact
  on quadratics for any spacing — that is the discriminating test here
  (linear fields cannot distinguish the schemes).
- `_limited_gradient`'s blend constant ½ was the uniform-mesh value; the
  distance-correct blend is β = d₀/(d₀+d₁).
- The structured-Jacobian transport operators mirror the residual's S_conv
  dTdx stencil; a strict per-entry structured-vs-DenseAD comparison on a
  *stretched* mesh (n = 8 — small n is where global normalization used to
  hide defects) gates the mirror.
"""

using Test
using Pyrolysis
using Pyrolysis.Geometry: opposite_distance_gradient, opposite_distance_gradient_coeffs
using Pyrolysis.Problem: StateCache
using Pyrolysis.Discretization: _limited_gradient, CENTRAL, VANLEER
using Pyrolysis.Jacobian:
    verify_jacobian_at, JacobianSpec, DenseAD_ForwardDiff, Structured, DirectSolve

# Quadratic test field: exact for any second-order gradient, first-order
# schemes miss by O(|d_R − d_L|·2c₂).
_quad(z) = 300.0 + 5.0e4 * z + 2.0e6 * z^2
_dquad(z) = 5.0e4 + 4.0e6 * z

@testset "nonuniform-mesh gradients" begin

    @testset "kernel exact on quadratics, any spacing" begin
        for (d_L, d_R) in ((1e-4, 3e-4), (5e-4, 1e-4), (2.5e-4, 2.5e-4))
            z_C = 0.004
            φ_L, φ_C, φ_R = _quad(z_C - d_L), _quad(z_C), _quad(z_C + d_R)
            g = opposite_distance_gradient(φ_L, φ_C, φ_R, d_L, d_R)
            @test g ≈ _dquad(z_C) rtol = 1e-12

            # coefficient form is the same operator, and ∂/∂φ entries sum to 0
            (c_L, c_C, c_R) = opposite_distance_gradient_coeffs(d_L, d_R)
            @test c_L * φ_L + c_C * φ_C + c_R * φ_R ≈ g rtol = 1e-12
            @test abs(c_L + c_C + c_R) < 1e-12 / d_L

            # the pre-fix forms are NOT exact off uniform spacing
            if d_L != d_R
                wide = (φ_R - φ_L) / (d_L + d_R)
                halfmean = ((φ_C - φ_L) / d_L + (φ_R - φ_C) / d_R) / 2
                @test !isapprox(wide, _dquad(z_C); rtol = 1e-6)
                @test !isapprox(halfmean, _dquad(z_C); rtol = 1e-6)
            end
        end
    end

    @testset "S_conv gradient second order on stretched mesh" begin
        solid = SolidComponent(:s, ρ = 500.0, c = 1500.0, k = 0.15)
        gas = GaseousComponent(:g, M = 0.044, c = 1200.0, k = 0.03, λ = 2e-5)
        mat = Material(name = :gradient_fixture, components = (solid, gas), reactions = ())

        n = 12
        mesh = create_stretched_mesh(0.01, n, Val(2), 1.5)
        for i = 1:n
            # Physical solid mass so k_eff > 0 (the boundary-face slope is
            # recovered as −q/k).
            mesh.cell_states[i] = Pyrolysis.Geometry.CellState{2}(
                _quad(mesh.cell_centers[i]),
                (500.0, 0.0),
                mesh.cell_volumes[i],
            )
        end
        state = StateCache(mesh)
        prop_cache = Pyrolysis.Properties.PropertyCache(n, n + 1)
        Pyrolysis.Properties.update_properties!(prop_cache, state, mesh, mat)

        J = 0.01  # constant gas flux at every face → J_avg = J exactly
        J_gas = zeros(2, n + 1)
        J_gas[2, :] .= J

        # Boundary-face fluxes consistent with the manufactured profile, in
        # the form the thermal BC writes them: q[f] = ∓k_i·(T_s − T_i)/Δz_half
        # with T_s = T(z_face). The kernel recovers the face slope as −q/k.
        f_bot = mesh.cells[1].face_ids[1]
        f_top = mesh.cells[n].face_ids[2]
        q_faces = zeros(n + 1)
        dzh_bot = state.volume[1] / state.A_face[f_bot] / 2
        dzh_top = state.volume[n] / state.A_face[f_top] / 2
        T_s_bot = _quad(state.face_positions[f_bot])
        T_s_top = _quad(state.face_positions[f_top])
        q_faces[f_bot] = -prop_cache.k[1] * (state.T[1] - T_s_bot) / dzh_bot
        q_faces[f_top] = -prop_cache.k[n] * (T_s_top - state.T[n]) / dzh_top

        S_conv = zeros(n)
        Pyrolysis.Physics.compute_gas_convective_source!(
            S_conv,
            J_gas,
            q_faces,
            prop_cache.k,
            state,
            mesh,
            mat,
        )
        cp_g = 1200.0
        for i = 2:(n-1)
            @test S_conv[i] ≈ -cp_g * J * _dquad(mesh.cell_centers[i]) rtol = 1e-10
        end

        # Boundary cells: the surface-pseudo-neighbor blend interpolates the
        # exact face slope (at z_c ± Δz_half/2) and the interior slope (at
        # z_c ∓ d_nb/2) back to z_c — exact on quadratics. The pre-fix
        # interior-one-sided gradient is provably NOT (first-order error
        # ∝ curvature on this stretched mesh).
        for (i, nb) in ((1, 2), (n, n - 1))
            exact = -cp_g * J * _dquad(mesh.cell_centers[i])
            @test S_conv[i] ≈ exact rtol = 1e-10
            d_nb = abs(mesh.cell_centers[nb] - mesh.cell_centers[i])
            s_int = (state.T[i] - state.T[nb]) / d_nb * (i == 1 ? -1.0 : 1.0)
            old_one_sided = -cp_g * J * s_int
            @test !isapprox(old_one_sided, exact; rtol = 1e-6)
        end
    end

    @testset "AMR indicator gradient second order on stretched mesh" begin
        n = 12
        mesh = create_stretched_mesh(0.01, n, Val(2), 1.5)
        φ = [_quad(mesh.cell_centers[i]) for i = 1:n]
        for i = 2:(n-1)
            g = Pyrolysis.Adaptivity.compute_gradient_1d(φ, mesh, i)
            @test g ≈ _dquad(mesh.cell_centers[i]) rtol = 1e-10
        end
    end

    @testset "_limited_gradient distance-correct blend (β = d₀/(d₀+d₁))" begin
        # Donor-side quadratic reconstruction with CENTRAL (ψ ≡ 1): the
        # distance-correct blend is exact on quadratics for any (d₀, d₁);
        # the old fixed-½ blend missed by c·(d₁−d₀)/2... ≠ 0.
        c2 = 2.0e6
        f(z) = 300.0 + 5.0e4 * z + c2 * z^2
        df(z) = 5.0e4 + 2 * c2 * z

        for (d0, d1) in ((1e-4, 3e-4), (4e-4, 1e-4))
            # donor on the LEFT of cell i at z = 0
            T_i, T_up, T_upup = f(0.0), f(-d0), f(-d0 - d1)
            g = _limited_gradient(T_i, T_up, T_upup, d0, d1, true, CENTRAL, true)
            @test g ≈ df(0.0) rtol = 1e-12

            # donor on the RIGHT
            T_up_r, T_upup_r = f(d0), f(d0 + d1)
            g_r = _limited_gradient(T_i, T_up_r, T_upup_r, d0, d1, true, CENTRAL, false)
            @test g_r ≈ df(0.0) rtol = 1e-12

            # the old ½ blend is provably NOT exact here
            grad_1 = (T_i - T_up) / d0
            grad_0 = (T_up - T_upup) / d1
            old = grad_1 + 0.5 * (grad_1 - grad_0)
            d0 == d1 || @test !isapprox(old, df(0.0); rtol = 1e-6)
        end

        # uniform spacing: β = ½ — identical to the historical scheme
        d = 2e-4
        T_i, T_up, T_upup = f(0.0), f(-d), f(-2d)
        g_u = _limited_gradient(T_i, T_up, T_upup, d, d, true, CENTRAL, true)
        grad_1 = (T_i - T_up) / d
        grad_0 = (T_up - T_upup) / d
        @test g_u ≈ grad_1 + 0.5 * (grad_1 - grad_0) rtol = 1e-14

        # TVD sanity on monotone data with VANLEER and skewed distances:
        # the limited gradient stays within the hull of the one-sided slopes.
        T_i, T_up, T_upup = 900.0, 700.0, 650.0
        d0, d1 = 1e-4, 4e-4
        g_tvd = _limited_gradient(T_i, T_up, T_upup, d0, d1, true, VANLEER, true)
        s_near = (T_i - T_up) / d0
        s_far = (T_up - T_upup) / d1
        lo, hi = min(s_near, s_far), max(s_near, s_far)
        @test lo - 1e-9 <= g_tvd
        @test g_tvd <= 2 * hi   # ψ ≤ 2 bound
    end
end

# ----------------------------------------------------------------------------
# Jacobian mirror: structured vs DenseAD on a STRETCHED mesh at small n.
# The transport operators hand-mirror the residual's S_conv dTdx stencil
# (apply_cell! AD path + analytical override coefficients); on a uniform mesh
# the new weights collapse to the old ½, so only a stretched mesh exercises
# the mirror.
# ----------------------------------------------------------------------------

function _build_stretched_problem(; darcy::Bool, use_ale::Bool = false)
    virgin = if darcy
        SolidComponent(
            :wood,
            ρ = 500.0,
            c = 1500.0,
            k = 0.15,
            ε = 0.9,
            κ = 1e-14,
            λ = 1e-5,
            ρ_intrinsic = 1500.0,
        )
    else
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86, λ = 1e-5)
    end
    gas = GaseousComponent(:gas, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = 870e3)
    mat = Material(name = :m21j, components = (virgin, gas), reactions = (rxn,))
    mesh = create_stretched_mesh(0.01, 8, Val(2), 1.6)
    bc = BoundaryConditionSet(
        thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()),
        mass = (top = DiffusiveMassBC(0.05), bottom = ImpermeableBC()),
    )
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    def, ws, u = initialize_state(problem; use_ale = use_ale)
    n = ws.layout.n_cells
    # Temperature ramp + gas everywhere: non-trivial dTdx and S_conv on the
    # stretched stencil.
    for i = 1:n
        u[i] = 600.0 + 400.0 * (i - 1) / (n - 1)
        u[n+i] = darcy ? 400.0 : 1100.0
        u[2n+i] = 8.0
    end
    return def, ws, u
end

@testset "structured Jacobian mirrors stretched-mesh dTdx" begin
    for (label, darcy, use_ale) in
        (("Fickian", false, false), ("Darcy", true, false), ("ALE Fickian", false, true))
        @testset "$label" begin
            def, ws, u = _build_stretched_problem(darcy = darcy, use_ale = use_ale)
            ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
            spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
            report = verify_jacobian_at(spec, def, ws, u, 0.0; reference_spec = ref_spec)
            @test report.mode == :strict
            @test report.passed
        end
    end
end
