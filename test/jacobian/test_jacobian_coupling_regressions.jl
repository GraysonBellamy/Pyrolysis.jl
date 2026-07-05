"""
Regression tests for three fixed Jacobian coupling defects:

- **Strict-verifier normalization** — `verify_jacobian_at`'s strict mode
  used a *global* max-|J_ref| normalization, so a wholly missing
  small-magnitude block (the ALE z-rows, ~10 orders below the Arrhenius
  ξ-rows) passed at large n and failed at small n: pass/fail depended on
  problem size, not correctness. The strict criterion is now per-entry
  with a row-magnitude floor. The unit test below drives
  `_strict_per_entry_compare` directly with a synthetic matrix that the
  legacy global criterion provably accepted.

- **z-row × (T, ξ) block** — the structured backend omitted this block
  entirely (dz_k/dt = w_node[k] depends on every upstream cell's strain
  rate via the prefix sum; the geometry AD block covers z *columns*
  only). Verified here at 6 cells — exactly the scale where the pre-fix
  strict verifier failed — by per-entry agreement of the z-row block
  against DenseAD. The WithChi probe additionally covers the
  θ_A → (T, ξ) back-reaction in every w_k.

- **Boundary mass flux on the energy row** — the transport operators
  zero their boundary faces but the residual's S_conv/S_gen include the
  boundary mass flux (filled by `apply_mass_bc_flux!` before the face
  average). Additionally, `BCThermalOperator` used to feed the
  surface Newton a film-law ṁ_g while the residual hard-zeroes it
  (transpiration is rejected at construction) — a Jacobian–residual
  inconsistency active under any permeable mass BC. Probed with
  DiffusiveMassBC at both boundaries, Fickian and Darcy, and the
  `:with_generation_sink` energy form.

  The boundary cell's whole `S_conv` energy row (face-averaged `J̄` ×
  the blended gradient through the surface-Newton `T_s`) is owned by
  `BCConvectionOperator`; `BCMassOperator` keeps only the linear
  `ΔS_gen` boundary-face slice. The same strict per-entry probes gate
  the split, plus a case with a nonlinear surface balance (convective +
  reradiative) where the `∂T_s` chain in the blended gradient is
  load-bearing.
"""

using Test
using SparseArrays
using Pyrolysis
using Pyrolysis: to_problem_def, build_workspace
using Pyrolysis.Jacobian:
    verify_jacobian_at,
    JacobianSpec,
    DenseAD_ForwardDiff,
    Structured,
    DirectSolve,
    _strict_per_entry_compare
using Pyrolysis.Problem: z_range, chi_range

# ----------------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------------

function _build_zrow_ale(; n_cells::Int = 6, with_chi::Bool = false)
    mat = if with_chi
        solid = SolidComponent(
            :wood,
            ρ = 500.0,
            c = 1500.0,
            k = 0.15,
            ε = 0.9,
            κ = 1e-14,
            λ = 1e-5,
        )
        gas = GaseousComponent(:pyrgas, M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
        rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
        Material(
            name = :wood_chi,
            components = (solid, gas),
            reactions = (rxn,),
            lateral_shrinkage_law = χ -> 1.0 - 0.2 * χ,
        )
    else
        virgin = SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
        rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
        Material(name = :PMMA, components = (virgin, gas), reactions = (rxn,))
    end
    mesh = create_uniform_mesh(0.01, n_cells, Val(2))
    bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
    def = to_problem_def(problem; use_ale = true)
    ws = build_workspace(problem, def)
    u = create_state_vector(ws.layout, ws.mesh)
    return def, ws, u
end

function _build_boundary_flux_eulerian(;
    darcy::Bool,
    thermal_top = HeatFluxBC(50e3),
    mass_top = DiffusiveMassBC(),
    mass_bottom = ImpermeableBC(),
    energy_form::Symbol = :standard,
)
    virgin = if darcy
        SolidComponent(:wood, ρ = 500.0, c = 1500.0, k = 0.15, ε = 0.9, κ = 1e-14, λ = 1e-5)
    else
        SolidComponent(:virgin, ρ = 1190, c = 1420, k = 0.21, ε = 0.86, λ = 1e-5)
    end
    gas = GaseousComponent(:gas, M = 0.100, c = 1500, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    mat = Material(name = :m, components = (virgin, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(0.01, 8, Val(2))
    bc = BoundaryConditionSet(
        thermal = (top = thermal_top, bottom = AdiabaticBC()),
        mass = (top = mass_top, bottom = mass_bottom),
    )
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
    def = to_problem_def(problem; use_ale = false, energy_form = energy_form)
    ws = build_workspace(problem, def)
    u = create_state_vector(ws.layout, ws.mesh)
    n = ws.layout.n_cells
    # Temperature ramp (dT/dx ≠ 0 at both boundaries), gas present everywhere
    # — the regime where the boundary flux dominates the energy row.
    for i = 1:n
        u[i] = 600.0 + 400.0 * (i - 1) / (n - 1)
        u[n+i] = darcy ? 400.0 : 1100.0
        u[2n+i] = 8.0
    end
    return def, ws, u
end

function _strict_report(def, ws, u)
    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
    return verify_jacobian_at(spec, def, ws, u, 0.0; reference_spec = ref_spec)
end

# ----------------------------------------------------------------------------
# Per-entry strict criterion
# ----------------------------------------------------------------------------

@testset "strict verifier uses per-entry tolerance" begin
    # Synthetic Jacobian spanning many orders of magnitude: a stiff
    # Arrhenius-like row (1e6) and a small z-like row (1e-2).
    J_ref = [
        1.0e6 2.0e3;
        1.0e-2 3.0e-2
    ]

    # Case 1: the small row's entries are entirely MISSING from the active
    # matrix — the missing-z-row shape. The legacy global criterion accepted
    # this (3e-2 / 1e6 = 3e-8 ≤ rtol 1e-6); per-entry must reject it.
    J_missing = copy(J_ref)
    J_missing[2, 1] = 0.0
    J_missing[2, 2] = 0.0
    @assert maximum(abs, J_missing .- J_ref) / maximum(abs, J_ref) <= 1e-6
    passed, max_abs_dev, _, worst = _strict_per_entry_compare(J_missing, J_ref, 1e-6, 1e-10)
    @test !passed
    @test max_abs_dev ≈ 3.0e-2
    @test worst[1] == 2

    # Case 2: honest roundoff at a true zero of the stiff row passes via the
    # row-magnitude floor (1e-9 ≤ rtol · 1e6).
    J_round = copy(J_ref)
    J_round[1, 2] += 1.0e-9
    passed, _, _, _ = _strict_per_entry_compare(J_round, J_ref, 1e-6, 1e-10)
    @test passed

    # Case 3: an in-support per-entry disagreement in the small row fails
    # even though it is invisible at the global scale.
    J_wrong = copy(J_ref)
    J_wrong[2, 2] = 2.0 * J_ref[2, 2]
    passed, _, _, _ = _strict_per_entry_compare(J_wrong, J_ref, 1e-6, 1e-10)
    @test !passed

    # Exact agreement passes.
    passed, max_abs_dev, _, _ = _strict_per_entry_compare(copy(J_ref), J_ref, 1e-6, 1e-10)
    @test passed
    @test max_abs_dev == 0.0
end

# ----------------------------------------------------------------------------
# z-row × (T, ξ) block
# ----------------------------------------------------------------------------

@testset "structured z-row × (T,ξ) block matches DenseAD at small n" begin
    def, ws, u = _build_zrow_ale()
    n = ws.layout.n_cells
    # Hot state — active Arrhenius so ∂θ/∂T, ∂θ/∂ξ are non-trivial.
    for i = 1:n
        u[i] = 1000.0
        u[n+i] = 1100.0
        u[2n+i] = 0.0
    end

    report = _strict_report(def, ws, u)
    @test report.mode == :strict
    @test report.passed

    # Explicit block check, independent of the verifier: the z-row × (T, ξ)
    # block must be non-trivially populated and agree with DenseAD per entry.
    # (Pre-fix the structured block was identically zero while DenseAD had
    # the exact lower-triangular prefix-sum shape, up to ~4.5e-2.)
    zr = z_range(ws.layout)
    J_s = Matrix(report.J_active)
    J_r = report.J_ref
    Zt_r = J_r[zr, 1:n]
    @test maximum(abs, Zt_r) > 1e-3                      # block is load-bearing
    @test maximum(abs, J_s[zr, 1:n] .- Zt_r) <= 1e-10 * max(1.0, maximum(abs, Zt_r))
    Zx_r = J_r[zr, (n+1):(3n)]
    @test maximum(abs, Zx_r) > 1e-5
    @test maximum(abs, J_s[zr, (n+1):(3n)] .- Zx_r) <= 1e-10 * max(1.0, maximum(abs, Zx_r))
end

@testset "WithChi θ_A back-reaction (z-rows and advection rows)" begin
    def, ws, u = _build_zrow_ale(with_chi = true)
    n = ws.layout.n_cells
    for i = 1:n
        u[i] = 900.0 + 100.0 * (i - 1) / (n - 1)
        u[n+i] = 400.0
        u[2n+i] = 2.0
    end
    # Non-trivial χ so the lateral-shrinkage law and its derivatives bite.
    for (k, idx) in enumerate(chi_range(ws.layout))
        u[idx] = 0.3 + 0.05 * k
    end

    report = _strict_report(def, ws, u)
    @test report.mode == :strict
    @test report.passed

    # The z-row × (T, ξ) reference block is non-trivial in this regime and
    # carries the rank-1 θ_A back-reaction on top of the prefix sum.
    zr = z_range(ws.layout)
    J_s = Matrix(report.J_active)
    J_r = report.J_ref
    @test maximum(abs, J_r[zr, 1:3n]) > 1e-6
    @test maximum(abs, J_s[zr, 1:3n] .- J_r[zr, 1:3n]) <=
          1e-10 * max(1.0, maximum(abs, J_r[zr, 1:3n]))
end

# ----------------------------------------------------------------------------
# Boundary mass flux on the energy row
# ----------------------------------------------------------------------------

@testset "boundary-flux energy row ($(case))" for (case, kwargs) in (
    ("Fickian, diffusive top", (darcy = false,)),
    ("Darcy, diffusive top", (darcy = true,)),
    (
        "Fickian, permeable both boundaries",
        (darcy = false, mass_bottom = DiffusiveMassBC()),
    ),
    ("Fickian, generation sink", (darcy = false, energy_form = :with_generation_sink)),
    # Deferred Minor 5: nonlinear surface balance — the blended boundary
    # gradient's ∂T_s chain (through h and εσT_s⁴) is load-bearing in the
    # boundary cell's S_conv row.
    (
        "Fickian, convective+reradiative top",
        (
            darcy = false,
            thermal_top = HeatFluxBC(50e3) +
                          ConvectiveBC(h = 25.0, T_inf = 300.0) +
                          RadiativeBC(T_inf = 300.0),
        ),
    ),
    (
        "Darcy, convective+reradiative top",
        (
            darcy = true,
            thermal_top = HeatFluxBC(50e3) +
                          ConvectiveBC(h = 25.0, T_inf = 300.0) +
                          RadiativeBC(T_inf = 300.0),
        ),
    ),
)
    def, ws, u = _build_boundary_flux_eulerian(; kwargs...)
    n = ws.layout.n_cells
    gascol_top = 2n + n

    report = _strict_report(def, ws, u)
    @test report.mode == :strict
    @test report.passed

    # Explicit slot check at the worst pre-fix entry: top-cell T-row ×
    # top-cell gas-ξ column (pre-fix: structured −4.05 vs DenseAD/FD −0.170,
    # 24× off). Per-entry agreement against the row scale.
    J_s = Matrix(report.J_active)
    J_r = report.J_ref
    row_scale = maximum(abs, J_r[n, :])
    @test abs(J_s[n, gascol_top] - J_r[n, gascol_top]) <= 1e-9 * row_scale
    @test abs(J_s[n, n] - J_r[n, n]) <= 1e-9 * row_scale
    # Neighbor columns — the blended ∂T/∂x and interior-face J̄ chains of the
    # boundary cell's S_conv row (owned by BCConvectionOperator).
    @test abs(J_s[n, n-1] - J_r[n, n-1]) <= 1e-9 * row_scale
    @test abs(J_s[n, gascol_top-1] - J_r[n, gascol_top-1]) <= 1e-9 * row_scale
    # Bottom boundary cell's T-row (impermeable: exact zero slice; permeable:
    # the mirrored σ-sign case).
    row_scale_b = maximum(abs, J_r[1, :])
    @test maximum(abs, J_s[1, :] .- J_r[1, :]) <= 1e-9 * row_scale_b
end
