using Test
using SparseArrays
using Pyrolysis
using Pyrolysis: to_problem_def, build_workspace
using Pyrolysis.Jacobian:
    verify_jacobian_at, JacobianSpec, DenseAD_ForwardDiff, Structured, DirectSolve

include(joinpath(@__DIR__, "..", "perf", "perf_helpers.jl"))

# ---------------------------------------------------------------------------
# Helpers for fix-specific Jacobian regression tests.
#
# These mutate the state vector to drive the system into regimes where the
# fix in question changes the structured Jacobian — verify_jacobian_at then
# checks that the `Structured` backend agrees with `DenseAD_ForwardDiff`
# (the verifier reference) there.
# ---------------------------------------------------------------------------

# Pack (T, ξ_1, ξ_2, ...) into u given the standard layout:
# u[1:n] = T, u[n+(j-1)n+1 : n+j*n] = ξ_j.
function _set_state!(u, n_cells, T_val, ξ_vals)
    @inbounds for i = 1:n_cells
        u[i] = T_val
    end
    @inbounds for (j, ξj) in enumerate(ξ_vals)
        off = n_cells + (j - 1) * n_cells
        for i = 1:n_cells
            u[off+i] = ξj
        end
    end
    return u
end

# Build a 2-component (solid + gas) Eulerian problem with a single reaction
# carrying user-specified `T_min` / `T_max`. The Arrhenius constants are
# picked so the reaction rate is non-negligible in the 700–1000 K range,
# letting us evaluate the Jacobian at temperatures within ~1 K of the
# window edge where the tanh-window correction term is largest
# (sech²(1)/(1+tanh(1)) ≈ 0.24).
function _build_window_problem(; T_min::Float64, T_max::Float64, n_cells::Int = 50)
    virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(
        :decomp,
        1 => 2,
        A = 1.0e10,
        E = 100e3,
        h = -300e3,
        T_min = T_min,
        T_max = T_max,
    )
    mat = Material(name = :PMMA_window, components = (virgin, gas), reactions = (rxn,))
    mesh = create_uniform_mesh(0.01, n_cells, Val(2))
    bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
    def = to_problem_def(problem; use_ale = false)
    ws = build_workspace(problem, def)
    u = create_state_vector(ws.layout, ws.mesh)
    return def, ws, u
end

# Build the wood-like Beer-Lambert + Darcy ALE problem (matches
# `build_ale_beer_lambert_darcy` in perf_helpers.jl) — exposed here so we
# can mutate `u` to a high-temperature, non-zero-gas state where the
# Sutherland viscosity contribution to ∂v_gas/∂T is meaningful.
_build_darcy_problem() = build_ale_beer_lambert_darcy()

function _build_generation_sink_problem(; n_cells::Int = 8)
    virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
    gas = GaseousComponent(:MMA, M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    mat = Material(
        name = :PMMA_generation_sink,
        components = (virgin, gas),
        reactions = (rxn,),
    )
    mesh = create_uniform_mesh(0.01, n_cells, Val(2))
    bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
    def = to_problem_def(problem; use_ale = false, energy_form = :with_generation_sink)
    ws = build_workspace(problem, def)
    u = create_state_vector(ws.layout, ws.mesh)
    return def, ws, u
end

# Generation-sink problem with a **temperature-dependent** gas cp. The S_gen
# source uses Δh_g(T,T_ref) = cp_g(T_m)·(T−T_ref) (midpoint rule, T_m=(T+T_ref)/2),
# whose ∂/∂T = cp_g(T_m) + ½(T−T_ref)·cp_g′(T_m). Constant-cp problems have
# cp_g′≡0 and so cannot exercise that half-slope chain term; a LinearProperty cp
# makes cp_g′ = b ≠ 0, validating the analytic S_gen Jacobian against AD.
function _build_generation_sink_tdep_cp_problem(; n_cells::Int = 8)
    virgin = SolidComponent(:virgin, ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
    gas = GaseousComponent(
        :MMA,
        M = 0.100,
        c = LinearProperty(1100.0, 0.8),
        k = 0.02,
        λ = 1e-5,
    )
    rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
    mat = Material(
        name = :PMMA_generation_sink_tdep,
        components = (virgin, gas),
        reactions = (rxn,),
    )
    mesh = create_uniform_mesh(0.01, n_cells, Val(2))
    bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
    problem =
        PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
    Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
    def = to_problem_def(problem; use_ale = false, energy_form = :with_generation_sink)
    ws = build_workspace(problem, def)
    u = create_state_vector(ws.layout, ws.mesh)
    return def, ws, u
end

function _set_generation_sink_state!(u, n_cells)
    @inbounds for i = 1:n_cells
        u[i] = 900.0 + 10.0 * i
        u[n_cells+i] = 700.0 - 20.0 * i
        u[2*n_cells+i] = 1.0 + 0.3 * i
    end
    return u
end

@testset "Jacobian regression fixes" begin
    @testset "Structured excludes gas from rho-cp derivatives" begin
        def, ws, u, _ = build_baseline_eulerian()
        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        report = verify_jacobian_at(
            spec,
            def,
            ws,
            u,
            0.0;
            reference_spec = ref_spec,
            rtol = 1e-6,
            atol = 1e-10,
        )
        @test report.passed
    end

    # -----------------------------------------------------------------------
    # Fix 1: T_min / T_max tanh-window derivative
    #
    # The Arrhenius rate is multiplied by 0.5(1+tanh(T-T_min)) and
    # 0.5(1+tanh(T_max-T)). Evaluating the analytical Jacobian within the
    # windows but close to either edge (sech² is non-trivial only there)
    # fails without the chain-rule additions in
    # `evaluate_reaction_derivatives!`.
    # -----------------------------------------------------------------------
    @testset "Structured includes T_min/T_max tanh-window derivative" begin
        # Lower-edge case: T = 750, T_min = 749 → tanh(1) regime, ~24% correction.
        let
            def, ws, u = _build_window_problem(T_min = 749.0, T_max = 1500.0)
            n = ws.layout.n_cells
            _set_state!(u, n, 750.0, (1100.0, 0.0))
            ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
            spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
            report = verify_jacobian_at(
                spec,
                def,
                ws,
                u,
                0.0;
                reference_spec = ref_spec,
                rtol = 1e-6,
                atol = 1e-10,
            )
            @test report.passed
        end

        # Upper-edge case: T = 1499, T_max = 1500 → symmetric tanh(1) regime.
        let
            def, ws, u = _build_window_problem(T_min = 400.0, T_max = 1500.0)
            n = ws.layout.n_cells
            _set_state!(u, n, 1499.0, (1100.0, 0.0))
            ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
            spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
            report = verify_jacobian_at(
                spec,
                def,
                ws,
                u,
                0.0;
                reference_spec = ref_spec,
                rtol = 1e-6,
                atol = 1e-10,
            )
            @test report.passed
        end
    end

    # -----------------------------------------------------------------------
    # Fix 2: Sutherland viscosity path in Darcy ∂v_gas/∂T
    #
    # Without the µ(T) chain-rule term, ∂v_gas/∂T misses ~60% of its
    # magnitude at T ≈ 1000 K. Drive the wood-like Darcy problem to a
    # high-T, non-zero-gas state so the pressure gradient (and hence the
    # viscosity contribution) is non-trivial.
    # -----------------------------------------------------------------------
    @testset "Structured includes Sutherland viscosity in Darcy ∂v_gas/∂T" begin
        def, ws, u, _ = _build_darcy_problem()
        n = ws.layout.n_cells
        # Set T = 1000 K (Sutherland correction ~60% of pressure path),
        # partially-decomposed virgin and a non-zero gas concentration so
        # Darcy flux and ∂P/∂{T,ξ} are both active.
        _set_state!(u, n, 1000.0, (400.0, 5.0))
        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        report = verify_jacobian_at(
            spec,
            def,
            ws,
            u,
            0.0;
            reference_spec = ref_spec,
            rtol = 1e-6,
            atol = 1e-10,
        )
        @test report.passed
    end

    @testset "Structured includes generation-sink derivative" begin
        def, ws, u = _build_generation_sink_problem()
        n = ws.layout.n_cells
        _set_generation_sink_state!(u, n)
        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        report = verify_jacobian_at(
            spec,
            def,
            ws,
            u,
            0.0;
            reference_spec = ref_spec,
            rtol = 1e-6,
            atol = 1e-10,
        )
        @test report.passed
        @test report.max_abs_dev_outside_support <= 1e-10
    end

    @testset "Structured generation-sink derivative (T-dependent cp)" begin
        # Validates the Δh_g = cp_g(T_m)·(T−T_ref) midpoint-rule S_gen Jacobian
        # (cp evaluated at T_m, explicit chain factor ½·cp_g′(T_m)) against AD.
        def, ws, u = _build_generation_sink_tdep_cp_problem()
        n = ws.layout.n_cells
        _set_generation_sink_state!(u, n)
        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        report = verify_jacobian_at(
            spec,
            def,
            ws,
            u,
            0.0;
            reference_spec = ref_spec,
            rtol = 1e-6,
            atol = 1e-10,
        )
        @test report.passed
        @test report.max_abs_dev_outside_support <= 1e-10
    end

    # -----------------------------------------------------------------------
    # Fix 3: Eulerian dS_dξ scratch slice consistency ([:,:,1] → [:,:,i])
    #
    # Spot-check at a non-default state where reaction rates differ
    # between cells. With a uniform u, every cell's slot would be written
    # identically and the slice-index bug would be invisible — so impose a
    # linear T gradient across the column and a non-zero gas concentration.
    # -----------------------------------------------------------------------
    @testset "Structured Eulerian dS_dξ uses per-cell scratch slice" begin
        def, ws, u, _ = build_baseline_eulerian()
        n = ws.layout.n_cells
        # Linear ramp in T across cells (650 K bottom → 800 K top) and a
        # small uniform gas concentration; reaction rates then vary cell
        # by cell, so the per-cell dS_dξ writes must hit the right slice.
        @inbounds for i = 1:n
            u[i] = 650.0 + (800.0 - 650.0) * (i - 1) / max(n - 1, 1)
        end
        off_v = n          # virgin block start
        off_g = 2 * n      # gas block start
        @inbounds for i = 1:n
            u[off_v+i] = 1100.0
            u[off_g+i] = 1.0
        end
        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
        report = verify_jacobian_at(
            spec,
            def,
            ws,
            u,
            0.0;
            reference_spec = ref_spec,
            rtol = 1e-6,
            atol = 1e-10,
        )
        @test report.passed
    end

    # -----------------------------------------------------------------------
    # Fix 6: apply_limiter type-stable for Real (was Float64-only)
    #
    # At default state ws.w_work ≈ 0 and gradients are zero, so the
    # `_limited_gradient` early-return on `abs(grad_1) < 1e-20` skips the
    # `apply_limiter` call. AD then never hits the type wall. With a
    # non-uniform T ramp + non-zero shrinkage rate, AD propagates Duals
    # through `apply_limiter`, which previously errored
    # (MethodError: apply_limiter(::Dual, ::Limiter)).
    #
    # This test fails to even compute the AD reference Jacobian without
    # the fix, so passing here is a regression guard.
    # -----------------------------------------------------------------------
    @testset "apply_limiter accepts Dual (AD-through-ALE-residual)" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        # Strong T gradient + non-zero gas → the ALE limited gradient is
        # fully active (no early return on uniform-state).
        @inbounds for i = 1:n
            u[i] = 600.0 + (1100.0 - 600.0) * (i - 1) / max(n - 1, 1)
        end
        off_v = n
        off_g = 2 * n
        @inbounds for i = 1:n
            u[off_v+i] = 800.0
            u[off_g+i] = 1.5
        end
        # Just computing the AD reference (which runs the residual on
        # Dual numbers) was the failure mode.
        ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
        J_ref = zeros(Float64, length(u), length(u))
        @test_nowarn Pyrolysis.Jacobian.compute_jacobian!(J_ref, ref_spec, u, def, 0.0, ws)
        @test all(isfinite, J_ref)
    end
end

# ============================================================================
# Off-default-state regression suite. Each canonical config is driven into
# a physical regime that exercises a different chain rule, and the single
# `Structured(strategy = DirectSolve())` backend is verified against
# `DenseAD_ForwardDiff` (verifier reference) at rtol = 1e-6 / atol = 1e-10.
# ============================================================================

function _assert_backends_agree(def, ws, u; rtol = 1e-6, atol = 1e-10)
    ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
    spec = JacobianSpec(Structured(strategy = DirectSolve()), def, ws)
    report = verify_jacobian_at(
        spec,
        def,
        ws,
        u,
        0.0;
        reference_spec = ref_spec,
        rtol = rtol,
        atol = atol,
    )
    @test report.passed
end

@testset "Off-default-state regressions" begin
    @testset "Eulerian PMMA at hot state (T=1000)" begin
        def, ws, u, _ = build_baseline_eulerian()
        n = ws.layout.n_cells
        _set_state!(u, n, 1000.0, (1100.0, 0.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "Eulerian PMMA at hot+gas-heavy state" begin
        def, ws, u, _ = build_baseline_eulerian()
        n = ws.layout.n_cells
        _set_state!(u, n, 1200.0, (400.0, 8.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "Eulerian PMMA partially decomposed" begin
        def, ws, u, _ = build_baseline_eulerian()
        n = ws.layout.n_cells
        _set_state!(u, n, 800.0, (700.0, 2.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "Eulerian surface-absorption at hot+gas state" begin
        def, ws, u, _ = build_eulerian_surface_absorption()
        n = ws.layout.n_cells
        _set_state!(u, n, 1000.0, (700.0, 5.0))
        _assert_backends_agree(def, ws, u)
    end

    # ALE configs at uniform off-default states (NOT ramp — the simple-upwind /
    # VANLEER divergence in compute_ale_advection_jacobian_*! shows up there).
    @testset "ALE no-radiation at hot state" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        _set_state!(u, n, 1000.0, (1100.0, 0.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "ALE no-radiation at hot+gas-heavy state" begin
        def, ws, u, _ = build_ale_no_radiation()
        n = ws.layout.n_cells
        _set_state!(u, n, 1200.0, (400.0, 8.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "ALE Beer-Lambert + Darcy at hot state" begin
        def, ws, u, _ = build_ale_beer_lambert_darcy()
        n = ws.layout.n_cells
        # Wood material: ρ_virgin = 500, so keep ξ_virgin < 500 (porosity ≥ 0).
        _set_state!(u, n, 1000.0, (400.0, 0.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "ALE Beer-Lambert + Darcy at hot+gas state (Sutherland active)" begin
        def, ws, u, _ = build_ale_beer_lambert_darcy()
        n = ws.layout.n_cells
        _set_state!(u, n, 1000.0, (400.0, 5.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "ALE Beer-Lambert + Darcy partially decomposed" begin
        def, ws, u, _ = build_ale_beer_lambert_darcy()
        n = ws.layout.n_cells
        _set_state!(u, n, 800.0, (200.0, 2.0))
        _assert_backends_agree(def, ws, u)
    end
end

# ============================================================================
# Limit-case regressions: small/large meshes, zero-rate, and 3-component
# materials. These probe edge cases the canonical 100-cell × 2-component
# CI matrix doesn't reach.
# ============================================================================

@testset "Limit-case regressions" begin
    @testset "Single-cell mesh" begin
        virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA; M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
        rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
        mat = Material(name = :PMMA_1cell, components = (virgin, gas), reactions = (rxn,))
        mesh = create_uniform_mesh(0.01, 1, Val(2))
        bc =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
        def = to_problem_def(problem; use_ale = false)
        ws = build_workspace(problem, def)
        u = create_state_vector(ws.layout, ws.mesh)
        _set_state!(u, 1, 1000.0, (1100.0, 0.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "200-cell mesh at hot+gas state" begin
        virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA; M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
        rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3)
        mat = Material(name = :PMMA_dense, components = (virgin, gas), reactions = (rxn,))
        mesh = create_uniform_mesh(0.01, 200, Val(2))
        bc =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
        def = to_problem_def(problem; use_ale = false)
        ws = build_workspace(problem, def)
        u = create_state_vector(ws.layout, ws.mesh)
        _set_state!(u, 200, 1000.0, (1100.0, 0.5))
        _assert_backends_agree(def, ws, u)
    end

    @testset "Zero-rate cold state (T well below T_min)" begin
        virgin = SolidComponent(:virgin; ρ = 1190.0, c = 1420.0, k = 0.21, ε = 0.86)
        gas = GaseousComponent(:MMA; M = 0.100, c = 1500.0, k = 0.02, λ = 1e-5)
        # T_min = 800 → at T=400 the rate is multiplied by 0.5(1+tanh(-400)) ≈ 0.
        rxn = Reaction(:decomp, 1 => 2, A = 8.5e12, E = 188e3, h = -870e3, T_min = 800.0)
        mat = Material(name = :PMMA_dorm, components = (virgin, gas), reactions = (rxn,))
        mesh = create_uniform_mesh(0.01, 30, Val(2))
        bc =
            BoundaryConditionSet(thermal = (top = HeatFluxBC(50e3), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
        def = to_problem_def(problem; use_ale = false)
        ws = build_workspace(problem, def)
        u = create_state_vector(ws.layout, ws.mesh)
        _set_state!(u, 30, 400.0, (1100.0, 1.0))
        _assert_backends_agree(def, ws, u)
    end

    @testset "3-component material (κ-mixing chain rule)" begin
        # Both solids carry permeability with different intrinsic κ → the
        # weighted κ-mixing rule fires non-trivially. ξ chosen to keep
        # porosity > 0: ξ_virgin/500 + ξ_char/200 = 0.65 → φ ≈ 0.35.
        virgin = SolidComponent(
            :wood;
            ρ = 500.0,
            c = 1500.0,
            k = 0.15,
            ε = 0.9,
            α = 1000.0,
            κ = 1e-14,
            λ = 1e-5,
        )
        char = SolidComponent(
            :char;
            ρ = 200.0,
            c = 1100.0,
            k = 0.10,
            ε = 0.95,
            α = 5000.0,
            κ = 1e-12,
            λ = 1e-5,
        )
        gas = GaseousComponent(:pyrgas; M = 0.06, c = 1800.0, k = 0.03, λ = 5e-5)
        rxn = Reaction(:pyrolysis, 1 => 2, A = 1e10, E = 150e3, h = -300e3)
        rxn2 = Reaction(:char_form, 1 => 3, A = 1e8, E = 120e3, h = -200e3)
        mat = Material(
            name = :wood3,
            components = (virgin, char, gas),
            reactions = (rxn, rxn2),
        )
        mesh = create_uniform_mesh(0.01, 50, Val(3))
        bc = BoundaryConditionSet(thermal = (top = HeatFluxBC(0.0), bottom = AdiabaticBC()))
        problem =
            PyrolysisProblem(mesh = mesh, material = mat, bc_set = bc, tspan = (0.0, 100.0))
        Pyrolysis.Problem.apply_initial_conditions!(problem.mesh, problem)
        def = to_problem_def(
            problem;
            use_ale = true,
            radiation_model = Pyrolysis.Materials.BEER_LAMBERT,
            radiation_intensity = 50e3,
        )
        ws = build_workspace(problem, def)
        u = create_state_vector(ws.layout, ws.mesh)
        _set_state!(u, 50, 900.0, (200.0, 50.0, 3.0))
        _assert_backends_agree(def, ws, u)
    end
end
