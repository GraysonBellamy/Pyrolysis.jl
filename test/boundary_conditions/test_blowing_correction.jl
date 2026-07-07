"""
Blowing (transpiration) correction tests.

Covers the Couette-film blowing correction `h → h·B/(e^B − 1)`,
`B = ṁ″·c_p,g/h` (opt-in via `use_blowing_correction`):

1. `blowing_factor` — limits, branch-seam continuity (value + derivative),
   large-B overflow safety under ForwardDiff, suction enhancement.
2. `evaluate_bc_flux_and_derivative` mcp dispatch — only convective BCs
   rescale; radiative/flux BCs ignore mcp; combined BCs forward it.
3. Surface Newton with blowing — converged balance, `blowing_mcp`
   diagnostic, off/on ordering, exact off-state identity.
4. AD cleanliness — ForwardDiff derivatives of `T_s` w.r.t. `T_cell` and
   the blowing mass flux match central finite differences.
5. Structured-backend Jacobian — strict support-complete agreement with
   `DenseAD_ForwardDiff` with blowing enabled (Fickian and Darcy–Fick,
   Eulerian and ALE).
6. Integration — `solve(use_blowing_correction = true)` raises the surface
   temperature under convective cooling, and the saved flux decomposition
   stays consistent with the boundary ledger.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.Geometry: BoundaryState
using Pyrolysis.Discretization:
    blowing_factor,
    solve_boundary_temperature,
    solve_boundary_temperature!,
    evaluate_bc_flux_and_derivative,
    compute_boundary_heat_fluxes
using Pyrolysis.Discretization: _blowing_conductance
using Pyrolysis.Materials: SURFACE_ABSORPTION, NO_RADIATION
using Pyrolysis.Jacobian: verify_jacobian_at, JacobianSpec, DenseAD_ForwardDiff, Structured

# ----------------------------------------------------------------------------
# Shared fixtures: 2-component solid + gas, constant cp_g so the blowing
# conductance is exactly J_gas · c_g and every balance is closed-form.
# ----------------------------------------------------------------------------
const _CP_GAS = 1000.0

function _blowing_material(; darcy::Bool = false)
    solid_kwargs = darcy ? (; κ = 1e-14) : (;)
    solid =
        SolidComponent(:solid; ρ = 1000.0, c = 1000.0, k = 0.2, ε = 0.9, solid_kwargs...)
    gas = GaseousComponent(:gas; M = 0.050, c = _CP_GAS, k = 0.02, λ = 1e-5)
    rxn = Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5)
    return Material(name = :blowing_test, components = (solid, gas), reactions = (rxn,))
end

# Reference factor straight from the definition (only valid away from 0/∞).
_f_exact(B) = B / (exp(B) - 1)

@testset "Blowing correction" begin

    # ========================================================================
    @testset "blowing_factor: values and limits" begin
        @test blowing_factor(0.0) === 1.0
        # Audit-verified arithmetic: B = 4 → h/h₀ = 0.0746…
        @test blowing_factor(4.0) ≈ _f_exact(4.0) rtol = 1e-13
        @test blowing_factor(4.0) ≈ 0.0746 atol = 1e-4
        @test blowing_factor(1.0) ≈ 1 / (ℯ - 1) rtol = 1e-13
        # Suction (B < 0) enhances.
        @test blowing_factor(-2.0) ≈ _f_exact(-2.0) rtol = 1e-13
        @test blowing_factor(-2.0) > 1.0
        # Monotone decreasing.
        grid = [-5.0, -1.0, -1e-6, 0.0, 1e-6, 0.5, 1.0, 4.0, 20.0, 200.0]
        vals = blowing_factor.(grid)
        @test all(diff(vals) .< 0)
        @test all(v -> v > 0, vals)
    end

    @testset "blowing_factor: branch seams (value + derivative)" begin
        for B0 in (1e-5, -1e-5)
            inner = blowing_factor(prevfloat(abs(B0)) * sign(B0))  # Taylor side
            outer = blowing_factor(nextfloat(abs(B0)) * sign(B0))  # exact side
            @test inner ≈ outer rtol = 1e-12
            d_inner = ForwardDiff.derivative(blowing_factor, prevfloat(abs(B0)) * sign(B0))
            d_outer = ForwardDiff.derivative(blowing_factor, nextfloat(abs(B0)) * sign(B0))
            @test d_inner ≈ d_outer rtol = 1e-8
        end
        # f′(0) = −1/2 from the series.
        @test ForwardDiff.derivative(blowing_factor, 0.0) ≈ -0.5 rtol = 1e-10
    end

    @testset "blowing_factor: large-B overflow safety" begin
        # Values: clean underflow to 0, never Inf/NaN.
        for B in (50.0, 700.0, 720.0, 1e4)
            v = blowing_factor(B)
            @test isfinite(v) && v >= 0.0
        end
        @test blowing_factor(700.0) ≈ 700.0 * exp(-700.0) rtol = 1e-12
        # ForwardDiff partials stay finite where the naive B/expm1(B)
        # produces Inf/Inf = NaN (expm1(720) = Inf).
        for B in (700.0, 720.0, 1e4)
            d = ForwardDiff.derivative(blowing_factor, B)
            @test isfinite(d)
        end
        # Strong suction is finite too: f(−B) → |B|.
        @test blowing_factor(-700.0) ≈ 700.0 rtol = 1e-12
        @test isfinite(ForwardDiff.derivative(blowing_factor, -700.0))
    end

    # ========================================================================
    @testset "evaluate_bc_flux_and_derivative: mcp dispatch" begin
        T_s, t, ε = 600.0, 0.0, 0.9
        h, T_inf = 20.0, 300.0
        conv = ConvectiveBC(h = h, T_inf = T_inf)
        mcp = 20.0                      # B = 1
        χ = blowing_factor(mcp / h)

        q, dq = evaluate_bc_flux_and_derivative(conv, T_s, t, ε, NO_RADIATION, mcp)
        @test q ≈ h * χ * (T_inf - T_s) rtol = 1e-14
        @test dq ≈ -h * χ rtol = 1e-14

        # mcp = 0 reproduces the 5-arg flux exactly.
        q0, dq0 = evaluate_bc_flux_and_derivative(conv, T_s, t, ε, NO_RADIATION, 0.0)
        q5, dq5 = evaluate_bc_flux_and_derivative(conv, T_s, t, ε, NO_RADIATION)
        @test q0 === q5
        @test dq0 === dq5

        # Non-convective BCs ignore mcp.
        rad = RadiativeBC(T_inf = T_inf)
        qr, dqr = evaluate_bc_flux_and_derivative(rad, T_s, t, ε, NO_RADIATION, mcp)
        qr5, dqr5 = evaluate_bc_flux_and_derivative(rad, T_s, t, ε, NO_RADIATION)
        @test qr == qr5 && dqr == dqr5

        flux = HeatFluxBC(2e4)
        qf, _ = evaluate_bc_flux_and_derivative(flux, T_s, t, ε, NO_RADIATION, mcp)
        @test qf == 2e4

        # Combined: convective child rescaled, others untouched.
        combo = flux + conv + rad
        qc, dqc = evaluate_bc_flux_and_derivative(combo, T_s, t, ε, NO_RADIATION, mcp)
        @test qc ≈ qf + q + qr rtol = 1e-14
        @test dqc ≈ dq + dqr rtol = 1e-14
    end

    # ========================================================================
    @testset "surface Newton with blowing: converged balance" begin
        mat = _blowing_material()
        h, T_inf = 20.0, 300.0
        bc = HeatFluxBC(2e4) + ConvectiveBC(h = h, T_inf = T_inf)
        T_cell, k_eff, Δz_half, t, ε = 600.0, 0.2, 5e-4, 0.0, 0.9
        J_gas = 0.02                                   # kg/(m²·s) outward
        J_zero = (0.0, 0.0)
        J_blow = (0.0, J_gas)

        state = solve_boundary_temperature(
            T_cell,
            k_eff,
            Δz_half,
            bc,
            mat,
            J_zero,
            T_inf,
            t,
            ε;
            use_blowing = true,
            J_blowing = J_blow,
        )
        @test state.converged

        # Constant cp_g → mcp is exactly J·cp, independent of T_s.
        mcp = J_gas * _CP_GAS
        @test state.blowing_mcp ≈ mcp rtol = 1e-12
        χ = blowing_factor(mcp / h)

        # Newton balance: q_applied = q_conducted, with the convective part
        # rescaled by χ.
        T_s = state.T_surface
        @test state.q_applied ≈ 2e4 + h * χ * (T_inf - T_s) rtol = 1e-10
        @test state.q_applied ≈ state.q_conducted rtol = 1e-8
        @test state.q_conducted ≈ k_eff * (T_s - T_cell) / Δz_half rtol = 1e-10
        # Transpiration stays off (S_conv-only convention).
        @test state.q_transpired == 0.0

        # Blowing suppresses convective cooling → hotter surface.
        state_off =
            solve_boundary_temperature(T_cell, k_eff, Δz_half, bc, mat, J_zero, T_inf, t, ε)
        @test state.T_surface > state_off.T_surface
        @test state_off.blowing_mcp == 0.0

        # Zero mass flux with blowing enabled is the exact no-blowing state.
        state_zero = solve_boundary_temperature(
            T_cell,
            k_eff,
            Δz_half,
            bc,
            mat,
            J_zero,
            T_inf,
            t,
            ε;
            use_blowing = true,
            J_blowing = J_zero,
        )
        @test state_zero.T_surface === state_off.T_surface
        @test state_zero.q_applied === state_off.q_applied

        # Net suction (inflow) enhances cooling → colder surface than off.
        state_suction = solve_boundary_temperature(
            T_cell,
            k_eff,
            Δz_half,
            bc,
            mat,
            J_zero,
            T_inf,
            t,
            ε;
            use_blowing = true,
            J_blowing = (0.0, -J_gas),
        )
        @test state_suction.T_surface < state_off.T_surface
        @test state_suction.blowing_mcp < 0.0
    end

    # ========================================================================
    @testset "AD: T_s derivatives match finite differences" begin
        mat = _blowing_material()
        h, T_inf = 20.0, 300.0
        bc = HeatFluxBC(2e4) + ConvectiveBC(h = h, T_inf = T_inf)
        k_eff, Δz_half, t, ε = 0.2, 5e-4, 0.0, 0.9

        function _Ts(T_cell::Real, J_gas::Real)
            st = BoundaryState()
            Jz = (zero(J_gas), zero(J_gas))
            result = solve_boundary_temperature!(
                st,
                T_cell,
                k_eff,
                Δz_half,
                bc,
                mat,
                Jz,
                T_inf,
                t,
                ε;
                use_blowing = true,
                J_blowing = (zero(J_gas), J_gas),
            )
            return result.T_s
        end

        # ∂T_s/∂T_cell
        d_ad = ForwardDiff.derivative(Tc -> _Ts(Tc, 0.02), 600.0)
        δ = 1e-3
        d_fd = (_Ts(600.0 + δ, 0.02) - _Ts(600.0 - δ, 0.02)) / (2δ)
        @test d_ad ≈ d_fd rtol = 1e-6

        # ∂T_s/∂ṁ″ (the blowing channel itself)
        dJ_ad = ForwardDiff.derivative(J -> _Ts(600.0, J), 0.02)
        δJ = 1e-6
        dJ_fd = (_Ts(600.0, 0.02 + δJ) - _Ts(600.0, 0.02 - δJ)) / (2δJ)
        @test abs(dJ_ad) > 0
        @test dJ_ad ≈ dJ_fd rtol = 1e-5
    end

    # ========================================================================
    @testset "Structured Jacobian: strict agreement with blowing" begin
        function _build(; darcy::Bool, use_ale::Bool)
            mat = _blowing_material(darcy = darcy)
            mesh = create_uniform_mesh(0.01, 12, Val(2))
            bcs = BoundaryConditionSet(
                thermal = (
                    top = HeatFluxBC(50e3) + ConvectiveBC(h = 10.0, T_inf = 300.0),
                    bottom = AdiabaticBC(),
                ),
                mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
            )
            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bcs,
                tspan = (0.0, 100.0),
            )
            def, ws, u =
                initialize_state(problem; use_ale = use_ale, use_blowing_correction = true)
            n = ws.layout.n_cells
            # Hot state with free gas so the boundary mass flux (and hence B)
            # is genuinely nonzero.
            for i = 1:n
                u[i] = 1000.0
                u[n+i] = 400.0          # ξ_solid
                u[2n+i] = 8.0           # ξ_gas
            end
            return def, ws, u
        end

        for (name, darcy, use_ale) in (
            ("eulerian Fickian", false, false),
            ("eulerian Darcy-Fick", true, false),
            ("ALE Fickian", false, true),
        )
            def, ws, u = _build(darcy = darcy, use_ale = use_ale)
            ref_spec = JacobianSpec(DenseAD_ForwardDiff(), def, ws)
            spec = JacobianSpec(Structured(), def, ws)
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
            @test report.mode == :strict
            @test report.passed
        end
    end

    # ========================================================================
    @testset "integration: solve() with blowing" begin
        function _solve(; blowing::Bool)
            mat = _blowing_material()
            mesh = create_uniform_mesh(
                0.005,
                8,
                Val(2);
                T_initial = 600.0,
                ξ_initial = (1000.0, 0.0),
            )
            bcs = BoundaryConditionSet(
                thermal = (
                    top = HeatFluxBC(50e3) + ConvectiveBC(h = 15.0, T_inf = 300.0),
                    bottom = AdiabaticBC(),
                ),
                mass = (top = DiffusiveMassBC(0.01), bottom = ImpermeableBC()),
            )
            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bcs,
                tspan = (0.0, 5.0),
            )
            return solve(
                problem;
                use_blowing_correction = blowing,
                saveat = [0.0, 2.5, 5.0],
                progress = false,
            )
        end

        sol_off = _solve(blowing = false)
        sol_on = _solve(blowing = true)
        @test sol_off.retcode == :Success
        @test sol_on.retcode == :Success

        # Outgassing under convective cooling: blowing suppresses the loss,
        # so the blowing-on surface ends hotter and loses mass faster.
        @test sol_on.surface_T[end] > sol_off.surface_T[end]
        @test sol_on.mass_loss_rate[end] > sol_off.mass_loss_rate[end]

        # Saved flux decomposition stays consistent with the boundary ledger:
        # at Newton convergence net = conduction (transpiration off), with the
        # convective component blowing-rescaled inside `net`.
        for sol in (sol_on, sol_off)
            net = sol.extras.surface_heat_flux_net[end]
            cond = sol.extras.surface_heat_flux_conduction[end]
            @test net ≈ cond rtol = 1e-6
            @test sol.extras.surface_heat_flux_transpiration[end] == 0.0
        end
        # Blowing-on convective loss is strictly smaller in magnitude than the
        # uncorrected h·(T_s − T_∞) at its own (hotter) surface temperature.
        h, T_inf = 15.0, 300.0
        T_s_on = sol_on.surface_T[end]
        conv_on = sol_on.extras.surface_heat_flux_convection[end]
        @test conv_on < h * (T_s_on - T_inf)
        @test conv_on > 0.0   # still a loss
    end
end
