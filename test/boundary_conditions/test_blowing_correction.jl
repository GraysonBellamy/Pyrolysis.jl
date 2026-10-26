"""
Blowing correction tests: the film factor, the face interface of the surface
Newton, and the `CouetteFilm` model end to end.

1. `blowing_factor` — limits, branch-seam continuity (value + derivative),
   large-φ overflow safety under ForwardDiff, suction enhancement.
2. `evaluate_bc_flux_and_derivative` face dispatch — only convective BCs
   with a blowing model rescale; `NoBlowing`, radiative and flux BCs leave the
   flux alone; combined BCs forward the face to each member.
3. `CouetteFilm` — φ = 2λ·ṁ″·c_p,g/h, each member its own λ, the diagnostic
   decomposition applies the same factor.
4. Surface Newton with blowing — converged balance, `blowing_mcp` and
   `blowing_ratio` diagnostics, off/on ordering, exact no-outflow identity.
5. AD cleanliness — ForwardDiff derivatives of `T_s` w.r.t. `T_cell` and the
   blowing mass flux match central finite differences.
6. Integration — the finite-volume core's `solve` with a film model raises
   the surface temperature under convective cooling, and the saved flux
   decomposition stays consistent with the boundary ledger.

The property-ratio and pooling models are checked against their reference
values in `test_blowing_models.jl`.
"""

using Test
using ForwardDiff
using Pyrolysis
using Pyrolysis.Materials
using Pyrolysis.BoundaryConditions
using Pyrolysis.BoundaryConditions:
    BoundaryState,
    blowing_factor,
    BlowingFace,
    blowing_state,
    blowing_ratio,
    solve_boundary_temperature!,
    evaluate_bc_flux_and_derivative,
    compute_boundary_heat_fluxes,
    _blowing_conductance
using Pyrolysis.Materials: SURFACE_ABSORPTION, NO_RADIATION
using Pyrolysis.States: trial_temperature

# ----------------------------------------------------------------------------
# Shared fixtures: 2-component solid + gas, constant cp_g so the blowing
# conductance is exactly J_gas · c_g and every balance is closed-form.
# ----------------------------------------------------------------------------
const _CP_GAS = 1000.0
function _blowing_material(; darcy::Bool = false)
    solid_kwargs = darcy ? (; κ = 1e-14) : (;)
    solid = SolidComponent(
        :solid;
        φ = 0.3,
        λ = 1e-5,
        ρ = 1000.0,
        c = 1000.0,
        k = 0.2,
        ε = 0.9,
        solid_kwargs...,
    )
    gas = GaseousComponent(:gas; M = 0.050, c = _CP_GAS, k = 0.02)
    rxn = Reaction(:decomp, 1 => 2, A = 1e5, E = 5e4, h = 1e5)
    return Material(name = :blowing_test, components = (solid, gas), reactions = (rxn,))
end

# Reference factor straight from the definition (only valid away from 0/∞).
_f_exact(B) = B / (exp(B) - 1)

# The surface balance at the trial clamp of each iterate, as the residual
# evaluates it; returns the diagnostic state.
function _blowing_surface(T_cell, k_eff, Δz_half, bc, mat, J, t, ε)
    state = BoundaryState()
    solve_boundary_temperature!(
        state,
        T_cell,
        k_eff,
        Δz_half,
        bc,
        mat,
        J,
        t,
        ε;
        trial = trial_temperature,
    )
    return state
end

@testset "Blowing correction" begin

    # ========================================================================
    @testset "blowing_factor: values and limits" begin
        @test blowing_factor(0.0) === 1.0
        # Hand-checked arithmetic: B = 4 → h/h₀ = 0.0746…
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
    @testset "evaluate_bc_flux_and_derivative: face dispatch" begin
        mat = _blowing_material()
        T_s, t, ε = 600.0, 0.0, 0.9
        h, T_inf = 20.0, 300.0
        J_gas = 0.02                                   # kg/(m²·s) outward
        face = BlowingFace((0.0, J_gas), T_s, mat)
        face0 = BlowingFace((0.0, 0.0), T_s, mat)
        mcp = J_gas * _CP_GAS                          # B = 1
        @test _blowing_conductance(mat, face.ṁ, T_s) == mcp
        χ = blowing_factor(mcp / h)

        conv = ConvectiveBC(h = h, T_inf = T_inf, blowing = CouetteFilm())
        q, dq = evaluate_bc_flux_and_derivative(conv, T_s, t, ε, NO_RADIATION, face)
        @test q ≈ h * χ * (T_inf - T_s) rtol = 1e-14
        @test dq ≈ -h * χ rtol = 1e-14

        # Zero outflow reproduces the uncorrected flux exactly.
        q0, dq0 = evaluate_bc_flux_and_derivative(conv, T_s, t, ε, NO_RADIATION, face0)
        @test q0 === h * (T_inf - T_s)
        @test dq0 === -h

        # The default model leaves the coefficient alone whatever the outflow.
        bare = ConvectiveBC(h = h, T_inf = T_inf)
        @test bare.blowing === NoBlowing()
        @test ConvectiveBC(h, T_inf).blowing === NoBlowing()
        qb, dqb = evaluate_bc_flux_and_derivative(bare, T_s, t, ε, NO_RADIATION, face)
        @test qb === q0
        @test dqb === dq0

        # Non-convective BCs ignore the face.
        rad = RadiativeBC(T_inf = T_inf)
        qr, dqr = evaluate_bc_flux_and_derivative(rad, T_s, t, ε, NO_RADIATION, face)
        qr5, dqr5 = evaluate_bc_flux_and_derivative(rad, T_s, t, ε, NO_RADIATION)
        @test qr == qr5 && dqr == dqr5

        flux = HeatFluxBC(2e4)
        qf, _ = evaluate_bc_flux_and_derivative(flux, T_s, t, ε, NO_RADIATION, face)
        @test qf == 2e4

        # Combined: convective child rescaled, others untouched.
        combo = flux + conv + rad
        qc, dqc = evaluate_bc_flux_and_derivative(combo, T_s, t, ε, NO_RADIATION, face)
        @test qc ≈ qf + q + qr rtol = 1e-14
        @test dqc ≈ dq + dqr rtol = 1e-14

        # A user rule is applied as given.
        custom = ConvectiveBC(
            h = h,
            T_inf = T_inf,
            blowing = FunctionBlowing((f, hh, Ti) -> 1 / (1 + sum(f.ṁ) / 0.01)),
        )
        qu, _ = evaluate_bc_flux_and_derivative(custom, T_s, t, ε, NO_RADIATION, face)
        @test qu ≈ h / 3 * (T_inf - T_s) rtol = 1e-14
    end

    # ========================================================================
    @testset "CouetteFilm: film constant" begin
        mat = _blowing_material()
        T_s, t, ε = 600.0, 0.0, 0.9
        h, T_inf = 20.0, 300.0
        face = BlowingFace((0.0, 0.02), T_s, mat)
        mcp = 0.02 * _CP_GAS

        @test CouetteFilm().λ == 0.5
        @test CouetteFilm(0.31).λ == 0.31
        @test CouetteFilm(λ = 0.31).λ == 0.31
        @test_throws ArgumentError CouetteFilm(0.0)
        @test_throws ArgumentError CouetteFilm(NaN)

        # The Couette default is bit-identical to the plain mcp/h argument.
        s = blowing_state(CouetteFilm(), face, h, T_inf)
        @test s.χ === blowing_factor(mcp / h)
        @test s.φ === mcp / h
        @test s.ψ == 1 && s.Y == 0 && s.Λ == 1 && s.states == 1

        # φ = 2λ·mcp/h; a smaller λ weakens the reduction of h.
        λ = 0.31
        conv = ConvectiveBC(h = h, T_inf = T_inf, blowing = CouetteFilm(λ))
        χ = blowing_factor(2λ * mcp / h)
        q, dq = evaluate_bc_flux_and_derivative(conv, T_s, t, ε, NO_RADIATION, face)
        @test q ≈ h * χ * (T_inf - T_s) rtol = 1e-14
        @test dq ≈ -h * χ rtol = 1e-14
        q_def, _ = evaluate_bc_flux_and_derivative(
            ConvectiveBC(h = h, T_inf = T_inf, blowing = CouetteFilm()),
            T_s,
            t,
            ε,
            NO_RADIATION,
            face,
        )
        @test abs(q) > abs(q_def)

        # The diagnostic decomposition applies the same factor.
        dec = compute_boundary_heat_fluxes(conv, T_s, t, ε, NO_RADIATION, face)
        @test dec.net ≈ q rtol = 1e-14
        @test dec.convection ≈ -q rtol = 1e-14

        # Inside a combined BC each convective member keeps its own model.
        combo = ConvectiveBC(h = h, T_inf = T_inf, blowing = CouetteFilm()) + conv
        qc, _ = evaluate_bc_flux_and_derivative(combo, T_s, t, ε, NO_RADIATION, face)
        @test qc ≈ q_def + q rtol = 1e-14
    end

    # ========================================================================
    @testset "surface Newton with blowing: converged balance" begin
        mat = _blowing_material()
        h, T_inf = 20.0, 300.0
        bc = HeatFluxBC(2e4) + ConvectiveBC(h = h, T_inf = T_inf, blowing = CouetteFilm())
        bc_bare = HeatFluxBC(2e4) + ConvectiveBC(h = h, T_inf = T_inf)
        T_cell, k_eff, Δz_half, t, ε = 600.0, 0.2, 5e-4, 0.0, 0.9
        J_gas = 0.02                                   # kg/(m²·s) outward
        J_zero = (0.0, 0.0)
        J_blow = (0.0, J_gas)

        state = _blowing_surface(T_cell, k_eff, Δz_half, bc, mat, J_blow, t, ε)
        @test state.converged

        # Constant cp_g → mcp is exactly J·cp, independent of T_s.
        mcp = J_gas * _CP_GAS
        @test state.blowing_mcp ≈ mcp rtol = 1e-12
        χ = blowing_factor(mcp / h)
        @test state.blowing_ratio ≈ χ rtol = 1e-12
        @test state.pooling_factor == 1.0
        @test state.face_fraction == 0.0

        # Newton balance: q_applied = q_conducted, with the convective part
        # rescaled by χ.
        T_s = state.T_surface
        @test state.q_applied ≈ 2e4 + h * χ * (T_inf - T_s) rtol = 1e-10
        @test state.q_applied ≈ state.q_conducted rtol = 1e-8
        @test state.q_conducted ≈ k_eff * (T_s - T_cell) / Δz_half rtol = 1e-10

        # Blowing suppresses convective cooling → hotter surface than the
        # bare coefficient under the same outflow.
        state_off = _blowing_surface(T_cell, k_eff, Δz_half, bc_bare, mat, J_blow, t, ε)
        @test state.T_surface > state_off.T_surface
        @test state_off.blowing_mcp ≈ mcp rtol = 1e-12   # the conductance is reported
        @test state_off.blowing_ratio == 1.0             # but not applied

        # Zero mass flux with a model is the exact no-blowing state.
        state_zero = _blowing_surface(T_cell, k_eff, Δz_half, bc, mat, J_zero, t, ε)
        state_bare0 = _blowing_surface(T_cell, k_eff, Δz_half, bc_bare, mat, J_zero, t, ε)
        @test state_zero.T_surface === state_bare0.T_surface
        @test state_zero.q_applied === state_bare0.q_applied
        @test state_zero.blowing_ratio == 1.0

        # Net suction (inflow) enhances cooling → colder surface than bare.
        state_suction =
            _blowing_surface(T_cell, k_eff, Δz_half, bc, mat, (0.0, -J_gas), t, ε)
        @test state_suction.T_surface < state_off.T_surface
        @test state_suction.blowing_mcp < 0.0
        @test state_suction.blowing_ratio > 1.0
    end

    # ========================================================================
    @testset "AD: T_s derivatives match finite differences" begin
        mat = _blowing_material()
        h, T_inf = 20.0, 300.0
        bc = HeatFluxBC(2e4) + ConvectiveBC(h = h, T_inf = T_inf, blowing = CouetteFilm())
        k_eff, Δz_half, t, ε = 0.2, 5e-4, 0.0, 0.9

        function _Ts(T_cell::Real, J_gas::Real)
            st = BoundaryState()
            result = solve_boundary_temperature!(
                st,
                T_cell,
                k_eff,
                Δz_half,
                bc,
                mat,
                (zero(J_gas), J_gas),
                t,
                ε;
                trial = trial_temperature,
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
    @testset "integration: solve() with blowing" begin
        function _solve(model)
            mat = _blowing_material()
            mesh = create_uniform_mesh(0.005, 8)
            bcs = BoundaryConditionSet(
                thermal = (
                    top = HeatFluxBC(50e3) +
                          ConvectiveBC(h = 15.0, T_inf = 300.0, blowing = model),
                    bottom = AdiabaticBC(),
                ),
                mass = (
                    top = AmbientGasBC(GasReservoir(P = 0.0, T = 300.0)),
                    bottom = ImpermeableBC(),
                ),
            )
            problem = PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bcs,
                T_initial = z -> 300.0,
                ξ_initial = [z -> 1000.0, z -> 0.0],
                tspan = (0.0, 5.0),
            )
            return Pyrolysis.Integration.solve(problem; saveat = [0.0, 2.5, 5.0])
        end

        sol_off = _solve(NoBlowing())
        sol_on = _solve(CouetteFilm())
        @test sol_off.retcode == :Success
        @test sol_on.retcode == :Success

        # Outgassing under convective cooling: blowing suppresses the loss,
        # so the blowing-on surface ends hotter and loses mass faster.
        @test sol_on.surface_T[end] > sol_off.surface_T[end]
        @test sol_on.mass_loss_rate[end] > sol_off.mass_loss_rate[end]
        @test all(sol_off.extras.surface_blowing_ratio .== 1.0)
        @test sol_on.extras.surface_blowing_ratio[end] < 1.0

        # Saved flux decomposition stays consistent with the boundary ledger:
        # at Newton convergence net = conduction, with the convective
        # component blowing-rescaled inside `net`.
        for sol in (sol_on, sol_off)
            net = sol.extras.surface_heat_flux_net[end]
            cond = sol.extras.surface_heat_flux_conduction[end]
            @test net ≈ cond rtol = 1e-6
        end
        # Blowing-on convective loss is strictly smaller in magnitude than the
        # uncorrected h·(T_s − T_∞) at its own (hotter) surface temperature.
        h, T_inf = 15.0, 300.0
        T_s_on = sol_on.surface_T[end]
        conv_on = sol_on.extras.surface_heat_flux_convection[end]
        @test conv_on < h * (T_s_on - T_inf)
        @test conv_on > 0.0   # still a loss

        # The pooling model runs through the integrator and records its state.
        sol_pool = _solve(PoolingFilm(NITROGEN; diffusion_volume = (gas = 60.0,)))
        @test sol_pool.retcode == :Success
        ψ = sol_pool.extras.surface_pooling_factor[end]
        @test 0.59 <= ψ <= 1.0
        @test 0.0 <= sol_pool.extras.surface_face_fraction[end] <= 1.0
        @test sol_pool.extras.surface_blowing_ratio[end] < 1.0
    end
end
