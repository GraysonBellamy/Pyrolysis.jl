# =============================================================================
# V-CPL-05 — TGA-to-gasification consistency ladder (0D→1D exactness chain)
# =============================================================================
#
# MaCFP-style internal-consistency chain: the same first-order kinetics run
# (a) 0D isothermal — covered exactly by the V-RX cases;
# (b) 0D constant heating rate — exact reference ξ(t) = ξ₀e^(−(A/β)I(T)) with
#     the temperature integral I computed to high precision (exact/tga.jl,
#     self-tested against the asymptotic expansion);
# (c) 1D slabs driven by a Dirichlet ramp T = T₀ + βt on BOTH faces, with a
#     thickness ladder L → 0: the volume-mean solid mass must converge to the
#     0D TGA curve at second order in L (the interior lags the faces by the
#     quasi-steady parabola ΔT ≈ β(L/2)²/2α ∝ L², which maps onto the mass
#     curve through dξ/dT). Catches BC/volume-coupling contamination of the
#     kinetics that no pure 0D test can.
#
# The depletion limiter is disabled: its tanh roll-off below ~1 kg/m³ holds a
# residual tail of order 0.1% of ξ₀ — the size of the thinnest rung's entire
# error — which would corrupt both the error metric and the observed order.
# h = 0 keeps the ramp exact (no reaction self-heating), everything is
# impermeable and Eulerian (held gas, no volume change), so temperature is
# the ONLY coupling between the 1D solve and the kinetics.

function _vcpl05_params(; β = 0.5)   # default 0.5 K/s (30 K/min)
    return (
        ξ0 = 1000.0,
        A = 1.0e10,
        E = 130.0e3,
        β = β,
        T0 = 300.0,
        t_end = 400.0 / β,  # T reaches 700 K; conversion complete by ~660 K
        k = 0.2,
        c = 1500.0,
    )
end

function _vcpl05_problem(L; n_cells, β = 0.5)
    P = _vcpl05_params(β = β)
    solid = SolidComponent(:A, ρ = P.ξ0, c = P.c, k = P.k)
    gas = GaseousComponent(:gas, M = 0.05, c = P.c, k = 0.02, λ = 1e-8)
    rxn = Reaction(:decomp, 1 => 2; A = P.A, E = P.E, h = 0.0)
    mat = Material(
        name = :VCPL05,
        components = (solid, gas),
        reactions = (rxn,),
        depletion_limiter = DepletionLimiter(enabled = false),
    )
    mesh =
        create_uniform_mesh(L, n_cells, Val(2); T_initial = P.T0, ξ_initial = (P.ξ0, 0.0))
    ramp = t -> P.T0 + P.β * t
    bc_set = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(ramp), bottom = DirichletThermalBC(ramp)),
    )
    return PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> P.T0,
        ξ_initial = [z -> P.ξ0, z -> 0.0],
        tspan = (0.0, P.t_end),
        output_times = collect(range(0.0, P.t_end, length = 81)),
    )
end

# Max-over-time deviation of the volume-mean solid mass from the 0D TGA curve
# for parameter set `P` — the ladder metric shared by the QoIs, the report
# curves, and the heating-rate sweep.
function _vcpl05_ladder_error(sol, P)
    errmax = 0.0
    for (j, t) in enumerate(sol.t)
        ξ̄ = sum(@view sol.ξ[1, :, j]) / size(sol.ξ, 2)
        ξex = tga_xi(t; ξ0 = P.ξ0, A = P.A, E = P.E, β = P.β, T0 = P.T0)
        errmax = max(errmax, abs(ξ̄ - ξex) / P.ξ0)
    end
    return errmax
end

function _verification_case_vcpl05()
    P = _vcpl05_params()
    L_ladder = [2e-3, 1e-3, 0.5e-3]   # build() provides the first rung

    build(n_cells) = _vcpl05_problem(L_ladder[1]; n_cells = n_cells)

    ladder_error(sol) = _vcpl05_ladder_error(sol, P)

    function qois(sol)
        out = QoIResult[]
        errs = [ladder_error(sol)]
        for L in L_ladder[2:end]
            s = solve(
                _vcpl05_problem(L; n_cells = 16);
                progress = false,
                abstol = 1e-10,
                reltol = 1e-8,
            )
            push!(errs, ladder_error(s))
        end

        # 1. Thinnest rung: 1D ≈ 0D to the quasi-steady-lag level.
        push!(
            out,
            QoIResult(
                name = "thin-limit deviation from 0D TGA (L = $(L_ladder[end]*1e3) mm)",
                value = errs[end],
                error = errs[end],
                tol = 3e-3,  # 3× observed: 1.0e-3 (n=16)
                comment = "max-over-t of |⟨ξ⟩ − ξ_0D|/ξ₀; lag ΔT ≈ β(L/2)²/2α ≈ 0.12 K",
            ),
        )

        # 2. Convergence order in L (the ladder itself): expect 2.
        rates = observed_orders(errs, L_ladder)
        for (i, r) in enumerate(rates)
            push!(
                out,
                QoIResult(
                    name = "L-convergence rate, rung $(i) → $(i + 1)",
                    value = r,
                    exact = 2.0,
                    error = abs(r - 2.0),
                    tol = 0.1,  # observed rates 2.021, 2.006
                    comment = "errors $(round.(errs, sigdigits = 3)) at L = $(L_ladder)",
                ),
            )
        end

        # 3. Completeness: conversion finished on every rung (guards against
        # asserting agreement between two curves that simply never reacted).
        ξ̄_end = sum(@view sol.ξ[1, :, length(sol.t)]) / size(sol.ξ, 2)
        push!(
            out,
            QoIResult(
                name = "complete conversion at t_end (coarsest rung)",
                value = ξ̄_end / P.ξ0,
                exact = 0.0,
                error = abs(ξ̄_end) / P.ξ0,
                tol = 1e-6,
                comment = "exact residual e^(−(A/β)I(700 K)) ≈ e⁻¹¹⁴",
            ),
        )
        return out
    end

    return VerificationCase(
        id = "V-CPL-05",
        title = "TGA-to-gasification consistency ladder (1D → 0D thin limit)",
        tags = [:kinetics, :consistency_ladder, :thin_limit, :bc_coupling],
        references = [
            "MaCFP condensed-phase guidelines (2020) — TGA target simulations",
            "Verification Suite Plan 2026-07-06, V-CPL-05",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        qois = qois,
        ci = CIConfig(n_cells = 16, budget_s = 25.0),
        report = (
            qoi_curves = (:mass_curves_overlay, :error_vs_L),
            # The report thickness ladder is in L, not n_cells (the CI ladder
            # extended by one rung on each side), so it lives in `curves`
            # rather than `n_cells_ladder`; the O(L²) reference line is
            # anchored at the thinnest rung.
            curves = _ -> begin
                Ls = [4e-3, 2e-3, 1e-3, 0.5e-3, 0.25e-3]
                sols = [
                    solve(
                        _vcpl05_problem(L; n_cells = 16);
                        progress = false,
                        abstol = 1e-10,
                        reltol = 1e-8,
                    ) for L in Ls
                ]
                errs = [_vcpl05_ladder_error(s, P) for s in sols]
                mass_series = [
                    CurveSeries(
                        label = "L = $(L * 1e3) mm",
                        x = s.t,
                        y = [
                            sum(@view s.ξ[1, :, j]) / size(s.ξ, 2) / P.ξ0 for
                            j in eachindex(s.t)
                        ],
                        exact = t ->
                            tga_xi(t; ξ0 = P.ξ0, A = P.A, E = P.E, β = P.β, T0 = P.T0) / P.ξ0,
                    ) for (L, s) in zip(Ls, sols)
                ]
                [
                    ReportCurve(
                        name = :mass_curves_overlay,
                        title = "V-CPL-05: volume-mean mass curves vs the 0D TGA " *
                                "reference, L ladder",
                        xlabel = "t (s)",
                        ylabel = "⟨ξ_s⟩/ξ₀",
                        series = mass_series,
                    ),
                    ReportCurve(
                        name = :error_vs_L,
                        title = "V-CPL-05: thin-limit convergence (line: O(L²) " *
                                "reference through the thinnest rung)",
                        xlabel = "L (m)",
                        ylabel = "max-over-t |⟨ξ⟩ − ξ_0D|/ξ₀",
                        xscale = :log10,
                        yscale = :log10,
                        series = [
                            CurveSeries(
                                label = "deviation from 0D",
                                x = Ls,
                                y = errs,
                                exact = L -> errs[end] * (L / Ls[end])^2,
                            ),
                        ],
                    ),
                ]
            end,
            # Heating-rate insensitivity at the thinnest CI rung: the 1D→0D
            # deviation must stay at the quasi-steady-lag level (∝ β·L²/α
            # through dξ/dT) across a 10× β range.
            param_sweep = (
                name = :beta,
                values = [10.0 / 60, 0.5, 100.0 / 60],
                build = (n, β) -> _vcpl05_problem(0.5e-3; n_cells = 16, β = β),
                measure = (sol, β) ->
                    (thin_limit_err = _vcpl05_ladder_error(sol, _vcpl05_params(β = β)),),
            ),
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl05)
