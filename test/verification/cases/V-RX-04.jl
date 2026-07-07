# =============================================================================
# V-RX-04 — Thermally-thin radiative cooling + steady radiative–convective
#           surface equilibrium
# =============================================================================
#
# The only closed-form transient involving the T⁴ nonlinearity: a thin slab
# cooling purely by radiation, dT/dt = −K(T⁴−T∞⁴) with K = εσ/(ρcL), exact
# implicit solution t = (G(Tᵢ)−G(T))/K (exact/radiative_cooling.jl, quadrature
# self-tested). Verifies `RadiativeBC` + the surface Newton solve against
# something other than steady algebra.
#
#   a) transient cooling — thin-limit setup: L = 0.1 mm, k = 100 W/mK makes
#      the internal ΔT ≈ q·L/2k ≈ 0.017 K at t = 0 (flux bias ~7e-5 relative),
#      far below the comparison tolerance; compared on the cell-mean
#      temperature, which is what the lumped ODE governs.
#   b) steady companion — HeatFluxBC + RadiativeBC + ConvectiveBC vs the exact
#      algebraic root of q_abs = εσ(Ts⁴−T∞⁴) + h(Ts−T∞): a direct check of
#      the surface-Newton residual (Ts ≈ 732 K, cone-like). Run to 36 thermal
#      time constants; at steady state with an adiabatic back the whole slab
#      must be uniform at Ts.
#
# Uses the code's own STEFAN_BOLTZMANN so the reference cannot drift from the
# implementation's constant.

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        ρ, c, k = 1000.0, 1000.0, 100.0
        L = 1e-4
        ε, T∞, Ti = 0.9, 300.0, 900.0
        K = ε * STEFAN_BOLTZMANN / (ρ * c * L)
        ts = [0.5, 1.0, 2.0, 4.0, 8.0]

        function build(n_cells)
            solid = SolidComponent(:solid, ρ = ρ, c = c, k = k, ε = ε)
            mat = Material(name = :VRX04a, components = (solid,), reactions = ())
            mesh = create_uniform_mesh(
                L,
                n_cells,
                Val(1);
                T_initial = Ti,
                ξ_initial = (ρ,),
            )
            bc_set = BoundaryConditionSet(
                thermal = (top = RadiativeBC(ε = ε, T_inf = T∞), bottom = AdiabaticBC()),
            )
            return PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> Ti,
                ξ_initial = [z -> ρ],
                tspan = (0.0, ts[end]),
                output_times = vcat(0.0, ts),
            )
        end

        function qois(sol)
            out = QoIResult[]
            err = maximum(
                abs(mean(sol.T[:, j+1]) - radiative_cooling_T(t; Ti = Ti, T∞ = T∞, K = K)) for (j, t) in enumerate(ts)
            )
            push!(
                out,
                QoIResult(
                    name = "max mean-T error vs implicit exact",
                    value = err,
                    error = err,
                    tol = 0.012,  # 3× observed 3.9e-3 K (n=4, 2026-07-06)
                    comment = "lumped T⁴ cooling; thin-limit bias ~7e-5 relative in flux",
                ),
            )
            ΔT_int = maximum(
                maximum(sol.T[:, j]) - minimum(sol.T[:, j]) for j = 2:length(sol.t)
            )
            push!(
                out,
                QoIResult(
                    name = "internal ΔT (thin-limit validity)",
                    value = ΔT_int,
                    error = ΔT_int,
                    tol = 0.05,  # a-priori estimate qL/2k ≈ 0.017 K at t=0
                    comment = "guards the lumped assumption of the reference solution",
                ),
            )
            return out
        end

        VerificationCase(
            id = "V-RX-04a",
            title = "thermally-thin radiative cooling (implicit exact)",
            tags = [:radiation, :radiative_bc, :surface_newton, :thermally_thin],
            references = [
                "partial-fraction antiderivative re-derived in exact/radiative_cooling.jl",
                "suite plan V-RX-04 (transient T⁴ closed form)",
            ],
            build = build,
            solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
            qois = qois,
            ci = CIConfig(n_cells = 4, budget_s = 15.0),
            report = (qoi_curves = (:T_history_vs_exact,),),
        )
    end,
)

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        ρ, c, k = 1000.0, 1000.0, 1.0
        L = 0.005
        q_abs, ε, h, T∞ = 20.0e3, 0.85, 15.0, 300.0
        T0 = 300.0
        t_end = 2000.0  # ≈ 36·τ with τ = ρcL/(h + 4εσTs³) ≈ 55 s
        Ts_exact = radiative_convective_equilibrium_Ts(
            q_abs = q_abs,
            ε = ε,
            h = h,
            T∞ = T∞,
            σ = STEFAN_BOLTZMANN,
        )

        function build(n_cells)
            solid = SolidComponent(:solid, ρ = ρ, c = c, k = k, ε = ε)
            mat = Material(name = :VRX04b, components = (solid,), reactions = ())
            mesh = create_uniform_mesh(
                L,
                n_cells,
                Val(1);
                T_initial = T0,
                ξ_initial = (ρ,),
            )
            bc_set = BoundaryConditionSet(
                thermal = (
                    top = HeatFluxBC(q_abs) +
                          RadiativeBC(ε = ε, T_inf = T∞) +
                          ConvectiveBC(h = h, T_inf = T∞),
                    bottom = AdiabaticBC(),
                ),
            )
            return PyrolysisProblem(
                mesh = mesh,
                material = mat,
                bc_set = bc_set,
                T_initial = z -> T0,
                ξ_initial = [z -> ρ],
                tspan = (0.0, t_end),
                output_times = [0.0, t_end / 2, t_end],
            )
        end

        function qois(sol)
            Ts_num = sol.surface_T[end]
            steadiness = abs(sol.surface_T[end] - sol.surface_T[end-1])
            uniformity = maximum(sol.T[:, end]) - minimum(sol.T[:, end])
            return [
                QoIResult(
                    name = "steady surface temperature vs algebraic root",
                    value = Ts_num,
                    exact = Ts_exact,
                    error = abs(Ts_num - Ts_exact),
                    tol = 1e-6,  # observed 3.2e-9 K (2026-07-06): Newton residual is exact
                    comment = "surface-Newton residual check; Ts_exact ≈ " *
                              "$(round(Ts_exact, digits = 2)) K",
                ),
                QoIResult(
                    name = "steadiness (last half-interval drift)",
                    value = steadiness,
                    error = steadiness,
                    tol = 1e-3,  # observed 1.1e-4 K (2026-07-06): integrator save noise, not transient
                    comment = "t_end ≈ 36 thermal time constants",
                ),
                QoIResult(
                    name = "steady uniformity (adiabatic back)",
                    value = uniformity,
                    error = uniformity,
                    tol = 1e-7,  # observed 7.6e-10 K (2026-07-06)
                    comment = "zero steady flux ⇒ uniform slab at Ts",
                ),
            ]
        end

        VerificationCase(
            id = "V-RX-04b",
            title = "steady radiative–convective surface equilibrium",
            tags = [:radiation, :radiative_bc, :convective_bc, :surface_newton, :steady],
            references = ["suite plan V-RX-04 companion (algebraic Ts root)"],
            build = build,
            solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
            qois = qois,
            ci = CIConfig(n_cells = 20, budget_s = 15.0),
            report = (qoi_curves = (:surface_T_history,),),
        )
    end,
)
