# =============================================================================
# V-HC-01 — Semi-infinite solid, convective surface BC (Biot-coupled erfc)
# =============================================================================
#
# The one verification case shared by all three reference codes
# (FDS `ht3d_slab`, Gpyro §3.4.1, ThermaKin Eq. 25), in Gpyro's variant:
# constant absorbed radiant flux q″ = 25 kW/m² plus convective exchange
# h = 20 W/m²K with T∞ = T₀ = 300 K; k = 0.2 W/(m·K), ρ = 1000 kg/m³,
# c = 1400 J/(kg·K), so profiles at t = 30/60/90/180 s are directly comparable
# to the dissertation figures.
#
# Because the combined surface BC is linear in T_s, absorbed-flux + convection
# is exactly convection toward T∞_eff = T∞ + q″/h = 1550 K, and the Carslaw &
# Jaeger §2.7 erfc solution (exact/semi_infinite.jl) applies unchanged.
#
# Coverage beyond the existing Dirichlet erf test: `HeatFluxBC + ConvectiveBC`
# enter the surface Newton solve, so the surface temperature is an unknown
# coupled to the interior — the place surface-BC bugs live.
#
# Domain: L = 5 cm ≫ thermal penetration (erfc(η) ≈ 3e-12 at the back face at
# t = 180 s), so the finite slab is semi-infinite to well below the QoI
# tolerances; the back-face-rise QoI enforces this assumption explicitly.

function _verification_case_vhc01()
    ρ, c, k = 1000.0, 1400.0, 0.2
    α = k / (ρ * c)
    q_abs, h, T0 = 25.0e3, 20.0, 300.0
    Tinf_eff = T0 + q_abs / h
    L = 0.05
    times = [30.0, 60.0, 90.0, 180.0]

    T_exact(x, t) =
        T_semiinf_convective(x, t; T0 = T0, Tinf = Tinf_eff, h = h, k = k, α = α)

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = ρ, c = c, k = k)
        mat = Material(name = :VHC01, components = (solid,), reactions = ())
        mesh = create_uniform_mesh(L, n_cells, Val(1); T_initial = T0, ξ_initial = (ρ,))
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(q_abs) + ConvectiveBC(h = h, T_inf = T0),
                bottom = AdiabaticBC(),
            ),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T0,
            ξ_initial = [z -> ρ],
            tspan = (0.0, times[end]),
            output_times = vcat(0.0, times),
        )
    end

    function qois(sol)
        out = QoIResult[]
        z = get_cell_positions(sol)
        for (j, t) in enumerate(times)
            Tn = @view sol.T[:, j+1]  # output_times = [0; times]
            err = maximum(abs(Tn[i] - T_exact(L - z[i], t)) for i in eachindex(z))
            push!(
                out,
                QoIResult(
                    name = "L∞ T-profile error, t=$(t) s",
                    value = err,
                    error = err,
                    tol = 0.75,  # 3× observed max: 0.24 K at t=30 s (n=100)
                    comment = "observed 0.24/0.21/0.18/0.14 K at t=30/60/90/180 s",
                ),
            )
        end
        for (j, t) in enumerate(times)
            Ts_num = sol.surface_T[j+1]
            Ts_exact = T_exact(0.0, t)
            push!(
                out,
                QoIResult(
                    name = "surface temperature, t=$(t) s",
                    value = Ts_num,
                    exact = Ts_exact,
                    error = abs(Ts_num - Ts_exact),
                    tol = 1.5,  # 3× observed max: 0.50 K at t=30 s (n=100)
                    comment = "surface Newton solve vs exact θ_s = 1 − erfcx(h√(αt)/k); " *
                              "observed 0.50/0.26/0.16/0.06 K at t=30/60/90/180 s",
                ),
            )
        end
        ΔT_back = sol.T[1, end] - T0
        push!(
            out,
            QoIResult(
                name = "back-face rise (semi-infinite validity)",
                value = ΔT_back,
                exact = 0.0,
                error = abs(ΔT_back),
                tol = 1e-6,  # observed 1.2e-9 K; exact rise at x=L, t=180 s is O(1e-9) K
                comment = "guards the semi-infinite assumption of the reference solution",
            ),
        )
        return out
    end

    # L2 error of the t = 60 s profile (output index 3) for the order check.
    function l2_error_t60(sol)
        z = get_cell_positions(sol)
        Tn = @view sol.T[:, 3]
        Δz = L / length(z)
        return sqrt(sum((Tn[i] - T_exact(L - z[i], 60.0))^2 * Δz for i in eachindex(z)) / L)
    end

    return VerificationCase(
        id = "V-HC-01",
        title = "semi-infinite solid, convective surface BC (erfc)",
        tags = [:heat_conduction, :convective_bc, :surface_newton],
        references = [
            "Carslaw & Jaeger (1959), Conduction of Heat in Solids, §2.7",
            "Lautenberger dissertation (2007) §3.4.1 (parameters)",
            "FDS Verification Guide `ht3d_slab`; ThermaKin (FSS 9:1141) Eq. 25",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        qois = qois,
        ci = CIConfig(
            n_cells = 100,
            budget_s = 20.0,
            convergence = ConvergenceSpec(
                # observed rates 1.991, 1.998
                n_cells = [40, 80, 160],
                error_fn = l2_error_t60,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
            ),
        ),
        report = (
            n_cells_ladder = [25, 50, 100, 133, 200, 400, 500],
            profile_times = times,
            qoi_curves = (:T_profiles, :surface_T_history),
            solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
            curves = sol -> begin
                x = L .- get_cell_positions(sol)
                # FDS 6.11.0 overlay (reference/fds/custom/): our deck
                # vhc01_semi_infinite.fds solves this exact case via the same
                # Robin-BC equivalence (pure convection toward T∞_eff = 1550 K,
                # ε = 0). devc depths from the exposed face; T in °C.
                fds = reference_columns(
                    "fds/custom/vhc01_semi_infinite_devc.csv",
                    "Time",
                    "Front",
                    "5 mm",
                    "10 mm",
                    "15 mm",
                    "20 mm",
                    "30 mm",
                )
                fds_depths = [0.0, 0.005, 0.010, 0.015, 0.020, 0.030]
                fds_T = [c .+ 273.15 for c in fds[2:end]]
                fds_row(t) = argmin(abs.(fds[1] .- t))
                # Gpyro 0.8200 overlay (reference/gpyro/): dissertation §3.4.1
                # deck (identical parameters — our case IS Gpyro's variant),
                # Δz = 0.1 mm, dt = 0.1 s. Long-format (t, x, T) profiles.
                gp_t, gp_x, gp_T = reference_columns(
                    "gpyro/vhc01_result.csv",
                    "t_s",
                    "x_m",
                    "T_K";
                    header_rows = 1,
                )
                gp_ts, gp_Ts = reference_columns(
                    "gpyro/vhc01_surfaceT.csv",
                    "t_s",
                    "T_surface_K";
                    header_rows = 1,
                )
                [
                    ReportCurve(
                        name = :T_profiles,
                        title = "V-HC-01: in-depth temperature profiles",
                        xlabel = "depth x (m)",
                        ylabel = "T (K)",
                        series = vcat(
                            [
                                CurveSeries(
                                    label = "t = $(t) s",
                                    x = x,
                                    y = sol.T[:, j+1],
                                    exact = xx -> T_exact(xx, t),
                                ) for (j, t) in enumerate(times)
                            ],
                            [
                                CurveSeries(
                                    label = "FDS 6.11.0, t = $(t) s",
                                    crosscode = true,
                                    x = fds_depths,
                                    y = [fds_T[i][fds_row(t)] for i = 1:6],
                                    exact = xx -> T_exact(xx, t),
                                ) for t in times
                            ],
                            [
                                let rows = findall(τ -> abs(τ - t) < 1e-6, gp_t)
                                    CurveSeries(
                                        label = "Gpyro 0.8200, t = $(t) s",
                                        crosscode = true,
                                        x = gp_x[rows],
                                        y = gp_T[rows],
                                        exact = xx -> T_exact(xx, t),
                                    )
                                end for t in times
                            ],
                        ),
                    ),
                    ReportCurve(
                        name = :surface_T_history,
                        title = "V-HC-01: surface temperature history",
                        xlabel = "t (s)",
                        ylabel = "T_s (K)",
                        series = [
                            CurveSeries(
                                label = "surface Newton solve",
                                x = sol.t[2:end],  # exact form is 0/0 at t = 0
                                y = sol.surface_T[2:end],
                                exact = t -> T_exact(0.0, t),
                            ),
                            # devc averages over the 0.2 s output interval
                            CurveSeries(
                                label = "FDS 6.11.0",
                                crosscode = true,
                                x = fds[1][2:end] .- 0.1,
                                y = fds_T[1][2:end],
                                exact = t -> T_exact(0.0, t),
                            ),
                            CurveSeries(
                                label = "Gpyro 0.8200",
                                crosscode = true,
                                x = gp_ts[2:end],
                                y = gp_Ts[2:end],
                                exact = t -> T_exact(0.0, t),
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol -> begin
                z = get_cell_positions(sol)
                e = [sol.T[i, 3] - T_exact(L - z[i], 60.0) for i in eachindex(z)]
                (L2 = sqrt(sum(abs2, e) / length(e)), Linf = maximum(abs, e))
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vhc01)
