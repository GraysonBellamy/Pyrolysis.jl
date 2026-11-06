# =============================================================================
# V-MMS-05 — manufactured solution with radiative transfer on a moving mesh
# =============================================================================
#
# The order-verification case of conduction coupled to radiative transfer in
# the solid. A gray, semi-transparent solid (ε = 0.8, κ = αξ ≈ 200–210 m⁻¹,
# κL ≈ 2) carries a manufactured temperature between 600 and 790 K, so the
# radiative conductivity 16σT³/(3κ) (0.32 W/m·K at 600 K) is of the order of
# the conductive k (0.38 W/m·K). Solved without the radiative term in the
# manufactured source, the temperature misses T* by 5 K in L2, 10⁵ times
# the 512-cell error. The mesh is compressed to 70 % of its thickness by
# V-MMS-03's prescribed uniform stretch
#
#   w*(z, t) = (ṡ/s)·z,   s(t) = 1 − 0.15·(1 − cos(2πt/20)),
#
# moving over a solid at rest, so every source is the Eulerian one,
# evaluated at the cell centers (exact/mms.jl header, "ALE"). The solid is
# uniform in z and decays in time, so κ changes with time; a z-varying solid
# would cross the moving faces, which the case avoids as V-MMS-03 does.
#
# The exposed face is a gray interface facing surroundings at T∞(t) through a
# `RadiativeBC`, beside the manufactured conductive flux (`HeatFluxBC`); the
# back carries only the manufactured conductive flux, so it is a mirror. The
# manufactured source removes the continuous two-stream source of T*
# (`mms_radiative_source`: the same face rows, the transfer solved exactly by
# quadrature over the moving column), while the solver solves its radiation
# rows for q± itself. The temperature error therefore carries the rows'
# discretization error: the emission is constant within each cell, which is
# second order for optically thin cells (κΔ = 0.033 at the CI mesh).
#
# Expected order: second. The L2 rates are 1.86, 1.92, 1.96, 1.98 from 64 to
# 1024 cells, the same in both energy forms; tightening the time tolerances
# by 100 moves the 512-cell error by 5×10⁻⁵ of itself, so the error is the
# spatial one.

function _verification_case_vmms05()
    L = 0.01
    t_end = 10.0
    τ = 8.0
    ω = 2π / 40.0
    ε = 0.8
    α = 0.33
    σ = STEFAN_BOLTZMANN

    s(t) = 1.0 - 0.15 * (1.0 - cos(2π * t / 20.0))
    sdot(t) = -0.15 * (2π / 20.0) * sin(2π * t / 20.0)
    w_star(z, t) = (sdot(t) / s(t)) * z
    z_top(t) = L * s(t)

    T_f(z, t) = 600.0 + 150.0 * sin(π * z / L) * cos(ω * t) + 40.0 * z / L
    ξs_f(z, t) = 600.0 + 40.0 * exp(-t / τ)
    T_inf(t) = 700.0 + 80.0 * sin(2π * t / 10.0)
    fields = MMSFields(T_f, (ξs_f,))

    # The continuous radiative source of the manufactured fields at the
    # position z: depth z_top − z below the moving exposed face.
    function radiation(z, t)
        top = z_top(t)
        return mms_radiative_source(
            x -> σ * T_f(top - x, t)^4,
            top - z,
            top,
            α * ξs_f(z, t);
            ε,
            Q = σ * T_inf(t)^4,
        )
    end

    make_material() = Material(
        name = :VMMS05,
        components = (
            SolidComponent(
                :s,
                ρ = 800.0,
                c = T -> 950.0 + 0.4 * T,
                k = T -> 0.20 + 3.0e-4 * T,
                ε = ε,
                α = α,
            ),
        ),
        reactions = (),
    )

    function build(n_cells)
        mat = make_material()
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(t -> -mms_heat_flux(mat, fields, z_top(t), t)) +
                      RadiativeBC(T_inf = T_inf),
                bottom = HeatFluxBC(t -> mms_heat_flux(mat, fields, 0.0, t)),
            ),
        )
        return PyrolysisProblem(
            mesh = uniform_mesh(L, n_cells),
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_f(z, 0.0),
            ξ_initial = [z -> ξs_f(z, 0.0)],
            tspan = (0.0, t_end),
            output_times = collect(range(0.0, t_end, length = 41)),
        )
    end

    solve_fn(problem; kwargs...) = solve_mms(
        problem,
        fields;
        use_ale = true,
        correction = PrescribedCorrection(w_star),
        energy_state = VERIFICATION_ENERGY_STATE,
        radiation_model = RADIATIVE_TRANSFER,
        radiation,
        kwargs...,
    )

    function qois(sol)
        err_T = mms_Linf_T(sol, fields)
        err_ξ = mms_Linf_xi(sol, fields, 1) / 600.0
        return [
            QoIResult(
                name = "L∞ T error on moving mesh, t=$(t_end) s",
                value = err_T,
                error = err_T,
                tol = 0.011,  # 3× observed: 3.56e-3 K (n=64), both energy forms
                comment = "manufactured T* with radiative transfer on a 30%-compressed mesh",
            ),
            QoIResult(
                name = "relative L∞ ξ_solid (uniform) error, t=$(t_end) s",
                value = err_ξ,
                error = err_ξ,
                tol = 1e-10,  # observed 2.5e-13 to 3.4e-13 (n=64): exact to the integrator
                comment = "uniform solid: its absorption coefficient follows the manufactured one",
            ),
        ]
    end

    l2_T_final(sol) = mms_L2_T(sol, fields)

    return VerificationCase(
        id = "V-MMS-05",
        title = "manufactured solution with radiative transfer on a moving mesh",
        tags = [:mms, :ale, :radiation, :order_verification],
        references = [
            "docs/guide/technical_reference/18_verification-validation.md §18.2 (MMS driver, ALE cases)",
            "Modest & Mazumder, Radiative Heat Transfer, 4th ed., §14.3 (two-flux method)",
            "exact/mms.jl (`mms_two_stream`, self-tested against the discrete rows)",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        solve_fn = solve_fn,
        qois = qois,
        ci = CIConfig(
            n_cells = 64,
            budget_s = 5.0,
            convergence = ConvergenceSpec(
                # observed rates 1.917, 1.957 in both energy forms
                n_cells = [128, 256, 512],
                error_fn = l2_T_final,
                expected_order = 2.0,
                min_order = 1.8,
                min_mean_order = 1.85,
                solve_kwargs = (abstol = 1e-11, reltol = 1e-9),
            ),
        ),
        report = (
            n_cells_ladder = [16, 32, 64, 128, 256, 512],
            qoi_curves = (
                :T_profiles,
                :radiative_source,
                :error_profiles,
                :thickness_history,
            ),
            curves = sol -> begin
                m = length(sol.t)
                jmid = argmin(abs.(sol.t .- t_end / 2))
                [
                    ReportCurve(
                        name = :T_profiles,
                        title = "V-MMS-05: manufactured T* with radiative transfer on the compressing mesh",
                        xlabel = "z (m)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "t = $(round(sol.t[j], digits = 1)) s",
                                x = sol.z[:, j],
                                y = sol.T[:, j],
                                exact = let tj = sol.t[j]
                                    z -> T_f(z, tj)
                                end,
                            ) for j in (1, jmid, m)
                        ],
                    ),
                    ReportCurve(
                        name = :radiative_source,
                        title = "V-MMS-05: continuous radiative source of T*",
                        xlabel = "z (m)",
                        ylabel = "S_rad (W/m³)",
                        series = [
                            CurveSeries(
                                label = "t = $(round(sol.t[j], digits = 1)) s",
                                x = sol.z[:, j],
                                y = [radiation(z, sol.t[j]) for z in sol.z[:, j]],
                            ) for j in (1, jmid, m)
                        ],
                    ),
                    ReportCurve(
                        name = :error_profiles,
                        title = "V-MMS-05: pointwise temperature error, t = $(t_end) s",
                        xlabel = "z (m)",
                        ylabel = "T − T* (K)",
                        series = [
                            CurveSeries(
                                label = "T − T*",
                                x = sol.z[:, m],
                                y = [
                                    sol.T[i, m] - T_f(sol.z[i, m], t_end) for
                                    i in axes(sol.T, 1)
                                ],
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :thickness_history,
                        title = "V-MMS-05: top-node trajectory under prescribed uniform stretch",
                        xlabel = "t (s)",
                        ylabel = "top node z (m)",
                        series = [
                            CurveSeries(
                                label = "top node (node ODE)",
                                x = sol.t,
                                y = [mms_top_node(sol, j) for j in eachindex(sol.t)],
                                exact = z_top,
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol ->
                (L2_T = mms_L2_T(sol, fields), Linf_T = mms_Linf_T(sol, fields)),
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vmms05)
