# =============================================================================
# V-MMS-03 — ALE manufactured solution with prescribed smooth mesh motion
# =============================================================================
#
# The moving-mesh order-verification companion to V-MMS-01 (and the "hardest
# and most valuable single test" of the suite plan's Tier 3): non-trivial
# manufactured fields solved on a mesh compressed to 70 % of its thickness by
# the prescribed uniform-stretch velocity
#
#   w*(z, t) = (ṡ/s)·z,   s(t) = 1 − 0.15·(1 − cos(2πt/20)),  z_i(t) = z_i(0)·s(t)
#
# through the `Workspace.w_override` hook. On a grid point moving with w, the
# grid time derivative is ∂tT* + w·∂zT*. For T the discrete mesh-relative
# advection kernel supplies the w·∂zT* part (w_physical = 0 without
# reactions), so the Eulerian manufactured T source — evaluated at the LIVE
# cell centers read from the z-block — makes T* exact on the moving mesh. The
# species ride the applied mesh (no advection operator), so their sources are
# assembled in the material frame with w = the override (see exact/mms.jl
# header, "ALE"). The T error norm therefore measures exactly the moving-mesh
# transport truncation: ALE advection (van-Leer-limited upwind
# reconstruction) + lab-frame operators on a nonuniform-in-time mesh; the gas
# error measures the diffusive flux and geometry on the moving mesh.
#
# No reactions: mesh motion is prescribed, not volume-change-driven, so the
# case isolates GCL/advection from kinetics (V-MMS-01 covers those). The
# solid field is spatially UNIFORM but time-varying (the material-frame
# source would admit a z-varying one; uniform keeps the T pins focused on
# the advection kernel); T and the gas field carry the full space–time
# structure.
#
# Expected order: the limiter's slope reconstruction is formally 2nd order on
# smooth data but clips at extrema (the sin peak crosses the domain), so the
# observed order sits between 1.5 and 2 — recorded at pin time; a drop below
# min_order flags a moving-mesh transport regression.
#
# Jacobian: dense FD (`jacobian = :fd`) — the structured geometry Jacobian
# assumes reaction-driven mesh velocity (see `Workspace.w_override` note).

function _verification_case_vmms03()
    L = 0.01
    t_end = 10.0
    τ = 8.0
    ω = 2π / 40.0

    s(t) = 1.0 - 0.15 * (1.0 - cos(2π * t / 20.0))
    sdot(t) = -0.15 * (2π / 20.0) * sin(2π * t / 20.0)
    w_star(z, t) = (sdot(t) / s(t)) * z

    # Lab-frame fields over the maximal domain [0, L] (compression-only
    # motion keeps every node inside it). Solid: uniform in z, decaying in t.
    T_f(z, t) = 480.0 + 120.0 * sin(π * z / L) * cos(ω * t) + 40.0 * z / L
    ξs_f(z, t) = 600.0 + 40.0 * exp(-t / τ)
    ξg_f(z, t) = 1.5 + 0.6 * sin(π * z / L) * exp(-t / τ)
    fields = MMSFields(T_f, (ξs_f, ξg_f))

    function make_material()
        solid = SolidComponent(
            :s,
            ρ = 800.0,
            c = T -> 950.0 + 0.4 * T,
            k = T -> 0.20 + 3.0e-4 * T,
            λ = 2.0e-5,
        )
        gas = GaseousComponent(:g, M = 0.028, c = 1000.0, k = 0.03, λ = 2.0e-5)
        return Material(name = :VMMS03, components = (solid, gas), reactions = ())
    end

    function build(n_cells)
        mat = make_material()
        mesh = create_uniform_mesh(
            L,
            n_cells,
            Val(2);
            T_initial = 500.0,
            ξ_initial = (600.0, 1.5),
        )
        # Manufactured BC fluxes must be evaluated at the MOVING boundary:
        # the top face sits at z_top(t) = L·s(t) (bottom pinned at 0), so the
        # fixed-L helpers do not apply here.
        z_top(t) = L * s(t)
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(t -> -mms_heat_flux(mat, fields, z_top(t), t)),
                bottom = HeatFluxBC(t -> mms_heat_flux(mat, fields, 0.0, t)),
            ),
            mass = (
                top = MassFluxBC(t -> -mms_gas_flux(mat, fields, 2, z_top(t), t)),
                bottom = MassFluxBC(t -> mms_gas_flux(mat, fields, 2, 0.0, t)),
            ),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_f(z, 0.0),
            ξ_initial = [z -> ξs_f(z, 0.0), z -> ξg_f(z, 0.0)],
            tspan = (0.0, t_end),
            output_times = collect(range(0.0, t_end, length = 41)),
        )
    end

    solve_fn(problem; kwargs...) = solve_mms(
        problem,
        fields;
        use_ale = true,
        w_override = w_star,
        # Dense FD accepted for offline report runs too (deferred decision,
        # 2026-07-06 audit): do NOT hand a physics+geometry union prototype to
        # colored FD without first proving the pattern covers the override's
        # ∂(dz/dt)/∂z coupling — silent-drop hazard.
        jacobian = :fd,
        kwargs...,
    )

    function qois(sol)
        out = QoIResult[]
        push!(
            out,
            QoIResult(
                name = "L∞ T error on moving mesh, t=$(t_end) s",
                value = mms_Linf_T(sol, fields),
                error = mms_Linf_T(sol, fields),
                tol = 0.25,  # 3× observed: 0.069 K (n=32, 2026-07-06)
                comment = "manufactured T* on 30%-compressed mesh",
            ),
        )
        for (name, idx, scale, tol, comment) in (
            (
                "ξ_solid (mesh-bound)",
                1,
                600.0,
                1e-10,  # observed 7.3e-12 (n=32, 2026-07-06): exact to integrator tol
                "uniform solid: verifies the mesh-bound (no-advection) species path + source injection",
            ),
            (
                "ξ_gas",
                2,
                1.5,
                6e-3,  # 3× observed: 1.8e-3 (n=32), second order
                "gas rides the applied mesh: diffusive flux + geometry on the moving mesh",
            ),
        )
            err = mms_Linf_xi(sol, fields, idx) / scale
            push!(
                out,
                QoIResult(
                    name = "relative L∞ $(name) error, t=$(t_end) s",
                    value = err,
                    error = err,
                    tol = tol,
                    comment = comment,
                ),
            )
        end
        return out
    end

    l2_T_final(sol) = mms_L2_T(sol, fields)

    return VerificationCase(
        id = "V-MMS-03",
        title = "ALE manufactured solution, prescribed smooth mesh motion",
        tags = [:mms, :ale, :gcl, :order_verification],
        references = [
            "docs/Verification_Suite_Plan_2026-07-06.md Tier 3 (ALE-MMS deliverable)",
            "Thomas & Lombard, AIAA J. 17 (1979) — geometric conservation law",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        solve_fn = solve_fn,
        qois = qois,
        ci = CIConfig(
            n_cells = 32,
            budget_s = 60.0,
            convergence = ConvergenceSpec(
                # observed rates 2.24, 2.18 (2026-07-06) — the anticipated
                # limiter clipping at the sin extremum did not depress the
                # global L2 order
                n_cells = [16, 32, 64],
                error_fn = l2_T_final,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (abstol = 1e-11, reltol = 1e-9),
            ),
        ),
        report = (
            n_cells_ladder = [16, 32, 64, 128],
            qoi_curves = (:T_profiles, :xi_profiles, :error_profiles, :thickness_history),
            curves = sol -> begin
                m = length(sol.t)
                jmid = argmin(abs.(sol.t .- t_end / 2))
                [
                    ReportCurve(
                        name = :T_profiles,
                        title = "V-MMS-03: manufactured T* on the compressing mesh",
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
                            ) for j in (jmid, m)
                        ],
                    ),
                    ReportCurve(
                        name = :xi_profiles,
                        title = "V-MMS-03: mesh-bound solid recovery (z-uniform, " *
                                "time-decaying), t = $(t_end) s",
                        xlabel = "z (m)",
                        ylabel = "ξ_s (kg/m³)",
                        series = [
                            CurveSeries(
                                label = "ξ_s",
                                x = sol.z[:, m],
                                y = sol.ξ[1, :, m],
                                exact = z -> ξs_f(z, t_end),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :xi_gas_profile,
                        title = "V-MMS-03: gas recovery on the moving mesh, " *
                                "t = $(t_end) s",
                        xlabel = "z (m)",
                        ylabel = "ξ_g (kg/m³)",
                        series = [
                            CurveSeries(
                                label = "ξ_g",
                                x = sol.z[:, m],
                                y = sol.ξ[2, :, m],
                                exact = z -> ξg_f(z, t_end),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :error_profiles,
                        title = "V-MMS-03: pointwise field errors, t = $(t_end) s",
                        xlabel = "z (m)",
                        ylabel = "error (K, kg/m³)",
                        series = [
                            CurveSeries(
                                label = "T − T*",
                                x = sol.z[:, m],
                                y = [
                                    sol.T[i, m] - T_f(sol.z[i, m], t_end) for
                                    i in axes(sol.T, 1)
                                ],
                            ),
                            CurveSeries(
                                label = "ξ_s − ξ_s*",
                                x = sol.z[:, m],
                                y = [
                                    sol.ξ[1, i, m] - ξs_f(sol.z[i, m], t_end) for
                                    i in axes(sol.T, 1)
                                ],
                            ),
                            CurveSeries(
                                label = "ξ_g − ξ_g*",
                                x = sol.z[:, m],
                                y = [
                                    sol.ξ[2, i, m] - ξg_f(sol.z[i, m], t_end) for
                                    i in axes(sol.T, 1)
                                ],
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :thickness_history,
                        title = "V-MMS-03: top-node trajectory under prescribed " *
                                "uniform stretch",
                        xlabel = "t (s)",
                        ylabel = "top node z (m)",
                        series = [
                            CurveSeries(
                                label = "top node (node ODE)",
                                x = sol.t,
                                y = [mms_top_node(sol, j) for j in eachindex(sol.t)],
                                exact = t -> L * s(t),
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

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vmms03)
