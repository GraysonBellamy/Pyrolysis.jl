# =============================================================================
# V-MMS-04 — induced-motion ALE manufactured solution (reaction-driven mesh)
# =============================================================================
#
# V-MMS-02/03 verify GCL + moving-mesh transport under PRESCRIBED motion
# (`w_override`), which moves the mesh over a solid at rest and deliberately
# bypasses the balance rows' mesh velocity (the integral of θ) and the
# cell-volume change with θ ≠ 0. Here the kinetics are manufactured so the
# rates are analytic on the manufactured fields and the motion is genuinely
# REACTION-DRIVEN — no `w_override` — closing the kinetics→θ→w→geometry chain.
#
# ## Manufactured design (see exact/mms.jl, "Induced mesh motion")
#
# Constant intrinsic densities, γ_solid = 1 / γ_gas = 0 (canonical), single
# first-order reaction s → g with the depletion limiter numerically inactive
# over the manufactured range (tanh(ξ_s/1.0) ≡ 1.0 in Float64 for
# ξ_s ≳ 20 kg/m³; here ξ_s ≥ 500). T* is spatially UNIFORM (time-varying):
# V-MMS-01/03 already order-verify the spatial T operators — this case's
# target is the kinetics→geometry chain. The solid and the gas are both
# z-VARYING and both ride the reacting mesh: every species source is
# assembled in the material frame, S = ∂tξ* + w*·∂zξ* + [∂zJ*]_gas − ω̇* +
# ξ*·θ*, so the moving-point chain rule and the growth of the cell volume are
# built in (`mms_species_source` with `w` from the `MMSMotion` companion).
#
#   θ*(z,t) = −rate*/ρ_s = −k(T*(t))·ξ_s*(z,t)/ρ_s   (volume_change_rate)
#   w*(z,t) = ∫₀ᶻ θ* dζ = −(k(t)/ρ_s)·[a(t)·z + b(t)·(L/π)(1 − cos(πz/L))]
#
# The case builds w* with `mms_mesh_velocity` (32-node Gauss–Legendre on the
# package kinetics — machine precision, requirement 4); the closed form above
# is the independent hand-check in exact/selftests.jl. θ* < 0 throughout ⇒
# monotone ~8 % compression over 20 s, every node inside the maximal domain
# [0, L].
#
# ## Moving-boundary BC positions
#
# Induced motion means the exact top-face trajectory solves the 1D node ODE
# dz*/dt = w*(z*, t), z*(0) = L — not analytic here, so `build` precomputes it
# with Vern9 at 1e-14 tolerances (requirement-4 compliant) and the BC closures
# read the dense interpolant. The bottom node is pinned (w*(0,t) = 0).

function _verification_case_vmms04()
    L = 0.01
    t_end = 20.0

    ρ_s = 800.0
    A_rxn = 5.0e4
    E_rxn = 80.0e3

    # Manufactured fields on the maximal domain [0, L].
    T_f(z, t) = 600.0 + 15.0 * sin(2π * t / 40.0)                 # z-uniform
    a_f(t) = 600.0 + 20.0 * cos(2π * t / 25.0)
    b_f(t) = 60.0 * exp(-t / 30.0)
    ξs_f(z, t) = a_f(t) + b_f(t) * sin(π * z / L)                 # z-varying solid
    ξg_f(z, t) = 2.0 + 0.8 * sin(π * z / L) * exp(-t / 15.0)
    fields = MMSFields(T_f, (ξs_f, ξg_f))

    function make_material()
        solid = SolidComponent(:s, ρ = ρ_s, c = 1200.0, k = 0.2, λ = 2.0e-5, γ = 1.0)
        gas = GaseousComponent(:g, M = 0.028, c = 1000.0, k = 0.03)  # γ = 0
        rxn = Reaction(:pyro, 1 => 2, A = A_rxn, E = E_rxn, h = 5.0e4)
        return Material(name = :VMMS04, components = (solid, gas), reactions = (rxn,))
    end

    mat = make_material()

    # Reaction-driven mesh velocity from the package kinetics (θ source of
    # truth = Physics.volume_change_rate; Gauss–Legendre to machine precision).
    w_star(z, t) = mms_mesh_velocity(mat, fields, z, t)
    motion = MMSMotion(w_star)

    # Exact top-face trajectory: dz*/dt = w*(z*, t), z*(0) = L, integrated far
    # below the truncation error under study (requirement 4).
    ztop_sol = OrdinaryDiffEq.solve(
        OrdinaryDiffEq.ODEProblem((z, p, t) -> w_star(z, t), L, (0.0, t_end)),
        OrdinaryDiffEq.Vern9();
        abstol = 1e-14,
        reltol = 1e-14,
    )
    z_top(t) = ztop_sol(t)

    function build(n_cells)
        mesh = create_uniform_mesh(L, n_cells)
        # Manufactured BC fluxes at the MOVING top face z*_top(t); bottom
        # pinned at 0. The thermal fluxes are identically zero for z-uniform
        # T*, but the closures stay field-derived so the case remains correct
        # if the fields are ever changed.
        bc_set = BoundaryConditionSet(
            thermal = (
                top = HeatFluxBC(t -> -mms_heat_flux(mat, fields, z_top(t), t)),
                bottom = HeatFluxBC(t -> mms_heat_flux(mat, fields, 0.0, t)),
            ),
            mass = (
                top = MassFluxBC(
                    fluxes = (g = t -> mms_gas_flux(mat, fields, 2, z_top(t), t),),
                    inflow_temperature = (g = t -> fields.T(z_top(t), t),),
                ),
                bottom = MassFluxBC(
                    fluxes = (g = t -> -mms_gas_flux(mat, fields, 2, 0.0, t),),
                    inflow_temperature = (g = t -> fields.T(0.0, t),),
                ),
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
        motion = motion,
        energy_state = VERIFICATION_ENERGY_STATE,
        kwargs...,
    )

    # Numeric top-node position from the saved cell centers (exact/mms.jl).
    numeric_top_node(sol) = mms_top_node(sol, length(sol.t))

    recession(t) = L - z_top(t)   # exact total surface recession [m]

    function qois(sol)
        out = QoIResult[]
        push!(
            out,
            QoIResult(
                name = "L∞ T error, t=$(t_end) s",
                value = mms_Linf_T(sol, fields),
                error = mms_Linf_T(sol, fields),
                tol = 6.0e-5,  # 3× observed: 2.0e-5 K (n=32). z-uniform T*, so the
                # error is the gas error fed back through ρcp
                # (ΔT ≈ ∫ Q_rxn/ρcp · δρcp/ρcp dt), not the spatial T operators
                comment = "z-uniform T* — kinetics/source channel through ρcp(ξ_g)",
            ),
        )
        for (name, idx, scale, tol, comment) in (
            (
                "ξ_solid (z-varying, material frame)",
                1,
                600.0,
                2.0e-6,  # 3× observed: 6.6e-7 (n=32)
                "solid rides the reaction-driven mesh: balance-row w + dilation of m/V + material-frame source",
            ),
            (
                "ξ_gas",
                2,
                2.0,
                7.0e-3,  # 3× observed: 2.3e-3 (n=32), second order (0.0184, 0.0046, 0.00115, 0.00029 at n = 16..128)
                "gas rides the reacting, compressing mesh: dilation of m/V + diffusion + material-frame source",
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

        # Surface-position error relative to the exact total recession — the
        # direct kinetics→θ→w→geometry chain measurement.
        z_exact = z_top(t_end)
        z_num = numeric_top_node(sol)
        rec = recession(t_end)
        pos_err = abs(z_num - z_exact) / rec
        push!(
            out,
            QoIResult(
                name = "relative surface-position error, t=$(t_end) s",
                value = pos_err,
                exact = z_exact,
                error = pos_err,
                tol = 5.0e-5,  # 3× observed: 1.6e-5 of the exact recession (n=32)
                comment = "|z_top − z*_top| / (L − z*_top); recession ≈ $(round(rec*1e3, sigdigits=3)) mm",
            ),
        )

        return out
    end

    # Convergence on the solid channel — the field that exercises the
    # induced-motion machinery (kinetics→geometry chain + material-frame
    # source).
    function l2_xi_s_final(sol)
        m = length(sol.t)
        t = sol.t[m]
        n = size(sol.T, 1)
        return sqrt(sum(abs2(sol.ξ[1, i, m] - ξs_f(sol.z[i, m], t)) for i = 1:n) / n)
    end

    return VerificationCase(
        id = "V-MMS-04",
        title = "Induced-motion ALE manufactured solution (reaction-driven mesh)",
        tags = [:mms, :ale, :kinetics, :order_verification],
        references = [
            "Roache, Verification and Validation in Computational Science (MMS methodology)",
            "docs/guide/technical_reference/18_verification-validation.md §18.2 (MMS driver, ALE cases)",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        solve_fn = solve_fn,
        qois = qois,
        ci = CIConfig(
            n_cells = 64,
            budget_s = 60.0,
            convergence = ConvergenceSpec(
                # observed rates 2.003, 2.007 (n = 16→32→64; L2 ξ_s errors
                # 7.6e-4 → 1.9e-4 → 4.7e-5)
                n_cells = [16, 32, 64],
                error_fn = l2_xi_s_final,
                expected_order = 2.0,
                min_order = 1.7,
                min_mean_order = 1.9,
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
                        title = "V-MMS-04: z-uniform manufactured T* on the " *
                                "reaction-driven mesh",
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
                        title = "V-MMS-04: z-varying solid recovery (material-frame " *
                                "source), t = $(t_end) s",
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
                        title = "V-MMS-04: gas recovery on the reacting, compressing " *
                                "mesh, t = $(t_end) s",
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
                        title = "V-MMS-04: pointwise field errors, t = $(t_end) s",
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
                        title = "V-MMS-04: reaction-driven top-face recession vs the " *
                                "exact node ODE",
                        xlabel = "t (s)",
                        ylabel = "top node z (m)",
                        series = [
                            CurveSeries(
                                label = "top node (induced motion)",
                                x = sol.t,
                                y = [mms_top_node(sol, j) for j in eachindex(sol.t)],
                                exact = t -> z_top(t),
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol ->
                (L2_xi_s = l2_xi_s_final(sol), Linf_xi_g = mms_Linf_xi(sol, fields, 2)),
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vmms04)
