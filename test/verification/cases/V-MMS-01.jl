# =============================================================================
# V-MMS-01 — Manufactured solution, fully coupled Eulerian system
# =============================================================================
#
# Method of Manufactured Solutions over the coupled energy + kinetics + gas
# transport operators — the only formal order-of-accuracy verification of the
# *coupled* system (every other case has zero residual source; tech ref §18.2
# documents this gap; FDS uses MMS for the same reason in the fluid phase).
#
# Active couplings, all with O(1) contributions over the space–time domain:
# - T-dependent k(T), c(T) on both condensed components (mixture rules live,
#   composition varying in z and t);
# - two-step kinetics A → 0.65 B + 0.35 gas (endothermic) and B → gas
#   (exothermic): intermediate B produced and consumed, reaction heat feeding
#   back into the energy equation, gas produced into the transport field;
# - partial-pressure gas transport J = −(Λ/T)∂z(ξT/φ), with condensed
#   matrix coefficients and the skeletal pore volume evaluated independently;
# - gas-advection energy source S_adv = −cp_g·J·∂zT;
# - manufactured Neumann BCs (HeatFluxBC/MassFluxBC callables of t) — the
#   Dirichlet face reconstruction is first-order and would cap the observed
#   order (see exact/mms.jl header, design decision 3).
#
# The manufactured sources are assembled by nested ForwardDiff through the
# package's own property/kinetics functions (exact/mms.jl) and injected at the
# state-derivative level as functions of (z, t), so the colored Jacobian stays
# exact.
#
# The depletion limiter is left at the material default: the source generator
# evaluates the identical rate law, and every reactant concentration stays
# > 20× the 1 kg/m³ threshold, where the tanh factor is exactly 1.0 in
# Float64. Temperature gates are at defaults (fully open).
#
# Fields (z̃ = z/L, smooth, strictly positive concentrations):
#   T*  = 500 + 150·sin(πz̃)·cos(ωt) + 50·z̃            [K]     (≈ 350–700 K)
#   ξ*_A = 350 + 80·z̃²·e^(−t/τ)                        [kg/m³]
#   ξ*_B = 220 + 60·(1−z̃)²·(1 − 0.5·e^(−t/τ))          [kg/m³]
#   ξ*_g = 1.2 + 0.5·sin(πz̃)·(1 + 0.3·sin(ωt))         [kg/m³]

function _verification_case_vmms01()
    L = 0.01
    τ = 10.0
    ω = 2π / 40.0
    t_end = 20.0
    times = [10.0, t_end]

    T_f(z, t) = 500.0 + 150.0 * sin(π * z / L) * cos(ω * t) + 50.0 * z / L
    ξA_f(z, t) = 350.0 + 80.0 * (z / L)^2 * exp(-t / τ)
    ξB_f(z, t) = 220.0 + 60.0 * (1.0 - z / L)^2 * (1.0 - 0.5 * exp(-t / τ))
    ξg_f(z, t) = 1.2 + 0.5 * sin(π * z / L) * (1.0 + 0.3 * sin(ω * t))
    fields = MMSFields(T_f, (ξA_f, ξB_f, ξg_f))

    function make_material()
        compA = SolidComponent(
            :A,
            ρ = 800.0,
            c = T -> 900.0 + 0.6 * T,
            k = T -> 0.10 + 6.0e-4 * T,   # ≈ 0.31 → 0.52 over 350–700 K
            λ = 1.5e-5,
        )
        compB = SolidComponent(
            :B,
            ρ = 550.0,
            c = T -> 1100.0 + 0.4 * T,
            k = T -> 0.32 + 3.0e-4 * T,
            λ = 2.5e-5,
        )
        gas = GaseousComponent(:g, M = 0.028, c = T -> 950.0 + 0.3 * T, k = 0.03)
        r1 = Reaction(
            :r1,
            1 => (2, 3);
            A = 200.0,
            E = 45e3,
            h = 1.2e5,
            yields = (0.65, 0.35),
        )
        r2 = Reaction(:r2, 2 => 3; A = 50.0, E = 40e3, h = -0.8e5)
        return Material(
            name = :VMMS01,
            components = (compA, compB, gas),
            reactions = (r1, r2),
        )
    end

    function build(n_cells)
        mat = make_material()
        mesh = create_uniform_mesh(L, n_cells)
        thermal = mms_thermal_bcs(mat, fields, L)
        mass = mms_mass_bcs(mat, fields, L)
        bc_set = BoundaryConditionSet(
            thermal = (top = thermal.top, bottom = thermal.bottom),
            mass = (top = mass.top, bottom = mass.bottom),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_f(z, 0.0),
            ξ_initial = [z -> ξA_f(z, 0.0), z -> ξB_f(z, 0.0), z -> ξg_f(z, 0.0)],
            tspan = (0.0, t_end),
            output_times = vcat(0.0, times),
        )
    end

    solve_fn(problem; kwargs...) =
        solve_mms(problem, fields; energy_state = VERIFICATION_ENERGY_STATE, kwargs...)

    function qois(sol)
        out = QoIResult[]
        for (j, t) in enumerate(times)
            m = j + 1  # output_times = [0; times]
            push!(
                out,
                QoIResult(
                    name = "L∞ T error, t=$(t) s",
                    value = mms_Linf_T(sol, fields; m = m),
                    error = mms_Linf_T(sol, fields; m = m),
                    tol = 0.05,  # K; prescribed physical accuracy bound
                    comment = "manufactured T* recovered on the coupled system; " *
                              "refined to 160 cells to resolve the pore-pressure gradients",
                ),
            )
        end
        for (name, idx, scale, tol, obs) in (
            ("ξ_A", 1, 350.0, 2e-4, "5.3e-6"),
            ("ξ_B", 2, 220.0, 1e-4, "1.7e-6"),
            ("ξ_gas", 3, 1.2, 3e-2, "1.34e-2"),
        )
            err = mms_Linf_xi(sol, fields, idx) / scale
            push!(
                out,
                QoIResult(
                    name = "relative L∞ $(name) error, t=$(t_end) s",
                    value = err,
                    error = err,
                    tol = tol,  # fixed concentration accuracy bounds
                    comment = "scale = manufactured baseline; observed $(obs); " *
                              "L∞ ξ order ≈ 2 on [20,40,80,160]",
                ),
            )
        end
        return out
    end

    l2_T_final(sol) = mms_L2_T(sol, fields)

    return VerificationCase(
        id = "V-MMS-01",
        title = "manufactured solution, coupled kinetics + T-props + gas transport",
        tags = [:mms, :coupled, :kinetics, :gas_transport, :order_verification],
        references = [
            "Roache, Verification and Validation in Computational Science (MMS methodology)",
            "docs/guide/technical_reference/18_verification-validation.md §18.2 (the documented gap)",
            "FDS Verification Guide (MMS usage precedent, fluid phase)",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        solve_fn = solve_fn,
        qois = qois,
        ci = CIConfig(
            n_cells = 160,
            budget_s = 30.0,
            convergence = ConvergenceSpec(
                # The independent pressure-potential source gives second-order
                # temperature and species convergence on this ladder.
                n_cells = [20, 40, 80],
                error_fn = l2_T_final,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (abstol = 1e-11, reltol = 1e-9),
            ),
        ),
        report = (
            n_cells_ladder = [20, 40, 80, 160, 320],
            profile_times = times,
            qoi_curves = (:T_profiles, :xi_profiles, :error_profiles),
            temporal_sweep = (n_cells = 80, dts = [0.5, 0.25, 0.125, 0.0625]),
            solve_kwargs = (abstol = 1e-11, reltol = 1e-9),
            curves = sol -> begin
                # MMSSolution: only t/T/ξ/z (no extras, no interpolation)
                species = (("ξ_A", 1, ξA_f), ("ξ_B", 2, ξB_f), ("ξ_g", 3, ξg_f))
                [
                    ReportCurve(
                        name = :T_profiles,
                        title = "V-MMS-01: manufactured T* recovery (coupled system)",
                        xlabel = "z (m)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "t = $(t) s",
                                x = sol.z[:, j+1],
                                y = sol.T[:, j+1],
                                exact = z -> T_f(z, t),
                            ) for (j, t) in enumerate(times)
                        ],
                    ),
                    ReportCurve(
                        name = :xi_profiles,
                        title = "V-MMS-01: manufactured condensed-species recovery, t = $(t_end) s",
                        xlabel = "z (m)",
                        ylabel = "ξ (kg/m³)",
                        series = [
                            CurveSeries(
                                label = name,
                                x = sol.z[:, end],
                                y = sol.ξ[idx, :, end],
                                exact = z -> f(z, t_end),
                            ) for (name, idx, f) in species[1:2]
                        ],
                    ),
                    ReportCurve(
                        name = :xi_gas_profile,
                        title = "V-MMS-01: manufactured gas-species recovery, t = $(t_end) s",
                        xlabel = "z (m)",
                        ylabel = "ξ_g (kg/m³)",
                        series = [
                            CurveSeries(
                                label = "ξ_g",
                                x = sol.z[:, end],
                                y = sol.ξ[3, :, end],
                                exact = z -> ξg_f(z, t_end),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :error_profiles,
                        title = "V-MMS-01: pointwise field errors, t = $(t_end) s",
                        xlabel = "z (m)",
                        ylabel = "error (K, kg/m³)",
                        series = vcat(
                            [
                                CurveSeries(
                                    label = "T − T*",
                                    x = sol.z[:, end],
                                    y = [
                                        sol.T[i, end] - T_f(sol.z[i, end], t_end) for
                                        i in axes(sol.T, 1)
                                    ],
                                ),
                            ],
                            [
                                CurveSeries(
                                    label = "$(name) − $(name)*",
                                    x = sol.z[:, end],
                                    y = [
                                        sol.ξ[idx, i, end] - f(sol.z[i, end], t_end) for
                                        i in axes(sol.T, 1)
                                    ],
                                ) for (name, idx, f) in species
                            ],
                        ),
                    ),
                ]
            end,
            error_norms = sol ->
                (L2_T = mms_L2_T(sol, fields), Linf_T = mms_Linf_T(sol, fields)),
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vmms01)
