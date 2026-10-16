# =============================================================================
# V-CPL-04 — Frank-Kamenetskii critical slab (reaction–conduction bifurcation)
# =============================================================================
#
# Zero-order exothermic Arrhenius source in a slab with Dirichlet walls
# (Frank-Kamenetskii). The only case where reaction heat release feeds back
# on the temperature field with a sharply known stability boundary
# δ_c = 0.878 (slab), δ = (E/RT₀²)·S₀L²/k, L = half-width,
# S₀ = |h|·A·e^(−E/RT₀).
#
# Realization: the volumetric heat source is a SELF-REGENERATING reaction
# A → A (order n = 0, h < 0) — product index = reactant index nets ξ̇ ≡ 0
# exactly (verified: drift 0.0) while depositing |h|·A·e^(−E/RT) per unit
# volume, i.e. a pure zero-order Arrhenius source with no reactant-depletion
# clock. The depletion limiter's tanh(ξ/1) factor is exactly 1.0 at
# ξ = 1000 kg/m³, and the default gates (T_min = 0) are exactly open, so the
# kinetics are pure Arrhenius. A T_max gate at T₀ + 100 K caps the
# SUPERCRITICAL branch so the runaway saturates instead of chasing e^θ to
# integrator death; the subcritical branch (rise ≤ 15 K) never feels it.
#
# References (exact/frank_kamenetskii.jl, self-tested): the FK closed-form
# profile at E/RT₀ → ∞, AND a full-Arrhenius RK4-shooting steady solution at
# the case's finite E/RT₀ = 40 — the primary comparison, since the FK-limit
# truncation error (~2–5% at E/RT₀ = 40) is physics, not discretization.
#
# Assert-mode protocol:
#   main solve  δ = 0.44 ≈ 0.50·δ_c → steady; center/interior T vs shooting
#   bracket     δ = 0.70 ≈ 0.80·δ_c → must stay bounded (near-critical, stable)
#               δ = 1.32 ≈ 1.50·δ_c → must run away past the detection rise
# The δ_c bisection itself (turning-point hunt) is report-mode work.

function _vcpl04_problem(δ; n_cells, t_end, n_outputs = 11)
    R = 8.314462618
    T0 = 500.0
    E = 40.0 * R * T0          # E/RT₀ = 40
    k = 0.2
    ρξ = 1000.0
    c = 1000.0
    L_half = 0.01
    εT = R * T0^2 / E          # 12.5 K per θ-unit
    S0 = δ * k * εT / L_half^2 # W/m³ at T₀
    h_rxn = -1.0e5             # exothermic (our h > 0 = endothermic)
    A_r = S0 / abs(h_rxn) * exp(E / (R * T0))

    comp = SolidComponent(:A, ρ = ρξ, c = c, k = k)
    rxn =
        Reaction(:fk_source, 1 => 1; A = A_r, E = E, h = h_rxn, n = 0.0, T_max = T0 + 100.0)
    mat = Material(name = :VCPL04, components = (comp,), reactions = (rxn,))
    mesh =
        create_uniform_mesh(2 * L_half, n_cells, Val(1); T_initial = T0, ξ_initial = (ρξ,))
    bc_set = BoundaryConditionSet(
        thermal = (top = DirichletThermalBC(T0), bottom = DirichletThermalBC(T0)),
    )
    problem = PyrolysisProblem(
        mesh = mesh,
        material = mat,
        bc_set = bc_set,
        T_initial = z -> T0,
        ξ_initial = [z -> ρξ],
        tspan = (0.0, t_end),
        output_times = collect(range(0.0, t_end, length = n_outputs)),
    )
    return (; problem, T0, E, k, L_half, εT, S0, R)
end

function _verification_case_vcpl04()
    δ_main = 0.44
    t_main = 2500.0   # diffusion time L²/α = 500 s; steady well before t_end

    build(n_cells) = _vcpl04_problem(δ_main; n_cells = n_cells, t_end = t_main).problem

    function qois(sol)
        P = _vcpl04_problem(δ_main; n_cells = 2, t_end = t_main) # params only
        out = QoIResult[]
        z = get_cell_positions(sol)
        L_half = P.L_half

        # High-precision references
        sh = fk_steady_shoot(S0 = P.S0, E = P.E, T0 = P.T0, k = P.k, L = L_half)
        θm_fk = fk_theta_max(δ_main)

        # 1. Steady center temperature vs full-Arrhenius shooting reference.
        jc = argmin(abs.(z .- L_half))  # cell nearest the slab center
        Tc = sol.T[jc, end]
        Tc_ref = sh.T_at(abs(z[jc] - L_half))
        err_c = abs(Tc - Tc_ref) / (Tc_ref - P.T0)
        push!(
            out,
            QoIResult(
                name = "steady center temperature vs full-Arrhenius shooting",
                value = Tc,
                exact = Tc_ref,
                error = err_c,
                tol = 2.5e-3,  # 3× observed: 7.4e-4 (n=40)
                comment = "δ = $(δ_main) ≈ 0.5 δ_c, E/RT₀ = 40, relative to the rise",
            ),
        )

        # 2. Same, vs the FK closed form — looser: FK-limit truncation is
        # physics (O(RT₀/E)), and its size is itself a documented number.
        θ_c_num = (Tc - P.T0) / P.εT
        err_fk = abs(θ_c_num - θm_fk) / θm_fk
        push!(
            out,
            QoIResult(
                name = "steady center rise vs FK closed form (limit check)",
                value = θ_c_num,
                exact = θm_fk,
                error = err_fk,
                tol = 6e-3,  # 3× observed: 1.7e-3 (n=40); FK truncation itself is only ~0.17% here
                comment = "quantifies the FK-limit error at finite E/RT₀",
            ),
        )

        # 3. Steady profile L∞ vs shooting (both half-slabs — symmetry included).
        err_prof =
            maximum(
                abs(
                    sol.T[i, end] -
                    sh.T_at(abs(z[i] - L_half) > L_half ? L_half : abs(z[i] - L_half)),
                ) for i in eachindex(z)
            ) / (sh.T_center - P.T0)
        push!(
            out,
            QoIResult(
                name = "steady profile L∞ vs shooting (relative to center rise)",
                value = err_prof,
                error = err_prof,
                tol = 2.5e-3,  # 3× observed: 7.4e-4 (n=40)
                comment = "includes the symmetry check",
            ),
        )

        # 4. Steadiness guard: last two outputs agree.
        drift = maximum(abs.(sol.T[:, end] - sol.T[:, end-1]))
        push!(
            out,
            QoIResult(
                name = "steady-state reached (last-output drift)",
                value = drift,
                exact = 0.0,
                error = drift,
                tol = 1e-3,  # K over 250 s
                comment = "t_end = 5 diffusion times",
            ),
        )

        # 5. Criticality bracket. Near-critical stable: δ = 0.8 δ_c must stay
        # bounded near its steady value; supercritical: δ = 1.5 δ_c must run
        # away beyond a 50 K rise (T_max gate caps it at ~100 K).
        r_sub = _vcpl04_problem(0.70; n_cells = 40, t_end = 5000.0)
        sol_sub = solve(r_sub.problem; progress = false, abstol = 1e-8, reltol = 1e-6)
        rise_sub = maximum(sol_sub.T[:, end]) - P.T0
        θ_sub = rise_sub / P.εT
        θ_sub_ref = fk_theta_max(0.70)
        push!(
            out,
            QoIResult(
                name = "near-critical δ = 0.80 δ_c stays on the stable branch",
                value = θ_sub,
                exact = θ_sub_ref,
                error = abs(θ_sub - θ_sub_ref) / θ_sub_ref,
                tol = 3e-2,  # 3× observed: 9.9e-3 (n=40); FK reference
                comment = "bounded, near the FK lower branch (θ ≈ $(round(θ_sub_ref, digits=3)))",
            ),
        )
        r_sup = _vcpl04_problem(1.32; n_cells = 40, t_end = 4000.0)
        sol_sup = solve(r_sup.problem; progress = false, abstol = 1e-8, reltol = 1e-6)
        rise_sup = maximum(sol_sup.T) - P.T0
        push!(
            out,
            QoIResult(
                name = "supercritical δ = 1.50 δ_c ignites (thermal runaway)",
                value = rise_sup,
                exact = NaN,
                error = rise_sup > 50.0 ? 0.0 : 1.0,
                tol = 0.5,
                comment = "max rise $(round(rise_sup, digits=1)) K must exceed 50 K " *
                          "(θ_c ≈ 1.19 ⇒ 14.9 K is the largest stable rise)",
            ),
        )
        return out
    end

    # Spatial convergence of the steady center temperature against the
    # (mesh-independent) shooting reference — smooth problem, clean order 2.
    function center_err(sol)
        P = _vcpl04_problem(δ_main; n_cells = 2, t_end = t_main)
        sh = fk_steady_shoot(S0 = P.S0, E = P.E, T0 = P.T0, k = P.k, L = P.L_half)
        z = get_cell_positions(sol)
        jc = argmin(abs.(z .- P.L_half))
        return abs(sol.T[jc, end] - sh.T_at(abs(z[jc] - P.L_half)))
    end

    return VerificationCase(
        id = "V-CPL-04",
        title = "Frank-Kamenetskii slab: subcritical steady state + criticality bracket",
        tags = [:reaction_conduction, :bifurcation, :exothermic, :coupled_jacobian],
        references = [
            "D.A. Frank-Kamenetskii, Diffusion and Heat Exchange in Chemical Kinetics",
            "Verification Suite Plan 2026-07-06, V-CPL-04",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-8, reltol = 1e-6),
        qois = qois,
        ci = CIConfig(
            n_cells = 40,
            budget_s = 30.0,
            convergence = ConvergenceSpec(
                n_cells = [20, 40, 80],
                error_fn = center_err,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
            ),
        ),
        report = (
            n_cells_ladder = [20, 40, 80, 160],
            qoi_curves = (:steady_profiles, :center_T_histories, :bifurcation_diagram),
            curves = sol -> begin
                P = _vcpl04_problem(δ_main; n_cells = 2, t_end = t_main) # params only
                δc = fk_delta_crit()
                z = get_cell_positions(sol)
                jc(zv) = argmin(abs.(zv .- P.L_half))
                shoot(δ) = fk_steady_shoot(
                    S0 = δ * P.k * P.εT / P.L_half^2,
                    E = P.E,
                    T0 = P.T0,
                    k = P.k,
                    L = P.L_half,
                )
                fksolve(δ; t_end) = solve(
                    _vcpl04_problem(δ; n_cells = 40, t_end = t_end, n_outputs = 41).problem;
                    progress = false,
                    abstol = 1e-8,
                    reltol = 1e-6,
                )
                # The assert-mode criticality bracket, re-solved with dense
                # outputs for the history figure.
                sol_sub = fksolve(0.70; t_end = 5000.0)
                sol_sup = fksolve(1.32; t_end = 4000.0)
                z_sub = get_cell_positions(sol_sub)
                z_sup = get_cell_positions(sol_sup)
                sh_main = shoot(δ_main)
                sh_sub = shoot(0.70)
                # δ_c bisection stand-in: the turning-point hunt as a swept
                # diagram. Long windows on every rung (settling slows near the
                # turning point on both sides; runaway rungs saturate at the
                # T_max gate instead of chasing e^θ). Points are classified by
                # OBSERVED outcome, not by δ vs the FK-limit δ_c: at finite
                # E/RT₀ = 40 the true (full-Arrhenius) turning point sits a
                # few percent above 0.8785 — δ = 0.90 settles to a genuine
                # steady state (shooting confirms a root at 0.90, none at
                # 0.92), which is exactly the FK-truncation physics the QoIs
                # quantify.
                δs = [0.2, 0.44, 0.70, 0.82, 0.86, 0.90, 1.00, 1.32]
                θs = map(δs) do δ
                    s = fksolve(δ; t_end = 25000.0)
                    (maximum(s.T[:, end]) - P.T0) / P.εT
                end
                θ_runaway = 3.0   # gate-capped runaway sits at θ ≈ 8
                set = findall(<=(θ_runaway), θs)
                ran = findall(>(θ_runaway), θs)
                [
                    ReportCurve(
                        name = :steady_profiles,
                        title = "V-CPL-04: steady T(z) vs full-Arrhenius shooting reference",
                        xlabel = "z (m)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "δ = 0.44 ≈ 0.50 δ_c",
                                x = z,
                                y = sol.T[:, end],
                                exact = zz ->
                                    sh_main.T_at(min(abs(zz - P.L_half), P.L_half)),
                            ),
                            CurveSeries(
                                label = "δ = 0.70 ≈ 0.80 δ_c",
                                x = z_sub,
                                y = sol_sub.T[:, end],
                                exact = zz ->
                                    sh_sub.T_at(min(abs(zz - P.L_half), P.L_half)),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :center_T_histories,
                        title = "V-CPL-04: center-temperature histories across the " *
                                "criticality bracket",
                        xlabel = "t (s)",
                        ylabel = "T at slab center (K)",
                        series = [
                            CurveSeries(
                                label = "δ = 0.44 (stable)",
                                x = sol.t,
                                y = sol.T[jc(z), :],
                            ),
                            CurveSeries(
                                label = "δ = 0.70 (near-critical, stable)",
                                x = sol_sub.t,
                                y = sol_sub.T[jc(z_sub), :],
                            ),
                            CurveSeries(
                                label = "δ = 1.32 (runaway, T_max-gated)",
                                x = sol_sup.t,
                                y = sol_sup.T[jc(z_sup), :],
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :bifurcation_diagram,
                        title = "V-CPL-04: bifurcation diagram (line: full-Arrhenius " *
                                "shooting; FK-limit δ_c ≈ $(round(δc, digits = 4)))",
                        xlabel = "δ = (E/RT₀²)·S₀L²/k",
                        ylabel = "θ_max = (T_max − T₀)·E/RT₀²",
                        series = [
                            CurveSeries(
                                label = "settled steady state",
                                x = δs[set],
                                y = θs[set],
                                # NaN past the full-Arrhenius turning point
                                # (shooting bisection loses its bracket there)
                                exact = δ -> try
                                    (shoot(δ).T_center - P.T0) / P.εT
                                catch
                                    NaN
                                end,
                            ),
                            CurveSeries(
                                label = "runaway (T_max-gate-capped)",
                                x = δs[ran],
                                y = θs[ran],
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl04)
