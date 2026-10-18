# =============================================================================
# V-RAD-01 — In-depth absorption with conduction (eigenfunction exact solution)
# =============================================================================
#
# Lautenberger dissertation §3.4.2 + Appendix A: finite slab δ = 1 cm,
# in-depth source q₀κe^(−κx) with κ = 1000 m⁻¹ (κδ = 10), convective loss at
# the irradiated face, insulated back. Exact eigenfunction solution re-derived
# in exact/indepth_absorption.jl (PDE-residual self-tested).
#
# This is the definitive BEER_LAMBERT + conduction + ConvectiveBC coupling
# case. The existing Beer–Lambert test checks the absorption profile and
# energy split WITHOUT conduction; the physics marker here is the SUB-SURFACE
# temperature peak (the peak sits at the front face at t → 0⁺ and migrates
# toward the insulated back as t → ∞), which only the coupled problem shows.
#
# Configuration pitfalls encoded below:
# - `RadiativeFluxBC(flux = q₀, absorptivity = 1.0)`: the EXPLICIT
#   absorptivity bypasses the Kirchhoff ε_eff rescaling at volumetric
#   injection, so exactly q₀ enters the Beer–Lambert sweep.
# - Under BEER_LAMBERT the incident flux contributes zero to the surface
#   Newton balance (deposited volumetrically); convective loss stays at the
#   surface — exactly the exact solution's BC structure. No RadiativeBC is
#   present, so there is no re-radiation term (the exact model has none).
# - The component absorption kwarg α is a MASS coefficient [m²/kg]:
#   α = κ/ρ = 1.0 m²/kg gives the target κ = 1000 m⁻¹.
# - The transmitted fraction e^(−κδ) ≈ 4.5e-5 of q₀ exits the insulated back
#   in both the exact model and the solver (never deposited): consistent.
#
# k, ρ, c, h, q₀ follow the Gpyro §3.4.1 property set for cross-code
# comparability of the dissertation figures (report mode).

function _verification_case_vrad01()
    ρ, c, k = 1000.0, 1400.0, 0.2
    α = k / (ρ * c)
    δ = 0.01
    κ = 1000.0
    q0, h, T0 = 25.0e3, 20.0, 300.0
    times = [30.0, 60.0, 90.0, 180.0]

    ex = make_indepth_absorption_exact(δ = δ, κ = κ, k = k, h = h, α = α, q0 = q0, T0 = T0)

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = ρ, c = c, k = k, α = κ / ρ)
        mat = Material(name = :VRAD01, components = (solid,), reactions = ())
        mesh = create_uniform_mesh(δ, n_cells, Val(1); T_initial = T0, ξ_initial = (ρ,))
        bc_set = BoundaryConditionSet(
            thermal = (
                top = RadiativeFluxBC(flux = q0, absorptivity = 1.0) +
                      ConvectiveBC(h = h, T_inf = T0),
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
        z = cell_positions(sol)
        for (j, t) in enumerate(times)
            Tn = @view sol.T[:, j+1]
            err = maximum(abs(Tn[i] - ex.T_exact(δ - z[i], t)) for i in eachindex(z))
            push!(
                out,
                QoIResult(
                    name = "L∞ T-profile error, t=$(t) s",
                    value = err,
                    error = err,
                    tol = 0.25,  # 3× observed max: 0.080 K at t=180 s (n=100)
                    comment = "observed 0.044/0.057/0.065/0.080 K at t=30/60/90/180 s",
                ),
            )
        end
        for (j, t) in enumerate(times)
            Ts_num = sol.surface_T[j+1]
            Ts_exact = ex.T_exact(0.0, t)
            push!(
                out,
                QoIResult(
                    name = "front-face temperature, t=$(t) s",
                    value = Ts_num,
                    exact = Ts_exact,
                    error = abs(Ts_num - Ts_exact),
                    tol = 0.25,  # 3× observed max: 0.080 K at t=30 s (n=100)
                    comment = "convective-only surface balance under BEER_LAMBERT",
                ),
            )
        end

        # Sub-surface peak at t = 90 s: value and location (location compared
        # against the cell-center quantization Δz).
        j90 = 4  # output_times = [0, 30, 60, 90, 180]
        Tn = @view sol.T[:, j90]
        i_peak = argmax(Tn)
        n = length(z)
        Δz = δ / n
        xs_fine = range(0.0, δ, length = 20001)
        T_fine = [ex.T_exact(x, times[3]) for x in xs_fine]
        i_ex = argmax(T_fine)
        push!(
            out,
            QoIResult(
                name = "sub-surface peak value, t=90 s",
                value = Tn[i_peak],
                exact = T_fine[i_ex],
                error = abs(Tn[i_peak] - T_fine[i_ex]),
                tol = 0.05,  # 3× observed 0.012 K (n=100)
                comment = "the coupled-physics marker no uncoupled test checks",
            ),
        )
        push!(
            out,
            QoIResult(
                name = "sub-surface peak location, t=90 s",
                value = δ - z[i_peak],           # depth below front face
                exact = xs_fine[i_ex],
                error = abs((δ - z[i_peak]) - xs_fine[i_ex]),
                tol = 1.5 * Δz,  # argmax is cell-center quantized
                comment = "peak must sit strictly below the surface",
            ),
        )
        # Interior peak sanity: not at either boundary cell.
        interior = (i_peak != 1 && i_peak != n) ? 0.0 : 1.0
        push!(
            out,
            QoIResult(
                name = "peak strictly interior, t=90 s",
                value = Float64(i_peak),
                error = interior,
                tol = 0.5,
                comment = "argmax at a boundary cell would mean no sub-surface peak",
            ),
        )
        return out
    end

    function l2_error_t90(sol)
        z = cell_positions(sol)
        Tn = @view sol.T[:, 4]
        Δz = δ / length(z)
        return sqrt(
            sum((Tn[i] - ex.T_exact(δ - z[i], 90.0))^2 * Δz for i in eachindex(z)) / δ,
        )
    end

    return VerificationCase(
        id = "V-RAD-01",
        title = "in-depth absorption + conduction (sub-surface peak)",
        tags = [:radiation, :beer_lambert, :convective_bc, :coupled],
        references = [
            "Lautenberger dissertation (2007) §3.4.2 + Appendix A",
            "exact/indepth_absorption.jl (re-derivation, PDE-residual self-test)",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8, radiation_model = BEER_LAMBERT),
        qois = qois,
        ci = CIConfig(
            n_cells = 100,
            budget_s = 20.0,
            convergence = ConvergenceSpec(
                # observed rates 1.996, 1.998
                n_cells = [50, 100, 200],
                error_fn = l2_error_t90,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (
                    abstol = 1e-12,
                    reltol = 1e-10,
                    radiation_model = BEER_LAMBERT,
                ),
            ),
        ),
        report = (
            n_cells_ladder = [25, 50, 100, 200, 400],
            profile_times = times,
            qoi_curves = (:T_profiles, :surface_T_history, :peak_trajectory),
            solve_kwargs = (abstol = 1e-12, reltol = 1e-10, radiation_model = BEER_LAMBERT),
            curves = sol -> begin
                z = cell_positions(sol)
                xs_fine = range(0.0, δ, length = 2001)
                peak_depth_exact(t) =
                    xs_fine[argmax([ex.T_exact(x, t) for x in xs_fine])]
                # Gpyro 0.8200 overlay (reference/gpyro/): dissertation §3.4.2
                # deck (identical parameters — the same eigenfunction case),
                # in-depth Beer–Lambert via KAPPA, Δz = 0.1 mm, dt = 0.1 s.
                gp_t, gp_x, gp_T = reference_columns(
                    "gpyro/vrad01_result.csv",
                    "t_s",
                    "x_m",
                    "T_K";
                    header_rows = 1,
                )
                gp_ts, gp_Ts = reference_columns(
                    "gpyro/vrad01_surfaceT.csv",
                    "t_s",
                    "T_front_K";
                    header_rows = 1,
                )
                [
                    ReportCurve(
                        name = :T_profiles,
                        title = "V-RAD-01: in-depth absorption profiles (sub-surface peak)",
                        xlabel = "depth x (m)",
                        ylabel = "T (K)",
                        series = vcat(
                            [
                                CurveSeries(
                                    label = "t = $(t) s",
                                    x = δ .- z,
                                    y = sol.T[:, j+1],
                                    exact = x -> ex.T_exact(x, t),
                                ) for (j, t) in enumerate(times)
                            ],
                            [
                                let rows = findall(τ -> abs(τ - t) < 1e-6, gp_t)
                                    CurveSeries(
                                        label = "Gpyro 0.8200, t = $(t) s",
                                        crosscode = true,
                                        x = gp_x[rows],
                                        y = gp_T[rows],
                                        exact = x -> ex.T_exact(x, t),
                                    )
                                end for t in times
                            ],
                        ),
                    ),
                    ReportCurve(
                        name = :surface_T_history,
                        title = "V-RAD-01: front-face temperature",
                        xlabel = "t (s)",
                        ylabel = "T_s (K)",
                        series = [
                            CurveSeries(
                                label = "front face",
                                x = sol.t[2:end],
                                y = sol.surface_T[2:end],
                                exact = t -> ex.T_exact(0.0, t),
                            ),
                            CurveSeries(
                                label = "Gpyro 0.8200",
                                crosscode = true,
                                x = gp_ts[2:end],
                                y = gp_Ts[2:end],
                                exact = t -> ex.T_exact(0.0, t),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :peak_trajectory,
                        title = "V-RAD-01: sub-surface peak migration toward the back",
                        xlabel = "t (s)",
                        ylabel = "peak depth below surface (m)",
                        series = [
                            CurveSeries(
                                label = "argmax of cell profile (Δz-quantized)",
                                x = sol.t[2:end],
                                y = [δ - z[argmax(sol.T[:, j])] for j = 2:length(sol.t)],
                                exact = peak_depth_exact,
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol -> begin
                z = cell_positions(sol)
                e = [sol.T[i, 4] - ex.T_exact(δ - z[i], 90.0) for i in eachindex(z)]
                (L2 = sqrt(sum(abs2, e) / length(e)), Linf = maximum(abs, e))
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad01)
