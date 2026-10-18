# =============================================================================
# V-CPL-01 — Neumann two-phase Stefan problem (latent-heat front)
# =============================================================================
#
# Gpyro's flagship moving-front verification (Lautenberger §3.4.3, Eqs.
# 3.137–3.138; Carslaw & Jaeger §11.2). Liquid initially at T₀ > T_m; wall
# (top face) dropped to T∞ < T_m at t = 0; a solidification front
# s(t) = 2λ√(αₛt) propagates inward releasing ΔH_m. Exact solution and
# eigenvalue in exact/stefan.jl (λ = 0.2317 for the dissertation parameters —
# self-tested against the published value and an independent quadrature
# energy balance).
#
# Two solver realizations:
#
# Route A (this case): apparent heat capacity — Gaussian c_p peak of variance
#   σ² = 0.1 K² (Gpyro Eq. 3.5c with σ_m² = 0.1) carrying the full ΔH_m,
#   exercising the T-dependent property path + property-derivative Jacobian
#   through a 38× spike in c. Uses the dissertation's polypropylene-like
#   parameter set verbatim, so profiles compare to dissertation Fig. 3.11.
#
# Route B (V-CPL-01b): latent heat via a T_max-gated reaction melt → solid
#   (h = −ΔH_m: solidification releases heat; our h > 0 = endothermic).
#   The temperature gate is a tanh ramp of configurable width
#   (`Reaction(...; gate_width)` [K], default 1.0; src/Physics/kinetics.jl).
#   At the default 1 K width the dissertation's 5 K driving ΔT shifts the
#   effective freezing point by ~20% of the driving force (measured: front
#   19% slow), so this CI case rescales to ±50 K driving ΔT (own exact λ
#   from the same library), making the gate a ~1–2% effect. The dissertation
#   parameter set with a narrow gate (e.g. 0.1 K) and the gate-width → 0
#   convergence sweep belong in report mode.
#
# Geometry: cold Dirichlet wall on TOP (z = L; depth x = L − z), adiabatic
# bottom, L = 4 cm ≫ liquid penetration √(4α_l t) ≈ 1.25 cm (back-face guard
# asserts undisturbedness). Geometrically stretched mesh, finest at the wall
# (Δz_top ≈ 14 μm, ≈ 103 μm at the deepest front position ≈ 3 mm — Gpyro
# used uniform 50 μm).
#
# No CI ConvergenceSpec: with fixed peak width σ (route A) or gate width
# (route B) the solver converges to the smeared problem, whose offset from
# the sharp-front Stefan solution dominates the mesh error — a meaningful
# order study is the σ → 0 + mesh ladder (a report-mode peak-width
# convergence study).

# Wall-depth frame shared by the QoIs and the report curves: depth x = L − z
# (cold wall on top), T interpolated at fixed depth, and the T_m-isotherm
# front position (linear interpolation on the bracketing cells).
function _vcpl01_frame(sol, L, T_m)
    z = cell_positions(sol)
    x = L .- z
    ord = sortperm(x)
    xs = x[ord]

    Tnum(depth, j) = begin
        i2i = findfirst(>=(depth), xs)
        x1, x2 = xs[i2i-1], xs[i2i]
        T1, T2 = sol.T[ord[i2i-1], j], sol.T[ord[i2i], j]
        T1 + (T2 - T1) * (depth - x1) / (x2 - x1)
    end

    front(j) = begin
        for i = 1:(length(xs)-1)
            Ta, Tb = sol.T[ord[i], j], sol.T[ord[i+1], j]
            if (Ta - T_m) * (Tb - T_m) <= 0 && Ta != Tb
                return xs[i] + (xs[i+1] - xs[i]) * (T_m - Ta) / (Tb - Ta)
            end
        end
        return NaN
    end

    return (; xs, ord, Tnum, front)
end

function _vcpl01_qois(p, λ, L, t_end, sol; tols_hist, tol_front, tol_slope, tol_back)
    αs = α_solid(p)
    F = _vcpl01_frame(sol, L, p.T_m)
    Tnum = F.Tnum
    ord = F.ord

    out = QoIResult[]

    # 1. Temperature histories at the four dissertation-figure depths
    # (t ≥ t_end/3 — earlier the front is inside the first few cells).
    depths = (0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3)
    for (d_i, depth) in enumerate(depths)
        err = maximum(
            abs(Tnum(depth, j) - stefan_T(depth, t, p, λ)) for
            (j, t) in enumerate(sol.t) if t >= t_end / 3
        )
        push!(
            out,
            QoIResult(
                name = "T-history L∞ error at x = $(depth*1e3) mm",
                value = err,
                error = err,
                tol = tols_hist[d_i],
                comment = "vs exact two-phase profile, t ≥ $(round(t_end/3)) s",
            ),
        )
    end

    # 2. Front position (T_m isotherm crossing, linear interpolation).
    front_pos = F.front
    for frac in (1 / 3, 2 / 3, 1.0)
        t_chk = frac * t_end
        j = findmin(abs.(sol.t .- t_chk))[2]
        s_num = front_pos(j)
        s_ex = stefan_front(sol.t[j]; λ = λ, αs = αs)
        err = abs(s_num - s_ex) / s_ex
        push!(
            out,
            QoIResult(
                name = "front position at t = $(round(sol.t[j])) s",
                value = s_num,
                exact = s_ex,
                error = err,
                tol = tol_front,
                comment = "s(t) = 2λ√(αₛt), λ = $(round(λ, digits = 5))",
            ),
        )
    end

    # 3. √t growth law: least-squares fit of s = m√t over the late window.
    js = [j for (j, t) in enumerate(sol.t) if t >= t_end / 3]
    ss = [front_pos(j) for j in js]
    ts = [sol.t[j] for j in js]
    m_fit = sum(ss .* sqrt.(ts)) / sum(ts)
    m_ex = 2 * λ * sqrt(αs)
    err_slope = abs(m_fit - m_ex) / m_ex
    push!(
        out,
        QoIResult(
            name = "√t front-growth coefficient",
            value = m_fit,
            exact = m_ex,
            error = err_slope,
            tol = tol_slope,
            comment = "fit over t ∈ [$(round(t_end/3)), $(t_end)] s",
        ),
    )

    # 4. Back face undisturbed (semi-infinite validity of the reference).
    dev_back = abs(sol.T[ord[end], end] - p.T0)
    push!(
        out,
        QoIResult(
            name = "back-face disturbance (semi-infinite validity)",
            value = dev_back,
            exact = 0.0,
            error = dev_back,
            tol = tol_back,
            comment = "L = 4 cm vs penetration ≈ 1.25 cm",
        ),
    )
    return out
end

# Report-mode overlay curves shared by both routes: the dissertation-figure
# depth histories vs the exact two-phase profile, and the T_m-isotherm front
# trajectory vs s(t) = 2λ√(αₛt). `id` labels the figure titles.
function _vcpl01_report_curves(p, λ, L, id)
    return sol -> begin
        F = _vcpl01_frame(sol, L, p.T_m)
        αs = α_solid(p)
        js = [j for (j, t) in enumerate(sol.t) if t > 0]
        ts = sol.t[js]
        front = [F.front(j) for j in js]
        keep = findall(!isnan, front)   # pre-detection outputs (front inside wall cell)
        # Gpyro 0.8200 overlay (reference/gpyro/): dissertation §3.4.3 deck —
        # the identical Stefan setup, latent heat via Gpyro's apparent-c_p
        # Gaussian (σₘ² = 0.1), Δz = 0.05 mm. Long-format depth histories +
        # T_m-isotherm front trajectory.
        gp_t, gp_x, gp_T = reference_columns(
            "gpyro/vcpl01_result.csv",
            "t_s",
            "x_m",
            "T_K";
            header_rows = 1,
        )
        gp_ft, gp_fs = reference_columns(
            "gpyro/vcpl01_front.csv",
            "t_s",
            "s_front_m";
            header_rows = 1,
        )
        [
            ReportCurve(
                name = :T_histories_at_depths,
                title = "$(id): T histories at the dissertation depths (cf. Fig. 3.11)",
                xlabel = "t (s)",
                ylabel = "T (K)",
                series = vcat(
                    [
                        CurveSeries(
                            label = "x = $(d * 1e3) mm",
                            x = ts,
                            y = [F.Tnum(d, j) for j in js],
                            exact = t -> stefan_T(d, t, p, λ),
                        ) for d in (0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3)
                    ],
                    [
                        let rows = findall(x -> abs(x - d) < 1e-8, gp_x)
                            CurveSeries(
                                label = "Gpyro 0.8200, x = $(d * 1e3) mm",
                                crosscode = true,
                                x = gp_t[rows],
                                y = gp_T[rows],
                                exact = t -> stefan_T(d, t, p, λ),
                            )
                        end for d in (0.5e-3, 1.0e-3, 1.5e-3, 2.0e-3)
                    ],
                ),
            ),
            ReportCurve(
                name = :front_trajectory,
                title = "$(id): solidification front vs s(t) = 2λ√(αₛt)",
                xlabel = "t (s)",
                ylabel = "front depth (m)",
                series = [
                    CurveSeries(
                        label = "T_m isotherm",
                        x = ts[keep],
                        y = front[keep],
                        exact = t -> stefan_front(t; λ = λ, αs = αs),
                    ),
                    CurveSeries(
                        label = "Gpyro 0.8200 (T_m isotherm)",
                        crosscode = true,
                        x = gp_ft,
                        y = gp_fs,
                        exact = t -> stefan_front(t; λ = λ, αs = αs),
                    ),
                ],
            ),
        ]
    end
end

# Relative error of the √t front-growth coefficient (the case's QoI 3) — the
# scalar measure for the report-mode parameter sweeps (σ² and A ladders).
function _vcpl01_front_slope_err(p, λ, L, t_end)
    return sol -> begin
        F = _vcpl01_frame(sol, L, p.T_m)
        js = [j for (j, t) in enumerate(sol.t) if t >= t_end / 3]
        m_fit = sum(F.front(j) * sqrt(sol.t[j]) for j in js) / sum(sol.t[j] for j in js)
        m_ex = 2 * λ * sqrt(α_solid(p))
        return (front_slope_relerr = abs(m_fit - m_ex) / m_ex,)
    end
end

function _verification_case_vcpl01()
    p = stefan_gpyro_problem()
    λ = stefan_lambda(p)
    L = 0.04
    t_end = 300.0
    σ2 = 0.1  # Gpyro Eq. 3.5c melting-peak variance [K²]

    function build_with_σ2(n_cells, σ2v)
        c_app = CallableProperty(
            T -> p.c_s + p.ΔH_m / sqrt(2π * σ2v) * exp(-(T - p.T_m)^2 / (2 * σ2v)),
        )
        comp = SolidComponent(:pp, ρ = p.ρ, c = c_app, k = p.k_s)
        mat = Material(name = :VCPL01, components = (comp,), reactions = ())
        # Matched-resolution rung (n ≥ 600): replicate Gpyro's uniform
        # Δz = 0.05 mm dissertation grid instead of stretching — the 1.03
        # geometric law at n = 800 puts the wall cell at ~6e-14 m (integrator
        # dies at t = 0). No convergence fit runs on this case, so the mixed
        # mesh family is harmless.
        mesh = if n_cells >= 600
            create_uniform_mesh(L, n_cells, Val(1); T_initial = p.T0, ξ_initial = (p.ρ,))
        else
            create_stretched_mesh(
                L,
                n_cells,
                Val(1),
                1.03;
                T_initial = p.T0,
                ξ_initial = (p.ρ,),
                stretch_from = :top,
            )
        end
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(p.T∞), bottom = AdiabaticBC()),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> p.T0,
            ξ_initial = [z -> p.ρ],
            tspan = (0.0, t_end),
            output_times = collect(range(0.0, t_end, length = 31)),
        )
    end
    build(n_cells) = build_with_σ2(n_cells, σ2)

    qois(sol) = _vcpl01_qois(
        p,
        λ,
        L,
        t_end,
        sol;
        # 3× observed: 0.014/0.027/0.038/0.134 K (n=150)
        tols_hist = (0.05, 0.1, 0.12, 0.4),
        tol_front = 0.03,   # 3× observed max 0.98%
        tol_slope = 0.03,   # 3× observed 0.95%
        tol_back = 0.01,    # observed 1e-4 K
    )

    return VerificationCase(
        id = "V-CPL-01",
        title = "Neumann Stefan problem via apparent heat capacity (Gpyro params)",
        tags = [:moving_front, :latent_heat, :state_dependent_properties, :coupled],
        references = [
            "Carslaw & Jaeger (1959) §11.2",
            "Lautenberger dissertation (2007) §3.4.3, Eqs. 3.137–3.138, Fig. 3.11",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-8, reltol = 1e-6),
        qois = qois,
        ci = CIConfig(n_cells = 150, budget_s = 35.0),
        report = (
            n_cells_ladder = [75, 150, 300, 800],
            qoi_curves = (:T_histories_at_depths, :front_trajectory),
            curves = _vcpl01_report_curves(p, λ, L, "V-CPL-01"),
            # Peak-width study the case header defers to report mode: the
            # smeared-problem offset from the sharp-front Stefan solution,
            # measured on the √t front-growth coefficient (σ² → 0 limit).
            param_sweep = (
                name = :sigma2,
                values = [0.4, 0.1, 0.025],
                build = (n, σ2v) -> build_with_σ2(n, σ2v),
                measure = _vcpl01_front_slope_err(p, λ, L, t_end),
            ),
        ),
    )
end

function _verification_case_vcpl01b()
    # Reaction-route variant, rescaled to ±50 K driving ΔT so the solver's
    # fixed ~1 K tanh gate ramp is a ~2% effect (see case header).
    p = StefanProblem(
        k_s = 0.20,
        k_l = 0.20,
        ρ = 900.0,
        c_s = 1700.0,
        c_l = 1700.0,
        T_m = 433.0,
        ΔH_m = 200.0e3,
        T0 = 483.0,
        T∞ = 383.0,
    )
    λ = stefan_lambda(p)
    L = 0.04
    t_end = 150.0

    function build_with_A(n_cells, A_rate)
        melt = SolidComponent(:melt, ρ = p.ρ, c = p.c_l, k = p.k_l)
        solidc = SolidComponent(:solid, ρ = p.ρ, c = p.c_s, k = p.k_s)
        # Freezing releases heat: h = −ΔH_m (our h > 0 = endothermic). Rate
        # A = 10 s⁻¹ converts within ≪ the front residence time per cell.
        rxn = Reaction(:freeze, 1 => 2; A = A_rate, E = 0.0, h = -p.ΔH_m, T_max = p.T_m)
        mat = Material(
            name = :VCPL01b,
            components = (melt, solidc),
            reactions = (rxn,),
            depletion_limiter = DepletionLimiter(enabled = false),
        )
        # Same matched-resolution rule as route A (uniform grid for n ≥ 600).
        mesh = if n_cells >= 600
            create_uniform_mesh(L, n_cells, Val(2); T_initial = p.T0, ξ_initial = (p.ρ, 0.0))
        else
            create_stretched_mesh(
                L,
                n_cells,
                Val(2),
                1.03;
                T_initial = p.T0,
                ξ_initial = (p.ρ, 0.0),
                stretch_from = :top,
            )
        end
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(p.T∞), bottom = AdiabaticBC()),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> p.T0,
            ξ_initial = [z -> p.ρ, z -> 0.0],
            tspan = (0.0, t_end),
            output_times = collect(range(0.0, t_end, length = 31)),
        )
    end
    build(n_cells) = build_with_A(n_cells, 10.0)

    qois(sol) = _vcpl01_qois(
        p,
        λ,
        L,
        t_end,
        sol;
        # 3× observed: 0.17/0.35/0.61/1.23 K (n=150)
        tols_hist = (0.55, 1.1, 1.9, 3.8),
        tol_front = 0.045,  # 3× observed max 1.4% (gate-ramp shift dominates)
        tol_slope = 0.035,  # 3× observed 1.2%
        tol_back = 0.01,    # observed ≈ 0 K
    )

    return VerificationCase(
        id = "V-CPL-01b",
        title = "Neumann Stefan problem via gated reaction (±50 K rescaled)",
        tags = [:moving_front, :latent_heat, :kinetics, :temperature_gate],
        references = ["Carslaw & Jaeger (1959) §11.2"],
        build = build,
        solve_kwargs = (abstol = 1e-8, reltol = 1e-6),
        qois = qois,
        ci = CIConfig(n_cells = 150, budget_s = 35.0),
        report = (
            n_cells_ladder = [75, 150, 300, 800],
            qoi_curves = (:T_histories_at_depths, :front_trajectory),
            curves = _vcpl01_report_curves(p, λ, L, "V-CPL-01b"),
            # Rate-insensitivity of the gated-reaction route: the front-growth
            # coefficient must be flat in A once conversion is much faster
            # than the per-cell front residence time.
            param_sweep = (
                name = :A_rate,
                values = [3.0, 10.0, 30.0],
                build = (n, A_rate) -> build_with_A(n, A_rate),
                measure = _vcpl01_front_slope_err(p, λ, L, t_end),
            ),
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl01)
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl01b)
