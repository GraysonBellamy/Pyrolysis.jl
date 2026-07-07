# =============================================================================
# V-CPL-02 — Landau ablation problem (exact steadily-regressing surface)
# =============================================================================
#
# H.G. Landau, Q. Appl. Math. 8:81 (1950). The flagship ALE case: none of
# FDS/Gpyro/ThermaKin verify against it. Semi-infinite solid, constant net
# surface flux q″, material removed at ablation temperature T_p with heat of
# ablation ΔH_p. Exact steady traveling wave (exact/landau.jl):
#
#   v_ss = q″/(ρ[c(T_p − T₀) + ΔH_p]),  T(ζ) = T₀ + (T_p − T₀)e^(−v_ss ζ/α)
#
# Realization: single reaction solid → gas with γ_gas = 0 (volume loss at
# intrinsic density ⇒ surface regression) under ALE with depletion-driven
# cell merging. The kinetics-controlled front approximates the fixed-T_p
# ablation limit; (A, E) are chosen on the leading-order Laplace locus
#
#   A = v² E (T_p − T₀) e^(E/RT_p) / (α R T_p²)
#
# so the designed regression rate is v_target at ablation temperature ≈ T_p.
# At finite E the front self-adjusts: the primary QoI therefore compares the
# measured regression rate v against the Landau balance evaluated at the
# MEASURED surface temperature — an identity that must hold for the steady
# wave regardless of kinetics — and the design locus is only a loose check.
#
# Gas heat capacity is set EQUAL to the solid's (c_gas = c_solid): then the
# heat of reaction ΔH(T) is T-independent (Kirchhoff) and the traveling-wave
# identity v = v_ss(measured T_s) is exact for ANY kinetics — the regression
# anchor. With c_gas ≈ 0 the identity is NOT exact at finite E: the true wave
# consumes only c(T̄_form − T₀) + ΔH per kg (conversion happens at T̄_form <
# T_s), so v runs ≈ +RT_s²/E·c/(cΔT+ΔH) fast against v_ss(T_s) — that
# configuration belongs to report-mode studies, not to this pinned QoI.
#
# The c_gas = c_solid anchor also guards the relative-velocity ALE energy
# formulation (tech-ref §11.2.2): a spurious full-mesh-velocity T-advection
# term (`+w·∂T/∂z` applied although the mesh rides the collapsing solid)
# would break this identity by ≈ c·(T_s − T̄_form) ≈ RT_s²/E·c per kg — a
# few-percent deficit. For the identity to be a sharp guard the QoI window
# must sit deep in the steady wave: from a uniform-T₀ start the true
# initial-value solution approaches the traveling wave as e^(−t/τ) with
# τ = α/v² ≈ 0.25 s — heat is still being banked into the growing thermal
# tail, so v(t) runs below v_ss(measured T_s) while T_s itself equilibrates
# early. That deficit is the correct transient of the initial-value problem,
# not solver error: the windowed first law (q″Δt = ΔE_stored + ablation
# enthalpy) closes throughout the approach, and the pre-ablation phase
# matches Landau's exact constant-flux heat-up T_s(t) = T₀ + 2q″√(t/(πρck)).
# The steady exact solution is simply not a valid reference before ~6τ. At
# the original window [0.55, 0.95] s (2–4τ) the residual deficit was ≈0.9%
# and mesh-INDEPENDENT (n=60→240 ladder, observed order 0.11 — diagnosed
# 2026-07-07). The window now sits at ≥6τ, where the remaining residual is
# the §12.3.2 half-cell T_s closure with the reaction zone inside the
# boundary half-cell — mesh-DEPENDENT as it should be (it shifts
# v_ss(measured T_s) but not the actual regression rate): ~1.3% at
# Δz ≈ 2δ_r, ~5e-4 at the CI mesh Δz ≈ δ_r — so the pinned tolerance
# actually arms the ALE guard.
#
# Blowing correction stays OFF (exact solution assumes prescribed net flux).
# Depletion limiter stays ON: the tanh roll-off is part of the regression
# physics here (front cells are consumed to zero) and its effect is contained
# in the pinned tolerances.
#
# CI runs no ConvergenceSpec: the case resolves a moving reaction front
# (δ_r ≈ δ·RT_s²/(E(T_s−T₀)) ≈ 53 μm ≈ Δz at the CI mesh) and a meaningful
# order study requires the report-mode Richardson ladder; the assert-mode
# QoIs below are all physics identities with pinned headroom instead.
# V-CPL-02b (broader front, E halved) doubles as the kinetics-limit
# insensitivity check of the ablation-temperature locus.

function _vcpl02_common(; E::Float64, v_target::Float64, n_default::Int)
    ρ, c, k = 1000.0, 1000.0, 1.0
    α = k / (ρ * c)
    T0, Tp, ΔH = 300.0, 700.0, 1.0e6
    q = v_target * ρ * (c * (Tp - T0) + ΔH)
    R = 8.314462618
    A = v_target^2 * E * (Tp - T0) * exp(E / (R * Tp)) / (α * R * Tp^2)
    L = 10e-3            # v·t_end = 4 mm consumed; 6 mm ≈ 12δ remains
    t_end = 2.0
    w1, w2 = 1.5, 1.9    # ≥6 development time constants (τ = α/v² = 0.25 s):
    # the e^(−t/τ) transient is <1e-3 of v_ss by w1 (see header)

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = ρ, c = c, k = k)
        # c_gas = c_solid: constant-ΔH(T) anchor — the Landau identity is
        # exact for any kinetics (see header).
        gas = GaseousComponent(:gas, M = 0.05, c = c, k = 0.02, λ = 1e-4)
        rxn = Reaction(:ablate, 1 => 2; A = A, E = E, h = ΔH)
        mat = Material(name = :VCPL02, components = (solid, gas), reactions = (rxn,))
        mesh = create_uniform_mesh(L, n_cells, Val(2); T_initial = T0, ξ_initial = (ρ, 0.0))
        bc_set = BoundaryConditionSet(
            thermal = (top = HeatFluxBC(q), bottom = AdiabaticBC()),
            mass = (
                top = ConvectiveMassBC(h_m = 10.0, ξ_inf = 0.0),
                bottom = ImpermeableBC(),
            ),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T0,
            ξ_initial = [z -> ρ, z -> 0.0],
            tspan = (0.0, t_end),
            output_times = collect(0.0:0.05:t_end),
        )
    end

    # Window means/differences shared by both cases' QoIs. sol.z stores NODE
    # positions, NaN-padded above the surface after depletion merges.
    function measures(sol)
        th = sol.extras.thickness_history
        i1 = findfirst(>=(w1), sol.t)
        i2 = findlast(<=(w2), sol.t)
        Δt = sol.t[i2] - sol.t[i1]
        v_thick = (th[i1] - th[i2]) / Δt
        v_ggr = sum(@view sol.extras.gas_generation_rate[i1:i2]) / (i2 - i1 + 1) / ρ
        m_escaped = sol.extras.total_mass[i1] - sol.extras.total_mass[i2]
        Ts = sol.surface_T[i2]
        v_landau = landau_vss(q = q, ρ = ρ, c = c, ΔH = ΔH, T0 = T0, Tp = Ts)
        return (; i1, i2, Δt, v_thick, v_ggr, m_escaped, Ts, v_landau)
    end

    return (; ρ, c, k, α, T0, Tp, ΔH, q, A, L, build, measures, v_target, n_default)
end

function _verification_case_vcpl02()
    # n = 160 puts Δz = 62.5 μm ≈ δ_r: the front-in-half-cell T_s bias that
    # dominates the identity residual at coarser meshes (~1.3% at n = 100)
    # is resolved away (~5e-4).
    C = _vcpl02_common(E = 1.0e5, v_target = 2.0e-3, n_default = 160)

    function qois(sol)
        m = C.measures(sol)
        out = QoIResult[]

        err_bal = abs(m.v_thick - m.v_landau) / m.v_landau
        push!(
            out,
            QoIResult(
                name = "steady energy balance: v vs Landau v_ss(measured T_s)",
                value = m.v_thick,
                exact = m.v_landau,
                error = err_bal,
                tol = 1.6e-3,  # 3× observed: 5.2e-4 (n=160, 2026-07-07, window
                # [1.5,1.9] s ≈ 6–7.6τ; residual = front-in-half-cell T_s
                # closure at Δz ≈ δ_r + holdup drift — see case header)
                comment = "traveling-wave identity, exact at any E for c_gas = c_solid",
            ),
        )

        err_int = abs(m.v_thick - m.v_ggr) / m.v_ggr
        push!(
            out,
            QoIResult(
                name = "regression consistency: thickness slope vs gas generation",
                value = m.v_thick,
                exact = m.v_ggr,
                error = err_int,
                tol = 2.1e-3,  # 3× the report-ladder worst: 7.0e-4 (n=240,
                # 2026-07-07); at the CI mesh the observed is 5.6e-6, but the
                # identity's floor is NOT monotone in Δz once Δz < δ_r (n=240
                # rung), so a CI-mesh 3× pin would fail finer meshes. Guard
                # target is O(1) γ_gas bookkeeping breakage — 2.1e-3 suffices.
                comment = "γ_gas = 0 ⇒ dL/dt = −ṁ‴_gas/ρ exactly",
            ),
        )

        # Mass ledger: window mass drop must equal the gas-generation integral
        # (holdup is steady over the plateau; trapezoid on output times).
        ggr = sol.extras.gas_generation_rate
        gen = sum(0.5 * (ggr[j] + ggr[j+1]) * (sol.t[j+1] - sol.t[j]) for j = m.i1:(m.i2-1))
        err_ledger = abs(m.m_escaped - gen) / gen
        push!(
            out,
            QoIResult(
                name = "mass ledger closure over plateau window",
                value = m.m_escaped,
                exact = gen,
                error = err_ledger,
                tol = 4e-3,  # 3× observed: 1.4e-3 (n=160, 2026-07-07, flux-form
                # gas transport; residual = holdup drift incl. surface sweep)
                comment = "total-mass drop vs ∫ ṁ″ dt; residual = holdup drift",
            ),
        )

        # Exponential in-depth profile: fit the decay rate of ln(T − T0) over
        # ζ/δ ∈ [0.75, 1.7] — below the reaction zone but inside the region the
        # tail has had time to equilibrate (a depth ζ needs ~ζ²/α; beyond
        # √(α·t) ≈ 2.7δ at the window end the tail is still developing and
        # reads as a steeper decay). Slope must equal −v/α.
        j = m.i2
        nodes = filter(!isnan, @view sol.z[:, j])
        δ = C.α / m.v_thick
        ζs = Float64[]
        ys = Float64[]
        for i = 1:(length(nodes)-1)
            ζ = nodes[end] - 0.5 * (nodes[i] + nodes[i+1])
            if 0.75 * δ <= ζ <= 1.7 * δ
                push!(ζs, ζ)
                push!(ys, log(sol.T[i, j] - C.T0))
            end
        end
        n = length(ζs)
        ζ̄ = sum(ζs) / n
        ȳ = sum(ys) / n
        slope = sum((ζs[i] - ζ̄) * (ys[i] - ȳ) for i = 1:n) / sum((ζ - ζ̄)^2 for ζ in ζs)
        v_slope = -slope * C.α
        err_prof = abs(v_slope - m.v_thick) / m.v_thick
        push!(
            out,
            QoIResult(
                name = "in-depth profile decay rate (surface-attached frame)",
                value = v_slope,
                exact = m.v_thick,
                error = err_prof,
                tol = 2.7e-2,  # 3× observed: 8.9e-3 (n=160, 2026-07-07):
                # reaction-zone edge + residual tail development
                comment = "ln(T−T₀) slope over ζ/δ ∈ [0.75, 1.7], $(n) cells",
            ),
        )

        err_locus = abs(m.v_thick - C.v_target) / C.v_target
        push!(
            out,
            QoIResult(
                name = "design-locus regression rate (loose)",
                value = m.v_thick,
                exact = C.v_target,
                error = err_locus,
                tol = 0.10,
                comment = "leading-order Laplace locus only pins T_p to O(RT_p²/E); " *
                          "informational, tightens as E ↑ (see V-CPL-02b)",
            ),
        )

        ΔT_back = sol.T[1, m.i2] - C.T0
        push!(
            out,
            QoIResult(
                name = "back-face rise (semi-infinite validity)",
                value = ΔT_back,
                exact = 0.0,
                error = abs(ΔT_back),
                tol = 0.01,  # observed 3.6e-4 K (n=160, 2026-07-07); loose vs
                # 3× — the guard is against gross domain-too-short violation
                comment = "guards the semi-infinite assumption",
            ),
        )
        return out
    end

    return VerificationCase(
        id = "V-CPL-02",
        title = "Landau ablation: steady surface regression under ALE",
        tags = [:ale, :moving_boundary, :kinetics, :surface_regression, :coupled],
        references = [
            "H.G. Landau, Q. Appl. Math. 8:81 (1950)",
            "Verification Suite Plan 2026-07-06, V-CPL-02",
        ],
        build = C.build,
        solve_kwargs = (
            use_ale = true,
            min_thickness = 1e-4,
            handle_depletion = true,
            depletion_threshold = 0.05,
            abstol = 1e-7,
            reltol = 1e-5,
        ),
        qois = qois,
        ci = CIConfig(n_cells = 160, budget_s = 110.0),
        report = (
            n_cells_ladder = [60, 100, 160, 240],
            E_ladder = [5e4, 1e5, 1.5e5, 2e5],
            qoi_curves = (:mlr_history, :thickness_history, :profile_surface_frame),
            curves = sol -> begin
                m = C.measures(sol)
                # surface-attached frame at the plateau-window end (see qois)
                j = m.i2
                nodes = filter(!isnan, sol.z[:, j])
                nc = length(nodes) - 1
                ζ = [nodes[end] - 0.5 * (nodes[i] + nodes[i+1]) for i = 1:nc]
                [
                    ReportCurve(
                        name = :mlr_history,
                        title = "V-CPL-02: MLR approach to the Landau steady wave",
                        xlabel = "t (s)",
                        ylabel = "ṁ″ (kg/m²s)",
                        series = [
                            CurveSeries(
                                label = "gas generation rate",
                                x = sol.t,
                                y = sol.extras.gas_generation_rate,
                                # steady-wave identity at the MEASURED T_s
                                exact = t -> C.ρ * m.v_landau,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :thickness_history,
                        title = "V-CPL-02: thickness regression (slope −v_ss)",
                        xlabel = "t (s)",
                        ylabel = "thickness (m)",
                        series = [
                            CurveSeries(
                                label = "thickness",
                                x = sol.t,
                                y = sol.extras.thickness_history,
                                # tangent through the plateau-window start with
                                # the Landau slope
                                exact = t ->
                                    sol.extras.thickness_history[m.i1] -
                                    m.v_landau * (t - sol.t[m.i1]),
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :profile_surface_frame,
                        title = "V-CPL-02: in-depth profile in the surface-attached frame",
                        xlabel = "ζ = depth below moving surface (m)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "t = $(round(sol.t[j], digits = 2)) s",
                                x = ζ,
                                y = sol.T[1:nc, j],
                                exact = zz -> landau_T(
                                    zz;
                                    T0 = C.T0,
                                    Tp = m.Ts,
                                    v = m.v_thick,
                                    α = C.α,
                                ),
                            ),
                        ],
                    ),
                ]
            end,
            # Richardson ladder over the steady-wave identities (the order
            # study the assert-mode header defers to report mode).
            error_norms = sol -> begin
                m = C.measures(sol)
                (
                    balance = abs(m.v_thick - m.v_landau) / m.v_landau,
                    consistency = abs(m.v_thick - m.v_ggr) / m.v_ggr,
                )
            end,
            # Kinetics-limit insensitivity: the Landau balance must close at
            # every E on the constant-ablation-temperature locus — as long as
            # the front stays resolved: δ_r ∝ 1/E, so at n = 160 the E = 2e5
            # point has δ_r ≈ 26 μm < Δz and its ~1% miss is the §12.3.2
            # front-in-half-cell T_s bias (resolution, not physics; 2026-07-07).
            param_sweep = (
                name = :E,
                values = [5e4, 1e5, 1.5e5, 2e5],
                build = (n, E) ->
                    _vcpl02_common(E = E, v_target = 2.0e-3, n_default = n).build(n),
                measure = sol -> begin
                    m = C.measures(sol)  # q, ρ, c, ΔH, T0 are E-independent
                    (
                        balance_err = abs(m.v_thick - m.v_landau) / m.v_landau,
                        consistency_err = abs(m.v_thick - m.v_ggr) / m.v_ggr,
                        v_thick = m.v_thick,
                    )
                end,
            ),
        ),
    )
end

function _verification_case_vcpl02b()
    # Kinetics-limit companion: E halved ⇒ reaction zone twice as wide
    # (δ_r/δ = RT_s²/(E·(T_s−T₀)) ≈ 0.09), resolvable on a coarser mesh. The
    # Landau balance with measured T_s must close at ANY E — comparing the
    # closure here and in V-CPL-02 is the ablation-limit insensitivity check.
    # n = 100 keeps Δz = 100 μm on the L = 10 mm domain (the pre-window-move
    # pin used n = 60 on L = 6 mm, the same spacing).
    C = _vcpl02_common(E = 5.0e4, v_target = 2.0e-3, n_default = 100)

    function qois(sol)
        m = C.measures(sol)
        err_bal = abs(m.v_thick - m.v_landau) / m.v_landau
        err_int = abs(m.v_thick - m.v_ggr) / m.v_ggr
        return [
            QoIResult(
                name = "steady energy balance: v vs Landau v_ss(measured T_s)",
                value = m.v_thick,
                exact = m.v_landau,
                error = err_bal,
                tol = 2e-3,  # 3× observed: 6.8e-4 (n=100, 2026-07-07; Δz ≈ δ_r
                # here since the E/2 front is twice as wide — matches the
                # ~5e-4 closure of V-CPL-02 at ITS δ_r-resolving mesh)
                comment = "identity must close at E = 50 kJ/mol as at 100 (locus check)",
            ),
            QoIResult(
                name = "regression consistency: thickness slope vs gas generation",
                value = m.v_thick,
                exact = m.v_ggr,
                error = err_int,
                tol = 2.5e-3,  # 3× observed: 8.3e-4 (n=100, 2026-07-07)
                comment = "γ_gas = 0 ⇒ dL/dt = −ṁ‴_gas/ρ exactly",
            ),
        ]
    end

    return VerificationCase(
        id = "V-CPL-02b",
        title = "Landau ablation: kinetics-limit insensitivity variant (E/2)",
        tags = [:ale, :moving_boundary, :kinetics, :surface_regression],
        references = ["H.G. Landau, Q. Appl. Math. 8:81 (1950)"],
        build = C.build,
        solve_kwargs = (
            use_ale = true,
            min_thickness = 1e-4,
            handle_depletion = true,
            depletion_threshold = 0.05,
            abstol = 1e-7,
            reltol = 1e-5,
        ),
        qois = qois,
        ci = CIConfig(n_cells = 100, budget_s = 15.0),
        report = (;),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl02)
push!(VERIFICATION_CASE_BUILDERS, _verification_case_vcpl02b)
