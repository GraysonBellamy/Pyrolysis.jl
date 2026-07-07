# =============================================================================
# V-RX-02 — nth-order and contracting-geometry isothermal decay
# =============================================================================
#
# Khawam & Flanagan (2006) exact isothermal integrations (FDS verifies its
# `part_baking_soda` case against the same forms). Sub-cases:
#
#   a) n = 0    — linear decay, exact finite burnout time. Compared only up to
#                 0.9·t_b with the depletion limiter DISABLED: the solver
#                 clamps ξ ≥ 0 inside ξⁿ, so at ξ = 0 the n = 0 rate stays A
#                 and would integrate ξ negative (kinetics.jl `_fast_pow`).
#   b) n = 2    — ThermaKin's second-order case as single-reactant n = 2:
#                 ξ = ξ₀/(1 + ξ₀kt).
#   c) n = 2/3  — contracting-sphere R3 model with FDS's NaHCO₃ parameters
#                 A_pub = 3.4e11 1/s, E = 103 kJ/mol at T = 450 K; the marker
#                 is m^(1/3) linear in t.
#   d) n = 0 with the limiter ENABLED (default) — non-negativity guard: run to
#                 2× burnout and require 0 ≤ ξ ≤ limiter threshold at the end.
#                 This is the case the limiter exists for.
#
# UNIT CONVERSION (pitfall): solver rate is A·e^(−E/RT)·ξⁿ with ξ in kg/m³;
# published constants for normalized m = ξ/ξ₀ convert as A_ours = k_pub·ξ₀^(1−n)
# (see exact/kinetics.jl header). Applied below for sub-case c.

function _vrx02_build(;
    n::Float64,
    A::Float64,
    E::Float64,
    ξ0::Float64,
    T_iso::Float64,
    t_end::Float64,
    output_times::Vector{Float64},
    limiter_enabled::Bool,
)
    function build(n_cells)
        compA = SolidComponent(:A, ρ = ξ0, c = 1000.0, k = 0.5)
        compB = SolidComponent(:B, ρ = ξ0, c = 1000.0, k = 0.5)
        rxn = Reaction(:decay, 1 => 2; A = A, E = E, h = 0.0, n = n)
        mat = Material(
            name = :VRX02,
            components = (compA, compB),
            reactions = (rxn,),
            depletion_limiter = DepletionLimiter(enabled = limiter_enabled),
        )
        mesh = create_uniform_mesh(
            1e-3,
            n_cells,
            Val(2);
            T_initial = T_iso,
            ξ_initial = (ξ0, 0.0),
        )
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(T_iso), bottom = DirichletThermalBC(T_iso)),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_iso,
            ξ_initial = [z -> ξ0, z -> 0.0],
            tspan = (0.0, t_end),
            output_times = output_times,
        )
    end
    return build
end

# --- a/b/c: exact-decay sub-cases -------------------------------------------

function _vrx02_exact_case(
    label,
    title,
    refs;
    n,
    k_conc,
    E = 0.0,
    A = k_conc,
    ξ0,
    T_iso = 400.0,
    ts,
    tol_xi,
    comment = "",
    overlay_series = () -> CurveSeries[],  # report-only third-party series
)
    build = _vrx02_build(
        n = n,
        A = A,
        E = E,
        ξ0 = ξ0,
        T_iso = T_iso,
        t_end = ts[end],
        output_times = vcat(0.0, ts),
        limiter_enabled = false,
    )

    function qois(sol)
        out = QoIResult[]
        err = maximum(
            abs(sol.ξ[1, 1, j+1] - nth_order_xi(t; ξ0 = ξ0, k = k_conc, n = n)) / ξ0 for
            (j, t) in enumerate(ts)
        )
        push!(
            out,
            QoIResult(
                name = "max relative ξ_A error (n=$(n))",
                value = err,
                error = err,
                tol = tol_xi,
                comment = comment,
            ),
        )
        mass_err = maximum(abs(sum(sol.ξ[:, 1, j]) - ξ0) / ξ0 for j = 1:length(sol.t))
        push!(
            out,
            QoIResult(
                name = "total-mass conservation",
                value = mass_err,
                error = mass_err,
                tol = 1e-13,  # observed ≤ 1.1e-15 across sub-cases (2026-07-06)
                comment = "A→B closed system",
            ),
        )
        return out
    end

    return VerificationCase(
        id = "V-RX-02$(label)",
        title = title,
        tags = [:kinetics, :reaction_order, :isothermal],
        references = refs,
        build = build,
        solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
        qois = qois,
        ci = CIConfig(n_cells = 1, budget_s = 10.0),
        report = (
            qoi_curves = (:xi_histories,),
            heating = :isothermal,
            curves = sol -> [
                ReportCurve(
                    name = :xi_histories,
                    title = "V-RX-02$(label): $(title)",
                    xlabel = "t (s)",
                    ylabel = "ξ_A (kg/m³)",
                    series = vcat(
                        [
                            CurveSeries(
                                label = "ξ_A (n = $(n))",
                                x = sol.t,
                                y = sol.ξ[1, 1, :],
                                exact = t ->
                                    nth_order_xi(t; ξ0 = ξ0, k = k_conc, n = n),
                            ),
                        ],
                        overlay_series(),
                    ),
                ),
            ],
        ),
    )
end

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        ξ0, k = 100.0, 5.0  # burnout t_b = 20 s
        _vrx02_exact_case(
            "a",
            "zero-order decay (exact burnout time)",
            ["Khawam & Flanagan (2006), J. Phys. Chem. B 110:17315"];
            n = 0.0,
            k_conc = k,
            ξ0 = ξ0,
            ts = [2.0, 6.0, 10.0, 14.0, 18.0],  # ≤ 0.9·t_b (see header)
            tol_xi = 1e-12,  # observed 2.8e-16 (2026-07-06); n=0 is exactly linear
            comment = "linear decay; compared only to 0.9·t_b, limiter off",
        )
    end,
)

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        ξ0, k = 100.0, 1e-3  # timescale 1/(ξ₀k) = 10 s
        _vrx02_exact_case(
            "b",
            "second-order decay (ThermaKin C+C analog)",
            [
                "Khawam & Flanagan (2006)",
                "Stoliarov & Lyon, FSS 9:1141 (2008), second-order verification case",
            ];
            n = 2.0,
            k_conc = k,
            ξ0 = ξ0,
            ts = [5.0, 10.0, 20.0, 50.0, 100.0],
            tol_xi = 5e-8,  # 3× observed 1.7e-8 (abstol=1e-12, 2026-07-06)
            comment = "ξ = ξ₀/(1+ξ₀kt); k in m³/(kg·s)",
        )
    end,
)

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        # FDS part_baking_soda parameters at T = 450 K, converted to
        # concentration units: k_pub = A_pub·exp(−E/RT); A_ours = A_pub·ξ0^(1/3).
        ξ0 = 1000.0
        A_pub, E, T_iso = 3.4e11, 103.0e3, 450.0
        A_ours = A_pub * ξ0^(1 / 3)
        R = 8.314462618
        k_pub = A_pub * exp(-E / (R * T_iso))          # normalized-units rate
        k_conc = A_ours * exp(-E / (R * T_iso))        # concentration-units rate
        tb = nth_order_burnout_time(ξ0 = ξ0, k = k_conc, n = 2 / 3)
        ts = collect(range(0.15, 0.9, length = 6)) .* tb
        _vrx02_exact_case(
            "c",
            "contracting-sphere R3 decay (FDS NaHCO₃, 450 K)",
            [
                "Khawam & Flanagan (2006), R3 model",
                "FDS Verification Guide `part_baking_soda` (A=3.4e11 1/s, E=103 kJ/mol)",
            ];
            n = 2 / 3,
            k_conc = k_conc,
            E = E,
            A = A_ours,
            ξ0 = ξ0,
            T_iso = T_iso,
            ts = ts,
            tol_xi = 1e-9,  # 3× observed 2.8e-10 (2026-07-06)
            comment = "m^(1/3) linear in t; k_pub=$(round(k_pub, sigdigits=4)) 1/s at 450 K; " *
                      "exercises live Arrhenius (E ≠ 0) at pinned T",
            # FDS 6.11.0 overlay (reference/fds/custom/): the span-matched
            # part_baking_soda_450K re-run (T_END = 2.4 s ≈ our window). Its
            # spherical variant follows the same R3 master curve, but FDS's
            # input A = 1.33e13 is a surface-regression convention — its
            # *effective* normalized rate is 3·k_first = 3·3.4e11·e^(−E/RT)
            # at T = 450.15 K (TMPA = 177 °C; measured 1.13983 s⁻¹ vs this
            # formula's 1.13997). Overlay in master-curve time:
            # t_ours = t_FDS · k_eff/k_ours, ξ = ξ0·RHO_S/2200.
            overlay_series = () -> begin
                k_eff = 3 * A_pub * exp(-E / (R * 450.15))
                t_scale = k_eff / k_pub
                t_fds, ρ_fds = reference_columns(
                    "fds/custom/part_baking_soda_450K_matched_devc.csv",
                    "Time",
                    "RHO_S spherical",
                )
                # devc values average over the output interval → −Δt_out/2;
                # trim to our comparison window (≤ 0.9·t_b): past burnout the
                # R3 form's fractional power is undefined
                t_fds = t_fds .- (t_fds[2] - t_fds[1]) / 2
                keep = findall(t -> 0 < t_scale * t ≤ ts[end], t_fds)
                [
                    CurveSeries(
                        label = "FDS 6.11.0 (rate-rescaled)",
                        crosscode = true,
                        x = t_scale .* t_fds[keep],
                        y = ξ0 .* ρ_fds[keep] ./ 2200.0,
                        exact = t -> nth_order_xi(t; ξ0 = ξ0, k = k_conc, n = 2 / 3),
                    ),
                ]
            end,
        )
    end,
)

# --- d: n = 0 with the limiter ON — non-negativity guard ---------------------

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        ξ0, k = 100.0, 5.0
        tb = nth_order_burnout_time(ξ0 = ξ0, k = k, n = 0.0)
        t_end = 2 * tb
        build = _vrx02_build(
            n = 0.0,
            A = k,
            E = 0.0,
            ξ0 = ξ0,
            T_iso = 400.0,
            t_end = t_end,
            output_times = collect(range(0.0, t_end, length = 21)),
            limiter_enabled = true,
        )

        function qois(sol)
            ξ_min = minimum(sol.ξ[1, 1, :])
            ξ_end = sol.ξ[1, 1, end]
            return [
                QoIResult(
                    name = "non-negativity through 2× burnout",
                    value = ξ_min,
                    exact = 0.0,
                    error = max(0.0, -ξ_min),
                    tol = 1e-9,  # observed 8.2e-13 undershoot (2026-07-06)
                    comment = "the depletion limiter's job for n = 0 kinetics",
                ),
                QoIResult(
                    name = "residual mass within limiter threshold",
                    value = ξ_end,
                    exact = 0.0,
                    error = abs(ξ_end),
                    tol = 1.0,  # limiter threshold [kg/m³]: residual must sit below it
                    comment = "tanh(ξ/th) ≈ ξ/th near 0 ⇒ exponential tail; observed ~8e-13",
                ),
            ]
        end

        VerificationCase(
            id = "V-RX-02d",
            title = "zero-order decay, depletion limiter ON (guard)",
            tags = [:kinetics, :depletion_limiter, :isothermal],
            references = ["src/Materials/reactions.jl DepletionLimiter"],
            build = build,
            solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
            qois = qois,
            ci = CIConfig(n_cells = 1, budget_s = 10.0),
            report = (
                qoi_curves = (:xi_histories,),
                note = "limiter roll-off characterization",
                curves = sol -> [
                    ReportCurve(
                        name = :xi_histories,
                        title = "V-RX-02d: zero-order decay through burnout, limiter ON",
                        xlabel = "t (s)",
                        ylabel = "ξ_A (kg/m³)",
                        series = [
                            CurveSeries(
                                label = "ξ_A (limiter roll-off near 0)",
                                x = sol.t,
                                y = sol.ξ[1, 1, :],
                                # unlimited kinetics reference: clamped linear decay
                                exact = t -> max(0.0, ξ0 - k * t),
                            ),
                        ],
                    ),
                ],
            ),
        )
    end,
)
