# =============================================================================
# V-RX-02 — nth-order and contracting-geometry isothermal decay
# =============================================================================
#
# Khawam & Flanagan (2006) exact isothermal integrations (FDS verifies its
# `part_baking_soda` case against the same forms). Sub-cases:
#
#   a) n = 0    — linear decay, exact finite burnout time. Compared only up to
#                 0.9·t_b with the depletion limiter DISABLED: the power is
#                 taken on the finite base of the conversion, so at ξ = 0 the
#                 n = 0 rate stays k ξ₀ and would integrate ξ negative
#                 (`DepletionLimiter(enabled = false)`).
#   b) n = 2    — single-species second-order decay (a C+C kinetic analog,
#                 not a reproduced ThermaKin run):
#                 ξ = ξ₀/(1 + kt).
#   c) n = 2/3  — contracting-sphere R3 model with FDS's NaHCO₃ parameters
#                 A = 3.4e11 1/s, E = 103 kJ/mol at T = 450 K; the marker
#                 is m^(1/3) linear in t.
#   d) n = 0 with the limiter ENABLED (default) — non-negativity guard: run to
#                 2× burnout and require 0 ≤ x ≤ limiter threshold at the end.
#                 This is the case the limiter exists for.
#
# The rate is k·ξ₀·xⁿ on the conversion x = ξ/ξ₀, k in s⁻¹ for every order, so
# published constants for normalized m = ξ/ξ₀ are the solver's as they stand
# (see exact/kinetics.jl header).

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
        mesh = create_uniform_mesh(1e-3, n_cells)
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
    k,
    E = 0.0,
    A = k,
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
            abs(sol.ξ[1, 1, j+1] - nth_order_xi(t; ξ0 = ξ0, k = k, n = n)) / ξ0 for
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
                tol = 1e-13,  # observed ≤ 1.1e-15 across sub-cases
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
                                exact = t -> nth_order_xi(t; ξ0 = ξ0, k = k, n = n),
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
        ξ0, k = 100.0, 0.05  # burnout t_b = 20 s
        _vrx02_exact_case(
            "a",
            "zero-order decay (exact burnout time)",
            ["Khawam & Flanagan (2006), J. Phys. Chem. B 110:17315"];
            n = 0.0,
            k = k,
            ξ0 = ξ0,
            ts = [2.0, 6.0, 10.0, 14.0, 18.0],  # ≤ 0.9·t_b (see header)
            tol_xi = 1e-12,  # observed 8.3e-15; n=0 is exactly linear
            comment = "linear decay; compared only to 0.9·t_b, limiter off",
        )
    end,
)

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        ξ0, k = 100.0, 0.1  # timescale 1/k = 10 s
        _vrx02_exact_case(
            "b",
            "second-order decay (ThermaKin C+C analog)",
            [
                "Khawam & Flanagan (2006)",
                "Stoliarov & Lyon, FSS 9:1141 (2008), second-order verification case",
            ];
            n = 2.0,
            k = k,
            ξ0 = ξ0,
            ts = [5.0, 10.0, 20.0, 50.0, 100.0],
            tol_xi = 1e-11,  # 3× observed 2.9e-12 (abstol=1e-12)
            comment = "ξ = ξ₀/(1+kt); k in s⁻¹",
        )
    end,
)

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        # FDS part_baking_soda parameters at T = 450 K: k = A·exp(−E/RT).
        ξ0 = 1000.0
        A, E, T_iso = 3.4e11, 103.0e3, 450.0
        R = 8.314462618
        k = A * exp(-E / (R * T_iso))
        tb = nth_order_burnout_time(k = k, n = 2 / 3)
        ts = collect(range(0.15, 0.9, length = 6)) .* tb
        _vrx02_exact_case(
            "c",
            "contracting-sphere R3 decay (FDS NaHCO₃, 450 K)",
            [
                "Khawam & Flanagan (2006), R3 model",
                "FDS Verification Guide `part_baking_soda` (A=3.4e11 1/s, E=103 kJ/mol)",
            ];
            n = 2 / 3,
            k = k,
            E = E,
            A = A,
            ξ0 = ξ0,
            T_iso = T_iso,
            ts = ts,
            tol_xi = 1e-9,  # 3× observed 2.8e-10
            comment = "m^(1/3) linear in t; k=$(round(k, sigdigits=4)) 1/s at 450 K; " *
                      "exercises live Arrhenius (E ≠ 0) at pinned T",
            # Native-input FDS/Gpyro runs are added by report mode from
            # reference/matched/. They use the same fixed-temperature rate.
        )
    end,
)

# --- d: n = 0 with the limiter ON — non-negativity guard ---------------------

push!(
    VERIFICATION_CASE_BUILDERS,
    () -> begin
        ξ0, k = 100.0, 0.05
        tb = nth_order_burnout_time(k = k, n = 0.0)
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
                    tol = 1e-9,  # no undershoot observed
                    comment = "the depletion limiter's job for n = 0 kinetics",
                ),
                QoIResult(
                    name = "residual mass within limiter threshold",
                    value = ξ_end,
                    exact = 0.0,
                    error = abs(ξ_end),
                    tol = 0.1,  # conversion threshold of ξ₀ [kg/m³]: residual must sit below it
                    comment = "quadratic onset near 0 ⇒ algebraic tail; observed 5.4e-5",
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
                                exact = t -> max(0.0, ξ0 * (1 - k * t)),
                            ),
                        ],
                    ),
                ],
            ),
        )
    end,
)
