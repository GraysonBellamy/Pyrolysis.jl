# =============================================================================
# V-PC-02a–c — A sealed, adiabatic cell relaxing to the isotherm, with its heat
# =============================================================================
#
# The cell of V-PC-01: a slab 4 mm thick of four cells of wood, 400 kg/m³
# oven-dry with a skeletal density of 1500 kg/m³, at a moisture content of
# 0.10 and 350 K, sealed and adiabatic, its pores holding air at 1 atm and no
# vapour. The phase change now takes its heat: the latent heat of water at
# 298.15 K, 2441.7 kJ/kg, with Kirchhoff's law from the heat capacities of
# water (4180 J/(kg K)) and of its vapour (1900 J/(kg K)), and the heat of
# sorption of the activity law. The water evaporates into the dry pores at
# the default exchange rate, 1000 per second, and the cell cools. Members:
# a, the CAST spruce isotherm; b, the CAST beech isotherm; c, free water.
#
# Exact answer: the two equations of the uniform cell, for the vapour and the
# temperature, dξ_v/dt = K (ξ_eq − ξ_v) and C dT/dt = −h K (ξ_eq − ξ_v), with
# h = L_ref + (c_v − c_w)(T − 298.15) + Q_s(X, T), integrated independently by
# classical Runge–Kutta at a step of a thousandth of 1/K
# (exact/phase_change.jl). For free water the total enthalpy, the components'
# sensible enthalpies with the latent heat of the vapour, is conserved, which
# with the water fixes the end state as well.
#
# Quantities: the temperature and the vapour at 0.5, 1, 2, 4 and 8 times 1/K;
# the closure of the solve's energy ledger; for free water, the temperature
# and the vapour at the end against the state the mass and the enthalpy fix.

function _vpc02_case(member, law, exact_law, pins)
    id = "V-PC-02$(member)"
    K = 1e3
    times = (0.5, 1.0, 2.0, 4.0, 8.0) ./ K
    t_end = 20 / K + 1
    build(n) = _vpc_sealed(
        n;
        activity = law,
        K,
        latent_heat = 2.4417e6,
        c_water = 4180.0,
        t_end,
        output_times = [0.0; collect(times); t_end],
    )
    cell = _vpc_cell(build(4), exact_law; K, L_ref = 2.4417e6, c_water = 4180.0)
    exact_vapour, exact_T = pc_thermal_history(cell, times)
    free = first(exact_law) === :free
    end_state = free ? pc_free_water_end(cell) : nothing
    last_run = Ref{Any}(nothing)
    function solve_with_ledger(problem; kwargs...)
        recorder = ConservationRecorder()
        sol = suite_solve(problem; callback = recorder, kwargs...)
        last_run[] = (; sol, audit = conservation_audit(sol, recorder))
        return sol
    end
    function qois(sol)
        run = last_run[]
        run !== nothing && run.sol === sol || error("$(id): no ledger for this solve")
        n = size(sol.T, 1)
        index = map(t -> _vpc_index(sol, t), times)
        temperature =
            maximum(abs(sol.T[i, k] - T) for (k, T) in zip(index, exact_T) for i = 1:n)
        vapour = maximum(
            abs(sol.ξ[3, i, k] / x - 1) for (k, x) in zip(index, exact_vapour) for i = 1:n
        )
        ledger = abs(run.audit.energy.relative)
        results = [
            QoIResult(
                name = "temperature at five times [K]",
                value = temperature,
                error = temperature,
                tol = pins.temperature[1],
                comment = pins.temperature[2],
            ),
            QoIResult(
                name = "vapour at five times / reference",
                value = vapour,
                error = vapour,
                tol = pins.vapour[1],
                comment = pins.vapour[2],
            ),
            QoIResult(
                name = "energy ledger closure (relative)",
                value = ledger,
                error = ledger,
                tol = pins.ledger[1],
                comment = pins.ledger[2],
            ),
        ]
        if free
            x∞, T∞ = end_state
            k = _vpc_index(sol, t_end)
            end_T = maximum(abs(sol.T[i, k] - T∞) for i = 1:n)
            end_x = maximum(abs(sol.ξ[3, i, k] / x∞ - 1) for i = 1:n)
            push!(
                results,
                QoIResult(
                    name = "temperature at the end / mass and enthalpy [K]",
                    value = sol.T[1, k],
                    exact = T∞,
                    error = end_T,
                    tol = pins.end_T[1],
                    comment = pins.end_T[2],
                ),
                QoIResult(
                    name = "vapour at the end / mass and enthalpy",
                    value = sol.ξ[3, 1, k],
                    exact = x∞,
                    error = end_x,
                    tol = pins.end_x[1],
                    comment = pins.end_x[2],
                ),
            )
        end
        return results
    end
    name =
        free ? "free water" :
        exact_law[2] === CAST_TABLE3.spruce ? "CAST spruce" : "CAST beech"
    return VerificationCase(
        id = id,
        title = "sealed adiabatic cell relaxing to the isotherm, $(name)",
        tags = [:phase_change, :sorption, :sealed, :adiabatic, :analytical],
        references = [
            "Dietenberger, Glass and Boardman (2024), the CAST isotherm, eqs. 18 and 19 and Table 3.",
            "Wagner and Pruß (2002), the saturation pressure of water, eq. 2.5.",
            "Independent Runge–Kutta integration of the two equations of the uniform sealed cell; see problem statement.",
        ],
        build = build,
        solve_fn = solve_with_ledger,
        solve_kwargs = (reltol = 1e-10, abstol = 1e-12),
        qois = qois,
        ci = CIConfig(n_cells = 4, budget_s = 10.0),
        report = (
            curves = sol -> [
                ReportCurve(
                    name = :temperature,
                    title = "$(id): temperature",
                    xlabel = "t (s)",
                    ylabel = "T (K)",
                    series = [
                        CurveSeries(
                            label = "temperature",
                            x = sol.t[1:(end-1)],
                            y = vec(sol.T[1, 1:(end-1)]),
                            exact = t -> last(pc_thermal_state(cell, t)),
                        ),
                    ],
                ),
            ],
        ),
    )
end

# Tolerances at three times the larger of the two energy forms' observations
# (reltol 1e-10).
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc02_case(
        "a",
        CAST_SPRUCE,
        (:cast, CAST_TABLE3.spruce),
        (
            temperature = (
                1.6e-9,
                "observed 5.13e-10 K (temperature) and 1.18e-10 K (enthalpy)",
            ),
            vapour = (5.4e-9, "observed 1.77e-9 (temperature) and 4.10e-10 (enthalpy)"),
            ledger = (7.6e-9, "observed 2.51e-9 (temperature) and 1.04e-9 (enthalpy)"),
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc02_case(
        "b",
        CAST_BEECH,
        (:cast, CAST_TABLE3.beech),
        (
            temperature = (
                1.7e-9,
                "observed 5.35e-10 K (temperature) and 3.38e-10 K (enthalpy)",
            ),
            vapour = (6.6e-9, "observed 1.90e-9 (temperature) and 2.17e-9 (enthalpy)"),
            ledger = (7.6e-9, "observed 2.51e-9 (temperature) and 1.09e-9 (enthalpy)"),
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc02_case(
        "c",
        FreeLiquid(),
        (:free,),
        (
            temperature = (
                1.6e-9,
                "observed 5.13e-10 K (temperature) and 1.81e-10 K (enthalpy)",
            ),
            vapour = (3.0e-9, "observed 9.73e-10 (temperature) and 4.18e-10 (enthalpy)"),
            ledger = (6.1e-9, "observed 2.03e-9 (temperature) and 9.11e-10 (enthalpy)"),
            end_T = (3.3e-11, "observed 5.9e-12 K (temperature) and 1.1e-11 K (enthalpy)"),
            end_x = (1.3e-12, "observed 2.3e-13 (temperature) and 4.2e-13 (enthalpy)"),
        ),
    ),
)
