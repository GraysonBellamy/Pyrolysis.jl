# =============================================================================
# V-PC-01a–c — A sealed cell relaxing to the sorption isotherm, athermal law
# =============================================================================
#
# A slab 4 mm thick of four cells of wood, 400 kg/m³ oven-dry with a skeletal
# density of 1500 kg/m³, holding water at a moisture content of 0.10 at
# 350 K. Both faces are sealed and adiabatic. The pores start with air at
# 1 atm and no vapour. The water and its vapour exchange through the phase
# change, K (ξ_eq − ξ_v), under an athermal law: the activity is the CAST
# spruce isotherm held at 298.15 K, a function of the moisture content alone
# with no heat of sorption, the latent heat is zero and the water takes the
# vapour's heat capacity, so the step takes no heat and the temperature stays
# where it starts. Every cell is the same, so nothing moves between them.
#
# Exact answer: the vapour follows dξ_v/dt = K (ξ_eq(ξ_w0 − ξ_v) − ξ_v), with
# the pore fraction and the moisture content following the water. Its end
# state is the root of ξ_eq(ξ_w0 − x) = x, and its history the quadrature of
# that equation (exact/phase_change.jl). Members a, b and c take
# K = 1, 100 and 10⁴ per second.
#
# Quantities: the vapour at 0.5, 1, 2, 4 and 8 times 1/K and at the end,
# 20/K + 1 s, relative to the quadrature; at the end the vapour pressure over
# the equilibrium pressure, both from the independent closures (at K = 1 the
# cell is then 21 relaxation times from its start, and e to the −21 of the way
# from equilibrium); the largest departure of any temperature from 350 K; and
# the accepted steps, which do not depend on K.

# The sealed cell of V-PC-01 and V-PC-02 on `n` cells: wood at moisture `X0`
# and `T0`, its pores holding air at 1 atm and no vapour, under the activity
# law `activity` with the latent heat `latent_heat`, the exchange rate `K`
# and the water's heat capacity `c_water`.
function _vpc_sealed(
    n;
    activity,
    K,
    latent_heat,
    c_water,
    t_end,
    output_times,
    T0 = 350.0,
    X0 = 0.10,
)
    wood = SolidComponent(
        :wood;
        ρ = 400.0,
        ρ_intrinsic = 1500.0,
        c = 1500.0,
        k = 0.12,
        λ = 5e-7,
        κ = 1e-15,
    )
    water = LiquidComponent(:water; ρ = 1000.0, c = c_water, k = 0.6)
    vapour = GaseousComponent(:vapour; M = 0.018015, c = 1900.0, k = 0.025)
    air = GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026)
    material = Material(
        name = :sealed_wood,
        components = (wood, water, vapour, air),
        phase_changes = (
            PhaseChange(
                :moisture,
                :water => :vapour;
                activity,
                sorbents = (:wood,),
                latent_heat,
                exchange_rate = K,
            ),
        ),
        mixing = MaterialMixing(
            k = MixingSpec(PARALLEL),
            λ = MixingSpec(PARALLEL),
            κ = MixingSpec(PARALLEL),
        ),
    )
    ξc = (400.0, X0 * 400.0, 0.0, 0.0)
    ξ_air = background_gas_concentration(material, :air, ξc, T0)
    return PyrolysisProblem(
        mesh = uniform_mesh(0.004, n),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> ξc[1], z -> ξc[2], z -> 0.0, z -> ξ_air],
        tspan = (0.0, t_end),
        dt_initial = 1e-6 / K,
        output_times = output_times,
    )
end

# The independent cell of `problem`'s initial state, under `law`.
function _vpc_cell(problem, law; K, L_ref, c_water, T0 = 350.0)
    ξ = [f(0.0) for f in problem.ξ_initial]
    return SealedCell(
        ξ_s = ξ[1],
        ρ_s = 1500.0,
        ξ_w0 = ξ[2],
        ξ_v0 = ξ[3],
        ξ_a = ξ[4],
        T0 = T0,
        law = law,
        K = K,
        L_ref = L_ref,
        c = (1500.0, c_water, 1900.0, 1006.0),
    )
end

# The saved index of time `t`.
_vpc_index(sol, t) = findfirst(s -> isapprox(s, t; rtol = 1e-12, atol = 1e-15), sol.t)

# The CAST spruce isotherm held at 298.15 K, a function of the moisture alone.
_vpc_frozen_spruce(X, T) = exp(-cast_gibbs(CAST_TABLE3.spruce, X, 298.15) / (PC_R * 298.15))

function _vpc01_case(member, K, steps_tol, comments)
    id = "V-PC-01$(member)"
    T0 = 350.0
    times = (0.5, 1.0, 2.0, 4.0, 8.0) ./ K
    t_end = 20 / K + 1
    build(n) = _vpc_sealed(
        n;
        activity = ActivityFunction(_vpc_frozen_spruce),
        K,
        latent_heat = 0.0,
        c_water = 1900.0,
        t_end,
        output_times = [0.0; collect(times); t_end],
    )
    cell = _vpc_cell(
        build(4),
        (:frozen, CAST_TABLE3.spruce, 298.15);
        K,
        L_ref = 0.0,
        c_water = 1900.0,
    )
    exact_vapour = map(t -> pc_athermal_vapour(cell, t), times)
    exact_end = pc_athermal_vapour(cell, t_end)
    function qois(sol)
        n = size(sol.T, 1)
        vapour = maximum(
            abs(sol.ξ[3, i, _vpc_index(sol, t)] / x - 1) for
            (t, x) in zip(times, exact_vapour) for i = 1:n
        )
        k = _vpc_index(sol, t_end)
        at_end = maximum(abs(sol.ξ[3, i, k] / exact_end - 1) for i = 1:n)
        pressure = maximum(
            abs(
                pc_vapour_pressure(cell, sol.ξ[2, i, k], sol.ξ[3, i, k], sol.T[i, k]) /
                pc_equilibrium_pressure(cell, sol.ξ[2, i, k], sol.T[i, k]) - 1,
            ) for i = 1:n
        )
        drift = maximum(abs, sol.T .- T0)
        steps = sol.statistics.accepted
        return [
            QoIResult(
                name = "vapour at five times / quadrature",
                value = vapour,
                error = vapour,
                tol = comments.vapour[1],
                comment = comments.vapour[2],
            ),
            QoIResult(
                name = "vapour at the end / quadrature",
                value = at_end,
                error = at_end,
                tol = comments.at_end[1],
                comment = comments.at_end[2],
            ),
            QoIResult(
                name = "|p_v/p_eq − 1| at the end",
                value = pressure,
                error = pressure,
                tol = comments.pressure[1],
                comment = comments.pressure[2],
            ),
            QoIResult(
                name = "temperature drift [K]",
                value = drift,
                error = drift,
                tol = 1e-10,
                comment = "the step takes no heat: the temperature stays at its start",
            ),
            QoIResult(
                name = "accepted steps (|log₂(steps/50)|)",
                value = steps,
                exact = 50.0,
                error = abs(log2(steps / 50)),
                tol = steps_tol,
                comment = comments.steps,
            ),
        ]
    end
    return VerificationCase(
        id = id,
        title = "sealed cell relaxing to the isotherm, athermal law, K = $(K) 1/s",
        tags = [:phase_change, :sorption, :sealed, :analytical],
        references = [
            "Dietenberger, Glass and Boardman (2024), the CAST isotherm, eq. 18 and Table 3 (spruce).",
            "Wagner and Pruß (2002), the saturation pressure of water, eq. 2.5.",
            "Independent quadrature of the one equation of the uniform sealed cell; see problem statement.",
        ],
        build = build,
        solve_kwargs = (reltol = 1e-10, abstol = 1e-12),
        qois = qois,
        ci = CIConfig(n_cells = 4, budget_s = 10.0),
        report = (
            curves = sol -> [
                ReportCurve(
                    name = :vapour,
                    title = "$(id): pore vapour",
                    xlabel = "t (s)",
                    ylabel = "vapour concentration (kg/m³)",
                    series = [
                        CurveSeries(
                            label = "vapour",
                            x = sol.t,
                            y = vec(sol.ξ[3, 1, :]),
                            exact = t -> pc_athermal_vapour(cell, t),
                        ),
                    ],
                ),
            ],
        ),
    )
end

# Tolerances at three times the larger of the two energy forms' observations
# (reltol 1e-10), and at 1e-13 where the observation is round-off.
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc01_case(
        "a",
        1.0,
        0.5,
        (
            vapour = (7.2e-9, "observed 1.16e-9 (temperature) and 2.39e-9 (enthalpy)"),
            at_end = (1e-11, "observed 3.2e-12 (temperature) and 9.0e-13 (enthalpy)"),
            pressure = (
                2.2e-9,
                "observed 7.19e-10 and 7.22e-10: the relaxation left after 21/K, e^(−21) of the way",
            ),
            steps = "observed 47 (temperature) and 48 (enthalpy); within a factor √2 of 50",
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc01_case(
        "b",
        1e2,
        0.5,
        (
            vapour = (7.2e-9, "observed 2.02e-9 (temperature) and 2.39e-9 (enthalpy)"),
            at_end = (1e-13, "observed 1.1e-15 and 1.3e-15, round-off"),
            pressure = (1e-13, "observed 1.1e-15 and 5.6e-16, round-off"),
            steps = "observed 51 (temperature) and 52 (enthalpy); within a factor √2 of 50",
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc01_case(
        "c",
        1e4,
        0.5,
        (
            vapour = (7.2e-9, "observed 1.85e-9 (temperature) and 2.39e-9 (enthalpy)"),
            at_end = (1e-13, "observed 8.9e-16 and 1.3e-15, round-off"),
            pressure = (1e-13, "observed 0 and 2.2e-16, round-off"),
            steps = "observed 54 in both forms; within a factor √2 of 50",
        ),
    ),
)
