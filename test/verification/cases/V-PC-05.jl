# =============================================================================
# V-PC-05a–b — Sorption into a slab from a humid reservoir, linear isotherm
# =============================================================================
#
# A slab 1 mm thick of wood, 400 kg/m³ oven-dry with a skeletal density of
# 1500 kg/m³, at 298.15 K, holds water in equilibrium with its pore vapour at
# an activity of 0.50, with air making up 1 atm. Its back is sealed and the
# wood has no permeability. At t = 0 its exposed face meets, in direct
# contact with no film, a reservoir of air at 1 atm whose humidity is an
# activity of 0.55. The vapour moves into the slab and the water takes it up.
# The water follows a linear law, a = X/0.2, through the phase change at an
# exchange rate of 10⁵ per second; the law does not depend on temperature,
# the latent heat is zero and the water has its vapour's heat capacity, so the
# step takes no heat and the slab stays at 298.15 K. The law's vanishing
# factor, X/(|X| + X₀), is taken with X₀ = 1e-7, where it moves the activity by
# 1e-6 and its slope by 1e-12. Members: a, the lumped closure with a transfer
# coefficient of 5e-7 m²/s; b, the dusty-gas closure, a tortuosity of 4, a
# pore diameter of 10 µm and a vapour–air diffusivity of 2.6e-5 m²/s at
# 298.15 K and 1 atm.
#
# Exact answer: in the stiff limit the vapour's density in the pores obeys a
# diffusion equation, its flux the transfer coefficient times its gradient,
# its storage the water's through the isotherm with the vapour's own beside it.
# The slab's gain is then the series of a slab held at a fixed value at its
# face (exact/phase_change.jl, the limit of exact/slab_series.jl's convective
# slab), with the gain at the end from the two equilibrium states and the
# diffusivity the transfer coefficient over the chord of the storage. Under
# the dusty-gas closure the transfer coefficient is that of vapour diffusing
# through air that does not move, evaluated at the mean of the two states'
# pore fractions and mole fractions; the pore fraction falls by 0.6 % as the
# water's volume grows, and the coefficient with it.
#
# Quantities: the water and vapour the slab has gained, per unit area, at
# Fourier numbers of 0.02, 0.05, 0.1, 0.2 and 0.5, relative to the series;
# and the largest departure of a temperature from 298.15 K. The exchange rate
# leaves the water a relaxation time of 7e-6 (a) and 7e-5 (b) of the
# diffusion time, and the water's lag behind the vapour costs about that over
# four times the Fourier number: 1e-4 and 8e-4 of the gain at the first time.
# The report's ladder in the step of the reservoir's activity, 0.05, 0.025 and
# 0.0125, shows how the error depends on the step: under the lumped closure
# the pore fraction enters only the vapour's own storage, and the step moves
# the error by 1e-7; under the dusty-gas closure the transfer coefficient
# follows the pore fraction and the mole fractions, and the error at the last
# time falls in proportion to the step.

const _VPC05_L = 1e-3
const _VPC05_T = 298.15
const _VPC05_P = 101325.0
const _VPC05_X_S = 0.2
const _VPC05_X0 = 1e-7
const _VPC05_A0 = 0.5
const _VPC05_S = 400.0 + PC_SORBENT_FLOOR
const _VPC05_FO = (0.02, 0.05, 0.1, 0.2, 0.5)
const _VPC05_FO_END = 1.0
const _VPC05_LAW = (:linear, _VPC05_X_S, _VPC05_X0)
const _VPC05_DUSTY = (τ = 4.0, d = 1e-5, D_ref = 2.6e-5)

_vpc05_linear(X, T) = X / _VPC05_X_S

function _vpc05_material(closure, K)
    wood = SolidComponent(
        :wood;
        ρ = 400.0,
        ρ_intrinsic = 1500.0,
        c = 1500.0,
        k = 0.12,
        λ = 5e-7,
        κ = 0.0,
    )
    water = LiquidComponent(:water; ρ = 1000.0, c = 1900.0, k = 0.6)
    vapour = GaseousComponent(:vapour; M = 0.018015, c = 1900.0, k = 0.025)
    air = GaseousComponent(:air; M = 0.02897, c = 1006.0, k = 0.026)
    transport =
        closure === :lumped ? LumpedTransfer() :
        DustyGas(
            tortuosity = _VPC05_DUSTY.τ,
            D_ref = _VPC05_DUSTY.D_ref,
            T_ref = _VPC05_T,
            P_ref = _VPC05_P,
            pore_diameter = _VPC05_DUSTY.d,
        )
    return Material(
        name = :sorbing_wood,
        components = (wood, water, vapour, air),
        phase_changes = (
            PhaseChange(
                :moisture,
                :water => :vapour;
                activity = ActivityFunction(_vpc05_linear; X0 = _VPC05_X0),
                sorbents = (:wood,),
                latent_heat = 0.0,
                exchange_rate = K,
            ),
        ),
        mixing = MaterialMixing(
            k = MixingSpec(PARALLEL),
            λ = MixingSpec(PARALLEL),
            κ = MixingSpec(PARALLEL),
        ),
        gas_transport = transport,
    )
end

# The uniform state in equilibrium at the activity `a`: wood, water, vapour
# and air making up 1 atm.
function _vpc05_state(material, a)
    X = pc_linear_moisture(a, _VPC05_X_S, _VPC05_X0)
    ξc = (400.0, X * _VPC05_S, 0.0, 0.0)
    ξ_v = equilibrium_vapour_concentration(material, ξc, _VPC05_T)
    p_v = equilibrium_vapour_pressure(material, ξc, _VPC05_T)
    ξ_air = background_gas_concentration(material, :air, ξc, _VPC05_T; P = _VPC05_P - p_v)
    return (ξc[1], ξc[2], ξ_v, ξ_air)
end

# The exact problem of the member: the gain at the end per unit area, the
# diffusivity and the times of the five Fourier numbers and the end.
function _vpc05_exact(closure, step)
    a∞ = _VPC05_A0 + step
    ρ_sat = 0.018015 * pc_saturation(_VPC05_T) / (PC_R * _VPC05_T)
    pore(a) = begin
        X = pc_linear_moisture(a, _VPC05_X_S, _VPC05_X0)
        1 - 400.0 / 1500.0 - X * _VPC05_S / 1000.0
    end
    W(a) = pc_linear_moisture(a, _VPC05_X_S, _VPC05_X0) * _VPC05_S + pore(a) * a * ρ_sat
    gain = (W(a∞) - W(_VPC05_A0)) * _VPC05_L
    Λ = if closure === :lumped
        5e-7
    else
        x(a) = a * pc_saturation(_VPC05_T) / _VPC05_P
        pc_dusty_gas_transfer(
            φ = (pore(_VPC05_A0) + pore(a∞)) / 2,
            τ = _VPC05_DUSTY.τ,
            d = _VPC05_DUSTY.d,
            D_ref = _VPC05_DUSTY.D_ref,
            T_ref = _VPC05_T,
            P_ref = _VPC05_P,
            T = _VPC05_T,
            P = _VPC05_P,
            M = 0.018015,
            x_a = 1 - (x(_VPC05_A0) + x(a∞)) / 2,
        )
    end
    D = Λ * step * ρ_sat * _VPC05_L / gain
    time(Fo) = Fo * _VPC05_L^2 / D
    return (; gain, D, Λ, times = map(time, _VPC05_FO), t_end = time(_VPC05_FO_END))
end

function _vpc05_problem(n; closure, K, step)
    material = _vpc05_material(closure, K)
    ξ0 = _vpc05_state(material, _VPC05_A0)
    x_v = (_VPC05_A0 + step) * WATER_SATURATION(_VPC05_T) / _VPC05_P
    reservoir = GasReservoir(
        P = _VPC05_P,
        T = _VPC05_T,
        mole_fractions = (vapour = x_v, air = 1 - x_v),
    )
    exact = _vpc05_exact(closure, step)
    times = sort(unique([collect(range(0.0, exact.t_end, 101)); collect(exact.times)]))
    return PyrolysisProblem(
        mesh = uniform_mesh(_VPC05_L, n),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = AdiabaticBC(), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(reservoir; film = nothing),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> _VPC05_T,
        ξ_initial = [z -> ξ0[j] for j = 1:4],
        tspan = (0.0, exact.t_end),
        dt_initial = 1e-3 / K,
        output_times = times,
    )
end

# The water and vapour gained per unit area at every saved time [kg/m²].
function _vpc05_gain(sol)
    n = size(sol.T, 1)
    V = [sol.m[1, i, 1] / sol.ξ[1, i, 1] for i = 1:n]
    area = sum(V) / _VPC05_L
    held(k) = sum(sol.m[2, i, k] + sol.m[3, i, k] for i = 1:n)
    return [(held(k) - held(1)) / area for k in eachindex(sol.t)]
end

# The relative errors of the gain at the five times, and the largest.
function _vpc05_errors(sol, exact)
    gain = _vpc05_gain(sol)
    return [
        gain[_vpc_index(sol, t)] / (exact.gain * pc_slab_uptake(Fo)) - 1 for
        (t, Fo) in zip(exact.times, _VPC05_FO)
    ]
end
_vpc05_error(sol, exact) = maximum(abs, _vpc05_errors(sol, exact))

function _vpc05_case(member, closure, n, pins)
    id = "V-PC-05$(member)"
    K = 1e5
    step = 0.05
    exact = _vpc05_exact(closure, step)
    build(n) = _vpc05_problem(n; closure, K, step)
    function qois(sol)
        gain = _vpc05_error(sol, exact)
        drift = maximum(abs, sol.T .- _VPC05_T)
        return [
            QoIResult(
                name = "gain at five times / series (relative)",
                value = gain,
                error = gain,
                tol = pins.gain[1],
                comment = pins.gain[2],
            ),
            QoIResult(
                name = "temperature drift [K]",
                value = drift,
                error = drift,
                tol = pins.drift[1],
                comment = pins.drift[2],
            ),
        ]
    end
    name = closure === :lumped ? "lumped transfer" : "dusty gas"
    return VerificationCase(
        id = id,
        title = "sorption into a slab from a humid reservoir, linear isotherm, $(name)",
        tags = [:phase_change, :sorption, :gas_transport, :analytical],
        references = [
            "Crank (1975), The Mathematics of Diffusion, eq. 4.18: uptake by a plane sheet.",
            "Mason and Malinauskas (1983), the dusty-gas model: diffusion through a stagnant gas.",
            "Wagner and Pruß (2002), the saturation pressure of water, eq. 2.5.",
        ],
        build = build,
        solve_kwargs = pins.tolerances,
        qois = qois,
        ci = CIConfig(n_cells = n, budget_s = pins.budget_s),
        report = (
            n_cells_ladder = pins.ladder,
            error_norms = sol -> (gain = _vpc05_error(sol, exact),),
            curves = sol -> [
                ReportCurve(
                    name = :gain,
                    title = "$(id): water and vapour gained",
                    xlabel = "t (s)",
                    ylabel = "gain (g/m²)",
                    series = [
                        CurveSeries(
                            label = "gain",
                            x = sol.t,
                            y = 1e3 .* _vpc05_gain(sol),
                            exact = t ->
                                1e3 * exact.gain * pc_slab_uptake(exact.D * t / _VPC05_L^2),
                        ),
                    ],
                ),
            ],
            param_sweep = (
                name = :step,
                values = [0.05, 0.025, 0.0125],
                build = (n, s) -> _vpc05_problem(n; closure, K, step = s),
                measure = (sol, s) -> begin
                    e = _vpc05_errors(sol, _vpc05_exact(closure, s))
                    (
                        largest_error_percent = 100 * maximum(abs, e),
                        error_at_Fo_0_5_percent = 100 * e[end],
                    )
                end,
            ),
        ),
    )
end

# Tolerances at three times the larger of the two energy forms' observations.
# The gain is the water's, not a trace: it is the same to 1e-8 from a relative
# tolerance of 1e-6 to 1e-12, and the case solves at the package's default.
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc05_case(
        "a",
        :lumped,
        50,
        (
            gain = (
                4.0e-3,
                "observed 1.343e-3 in both forms, at the first time; the mesh's: 4.0e-4 at 100 cells, 1.7e-4 at 200",
            ),
            drift = (
                1e-10,
                "the step takes no heat; observed 5.7e-14 K (temperature) and 0 (enthalpy)",
            ),
            tolerances = (reltol = 1e-6, abstol = 1e-8),
            ladder = [25, 50, 100, 200],
            budget_s = 5.0,
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc05_case(
        "b",
        :dusty,
        50,
        (
            gain = (
                6.9e-3,
                "observed 2.309e-3 in both forms, at the first time; with the coefficient's trace limit, x_a = 1, it would be 7.6e-3",
            ),
            drift = (
                1e-10,
                "the step takes no heat; observed 5.7e-14 K (temperature) and 0 (enthalpy)",
            ),
            tolerances = (reltol = 1e-6, abstol = 1e-8),
            ladder = [25, 50, 100, 200],
            budget_s = 5.0,
        ),
    ),
)
