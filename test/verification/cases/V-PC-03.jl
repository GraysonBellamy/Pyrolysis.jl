# =============================================================================
# V-PC-03a–b — A sealed slab heated and cooled through one face: the latent
# heat over a closed cycle
# =============================================================================
#
# A slab 2 mm thick of ten cells of wood, 400 kg/m³ oven-dry with a skeletal
# density of 1500 kg/m³, holding water at a moisture content of 0.10 at
# 300 K, its pore vapour in equilibrium with the water and air making up
# 1 atm. Both faces are sealed; the back is adiabatic and the exposed face
# convects, at 50 W/(m² K), to an ambient that rises smoothly from 300 to
# 400 K and back over 300 s (300 + 100 sin² of π t/300) and is then held at
# 300 K for 900 s, until the slab is at rest. The water evaporates as the
# slab heats, the vapour moves inside it, and the water condenses again as it
# cools. The phase change takes the latent heat of water at 298.15 K,
# 2441.7 kJ/kg, with Kirchhoff's law from the heat capacities of water
# (4180 J/(kg K)) and of its vapour (1900 J/(kg K)), and the activity law's
# heat of sorption. Members: a, free water; b, the CAST spruce isotherm.
#
# Exact answer: with the total enthalpy H, the sensible enthalpies of all the
# components with the latent heat at 298.15 K of the vapour, the heat through
# the face equals the change of H at every time when the latent heat is the
# components' own enthalpy difference: for free water exactly, and the slab
# returns to its start. Under an isotherm the heat of sorption is not the
# derivative of a state function, so over the cycle the heat through the face
# and the change of H differ by the heat of sorption exchanged, its time
# integral over the cells of Q_s times the evaporation rate. A latent heat
# that ignored the heat capacities, the latent heat at 298.15 K at every
# temperature, would leave 8e-4 of the peak rise of H in that difference over
# this cycle of free water.
#
# The heat through the face is the time integral of the flux the convective
# condition applies, by Simpson's rule on the saved times (every 0.25 s); the
# total enthalpy is evaluated from the saved masses and temperatures with
# the components' heat capacities, independently of the solver's ledger.
#
# Quantities, a: the heat through the face less the change of H, over the
# peak rise of H, at 75, 150, 225, 300 and 600 s and at the end, 1200 s; the
# largest departure of a temperature from 300 K and of the vapour from its
# start at the end. b: the same difference less the heat of sorption
# exchanged, over the peak rise of H, at the same times; and at the end the
# difference itself over the heat of sorption exchanged in either direction,
# which the excess heat capacity that the model leaves out bounds at 5 %.

const _VPC03_C = (1500.0, 4180.0, 1900.0, 1006.0)
const _VPC03_L_REF = 2.4417e6
const _VPC03_DT = 0.25
const _VPC03_CYCLE = 300.0
const _VPC03_END = 1200.0
const _VPC03_TIMES = (75.0, 150.0, 225.0, 300.0, 600.0, _VPC03_END)

_vpc03_ambient(t) = t < _VPC03_CYCLE ? 300.0 + 100.0 * sinpi(t / _VPC03_CYCLE)^2 : 300.0

function _vpc03_problem(n, activity)
    wood = SolidComponent(
        :wood;
        ρ = 400.0,
        ρ_intrinsic = 1500.0,
        c = _VPC03_C[1],
        k = 0.12,
        λ = 5e-7,
        κ = 1e-15,
    )
    water = LiquidComponent(:water; ρ = 1000.0, c = _VPC03_C[2], k = 0.6)
    vapour = GaseousComponent(:vapour; M = 0.018015, c = _VPC03_C[3], k = 0.025)
    air = GaseousComponent(:air; M = 0.02897, c = _VPC03_C[4], k = 0.026)
    material = Material(
        name = :cycled_wood,
        components = (wood, water, vapour, air),
        phase_changes = (
            PhaseChange(:moisture, :water => :vapour; activity, sorbents = (:wood,)),
        ),
        mixing = MaterialMixing(
            k = MixingSpec(PARALLEL),
            λ = MixingSpec(PARALLEL),
            κ = MixingSpec(PARALLEL),
        ),
    )
    T0 = 300.0
    ξc = (400.0, 40.0, 0.0, 0.0)
    ξ_v = equilibrium_vapour_concentration(material, ξc, T0)
    p_v = equilibrium_vapour_pressure(material, ξc, T0)
    ξ_air = background_gas_concentration(material, :air, ξc, T0; P = 101325.0 - p_v)
    return PyrolysisProblem(
        mesh = uniform_mesh(0.002, n),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (
                top = ConvectiveBC(h = 50.0, T_inf = _vpc03_ambient),
                bottom = AdiabaticBC(),
            ),
            mass = (top = ImpermeableBC(), bottom = ImpermeableBC()),
        ),
        T_initial = z -> T0,
        ξ_initial = [z -> ξc[1], z -> ξc[2], z -> ξ_v, z -> ξ_air],
        tspan = (0.0, _VPC03_END),
        output_times = collect(0.0:_VPC03_DT:_VPC03_END),
    )
end

# The integral of the samples `y` at spacing `h` from the first to the `k`-th,
# `k − 1` even, by Simpson's rule.
function _vpc03_simpson(y, h, k)
    iseven(k - 1) ||
        throw(ArgumentError("Simpson's rule needs an even number of intervals"))
    k == 1 && return 0.0
    return h / 3 * (y[1] + y[k] + 4 * sum(y[2:2:(k-1)]) + 2 * sum(y[3:2:(k-2)]; init = 0.0))
end

# The saved series of the cycle: per saved time the heat flux through the
# face [W] and the total enthalpy [J], and per cell and saved time the heat of
# sorption times the evaporation rate [W], from the independent closures.
function _vpc03_series(sol, law)
    nt = length(sol.t)
    maximum(abs.(diff(sol.t) .- _VPC03_DT)) <= 1e-9 ||
        error("V-PC-03: the saved times are not every $(_VPC03_DT) s")
    n = size(sol.T, 1)
    V = [sol.m[1, i, 1] / sol.ξ[1, i, 1] for i = 1:n]
    area = sum(V) / 0.002
    q = sol.extras.surface_heat_flux_net .* area
    H = map(1:nt) do k
        sum(
            sum(sol.m[j, i, k] * _VPC03_C[j] for j = 1:4) * (sol.T[i, k] - PC_T_REF) +
            _VPC03_L_REF * sol.m[3, i, k] for i = 1:n
        )
    end
    cell = SealedCell(
        ξ_s = 400.0,
        ρ_s = 1500.0,
        ξ_w0 = 40.0,
        ξ_a = 0.0,
        T0 = 300.0,
        law = law,
        K = 1e3,
        c = _VPC03_C,
    )
    sorption = map(1:nt) do k
        sum(1:n) do i
            ξ_w, ξ_v, T = sol.ξ[2, i, k], sol.ξ[3, i, k], sol.T[i, k]
            ω = cell.K * (pc_equilibrium_vapour(cell, ξ_w, T) - ξ_v)
            pc_sorption_heat(law, pc_moisture(cell, ξ_w), T, cell.M) * ω * V[i]
        end
    end
    return (; q, H, sorption, absolute = abs.(sorption), area)
end

function _vpc03_case(member, activity, law, pins)
    id = "V-PC-03$(member)"
    free = first(law) === :free
    build(n) = _vpc03_problem(n, activity)
    index(sol, t) = _vpc_index(sol, t)
    function qois(sol)
        s = _vpc03_series(sol, law)
        peak = maximum(s.H .- s.H[1])
        residual(k) = _vpc03_simpson(s.q, _VPC03_DT, k) - (s.H[k] - s.H[1])
        last = index(sol, _VPC03_END)
        if free
            cycle = maximum(abs(residual(index(sol, t))) for t in _VPC03_TIMES) / peak
            temperature = maximum(abs, sol.T[:, last] .- 300.0)
            vapour = abs(sum(sol.m[3, :, last]) / sum(sol.m[3, :, 1]) - 1)
            return [
                QoIResult(
                    name = "|∫q dt − ΔH| / peak rise of H, five times and the end",
                    value = cycle,
                    error = cycle,
                    tol = pins.cycle[1],
                    comment = pins.cycle[2],
                ),
                QoIResult(
                    name = "temperature back at the start [K]",
                    value = temperature,
                    error = temperature,
                    tol = pins.temperature[1],
                    comment = pins.temperature[2],
                ),
                QoIResult(
                    name = "vapour back at the start (relative)",
                    value = vapour,
                    error = vapour,
                    tol = pins.vapour[1],
                    comment = pins.vapour[2],
                ),
            ]
        end
        sorbed(k) = _vpc03_simpson(s.sorption, _VPC03_DT, k)
        identity =
            maximum(
                abs(residual(index(sol, t)) - sorbed(index(sol, t))) for t in _VPC03_TIMES
            ) / peak
        exchanged = _vpc03_simpson(s.absolute, _VPC03_DT, last)
        excess = abs(residual(last)) / exchanged
        return [
            QoIResult(
                name = "|∫q dt − ΔH − ∫Q_s ω V dt| / peak rise of H, five times and the end",
                value = identity,
                error = identity,
                tol = pins.identity[1],
                comment = pins.identity[2],
            ),
            QoIResult(
                name = "cycle residual / heat of sorption exchanged",
                value = excess,
                error = excess,
                tol = 0.05,
                comment = pins.excess,
            ),
        ]
    end
    return VerificationCase(
        id = id,
        title = "sealed slab heated and cooled through one face, $(free ? "free water" : "CAST spruce")",
        tags = [:phase_change, :sorption, :energy_conservation, :sealed],
        references = [
            "Total enthalpy of the components with the latent heat at 298.15 K: Kirchhoff's law; see problem statement.",
            "Dietenberger, Glass and Boardman (2024), the CAST isotherm, eqs. 18 and 19 and Table 3.",
            "Wagner and Pruß (2002), the saturation pressure of water, eq. 2.5.",
        ],
        build = build,
        solve_kwargs = (reltol = 1e-10, abstol = 1e-12),
        qois = qois,
        ci = CIConfig(n_cells = 10, budget_s = 60.0),
        report = (
            curves = sol -> begin
                s = _vpc03_series(sol, law)
                k = length(sol.t)
                heat = [_vpc03_simpson(s.q, _VPC03_DT, j) for j = 1:2:k]
                [
                    ReportCurve(
                        name = :ledger,
                        title = "$(id): heat through the face and change of the total enthalpy",
                        xlabel = "t (s)",
                        ylabel = "energy (J/m²)",
                        series = [
                            CurveSeries(
                                label = "∫q dt",
                                x = sol.t[1:2:k],
                                y = heat ./ s.area,
                            ),
                            CurveSeries(
                                label = "ΔH",
                                x = sol.t[1:2:k],
                                y = (s.H[1:2:k] .- s.H[1]) ./ s.area,
                            ),
                        ],
                    ),
                ]
            end,
        ),
    )
end

# Tolerances at three times the larger of the two energy forms' observations
# (reltol 1e-10). A latent heat that ignored the heat capacities would leave
# 8.1e-4 of the peak rise of H in the cycle residual of member a: a value above
# 1e-4 is a defect in the heat, not a tolerance to widen. At the end the cycle has
# moved water toward the back, and with no capillary flow it diffuses back as
# vapour over thousands of seconds: the evaporation and condensation that
# drives hold the cells about 2e-5 K from 300 K (2.31e-5 K at 3000 s).
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc03_case(
        "a",
        FreeLiquid(),
        (:free,),
        (
            cycle = (4.9e-10, "observed 1.64e-10 (temperature) and 5.96e-11 (enthalpy)"),
            temperature = (
                7.1e-5,
                "observed 2.34e-5 K in both forms: the water the cycle moved still flows back",
            ),
            vapour = (6.7e-6, "observed 2.21e-6 in both forms, as the temperature"),
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc03_case(
        "b",
        CAST_SPRUCE,
        (:cast, CAST_TABLE3.spruce),
        (
            identity = (3.9e-8, "observed 1.27e-8 (temperature) and 2.37e-9 (enthalpy)"),
            excess = "the limit of the excess heat capacity; observed 0.0171 in both forms",
        ),
    ),
)
