# =============================================================================
# V-PC-04a–c — A drying front of free water against the Neumann solution
# =============================================================================
#
# A slab 40 mm thick of wood, 400 kg/m³ oven-dry with a skeletal density of
# 1500 kg/m³ and a heat capacity of 1500 J/(kg K), holds free water at a
# moisture content of 0.25, its pore vapour in equilibrium with the water and
# air making up 1 atm. At t = 0 its exposed face is raised to 473.15 K; the
# back is adiabatic and sealed. The face is in direct contact with dry air at
# 1 atm, and the wood's permeability is large (1e-10 m²) with no partial
# pressure transfer, so the vapour leaves through the face at about 1 atm and
# nothing condenses ahead of the front. The vapour and the air have a heat
# capacity of 1 J/(kg K): the vapour carries no sensible heat out. The water
# leaves through the phase change of free water, a = X/(|X| + 1e-4), at the
# exchange rate K, with the latent heat of water at 298.15 K, 2441.7 kJ/kg,
# and Kirchhoff's law from the heat capacities of water (4180 J/(kg K)) and
# vapour.
#
# Exact answer: in the stiff limit the water leaves at a front held at the
# boiling temperature T_b, where a p_sat(T_b) = 1 atm, 373.136 K at this
# moisture, and takes the latent heat L(T_b) there. The front then moves as
# the Neumann solution of the two-phase Stefan problem, s(t) = 2 λ √(α_d t),
# with the dry and wet zones' conductivities from the model's mixing law and
# their heat capacities from their contents (exact/stefan.jl,
# exact/phase_change.jl). Members, each on 400 cells: a, from 323.15 K at
# K = 1000 per second; b, the same at K = 10⁴ per second; c, from the boiling
# temperature at K = 1000 per second, which leaves the wet zone at
# T_b and the second phase of the problem out. In c the pores start with
# vapour alone at 1 atm, and the air, absent but for round-off, is held only
# to the solve's tolerance on the cell's mass: at a relative tolerance of 1e-6
# it sits near −0.2 Pa in the cell at the face, close to the −0.5 Pa the
# domain test allows a gas, and in the enthalpy form a step takes it past
# that. Just behind the front, where a cell's last water passes the knee of
# the law's vanishing factor, the temperature follows the solve's tolerance
# closely only from 1e-9. Member c therefore solves at 1e-9, the others at
# 1e-6.
#
# The front is read as the depth at which half the initial water is left,
# interpolated between cell centres, at the saved times from 200 to 600 s, every
# 20 s. Quantities: the mean of its error relative to the exact front; the
# error of the slope of a fit of the front against √t through the origin; the
# largest temperature error over the cells and those times; the energy ledger's
# closure; and the accepted steps, against 1.5 times the count a separate
# implementation of the model took for the same deck on the same mesh. The front's error is set by the exchange
# rate, not by the mesh: a finite rate spreads the front over a depth that does
# not shrink with the cells, which the report's ladders in the cells and in K
# show.

const _VPC04_L = 0.04
const _VPC04_X0 = 0.25
const _VPC04_T_FACE = 473.15
const _VPC04_P = 101325.0
const _VPC04_C = (1500.0, 4180.0, 1.0, 1.0)
const _VPC04_L_REF = 2.4417e6
const _VPC04_END = 600.0
const _VPC04_TIMES = collect(0.0:20.0:_VPC04_END)
const _VPC04_FIT = filter(t -> t >= _VPC04_END / 3, _VPC04_TIMES)

function _vpc04_material(K)
    wood = SolidComponent(
        :wood;
        ρ = 400.0,
        ρ_intrinsic = 1500.0,
        c = _VPC04_C[1],
        k = 0.15,
        λ = 0.0,
        κ = 1e-10,
    )
    water = LiquidComponent(:water; ρ = 1000.0, c = _VPC04_C[2], k = 0.6)
    vapour = GaseousComponent(:vapour; M = 0.018015, c = _VPC04_C[3], k = 1e-6)
    air = GaseousComponent(:air; M = 0.02897, c = _VPC04_C[4], k = 1e-6)
    return Material(
        name = :drying_wood,
        components = (wood, water, vapour, air),
        phase_changes = (
            PhaseChange(
                :moisture,
                :water => :vapour;
                activity = FreeLiquid(),
                sorbents = (:wood,),
                exchange_rate = K,
            ),
        ),
        mixing = MaterialMixing(
            k = MixingSpec(PARALLEL),
            λ = MixingSpec(PARALLEL),
            κ = MixingSpec(PARALLEL),
        ),
    )
end

# The wet start at `T_i`: wood, water, the vapour in equilibrium and air
# making up 1 atm (none at the boiling temperature).
function _vpc04_start(material, T_i)
    ξc = (400.0, _VPC04_X0 * 400.0, 0.0, 0.0)
    ξ_v = equilibrium_vapour_concentration(material, ξc, T_i)
    p_v = equilibrium_vapour_pressure(material, ξc, T_i)
    ξ_air =
        background_gas_concentration(material, :air, ξc, T_i; P = max(_VPC04_P - p_v, 0.0))
    return (ξc[1], ξc[2], ξ_v, ξ_air)
end

function _vpc04_problem(n; K, T_i)
    material = _vpc04_material(K)
    ξ0 = _vpc04_start(material, T_i)
    reservoir = GasReservoir(
        P = _VPC04_P,
        T = _VPC04_T_FACE,
        mole_fractions = (air = 1.0, vapour = 0.0),
    )
    return PyrolysisProblem(
        mesh = uniform_mesh(_VPC04_L, n),
        material = material,
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(_VPC04_T_FACE), bottom = AdiabaticBC()),
            mass = (
                top = AmbientGasBC(reservoir; film = nothing),
                bottom = ImpermeableBC(),
            ),
        ),
        T_initial = z -> T_i,
        ξ_initial = [z -> ξ0[j] for j = 1:4],
        tspan = (0.0, _VPC04_END),
        dt_initial = 1e-6,
        output_times = _VPC04_TIMES,
    )
end

# The boiling temperature of the water at the start's moisture, from the
# independent closures.
const _VPC04_T_B = pc_boiling_temperature((:free,), _VPC04_X0, _VPC04_P)

# The Neumann problem of the deck started at `T_i`, its λ and α_d.
function _vpc04_exact(T_i)
    material = _vpc04_material(1e3)
    ξ0 = _vpc04_start(material, T_i)
    ξ_dry = (400.0, 0.0, 0.0, 0.0)
    k_dry = Pyrolysis.Properties.effective_conductivity(ξ_dry, _VPC04_T_FACE, material)
    k_wet = Pyrolysis.Properties.effective_conductivity(ξ0, T_i, material)
    p = pc_drying_front(
        k_dry = k_dry,
        k_wet = k_wet,
        C_dry = 400.0 * _VPC04_C[1],
        C_wet = 400.0 * _VPC04_C[1] + ξ0[2] * _VPC04_C[2],
        ξ_w = ξ0[2],
        L = pc_latent_heat(_VPC04_L_REF, _VPC04_C[2], _VPC04_C[3], _VPC04_T_B),
        T_s = _VPC04_T_FACE,
        T_b = _VPC04_T_B,
        T_i = T_i,
    )
    λ = stefan_lambda(p)
    return (; p, λ, α_d = α_solid(p))
end

# Depths of the cell centres below the exposed face [m].
_vpc04_depths(n) = [_VPC04_L - (i - 0.5) * _VPC04_L / n for i = 1:n]

# The depth at which half the initial water is left at saved index `k`,
# interpolated between cell centres; NaN before the first cell has dried.
function _vpc04_front(sol, k)
    n = size(sol.T, 1)
    depth = _vpc04_depths(n)
    order = sortperm(depth)
    half = 0.5 * _VPC04_X0 * 400.0
    w = sol.ξ[2, order, k]
    j = findfirst(>=(half), w)
    (j === nothing || j == 1) && return NaN
    d = depth[order]
    return d[j-1] + (d[j] - d[j-1]) * (half - w[j-1]) / (w[j] - w[j-1])
end

# The errors of the front and the temperatures against the exact solution.
function _vpc04_errors(sol, exact)
    (; p, λ, α_d) = exact
    ks = [_vpc_index(sol, t) for t in _VPC04_FIT]
    front = [_vpc04_front(sol, k) for k in ks]
    t = sol.t[ks]
    relative = front ./ (2λ .* sqrt.(α_d .* t)) .- 1
    slope = sum(front .* sqrt.(t)) / sum(t)
    depth = _vpc04_depths(size(sol.T, 1))
    temperature = maximum(
        abs(sol.T[i, k] - stefan_T(depth[i], sol.t[k], p, λ)) for k in ks for
        i in eachindex(depth)
    )
    return (;
        front = sum(relative) / length(relative),
        front_max = maximum(abs, relative),
        slope = slope / (2λ * sqrt(α_d)) - 1,
        temperature,
    )
end

function _vpc04_case(member, K, T_i, n, pins)
    id = "V-PC-04$(member)"
    exact = _vpc04_exact(T_i)
    build(n) = _vpc04_problem(n; K, T_i)
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
        e = _vpc04_errors(sol, exact)
        steps = sol.statistics.accepted
        reference = get(pins.steps[1], size(sol.T, 1), nothing)
        ledger = abs(run.audit.energy.relative)
        results = [
            QoIResult(
                name = "front position, mean error (relative)",
                value = e.front,
                error = abs(e.front),
                tol = pins.front[1],
                comment = pins.front[2],
            ),
            QoIResult(
                name = "slope of the front against √t (relative)",
                value = e.slope,
                error = abs(e.slope),
                tol = pins.slope[1],
                comment = pins.slope[2],
            ),
            QoIResult(
                name = "temperature, largest error [K]",
                value = e.temperature,
                error = e.temperature,
                tol = pins.temperature[1],
                comment = pins.temperature[2],
            ),
            QoIResult(
                name = "energy ledger closure",
                value = ledger,
                error = ledger,
                tol = pins.ledger[1],
                comment = pins.ledger[2],
            ),
        ]
        # The steps against the count for the same deck on the same mesh, where
        # one exists.
        reference === nothing || push!(
            results,
            QoIResult(
                name = "accepted steps / the reference count",
                value = steps,
                exact = reference,
                error = steps / reference,
                tol = 1.5,
                comment = pins.steps[2],
            ),
        )
        return results
    end
    start = T_i == _VPC04_T_B ? "from the boiling temperature" : "from $(T_i) K"
    report = (
        n_cells_ladder = pins.ladder,
        error_norms = sol -> begin
            e = _vpc04_errors(sol, exact)
            (front = abs(e.front), slope = abs(e.slope), temperature = e.temperature)
        end,
        curves = sol -> begin
            ks = filter(k -> sol.t[k] > 0, eachindex(sol.t))
            depth = _vpc04_depths(size(sol.T, 1))
            [
                ReportCurve(
                    name = :front,
                    title = "$(id): drying front",
                    xlabel = "t (s)",
                    ylabel = "depth of the front (mm)",
                    series = [
                        CurveSeries(
                            label = "front",
                            x = sol.t[ks],
                            y = [1e3 * _vpc04_front(sol, k) for k in ks],
                            exact = t ->
                                1e3 * stefan_front(t; λ = exact.λ, αs = exact.α_d),
                        ),
                    ],
                ),
                ReportCurve(
                    name = :temperature,
                    title = "$(id): temperature at $(_VPC04_END) s",
                    xlabel = "depth (mm)",
                    ylabel = "T (K)",
                    series = [
                        CurveSeries(
                            label = "temperature",
                            x = 1e3 .* depth,
                            y = vec(sol.T[:, end]),
                            exact = x ->
                                stefan_T(x / 1e3, _VPC04_END, exact.p, exact.λ),
                        ),
                    ],
                ),
            ]
        end,
    )
    if member == "a"
        report = merge(
            report,
            (
                param_sweep = (
                    name = :K,
                    values = [1e2, 1e3, 1e4],
                    build = (n, K) -> _vpc04_problem(n; K, T_i),
                    measure = sol -> begin
                        e = _vpc04_errors(sol, exact)
                        (
                            front_mean_error_percent = 100 * e.front,
                            slope_error_percent = 100 * e.slope,
                            temperature_error_K = e.temperature,
                        )
                    end,
                ),
            ),
        )
    end
    return VerificationCase(
        id = id,
        title = "drying front of free water $(start), K = $(K) 1/s",
        tags = [:phase_change, :free_water, :moving_front, :analytical],
        references = [
            "Carslaw and Jaeger (1959), §11.2: the Neumann solution of the two-phase Stefan problem.",
            "Wagner and Pruß (2002), the saturation pressure of water, eq. 2.5.",
            "Kirchhoff's law for the latent heat from the heat capacities of water and its vapour; see problem statement.",
        ],
        build = build,
        solve_fn = solve_with_ledger,
        solve_kwargs = (reltol = pins.reltol, abstol = pins.reltol / 100, dtmax = 5.0),
        qois = qois,
        ci = CIConfig(n_cells = n, budget_s = pins.budget_s),
        report = report,
    )
end

# Tolerances at three times the larger of the two energy forms' observations,
# and at the targets where three times is above them: the front and its slope
# 1 %, the temperatures 1.5 K. The step budget is 1.5 times the count a
# separate implementation of the model took for the same deck on the same mesh
# at the same tolerance.
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc04_case(
        "a",
        1e3,
        323.15,
        400,
        (
            front = (
                4.2e-3,
                "observed −1.413e-3 in both forms; +2.17e-4 at 200 cells and −2.26e-3 at 800, toward the exchange rate's error",
            ),
            slope = (4.0e-3, "observed −1.328e-3 in both forms"),
            temperature = (
                1.5,
                "observed 1.015 K in both forms, at the front; three times is above the target, 1.5 K",
            ),
            steps = (
                Dict(200 => 4661.0, 400 => 8193.0, 800 => 12661.0),
                "a separate implementation of the model took 4661, 8193 and 12 661 at 200, 400 and 800 cells; observed 8193 (temperature) and 10 949 (enthalpy) at 400",
            ),
            ledger = (2.0e-6, "observed 2.6e-7 (temperature) and 6.6e-7 (enthalpy)"),
            reltol = 1e-6,
            ladder = [200, 400, 800],
            budget_s = 30.0,
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc04_case(
        "b",
        1e4,
        323.15,
        400,
        (
            front = (1.8e-3, "observed +6.08e-4 in both forms; +2.38e-3 at 200 cells"),
            slope = (1.9e-3, "observed +6.47e-4 in both forms"),
            temperature = (
                1.5,
                "observed 1.022 K in both forms, at the front; three times is above the target, 1.5 K",
            ),
            steps = (
                Dict(400 => 6957.0),
                "a separate implementation of the model took 6957; observed 6955 (temperature) and 9557 (enthalpy)",
            ),
            ledger = (9.2e-6, "observed 3.08e-6 (temperature) and 1.42e-6 (enthalpy)"),
            reltol = 1e-6,
            ladder = [200, 400],
            budget_s = 30.0,
        ),
    ),
)
push!(
    VERIFICATION_CASE_BUILDERS,
    () -> _vpc04_case(
        "c",
        1e3,
        _VPC04_T_B,
        400,
        (
            front = (1.5e-3, "observed −4.885e-4 in both forms"),
            slope = (1.4e-3, "observed −4.594e-4 in both forms"),
            temperature = (1.43, "observed 0.478 K in both forms"),
            steps = (
                Dict(400 => 10059.0),
                "a separate implementation of the model took 10 059; observed 10 059 (temperature) and 10 221 (enthalpy)",
            ),
            ledger = (
                1.9e-5,
                "observed 6.2e-6 (temperature) and 4.7e-6 (enthalpy); 2.7e-3 at a relative tolerance of 1e-6",
            ),
            reltol = 1e-9,
            ladder = [200, 400],
            budget_s = 40.0,
        ),
    ),
)
