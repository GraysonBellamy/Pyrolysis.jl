# =============================================================================
# V-RAD-06 — the radiative energy ledger through burnout
# =============================================================================
#
# A 4 mm GIDAZE+ PMMA sample (V-RAD-05's deck, 40 cells, moving mesh) burns
# out completely under a collimated heater of 100 kW/m², facing surroundings
# at 293 K through its exposed face and through an open back (a `RadiativeBC`
# on `:bottom`), so radiation leaves through both faces. The heater is
# absorbed in depth and the sample pyrolyses through its volume, so its cells
# empty together, in three deactivation invocations between 86 and 123 s,
# and the exposed face, with the radiation rows' boundary, moves down through
# the sample until none is left. The pore gas runs under a pressure tolerance
# of 1e-3: without it the gas the emptying cells release closes the mass
# ledger to 1e-4 instead of 1e-9.
#
# At every accepted step, through every deactivation event:
#
# - the radiation crossing the faces balances what the cells absorb: what
#   enters through the exposed face (`radiation_absorbed`), less what leaves
#   through it (`radiation_emitted`) and through the back (`radiation_back`),
#   equals the summed radiative sources of the live cells per unit area, to
#   round-off (an identity of the rows at any state);
# - what enters is the heater and the surroundings at the interface's
#   transmissivity, ε(q + σT∞⁴), ε = 0.96 of PMMA, whichever cell is exposed;
# - the live region's energy and mass close (`ConservationRecorder`,
#   `conservation_audit`): the change of the inventories equals what entered
#   less what the events took out, at the solve's tolerance.
#
# Reported beside: the radiation absorbed, emitted and lost through the back
# over the burn [J/m²].

# A ledger of the radiation at the accepted steps of one solve: the largest
# defect of the face balance against the cells' sources, the largest
# departure of the entering radiation from ε(q + σT∞⁴), and the trapezoidal
# integrals of the three crossings over time [J/m²].
if !isdefined(@__MODULE__, :_VRAD06Ledger)
    mutable struct _VRAD06Ledger
        worst_balance::Float64
        worst_entering::Float64
        absorbed::Float64
        emitted::Float64
        back::Float64
        steps::Int
        last::Tuple{Float64,Float64,Float64,Float64}
    end
    _VRAD06Ledger() = _VRAD06Ledger(0.0, 0.0, 0.0, 0.0, 0.0, 0, (NaN, NaN, NaN, NaN))
end

function _verification_case_vrad06()
    L, n = 4e-3, 40
    q = 100e3
    ε = 0.96
    entering = ε * (q + STEFAN_BOLTZMANN * _VRAD_GIDAZE_T0^4)
    last_run = Ref{Any}(nothing)   # the ledger and audit of the latest solve

    build(n_cells) = _vrad_gidaze_problem(
        collect(range(0.0, L; length = n_cells + 1)),
        z -> _VRAD_GIDAZE_T0,
        z -> _VRAD_GIDAZE_ρ;
        q,
        bottom = RadiativeBC(T_inf = _VRAD_GIDAZE_T0),
        tspan = (0.0, 300.0),
        output_times = collect(0.0:1.0:300.0),
    )

    # The ledger's entry at the accepted state u at time t, read with the
    # definition as it stands at that step (deactivation moves the kernel's
    # tables with the cells) into a private workspace.
    function record!(ledger, def, ws, u, t)
        kernel = find_kernel(RadiativeTransferKernel, def.kernels)
        topo = def.setup.topology
        ds, c = Pyrolysis.Assembly.buffers(ws, def, Float64)
        derive!(ds, u, def.layout, def.setup)
        any(ds.active) || return nothing
        Pyrolysis.Kernels.reset!(c)
        Pyrolysis.Kernels.contribute!(c, kernel, ds, topo, def.setup.material, t)
        A = ds.geometry.face_area[1]
        sources = sum(c.S_energy[i] for i in eachindex(ds.active) if ds.active[i]) / A
        x = radiation_crossing(kernel, ds, t)
        ledger.worst_balance =
            max(ledger.worst_balance, abs(sources - (x.absorbed - x.emitted - x.back)))
        ledger.worst_entering = max(ledger.worst_entering, abs(x.absorbed - entering))
        t0, a0, e0, b0 = ledger.last
        if isfinite(t0)
            ledger.absorbed += A * (a0 + x.absorbed) * (t - t0) / 2
            ledger.emitted += A * (e0 + x.emitted) * (t - t0) / 2
            ledger.back += A * (b0 + x.back) * (t - t0) / 2
        end
        ledger.last = (t, x.absorbed, x.emitted, x.back)
        ledger.steps += 1
        return nothing
    end

    function solve_fn(problem; kwargs...)
        prepared = setup(problem; energy_state = VERIFICATION_ENERGY_STATE, kwargs...)
        private = deepcopy(prepared.ws)
        ledger = _VRAD06Ledger()
        record!(ledger, prepared.def, private, prepared.u0, first(problem.tspan))
        step = DiscreteCallback(
            (u, t, integrator) -> true,
            integrator ->
                record!(ledger, prepared.def, private, integrator.u, integrator.t);
            save_positions = (false, false),
        )
        recorder = ConservationRecorder()
        # The cells reach the deactivation threshold, 1 % of their mass, at
        # 1 µm, the default thickness stop: below it, so the last cell is
        # deactivated rather than stopped.
        sol = integrate(prepared; callback = (recorder, step), min_thickness = 1e-8)
        last_run[] = (; sol, ledger, audit = conservation_audit(sol, recorder))
        return sol
    end

    function qois(sol)
        run = last_run[]
        run !== nothing && run.sol === sol || error("V-RAD-06: no ledger for this solve")
        ledger, audit = run.ledger, run.audit
        scale = ε * q
        return [
            QoIResult(
                name = "the sample burns out, every cell deactivated",
                value = length(sol.events),
                exact = n,
                error = (sol.outcome === :burned_out && length(sol.events) == n) ? 0.0 :
                        1.0,
                tol = 0.5,
                comment = "outcome $(sol.outcome), $(length(sol.events)) events",
            ),
            QoIResult(
                name = "radiation crossing the faces against the cells' sources, every accepted step",
                value = ledger.worst_balance / scale,
                error = ledger.worst_balance / scale,
                tol = 1e-12,
                comment = "absorbed − emitted − back against Σ S_rad / A over $(ledger.steps) states, over εq",
            ),
            QoIResult(
                name = "radiation entering, every accepted step",
                value = ledger.worst_entering / scale,
                error = ledger.worst_entering / scale,
                tol = 1e-12,
                comment = "ε(q + σT∞⁴) through whichever cell is exposed, over εq",
            ),
            QoIResult(
                name = "energy closure of the live region",
                value = audit.energy.relative,
                error = abs(audit.energy.relative),
                tol = 2.3e-7,  # 3× observed: 7.6e-8 (temperature form), 1.1e-8 (enthalpy)
                comment = "change − inflow + lost over the largest of them; inflow $(round(audit.energy.inflow, sigdigits = 4)) J",
            ),
            QoIResult(
                name = "mass closure of the live region",
                value = audit.mass.relative,
                error = abs(audit.mass.relative),
                tol = 4.6e-9,  # 3× observed: 1.5e-9 (temperature form), 1.0e-9 (enthalpy)
                comment = "lost $(round(audit.mass.lost, sigdigits = 4)) kg of $(round(audit.mass.initial, sigdigits = 4)) kg",
            ),
            QoIResult(
                name = "radiation absorbed over the burn [J/m²]",
                value = ledger.absorbed,
                error = 0.0,
                tol = Inf,
                comment = "reported",
            ),
            QoIResult(
                name = "radiation emitted through the exposed face over the burn [J/m²]",
                value = ledger.emitted,
                error = 0.0,
                tol = Inf,
                comment = "reported",
            ),
            QoIResult(
                name = "radiation leaving through the back over the burn [J/m²]",
                value = ledger.back,
                error = 0.0,
                tol = Inf,
                comment = "reported",
            ),
        ]
    end

    return VerificationCase(
        id = "V-RAD-06",
        title = "the radiative energy ledger through burnout",
        tags = [:radiation, :radiative_transfer, :burnout, :conservation],
        references = [
            "the rows' telescoping identity (exact at any state)",
            "conservation_audit with a ConservationRecorder (the live region's ledger)",
        ],
        build = build,
        solve_kwargs = (
            abstol = 1e-8,
            reltol = 1e-6,
            radiation_model = RADIATIVE_TRANSFER,
            use_ale = true,
            deactivation = true,
            pressure_tolerance = 1e-3,
        ),
        solve_fn = solve_fn,
        expected_retcode = :Terminated,
        qois = qois,
        ci = CIConfig(n_cells = n, budget_s = 5.0),
        report = (
            qoi_curves = (:radiation_crossing,),
            curves = sol -> begin
                x = sol.extras
                live = x.live_cells .> 0
                series(label, y) =
                    CurveSeries(label = label, x = sol.t[live], y = y[live])
                [
                    ReportCurve(
                        name = :radiation_crossing,
                        title = "V-RAD-06: the radiation crossing the faces through burnout",
                        xlabel = "t (s)",
                        ylabel = "flux (W/m²)",
                        series = [
                            series("entering the exposed face", x.radiation_absorbed),
                            series("leaving the exposed face", x.radiation_emitted),
                            series("leaving through the back", x.radiation_back),
                        ],
                    ),
                ]
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vrad06)
