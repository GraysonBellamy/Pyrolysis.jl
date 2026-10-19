# The prototype's experiments, by name (Appendix A of the Phase 2 plan):
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase2_prototype/run.jl <experiment>
# with <experiment> one of check, sweep, band, repeat, ladder, order, continuous, scaling.
# Each run is a fresh session; the overrides of common.jl stay in it.

include(joinpath(@__DIR__, "common.jl"))
const EXPERIMENT = isempty(ARGS) ? "check" : ARGS[1]

if EXPERIMENT == "check"
    # containment of the widened pattern with every node partially activated (r_hi = 1.5 > 1 > r_lo = 0.5,
    # so the activation and its derivatives are exercised) and the monitor's mass term on
    cp2, prob = build_vbo(4; τ = 1.0, r_hi = 1.5, r_lo = 0.5, β = 1.0)
    println(
        "cells = ",
        St.n_cells(cp2.def.layout),
        " colors = ",
        cp2.jacobian.n_colors,
        " state length = ",
        length(cp2.u0),
    )
    u = copy(cp2.u0)
    # perturb the state so the activation and fluxes are generic
    layout = cp2.def.layout
    for i = 1:St.n_cells(layout)
        u[St.energy_index(layout, i)] += 50.0 * sin(i)
        u[St.position_index(layout, i + 1)] *= 1 + 0.05 * cos(i)
        for j = 1:St.n_components(layout)
            u[St.mass_index(layout, i, j)] *= 1 + 0.1 * sin(3i + j)
        end
    end
    Ig.consistent_state!(u, cp2.def, cp2.ws, 0.0)
    Jd = CJ.dense_jacobian(cp2.def, cp2.ws, u, 0.0)
    miss = CJ.missing_entries(Jd, cp2.jacobian.pattern)
    println("missing entries: ", length(miss), " ", first(miss, 10))
    J = CJ.jac_prototype(cp2.jacobian)
    CJ.jacobian!(J, cp2.jacobian, u, 0.0)
    println("colored vs dense: ", CJ.compare_entries(Matrix(J), Jd; rtol = 1e-10))
    # sanity: with the correction off the pattern equals the package's
    CORR.on = false
    println(
        "c at node 5 with correction off: ",
        Pyrolysis.Closures.corr_c(cp2.ws.derived, cp2.def.setup.topology, 5),
    )
end

if EXPERIMENT == "sweep"
    BLOCK[] = true
    const BL = (; radiation_model = BEER_LAMBERT)
    for (tag, kw) in (
        ("c=0, block rule", (; on = false, case_kw = BL)),
        ("tau=1", (; τ = 1.0, case_kw = BL)),
        ("tau=0.2", (; τ = 0.2, case_kw = BL)),
        ("tau=5", (; τ = 5.0, case_kw = BL)),
    )
        r = Base.invokelatest(run_vbo, 60; kw...)
        Base.invokelatest(report, "V-BO-01 BL n=60 $(tag)", r)
        Base.invokelatest(vbo_qois, r.sol)
        flush(stdout)
    end
end

if EXPERIMENT == "band"
    BLOCK[] = true
    const BL = (; radiation_model = BEER_LAMBERT)
    for (κ, τ, n) in (
        (3.0, 0.2, 60),
        (3.0, 1.0, 60),
        (10.0, 0.2, 60),
        (10.0, 1.0, 60),
        (3.0, 0.2, 30),
        (3.0, 0.2, 120),
        (3.0, 1.0, 120),
    )
        BAND[] = κ
        r = Base.invokelatest(run_vbo, n; τ, case_kw = BL)
        Base.invokelatest(report, "V-BO-01 BL n=$(n) tau=$(τ) band=$(κ)", r)
        Base.invokelatest(vbo_qois, r.sol)
        flush(stdout)
    end
end

if EXPERIMENT == "repeat"
    BLOCK[] = true
    const BL = (; radiation_model = BEER_LAMBERT)
    for (κ, τ, n, nrep) in ((10.0, 1.0, 60, 3), (10.0, 0.2, 60, 3))
        BAND[] = κ
        r = Base.invokelatest(run_vbo, n; τ, case_kw = BL, nrep)
        Base.invokelatest(report, "V-BO-01 BL n=$(n) tau=$(τ) band=$(κ) nrep=$(nrep)", r)
        Base.invokelatest(vbo_qois, r.sol)
        flush(stdout)
    end
end

if EXPERIMENT == "ladder"
    BLOCK[] = true
    const BL = (; radiation_model = BEER_LAMBERT)
    for (tag, κ, n, kw) in (
        ("tau=0.2 band=10", 10.0, 120, (; τ = 0.2)),
        ("tau=0.2 band=10", 10.0, 30, (; τ = 0.2)),
        ("tau=0.1 band=10", 10.0, 60, (; τ = 0.1)),
        ("tau=0.2 band=10 beta=1", 10.0, 60, (; τ = 0.2, β = 1.0)),
        ("tau=0.2 band=10 and-activation", 10.0, 60, (; τ = 0.2, mode = :and)),
        ("tau=0.2 band=10 eps_m=3e-3", 10.0, 60, (; τ = 0.2, threshold = 3e-3)),
        ("tau=0.2 band=10 eps_m=3e-2", 10.0, 60, (; τ = 0.2, threshold = 3e-2)),
        (
            "tau=0.2 band=10 activation 0.5-0.2",
            10.0,
            60,
            (; τ = 0.2, r_hi = 0.5, r_lo = 0.2),
        ),
        (
            "tau=0.2 band=10 enthalpy",
            10.0,
            60,
            (;
                τ = 0.2,
                case_kw = (; radiation_model = BEER_LAMBERT, energy_state = :enthalpy),
            ),
        ),
    )
        BAND[] = κ
        ckw = haskey(kw, :case_kw) ? kw : merge(kw, (; case_kw = BL))
        r = try
            Base.invokelatest(run_vbo, n; ckw...)
        catch err
            println(
                "[V-BO-01 BL n=$(n) $(tag)] ERROR: ",
                sprint(showerror, err)[1:min(end, 400)],
            )
            flush(stdout)
            continue
        end
        Base.invokelatest(report, "V-BO-01 BL n=$(n) $(tag)", r)
        Base.invokelatest(vbo_qois, r.sol)
        flush(stdout)
    end
end

if EXPERIMENT == "order"
    BLOCK[] = true
    const BL = (; radiation_model = BEER_LAMBERT)
    for (κ, τ) in ((3.0, 0.2), (3.0, 1.0), (10.0, 1.0))
        BAND[] = κ
        r = Base.invokelatest(run_vbo, 60; τ, case_kw = BL, nrep = 2)
        Base.invokelatest(report, "order check: V-BO-01 BL n=60 tau=$(τ) band=$(κ)", r)
        flush(stdout)
    end
end

if EXPERIMENT == "continuous"
    isdefined(Pyrolysis.Integration, :CONT) ||
        @eval Pyrolysis.Integration const CONT = Ref(false)
    const CONT = Pyrolysis.Integration.CONT
    CORR.on = false
    BLOCK[] = true
    @eval Pyrolysis.Integration function _deactivation_callback(
        cp::CoreProblem,
        check::DomainCheck,
        policy,
    )
        cells = cp.deactivation
        cells === nothing && return ()
        def, ws = cp.def, cp.ws
        layout = def.layout
        has_event =
            FunctionWrapper{Bool,Tuple{_StateVector}}(u -> any_below(cells, layout, u))
        deactivate = FunctionWrapper{Symbol,Tuple{_StateVector,Float64}}(
            (u, t) -> begin
                status = deactivate!(u, cells, def, t)
                status === :transferred || return status
                _consistent_intensities!(u, def, ws, t)
                _consistent_motion!(u, def, ws, t) === nothing || return :inconsistent
                inadmissible(u, layout, check.abstol) && return :inadmissible
                _porosity_out(check, _derived(check, u, t), def.setup.gas) &&
                    return :inadmissible
                return status
            end,
        )
        # at a located root the crossing cell sits at the threshold: nudge it below
        # by four ulps so the strict predicate sees it, then deactivate
        nudge = FunctionWrapper{Nothing,Tuple{_StateVector}}(
            u -> begin
                best, gbest = 0, Inf
                for i = 1:n_cells(layout)
                    layout.mask[energy_index(layout, i)] || continue
                    cells.reference[i] > 0 || continue
                    g =
                        condensed_mass(cells, layout, u, i) /
                        (cells.threshold * cells.reference[i]) - 1
                    g < gbest && (best, gbest = i, g)
                end
                if best > 0 && !below_threshold(cells, layout, u, best)
                    for j = 1:n_components(layout)
                        cells.solid[j] && (u[mass_index(layout, best, j)] *= (1 - 4eps()))
                    end
                end
                nothing
            end,
        )
        function affect!(integrator)
            CONT[] && nudge(integrator.u)
            status = deactivate(integrator.u, integrator.t)
            if status === :inadmissible
                @warn "solve(): the state after deactivating cell(s) at t = $(integrator.t) is inadmissible"
                terminate!(integrator, ReturnCode.Failure)
            elseif status === :inconsistent
                @warn "solve(): no consistent node velocities after deactivating cell(s) at t = $(integrator.t)"
                terminate!(integrator, ReturnCode.Failure)
            elseif status !== :transferred
                terminate!(integrator)
            elseif policy.restart === :forced_dt
                set_proposed_dt!(integrator, policy.restart_dt)
            end
            return nothing
        end
        if CONT[]
            # the smallest margin above the threshold over the live cells, scaled
            margin = FunctionWrapper{Float64,Tuple{_StateVector}}(
                u -> begin
                    g = Inf
                    for i = 1:n_cells(layout)
                        layout.mask[energy_index(layout, i)] || continue
                        cells.reference[i] > 0 || continue
                        g = min(
                            g,
                            condensed_mass(cells, layout, u, i) /
                            (cells.threshold * cells.reference[i]) - 1,
                        )
                    end
                    g
                end,
            )
            condition_d = (u, t, integrator) -> has_event(u)
            return (
                ContinuousCallback(
                    (u, t, integrator) -> margin(u),
                    nothing;
                    affect_neg! = affect!,
                    rootfind = SciMLBase.RightRootFind,
                    save_positions = (false, false),
                ),
                DiscreteCallback(condition_d, affect!; save_positions = (false, false)),
            )
        end
        condition = (u, t, integrator) -> has_event(u)
        return (DiscreteCallback(condition, affect!; save_positions = (false, false)),)
    end

    prob = Base.invokelatest(_vbo_problem, 60; substrate_alpha = 0.0, pmma_alpha = 0.0)
    results = Dict()
    for cont in (false, true), reltol in (1e-6, 1e-8)
        CONT[] = cont
        kw = (
            deactivation = true,
            use_ale = false,
            radiation_model = SURFACE_ABSORPTION,
            abstol = reltol * 1e-2,
            reltol,
        )
        sol = Base.invokelatest(Ig.solve, prob; kw...)
        t = minimum([(@elapsed sol = Base.invokelatest(Ig.solve, prob; kw...)) for _ = 1:2])
        st = sol.statistics
        results[(cont, reltol)] = [e.t for e in sol.events]
        @printf(
            "[V-BO-02 %s reltol=%.0e] retcode=%s events=%d wall=%.3f s steps=%d rej=%d jac=%d nf=%d first=%.3f last=%.3f\n",
            cont ? "continuous" : "discrete",
            reltol,
            sol.retcode,
            length(sol.events),
            t,
            st.accepted,
            st.rejected,
            st.jacobians,
            st.rhs_evaluations,
            sol.events[1].t,
            sol.events[end].t
        )
        m = Base.invokelatest(_vbo_measures, sol)
        @printf(
            "    peak=%.3f t_peak=%.1f mass_loss=%.4f ledger=%.2e\n",
            m.peak,
            m.t_peak,
            m.mass_loss,
            m.ledger
        )
        flush(stdout)
    end
    for cont in (false, true)
        a, b = results[(cont, 1e-6)], results[(cont, 1e-8)]
        n = min(length(a), length(b))
        a = a[1:n]
        b = b[1:n]
        interval = median(diff(b))
        @printf(
            "[event-time shift %s] max |t(1e-6) − t(1e-8)| = %.4f s = %.3f of the median interval %.3f s\n",
            cont ? "continuous" : "discrete",
            maximum(abs.(a[1:n] .- b[1:n])),
            maximum(abs.(a[1:n] .- b[1:n])) / interval,
            interval
        )
    end
    a, b = results[(false, 1e-8)], results[(true, 1e-8)]
    n = min(length(a), length(b))
    a = a[1:n]
    b = b[1:n]
    @printf(
        "[discrete vs continuous at 1e-8] max |Δt_event| = %.4f s\n",
        maximum(abs.(a .- b))
    )
end

if EXPERIMENT == "scaling"
    # the N² scaling of the relaxation: τ scaled by (60/N)² should give the 60-cell behaviour at 30 and 120 cells
    BLOCK[] = true
    BAND[] = 10.0
    const BL = (; radiation_model = BEER_LAMBERT)
    for (n, τ) in ((120, 0.05), (30, 0.8), (120, 0.1))
        r = Base.invokelatest(run_vbo, n; τ, case_kw = BL)
        Base.invokelatest(report, "V-BO-01 BL n=$(n) tau=$(τ) band=10 (N² scaling)", r)
        Base.invokelatest(vbo_qois, r.sol)
        flush(stdout)
    end
end
