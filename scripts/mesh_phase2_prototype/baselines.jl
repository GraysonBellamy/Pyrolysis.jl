# The baseline measurements of the Phase 2 plan (Appendix A.1, A.2, A.6 to A.9):
#   julia --project=test --startup-file=no --threads=1 scripts/mesh_phase2_prototype/baselines.jl <name>
# with <name> one of vbo_old, vbo_new, vbo_ladder_old, vbo01_qois_old, vbo02_restart,
# df_cores, integrators, tolweights. These override nothing and may share a session.

const NAME = isempty(ARGS) ? "vbo_old" : ARGS[1]

# The Douglas fir deck of examples/07_douglas_fir_cone.jl: the material, mesh and
# boundary conditions, without the example's printing, plotting and solve.
function douglas_fir_deck_source()
    lines = readlines(joinpath(@__DIR__, "..", "..", "examples", "07_douglas_fir_cone.jl"))
    a = findfirst(startswith("const T_HFM_MAX"), lines)
    b = findfirst(startswith("solution = solve("), lines)
    out = String[]
    skip = false
    for line in lines[a:(b-1)]
        if skip
            startswith(line, ")") && (skip = false)
            continue
        end
        if startswith(line, "println(")
            skip = !endswith(rstrip(line), ")")
            continue
        end
        startswith(line, "flush(") && continue
        push!(out, line)
    end
    return join(out, "\n")
end

if NAME == "vbo_old"
    ENV["PYROLYSIS_CORE"] = "old"
    ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
    empty!(ARGS)
    push!(ARGS, "V-BO-01,V-BO-02,V-BO-03,V-BO-04,V-BO-05")
    include(joinpath(@__DIR__, "cases.jl"))
end

if NAME == "vbo_new"
    ENV["PYROLYSIS_CORE"] = "new"
    ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
    empty!(ARGS)
    push!(ARGS, "V-BO-01,V-BO-02,V-BO-03,V-BO-04,V-BO-05")
    include(joinpath(@__DIR__, "cases.jl"))
end

if NAME == "vbo_ladder_old"
    ENV["PYROLYSIS_CORE"] = "old"
    ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
    using Pyrolysis, Statistics
    const REPO = normpath(joinpath(@__DIR__, "..", ".."))
    include(joinpath(REPO, "test", "verification", "casedefs.jl"))
    cases =
        filter(c -> c.id in ("V-BO-01", "V-BO-02", "V-BO-03"), load_verification_cases())
    for case in cases, n in (30, 60, 120)
        prob = case.build(n)
        _case_solve(case, prob; case.solve_kwargs...)
        ts = [
            (@elapsed sol = _case_solve(case, case.build(n); case.solve_kwargs...)) for
            _ = 1:3
        ]
        sol = _case_solve(case, case.build(n); case.solve_kwargs...)
        st = sol.raw_solution.stats
        println(
            "[$(case.id) old n_pmma=$n cells=$(length(prob.mesh.cells))] wall_median=$(round(median(ts), digits=3)) s of $(round.(ts, digits=3)) steps=$(st.naccept) nlfail=$(st.nnonlinconvfail) parks=$(length(sol.extras.park_times)) retcode=$(sol.retcode)",
        )
    end
end

if NAME == "vbo01_qois_old"
    ENV["PYROLYSIS_CORE"] = "old"
    ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
    using Pyrolysis, Statistics, Printf
    const REPO = normpath(joinpath(@__DIR__, "..", ".."))
    include(joinpath(REPO, "test", "verification", "casedefs.jl"))
    case = only(filter(c -> c.id == "V-BO-01", load_verification_cases()))
    for n in (30, 60, 120)
        sol = _case_solve(case, case.build(n); case.solve_kwargs...)
        m = Base.invokelatest(_vbo_measures, sol)
        @printf(
            "[V-BO-01 old n=%d] peak=%.3f t_peak=%.1f t_bo=%.1f mass_loss=%.4f ledger=%.2e L_end=%.6f parks=%d first_park=%.2f last_park=%.2f\n",
            n,
            m.peak,
            m.t_peak,
            m.t_bo,
            m.mass_loss,
            m.ledger,
            m.L_end,
            length(sol.extras.park_times),
            sol.extras.park_times[1],
            sol.extras.park_times[end]
        )
    end
end

if NAME == "vbo02_restart"
    ENV["PYROLYSIS_CORE"] = "new"
    ENV["PYROLYSIS_ENERGY_STATE"] = "temperature"
    using Pyrolysis, Statistics, Printf
    const Ig = Pyrolysis.Integration
    const REPO = normpath(joinpath(@__DIR__, "..", ".."))
    include(joinpath(REPO, "test", "verification", "casedefs.jl"))
    load_verification_cases()
    prob = Base.invokelatest(_vbo_problem, 60; substrate_alpha = 0.0, pmma_alpha = 0.0)
    for (tag, restart) in (("solver-chosen", :solver_chosen), ("forced 1e-6", :forced_dt))
        cp = Ig.setup(
            prob;
            use_ale = false,
            deactivation = true,
            radiation_model = SURFACE_ABSORPTION,
            abstol = 1e-8,
            reltol = 1e-6,
            time_scale = 1e-3,
        )
        pol = Ig.rodas5p_policy(Ig.default_linear_solver(); restart)
        sol = Ig.integrate(cp; integrator_policy = pol)
        t = minimum([
            (@elapsed sol = Ig.integrate(cp; integrator_policy = pol)) for _ = 1:3
        ])
        st = sol.statistics
        m = Base.invokelatest(_vbo_measures, sol)
        @printf(
            "[V-BO-02 new Rodas5P restart=%s] retcode=%s events=%d wall=%.3f s steps=%d rej=%d jac=%d nf=%d lin=%d peak=%.3f ledger=%.2e\n",
            tag,
            sol.retcode,
            length(sol.events),
            t,
            st.accepted,
            st.rejected,
            st.jacobians,
            st.rhs_evaluations,
            st.linear_solves,
            m.peak,
            m.ledger
        )
    end
    # KenCarp4 on the new core for the event-dense fixed-mesh deck, for the restart comparison
    cp = Ig.setup(
        prob;
        use_ale = false,
        deactivation = true,
        radiation_model = SURFACE_ABSORPTION,
        abstol = 1e-8,
        reltol = 1e-6,
        time_scale = 1e-3,
    )
    pol = Ig.kencarp4_policy(Ig.default_linear_solver())
    sol = Ig.integrate(cp; integrator_policy = pol)
    t = minimum([(@elapsed sol = Ig.integrate(cp; integrator_policy = pol)) for _ = 1:2])
    st = sol.statistics
    @printf(
        "[V-BO-02 new KenCarp4 raw] retcode=%s events=%d wall=%.3f s steps=%d rej=%d jac=%d nf=%d nlfail=%d\n",
        sol.retcode,
        length(sol.events),
        t,
        st.accepted,
        st.rejected,
        st.jacobians,
        st.rhs_evaluations,
        st.nonlinear_failures
    )
end

if NAME == "df_cores"
    using Pyrolysis, Statistics, Printf
    const T_FINAL = 60.0
    include_string(@__MODULE__, douglas_fir_deck_source())
    problem = PyrolysisProblem(
        mesh = mesh,
        material = douglas_fir,
        bc_set = bc_set,
        T_initial = T_initial_func,
        ξ_initial = ξ_initial_funcs,
        tspan = (0.0, T_FINAL),
        dt_initial = dt_initial,
        output_times = collect(0.0:1.0:T_FINAL),
    )
    println(
        "DF deck: ",
        length(problem.mesh.cells),
        " cells, ",
        Pyrolysis.n_components(douglas_fir),
        " components, ",
        Pyrolysis.n_reactions(douglas_fir),
        " reactions",
    )
    kw = (
        abstol = 1e-9,
        reltol = 1e-7,
        dtmin = 1e-14,
        use_ale = true,
        radiation_model = SURFACE_ABSORPTION,
    )
    # old core
    sol_old = solve(problem; progress = false, handle_depletion = false, kw...)
    t_old =
        @elapsed sol_old = solve(problem; progress = false, handle_depletion = false, kw...)
    st = sol_old.raw_solution.stats
    @printf(
        "[DF old] retcode=%s wall=%.2f s steps=%d nlfail=%d jacs=%d T_s(end)=%.2f thickness(end)=%.5f mlr(end)=%.4e\n",
        sol_old.retcode,
        t_old,
        st.naccept,
        st.nnonlinconvfail,
        st.njacs,
        sol_old.surface_T[end],
        sol_old.extras.thickness_history[end],
        sol_old.mass_loss_rate[end]
    )
    # new core
    try
        cp = Pyrolysis.Integration.setup(
            problem;
            use_ale = true,
            radiation_model = SURFACE_ABSORPTION,
            abstol = 1e-9,
            reltol = 1e-7,
            time_scale = dt_initial,
        )
        println(
            "[DF new] colors = ",
            cp.jacobian.n_colors,
            " state length = ",
            length(cp.u0),
        )
        sol_new = Pyrolysis.Integration.integrate(cp; dtmin = 1e-14)
        t_new = @elapsed sol_new = Pyrolysis.Integration.integrate(cp; dtmin = 1e-14)
        s = sol_new.statistics
        @printf(
            "[DF new] retcode=%s wall=%.2f s steps=%d rej=%d jacs=%d cres=%.2e T_s(end)=%.2f thickness(end)=%.5f mlr(end)=%.4e\n",
            sol_new.retcode,
            t_new,
            s.accepted,
            s.rejected,
            s.jacobians,
            s.largest_constraint_residual,
            sol_new.surface_T[end],
            sol_new.extras.thickness_history[end],
            sol_new.mass_loss_rate[end]
        )
        dT = maximum(abs.(sol_new.surface_T .- sol_old.surface_T))
        println(
            "[DF] max |ΔT_s| between cores over outputs = ",
            dT,
            " K; thickness diff = ",
            abs(
                sol_new.extras.thickness_history[end] -
                sol_old.extras.thickness_history[end],
            ),
        )
    catch err
        println("[DF new] ERROR: ", sprint(showerror, err)[1:min(end, 3000)])
    end
end

if NAME == "integrators"
    ENV["PYROLYSIS_CORE"] = "new"
    using Pyrolysis, Statistics, Printf, OrdinaryDiffEq
    const AutoFiniteDiff = Pyrolysis.Integration.AutoFiniteDiff
    const Ig = Pyrolysis.Integration
    const REPO = normpath(joinpath(@__DIR__, "..", ".."))
    include(joinpath(REPO, "scripts", "mesh_phase0", "baselines", "jacobian_common.jl"))
    _mod(name) = first(m for (k, m) in Base.loaded_modules if k.name == name)
    RB = _mod("OrdinaryDiffEqRosenbrock")
    SD = _mod("OrdinaryDiffEqSDIRK")
    println(
        "Rosenbrock methods: ",
        filter(n -> occursin(r"^(Rodas|ROS|Ros|Velds|GRK|Scholz)", string(n)), names(RB)),
    )
    candidates = Any[
        ("Rodas5P", ls -> RB.Rodas5P(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        ("Rodas5Pe", ls -> RB.Rodas5Pe(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        ("Rodas5Pr", ls -> RB.Rodas5Pr(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        ("Rodas4P2", ls -> RB.Rodas4P2(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        ("Rodas4P", ls -> RB.Rodas4P(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        ("Rodas3P", ls -> RB.Rodas3P(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        ("ROS34PW2", ls -> RB.ROS34PW2(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        ("Rodas23W", ls -> RB.Rodas23W(autodiff = AutoFiniteDiff(), linsolve = ls), :clear),
        (
            "KenCarp4 raw",
            ls -> SD.KenCarp4(
                autodiff = AutoFiniteDiff(),
                linsolve = ls,
                smooth_est = false,
            ),
            :keep,
        ),
        (
            "KenCarp4 smooth",
            ls -> SD.KenCarp4(autodiff = AutoFiniteDiff(), linsolve = ls),
            :keep,
        ),
        ("Kvaerno5", ls -> SD.Kvaerno5(autodiff = AutoFiniteDiff(), linsolve = ls), :keep),
        (
            "KenCarp47",
            ls -> SD.KenCarp47(autodiff = AutoFiniteDiff(), linsolve = ls),
            :keep,
        ),
    ]
    function new_solve(problem, kw, alg, flag, reltol; policy_restart = :solver_chosen)
        cp = Ig.setup(
            problem;
            use_ale = kw.use_ale,
            radiation_model = kw.radiation_model,
            reltol,
            abstol = reltol / 100,
            time_scale = 1e-3,
        )
        ls = Ig.default_linear_solver()
        pol = Ig.IntegratorPolicy(alg(ls), policy_restart, flag, 1e-6)
        return Ig.integrate(cp; integrator_policy = pol)
    end
    for (kind, n) in ((:ale, 100), (:eulerian, 100), (:ale_darcy, 100))
        problem, kw = problem_for(kind, n; solid_λ = kind == :ale_darcy ? nothing : 1e-6)
        ref = new_solve(problem, kw, candidates[1][2], :clear, 1e-10)
        old_ref = solve(problem; kw..., reltol = 1e-10, abstol = 1e-12)
        old = solve(problem; kw..., reltol = 1e-6, abstol = 1e-8)
        t_old = minimum([
            (@elapsed solve(problem; kw..., reltol = 1e-6, abstol = 1e-8)) for _ = 1:3
        ])
        target = maximum(abs, old.T .- old_ref.T)
        @printf(
            "[%s %d] old core default: error %.3e K, wall %.4f s, steps %d\n",
            kind,
            n,
            target,
            t_old,
            old.raw_solution.stats.naccept
        )
        for (name, alg, flag) in candidates
            for reltol in (1e-6, 3e-7, 1e-7)
                try
                    sol = new_solve(problem, kw, alg, flag, reltol)
                    t = minimum([
                        (@elapsed new_solve(problem, kw, alg, flag, reltol)) for _ = 1:2
                    ])
                    err = maximum(abs, sol.T .- ref.T)
                    st = sol.statistics
                    @printf(
                        "[%s %d] %-16s reltol=%.0e retcode=%s error=%.3e K (%.2f× target) wall=%.4f s (%.2f× old) steps=%d rej=%d jac=%d nf=%d nlfail=%d\n",
                        kind,
                        n,
                        name,
                        reltol,
                        sol.retcode,
                        err,
                        err / target,
                        t,
                        t / t_old,
                        st.accepted,
                        st.rejected,
                        st.jacobians,
                        st.rhs_evaluations,
                        st.nonlinear_failures
                    )
                catch e
                    println(
                        "[$kind $n] $name reltol=$reltol ERROR: ",
                        sprint(showerror, e)[1:min(end, 300)],
                    )
                end
                flush(stdout)
            end
        end
    end
end

if NAME == "tolweights"
    ENV["PYROLYSIS_CORE"] = "new"
    using Pyrolysis, Statistics, Printf
    const Ig = Pyrolysis.Integration
    const St = Pyrolysis.States
    const REPO = normpath(joinpath(@__DIR__, "..", ".."))
    include(joinpath(REPO, "scripts", "mesh_phase0", "baselines", "jacobian_common.jl"))
    function run(problem, kw, reltol; vel = 1.0, pos = 1.0, mass = 1.0, energy = 1.0)
        cp = Ig.setup(
            problem;
            use_ale = kw.use_ale,
            radiation_model = kw.radiation_model,
            reltol,
            abstol = reltol / 100,
            time_scale = 1e-3,
        )
        layout = cp.def.layout
        if St.is_moving(layout)
            for nd = 1:St.n_nodes(layout)
                cp.abstol[St.velocity_index(layout, nd)] *= vel
                cp.abstol[St.position_index(layout, nd)] *= pos
            end
        end
        for i = 1:St.n_cells(layout)
            cp.abstol[St.energy_index(layout, i)] *= energy
            for j = 1:St.n_components(layout)
                cp.abstol[St.mass_index(layout, i, j)] *= mass
            end
        end
        return Ig.integrate(cp)
    end
    for kind in (:ale, :eulerian)
        problem, kw = problem_for(kind, 100; solid_λ = 1e-6)
        ref = run(problem, kw, 1e-10)
        old_ref = solve(problem; kw..., reltol = 1e-10, abstol = 1e-12)
        old = solve(problem; kw..., reltol = 1e-6, abstol = 1e-8)
        target = maximum(abs, old.T .- old_ref.T)
        @printf(
            "[%s 100] old default error %.3e K steps %d\n",
            kind,
            target,
            old.raw_solution.stats.naccept
        )
        for (tag, w) in (
            ("baseline", (;)),
            ("vel×1e2", (; vel = 1e2)),
            ("vel×1e4", (; vel = 1e4)),
            ("pos×1e2", (; pos = 1e2)),
            ("mass×1e2", (; mass = 1e2)),
            ("mass×1e4", (; mass = 1e4)),
            ("energy×1e2", (; energy = 1e2)),
            ("energy×1e-2", (; energy = 1e-2)),
            ("mass×1e4 vel×1e4", (; mass = 1e4, vel = 1e4)),
        )
            for reltol in (1e-6, 3e-7)
                sol = run(problem, kw, reltol; w...)
                err = maximum(abs, sol.T .- ref.T)
                @printf(
                    "[%s 100] %-18s reltol=%.0e error=%.3e K (%.2f× target) steps=%d rej=%d\n",
                    kind,
                    tag,
                    reltol,
                    err,
                    err / target,
                    sol.statistics.accepted,
                    sol.statistics.rejected
                )
            end
        end
        flush(stdout)
    end
end
