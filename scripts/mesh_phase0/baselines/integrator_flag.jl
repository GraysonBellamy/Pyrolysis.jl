# Step-update flag × integrator × tolerance matrix on the finite-slab
# heat-conduction cases.
#
# The package's step-update callback mirrors the accepted state into the mesh
# and flags a derivative discontinuity after every accepted step. This driver
# replaces that callback in-process with one of three variants and measures,
# for each integrator and tolerance pair, the integrator counters, the warm
# wall time and the profile error against the exact series:
#
#   current   mirror + derivative_discontinuity!(integrator, true)   (package)
#   flag_off  mirror, no flag call (the integrator's own default applies)
#   observer  mirror + derivative_discontinuity!(integrator, false)
#
# `u_modified!` is a deprecated alias of `derivative_discontinuity!` in the
# pinned SciMLBase, so the observer variant calls the latter directly.
#
# Output: scripts/mesh_phase0/results/integrator_flag.csv

include(joinpath(@__DIR__, "common.jl"))

include(joinpath(ROOT, "test", "verification", "cases", "V-HC-03.jl"))
include(joinpath(ROOT, "test", "verification", "cases", "V-HC-02.jl"))

const STEP_VARIANT = Ref(:current)
const MAXITERS = 20_000   # successful runs need < 300 steps; a stalled run stops here

# In-process replacement of the step-update callback (same signature as the
# package method; `solve` builds it at every call, so the variant is read at
# solve time). The only addition to the package code is the wall-clock guard.
@eval Pyrolysis.Solver function create_step_update_callback(ws, material)
    mesh = ws.mesh
    sync_mesh! = FunctionWrapper{Nothing,Tuple{_StateVector}}(
        u -> begin
            layout = ws.layout
            mesh.cross_section_area =
                cross_section_area_from_state(mesh, material, u, layout)
            update_mesh_from_state!(mesh, u, layout)
        end,
    )
    condition = (u, t, integrator) -> true
    variant = Main.STEP_VARIANT[]
    deadline = Main.DEADLINE
    timed_out = Main.TIMED_OUT
    last_t = Main.LAST_T
    function affect!(integrator)
        sync_mesh!(integrator.u)
        last_t[] = integrator.t
        if time() > deadline[]
            timed_out[] = true
            terminate!(integrator)
        end
        if variant === :current
            derivative_discontinuity!(integrator, true)
        elseif variant === :observer
            derivative_discontinuity!(integrator, false)
        end
        return nothing
    end
    return DiscreteCallback(condition, affect!; save_positions = (false, false))
end

# ----------------------------------------------------------------------------
# Cases and their exact profiles T(z, t) (constants as in the case files).
# ----------------------------------------------------------------------------
function exact_profile(case_id::String)
    if case_id == "V-HC-03"
        ρ, c, k, L, T0, q = 1000.0, 1000.0, 0.5, 0.02, 300.0, 5000.0
        α = k / (ρ * c)
        return (z, t) -> T0 + (q * L / k) * slab_constflux_phi(z / L, α * t / L^2)
    elseif case_id == "V-HC-02c"
        ρ, c, k, L, T_i, T_inf, Bi = 1000.0, 1000.0, 0.5, 0.1, 300.0, 500.0, 10.0
        α = k / (ρ * c)
        λs = biot_eigenvalues(Bi, 60)
        return (z, t) ->
            T_inf + (T_i - T_inf) * slab_convective_theta(z / L, α * t / L^2, λs)
    end
    error("no exact profile for $case_id")
end

# RMS and maximum profile error at every output time after t = 0.
function profile_errors(sol, T_exact)
    z = get_cell_positions(sol)
    L2 = Float64[]
    Linf = Float64[]
    for j = 2:length(sol.t)
        e = [sol.T[i, j] - T_exact(z[i], sol.t[j]) for i in eachindex(z)]
        push!(L2, sqrt(sum(abs2, e) / length(e)))
        push!(Linf, maximum(abs, e))
    end
    return L2, Linf
end

const CASES = Dict(
    "V-HC-03" => _verification_case_vhc03(),
    "V-HC-02c" => _vhc02_case("c", 10.0; with_convergence = true),
)

const INTEGRATORS = (
    (
        "KenCarp4",
        () -> KenCarp4(autodiff = AutoFiniteDiff(), linsolve = KLUFactorization()),
    ),
    (
        "KenCarp4_raw",
        () -> KenCarp4(
            smooth_est = false,
            autodiff = AutoFiniteDiff(),
            linsolve = KLUFactorization(),
        ),
    ),
    ("FBDF", () -> FBDF(autodiff = AutoFiniteDiff(), linsolve = KLUFactorization())),
    ("Rodas5P", () -> Rodas5P(autodiff = AutoFiniteDiff(), linsolve = KLUFactorization())),
)
const VARIANTS = (:current, :flag_off, :observer)
const TOLERANCES = ((1e-6, 1e-8), (1e-8, 1e-10))   # (reltol, abstol)
const RUNS = (("V-HC-03", 320), ("V-HC-03", 80), ("V-HC-02c", 320))

const MAX_OUTPUT_TIMES = 4
const HEADER = (
    "case",
    "n_cells",
    "integrator",
    "callback_variant",
    "reltol",
    "abstol",
    "retcode",
    "timed_out",
    "message",
    "steps",
    "rejects",
    "nf",
    "njacs",
    "nw",
    "nsolve",
    "nonliniter",
    "nonlinfail",
    "first_s",
    "warm_median_s",
    "warm_min_s",
    "warm_max_s",
    "repeats_consistent",
    "output_times_or_t_reached",
    ("L2_t$i" for i = 1:MAX_OUTPUT_TIMES)...,
    ("Linf_t$i" for i = 1:MAX_OUTPUT_TIMES)...,
    "L2_max",
    "Linf_max",
    "L2_final",
    "Linf_final",
    "qoi_failures",
    "qoi_max_error_over_tol",
)

padded(v) = vcat(v, fill(NaN, MAX_OUTPUT_TIMES - length(v)))

function run_matrix(io)
    checked = Set{String}()
    for run in RUNS, intg in INTEGRATORS, variant in VARIANTS, tol in TOLERANCES
        (case_id, n), (int_name, make_alg), (reltol, abstol) = run, intg, tol
        case = CASES[case_id]
        T_exact = exact_profile(case_id)
        STEP_VARIANT[] = variant
        say("RUN ", (case_id, n, int_name, variant, reltol, abstol))
        r = timed_runs(
            case,
            n;
            integrator = make_alg(),
            reltol = reltol,
            abstol = abstol,
            maxiters = MAXITERS,
        )
        st = solve_stats(r.sol)
        finished = r.sol !== nothing && r.retcode == case.expected_retcode && !r.timed_out
        if !finished
            # A solve that stopped early has no comparable profiles; record
            # the counters, the return code and the time it reached.
            L2 = Float64[]
            Linf = Float64[]
            times = "reached t=" * string(LAST_T[])
            qoi_fail = -1
            qoi_ratio = NaN
        else
            L2, Linf = profile_errors(r.sol, T_exact)
            times = join(string.(r.sol.t[2:end]), ";")
            qois = Base.invokelatest(case.qois, r.sol)
            qoi_fail = count(q -> !(q.error < q.tol), qois)
            qoi_ratio = maximum(q.error / q.tol for q in qois)
            # Guard against constant drift between this file and the case
            # file: the case's own L∞ QoI at the first output time must match.
            if !(case_id in checked)
                q1 = first(qois)
                abs(q1.value - Linf[1]) <= 1e-9 * max(1.0, q1.value) || error(
                    "exact-profile mismatch for $case_id: case QoI $(q1.value) vs " *
                    "driver $(Linf[1])",
                )
                push!(checked, case_id)
            end
        end
        csv_row(
            io,
            case_id,
            n,
            int_name,
            string(variant),
            reltol,
            abstol,
            string(r.retcode),
            r.timed_out,
            r.message,
            st.steps,
            st.rejects,
            st.nf,
            st.njacs,
            st.nw,
            st.nsolve,
            st.nonliniter,
            st.nonlinfail,
            r.first_s,
            r.warm_median_s,
            r.warm_min_s,
            r.warm_max_s,
            r.repeats_consistent,
            times,
            padded(L2)...,
            padded(Linf)...,
            isempty(L2) ? NaN : maximum(L2),
            isempty(Linf) ? NaN : maximum(Linf),
            isempty(L2) ? NaN : L2[end],
            isempty(Linf) ? NaN : Linf[end],
            qoi_fail,
            qoi_ratio,
        )
        say(
            "   -> ",
            (r.retcode, st.steps, st.rejects, st.nonlinfail, r.warm_median_s),
            " L2max=",
            isempty(L2) ? NaN : maximum(L2),
            " Linf_final=",
            isempty(Linf) ? NaN : Linf[end],
        )
    end
end

t_total = @elapsed open(joinpath(RESULTS_DIR, "integrator_flag.csv"), "w") do io
    csv_row(io, HEADER...)
    Base.invokelatest(run_matrix, io)
end
say("TOTAL wall time (s): ", t_total)
