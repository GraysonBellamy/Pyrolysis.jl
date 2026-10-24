# =============================================================================
# V-MMS-02 — ALE free-stream preservation under prescribed mesh motion
# =============================================================================
#
# The sharpest possible discrete-GCL assert: spatially uniform T*, ξ* with no
# reactions and zero-flux manufactured BCs, while the mesh is compressed to
# 70 % of its initial thickness through the `Workspace.w_override` hook with
# the uniform-stretch velocity
#
#   w*(z, t) = (ṡ/s)·z,   s(t) = 1 − 0.15·(1 − cos(2πt/20)),   z_bottom pinned
#
# (node trajectories z_i(t) = z_i(0)·s(t) exactly). A uniform state must stay
# uniform to machine precision on the moving mesh: the ALE advection of a
# constant field must vanish discretely (limiter slope ≡ 0), the dilation
# correction is zero (θ reaction-driven, no reactions), and the conductive/
# species fluxes are zero. Any drift is a moving-mesh transport defect that
# smooth-field order studies (V-MMS-03) can hide inside their truncation
# error.
#
# The node-motion ODE ż_i = w*(z_i,t) is integrated by the same stiff
# integrator, so the top-cell-center trajectory doubles as an exactness check
# on the prescribed-motion path (temporal tolerance only).
#
# Jacobian: colored/dense FD (`jacobian = :fd`) — the structured geometry
# Jacobian assumes reaction-driven mesh velocity and does not model the
# override (documented on `Workspace.w_override`).

function _verification_case_vmms02()
    L = 0.01
    t_end = 10.0
    T0 = 500.0
    ξs0, ξg0 = 600.0, 1.5

    s(t) = 1.0 - 0.15 * (1.0 - cos(2π * t / 20.0))
    sdot(t) = -0.15 * (2π / 20.0) * sin(2π * t / 20.0)
    w_star(z, t) = (sdot(t) / s(t)) * z

    fields = MMSFields((z, t) -> T0, ((z, t) -> ξs0, (z, t) -> ξg0))

    function make_material()
        solid = SolidComponent(:s, ρ = 800.0, c = 1000.0, k = 0.3, λ = 2.0e-5)
        gas = GaseousComponent(:g, M = 0.028, c = 1000.0, k = 0.03, λ = 2.0e-5)
        return Material(name = :VMMS02, components = (solid, gas), reactions = ())
    end

    function build(n_cells)
        mat = make_material()
        mesh = core_mesh(
            create_uniform_mesh,
            L,
            n_cells,
            Val(2);
            T_initial = T0,
            ξ_initial = (ξs0, ξg0),
        )
        # Uniform fields ⇒ every manufactured flux is identically zero; the
        # BCs still go through the manufactured-flux callables so the whole
        # path is the one V-MMS-03 uses.
        thermal = mms_thermal_bcs(mat, fields, L)
        mass = mms_mass_bcs(mat, fields, L)
        bc_set = BoundaryConditionSet(
            thermal = (top = thermal.top, bottom = thermal.bottom),
            mass = (top = mass.top, bottom = mass.bottom),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T0,
            ξ_initial = [z -> ξs0, z -> ξg0],
            tspan = (0.0, t_end),
            output_times = collect(range(0.0, t_end, length = 41)),
        )
    end

    solve_fn(problem; kwargs...) = solve_mms(
        problem,
        fields;
        use_ale = true,
        w_override = w_star,
        jacobian = :fd,
        core = Symbol(VERIFICATION_CORE),
        energy_state = VERIFICATION_ENERGY_STATE,
        kwargs...,
    )

    function qois(sol)
        out = QoIResult[]
        T_drift = maximum(abs.(sol.T .- T0))
        push!(
            out,
            QoIResult(
                name = "free-stream T drift (all cells, all snapshots)",
                value = T_drift,
                exact = 0.0,
                error = T_drift,
                tol = 1e-10,  # observed exactly 0.0 (n=20)
                comment = "uniform field must survive 30% mesh compression",
            ),
        )
        for (name, idx, scale) in (("ξ_solid", 1, ξs0), ("ξ_gas", 2, ξg0))
            drift = maximum(abs.(sol.ξ[idx, :, :] .- scale)) / scale
            push!(
                out,
                QoIResult(
                    name = "free-stream $(name) relative drift",
                    value = drift,
                    exact = 0.0,
                    error = drift,
                    tol = 1e-12,  # observed 0.0 (solid) / 1.5e-16 (gas) at n=20
                    comment = "ALE advection + dilation must vanish for a uniform state",
                ),
            )
        end
        # Prescribed-motion exactness: top cell center must track
        # z(0)·s(t) (temporal integration error only).
        n = size(sol.z, 1)
        z_top_exact = sol.z[n, 1] * s(sol.t[end])
        z_err = abs(sol.z[n, end] - z_top_exact) / L
        push!(
            out,
            QoIResult(
                name = "top cell-center trajectory vs z₀·s(t_end)",
                value = sol.z[n, end],
                exact = z_top_exact,
                error = z_err,
                tol = 1.5e-6,  # 3× observed: 4.0e-7 rel. (n=20, reltol=1e-8)
                comment = "node ODE ż = (ṡ/s)z integrated at solver tolerance",
            ),
        )
        return out
    end

    return VerificationCase(
        id = "V-MMS-02",
        title = "ALE free-stream preservation under prescribed mesh motion",
        tags = [:mms, :ale, :gcl, :free_stream],
        references = [
            "Thomas & Lombard, AIAA J. 17 (1979) — geometric conservation law",
            "docs/guide/technical_reference/18_verification-validation.md §18.2 (MMS driver, ALE cases)",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        solve_fn = solve_fn,
        qois = qois,
        ci = CIConfig(n_cells = 20, budget_s = 30.0),
        report = (
            qoi_curves = (:thickness_history,),
            # The free-stream drifts are exactly zero (the QoIs), so the one
            # meaningful overlay is the prescribed-motion exactness check: the
            # top node against the exact uniform-stretch trajectory L·s(t).
            curves = sol -> [
                ReportCurve(
                    name = :thickness_history,
                    title = "V-MMS-02: top-node trajectory under prescribed " *
                            "uniform stretch (30% compression)",
                    xlabel = "t (s)",
                    ylabel = "top node z (m)",
                    series = [
                        CurveSeries(
                            label = "top node (node ODE)",
                            x = sol.t,
                            y = [mms_top_node(sol, j) for j in eachindex(sol.t)],
                            exact = t -> L * s(t),
                        ),
                    ],
                ),
            ],
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vmms02)
