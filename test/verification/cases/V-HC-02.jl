# =============================================================================
# V-HC-02 — Finite slab, convective BC, Biot sweep (FDS heat_conduction_a–d)
# =============================================================================
#
# 0.1 m slab initially at Tᵢ = 300 K, convective front (T∞ = 500 K), adiabatic
# back, Bi = hL/k ∈ {0.1, 1, 10, 100} — four regimes from one definition:
# Bi = 0.1 is lumped-like (BC-dominated), Bi = 100 is surface-Dirichlet-like
# (interior-dominated). Exact eigenfunction series with eigencondition
# λₙ tan λₙ = Bi (exact/slab_series.jl).
#
# Mesh mapping: adiabatic face is z = 0, so the series coordinate is x̃ = z/L.
# Comparison times are Fo = αt/L² ∈ {0.05, 0.1, 0.2, 0.5}, all well above the
# 60-root series validity floor (Fo ≳ 1.1e-3).
#
# Convergence spot check runs on the Bi = 10 member only (the mixed regime);
# the other three assert QoIs at the CI mesh.

function _vhc02_case(label::String, Bi::Float64; with_convergence::Bool)
    ρ, c, k = 1000.0, 1000.0, 0.5
    α = k / (ρ * c)
    L = 0.1
    T_i, T_inf = 300.0, 500.0
    h = Bi * k / L
    Fos = [0.05, 0.1, 0.2, 0.5]
    times = [Fo * L^2 / α for Fo in Fos]
    λs = biot_eigenvalues(Bi, 60)

    T_exact(x̃, Fo) = T_inf + (T_i - T_inf) * slab_convective_theta(x̃, Fo, λs)

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = ρ, c = c, k = k)
        mat = Material(name = :VHC02, components = (solid,), reactions = ())
        mesh = create_uniform_mesh(L, n_cells, Val(1); T_initial = T_i, ξ_initial = (ρ,))
        bc_set = BoundaryConditionSet(
            thermal = (top = ConvectiveBC(h = h, T_inf = T_inf), bottom = AdiabaticBC()),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_i,
            ξ_initial = [z -> ρ],
            tspan = (0.0, times[end]),
            output_times = vcat(0.0, times),
        )
    end

    function qois(sol)
        out = QoIResult[]
        z = get_cell_positions(sol)
        for (j, Fo) in enumerate(Fos)
            Tn = @view sol.T[:, j+1]
            err = maximum(abs(Tn[i] - T_exact(z[i] / L, Fo)) for i in eachindex(z))
            push!(
                out,
                QoIResult(
                    name = "L∞ T-profile error, Fo=$(Fo)",
                    value = err,
                    error = err,
                    tol = 0.06,  # 3× worst across the sweep: 0.020 K (Bi=100, Fo=0.05; n=80, 2026-07-06)
                    comment = "observed max per Bi: 2.3e-4 (0.1), 2.4e-3 (1), 1.8e-2 (10), 2.0e-2 (100) K",
                ),
            )
            Ts_num = sol.surface_T[j+1]
            Ts_exact = T_exact(1.0, Fo)
            push!(
                out,
                QoIResult(
                    name = "surface temperature, Fo=$(Fo)",
                    value = Ts_num,
                    exact = Ts_exact,
                    error = abs(Ts_num - Ts_exact),
                    tol = 0.015,  # 3× worst across the sweep: 4.9e-3 K (Bi=10, Fo=0.05; n=80, 2026-07-06)
                    comment = "front face from the surface Newton solve vs series at x̃=1; " *
                              "observed max per Bi: 9.2e-4 (0.1), 4.4e-3 (1), 4.9e-3 (10), 9.8e-4 (100) K",
                ),
            )
        end
        return out
    end

    # L2 profile error at Fo = 0.2 (output index 4).
    function l2_error_fo02(sol)
        z = get_cell_positions(sol)
        Tn = @view sol.T[:, 4]
        Δz = L / length(z)
        return sqrt(sum((Tn[i] - T_exact(z[i] / L, 0.2))^2 * Δz for i in eachindex(z)) / L)
    end

    return VerificationCase(
        id = "V-HC-02$(label)",
        title = "finite slab, convective BC, Bi = $(Bi)",
        tags = [:heat_conduction, :convective_bc, :biot_sweep],
        references = [
            "Carslaw & Jaeger (1959) §3.11 (eigenfunction series)",
            "FDS Verification Guide `heat_conduction_a–d`",
        ],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        qois = qois,
        ci = CIConfig(
            n_cells = 80,
            budget_s = 20.0,
            convergence = if with_convergence
                ConvergenceSpec(
                    # observed rates 2.0005, 2.0001 at Bi=10 (2026-07-06)
                    n_cells = [40, 80, 160],
                    error_fn = l2_error_fo02,
                    expected_order = 2.0,
                    min_order = 1.5,
                    min_mean_order = 1.8,
                    solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
                )
            else
                nothing
            end,
        ),
        report = (
            n_cells_ladder = [20, 40, 80, 160, 320],
            profile_times = times,
            qoi_curves = (:T_profiles, :surface_T_history, :back_T_history),
        ),
    )
end

for (label, Bi) in (("a", 0.1), ("b", 1.0), ("c", 10.0), ("d", 100.0))
    push!(
        VERIFICATION_CASE_BUILDERS,
        () -> _vhc02_case(label, Bi; with_convergence = Bi == 10.0),
    )
end
