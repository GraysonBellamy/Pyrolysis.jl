# =============================================================================
# V-HC-05 — Temperature-dependent k(T): steady Kirchhoff-exact profile
# =============================================================================
#
# Steady conduction across a slab with k(T) = k₀(1 + βT), fixed T at both
# faces. The Kirchhoff transform gives the exact flux and implicit profile:
#
#   q = (k₀/L)·[(T_top − T_bot) + β(T_top² − T_bot²)/2]
#   z(T) = k₀[(T − T_bot) + β(T² − T_bot²)/2] / q          (bottom cold face)
#
# The implicit form is compared invert-free: evaluate z(T_num[i]) and compare
# to the cell position. k varies 2.25× across 300–800 K with β = 0.01 (the
# linear form caps the ratio at T_top/T_bot; 2.25× is char-vs-virgin scale).
#
# What is verified: the T-dependent property path (LinearProperty) and the
# distance-weighted harmonic-mean face interpolation of cell-centered k(T)
# under strong nonlinearity. Note k(T) is evaluated at cell centers and
# harmonic-averaged to faces — NOT the integrated Kirchhoff face conductance —
# so unlike V-HC-04 this is a genuine discretization-error measurement, with
# a convergence check. Transient T-dependent {k,c} has no exact solution
# (that regime is B-02, code-to-code).

function _verification_case_vhc05()
    ρ, c = 1000.0, 1000.0
    k0, β = 0.2, 0.01
    L = 0.02
    T_bot, T_top = 300.0, 800.0
    q_exact = (k0 / L) * ((T_top - T_bot) + β * (T_top^2 - T_bot^2) / 2)
    z_of_T(T) = k0 * ((T - T_bot) + β * (T^2 - T_bot^2) / 2) / q_exact
    α_min = k0 * (1 + β * T_bot) / (ρ * c)
    t_end = 10 * L^2 / α_min

    function build(n_cells)
        solid = SolidComponent(:solid, ρ = ρ, c = c, k = LinearProperty(k0, k0 * β))
        mat = Material(name = :VHC05, components = (solid,), reactions = ())
        mesh = create_uniform_mesh(L, n_cells, Val(1); T_initial = T_bot, ξ_initial = (ρ,))
        bc_set = BoundaryConditionSet(
            thermal = (top = DirichletThermalBC(T_top), bottom = DirichletThermalBC(T_bot)),
        )
        return PyrolysisProblem(
            mesh = mesh,
            material = mat,
            bc_set = bc_set,
            T_initial = z -> T_bot,
            ξ_initial = [z -> ρ],
            tspan = (0.0, t_end),
            output_times = [0.0, t_end / 2, t_end],
        )
    end

    function qois(sol)
        out = QoIResult[]
        z = cell_positions(sol)
        Tn = @view sol.T[:, end]

        err = maximum(abs(z_of_T(Tn[i]) - z[i]) for i in eachindex(z)) / L
        push!(
            out,
            QoIResult(
                name = "L∞ Kirchhoff-profile error (z(T) map, /L)",
                value = err,
                error = err,
                tol = 1.2e-4,  # 3× observed 3.9e-5 (n=80)
                comment = "implicit profile compared invert-free; n=80 uniform",
            ),
        )

        # Exact flux: reconstructed from the Kirchhoff integral between the
        # two cells adjacent to the domain center (pure interior estimate,
        # avoids the O(Δz) boundary half-cell k-evaluation bias).
        n = length(z)
        i1, i2 = n ÷ 2, n ÷ 2 + 1
        q_num = k0 * ((Tn[i2] - Tn[i1]) + β * (Tn[i2]^2 - Tn[i1]^2) / 2) / (z[i2] - z[i1])
        push!(
            out,
            QoIResult(
                name = "steady flux vs Kirchhoff-exact",
                value = q_num,
                exact = q_exact,
                error = abs(q_num - q_exact) / q_exact,
                tol = 6e-5,  # 3× observed 1.8e-5
                comment = "q = (k₀/L)[ΔT + βΔ(T²)/2] = $(round(q_exact, digits = 1)) W/m²",
            ),
        )

        steadiness = maximum(abs.(sol.T[:, end] .- sol.T[:, end-1]))
        push!(
            out,
            QoIResult(
                name = "steadiness (last half-interval drift)",
                value = steadiness,
                error = steadiness,
                tol = 1e-6,  # observed 4.7e-8 K
                comment = "t_end = 10·L²/α_min",
            ),
        )
        return out
    end

    function profile_error(sol)
        z = cell_positions(sol)
        Tn = @view sol.T[:, end]
        return maximum(abs(z_of_T(Tn[i]) - z[i]) for i in eachindex(z)) / L
    end

    return VerificationCase(
        id = "V-HC-05",
        title = "steady Kirchhoff profile with k(T) = k₀(1+βT)",
        tags = [:heat_conduction, :temperature_dependent_k, :kirchhoff, :steady],
        references = ["Kirchhoff transform (exact steady flux/profile), standard texts"],
        build = build,
        solve_kwargs = (abstol = 1e-10, reltol = 1e-8),
        qois = qois,
        ci = CIConfig(
            n_cells = 80,
            budget_s = 20.0,
            convergence = ConvergenceSpec(
                # observed rates 1.965, 1.982
                n_cells = [40, 80, 160],
                error_fn = profile_error,
                expected_order = 2.0,
                min_order = 1.5,
                min_mean_order = 1.8,
                solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
            ),
        ),
        report = (
            n_cells_ladder = [20, 40, 80, 160, 320],
            profile_times = [t_end],
            qoi_curves = (:T_profile_steady, :flux_history),
            solve_kwargs = (abstol = 1e-12, reltol = 1e-10),
            curves = sol -> begin
                z = cell_positions(sol)
                n = length(z)
                i1, i2 = n ÷ 2, n ÷ 2 + 1
                # invert the implicit Kirchhoff profile z(T) for the overlay
                T_of_z(zz) = bisect_root(T -> z_of_T(T) - zz, T_bot, T_top)
                q_interior(j) =
                    k0 * (
                        (sol.T[i2, j] - sol.T[i1, j]) +
                        β * (sol.T[i2, j]^2 - sol.T[i1, j]^2) / 2
                    ) / (z[i2] - z[i1])
                [
                    ReportCurve(
                        name = :T_profile_steady,
                        title = "V-HC-05: steady Kirchhoff profile, k(T) = k₀(1+βT)",
                        xlabel = "z (m)",
                        ylabel = "T (K)",
                        series = [
                            CurveSeries(
                                label = "steady solver profile",
                                x = z,
                                y = sol.T[:, end],
                                exact = T_of_z,
                            ),
                        ],
                    ),
                    ReportCurve(
                        name = :flux_history,
                        title = "V-HC-05: interior flux approach to Kirchhoff-exact",
                        xlabel = "t (s)",
                        ylabel = "q″ (W/m²)",
                        series = [
                            CurveSeries(
                                label = "center-pair Kirchhoff flux",
                                x = sol.t,
                                y = [q_interior(j) for j in eachindex(sol.t)],
                                exact = t -> q_exact,
                            ),
                        ],
                    ),
                ]
            end,
            error_norms = sol -> begin
                z = cell_positions(sol)
                e = [z_of_T(sol.T[i, end]) - z[i] for i in eachindex(z)]
                (L2 = sqrt(sum(abs2, e) / length(e)) / L, Linf = maximum(abs, e) / L)
            end,
        ),
    )
end

push!(VERIFICATION_CASE_BUILDERS, _verification_case_vhc05)
